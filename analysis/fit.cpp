/*
 * =====================================================================================
 *
 *       Filename:  fit.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  19.10.2016 23:27:58
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */

#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TMath.h>

#include <boost/format.hpp>

#include <RooSimultaneous.h>
#include <RooFitResult.h>
#include <RooChi2Var.h>
#include <RooNLLVar.h>
#include <RooRealVar.h>
#include <RooPolynomial.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooCategory.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooBukinPdf.h>
using namespace RooFit;

#include "RooMcbPdf.h"
#include "fit.h"

bool OPT_NOBGSLOPE=false; //no slope for the background
bool OPT_NOBG=false; //no background
bool OPT_NOGAUSRAD=false; //no gaus rad
std::string OPT_PARAM_CONFIG_FILE=""; 

void fit(TH1 * his)
{
  std::list<TH1*> hlst = {his};
  fit(hlst);
}

void fit(TH1 * hisKK, TH1 * hisUU)
{
  hisKK->SetName("KK");
  hisUU->SetName("UU");
  std::list<TH1*> hlst = {hisKK,hisUU};
  fit(hlst);
}




void fit(std::list<TH1*> & hlst)
{
  //define maximum and minimum recoil mass
  double Mmin  = (*std::max_element(std::begin(hlst), std::end(hlst), 
      [](const TH1 * h1, const TH1 *h2) { return h1->GetXaxis()->GetXmin() < h2->GetXaxis()->GetXmin(); } ))->GetXaxis()->GetXmin();

  double Mmax  = (*std::min_element(std::begin(hlst), std::end(hlst), 
      [](const TH1 * h1, const TH1 *h2) { return h1->GetXaxis()->GetXmax() < h2->GetXaxis()->GetXmax(); } ))->GetXaxis()->GetXmax();

  std::list < std::string > name_lst;
  std::map<std::string,TH1*> hisMap;
  for(auto h : hlst) 
  {
    name_lst.push_back(h->GetName());
    hisMap[h->GetName()] = h;
  }
  //i decided to make own variable for mmrec
  std::map<std::string, RooRealVar*> MMrec;

  boost::format mrec_title_fmt("M_{rec}^{%s}(#pi^{+}#pi^{-}) - 3097 MeV");
  for(auto name : name_lst)
  {
    MMrec[name] = new RooRealVar(
        ("Mrec"+name).c_str(), 
         (mrec_title_fmt % name).str().c_str(), 
         hisMap[name]->GetXaxis()->GetXmin(),
        hisMap[name]->GetXaxis()->GetXmax(), "MeV"
        );
    MMrec[name]->setBins(hisMap[name]->GetNbinsX());
  }



  RooRealVar  Mrec("Mrec",(mrec_title_fmt % "").str().c_str(), Mmin,Mmax, "MeV");
  RooRealVar  mean( "mean",  "mean",   -0.1 + 0.5*(Mmin+Mmax) ,  Mmin,  Mmax, "MeV") ;
  RooRealVar sigma("sigma","sigma",1.4,0,10, "MeV") ;
  RooRealVar n1("n1","n1", 3,  1,100) ;
  RooRealVar n2("n2","n2", 2,  1,100) ;

  std::vector<RooRealVar> staple =
  {
    RooRealVar("L1", "Left Gaus range"   , 2        , "MeV"),
    RooRealVar("L2", "Left Exp range"    , 1        , 0      , (Mmax-Mmin)*0.5 , "MeV"),
    RooRealVar("L3", "Left Power range"  , 10       , 0      , (Mmax-Mmin)*0.5 , "MeV"),
    RooRealVar("L4", "Left Exp range 2"  , 20       , 0      , (Mmax-Mmin)*0.5 , "MeV"),
    RooRealVar("R1", "Right Gaus range"  , 2        , "MeV"),
    RooRealVar("R2", "Right Exp range"   , 1        , 0      , (Mmax-Mmin)*0.5 , "MeV"), 
    RooRealVar("R3", "Right Power range" , 30       , 0      , (Mmax-Mmin)*0.5 , "MeV"), 
    RooRealVar("R4", "Right Exp range 2" , 1        , 0      , (Mmax-Mmin)*0.5 , "MeV")
  };

  std::map<std::string, RooAbsPdf*> McbPdfMap;
  std::map<std::string, RooAbsPdf*> SignalPdfMap;
  for(auto name : name_lst)
  {
    McbPdfMap[name]=new RooMcb2Pdf(("Mcb2Pdf_"+name).c_str(), ("My modified Crystal Bal function for " + name).c_str(),
        *MMrec[name], mean,sigma,staple, n1, n2, hisMap[name]->GetXaxis()->GetXmin(), hisMap[name]->GetXaxis()->GetXmax());
  }
	//RooMcb2Pdf *mcbPdf=0;
  //mcbPdf = new RooMcb2Pdf("mcb2","mcb2",Mrec,mean,sigma, staple, n1,n2,Mmin,Mmax);


  std::vector<RooRealVar*> meanRad;
  if(!OPT_NOGAUSRAD) meanRad.resize(2);
  std::vector<RooRealVar*> sigmaRad(meanRad.size());
  std::vector<RooGaussian*> radPdf(meanRad.size());
  std::vector<RooRealVar*> radFrac(meanRad.size());

  std::map<std::string, std::vector<RooGaussian*>> radPdfMap;

  //create parameters for radiative gaussian (mean and sigma)
  for(int i=0;i<meanRad.size();i++)
  {
    std::string istr = to_string(i);
    meanRad[i] = new RooRealVar(("mean_rad" + istr).c_str(), ("mean_rad" + istr).c_str(), 20,
        (Mmax+Mmin)*0.5+5, Mmax, "MeV");
    sigmaRad[i] = new RooRealVar(("sigma_rad" + istr).c_str(), ("sigma_rad" + istr).c_str(), 10,
        2, 100, "MeV");
  }

  //name create the SignalPdfs for each channel
  for(auto name : name_lst)
  {
    RooArgList PdfList;
    RooArgList RadFracList;
    radPdfMap[name].resize(meanRad.size());

    for(int i=0;i<meanRad.size();i++)
    {
      std::string istr = to_string(i);
      auto pdf = new RooGaussian(("gaus_radPdf" + name + istr).c_str(),
          ("Radiative gauss " + istr + " for " + name).c_str(), *MMrec[name], *meanRad[i], *sigmaRad[i]);
      radPdfMap[name][i] = pdf;
      PdfList.add(*pdf);
      RadFracList.add(*radFrac[i]);
    }
    //PdfList.add(*mcbPdf);
    PdfList.add(*McbPdfMap[name]);
    //RooAbsPdf * signalPdf =  new RooAddPdf("signalPdf","Signal model", PdfList, RadFracList);
    SignalPdfMap[name] = new RooAddPdf(("signalPdf"+name).c_str(),("Signal model for "+name).c_str(), PdfList, RadFracList);
  }

  //Mixing with the background
  std::map<std::string, RooAbsPdf * > SamplePdf; //here will be Pdf with signal and background mixed

  std::map<std::string, RooRealVar * > Nsig; //number of signal events
  std::map<std::string, RooRealVar * > Nbg; //number of background events
  std::map<std::string, RooPolynomial*> bgPdf; //background pdf function
  std::map<std::string, RooRealVar*> bg_c1; //slope of the background
  std::map<std::string, long> Nhis;
  std::map<std::string, RooPlot*> frame;

 	RooCategory sample("sample","sample");
  RooArgList MrecArgList;
  std::map<std::string, RooDataHist *> dataMap;

  for ( auto his : hlst)
  {
    std::string name = his->GetName();
    Nhis[name] = his->GetEntries();

    bg_c1[name] = new RooRealVar(("bg"+name+"_c1").c_str(),("background slope for "+name).c_str(),0,- 1.0/Mmax, - 1.0/Mmin);
    if(!OPT_NOBGSLOPE)
    {
      bgPdf[name] = new RooPolynomial(("bg"+name).c_str(), ( name + " background pdf").c_str(), *MMrec[name], *bg_c1[name]);	
    }
    else
    {
      bgPdf[name] = new RooPolynomial(("bg"+name).c_str(), ( name + " background pdf").c_str(), *MMrec[name]);	
    }
    Nsig[name] = new RooRealVar(("Nsig"+name).c_str(), ("Number of " + name +" events").c_str(), Nhis[name],0, Nhis[name]*100);
    Nbg[name] = new RooRealVar(("Nbg"+name).c_str(), ("Number of background events for " + name + " channel").c_str(), 0, 0, Nhis[name]*100);
    if(!OPT_NOBG)
    {
      SamplePdf[name] = new RooAddPdf((name+"Pdf").c_str(), (name+" signal + background p.d.f.").c_str(), 
          //RooArgList(*bgPdf[name], *signalPdf), 
          RooArgList(*bgPdf[name], *SignalPdfMap[name]), 
          RooArgList(*Nbg[name],  *Nsig[name]));
    }
    else
    {
      SamplePdf[name] = new RooAddPdf((name+"Pdf").c_str(), (name+" signal + background p.d.f.").c_str(), 
          RooArgList(*SignalPdfMap[name]), RooArgList(*Nsig[name]));
          //RooArgList(*signalPdf), RooArgList(*Nsig[name]));
    }
    dataMap[name]=new RooDataHist(name.c_str(), name.c_str(), *MMrec[name], Import(*hisMap[name]));
    sample.defineType(name.c_str());
    MrecArgList.add(*MMrec[name]);
    frame[name] = MMrec[name]->frame(Title(("#pi^{+}#pi^{-} recoil mass for " + name).c_str()));
  }

  RooDataHist * data = new RooDataHist("combData","Combined data", MrecArgList, sample, dataMap);
  //RooDataHist * data = new RooDataHist("combData","Combined data", MrecArgList, sample, hisMap);
  RooSimultaneous simPdf("simPdf","simultaneous pdf",SamplePdf,sample) ;

  auto p = simPdf.getParameters(*MMrec[name_lst.front()]);
  if(OPT_PARAM_CONFIG_FILE!="")
  {
    p->readFromFile(OPT_PARAM_CONFIG_FILE.c_str());
    p->Print("v");
  }

	auto theFitResult = simPdf.fitTo(*data, Extended(), Strategy(2), Minos(), Save());

  p->writeToFile("tmp_fit_result.txt");

	RooChi2Var chi2Var("chi2", "chi2", simPdf, *data);
  //number of free parameters
  int nfp = theFitResult->floatParsFinal().getSize() ;
  //chi square
	double chi2 =  chi2Var.getValV();
  //number of degree of freedom
  int ndf = -nfp;
  for(auto his : hisMap) ndf+=his.second->GetNbinsX();
  //normalized chi2
  double chi2ndf = chi2/ndf;
  //probability for this hypotisis to get chi2 higher then for this data
  double chi2prob = TMath::Prob(chi2,ndf);



	RooNLLVar nll("nll","nll",simPdf,*data) ;

  RooPlot* frame_sigma = sigma.frame(Range(0.8, 5.0)) ;
  nll.plotOn(frame_sigma,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kRed)) ;

  RooPlot* frame_mean = mean.frame(Range(-2.0, 2.0)) ;
  nll.plotOn(frame_mean,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kRed)) ;

  std::vector<RooPlot*> frame_staple(staple.size());
  std::vector<int> colors = {kBlack, kBlue,kCyan+2,kGreen+2, kGreen, kRed, kMagenta, kRed+2}; 
  std::vector<int> line_styles = {kSolid, kSolid, kDashed, kDotted, kSolid, kSolid, kDashed, kDotted};
  for(int i=1;i<4;i++)
  {
    frame_staple[i] = staple[i].frame();
    frame_staple[i+4] = staple[i+4].frame();
    nll.plotOn(frame_staple[i],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(colors[i]), LineStyle(line_styles[i])) ;
    nll.plotOn(frame_staple[i+4],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(colors[i+4]),LineStyle(line_styles[i+4])) ;
  }


  RooPlot* frame_n1 = n1.frame(Range(1,10)) ;
  nll.plotOn(frame_n1,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kBlue)) ;

  RooPlot* frame_n2 = n2.frame(Range(1, 10)) ;
  nll.plotOn(frame_n2,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kRed));


	TCanvas * cnll =new TCanvas;
	cnll->Divide(2, 2);
	cnll->cd(1);
	frame_mean->Draw();
	cnll->cd(2);
	frame_sigma->Draw();
	cnll->cd(3);
	frame_n1->Draw();
	frame_n2->Draw("same");
	cnll->cd(4);
	frame_staple[1]->Draw();
	frame_staple[5]->Draw("same");
	frame_staple[2]->Draw("same");
	frame_staple[6]->Draw("same");
	frame_staple[3]->Draw("same");
	frame_staple[7]->Draw("same");

  std::vector<RooPlot*> Nframe;
  TCanvas * Ncanvas = new TCanvas("Ncanvas","Log likelihood for number of events");
  Ncanvas->Divide(2,hlst.size());

  for(auto his : hlst)
  {
    std::string name = his->GetName();
    Nframe.push_back(Nsig[name]->frame(Range(Nsig[name]->getValV()*0.8,Nsig[name]->getValV()*1.2)));
    Nframe.push_back( Nbg[name]->frame(Range(Nbg[name]->getValV()*0.8, Nbg[name]->getValV()*1.2)));
  }
  for(int i=0;i<Nframe.size();i++)
  {
    nll.plotOn(Nframe[i],PrintEvalErrors(-1),ShiftToZero(), EvalErrorValue(nll.getVal()+10),LineColor(kRed)) ;
    Ncanvas->cd(i+1);
    Nframe[i]->Draw();
  }

  //RooPlot * Frame = Mrec.frame(Title("#pi^{+}#pi^{-} recoil mass"));
  RooArgSet viewArgSet;
  int i=0;
  for( auto & his : hlst)
  {
    std::string name = his->GetName();
    data->plotOn(frame[name], XErrorSize(0), MarkerSize(0.5),  Cut(("sample==sample::"+name).c_str()),LineColor(colors[i]), MarkerColor(colors[i]));
    simPdf.plotOn(frame[name], Slice(sample,name.c_str()),ProjWData(sample,*data), LineWidth(1),LineColor(colors[i]), MarkerColor(colors[i]));
    simPdf.plotOn(frame[name], Slice(sample,name.c_str()),ProjWData(sample,*data), Components(*bgPdf[name]),LineStyle(kDashed),LineWidth(1),LineColor(colors[i]), MarkerColor(colors[i]));
    frame[name]->SetMinimum(0.1);
    frame[name]->GetYaxis()->SetTitleOffset(1.6) ; 
    viewArgSet.add(*Nsig[name]);
    viewArgSet.add(*Nbg[name]);
    i++;
  }
	//Frame->SetMinimum(0.1);


  TCanvas* c = new TCanvas("mcb","fit") ;
  gPad->SetLeftMargin(0.15) ; 
  gPad->SetLogy();
  //now find frame with maximum number of events
  std::string most_entries_his_name;
  long long tmpN = 0;
  for(auto name : name_lst)
  {
    if(Nhis[name]>tmpN)
    {
      tmpN=Nhis[name];
      most_entries_his_name = name;
    }
  }

  auto chi_fmt = boost::format("#chi^{2}/ndf = %.1f/%d = %.2f") % chi2 % ndf % chi2ndf;
  auto prob_fmt = boost::format("prob = %.1f%%") % chi2prob;
  simPdf.paramOn(frame[most_entries_his_name], Label(chi_fmt.str().c_str()), Parameters(viewArgSet));

  //now draw this frame first
  frame[most_entries_his_name]->Draw();
  //other should be draw with the option "same"
  for(auto name :name_lst)
  {
    if(name != most_entries_his_name)
    {
      frame[name]->Draw("same");
    }
  }


  //Frame->GetYaxis()->SetTitleOffset(1.6) ; 
  //
  //Frame->Draw() ;


  mean.Print();
  sigma.Print();
  for(auto & his : hlst)
  {
    std::string name = his->GetName();
    Nsig[name]->Print();
    Nbg[name]->Print();
  }
  std::cout << boost::format("chi2/ndf = %f/(%d-%d) = %f,  prob = %f") % chi2 % ndf % nfp % chi2ndf %  chi2prob<< std::endl;

  //No print the ration of signal event to number of signal event for last data
  //sample
  for(auto his : hlst)
  {
    if(his == hlst.back()) continue;
    std::string name = his->GetName();
    std::string name0 = hlst.back()->GetName();
    auto  N = Nsig[name];
    auto  N0 = Nsig[name0];
    double n   =  N -> getValV();
    double dn  =  N -> getError();
    double n0  = N0 -> getValV();
    double dn0 = N0 -> getError();
    double theRatio = n/n0;
    double theRelError = sqrt( pow(dn/n, 2)  +  pow(dn0/n0, 2) );
    double theError = theRatio * theRelError;
    boost::format fmt("%s/%s = (%6.3f ± %-4.3f) * 1e-3 ( ± %.2f%%)");
    cout << fmt % N->GetName() % N0->GetName() % (theRatio*1000) %  (theError*1000) %  (theRelError*100)  <<  endl;
  }
}

