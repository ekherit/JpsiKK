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

  RooRealVar Mrec("Mrec","M_{rec}(#pi^{+}#pi^{-}) - 3097 ", Mmin,Mmax, "MeV");
  RooRealVar  mean( "mean",  "mean",   -0.1 + 0.5*(Mmin+Mmax) ,  Mmin,  Mmax, "MeV") ;
  RooRealVar sigma("sigma","sigma",1.4,0,10, "MeV") ;
  RooRealVar n1("n1","n1", 3,  1,100) ;
  RooRealVar n2("n2","n2", 2,  1,100) ;

  std::vector<RooRealVar> staple =
  {
    RooRealVar("L1" , "Left Gaus range"   , 2        , "MeV"),
    RooRealVar("L2" , "Left Exp range"    , 0.853    , 0      , (Mmax-Mmin)*0.5 , "MeV"),
    RooRealVar("L3" , "Left Power range"  , 12.8861  , 0      , (Mmax-Mmin)*0.5 , "MeV"),
    RooRealVar("L4" , "Left Exp range 2"  , 22.7831  , 0      , (Mmax-Mmin)*0.5 , "MeV"),
    RooRealVar("R1" , "Right Gaus range"  , 2        , "MeV"),
    RooRealVar("R2" , "Right Exp range"   , 2.6e-7   , 0      , (Mmax-Mmin)*0.5 , "MeV"), 
    RooRealVar("R3" , "Right Power range" , 31.9138  , 0      , (Mmax-Mmin)*0.5 , "MeV"), 
    RooRealVar("R4" , "Right Exp range 2" , 0.125333 , 0      , (Mmax-Mmin)*0.5 , "MeV")
  };

	RooMcb2Pdf *mcbPdf=0;
  mcbPdf = new RooMcb2Pdf("mcb2","mcb2",Mrec,mean,sigma,
      staple,
      n1,n2,Mmin,Mmax);



	RooRealVar radMean("radMean", "rad mean",20,  1, 40, "MeV");
	RooRealVar radSigma("radSigma", "rad sigma", 10,  5, 30, "MeV");
	RooGaussian radGausPdf("smals_gauss", "small_gauss",  Mrec,  radMean,  radSigma);

	RooAbsPdf * signalPdf =0;
	signalPdf= mcbPdf;

	RooRealVar fRad("fRad", "Fraction of rad", 1e-3, 0, 0.1);
	RooAddPdf  signalRadPdf("signalPdf", "Signal with radiative photons",  RooArgList(radGausPdf,  *mcbPdf),  fRad);

  std::map<std::string, RooRealVar * > Nsig;
  std::map<std::string, RooRealVar * > Nbg;
  std::map<std::string, RooAbsPdf * > SamplePdf;
  std::map<std::string, RooPolynomial*> bgPdf;
  std::map<std::string, RooRealVar*> bg_c1;
  std::map<std::string, long> Nhis;
  std::map<std::string, RooDataHist*> dataMap;
  std::map<std::string,TH1*> hisMap;

  std::map<std::string, RooPlot*> frame;

 	RooCategory sample("sample","sample") ;

  for ( auto his : hlst)
  {
    std::string name = his->GetName();
    hisMap[name] = his;
    Nhis[name] = his->GetEntries();

    bg_c1[name] = new RooRealVar(("bg"+name+"_c1").c_str(),("background slope for "+name).c_str(),0,- 1.0/Mmax, - 1.0/Mmin);
    bgPdf[name] = new RooPolynomial(("bg"+name).c_str(), ( name + " background pdf").c_str(), Mrec, *bg_c1[name]);	
    Nsig[name] = new RooRealVar(("Nsig"+name).c_str(), ("Number of " + name +" events").c_str(), Nhis[name],0, Nhis[name]*100);
    Nbg[name] = new RooRealVar(("Nbg"+name).c_str(), ("Number of background events for " + name + " channel").c_str(), 0, 0, Nhis[name]*100);
    SamplePdf[name] = new RooAddPdf((name+"Pdf").c_str(), (name+" signal + background p.d.f.").c_str(), RooArgList(*bgPdf[name], *signalPdf), RooArgList(*Nbg[name],  *Nsig[name]));
    dataMap[name]=new RooDataHist(name.c_str(), name.c_str(), Mrec, Import(*hisMap[name]));
    sample.defineType(name.c_str());
    frame[name] = Mrec.frame(Title("#pi^{+}#pi^{-} recoil mass"));
  }

  RooDataHist * data = new RooDataHist("combData","Combined data", Mrec, sample, dataMap);
  RooSimultaneous simPdf("simPdf","simultaneous pdf",SamplePdf,sample) ;

	simPdf.fitTo(*data, Extended(), Strategy(2), Minos());

	RooChi2Var chi2Var("chi2", "chi2", simPdf, *data);
	cout << "chi2 = " << chi2Var.getVal() << endl;
  RooArgSet args(Mrec);
  std::cout << "exp. Events = " << simPdf.expectedEvents(&args) << std::endl;



	RooNLLVar nll("nll","nll",simPdf,*data) ;

  RooPlot* frame_sigma = sigma.frame(Range(0.2, 2.0)) ;
  nll.plotOn(frame_sigma,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kRed)) ;

  RooPlot* frame_mean = mean.frame(Range(-2.0, 2.0)) ;
  nll.plotOn(frame_mean,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kRed)) ;

  std::vector<RooPlot*> frame_staple(staple.size());
  std::vector<int> colors = {kBlack,kBlue,kCyan+2,kGreen+2, kBlack, kRed, kMagenta+2, kRed+2}; 
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

  for( auto & his : hlst)
  {
    std::string name = his->GetName();
    data->plotOn(frame[name], MarkerSize(0.5),  Cut(("sample==sample::"+name).c_str()));
    simPdf.plotOn(frame[name], Slice(sample,name.c_str()),ProjWData(sample,*data), LineWidth(1));
    simPdf.plotOn(frame[name], Slice(sample,name.c_str()),ProjWData(sample,*data), Components(*bgPdf[name]),LineStyle(kDashed),LineWidth(1));
    frame[name]->SetMinimum(0.1);
  }

  TCanvas* c = new TCanvas("mcb","fit") ;
	c->Divide(1, hlst.size());
  int i=1;
  for(auto & his : hlst)
  {
    std::string name = his->GetName();
    c->cd(i);
    gPad->SetLeftMargin(0.15) ; 
    gPad->SetLogy();
    frame[name]->GetYaxis()->SetTitleOffset(1.6) ; 
    frame[name]->Draw() ;
    i++;
  }

  mean.Print();
  sigma.Print();
  for(auto & his : hlst)
  {
    std::string name = his->GetName();
    Nsig[name]->Print();
    Nbg[name]->Print();
  }

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
    cout << fmt % name % name0 % (theRatio*1000) %  (theError*1000) %  (theRelError*100)  <<  endl;
  }

}

/*  
void fit2(std::list<TH1*> & hlst)
{
  //define maximum and minimum recoil mass
  double Mmin  = (*std::max_element(std::begin(hlst), std::end(hlst), 
      [](const TH1 * h1, const TH1 *h2) { return h1->GetXaxis()->GetXmin() < h2->GetXaxis()->GetXmin(); } ))->GetXaxis()->GetXmin();

  double Mmax  = (*std::min_element(std::begin(hlst), std::end(hlst), 
      [](const TH1 * h1, const TH1 *h2) { return h1->GetXaxis()->GetXmax() < h2->GetXaxis()->GetXmax(); } ))->GetXaxis()->GetXmax();

  //define variables
  RooRealVar Mrec("Mrec","M_{rec}(#pi^{+}#pi^{-}) - 3097 ", Mmin,Mmax, "MeV");
  RooRealVar mean( "mean",  "mean",   -0.1 + 0.5*(Mmin+Mmax) ,  Mmin,  Mmax, "MeV") ;
  RooRealVar sigma("sigma","sigma",1.4,0,10, "MeV") ;
  RooRealVar n1("n1","n1", 3,  1,100) ;
  RooRealVar n2("n2","n2", 2,  1,100) ;

  RooRealVar staple1("L1" , "Left Gaus range"   , 2        , "MeV");
  RooRealVar staple2("L2" , "Left Exp range"    , 0.853    , 0      , (Mmax-Mmin)*0.5 , "MeV");
  RooRealVar staple3("L3" , "Left Power range"  , 12.8861  , 0      , (Mmax-Mmin)*0.5 , "MeV");
  RooRealVar staple4("L4" , "Left Exp range 2"  , 22.7831  , 0      , (Mmax-Mmin)*0.5 , "MeV");
  RooRealVar staple5("R1" , "Right Gaus range"  , 2        , "MeV");
  RooRealVar staple6("R2" , "Right Exp range"   , 2.6e-7   , 0      , (Mmax-Mmin)*0.5 , "MeV"); 
  RooRealVar staple7("R3" , "Right Power range" , 31.9138  , 0      , (Mmax-Mmin)*0.5 , "MeV"); 
  RooRealVar staple8("R4" , "Right Exp range 2" , 0.125333 , 0      , (Mmax-Mmin)*0.5 , "MeV");

  //create McbPdf
	RooMcb2Pdf * mcbPdf = new RooMcb2Pdf("mcb2","mcb2",Mrec,mean,sigma,
      staple1,staple2,staple3,staple4,staple5,staple6,staple7,staple8,
      n1,n2,Mmin,Mmax);

	RooRealVar radMean("radMean", "rad mean",20,  1, 40, "MeV");
	RooRealVar radSigma("radSigma", "rad sigma", 10,  5, 30, "MeV");
	RooGaussian radGausPdf("smals_gauss", "small_gauss",  Mrec,  radMean,  radSigma);

	RooAbsPdf * signalPdf =0;
	signalPdf= mcbPdf;

	RooRealVar fRad("fRad", "Fraction of rad", 1e-3, 0, 0.1);
	//RooRealVar fRad("fRad", "Fraction of rad", 0);
	RooAddPdf  signalRadPdf("signalPdf", "Signal with radiative photons",  RooArgList(radGausPdf,  *mcbPdf),  fRad);


  std::map<std::string, RooDataHist*> dataMap;

  std::map<std::string, RooRealVar*> Nsig;
  std::map<std::string, RooRealVar*> Nbg;
  std::map<std::string, RooRealVar*> bg_c1;
  std::map<std::string, RooPolynomial*> bgPdf;
  std::map<std::string, RooAbsPdf*> samplePdf;
  std::map<std::string, RooPlot*> frame;

 	RooCategory sample("sample","sample") ;

  std::cout <<"Before init Nsig Nbg etc..." <<std::endl;
  for(auto & his : hlst)
  {
    std::string name = his->GetName();
    dataMap[name]= new RooDataHist(name.c_str(), name.c_str(), Mrec, Import(*his));

    bg_c1[name] = new RooRealVar(("bg"+name+"_c1").c_str(),("background slope for " + name).c_str(), 0, -10,10);
    bgPdf[name] = new RooPolynomial(("bg"+name).c_str(),("background Pdf for " + name).c_str(),Mrec);

    Nsig[name] = new RooRealVar( ("Nsig"+name).c_str(), ("Number of signal events for " + name).c_str(), his->GetEntries(), 0, his->GetEntries()*100);

    Nbg[name] = new RooRealVar( ("Nbg"+name).c_str(), ("Number of background events for " + name).c_str(), 0, 0, his->GetEntries()*100);

    samplePdf[name] = new RooAddPdf(name.c_str(), ("Pdf for " + name).c_str(), 
        RooArgList(*bgPdf[name], *signalPdf), RooArgList(*Nbg[name],  *Nsig[name]));

    sample.defineType(name.c_str());
    frame[name] = Mrec.frame(Title("#pi^{+}#pi^{-} recoil mass"));
  }
  std::cout <<"After init Nsig Nbg etc..." <<std::endl;

  RooDataHist data("combData","Combined data", Mrec, sample, dataMap);

  std::cout << "Before RooSimultaneous"<<std::endl;
  RooSimultaneous simPdf("simPdf","simultaneous pdf",samplePdf,sample) ;

  
  std::cout << "Before fit" <<std::endl;
	simPdf.fitTo(data, Extended(), Strategy(2), Minos());
  std::cout << "After fit" <<std::endl;

	RooChi2Var chi2Var("chi2", "chi2", simPdf, data);
	cout << "chi2 = " << chi2Var.getVal() << endl;
  RooArgSet args(Mrec);
  std::cout << "exp. Events = " << simPdf.expectedEvents(&args) << std::endl;


  for( auto & his : hlst)
  {
    std::string name = his->GetName();
    data.plotOn(frame[name], MarkerSize(0.5),  Cut(("sample==sample::"+name).c_str()));
    simPdf.plotOn(frame[name], Slice(sample,name.c_str()),ProjWData(sample,data), LineWidth(1));
    simPdf.plotOn(frame[name], Slice(sample,name.c_str()),ProjWData(sample,data), Components(*bgPdf[name]),LineStyle(kDashed),LineWidth(1));
    frame[name]->SetMinimum(0.1);
  }

  mean.Print();
  sigma.Print();
  for(auto & his : hlst)
  {
    std::string name = his->GetName();
    Nsig[name]->Print();
    Nbg[name]->Print();
  }

  TCanvas* c = new TCanvas("mcb","Test for Modified CrystalBall Fit") ;
	c->Divide(1, hlst.size());
  int i=1;
  for(auto & his : hlst)
  {
    std::string name = his->GetName();
    c->cd(i);
    gPad->SetLeftMargin(0.15) ; 
    gPad->SetLogy();
    frame[name]->GetYaxis()->SetTitleOffset(1.6) ; 
    frame[name]->Draw() ;
    i++;
  }
  


	RooNLLVar nll("nll","nll",simPdf,data) ;

  RooPlot* frame_sigma = sigma.frame(Range(0.2, 2.0), Title("-log(L) scan vs sigma")) ;
  nll.plotOn(frame_sigma,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kRed)) ;

  RooPlot* frame_mean = mean.frame(Range(-2.0, 2.0), Title("-log(L) scan vs mean")) ;
  nll.plotOn(frame_mean,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kRed)) ;

  RooPlot* frame_staple2 = staple2.frame(Range(0.,40.),Title("staple2 -log(L)")) ;
  nll.plotOn(frame_staple2,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kBlue)) ;

  RooPlot* frame_staple6 = staple6.frame(Range(0.,40.),Title("staple6 -log(L)")) ;
  nll.plotOn(frame_staple6,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kRed)) ;

  RooPlot* frame_staple3 = staple3.frame(Range(0.,40.),Title("staple3 -log(L)")) ;
  nll.plotOn(frame_staple3,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kBlue)) ;

  RooPlot* frame_staple7 = staple7.frame(Range(0.,40.),Title("staple7 -log(L)")) ;
  nll.plotOn(frame_staple7,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kRed)) ;

  RooPlot* frame_staple4 = staple4.frame(Range(0.,40.),Title("staple4 -log(L)")) ;
  nll.plotOn(frame_staple4,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kBlue)) ;

  RooPlot* frame_staple8 = staple8.frame(Range(0.,40.),Title("staple8 -log(L)")) ;
  nll.plotOn(frame_staple8,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kRed)) ;

  RooPlot* frame_n1 = n1.frame(Title("n1 -log(L)"), Range(1, 10)) ;
  nll.plotOn(frame_n1,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kBlue)) ;

  RooPlot* frame_n2 = n2.frame(Title("n2 -log(L)"), Range(1, 10)) ;
  nll.plotOn(frame_n2,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kRed)) ;


	TCanvas * cnll =new TCanvas;
	cnll->Divide(3, 3);
	cnll->cd(1);
	frame_mean->Draw();
	cnll->cd(2);
	frame_sigma->Draw();
	cnll->cd(3);
	frame_staple2->Draw();
	frame_staple6->Draw("same");
	cnll->cd(4);
	frame_staple3->Draw();
	frame_staple7->Draw("same");
	cnll->cd(5);
	frame_staple4->Draw();
	frame_staple8->Draw("same");
	cnll->cd(6);
	frame_n1->Draw();
	frame_n2->Draw("same");

}
*/

