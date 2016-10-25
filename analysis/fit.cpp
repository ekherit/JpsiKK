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

#include <TSystem.h>

#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TMath.h>
#include <TTree.h>

#include <boost/format.hpp>

#include <RooSimultaneous.h>
#include <RooFitResult.h>
#include <RooChi2Var.h>
#include <RooNLLVar.h>
#include <RooAddition.h>
#include <RooRealVar.h>
#include <RooPolynomial.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooCategory.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooBukinPdf.h>
#include <RooMinuit.h>

#include <RooHist.h>
#include <RooCurve.h>
using namespace RooFit;

#include "RooMcbPdf.h"
#include "fit.h"

bool OPT_NOBGSLOPE=false; //no slope for the background
bool OPT_NOBG=false; //no background
bool OPT_NOGAUSRAD=false; //no gaus rad
std::string OPT_PARAM_CONFIG_FILE=""; 
bool OPT_SEPARATE_MREC=false; //separate Mrec for each channel
std::string OPT_FIT_METHOD="lh";
std::string OPT_FIT_RESULT_FILE_NAME="fit_result.txt";

void fit(TH1 * his)
{
  std::list<TH1*> hlst = {his};
  std::list<TTree*> tlst;
  fit(hlst,tlst);
}

void fit(TH1 * hisKK, TH1 * hisUU)
{
  hisKK->SetName("KK");
  hisUU->SetName("UU");
  std::list<TH1*> hlst = {hisKK,hisUU};
  std::list<TTree*> tlst;
  fit(hlst,tlst);
}


RooArgSet * createMcbVars(double Mmin, double Mmax)
{
  RooArgSet * args = new RooArgSet;
  args->add( * new RooRealVar("mean",  "mean",   -0.1 + 0.5*(Mmin+Mmax) ,  Mmin,  Mmax, "MeV"));
  args->add( * new RooRealVar("sigma", "sigma",1.4, 0, 10, "MeV") ) ;
  args->add( * new RooRealVar("n1","n1", 3,  1,100));
  args->add( * new RooRealVar("n2","n2", 2,  1,100));
  args->add( * new  RooRealVar("L1", "Left Gaus range"   , 2        , "MeV"));
  args->add( * new  RooRealVar("L2", "Left Exp range"    , 1        , 0      , (Mmax-Mmin)*0.5 , "MeV"));
  args->add( * new  RooRealVar("L3", "Left Power range"  , 10       , 0      , (Mmax-Mmin)*0.5 , "MeV"));
  args->add( * new  RooRealVar("L4", "Left Exp range 2"  , 20       , 0      , (Mmax-Mmin)*0.5 , "MeV"));
  args->add( * new  RooRealVar("R1", "Right Gaus range"  , 2        , "MeV") );
  args->add( * new  RooRealVar("R2", "Right Exp range"   , 1        , 0      , (Mmax-Mmin)*0.5 , "MeV")); 
  args->add( * new  RooRealVar("R3", "Right Power range" , 30       , 0      , (Mmax-Mmin)*0.5 , "MeV")); 
  args->add( * new  RooRealVar("R4", "Right Exp range 2" , 1        , 0      , (Mmax-Mmin)*0.5 , "MeV"));

  for(int i=0;i<4;i++)
  {
    std::string istr = to_string(i);
    double Mean, Min,Max;
    if(i<2)
    {
      Mean = 20;
      Min = (Mmax+Mmin)*0.5+5;
      Max = Mmax;
    }
    else
    {
      Mean = -20;
      Min = Mmin;
      Max = (Mmax+Mmin)*0.5-5;
    }
    Min = Mmin;
    Max = Mmax;
    auto mean = new RooRealVar(("rad_mean" + istr).c_str(), ("Radiative gauss mean " + istr).c_str(), Mean, Min, Max, "MeV"); 
    auto sigma = new RooRealVar(("rad_sigma" + istr).c_str(), ("Radiative gaus sigma" + istr).c_str(), 10, 2, 100, "MeV");
    auto frac = new RooRealVar(("rad_frac"+istr).c_str(),("fraction " + istr + "of radiative gauss").c_str(), 0, 0.2);
    if(OPT_NOGAUSRAD)
    {
      mean->setConstant();
      sigma->setConstant();
      frac->setConstant();
    }
    args->add(*mean);
    args->add(*sigma);
    args->add(*frac);
  }
  return args;
}

RooAbsPdf * createMcbPdf(std::string name, RooRealVar & Mrec, RooArgSet & pars, double Mmin, double Mmax)
{
  std::string mcb_name = "Mcb2Pdf_"+name;
  std::string mcb_title = "My modified Crystal Bal function for " + name;
  RooMcb2Pdf * mcbPdf =  new RooMcb2Pdf(mcb_name.c_str(), mcb_title.c_str(), 
      Mrec , 
      (RooRealVar&)pars["mean"],
      (RooRealVar&)pars["sigma"],
      (RooRealVar&)pars["L1"],
      (RooRealVar&)pars["L2"],
      (RooRealVar&)pars["L3"],
      (RooRealVar&)pars["L4"],
      (RooRealVar&)pars["R1"],
      (RooRealVar&)pars["R2"],
      (RooRealVar&)pars["R3"],
      (RooRealVar&)pars["R4"],
      (RooRealVar&)pars["n1"],
      (RooRealVar&)pars["n2"]
      );
  mcbPdf->setRange(Mmin,Mmax);
  return mcbPdf;
}


RooAbsPdf * addRad(RooAbsPdf * mcbPdf, RooRealVar & Mrec, RooArgSet & pars)
{
  std::vector<RooGaussian*> radPdf(2);
  double Mmin = Mrec.getMin();
  double Mmax = Mrec.getMax();
  std::string name = mcbPdf->GetName();
  //name create the SignalPdfs for each channel
  RooArgList PdfList;
  RooArgList RadFracList;
  for(int i=0;i<radPdf.size();i++)
  {
    std::string istr = to_string(i);
    radPdf[i] = new RooGaussian(
        ("gaus_radPdf" + name + istr).c_str(),
        ("Radiative gauss " + istr + " for " + name).c_str(), 
       Mrec, 
       (RooRealVar&)pars[("rad_mean"+istr).c_str()],
       (RooRealVar&)pars[("rad_sigma"+istr).c_str()]
       );
    PdfList.add(*radPdf[i]);
    RadFracList.add((RooRealVar&)pars[("rad_frac"+istr).c_str()]);
  }
  PdfList.add(*mcbPdf);
  return  new RooAddPdf(("signal"+name).c_str(),("Signal model for "+name).c_str(), PdfList, RadFracList);
}

RooAbsPdf *addBg(std::string name, RooAbsPdf * signalPdf, RooAbsPdf * bgPdf, RooRealVar & Mrec, RooArgSet & pars)
{
  double Mmin = Mrec.getMin();
  double Mmax = Mrec.getMax();
  RooRealVar * bg_c1 = new RooRealVar(("bg"+name+"_c1").c_str(),("background slope for "+name).c_str(),0,- 1.0/Mmax, - 1.0/Mmin);
  if(OPT_NOBGSLOPE) bg_c1->setConstant();
  pars.add(*bg_c1);
  bgPdf = new RooPolynomial(("bg"+name).c_str(), ( name + " background pdf").c_str(), Mrec, *bg_c1);	
  RooRealVar * Nsig = new RooRealVar(("Nsig"+name).c_str(), ("Number of " + name +" events").c_str(), 0,0, 1e8);
  RooRealVar * Nbg  = new RooRealVar(("Nbg"+name).c_str(), ("Number of background events for " + name + " channel").c_str(), 0, 0, 1e8);
  pars.add(*Nsig);
  pars.add(*Nbg);
  if(OPT_NOBG)
  {
    Nbg->setConstant();
    bg_c1->setConstant();
  }
  return new RooAddPdf
    (
     (name+"Pdf").c_str(), 
     (name+" signal + background p.d.f.").c_str(), 
     RooArgList(*bgPdf, *signalPdf), 
     RooArgList(*Nbg,  *Nsig)
    );
}




void draw_nll_scan(RooAbsReal * nll, RooArgSet & args, std::list<std::string> & name_lst)
{
  std::cout << "Scan nll variable over parametes" << std::endl;
  args.Print();
  RooRealVar & sigma = (RooRealVar&)args["sigma"];
  RooRealVar &  mean =  (RooRealVar&)args["mean"];
  std::vector<RooRealVar*> staple(8);
  for(int i=0;i<4;i++)
  {
    std::string istr = to_string(i+1);
    std::string L = "L"+istr;
    std::string R = "R"+istr;
    staple[i] =  (RooRealVar*)(& args[L.c_str()]);
    staple[i+4] =  (RooRealVar*)(& args[R.c_str()]);
  }
  RooRealVar &  n1 = (RooRealVar&) args["n1"];
  RooRealVar &  n2 =  (RooRealVar&)args["n2"];

  std::map<std::string, RooRealVar *> Nsig;
  std::map<std::string, RooRealVar *> Nbg;

  for(auto name : name_lst)
  {
    Nsig[name] = (RooRealVar*)(&args[("Nsig"+name).c_str()]);
    Nbg[name] = (RooRealVar*)(&args[("Nbg"+name).c_str()]);
  }


  RooPlot* frame_sigma = sigma.frame(Range(0.8, 5.0)) ;
  nll->plotOn(frame_sigma,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll->getVal()+10),LineColor(kRed)) ;

  RooPlot* frame_mean = mean.frame(Range(-2.0, 2.0)) ;
  nll->plotOn(frame_mean,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll->getVal()+10),LineColor(kRed)) ;

  std::vector<RooPlot*> frame_staple(staple.size());
  std::vector<int> colors = {kBlack, kBlue,kCyan+2,kGreen+2, kGreen, kRed, kMagenta, kRed+2}; 
  std::vector<int> line_styles = {kSolid, kSolid, kDashed, kDotted, kSolid, kSolid, kDashed, kDotted};
  for(int i=1;i<4;i++)
  {
    frame_staple[i] = staple[i]->frame();
    frame_staple[i+4] = staple[i+4]->frame();
    nll->plotOn(frame_staple[i],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll->getVal()+10),LineColor(colors[i]), LineStyle(line_styles[i])) ;
    nll->plotOn(frame_staple[i+4],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll->getVal()+10),LineColor(colors[i+4]),LineStyle(line_styles[i+4])) ;
  }


  RooPlot* frame_n1 = n1.frame(Range(1,10)) ;
  nll->plotOn(frame_n1,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll->getVal()+10),LineColor(kBlue)) ;

  RooPlot* frame_n2 = n2.frame(Range(1, 10)) ;
  nll->plotOn(frame_n2,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll->getVal()+10),LineColor(kRed));


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
  Ncanvas->Divide(2,name_lst.size());

  for(auto name : name_lst)
  {
    Nframe.push_back(Nsig[name]->frame(Range(Nsig[name]->getValV()*0.8,Nsig[name]->getValV()*1.2)));
    Nframe.push_back( Nbg[name]->frame(Range(Nbg[name]->getValV()*0.8, Nbg[name]->getValV()*1.2)));
  }
  for(int i=0;i<Nframe.size();i++)
  {
    nll->plotOn(Nframe[i],PrintEvalErrors(-1),ShiftToZero(), EvalErrorValue(nll->getVal()+10),LineColor(kRed)) ;
    Ncanvas->cd(i+1);
    Nframe[i]->Draw();
  }
  std::cout << "Done" << std::endl;
}

void fit(std::list<TH1*> & hlst, std::list<TTree*> & tree_list, bool use_tree)
{
  //define maximum and minimum recoil mass
  double Mmin  = (*std::max_element(std::begin(hlst), std::end(hlst), 
      [](const TH1 * h1, const TH1 *h2) { return h1->GetXaxis()->GetXmin() < h2->GetXaxis()->GetXmin(); } ))->GetXaxis()->GetXmin();

  double Mmax  = (*std::min_element(std::begin(hlst), std::end(hlst), 
      [](const TH1 * h1, const TH1 *h2) { return h1->GetXaxis()->GetXmax() < h2->GetXaxis()->GetXmax(); } ))->GetXaxis()->GetXmax();


  std::list < std::string > name_lst;
  std::map<std::string,TH1*> hisMap;
  std::map<std::string, Long64_t> Nhis;
  std::map<std::string, TTree*> treeMap;
  for(auto h : hlst) 
  {
    name_lst.push_back(h->GetName());
    hisMap[h->GetName()] = h;
    Nhis[h->GetName()] = h->GetEntries();
  }
  if(use_tree)
  {
    name_lst.clear();
    for(auto t : tree_list)
    {
      if(t==nullptr) continue;
      name_lst.push_back(t->GetName());
      treeMap[t->GetName()] = t;
      Nhis[t->GetName()] = t->GetEntries();
    }
  }
  //i decided to make own variable for mmrec
  std::map<std::string, RooAbsReal*> MMrec;

  boost::format mrec_title_fmt("M_{rec}^{%s}(#pi^{+}#pi^{-}) - 3097 MeV");
  RooRealVar  MrecInTree("Mrec","Mrec in tree",3.097-0.045, 3.097+0.045, "GeV");
  RooRealVar  Mrec("M",(mrec_title_fmt % "").str().c_str(), Mmin,Mmax, "MeV");
  RooFormulaVar DM ("DM","(Mrec-3.097)*1000",MrecInTree);
  std::cout << "Before Mrec asignment" << std::endl;
  for(auto name : name_lst)
  {
    std::cout << "MrecAsign 1" << std::endl;
    MMrec[name]=&Mrec; //default observable
    std::cout << "MrecAsign 2" << std::endl;
    if(use_tree) MMrec[name]=&DM;
    std::cout << "MrecAsign 3" << std::endl;
    if(OPT_SEPARATE_MREC)
    {
      if(!use_tree)
      {
        MMrec[name] = new RooRealVar(
            ("M"+name).c_str(), 
            (mrec_title_fmt % name).str().c_str(), 
            hisMap[name]->GetXaxis()->GetXmin(),
            hisMap[name]->GetXaxis()->GetXmax(), 
            "MeV"
            );
      }
      else
      {
        MMrec[name] = new RooFormulaVar(("DM"+name).c_str(),"(Mrec-3.097)*1000", MrecInTree);
      }
    }
    //MMrec[name]->setBins(hisMap[name]->GetNbinsX());
    //Mrec.setBins(hisMap[name]->GetNbinsX(),name.c_str());
    //Mrec.setBins(10000);
  }
  std::cout << "After MrecAsign" << std::endl;


  RooArgSet & args = * createMcbVars(Mmin,Mmax);

  /*
  RooRealVar  mean( "mean",  "mean",   -0.1 + 0.5*(Mmin+Mmax) ,  Mmin,  Mmax, "MeV") ;
  RooRealVar sigma("sigma","sigma",1.4,0.5,10, "MeV") ;
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
  */

  /*  
  RooRealVar & mean  = (RooRealVar&)args["mean"];
  RooRealVar & sigma = (RooRealVar&)args["sigma"];
  RooRealVar & n1 =    (RooRealVar&)args["n1"];
  RooRealVar & n2 =    (RooRealVar&)args["n2"];
  std::vector<RooRealVar*> staple =
  {
    (RooRealVar*)&args["L1"],
    (RooRealVar*)&args["L2"],
    (RooRealVar*)&args["L3"],
    (RooRealVar*)&args["L4"],
    (RooRealVar*)&args["R1"],
    (RooRealVar*)&args["R2"],
    (RooRealVar*)&args["R3"],
    (RooRealVar*)&args["R4"],
  };
  */

  std::map<std::string, RooAbsPdf*> McbPdfMap;
  std::map<std::string, RooAbsPdf*> SignalPdfMap;
  /*  
  std::cout << "Before mcb pdf construction" << std::endl;
  for(auto name : name_lst)
  {
    std::cout << "MMrec[name]->GetName() = " << MMrec[name]->GetName() << std::endl;
    //std::cout << "min = " << ((RooRealVar*)MMrec[name])->getMin() << " max = " << ((RooRealVar*)MMrec[name])->getMax() << std::endl;
    std::cout << "In mcb pdf construction " << name <<  std::endl;
    //McbPdfMap[name]=new RooMcb2Pdf(("Mcb2Pdf_"+name).c_str(), ("My modified Crystal Bal function for " + name).c_str(),
    //    *MMrec[name], mean,sigma,staple, n1, n2);
    McbPdfMap[name]=new RooMcb2Pdf(("Mcb2Pdf_"+name).c_str(), ("My modified Crystal Bal function for " + name).c_str(),
        *MMrec[name], mean,sigma,
        *staple[0], 
        *staple[1], 
        *staple[2], 
        *staple[3], 
        *staple[4], 
        *staple[5], 
        *staple[6], 
        *staple[7], 
        n1, n2);
    ((RooMcb2Pdf*)McbPdfMap[name])->setRange(Mmin,Mmax);
  }
  std::cout << "After mcb construction" << std::endl;
  */
	//RooMcb2Pdf *mcbPdf=0;
  //mcbPdf = new RooMcb2Pdf("mcb2","mcb2",Mrec,mean,sigma, staple, n1,n2,Mmin,Mmax);


  for(auto name: name_lst)
  {
    RooRealVar & Mobs = *(RooRealVar*)MMrec[name];
    McbPdfMap[name] = createMcbPdf(name, Mobs, args,Mmin,Mmax);
    SignalPdfMap[name] = OPT_NOGAUSRAD ? McbPdfMap[name] : addRad(McbPdfMap[name], Mobs, args);
  }


  //std::vector<RooRealVar*> meanRad;
  //if(!OPT_NOGAUSRAD) meanRad.resize(4);
  //std::vector<RooRealVar*> sigmaRad(meanRad.size());
  //std::vector<RooGaussian*> radPdf(meanRad.size());
  //std::vector<RooRealVar*> radFrac(meanRad.size());

  //std::map<std::string, std::vector<RooGaussian*>> radPdfMap;

  //create parameters for radiative gaussian (mean and sigma)
  /*
  for(int i=0;i<meanRad.size();i++)
  {
    std::string istr = to_string(i);
    double Mean, Min,Max;
    if(i<2)
    {
      Mean = 20;
      Min = (Mmax+Mmin)*0.5+5;
      Max = Mmax;
    }
    else
    {
      Mean = -20;
      Min = Mmin;
      Max = (Mmax+Mmin)*0.5-5;
    }
    Min = Mmin;
    Max = Mmax;
    meanRad[i]  = new RooRealVar(("rad_mean" + istr).c_str(),   ("mean_rad" + istr).c_str(), Mean, Min, Max, "MeV");
    sigmaRad[i] = new RooRealVar(("rad_sigma" + istr).c_str(), ("sigma_rad" + istr).c_str(), 10, 2, 100, "MeV");
    radFrac[i] = new RooRealVar(("rad_frac"+istr).c_str(),("fraction " + istr + " radiative gauss").c_str(), 0, 0.2);
  }
  */

  //for(int i=0;i<meanRad.size();i++)
  //{
  //  std::string istr = to_string(i);
  //  meanRad[i] = (RooRealVar*)(&args[("rad_mean"+istr).c_str()]);
  //  sigmaRad[i] = (RooRealVar*)(&args[("rad_sigma"+istr).c_str()]);
  //  radFrac[i] = (RooRealVar*)(&args[("rad_frac"+istr).c_str()]);
  //}

  ////name create the SignalPdfs for each channel
  //for(auto name : name_lst)
  //{
  //  RooArgList PdfList;
  //  RooArgList RadFracList;
  //  radPdfMap[name].resize(meanRad.size());

  //  for(int i=0;i<meanRad.size();i++)
  //  {
  //    std::string istr = to_string(i);
  //    auto pdf = new RooGaussian(("radPdf" + name + istr).c_str(),
  //        ("Radiative gauss " + istr + " for " + name).c_str(), *MMrec[name], *meanRad[i], *sigmaRad[i]);
  //    radPdfMap[name][i] = pdf;
  //    PdfList.add(*pdf);
  //    RadFracList.add(*radFrac[i]);
  //  }
  //  //PdfList.add(*mcbPdf);
  //  PdfList.add(*McbPdfMap[name]);
  //  //RooAbsPdf * signalPdf =  new RooAddPdf("signalPdf","Signal model", PdfList, RadFracList);
  //  SignalPdfMap[name] = new RooAddPdf(("signalPdf"+name).c_str(),("Signal model for "+name).c_str(), PdfList, RadFracList);
  //}



  //Mixing with the background
  std::map<std::string, RooAbsPdf * > SamplePdf; //here will be Pdf with signal and background mixed

  std::map<std::string, RooRealVar * > Nsig; //number of signal events
  std::map<std::string, RooRealVar * > Nbg; //number of background events
  std::map<std::string, RooPolynomial*> bgPdf; //background pdf function
  std::map<std::string, RooRealVar*> bg_c1; //slope of the background
  std::map<std::string, RooPlot*> frame;

 	RooCategory sample("sample","sample");
  RooArgList MrecArgList;
  std::map<std::string, RooDataHist *> dataMap;
  std::map<std::string, RooDataSet *> dataSetMap;
  RooPlot * Frame=nullptr;
  if(!use_tree) Frame= Mrec.frame(Title("#pi^{+}#pi^{-} recoil mass" ));
  else Frame= MrecInTree.frame(Title("#pi^{+}#pi^{-} recoil mass" ));

  std::cout << "Before data initialization: " << std::endl;
  for ( auto name : name_lst)
  {
    bg_c1[name] = new RooRealVar(("bg"+name+"_c1").c_str(),("background slope for "+name).c_str(),0,- 1.0/Mmax, - 1.0/Mmin);
    if(OPT_NOBGSLOPE) bg_c1[name]->setConstant();
    bgPdf[name] = new RooPolynomial(("bg"+name).c_str(), ( name + " background pdf").c_str(), *MMrec[name], *bg_c1[name]);	
    Nsig[name] = new RooRealVar(("Nsig"+name).c_str(), ("Number of " + name +" events").c_str(), Nhis[name],0, Nhis[name]*100);
    Nbg[name] = new RooRealVar(("Nbg"+name).c_str(), ("Number of background events for " + name + " channel").c_str(), 0, 0, Nhis[name]*100);
    if(OPT_NOBG)
    {
      Nbg[name]->setConstant();
      bg_c1[name]->setConstant();
    }
    SamplePdf[name] = new RooAddPdf((name+"Pdf").c_str(), (name+" signal + background p.d.f.").c_str(), 
        //RooArgList(*bgPdf[name], *signalPdf), 
        RooArgList(*bgPdf[name], *SignalPdfMap[name]), 
        RooArgList(*Nbg[name],  *Nsig[name]));

    if(!use_tree)
    {
      dataMap[name] = new RooDataHist(name.c_str(), name.c_str(), *MMrec[name], Import(*hisMap[name]));
    }
    else
    {
      dataMap[name] = new RooDataHist(name.c_str(), name.c_str(), Mrec, Import(*hisMap[name]));
      std::cout << "Initializing RooDataSet: name = " << name << " Mrec name = " << MMrec[name]->GetName() << std::endl;
      dataSetMap[name] = new RooDataSet(name.c_str(), name.c_str(), MrecInTree, Import(*treeMap[name]));
      dataSetMap[name]->addColumn(*MMrec[name]);
    }
    sample.defineType(name.c_str());
    MrecArgList.add(*MMrec[name]);
    if(OPT_SEPARATE_MREC) frame[name] = MrecInTree.frame(Title(("#pi^{+}#pi^{-} recoil mass for " + name).c_str()));
    else frame[name] = Frame;
  }
  std::cout << "After data initialization" << std::endl;

  RooAbsData * data = nullptr;
  if(!use_tree) data = new RooDataHist("combData","Combined data", MrecArgList, sample, dataMap);
  else data = new RooDataSet("comb","Combined data", MrecInTree, Index(sample), Import(dataSetMap));
  std::cout << "After combined data set " << std::endl;

  RooSimultaneous simPdf("simPdf","simultaneous pdf",SamplePdf,sample) ;

  std::cout << "Before chi2 nll initialization" << std::endl;
  std::map<std::string, RooNLLVar*> nllMap;
  std::map<std::string, RooChi2Var*> chi2Map;
  RooArgSet  nll_arg_set;
  RooArgSet  chi2_arg_set;
  for(auto name : name_lst)
  {
    nllMap[name] = new RooNLLVar(("nll"+name).c_str(),("Nll for " + name).c_str(),*SamplePdf[name],*dataMap[name]);
    nll_arg_set.add(*nllMap[name]);
    chi2Map[name] = new RooChi2Var(("chi2"+name).c_str(),("chi2 for " + name).c_str(),*SamplePdf[name],*dataMap[name],
        DataError(RooAbsData::Expected));
    chi2_arg_set.add(*chi2Map[name]);
  }

  RooAddition nll("nll","nll", nll_arg_set);
	RooAddition chi2Var("chi2", "chi2", chi2_arg_set);
	//RooNLLVar nll("nll","nll",simPdf,*data) ;

  //RooDataHist * data = new RooDataHist("combData","Combined data", MrecArgList, sample, hisMap);
  std::cout << "After chi2 nll initialization" << std::endl;

  auto p = simPdf.getParameters(*MMrec[name_lst.front()]);
  if(OPT_PARAM_CONFIG_FILE!="")
  {
    p->readFromFile(OPT_PARAM_CONFIG_FILE.c_str());
    p->Print("v");
  }

  
  if(false)
  {
    for(auto name : name_lst)
    {
      new TCanvas;
      RooPlot * f;
      if(use_tree) f = MrecInTree.frame(Title(("#pi^{+}#pi^{-} recoil mass for " + name).c_str()));
      else f = ((RooRealVar*)MMrec[name])->frame(Title(("#pi^{+}#pi^{-} recoil mass for " + name).c_str()));
      data->plotOn(f, XErrorSize(0), MarkerSize(0.5),  Cut(("sample==sample::"+name).c_str()), AutoBinning());
      simPdf.plotOn(f, Slice(sample,name.c_str()),ProjWData(sample,*data), LineWidth(2),LineColor(kRed));
      f->SetMinimum(0.1);
      f->GetYaxis()->SetTitleOffset(1.6) ; 
      f->Draw();
    }
    while(true)
    {
      gSystem->ProcessEvents();
      usleep(100000); 
    }
  }

  RooFitResult * theFitResult=nullptr;
  if(OPT_FIT_METHOD == "lh")
  {
    std::cout << "Starting standart likelihood fit" << std::endl;
    //standart fit method it is used to mach calculation for simpultaneuse fit on separate Mrec
    theFitResult = simPdf.fitTo(*data, Extended(), Strategy(2), Minos(), Save());
    std::cout << "Fit done" << std::endl;
  }
  if(OPT_FIT_METHOD == "chi2")
  {
    std::cout << "Starting chi2 fit " << std::endl;
    RooMinuit minuit(chi2Var);
    minuit.migrad();
    minuit.hesse();
    minuit.minos();
    theFitResult = minuit.save();
    std::cout << "Fit done" << std::endl;
  }
  if(OPT_FIT_METHOD == "lh2")
  {
    std::cout << "Starting my likelihood fit " << std::endl;
    RooMinuit minuit(nll);
    minuit.migrad();
    minuit.hesse();
    minuit.minos();
    theFitResult = minuit.save();
    std::cout << "Fit done" << std::endl;
  }

  std::cout << "Writing result of the fit into file: " << OPT_FIT_RESULT_FILE_NAME << std::endl;
  p->writeToFile(OPT_FIT_RESULT_FILE_NAME.c_str());
  std::cout << "Writing done" << std::endl;

  std::cout << "Calculating chi2 and number of degree of freedom" << std::endl;
  //number of free parameters
  int nfp = theFitResult->floatParsFinal().getSize() ;
  //chi square
	double chi2 =  chi2Var.getVal();
  int Ndata_entries=0;
  for(auto his : hisMap) Ndata_entries+=his.second->GetNbinsX();
  //number of degree of freedom
  int ndf = Ndata_entries - nfp;
  //normalized chi2
  double chi2ndf = chi2/ndf;
  double chi2prob = TMath::Prob(chi2,ndf);

  draw_nll_scan(&nll, *simPdf.getParameters(data), name_lst);




  std::cout << "Plotting data and pdf" << std::endl;
  //RooPlot * Frame = Mrec.frame(Title("#pi^{+}#pi^{-} recoil mass"));
  std::vector<int> colors = {kBlack, kBlue,kCyan+2,kGreen+2, kGreen, kRed, kMagenta, kRed+2}; 
  std::vector<int> line_styles = {kSolid, kSolid, kDashed, kDotted, kSolid, kSolid, kDashed, kDotted};
  RooArgSet viewArgSet;
  int i=0;
  for(auto name : name_lst)
  {
    std::cout << "     plotting " << name  << std::endl;
    data->plotOn(frame[name], XErrorSize(0), MarkerSize(0.5),  Cut(("sample==sample::"+name).c_str()),LineColor(colors[i]), MarkerColor(colors[i]), AutoBinning());
    std::cout << "     plotting Pdf for" << name  << std::endl;
    simPdf.plotOn(frame[name], Slice(sample,name.c_str()),ProjWData(sample,*data), LineWidth(1),LineColor(colors[i]));
    simPdf.plotOn(frame[name], Slice(sample,name.c_str()),ProjWData(sample,*data), Components(*bgPdf[name]),LineStyle(kDashed),LineWidth(1),LineColor(colors[i]));
    frame[name]->SetMinimum(0.1);
    frame[name]->GetYaxis()->SetTitleOffset(1.6) ; 
    viewArgSet.add(*Nsig[name]);
    viewArgSet.add(*Nbg[name]);
    std::cout << "chi^2 from frame  = " << frame[name]->chiSquare() << endl ;
    i++;
  }
  std::cout << "Done" << std::endl;
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
  if(OPT_SEPARATE_MREC)
  {
    for(auto name :name_lst)
    {
      if(name != most_entries_his_name)
      {
        frame[name]->Draw("same");
      }
    }
  }


  //Frame->GetYaxis()->SetTitleOffset(1.6) ; 
  //
  //Frame->Draw() ;


  args.Print();
  /*  
  mean.Print();
  sigma.Print();
  for(auto name : name_lst)
  {
    Nsig[name]->Print();
    Nbg[name]->Print();
  }
  */
  std::cout << boost::format("chi2/ndf = %.15f/(%d-%d) = %f,  prob = %f") % chi2 % Ndata_entries % nfp % chi2ndf %  chi2prob<< std::endl;


  //my own chi2 calculation
  /*
  double my_sum_chi2=0;
  for(auto name : name_lst)
  {
    //RooHist* histo = (RooHist*) frame[name]->findObject(name.c_str()) ;
    //RooCurve* func = (RooCurve*) frame[name]->findObject((name+"Pdf").c_str()) ;
    RooHist* histo = (RooHist*) frame[name]->getHist() ;
    RooCurve* func = (RooCurve*) frame[name]->getCurve() ;
    for (Int_t i=0 ; i<histo->GetN() ; i++) 
    {
      Double_t xdata,ydata ;
      histo->GetPoint(i,xdata,ydata) ;
      Double_t yfunc = func->interpolate(xdata) ;
      std::cout << "n=" << ydata << " n_obs  - n_exp" << (ydata-yfunc) << "  (n_obs-n_exp)^2/n_obs =" << (ydata-yfunc)*(ydata-yfunc)/ydata << std::endl;
      my_sum_chi2 += (ydata-yfunc)*(ydata-yfunc)/ydata;
    }
  }
  std::cout << "mysum_chi2 = " << my_sum_chi2 <<std::endl;
  */

  //No print the ration of signal event to number of signal event for last data
  //sample
  for(auto name : name_lst)
  {
    std::string name0 = name_lst.back();
    if(name == name0) continue;
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

void fit2(std::list<TH1*> & hlst)
{
  //define maximum and minimum recoil mass
  double Mmin  = (*std::max_element(std::begin(hlst), std::end(hlst), 
      [](const TH1 * h1, const TH1 *h2) { return h1->GetXaxis()->GetXmin() < h2->GetXaxis()->GetXmin(); } ))->GetXaxis()->GetXmin();

  double Mmax  = (*std::min_element(std::begin(hlst), std::end(hlst), 
      [](const TH1 * h1, const TH1 *h2) { return h1->GetXaxis()->GetXmax() < h2->GetXaxis()->GetXmax(); } ))->GetXaxis()->GetXmax();

  std::list < std::string > name_lst;
  std::map<std::string,TH1*> hisMap;
  std::map<std::string, Long64_t> Nhis;
  for(auto h : hlst) 
  {
    name_lst.push_back(h->GetName());
    hisMap[h->GetName()] = h;
    Nhis[h->GetName()] = h->GetEntries();
  }
  //i decided to make own variable for mmrec
  std::map<std::string, RooRealVar*> MMrec;

  boost::format mrec_title_fmt("M_{rec}^{%s}(#pi^{+}#pi^{-}) - 3097 MeV");
  RooRealVar  Mrec("Mrec",(mrec_title_fmt % "").str().c_str(), Mmin,Mmax, "MeV");

  std::map<std::string, RooAbsPdf * > SamplePdf; //here will be Pdf with signal and background mixed

  std::map<std::string, RooRealVar * > Nsig; //number of signal events
  std::map<std::string, RooRealVar * > Nbg; //number of background events
  std::map<std::string, RooAbsPdf * > bgPdf; //number of background events
  std::map<std::string, RooPlot*> frame;

  RooArgSet & pars = *createMcbVars(Mmin,Mmax);
 	RooCategory sample("sample","sample");
  for(auto name : name_lst)
  {
    MMrec[name]=&Mrec; //default observable
    //MMrec[name] = new RooRealVar(
    //    ("Mrec"+name).c_str(), 
    //     (mrec_title_fmt % name).str().c_str(), 
    //     hisMap[name]->GetXaxis()->GetXmin(),
    //    hisMap[name]->GetXaxis()->GetXmax(), "MeV"
    //    );
    MMrec[name]->setBins(hisMap[name]->GetNbinsX());
    //Mrec.setBins(hisMap[name]->GetNbinsX(),name.c_str());
    //Mrec.setBins(10000);
    auto mcb = createMcbPdf(name, *MMrec[name], pars, Mmin,Mmax);
    auto mcb_rad = addRad(mcb, *MMrec[name], pars);
    RooAbsPdf * bgpdf;
    SamplePdf[name] = addBg(name, mcb_rad ,bgpdf ,*MMrec[name],pars);
    bgPdf[name] = bgpdf;
    sample.defineType(name.c_str());

    Nsig[name]= &((RooRealVar&)pars[("Nsig"+name).c_str()]);
    Nbg[name]= &((RooRealVar&)pars[("Nbg"+name).c_str()]);
    Nsig[name]->setVal(hisMap[name]->GetEntries());
  }



  //data
  std::map<std::string, RooDataHist *> dataMap;
  RooArgList MrecArgList;
  for ( auto name : name_lst)
  {
    dataMap[name]=new RooDataHist(name.c_str(), name.c_str(), *MMrec[name], Import(*hisMap[name]));
    MrecArgList.add(*MMrec[name]);
    frame[name] = MMrec[name]->frame(Title(("#pi^{+}#pi^{-} recoil mass for " + name).c_str()));
  }

  RooDataHist * data = new RooDataHist("combData","Combined data", MrecArgList, sample, dataMap);
  RooSimultaneous simPdf("simPdf","simultaneous pdf",SamplePdf,sample) ;

  std::map<std::string, RooNLLVar*> nllMap;
  std::map<std::string, RooChi2Var*> chi2Map;
  RooArgSet  nll_arg_set;
  RooArgSet  chi2_arg_set;
  for(auto name : name_lst)
  {
    nllMap[name] = new RooNLLVar(("nll"+name).c_str(),("Nll for " + name).c_str(),*SamplePdf[name],*dataMap[name]);
    nll_arg_set.add(*nllMap[name]);
    chi2Map[name] = new RooChi2Var(("chi2"+name).c_str(),("chi2 for " + name).c_str(),*SamplePdf[name],*dataMap[name], DataError(RooAbsData::Expected));
    chi2_arg_set.add(*chi2Map[name]);
  }

  RooAddition nll("nll","nll", nll_arg_set);
	//RooAddition chi2Var("chi2", "chi2", chi2_arg_set);
	RooChi2Var chi2Var("chi2", "chi2", simPdf, *data, DataError(RooAbsData::Expected));
	//RooNLLVar nll("nll","nll",simPdf,*data) ;

  //RooDataHist * data = new RooDataHist("combData","Combined data", MrecArgList, sample, hisMap);

  auto p = simPdf.getParameters(*MMrec[name_lst.front()]);
  if(OPT_PARAM_CONFIG_FILE!="")
  {
    p->readFromFile(OPT_PARAM_CONFIG_FILE.c_str());
    p->Print("v");
  }


  //RooMinuit minuit(chi2Var);
  //minuit.migrad();
  //minuit.hesse();
  //minuit.minos();
  //auto theFitResult = minuit.save();

	auto theFitResult = simPdf.fitTo(*data, Extended(), Strategy(2), Minos(), Save());
  theFitResult->Print();

  p->writeToFile("tmp_fit_result.txt");

  //number of free parameters
  int nfp = theFitResult->floatParsFinal().getSize() ;
  //chi square
	double chi2 =  chi2Var.getVal();
  //number of degree of freedom
  int Ndata_entries=0;
  for(auto his : hisMap) Ndata_entries+=his.second->GetNbinsX();
  //number of degree of freedom
  int ndf = Ndata_entries - nfp;
  //normalized chi2
  double chi2ndf = (Ndata_entries-nfp)/ndf;
  //probability for this hypotisis to get chi2 higher then for this data
  double chi2prob = TMath::Prob(chi2,ndf);

  int i=0;
  std::vector<int> colors = {kBlack, kBlue,kCyan+2,kGreen+2, kGreen, kRed, kMagenta, kRed+2}; 
  std::vector<int> line_styles = {kSolid, kSolid, kDashed, kDotted, kSolid, kSolid, kDashed, kDotted};
  for(auto name : name_lst)
  {
    std::cout << "plotting data" << name << std::endl;
    data->plotOn(frame[name], XErrorSize(0), MarkerSize(0.5),  Cut(("sample==sample::"+name).c_str()),LineColor(colors[i]), MarkerColor(colors[i]));
    std::cout << "plotting simPdf " << name << std::endl;
    simPdf.plotOn(frame[name], Slice(sample,name.c_str()),ProjWData(sample,*data), LineWidth(1),LineColor(colors[i]));
    std::cout << "plotting simPdf with bgcomptonents " << name << std::endl;
    /*  
    simPdf.plotOn(frame[name], 
        Slice(sample,name.c_str()),
        ProjWData(sample,*data), 
        Components(* bgPdf[name]),
        LineStyle(kDashed),
        LineWidth(1),
        LineColor(colors[i]));
        */
    std::cout << "frame sete minimum " << name << std::endl;
    frame[name]->SetMinimum(0.1);
    std::cout << "frame set title " << name << std::endl;
    frame[name]->GetYaxis()->SetTitleOffset(1.6) ; 
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
  //simPdf.paramOn(frame[most_entries_his_name], Label(chi_fmt.str().c_str()), Parameters(viewArgSet));

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
}


