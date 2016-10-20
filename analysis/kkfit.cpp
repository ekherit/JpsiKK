// =====================================================================================
//
//       Filename:  test.cpp
//
//    Description:  
//
//        Version:  1.0
//        Created:  03.07.2015 15:26:30
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================
#include <iostream>
#include <limits>

#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>
extern void InitGui();
VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
TROOT root("kkfit","kkfit", initfuncs);


#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TTree.h>
#include <TCut.h>
#include <TLegend.h>
#include <TMath.h>

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
#include <TCanvas.h>
#include <RooPlot.h>
#include <TAxis.h>
#include <RooBukinPdf.h>
#include <RooSimultaneous.h>

#include "RooMcbPdf.h"

#include "hashlist.h"

#include "fit.h"

bool BG_NOSLOPE=false;
bool BG_FIX_TO_ZERO=false;
bool GAUS_RAD=false;

using namespace RooFit;
using namespace std;

struct RooFitItem_t
{
  TH1           * his;
  RooDataHist   * data;
  RooRealVar    * Nsig;
  RooRealVar    * Nbg;
  RooRealVar    * bg_c1;
  RooRealVar    * bgB;
  RooPolynomial * bgPdf;
  RooAddPdf     * addPdf;
  RooPlot       * frame;
  std::string fName;
  double Mmin;
  double Mmax;

  RooFitItem_t(TH1 * h, RooAbsPdf* mcbPdf, RooRealVar & Mrec, double mmin, double mmax) 
    : RooFitItem_t(h->GetName(), h, mcbPdf, Mrec, mmin, mmax)
  {
  }

  RooFitItem_t(std::string name, TH1 * h, RooAbsPdf* mcbPdf, RooRealVar & Mrec, double mmin, double mmax)
  {
    fName = name;
    Mmin = mmin;
    Mmax = mmax;
    his = h;
    data = new RooDataHist(h->GetName(),  h->GetTitle(), Mrec, Import(*h));
    Nsig = new RooRealVar(
        (std::string("Nsig")+h->GetName()).c_str(), 
        (std::string("Number of signal events for  ") + h->GetTitle()).c_str(), 
        h->GetEntries(), 
        0, 
        h->GetEntries()*100
        );
    if(!BG_FIX_TO_ZERO)
    {
      Nbg = new RooRealVar(
          (std::string("Nbg")+h->GetName()).c_str(), 
          (std::string("Number of substrate events for  ") + h->GetTitle()).c_str(), 
          0, 
          0, 
          h->GetEntries()*100
          );
    }
    else
    {
      Nbg = new RooRealVar(
          (std::string("Nbg")+h->GetName()).c_str(), 
          (std::string("Number of substrate events for  ") + h->GetTitle()).c_str(), 
          0);
    }
    std::string bg_name = std::string("bg")+h->GetName();
    double bgc1min=-10; 
    double bgc1max=+10;
    if(Mmin < 0 && Mmax > 0) 
    {
      bgc1min = - 1.0/Mmax;  
      bgc1max = - 1.0/Mmin;
    }
    if(Mmin > 0 && Mmax > 0)
    {
      bgc1min = - 1.0/Mmax;
      bgc1max = std::numeric_limits<double>::max();
    }
    if(Mmin > 0 && Mmax < 0)
    {
      std::cerr << "ERROR: wrong fit range: Mmin (" << Mmin << ") > Mmax ("<<Mmax<<"): "  << std::endl;
      exit(1);
    }
    if(Mmin < 0 && Mmax < 0)
    {
      bgc1min = -std::numeric_limits<double>::max();
      bgc1max = - 1.0/Mmin;
    }
    
    bg_c1 = new RooRealVar((bg_name+"c1").c_str(),(std::string("Linear coefficient for background for ") + h->GetName()).c_str(),0, bgc1min, bgc1max);
    //bgB = new RooRealVar((bg_name+"B").c_str(),(std::string("Background shift") + h->GetName()).c_str(),0, -10000, 10000);
    std::string bg_title = std::string("Background Pdf for ") + h->GetName();
    //when proceed data with muond (U channel) not using slope of the background
    if(BG_NOSLOPE)
    {
      bgPdf = new RooPolynomial(bg_name.c_str(),bg_title.c_str(), Mrec);
    }
    else
    {
      if(fName.find("U") == std::string::npos)
      {
        bgPdf = new RooPolynomial(bg_name.c_str(),bg_title.c_str(), Mrec, *bg_c1);
      }
      else
      {
        bgPdf = new RooPolynomial(bg_name.c_str(),bg_title.c_str(), Mrec, *bg_c1);
        //bgPdf = new RooPolynomial(bg_name.c_str(),bg_title.c_str(), Mrec);
      }
    }
    //bgPdf = new RooBgPdf(bg_name.c_str(), bg_title.c_str(), Mrec, *bgB, his->GetXaxis()->GetXmin(), his->GetXaxis()->GetXmax());
    addPdf = new RooAddPdf(name.c_str(),h->GetTitle(), RooArgList(*bgPdf,*mcbPdf), RooArgList(*Nbg, *Nsig));
  }

  void plotData(RooDataHist * d)
  {
    d->plotOn(frame, MarkerSize(0.5),  Cut((std::string("sample==sample::")+his->GetName()).c_str())) ;
  }

  void plotPdf(RooDataHist *d, RooSimultaneous & simPdf,  RooCategory & sample)
  {
    simPdf.plotOn(frame, Slice(sample,his->GetName()),ProjWData(sample,*d), LineWidth(1)) ;
    simPdf.plotOn(frame, Slice(sample,his->GetName()),ProjWData(sample,*d), Components(*bgPdf),LineStyle(kDashed),LineWidth(1)) ;
  }
  std::string name(void) const 
  {
    return fName;
  }
};

/*  

void combfit(TH1 * h1, TH1 *h2)
{
  double Mmin  = h1->GetXaxis()->GetXmin();
  double Mmax  = h1->GetXaxis()->GetXmax();
  assert(Mmin == h2->GetXaxis()->GetXmin() && Mmax == h2->GetXaxis()->GetXmax());
  RooRealVar Mrec("Mrec","M_{rec}(#pi^{+}#pi^{-})", Mmin,Mmax, "MeV");
  RooRealVar sigma("sigma","sigma",1.4,0,  10, "MeV") ;
  RooRealVar staple1("staple1","staple1",  15,   Mmin,Mmax, "MeV") ;
  RooRealVar staple2("staple2","staple2",  3,   Mmin,Mmax, "MeV") ;
  RooRealVar staple3("staple3","staple3",  2,   Mmin,Mmax, "MeV") ;
  RooRealVar staple4("staple4","staple4",  0,   2, 2, "MeV") ;
  RooRealVar staple5("staple5","staple5",  2,   Mmin,Mmax, "MeV") ;
  RooRealVar staple6("staple6","staple6",  4,   Mmin,Mmax, "MeV") ;
  RooRealVar staple7("staple7","staple7",  15,   Mmin,Mmax, "MeV") ;
  RooRealVar n1("n1","n1", 3,  1,100) ;
  RooRealVar n2("n2","n2", 2,  1,100) ;

  RooMcbPdf * mcbPdf =  new RooMcbPdf("ModCB", "Simple Modified CrystalBall: gaus + power + exp",  
      Mrec,  
      sigma,  
      staple1, 
      staple2, 
      staple4, 
      staple6, 
      staple7, 
      n1, 
      n2);

 	RooCategory sample("sample","sample") ;
  sample.defineType(h1->GetName()) ;
  sample.defineType(h2->GetName()) ;


  std::clog << "Before FitItem_t f1" << std::endl;
  RooFitItem_t f1(h1->GetName(),h1, mcbPdf, Mrec, Mmin,Mmax);
  std::clog << "Before FitItem_t f2" << std::endl;
  RooFitItem_t f2(h2->GetName(),h2, mcbPdf, Mrec, Mmin,Mmax);

  std::clog << "Before RooDataHist" << std::endl;
  RooDataHist * data = new RooDataHist("data","combined h1 + h2",Mrec,Index(sample),Import(h1->GetName(),*h1),Import(h2->GetName(),*h2)) ;

  std::clog << "Before simPdf" << std::endl;
  RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;

  // Associate model with the physics state and model_ctl with the control state
  std::clog << "Before simPdf.addPdf(f1)" << std::endl;
  simPdf.addPdf(*(f1.addPdf),f1.his->GetName()) ;
  std::clog << "Before simPdf.addPdf(f2)" << std::endl;
  simPdf.addPdf(*(f2.addPdf),f2.his->GetName()) ;
  std::clog << "Before fitTo"<< std::endl;
	simPdf.fitTo(*data, Extended(), Strategy(2), Minos());

  f1.plotData(data);
  f2.plotData(data);
  f1.plotPdf(data, simPdf,sample);
  f2.plotPdf(data, simPdf,sample);


  TCanvas * c = new TCanvas("combined_fit", (std::string("Combined fit of the ") + h1->GetTitle() + " and " + h2->GetTitle()).c_str()) ;
  c->Divide(2,1);
  c->cd(1);
	gPad->SetLeftMargin(0.15) ; 
	gPad->SetLogy();
  f1.frame->Draw();
  c->cd(2);
	gPad->SetLeftMargin(0.15) ; 
	gPad->SetLogy();
  f2.frame->Draw();

  new TCanvas;
  f1.frame->Draw();
  f2.frame->Draw("same");

	RooNLLVar nll("nll","nll",simPdf,*data) ;
  RooPlot* frame_sigma = sigma.frame(Range(0.2, 5.0), Title("-log(L) scan vs sigma")) ;
  nll.plotOn(frame_sigma,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kRed)) ;
  new TCanvas;
  frame_sigma->Draw();

}


void combfit( std::list<TH1*>  & his_list)
{
  //calculate the range

  double Mmin  = (*std::max_element(std::begin(his_list), std::end(his_list), 
      [](const TH1 * h1, const TH1 *h2) { return h1->GetXaxis()->GetXmin() < h2->GetXaxis()->GetXmin(); } ))->GetXaxis()->GetXmin();

  double Mmax  = (*std::min_element(std::begin(his_list), std::end(his_list), 
      [](const TH1 * h1, const TH1 *h2) { return h1->GetXaxis()->GetXmax() < h2->GetXaxis()->GetXmax(); } ))->GetXaxis()->GetXmax();

  std::cout << "Mmin = " << Mmin << "  Mmax = " << Mmax << std::endl;

  RooRealVar Mrec("Mrec","M_{rec}(#pi^{+}#pi^{-})", Mmin,Mmax, "MeV");

  RooRealVar sigma("sigma","sigma",1.4,0,  10, "MeV") ;
  RooRealVar staple[] = 
  {
    RooRealVar("staple1","staple1",  -15,   Mmin,Mmax, "MeV"),
    RooRealVar("staple2","staple2",  -3,   Mmin,Mmax, "MeV") ,
    RooRealVar("staple3","staple3",  -2,   Mmin,Mmax, "MeV") ,
    RooRealVar("staple4","staple4",  0,    2, 2, "MeV")      ,
    RooRealVar("staple5","staple5",  +2,   Mmin,Mmax, "MeV") ,
    RooRealVar("staple6","staple6",  +4,   Mmin,Mmax, "MeV") ,
    RooRealVar("staple7","staple7",  +15,   Mmin,Mmax, "MeV")
  };
  RooRealVar n1("n1","n1", 3,  1,100) ;
  RooRealVar n2("n2","n2", 2,  1,100) ;

  //this function describes signal
  RooMcbPdf * mcbPdf =  new RooMcbPdf("ModCB", "Simple Modified CrystalBall: gaus + power + exp",  
      Mrec,  
      sigma,  
      staple[0], 
      staple[1], 
      staple[2], 
      staple[3], 
      staple[4], 
      staple[5], 
      staple[6], 
      n1, 
      n2);
  RooCategory sample("sample","sample");
  std::vector< RooFitItem_t *> fi_lst;
  std::map<std::string, RooDataHist*> dataMap;
  for(auto h : his_list)
  {
    RooFitItem_t * f = new RooFitItem_t(h, mcbPdf, Mrec,Mmin,Mmax);
    fi_lst.push_back(f);
    auto name = f->name();
    sample.defineType(name.c_str());
  }
  RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;
  for(auto f : fi_lst)
  {
    auto name = f->his->GetName();
    simPdf.addPdf(*(f->addPdf),name) ;
    dataMap[name] = f->data;
  }
  RooDataHist * data = new RooDataHist("data","Combined data", Mrec, sample, dataMap);
	simPdf.fitTo(*data, Extended(), Strategy(2), Minos());
  auto frame = Mrec.frame(Title("title"));
  std::vector<int> colors ={kBlack, kBlue, kRed, kGreen};
  //TLegend * legend = new TLegend(0.8,0.8,1.0,1.0);
  for(int i =0;i< fi_lst.size(); i++)
  {
    auto f = fi_lst[i];
    data->plotOn(frame, MarkerSize(0.5),  Cut((std::string("sample==sample::")+f->his->GetName()).c_str())) ;
    simPdf.plotOn(frame, Slice(sample,f->his->GetName()),ProjWData(sample,*data), LineWidth(1),LineColor(colors[i])) ;
    simPdf.plotOn(frame, Slice(sample,f->his->GetName()),ProjWData(sample,*data), Components(*f->bgPdf),LineStyle(kDashed),LineWidth(1),LineColor(colors[i])) ;
  }
	frame->SetMinimum(0.1);
  TCanvas * c = new TCanvas;
  c->SetLogy();
  frame->Draw();
  sigma.Print() ;
  n1.Print();
  n2.Print();
  for(int i=0;i<7;i++) staple[i].Print();
  for(auto f: fi_lst)
  {
    f->Nsig->Print();
    f->Nbg->Print();
  } 
	RooNLLVar nll("nll","nll",simPdf,*data) ;
  RooPlot* frame_sigma = sigma.frame(Range(0.2, 5.0), Title("-log(L) scan vs sigma")) ;
  nll.plotOn(frame_sigma,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kRed)) ;
  new TCanvas;
  frame_sigma->Draw();

}
*/

void combfit2( std::list<TH1*>  & his_list)
{
  //calculate the range
  double Mmin  = (*std::max_element(std::begin(his_list), std::end(his_list), 
      [](const TH1 * h1, const TH1 *h2) { return h1->GetXaxis()->GetXmin() < h2->GetXaxis()->GetXmin(); } ))->GetXaxis()->GetXmin();

  double Mmax  = (*std::min_element(std::begin(his_list), std::end(his_list), 
      [](const TH1 * h1, const TH1 *h2) { return h1->GetXaxis()->GetXmax() < h2->GetXaxis()->GetXmax(); } ))->GetXaxis()->GetXmax();

  std::cout << "Mmin = " << Mmin << "  Mmax = " << Mmax << std::endl;

  RooRealVar Mrec("Mrec", "M_{rec}(#pi^{+}#pi^{-})", Mmin, Mmax, "MeV");
  RooPlot *  frame = Mrec.frame(Title("#pi^{+}#pi^{-} recoil mass")) ;

  RooRealVar sigma("sigma", "sigma",              1.39,    0.5,    10, "MeV") ;
  RooRealVar  mean( "mean",  "mean",   -0.1 + 0.5*(Mmin+Mmax) ,  Mmin,  Mmax, "MeV") ;
    RooRealVar staple1("L1" , "Left Gaus range"   , 2        , "MeV");
    RooRealVar staple2("L2" , "Left Exp range"    , 0.853    , 0      , (Mmax-Mmin)*0.5 , "MeV");
    RooRealVar staple3("L3" , "Left Power range"  , 12.8861  , 0      , (Mmax-Mmin)*0.5 , "MeV");
    RooRealVar staple4("L4" , "Left Exp range 2"  , 22.7831  , 0      , (Mmax-Mmin)*0.5 , "MeV");
    RooRealVar staple5("R1" , "Right Gaus range"  , 2        , "MeV");
    RooRealVar staple6("R2" , "Right Exp range"   , 2.6e-7   , 0      , (Mmax-Mmin)*0.5 , "MeV"); 
    RooRealVar staple7("R3" , "Right Power range" , 31.9138  , 0      , (Mmax-Mmin)*0.5 , "MeV"); 
    RooRealVar staple8("R4" , "Right Exp range 2" , 0.125333 , 0      , (Mmax-Mmin)*0.5 , "MeV");

    std::vector<RooRealVar*> staple = {&staple1,&staple2,&staple2,&staple3,&staple4,&staple5,&staple6,&staple7,&staple8};

  RooRealVar n1("Ln", "Left power",  1.326,  1, 100) ;
  RooRealVar n2("Rn", "Right power", 1.507,  1, 100) ;

  RooMcb2Pdf * mcbPdf =  new RooMcb2Pdf("ModCB", "Modified CrystalBall: gaus + exp + power + exp",  
      Mrec,  
      mean,  
      sigma,  
      staple1, 
      staple2, 
      staple3, 
      staple4, 
      staple5, 
      staple6, 
      staple7, 
      staple8, 
      n1, 
      n2,
      Mmin,
      Mmax);

  //this is needed for fit of the monte carlo UU
  RooRealVar sigma2("sigma2", "sigma2",   5,    2,    20, "MeV") ;
  RooRealVar  mean2( "mean2",  "mean2",   20,   (Mmax+Mmin)*0.5+10,  Mmax, "MeV") ;

  RooRealVar sigma3("sigma3", "sigma3",   5,    2,    20, "MeV") ;
  RooRealVar  mean3( "mean3",  "mean3",   10,   (Mmax+Mmin)*0.5+2,  Mmax, "MeV") ;

  RooGaussian rad_gausPdf("rad_gaus","Radiative gauss for MonteCarlo", Mrec, mean2, sigma2);
  RooGaussian rad_gausPdf3("rad_gaus3","Radiative gauss for MonteCarlo", Mrec, mean3, sigma3);

  RooRealVar  rad_gaus_fraction("rgfrac","fraction of radiative gaus",0,0.,1.);
  RooRealVar  rad_gaus_fraction2("rgfrac2","fraction of radiative gaus",0,0.,1.);

  RooAddPdf * mcbRadPdf = new RooAddPdf("mcb_rad","Modified crystal ball + radiative Gaus",
      RooArgList(*mcbPdf, rad_gausPdf, rad_gausPdf3), RooArgList(rad_gaus_fraction, rad_gaus_fraction2)); 


  RooAbsPdf * modelPdf = mcbPdf;
  if ( GAUS_RAD ) modelPdf = mcbRadPdf;

  RooCategory sample("sample","sample");
  std::vector< RooFitItem_t *> fi_lst;
  std::map<std::string, RooDataHist*> dataMap;
  for(auto h : his_list)
  {
    RooFitItem_t * f = new RooFitItem_t(h, modelPdf, Mrec,Mmin,Mmax);
    f->frame = frame;
    fi_lst.push_back(f);
    auto name = f->name();
    sample.defineType(name.c_str());
  }
  RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;
  for(auto f : fi_lst)
  {
    auto name = f->his->GetName();
    simPdf.addPdf(*(f->addPdf),name) ;
    dataMap[name] = f->data;
  }
  RooDataHist * data = new RooDataHist("data","Combined data", Mrec, sample, dataMap);

	auto theFitResult = simPdf.fitTo(*data, Extended(), Strategy(2), Minos(), Save());
  auto p = simPdf.getParameters(Mrec);
  p->writeToFile("tmp_fit_result.txt");
	//simPdf.fitTo(*data, Extended(), Strategy(2), Minos());

  int nfp = theFitResult->floatParsFinal().getSize() ;
  std::cout << "nfp = " << nfp << std::endl;

	RooChi2Var chi2Var("chi2", "chi2", simPdf, *data);
  double chi2 = chi2Var.getValV();
  //calculate number of degree of freadom
  int ndf=0;
  int nfree_param = 10; //mean,sigma, n1,n2, 6 staple intervals
  for(auto h : his_list)
  {
    std::cout << h->GetName() << " " << h->GetNbinsX() << std::endl;
    ndf+=h->GetNbinsX();
    nfree_param += 3; //Nsig, Nbg, bg_slope
  }
  double chi2ndf = chi2/(ndf-nfree_param);
  double chi2prob = TMath::Prob(chi2,ndf - nfree_param);
  std::cout << boost::format("chi2/ndf = %f/(%d-%d) = %f,  prob = %f") % chi2 % ndf % nfree_param % chi2ndf %  chi2prob<< std::endl;

  std::vector<int> colors ={kBlack, kBlue, kRed, kGreen};
  TLegend * legend = new TLegend(0.8,0.8,1.0,1.0);
  RooArgSet viewArgSet;
  ///viewArgSet.add(myChi2);
  for(int i =0;i< fi_lst.size(); i++)
  {
    auto f = fi_lst[i];
    data->plotOn(frame, XErrorSize(0),MarkerSize(0.5),  Cut((std::string("sample==sample::")+f->his->GetName()).c_str()), LineColor(colors[i]), MarkerColor(colors[i])) ;
    //example of plotting statistics
    //data->statOn(frame,Layout(0.55,0.99,0.7));
    simPdf.plotOn(frame, Slice(sample,f->his->GetName()),ProjWData(sample,*data), LineWidth(1),LineColor(colors[i])) ;
    simPdf.plotOn(frame, Slice(sample,f->his->GetName()),ProjWData(sample,*data), Components(*f->bgPdf),LineStyle(kDashed),LineWidth(1),LineColor(colors[i])) ;
    viewArgSet.add(*(f->Nsig));
    viewArgSet.add(*(f->Nbg));
  }
	frame->SetMinimum(0.1);
  //chi2 = frame->chiSquare("simPdf","data",nfree_param);
  //double ndoff = frame->GetNbinsX();
  //chi2prob = TMath::Prob(chi2,ndoff);
  //std::cout << " chi2 = " << chi2 << " ndoff = " << ndoff  << "  chi2prob = " << chi2prob << std::endl;
  auto chi_fmt = boost::format("#chi^{2}/ndf = %.1f/%d = %.2f") % chi2 % (ndf-nfree_param) % chi2ndf;
  auto prob_fmt = boost::format("prob = %.1f%%") % chi2prob;
  simPdf.paramOn(frame, Label(chi_fmt.str().c_str()), Parameters(viewArgSet));
  //Mrec.plotOn(frame);
  TCanvas * c = new TCanvas;
  c->SetLogy();
  frame->Draw();
  sigma.Print();
  mean.Print();
  n1.Print();
  n2.Print();
  for(auto & s: staple) s->Print();
  for(auto f: fi_lst)
  {
    f->Nsig->Print();
    f->Nbg->Print();
  } 


  //print Nsig's ratio for every channel
  for(auto f: fi_lst)
  {
    auto fm = fi_lst.back();
    if(fm == f) continue;
    double theRatio = f->Nsig->getValV()/fm->Nsig->getValV();
    double theRelError = sqrt( pow(f->Nsig->getError()/f->Nsig->getValV(), 2)  +  pow(fm->Nsig->getError()/fm->Nsig->getValV(), 2) );
    double theError = theRatio * theRelError;
    boost::format fmt("%s/%s = (%6.3f ± %-4.3f) * 1e-3 ( ± %.2f%%)");
    cout << fmt % f->Nsig->GetName() % fm->Nsig->GetName()% (theRatio*1000) %  (theError*1000) %  (theRelError*100)  <<  endl;
  }

	RooNLLVar nll("nll","nll",simPdf,*data) ;
  std::cout << "nllvar = " << nll.getValV() << std::endl;
  RooPlot* frame_mean = mean.frame(Range(Mmin, Mmax), Title("-log(L) scan vs mean")) ;
  nll.plotOn(frame_mean,LineColor(kRed)) ;
  RooPlot* frame_n1 = n1.frame(Title("n1 -log(L)"), Range(1, 5)) ;
  nll.plotOn(frame_n1,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kBlue)) ;

	TCanvas * cnll =new TCanvas("cnll", "- Log likelihood");
  frame_mean->Draw();
}

int main(int argc,  char ** argv)
{
  std::ios_base::sync_with_stdio(false);
  namespace po=boost::program_options;
  po::options_description opt_desc("Allowed options");
  std::string tree_name;
  std::string mctree_name;
  std::string file_name;
  std::string his_name;
  std::string cut;
  std::string hash_list_str;
  std::string suffix;
  std::string combine_str;
  double nbg=0;
  opt_desc.add_options()
    ("help,h","Print this help")
    ("input", po::value<std::string>(&file_name), "Root file (.root) with the data")
    ("his", po::value<std::string>(&his_name), "Histogram name")
    ("tree", po::value<std::string>(&tree_name), "Tree name ")
    ("mctree", po::value<std::string>(&mctree_name), "MonteCarlo Tree name")
    ("cut", po::value<std::string>(&cut), "Cut")
    ("hash",po::value<std::string>(&hash_list_str),"List of hashes to draw")
    ("skip","skip hashes")
    ("Nbg",po::value<double>(&nbg),"Background events")
    ("suffix,s",po::value<std::string>(&suffix),"Suffix")
    ("mc", "Use MonteCarlo information mctopo")
    ("trk","Tracking efficiency fit")
    ("combine,c",po::value<std::string>(&combine_str), "combined fit")
    ("nobgslope","No bg slope")
    ("nobg","No background events")
    ("nograd","No radiative gaus")
    ("gaus_rad","Gaus rad")
    ;
  po::positional_options_description pos;
  pos.add("input", 1);
  po::variables_map opt; //options container
  try
  {
    po::store(po::command_line_parser(argc, argv).options(opt_desc).positional(pos).run(), opt);
    po::notify(opt);
  } 
  catch (boost::program_options::error & po_error)
  {
    std::cerr << "WARGNING: configuration: "<< po_error.what() << std::endl;
  }

  if(opt.count("help") && !opt.count("input"))
  {
    std::cout << "Usage: testMcbFit <root_file>" << std::endl;
    std::clog << opt_desc;
    return 0;
  }

  OPT_NOBGSLOPE=opt.count("nobgslope");
  OPT_NOBG = opt.count("nobg");
  OPT_NOGAUSRAD = opt.count("nograd");


  BG_NOSLOPE = opt.count("nobgslope");
  BG_FIX_TO_ZERO = opt.count("nobg");
  GAUS_RAD = opt.count("gaus_rad");

  std::cout << " " << suffix << " " << file_name << std::endl;
  
  TCut hash_cut = "";

  if(opt.count("hash"))
  {
    hash_cut = make_cut (hash_list_str);
    if(opt.count("skip"))
    {
      hash_cut = ! hash_cut;
    }
    std::cout << "Using hash cut =" << hash_cut << std::endl;
  }

  TFile file(file_name.c_str());
  TH1    * his=nullptr;
  TTree * tree=nullptr;

	TApplication theApp("root_app", &argc, argv);

  auto create_histogram = [&](std::string suffix) -> TH1*
  {
    tree_name   = "event"+suffix;
    his_name    = "hMrec"+suffix;
    mctree_name = "mctopo"+suffix;
    std::cout << tree_name << "  " << his_name << "  " << mctree_name << std::endl;
    auto h = (TH1*)file.Get(his_name.c_str());
    std::string title = h->GetTitle();
    if( h == nullptr ) 
    {
      std::cerr << "ERROR: Unable to find histogram " << his_name << std::endl;
      exit(1);
    }
    tree = (TTree*) file.Get(tree_name.c_str());
    if ( tree == nullptr )
    {
      std::cerr << "ERROR: Unable to find histogram " << tree_name << std::endl;
      exit(1);
    }

    TTree *mctree = (TTree*) file.Get(mctree_name.c_str());
    if ( mctree == nullptr )
    {
      std::cerr << "ERROR: Unable to find histogram " << mctree_name << std::endl;
      exit(1);
    }
    if(mctree->GetEntries() > 0)
    {
      std::cout << "MonteCarlo data for " << mctree_name << "  exists" << std::endl;
      if(mctree->GetEntries() != tree->GetEntries())
      {
        std::cerr << "ERROR: number of events in " << mctree_name <<  " (" << mctree->GetEntries() << ") " << " doesnt fit to number of entries in main tree " << tree_name  << " (" << tree->GetEntries() << ")" << std::endl;
      }
      else 
      {
        if(mctree) tree->AddFriend(mctree);
      }
    }
    int    Nbin = h->GetNbinsX();
    double Mmin = h->GetXaxis()->GetXmin();
    double Mmax = h->GetXaxis()->GetXmax();
    auto varexp = boost::format("(Mrec-3.097)*1000 >> Mrec%s(%d,%f,%f)") % suffix % Nbin % Mmin % Mmax; 
    tree->Draw(varexp.str().c_str(),hash_cut && cut.c_str(),"goff");
    h = (TH1F*) tree->GetHistogram();
    h->SetName(suffix.c_str());
    h->SetTitle(title.c_str());
    return h;
  };

  if(opt.count("suffix"))
  {
    std::vector<std::string> suffix_list;
    boost::split(suffix_list,suffix,boost::is_any_of(", "));
    std::list<TH1*> his_lst;
    for(auto str : suffix_list)
    {
      his_lst.push_back(create_histogram(str));
    }
    //combfit2(his_lst);
    //combfit(his_lst.front(), his_lst.back());
    fit(his_lst);
    theApp.Run();
  }
	theApp.Run();
}
