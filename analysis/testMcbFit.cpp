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


#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>
extern void InitGui();
VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
TROOT root("polarimeter","polarimeter", initfuncs);


#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>

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

#include "RooMcbPdf.h"

using namespace RooFit;
using namespace std;
int main(int argc,  char ** argv)
{
  if(argc<2) return 1;
	TFile file(argv[2]);
	TH1 * his = (TH1*)file.Get(argv[1]);
	TApplication theApp("root_app", &argc, argv);

	double Mmin  = his->GetXaxis()->GetXmin();
	double Mmax  = his->GetXaxis()->GetXmax();
  RooRealVar Mrec("Mrec","M_{rec}(#pi^{+}#pi^{-})", Mmin,Mmax, "MeV");
  RooRealVar sigma("sigma","sigma",1.4,0,10, "MeV") ;
  RooRealVar staple1("staple1","staple1",-15,   Mmin,Mmax, "MeV") ;
  RooRealVar staple2("staple2","staple2", -3,   Mmin,Mmax, "MeV") ;
  RooRealVar staple3("staple3","staple3", -2,   Mmin,Mmax, "MeV") ;
  RooRealVar staple4("staple4","staple4",  0,   -2,2, "MeV") ;
  RooRealVar staple5("staple5","staple5",  2,   Mmin,Mmax, "MeV") ;
  RooRealVar staple6("staple6","staple6",  4,   Mmin,Mmax, "MeV") ;
  RooRealVar staple7("staple7","staple7", 15,   Mmin,Mmax, "MeV") ;
  RooRealVar n1("n1","n1", 3,  1,100) ;
  RooRealVar n2("n2","n2", 2,  1,100) ;


	std::vector<RooRealVar> staple(7);
	for(int i=0;i<staple.size();i++)
	{
		char buf[128];
		sprintf(buf, "staple-%d",i );
		staple[i] = RooRealVar(buf, buf, (i-3.5)*5,  -40, +40, "MeV");
	}
	std::vector<RooRealVar> N(2);
	N[0] = RooRealVar("nl", "Left power",  2,  0.1, 100);
	N[1] = RooRealVar("nr", "Right power", 2,  0.1, 100);

	RooMcbPdf mcb("ModCB", "Modified CrystalBall",  Mrec,  sigma,  
			staple1, 
			staple2, 
			staple3, 
			staple4, 
			staple5, 
			staple6, 
			staple7, 
			n1, 
			n2);

	RooRealVar mean2("mean2", "mean2",0,  -5, 5, "MeV");
	RooRealVar sigma2("sigma2", "sigma2", 0.5,  0, 3, "MeV");
	RooGaussian small_gaus("smals_gauss", "small_gauss",  Mrec,  mean2,  sigma2);


	RooRealVar mcb_yield("mcb_yield", "MCB yield", his->GetEntries(), 0, his->GetEntries());
	RooRealVar sg_yield("sg_yield", "small gauss yield", 0, 0, his->GetEntries());
	RooArgList sig_shapes;
	RooArgList sig_yields;
	sig_shapes.add(small_gaus);      sig_yields.add(sg_yield);
	sig_shapes.add(mcb);             sig_yields.add(mcb_yield);
	RooAddPdf complexSignal("complexSignal", "signal with small gaus",  sig_shapes,  sig_yields);

	RooRealVar poly_c1("poly_c1","coefficient of x^1 term",0, -10, 10);
	//RooRealVar poly_c1("poly_c1","coefficient of x^1 term",0);
	RooPolynomial bkgd_poly("bkgd_poly", "linear function for background", Mrec, RooArgList(poly_c1));	


	RooRealVar Nsig("Nsig", "Number of signal events", his->GetEntries(), 0, his->GetEntries()*100);
	RooRealVar Nbg("Nbg", "Number of background events", 0, 0,  his->GetEntries());
	//RooRealVar bkgd_yield("bkgd_yield", "yield of background", 0);
	RooArgList shapes;
	RooArgList yields;
	shapes.add(bkgd_poly);      yields.add(Nbg);
	shapes.add(mcb);  					yields.add(Nsig);
	RooAddPdf  totalPdf("totalPdf", "sum of signal and background PDF's", shapes, yields);


	RooAbsReal* igx = mcb.createIntegral(Mrec) ;
	cout << "gx_Int[x] = " << igx->getVal() << endl;


	RooArgSet ntupleVarSet(Mrec);
	//RooDataSet * data = new RooDataSet("tree", "tree",  tree,  ntupleVarSet);
	//RooDataSet * data = new RooDataSet("dh", "dh", Mrec, Index(rooCategory),  Import("his",*his));
	RooDataHist * data = new RooDataHist("dh", "dh", Mrec, Import(*his));
  
  RooPlot* xframe3 = Mrec.frame(Title("Show my cb function")) ;
	mcb.plotOn(xframe3);
	TCanvas * c_mcb = new TCanvas;
	xframe3->Draw();

	//theApp.Run();
	//totalPdf.fitTo(*data,  Extended(), FitOptions("qmh"));
	totalPdf.fitTo(*data,  Extended(), Strategy(2));

	RooChi2Var chi2Var("chi2", "chi2", totalPdf, *data);
	cout << "chi2/ndf = " << chi2Var.getVal() << endl;

  RooPlot* xframe2 = Mrec.frame(Title("Fit by Modified CrystalBall")) ;
  data->plotOn(xframe2, MarkerSize(0.5)) ;
  totalPdf.plotOn(xframe2,  LineWidth(1)) ;
	totalPdf.plotOn(xframe2,Components(bkgd_poly),LineStyle(kDashed),  LineWidth(1)) ;
  

  sigma.Print() ;
	Nsig.Print();
	Nbg.Print();
	//for(auto s : staple)
	//{
	//	s.Print();
	//}
	//for(auto n : N)
	//{
	//	n.Print();
	//}

  // Draw all frames on a canvas
  TCanvas* c = new TCanvas("mcb","Test for Modified CrystalBall Fit") ;
	//TH2F * h2=new TH2F("h2", "h2", 100, Mmin, Mmax, 100, -0.2*his->GetMaximum(), his->GetMaximum()*1.2);
	//h2->Draw();
	c->SetLogy();
  c->cd(2) ; 
	gPad->SetLeftMargin(0.15) ; 
	xframe2->GetYaxis()->SetTitleOffset(1.6) ; 
	xframe2->Draw() ;



	RooNLLVar nll("nll","nll",totalPdf,*data) ;

  RooPlot* frame_sigma = sigma.frame(Range(0.2, 2.0), Title("-log(L) scan vs sigma")) ;
  nll.plotOn(frame_sigma,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kRed)) ;
  RooPlot* frame_staple4 = staple4.frame(Title("staple4 -log(L)")) ;
  nll.plotOn(frame_staple4,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kRed)) ;

  RooPlot* frame_staple3 = staple3.frame(Title("staple3 -log(L)")) ;
  nll.plotOn(frame_staple3,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kRed)) ;
  RooPlot* frame_staple5 = staple5.frame(Title("staple5 -log(L)")) ;
  nll.plotOn(frame_staple5,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kBlue)) ;

  RooPlot* frame_staple2 = staple2.frame(Title("staple2 -log(L)")) ;
  nll.plotOn(frame_staple2,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kBlue)) ;
  RooPlot* frame_staple6 = staple6.frame(Title("staple6 -log(L)")) ;
  nll.plotOn(frame_staple6,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kRed)) ;

  RooPlot* frame_staple1 = staple1.frame(Title("staple1 -log(L)")) ;
  nll.plotOn(frame_staple1,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kBlue)) ;
  RooPlot* frame_staple7 = staple7.frame(Title("staple7 -log(L)")) ;
  nll.plotOn(frame_staple7,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kRed)) ;


  RooPlot* frame_n1 = n1.frame(Title("n1 -log(L)"), Range(1, 5)) ;
  nll.plotOn(frame_n1,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kBlue)) ;

  RooPlot* frame_n2 = n2.frame(Title("n2 -log(L)"), Range(1, 5)) ;
  nll.plotOn(frame_n2,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll.getVal()+10),LineColor(kRed)) ;

  RooPlot* frame_Nsig = Nsig.frame(Range(0, his->GetEntries()*1.2), Title("Nsig -log(L)")) ;
  nll.plotOn(frame_Nsig,PrintEvalErrors(-1),ShiftToZero(), EvalErrorValue(nll.getVal()+10),LineColor(kRed)) ;
	frame_Nsig->SetMinimum(0);
	frame_Nsig->SetMaximum(40);

  RooPlot* frame_Nbg = Nbg.frame(Range(0, Nbg.getValV()*2), Title("Nbg -log(L)")) ;
  nll.plotOn(frame_Nbg,PrintEvalErrors(-1),ShiftToZero(), EvalErrorValue(nll.getVal()+10),LineColor(kRed)) ;

	TCanvas * cnll =new TCanvas;
	cnll->Divide(3, 3);
	cnll->cd(1);
	frame_sigma->Draw();
	cnll->cd(2);
	frame_staple4->Draw();
	cnll->cd(3);
	frame_staple3->Draw();
	frame_staple5->Draw("same");
	cnll->cd(5);
	frame_staple2->Draw();
	frame_staple6->Draw("same");
	cnll->cd(6);
	frame_staple1->Draw();
	frame_staple7->Draw("same");
	cnll->cd(7);
	frame_n1->Draw();
	frame_n2->Draw("same");
	cnll->cd(8);
	frame_Nsig->Draw();
	cnll->cd(9);
	frame_Nbg->Draw();

	for(int i=0;i<his->GetNbinsX();i++)
	{
		double x  = his->GetBinCenter(i);
		double N  = his->GetBinContent(i);
		double Nth;
		//double Nth = mcb.evaluate(); 
		//cout <<  i << " " << x << "  " << N <<   " " << Nth << endl;
	}
	
	theApp.Run();
}
