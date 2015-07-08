// =====================================================================================
//
//       Filename:  CombMcbFit
//
//    Description:  Simultaneous fit of KK and uu events with the same signal shape
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
#include <TH2F.h>

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
	TFile file(argv[1]);
	TH1 * hisKK = (TH1*)file.Get("hMrecKK");
	TH1 * hisUU = (TH1*)file.Get("hMrecUU");
	//TTree * treeKK = (TTree*)file.Get("eventKK");
	//TTree * treeUU = (TTree*)file.Get("eventUU");
	TApplication theApp("root_app", &argc, argv);

	double Mmin  = hisKK->GetXaxis()->GetXmin();
	double Mmax  = hisKK->GetXaxis()->GetXmax();
  RooRealVar Mrec("Mrec","M_{rec}(#pi^{+}#pi^{-}) - 3097 ", Mmin,Mmax, "MeV");
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


	//RooRealVar poly_c1("poly_c1","coefficient of x^1 term",0, -10, 10);
	//RooRealVar poly_c1("poly_c1","coefficient of x^1 term",0);
	RooRealVar bgKK_c1("bgKK_c1","coefficient of Mrec term for KK channel",0, -10, 10);
	//RooPolynomial bgKK("bgKK", "linear function for KK background", Mrec, RooArgList(bgKK_c1));	
	RooPolynomial bgKK("bgKK", "linear function for KK background", Mrec);	
	RooPolynomial bgUU("bgKK", "constant background for UU channel", Mrec);	


	RooRealVar NKK("NKK", "Number of KK events", hisKK->GetEntries(), 0, hisKK->GetEntries()*100);
	RooRealVar NbgKK("NbgKK", "Number of background events for KK channel", 0, 0,  hisKK->GetEntries()*100);
	//RooRealVar bkgd_yield("bkgd_yield", "yield of background", 0);
	RooArgList shapesKK;
	RooArgList yieldsKK;
	shapesKK.add(bgKK);         yieldsKK.add(NbgKK);
	shapesKK.add(mcb);  					yieldsKK.add(NKK);
	RooAddPdf  KKPdf("KKPdf", "KK + background p.d.f.", shapesKK, yieldsKK);

	RooRealVar NUU("NUU", "Number of UU events", hisUU->GetEntries(), 0, hisUU->GetEntries()*100);
	RooRealVar NbgUU("NbgUU", "Number of background events for UU channel", 0, 0,  hisUU->GetEntries()*100);
	RooArgList shapesUU;
	RooArgList yieldsUU;
	shapesUU.add(bgUU);      yieldsUU.add(NbgUU);
	shapesUU.add(mcb);  					yieldsUU.add(NUU);
	RooAddPdf  UUPdf("UUPdf", "UU + background p.d.f.", shapesUU, yieldsUU);




	RooAbsReal* igx = mcb.createIntegral(Mrec) ;
	cout << "gx_Int[x] = " << igx->getVal() << endl;


	//RooArgSet ntupleVarSet(Mrec);
	//RooDataSet * dataSetKK = new RooDataSet("KK", "KK",  treeKK,  ntupleVarSet);
	//RooDataSet * dataSetUU = new RooDataSet("UU", "UU",  treeUU,  ntupleVarSet);


	RooDataHist * dataKK = new RooDataHist("KK", "KK", Mrec, Import(*hisKK));
	RooDataHist * dataUU = new RooDataHist("UU", "UU", Mrec, Import(*hisUU));

 	RooCategory sample("sample","sample") ;
  sample.defineType("KK") ;
  sample.defineType("UU") ;

  RooDataHist combData("combData","combined data KK + uu channel",Mrec,Index(sample),Import("KK",*dataKK),Import("UU",*dataUU)) ;
	RooDataHist * data = & combData;

  //RooDataSet combData("combData","combined data KK + uu channel",Mrec,Index(sample),Import("KK",*dataSetKK),Import("UU",*dataSetUU)) ;
	//RooDataSet * data = & combData;


  // C o n s t r u c t   a   s i m u l t a n e o u s   p d f   i n   ( x , s a m p l e )
  // -----------------------------------------------------------------------------------

  // Construct a simultaneous pdf using category sample as index
  RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;

  // Associate model with the physics state and model_ctl with the control state
  simPdf.addPdf(KKPdf,"KK") ;
  simPdf.addPdf(UUPdf,"UU") ;
  
  RooPlot* xframe3 = Mrec.frame(Title("Show my cb function")) ;
	mcb.plotOn(xframe3);
	TCanvas * c_mcb = new TCanvas;
	xframe3->Draw();

	simPdf.fitTo(combData, Extended(), Strategy(2), Minos());

	//RooChi2Var chi2Var("chi2", "chi2", totalPdf, *data);
	//cout << "chi2/ndf = " << chi2Var.getVal() << endl;

  RooPlot* frameKK = Mrec.frame(Title("K^{+}K^{-} events")) ;
  RooPlot* frameUU = Mrec.frame(Title("#mu^{+}mu^{-} events")) ;
  combData.plotOn(frameKK, MarkerSize(0.5),  Cut("sample==sample::KK")) ;
  simPdf.plotOn(frameKK, Slice(sample,"KK"),ProjWData(sample,combData),   LineWidth(1)) ;
  simPdf.plotOn(frameKK, Slice(sample,"KK"),ProjWData(sample,combData), Components(bgKK),LineStyle(kDashed),LineWidth(1)) ;
	frameKK->SetMinimum(1);
  combData.plotOn(frameUU, MarkerSize(0.5),  Cut("sample==sample::UU")) ;
  simPdf.plotOn(frameUU, Slice(sample,"UU"),ProjWData(sample,combData),   LineWidth(1)) ;
  simPdf.plotOn(frameUU, Slice(sample,"UU"),ProjWData(sample,combData), Components(bgUU),LineStyle(kDashed),LineWidth(1)) ;
	frameUU->SetMinimum(1);
	//frameUU->SetMaximum();

  TCanvas* c = new TCanvas("mcb","Test for Modified CrystalBall Fit") ;
	c->Divide(1, 2);
	c->SetLogy();
  c->cd(1) ; 
	gPad->SetLeftMargin(0.15) ; 
	gPad->SetLogy();
	frameKK->GetYaxis()->SetTitleOffset(1.6) ; 
	frameKK->Draw() ;
  c->cd(2) ; 
	gPad->SetLogy();
	gPad->SetLeftMargin(0.15) ; 
	frameUU->GetYaxis()->SetTitleOffset(1.6) ; 
	frameUU->Draw() ;
  

  sigma.Print() ;
	NKK.Print();
	NUU.Print();
	NbgKK.Print();
	NbgUU.Print();


	RooNLLVar nll("nll","nll",simPdf,*data) ;

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

  RooPlot* frame_NKK = NKK.frame(Range(0, hisKK->GetEntries()*1.2), Title("NKK -log(L)")) ;
  nll.plotOn(frame_NKK,PrintEvalErrors(-1),ShiftToZero(), EvalErrorValue(nll.getVal()+10),LineColor(kRed)) ;
	frame_NKK->SetMinimum(0);
	frame_NKK->SetMaximum(40);

  RooPlot* frame_NbgKK = NbgKK.frame(Range(0, NbgKK.getValV()*2), Title("NbgKK -log(L)")) ;
  nll.plotOn(frame_NbgKK,PrintEvalErrors(-1),ShiftToZero(), EvalErrorValue(nll.getVal()+10),LineColor(kRed)) ;

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
	frame_NKK->Draw();
	cnll->cd(9);
	frame_NbgKK->Draw();


	double theRatio = NKK.getValV()/NUU.getValV();
	double theRelError = sqrt( pow(NKK.getError()/NKK.getValV(), 2)  +  pow(NUU.getError()/NUU.getValV(), 2) );
	double theError = theRatio * theRelError;
	cout << "Nkk/Nuu = " << theRatio << " +/- " << theError  <<  " ( +/- " << theRelError*100 << "% )" <<  endl;
	theApp.Run();
}
