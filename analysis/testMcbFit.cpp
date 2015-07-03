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
#include "RooMcbPdf.h"

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

using namespace RooFit;
using namespace std;
int main(int argc,  char ** argv)
{
	TApplication theApp("root_app", &argc, argv);
	TFile file("data12-Mrec90MeV.root");
	TH1 * his = (TH1*)file.Get("hMrecUU");
	his->Draw();
	//theApp.Run();
  // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
	double Mmin  = his->GetXaxis()->GetXmin();
	double Mmax  = his->GetXaxis()->GetXmax();
  RooRealVar Mrec("Mrec","Mrec",Mmin,Mmax);
  RooRealVar sigma("sigma","sigma",his->GetRMS(),0,100, "MeV") ;
  RooRealVar staple1("staple1","staple1",-20,   Mmin,Mmax, "MeV") ;
  RooRealVar staple2("staple2","staple2", -2,  Mmin,Mmax, "MeV") ;
  RooRealVar staple3("staple3","staple3", -1.75,  Mmin,Mmax, "MeV") ;
  RooRealVar staple4("staple4","staple4",  0,   Mmin,Mmax, "MeV") ;
  RooRealVar staple5("staple5","staple5",  1.2,  Mmin,Mmax, "MeV") ;
  RooRealVar staple6("staple6","staple6",  2,  Mmin,Mmax, "MeV") ;
  RooRealVar staple7("staple7","staple7", 20,   Mmin,Mmax, "MeV") ;
  RooRealVar n1("n1","n1", 3.1,  1,100) ;
  RooRealVar n2("n2","n2", 2.4,  1,100) ;

	std::vector<RooRealVar> staple(7);
	for(int i=0;i<staple.size();i++)
	{
		char buf[128];
		sprintf(buf, "staple-%d",i );
		staple[i] = RooRealVar(buf, buf, (i-3.5)*5,  -40, +40, "MeV");
	}
	std::vector<RooRealVar> N(2);
	N[0] = RooRealVar("nl", "Left power", 2,  0.1, 100);
	N[1] = RooRealVar("nr", "Right power", 2,  0.1, 100);

	//ModifiedCrystalBall mcb("ModCB", "Modified CrystalBall",  Mrec,  sigma,  
	//		staple1, 
	//		staple2, 
	//		staple3, 
	//		staple4, 
	//		staple5, 
	//		staple6, 
	//		staple7, 
	//		n1, 
	//		n2);
	//RooGaussian mcb("test_gaus", "test_gaus",  Mrec,  mean,  sigma);
	//RooMcbPdf mcb("argus", "argus", Mrec, staple1, sigma);
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
	cout << "After mcb" << endl;

	RooRealVar mean2("mean2", "mean2",0,  -5, 5, "MeV");
	RooRealVar sigma2("sigma2", "sigma2", 0.5,  0, 3, "MeV");
	RooGaussian small_gaus("smals_gauss", "small_gauss",  Mrec,  mean2,  sigma2);


	RooRealVar mcb_yield("mcb_yield", "MCB yield", his->GetEntries(), 0, 1000000);
	RooRealVar sg_yield("sg_yield", "small gauss yield", 0, 0, 1000000);
	RooArgList sig_shapes;
	RooArgList sig_yields;
	sig_shapes.add(small_gaus);      sig_yields.add(sg_yield);
	sig_shapes.add(mcb);             sig_yields.add(mcb_yield);
	RooAddPdf complexSignal("complexSignal", "signal with small gaus",  sig_shapes,  sig_yields);

	//RooRealVar poly_c1("poly_c1","coefficient of x^1 term",0, -10, 10);
	RooRealVar poly_c1("poly_c1","coefficient of x^1 term",0);
	RooPolynomial bkgd_poly("bkgd_poly", "linear function for background", Mrec, RooArgList(poly_c1));	

  // Construct plot frame in 'x'
  //RooPlot* xframe = Mrec.frame(Title("Gaussian p.d.f.")) ;

	RooRealVar peak_yield("peak_yield", "yield signal peak", his->GetEntries(), 0, 1e100);
	RooRealVar bkgd_yield("bkgd_yield", "yield of background", 0, 0, 1e100);
	RooArgList shapes;
	RooArgList yields;
	shapes.add(bkgd_poly);      yields.add(bkgd_yield);
	shapes.add(mcb);  yields.add(peak_yield);
	RooAddPdf  totalPdf("totalPdf", "sum of signal and background PDF's", shapes, yields);
	cout << "After total pdf" << endl;


	RooAbsReal* igx = mcb.createIntegral(Mrec) ;
	cout << "gx_Int[x] = " << igx->getVal() << endl;


  // P l o t   m o d e l   a n d   c h a n g e   p a r a m e t e r   v a l u e s
  // ---------------------------------------------------------------------------

  // Plot gauss in frame (i.e. in x) 
  //gauss.plotOn(xframe) ;

  // Change the value of sigma to 3
  //sigma.setVal(3) ;

  // Plot gauss in frame (i.e. in x) and draw frame on canvas
  //gauss.plotOn(xframe,LineColor(kRed)) ;
  

	// Create category observable c that serves as index for the ROOT histograms
  //RooCategory rooCategory("c","c") ;
  //rooCategory.defineType("his") ;
  // G e n e r a t e   e v e n t s 
  // -----------------------------

  // Generate a dataset of 1000 events in x from gauss
  //RooDataSet* data = gauss.generate(x,10000) ;  
	RooArgSet ntupleVarSet(Mrec);
	//RooDataSet * data = new RooDataSet("tree", "tree",  tree,  ntupleVarSet);
	//RooDataSet * data = new RooDataSet("dh", "dh", Mrec, Index(rooCategory),  Import("his",*his));
	RooDataHist * data = new RooDataHist("dh", "dh", Mrec, Import(*his));
	cout << "After import his" << endl;
  
	totalPdf.fitTo(*data,  FitOptions("qmh"));
  // Make a second plot frame in x and draw both the 
  // data and the p.d.f in the frame
  RooPlot* xframe2 = Mrec.frame(Title("Gaussian p.d.f. with data")) ;
  data->plotOn(xframe2, MarkerSize(0.5)) ;
	cout << "After plot on frame" << endl;;
  totalPdf.plotOn(xframe2,  LineWidth(1)) ;
	cout << "After total pdf plot" << endl;;
	totalPdf.plotOn(xframe2,Components(bkgd_poly),LineStyle(kDashed),  LineWidth(1)) ;
	//mcb.plotOn(xframe2,Components(bkgd_poly),LineStyle(kDashed)) ;
	//mcb.plotOn(xframe2) ;
	cout << "Components plot" << endl;;
  

  // F i t   m o d e l   t o   d a t a
  // -----------------------------

  // Fit pdf to data
  //gauss.fitTo(*data) ;
	//totalPdf.fitTo(*data, Extended());
	//cout << "Before fit" << endl;
	//totalPdf.fitTo(*data);
	////mcb.fitTo(*data, FitOptions("qmh"));
	//cout << "After fit" << endl;

  // Print values of mean and sigma (that now reflect fitted values and errors)
  //mean.Print() ;
  sigma.Print() ;
	//for(auto s : staple)
	//{
	//	s.Print();
	//}
	//for(auto n : N)
	//{
	//	n.Print();
	//}

  // Draw all frames on a canvas
  TCanvas* c = new TCanvas("rf101_basics","rf101_basics",800,400) ;
	c->SetLogy();
  //c->Divide(2) ;
  //c->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe->GetYaxis()->SetTitleOffset(1.6) ; xframe->Draw() ;
  c->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe2->GetYaxis()->SetTitleOffset(1.6) ; xframe2->Draw() ;





	theApp.Run();
}
