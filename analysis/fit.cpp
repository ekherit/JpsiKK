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


#include "fit.h"

void fit(TH1 * hisKK, TH1 * hisUU)
{
	double Mmin  = hisKK->GetXaxis()->GetXmin();
	double Mmax  = hisKK->GetXaxis()->GetXmax();
  RooRealVar Mrec("Mrec","M_{rec}(#pi^{+}#pi^{-}) - 3097 ", Mmin,Mmax, "MeV");
  RooRealVar  mean( "mean",  "mean",   -0.1 + 0.5*(Mmin+Mmax) ,  Mmin,  Mmax, "MeV") ;
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

	RooMcb2Pdf *mcbPdf=0;
  mcbPdf = new RooMcb2Pdf("mcb2","mcb2",Mrec,mean,sigma,
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
	//if(opt.count("rad"))
	//{
	//	signalPdf = &signalRadPdf;
	//}

	//RooRealVar poly_c1("poly_c1","coefficient of x^1 term",0, -10, 10);
	//RooRealVar poly_c1("poly_c1","coefficient of x^1 term",0);
	RooRealVar bgKK_c1("bgKK_c1","coefficient of Mrec term for KK channel",0, -10, 10);
	//RooPolynomial bgKK("bgKK", "linear function for KK background", Mrec, RooArgList(bgKK_c1));	
	RooPolynomial bgKK("bgKK", "linear function for KK background", Mrec);	
	RooPolynomial bgUU("bgKK", "constant background for UU channel", Mrec);	


	RooRealVar NKK("NKK", "Number of KK events", hisKK->GetEntries(), hisKK->GetEntries()/2.0, hisKK->GetEntries()*100);
	RooRealVar NbgKK("NbgKK", "Number of background events for KK channel", 0, 0,  hisKK->GetEntries()*0.5);
	RooAddPdf  KKPdf("KKPdf", "KK + background p.d.f.", RooArgList(bgKK, *signalPdf), RooArgList(NbgKK,  NKK));

	RooRealVar NUU("NUU", "Number of UU events", hisUU->GetEntries(), hisUU->GetEntries()/2.0, hisUU->GetEntries()*100);
	RooRealVar NbgUU("NbgUU", "Number of background events for UU channel", 0, 0,  hisUU->GetEntries()*0.5);
	RooAddPdf  UUPdf("UUPdf", "UU + background p.d.f.", RooArgList(bgUU, *signalPdf), RooArgList(NbgUU,  NUU));




	RooAbsReal* igx = mcbPdf->createIntegral(Mrec) ;
	cout << "gx_Int[x] = " << igx->getVal() << endl;


	//RooArgSet ntupleVarSet(Mrec);
	//RooDataSet * dataSetKK = new RooDataSet("KK", "KK",  treeKK,  ntupleVarSet);
	//RooDataSet * dataSetUU = new RooDataSet("UU", "UU",  treeUU,  ntupleVarSet);


	RooDataHist * dataKK = new RooDataHist("KK", "KK", Mrec, Import(*hisKK));
	RooDataHist * dataUU = new RooDataHist("UU", "UU", Mrec, Import(*hisUU));

 	RooCategory sample("sample","sample") ;
  sample.defineType("KK") ;
  sample.defineType("UU") ;

  std::map<std::string, RooDataHist*> dataMap;
  dataMap["KK"]=dataKK;
  dataMap["UU"]=dataUU;

  //RooDataHist combData("combData","combined data KK + uu channel",Mrec,Index(sample),Import("KK",*dataKK),Import("UU",*dataUU)) ;
  //
  RooDataHist combData("combData","Combined data", Mrec, sample, dataMap);
	RooDataHist * data = & combData;


  //RooDataSet combData("combData","combined data KK + uu channel",Mrec,Index(sample),Import("KK",*dataSetKK),Import("UU",*dataSetUU)) ;
	//RooDataSet * data = & combData;


  RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;

  // Associate model with the physics state and model_ctl with the control state
  simPdf.addPdf(KKPdf,"KK") ;
  simPdf.addPdf(UUPdf,"UU") ;
  
  RooPlot* xframe3 = Mrec.frame(Title("Show my cb function")) ;
	mcbPdf->plotOn(xframe3);
	TCanvas * c_mcb = new TCanvas;
	xframe3->Draw();

	simPdf.fitTo(combData, Extended(), Strategy(2), Minos());

	RooChi2Var chi2Var("chi2", "chi2", simPdf, *data);
	cout << "chi2 = " << chi2Var.getVal() << endl;
  RooArgSet args(Mrec);
  std::cout << "exp. Events = " << simPdf.expectedEvents(&args) << std::endl;

  RooPlot* frameKK = Mrec.frame(Title("K^{+}K^{-} events")) ;
  RooPlot* frameUU = Mrec.frame(Title("#mu^{+}#mu^{-} events")) ;
  combData.plotOn(frameKK, MarkerSize(0.5),  Cut("sample==sample::KK")) ;
  simPdf.plotOn(frameKK, Slice(sample,"KK"),ProjWData(sample,combData),   LineWidth(1)) ;
  simPdf.plotOn(frameKK, Slice(sample,"KK"),ProjWData(sample,combData), Components(bgKK),LineStyle(kDashed),LineWidth(1)) ;
	frameKK->SetMinimum(0.1);
  combData.plotOn(frameUU, MarkerSize(0.5),  Cut("sample==sample::UU"), "X0") ;
  simPdf.plotOn(frameUU, Slice(sample,"UU"),ProjWData(sample,combData),   LineWidth(1)) ;
  simPdf.plotOn(frameUU, Slice(sample,"UU"),ProjWData(sample,combData), Components(bgUU),LineStyle(kDashed),LineWidth(1)) ;
	frameUU->SetMinimum(0.1);
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

  RooPlot* frame_NKK = NKK.frame(Range(hisKK->GetEntries()*0.75, hisKK->GetEntries()*1.2), Title("NKK -log(L)")) ;
  nll.plotOn(frame_NKK,PrintEvalErrors(-1),ShiftToZero(), EvalErrorValue(nll.getVal()+10),LineColor(kRed)) ;
	frame_NKK->SetMinimum(0);
	frame_NKK->SetMaximum(40);

  RooPlot* frame_NbgKK = NbgKK.frame(Range(0, NbgKK.getValV()*2), Title("NbgKK -log(L)")) ;
  nll.plotOn(frame_NbgKK,PrintEvalErrors(-1),ShiftToZero(), EvalErrorValue(nll.getVal()+10),LineColor(kRed)) ;

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
	cnll->cd(7);
	frame_NKK->Draw();
	cnll->cd(8);
	frame_NbgKK->Draw();


	double theRatio = NKK.getValV()/NUU.getValV();
	double theRelError = sqrt( pow(NKK.getError()/NKK.getValV(), 2)  +  pow(NUU.getError()/NUU.getValV(), 2) );
	double theError = theRatio * theRelError;
	boost::format fmt("%s = (%6.3f ± %-4.3f) * 1e-3 ( ± %.2f%%)");
	//cout << fmt % "Nkk/Nuu = (" << theRatio*1000 << " +/- " << theError*1000  <<  ") * 1e-3 ( +/- " << theRelError*100 << "% )" <<  endl;
	cout << fmt % "Nkk/Nuu" % (theRatio*1000) %  (theError*1000) %  (theRelError*100)  <<  endl;

}
