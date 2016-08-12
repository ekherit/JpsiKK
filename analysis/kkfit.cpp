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

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

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
#include <TTree.h>
#include <TCut.h>

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
  namespace po=boost::program_options;
  po::options_description opt_desc("Allowed options");
  std::string tree_name;
  std::string mctree_name;
  std::string file_name;
  std::string his_name;
  std::string cut;
  std::string hash_list_str;
  std::string suffix;
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
    ("suffix",po::value<std::string>(&suffix),"Suffix")
    ("mc", "Use MonteCarlo information mctopo")
    ;
  po::positional_options_description pos;
  pos.add("input",-1);
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
  std::string bg_from_jpsi = "cfe7a549,4d7eec07,59476175,8a4d7453,4d7eec07,cfe7a549,1ef1d66d,79d29da3,4f68a152,d98afca4,d1310577,84fce654,5d977b97,df0e32d9,57600321,47b3a057,e4ebc641,acaee554,6d47f0a5,9089b143,2555a9e3,7181a332,f5bb9d48,e990b7b9,78a768c7,8fa386fb,5721f99e,6cceb77a,270c3648,6925d619,19f87eb1,167c1bc8";

  std::string nojpsi_bg = "183789a,d0c6c2a9,35c7a381,5f136e0a,708b1663,6df8df1e,59b99e9a,6f6bf2a3,2eb71455,859ea6c1,b1dfe745,a7ee33a3,7025be96,a8b93aa8,5fbdd494,bcdd2654,4200e40e,ccf8f4be,3bfc1a82,cfa021dd,e03859b4,d6d8305c,21dcde60,e6f0a031,e7aec6e2,861f5eda,f46414d5,110dcd9a,3e95b5f3,a5f18375,2ea2390e,54ca98db,7a564e5a,fac211ae,27940cb5,13d54d31,32d836e2,e39d8cd1,552c11f1,d7b558bf,3f291a5a,e01958bf,329ac69e,f922962d,16d4eedf,e8092c85,163b9509,c77e2f3a,4854d053,7c5b6e7b,99116a60";
    
  std::string uu_str = "1397e7e7,8ac60398,aca004d7,cf7de4a7,daa0af44";
  std::string KK_str = "3031ea55,cdffabb3,ecdbe915";
  std::string pipi_str = "67c1bc8,19f87eb1,270c3648,6925d619,5f136e0a,708b1663,a8b93aa8,5fbdd494,cfa021dd,e03859b4,e6f0a031,f46414d5,f5bb9d48,110dcd9a,3e95b5f3,6d47f0a5,9089b143,e4ebc641,7a564e5a,fac211ae,d1310577,d98afca4,27940cb5,13d54d31,79d29da3,4f68a152,1ef1d66d,cfe7a549,4d7eec07,59476175,8a4d7453,4d7eec07,cfe7a549";
  std::string pipiKK_str = "183789a,d0c6c2a9,35c7a381,6df8df1e,59b99e9a,6f6bf2a3,2eb71455,859ea6c1,b1dfe745,bcdd2654,4200e40e,ccf8f4be,3bfc1a82,d6d8305c,21dcde60,a5f18375,5d977b97,df0e32d9,84fce654,32d836e2,e39d8cd1,552c11f1,d7b558bf,3f291a5a,e01958bf,329ac69e,f922962d,16d4eedf,e8092c85,163b9509,c77e2f3a,4854d053,7c5b6e7b,99116a60";
  
  TCut hash_cut = "";

  //std::vector<std::string> hash_list_str;
  std::list<std::string> hash_list;
  if(opt.count("hash"))
  {
    boost::split(hash_list,hash_list_str,boost::is_any_of(", "));
    std::cout << " Initial hash_list: ";
    for ( auto & x : hash_list)
    {
      std::cout << x << ",";
    }
    std::cout << std::endl;
    std::list<std::string> tmp_hash_list;
    for( auto & x : hash_list)
    {
      std::string s=x;
      if ( x == "jpsi" )  s = bg_from_jpsi;
      if ( x == "nojpsi" ) s = nojpsi_bg;
      if ( x == "uu" )  s = uu_str;
      if ( x == "KK" )  s = KK_str;
      if ( x == "pipi" )  s = pipi_str;
      if ( x == "pipiKK" )  s = pipiKK_str;
      std::list<std::string> v;
      boost::split(v, s, boost::is_any_of(", "));
      tmp_hash_list.merge(v);
    }
    hash_list = tmp_hash_list;
    std::cout << "New hash_list: ";
    for ( auto & x : hash_list )
    {
      std::cout << x << ",";
    }
    std::cout << std::endl;
    for(auto & x : hash_list)
    {
      hash_cut = hash_cut || TCut(("hash==0x"+x).c_str());
    }
    if(opt.count("skip"))
    {
      hash_cut = !hash_cut;
    }
  }
  std::cout << "hash cut =" << hash_cut << std::endl;
  //if(argc<2) return 1;
	//TFile file(argv[2]);
  TFile file(file_name.c_str());
  TH1 * his=nullptr;
  TTree * tree=nullptr;
	TApplication theApp("root_app", &argc, argv);

  if(opt.count("his"))
  {
    his = (TH1*)file.Get(his_name.c_str());
  }

  /*
  if(opt.count("tree"))
  {
    tree = (TTree*) file.Get(tree_name.c_str());
    if(opt.count("mctree"))
    {
      TTree *mctree = (TTree*) file.Get(mctree_name.c_str());
      if(mctree) tree->AddFriend(mctree);
    }

    int Nbin = his->GetNbinsX();
    double Mmin = his->GetXaxis()->GetXmin();
    double Mmax  = his->GetXaxis()->GetXmax();

    char varexp[1024];
    sprintf(varexp,"(Mrec-3.097)*1000 >> hMrec(%d,%f,%f)", Nbin, Mmin, Mmax);
    tree->Draw(varexp,hash_cut && cut.c_str(),"goff");
    his = (TH1F*) tree->GetHistogram();
    std::cout << "New entries = " << his->GetEntries() << std::endl;
  }
  */

  if(opt.count("suffix"))
  {
    tree_name   = "event"+suffix;
    his_name    = "hMrec"+suffix;
    mctree_name = "mctopo"+suffix;
    std::cout << tree_name << "  " << his_name << "  " << mctree_name << std::endl;
    his = (TH1*)file.Get(his_name.c_str());
    tree = (TTree*) file.Get(tree_name.c_str());
    if(opt.count("mc"))
    {
      TTree *mctree = (TTree*) file.Get(mctree_name.c_str());
      if(mctree) tree->AddFriend(mctree);
    }
    int Nbin = his->GetNbinsX();
    double Mmin = his->GetXaxis()->GetXmin();
    double Mmax  = his->GetXaxis()->GetXmax();
    char varexp[1024];
    sprintf(varexp,"(Mrec-3.097)*1000 >> hMrec(%d,%f,%f)", Nbin, Mmin, Mmax);
    tree->Draw(varexp,hash_cut && cut.c_str(),"goff");
    his = (TH1F*) tree->GetHistogram();
  }

  double Mmin  = his->GetXaxis()->GetXmin();
  double Mmax  = his->GetXaxis()->GetXmax();
  Long64_t nEntries=his->GetEntries();
  double SigmaInitial=1.4;
  double Scale = 1;
  double M0 = 0;
  
  //if(his)
  //{
  //  nEntries = his->GetEntries();
  //}
  //if(tree)
  //{
  //  nEntries = tree->GetEntries();
  //  Scale = 1000;
  //  M0=3.097;
  //  Mmin = M0 + Mmin/Scale;
  //  Mmax = M0 + Mmax/Scale;
  //  std::cout << nEntries << std::endl;
  //}
  RooRealVar Mrec("Mrec","M_{rec}(#pi^{+}#pi^{-})", Mmin,Mmax, "MeV");
  RooRealVar sigma("sigma","sigma",1.4/Scale,0,10/Scale, "MeV") ;
  RooRealVar staple1("staple1","staple1",M0-15/Scale,   Mmin,Mmax, "MeV") ;
  RooRealVar staple2("staple2","staple2",M0-3/Scale,   Mmin,Mmax, "MeV") ;
  RooRealVar staple3("staple3","staple3",M0 -2/Scale,   Mmin,Mmax, "MeV") ;
  RooRealVar staple4("staple4","staple4",M0+0,   M0-2/Scale,M0+2/Scale, "MeV") ;
  RooRealVar staple5("staple5","staple5",M0+2/Scale,   Mmin,Mmax, "MeV") ;
  RooRealVar staple6("staple6","staple6",M0+4/Scale,   Mmin,Mmax, "MeV") ;
  RooRealVar staple7("staple7","staple7",M0+15/Scale,   Mmin,Mmax, "MeV") ;
  RooRealVar n1("n1","n1", 3,  1,100) ;
  RooRealVar n2("n2","n2", 2,  1,100) ;



	//std::vector<RooRealVar> staple(7);
	//for(int i=0;i<staple.size();i++)
	//{
	//	char buf[128];
	//	sprintf(buf, "staple-%d",i );
	//	staple[i] = RooRealVar(buf, buf, (i-3.5)*5,  -40, +40, "MeV");
	//}
	//std::vector<RooRealVar> N(2);
	//N[0] = RooRealVar("nl", "Left power",  2,  0.1, 100);
	//N[1] = RooRealVar("nr", "Right power", 2,  0.1, 100);

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


	RooRealVar mcb_yield("mcb_yield", "MCB yield", nEntries, 0, nEntries);
	RooRealVar sg_yield("sg_yield", "small gauss yield", 0, 0, nEntries);
	RooArgList sig_shapes;
	RooArgList sig_yields;
	sig_shapes.add(small_gaus);      sig_yields.add(sg_yield);
	sig_shapes.add(mcb);             sig_yields.add(mcb_yield);
	RooAddPdf complexSignal("complexSignal", "signal with small gaus",  sig_shapes,  sig_yields);

	RooRealVar poly_c1("poly_c1","coefficient of x^1 term",0, -10, 10);
	//RooRealVar poly_c1("poly_c1","coefficient of x^1 term",0);
	RooPolynomial bkgd_poly("bkgd_poly", "linear function for background", Mrec, RooArgList(poly_c1));	


	RooRealVar Nsig("Nsig", "Number of signal events", nEntries, 0, nEntries*100);
  RooRealVar Nbg("Nbg", "Number of background events", nbg, 0,  nEntries);
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
  //RooAbsData * data;
  RooDataHist * data = new RooDataHist("dh", "dh", Mrec, Import(*his));
  //if(opt.count("his"))
  //{
  //  data =   }
  //if(opt.count("tree"))
  //{
  //  data = new RooDataSet("dh","dh", Mrec, Import(*tree), Cut(cut.c_str()));
  //  //data = new RooDataSet("dh","dh", Mrec, Import(*tree));
  //  //new RooDataSet("dh","dh", Mrec, Import(*tree));
  //}
  
  RooPlot* xframe3 = Mrec.frame(Title("Show my cb function")) ;
	mcb.plotOn(xframe3);
	TCanvas * c_mcb = new TCanvas;
	xframe3->Draw();

	//theApp.Run();
	//totalPdf.fitTo(*data,  Extended(), FitOptions("qmh"));
	totalPdf.fitTo(*data,  Extended(), Strategy(2));

	//RooChi2Var chi2Var("chi2", "chi2", totalPdf, *data);
	//cout << "chi2/ndf = " << chi2Var.getVal() << endl;

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

  RooPlot* frame_sigma = sigma.frame(Range(0.2/Scale, 5.0/Scale), Title("-log(L) scan vs sigma")) ;
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

  RooPlot* frame_Nsig = Nsig.frame(Range(0, nEntries*1.2), Title("Nsig -log(L)")) ;
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

  /*
	for(int i=0;i<his->GetNbinsX();i++)
	{
		double x  = his->GetBinCenter(i);
		double N  = his->GetBinContent(i);
		double Nth;
		//double Nth = mcb.evaluate(); 
		//cout <<  i << " " << x << "  " << N <<   " " << Nth << endl;
	}
  */
	
	theApp.Run();
}
