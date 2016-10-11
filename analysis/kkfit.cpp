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

#include <boost/format.hpp>
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
#include <RooSimultaneous.h>

#include "RooMcbPdf.h"

#include "hashlist.h"

using namespace RooFit;
using namespace std;

struct RooFitItem_t
{
  TH1           * his;
  RooDataHist   * data;
  RooRealVar    * Nsig;
  RooRealVar    * Nbg;
  RooRealVar    * bg_c1;
  RooPolynomial * bgPdf;
  RooAddPdf     * addPdf;
  RooPlot       * frame;
  std::string fName;

  RooFitItem_t(TH1 * h, RooMcbPdf* mcbPdf, RooRealVar & Mrec)
  {
    RooFitItem_t(h->GetName(), h, mcbPdf, Mrec);
  }

  RooFitItem_t(std::string name, TH1 * h, RooMcbPdf* mcbPdf, RooRealVar & Mrec)
  {
    fName = name;
    his = h;
    data = new RooDataHist(h->GetName(),  h->GetTitle(), Mrec, Import(*h));
    Nsig = new RooRealVar(
        (std::string("Nsig")+h->GetName()).c_str(), 
        (std::string("Number of signal events for  ") + h->GetTitle()).c_str(), 
        h->GetEntries(), 
        0, 
        h->GetEntries()*100
        );
    Nbg = new RooRealVar(
        (std::string("Nbg")+h->GetName()).c_str(), 
        (std::string("Number of substrate events for  ") + h->GetTitle()).c_str(), 
        h->GetEntries()*0.001, 
        0, 
        h->GetEntries()*100
        );
    std::string bg_name = std::string("bg")+h->GetName();
    bg_c1 = new RooRealVar((bg_name+"c1").c_str(),(std::string("Linear coefficient for background for ") + h->GetName()).c_str(),0, -10, 10);
    std::string bg_title = std::string("Background Pdf for ") + h->GetName();
    bgPdf = new RooPolynomial(bg_name.c_str(),bg_title.c_str(), Mrec, *bg_c1);
    addPdf = new RooAddPdf(name.c_str(),h->GetTitle(), RooArgList(*bgPdf,*mcbPdf), RooArgList(*Nbg, *Nsig));
    frame = Mrec.frame(Title(his->GetTitle())) ;
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
  RooFitItem_t f1(h1->GetName(),h1, mcbPdf, Mrec);
  std::clog << "Before FitItem_t f2" << std::endl;
  RooFitItem_t f2(h2->GetName(),h2, mcbPdf, Mrec);

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
  RooRealVar staple1("staple1","staple1",  15,   Mmin,Mmax, "MeV") ;
  RooRealVar staple2("staple2","staple2",  3,   Mmin,Mmax, "MeV") ;
  RooRealVar staple3("staple3","staple3",  2,   Mmin,Mmax, "MeV") ;
  RooRealVar staple4("staple4","staple4",  0,   2, 2, "MeV") ;
  RooRealVar staple5("staple5","staple5",  2,   Mmin,Mmax, "MeV") ;
  RooRealVar staple6("staple6","staple6",  4,   Mmin,Mmax, "MeV") ;
  RooRealVar staple7("staple7","staple7",  15,   Mmin,Mmax, "MeV") ;
  RooRealVar n1("n1","n1", 3,  1,100) ;
  RooRealVar n2("n2","n2", 2,  1,100) ;

  //this function describes signal
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
  std::cout << "Before sample" << std::endl;
  RooCategory sample("sample","sample");
  std::vector< RooFitItem_t *> fi_lst;
  std::map<std::string, RooDataHist*> dataMap;
  std::cout << "Before his loop" << std::endl;
  for(auto h : his_list)
  {
    std::cout << "Creating fitItem: " << h->GetName();
    RooFitItem_t * f = new RooFitItem_t(h, mcbPdf, Mrec);
    fi_lst.push_back(f);
    auto name = f->name();
    sample.defineType(name.c_str());
    std::cout << name << std::endl;
  }
  std::cout << "Before simPdf" << std::endl;
  RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;
  std::cout << "Before fi_lst loop" << std::endl;
  for(auto f : fi_lst)
  {
    std::cout << "Inside fi_lst loop" << std::endl;
    std::cout << f->his->GetName() << std::endl;
    auto name = f->his->GetName();
    std::cout << "Adding to simPdf: " << f->his->GetName();
    simPdf.addPdf(*(f->addPdf),name) ;
    dataMap[name] = f->data;
  }
  std::cout << "Before RooDataHist" << std::endl;
  RooDataHist * data = new RooDataHist("data","Combined data", Mrec, sample, dataMap);
  std::clog << "Before fitTo " << std::endl;
	simPdf.fitTo(*data, Extended(), Strategy(2), Minos());
  auto frame = Mrec.frame(Title("title"));
  for(auto f : fi_lst )
  {
    data->plotOn(frame, MarkerSize(0.5),  Cut((std::string("sample==sample::")+f->his->GetName()).c_str())) ;
    simPdf.plotOn(frame, Slice(sample,f->his->GetName()),ProjWData(sample,*data), LineWidth(1)) ;
    simPdf.plotOn(frame, Slice(sample,f->his->GetName()),ProjWData(sample,*data), Components(*f->bgPdf),LineStyle(kDashed),LineWidth(1)) ;
  }
  new TCanvas;
  frame->Draw();
}


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
    ;
  po::positional_options_description pos;
  pos.add("suffix", 1);
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

  std::cout << " " << suffix << " " << file_name << std::endl;
  //exit(1);
  
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
    h->SetTitle(title.c_str());
    return h;
  };

  if(opt.count("combine"))
  {
    std::vector<std::string> suffix_list;
    boost::split(suffix_list,combine_str,boost::is_any_of(", "));
    std::list<TH1*> his_lst;
    for(auto str : suffix_list)
    {
      his_lst.push_back(create_histogram(str));
    }
    //auto h1 = create_histogram(suffix_list[0]);
    //auto h2 = create_histogram(suffix_list[1]);
    //combfit(h1,h2);
    combfit(his_lst);
    theApp.Run();
  }

  if(opt.count("his"))
  {
    his = (TH1*)file.Get(his_name.c_str());
  }

  const double Scale = 1;
  const double M0 = 0;


  if(opt.count("suffix"))
  {
    his = create_histogram(suffix);
  }

  double Mmin  = his->GetXaxis()->GetXmin();
  double Mmax  = his->GetXaxis()->GetXmax();
  Long64_t nEntries=his->GetEntries();
  double SigmaInitial=1.4;
  
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
  RooDataHist * data = new RooDataHist("dh", "dh", Mrec, Import(*his));
	totalPdf.fitTo(*data,  Extended(), Strategy(2));

  RooPlot* xframe2 = Mrec.frame(Title("Fit by Modified CrystalBall")) ;
  data->plotOn(xframe2, MarkerSize(0.5)) ;
  totalPdf.plotOn(xframe2,  LineWidth(1)) ;
	totalPdf.plotOn(xframe2,Components(bkgd_poly),LineStyle(kDashed),  LineWidth(1)) ;
  

  sigma.Print() ;
	Nsig.Print();
	Nbg.Print();


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

	TCanvas * cnll =new TCanvas("cnll", "- Log likelihood");
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

  // Draw all frames on a canvas
  TCanvas* c = new TCanvas("mcb","Test for Modified CrystalBall Fit") ;
	c->SetLogy();
  c->cd(2) ; 
	gPad->SetLeftMargin(0.15) ; 
	xframe2->GetYaxis()->SetTitleOffset(1.6) ; 
	xframe2->Draw() ;
	
	theApp.Run();
}
