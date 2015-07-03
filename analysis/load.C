#include <regex>
#include <math.h>

#include <list>
#include <TRandom.h>

//#include "mctop/McTopo.h"
//#include "mctop/MyEvent.h"
#include "CrystalBall.h"

#include "analize.h"
#include <TSystemDirectory.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCut.h>
#include <TTree.h>
#include <TH2F.h>


#include "libFit.h"

list<string> list_files(const char *dirname=".", const char *ext=".root")
{

	list<string> flist;
	cout << dirname << endl;
	TSystemDirectory dir(dirname, dirname);
	TList *files = dir.GetListOfFiles();
	if (files) 
	{
		TSystemFile *file;
		TString fname;
		TIter next(files);
		while ((file=(TSystemFile*)next()))
		{
			fname = file->GetName();
			if (!file->IsDirectory() && fname.EndsWith(ext))
			{
				flist.push_back(fname.Data());
			}
		}
	}
	flist.sort();
	return flist;
}
//TChain * load_tree(const char *dirname=".", const char *ext=".root")
TChain * load_tree(string dirname=".", string ext=".root",const char * treename="event")
{
	TChain * event = new TChain(treename, treename);
	TString pwd(gSystem->pwd());
	TSystemDirectory dir(dirname.c_str(), dirname.c_str());
	TList *files = dir.GetListOfFiles();
	gSystem->cd(pwd.Data());
	if (files) 
	{
		TSystemFile *file;
		TString fname;
		TIter next(files);
		while ((file=(TSystemFile*)next()))
		{
			fname = file->GetName();
			if (!file->IsDirectory() && fname.EndsWith(ext.c_str()))
			{
				string f = dirname + "/" + fname.Data();
				event->AddFile(f.c_str());
				//flist.push_back(fname.Data());
			}
		}
	}
	return event;
}

TChain * load_tree2(string dirname=".", string ext=".root",const char * treename="event")
{
	TChain * event = new TChain(treename, treename);
	TString pwd(gSystem->pwd());
	TSystemDirectory dir(dirname.c_str(), dirname.c_str());
	TList *files = dir.GetListOfFiles();
	gSystem->cd(pwd.Data());
	if (files) 
	{
		TSystemFile *file;
		TString fname;
		TIter next(files);
		while ((file=(TSystemFile*)next()))
		{
			fname = file->GetName();
			if (!file->IsDirectory() && fname.EndsWith(ext.c_str()))
			{
				string f = dirname + "/" + fname.Data();
				event->AddFile(f.c_str());
				//flist.push_back(fname.Data());
			}
		}
	}
	return event;
}


void Analize(const char * dir, const char * file="analize_result.root", Long64_t N=0)
{
	gSystem->CompileMacro("analize.C","kO","","/tmp");
	TChain * event = load_tree(dir,".root","event");
	TChain * mctopo = load_tree(dir,".root","mctopo");
  cout << "Total number of events: " << event->GetEntriesFast() << std::endl;
  event->AddFriend(mctopo);
  analize an;
  //an.output_file_name = file;
  if(N==0) event->Process(&an,file);
  else event->Process(&an,file,N);
	cout << "After process" << endl;
}


void fit(TH1F * his)
{
  vector<double> res(3);
  res = Fit(his);
  //cout << "Number of selected events" << res[0] << " " << -res[1] << " +" << res[2] << endl;
}

void fit(const char * rootfile, const char * hisname)
{
  TFile * f = new TFile(rootfile);
  f->ls();
  TH1F * h = (TH1F*)f->Get(hisname);
  h->Draw("E");
  fit(h);
}


void draw_by_topo(TTree *t, vector<double> pdgid)
{
  for(auto v: pdgid) cout << v << endl;
}


void show_cuts(void)
{
  TChain * mc_event   =  load_tree("../data/mc09",".root","event");
  TChain * mc_mdc     =  load_tree("../data/mc09",".root","mdc");
  TChain * mc_emc     =  load_tree("../data/mc09",".root","emc");

  TChain * data_event =  load_tree("../data/data09",".root","event");
  TChain * data_mdc   =  load_tree("../data/data09",".root","mdc");
  TChain * data_emc   =  load_tree("../data/data09",".root","emc");
  TCanvas * cp = new TCanvas("pion_cuts_canvas","Pion momentum cuts");
  TCut main_cut = "kin_chi2<40 && pid_chi2<20 && 3.08 < Mrec && Mrec<3.114";
  mc_event->SetLineColor(kRed);
  data_event->SetLineColor(kRed);
  data_event->SetMarkerStyle(8);
  data_event->SetMarkerColor(kRed);
  mc_event->Draw("p[0]>>hmcKKp(100,0,0.45)","KK" && main_cut ,"NORM");
	TH1* hmcKKp = (TH1*)gDirectory->Get("hmcKKp");
  hmcKKp->SetTitle("Pion momentum");
  mc_event->GetHistogram()->GetXaxis()->SetTitle("p, GeV/c");
  data_event->Draw("p[0]>>hdataKKp(100,0,0.45)","KK" && main_cut ,"ENORMSAME");
	TH1 * hdataKKp = (TH1*)gDirectory->Get("hdataKKp");
  mc_event->SetLineColor(kBlue);
  data_event->SetLineColor(kBlue);
  data_event->SetMarkerStyle(8);
  data_event->SetMarkerColor(kBlue);
  mc_event->Draw("p[0]>>hmcUUp(100,0,0.45)","uu" && main_cut ,"NORMSAM");
	TH1 * hmcUUp = (TH1*)gDirectory->Get("hmcUUp");
  mc_event->GetHistogram()->GetXaxis()->SetTitle("p, GeV/c");
  data_event->Draw("p[0]>>hdataUUp(100,0,0.45)","uu" && main_cut ,"ENORMSAME");
	TH1 * hdataUUp = (TH1*)gDirectory->Get("hdataUUp");
  TLegend *l = new TLegend(0.8,0.8,1.0,1.0);
  l->AddEntry(hmcKKp,"MC2009 K^{+}K^{-}","lp");
  l->AddEntry(hdataKKp,"data2009 K^{+}K^{-}","l");
  l->AddEntry(hmcUUp,"MC2009 #mu^{+}#mu^{-}","lp");
  l->AddEntry(hdataUUp,"data2009 #mu^{+}#mu^{-}","lp");
  l->Draw();
  TCanvas * Mrec_canvas = new TCanvas("canvas_pion_Mrec","Pion recoil invariant mass");
  //main_cut = "kin_chi2<40 && pid_chi2<20";
  main_cut = "";
  mc_event->SetLineColor(kRed);
  data_event->SetLineColor(kRed);
  data_event->SetMarkerStyle(8);
  data_event->SetMarkerColor(kRed);
  mc_event->Draw("Mrec>>hmcKKMrec(100,3.0,3.2)","KK" && main_cut ,"NORM");
	TH1 * hmcKKMrec = (TH1*)gDirectory->Get("hmcKKMrec");
  hmcKKp->SetTitle("Pion momentum");
  mc_event->GetHistogram()->GetXaxis()->SetTitle("p, GeV/c");
  data_event->Draw("Mrec>>hdataKKMrec(100,3.0,3.2)","KK" && main_cut ,"ENORMSAME");
	TH1 * hdataKKMrec = (TH1*)gDirectory->Get("hdataKKMrec");
  mc_event->SetLineColor(kBlue);
  data_event->SetLineColor(kBlue);
  data_event->SetMarkerStyle(8);
  data_event->SetMarkerColor(kBlue);
  mc_event->Draw("Mrec>>hmcUUMrec(100,3.0,3.2)","uu" && main_cut ,"NORMSAME");
	TH1 * hmcUUMrec = (TH1*)gDirectory->Get("hmcUUMrec");
  mc_event->GetHistogram()->GetXaxis()->SetTitle("p, GeV/c");
  data_event->Draw("Mrec>>hdataUUMrec(100,3.0,3.2)","uu" && main_cut ,"ENORMSAME");
	TH1 * hdataUUMrec = (TH1*)gDirectory->Get("hdataUUMrec");
  TLegend *lMrec = new TLegend(0.8,0.8,1.0,1.0);
  lMrec->AddEntry(hmcKKMrec,"MC2009 K^{+}K^{-}","lp");
  lMrec->AddEntry(hdataKKMrec,"data2009 K^{+}K^{-}","lp");
  lMrec->AddEntry(hmcUUMrec,"MC2009 #mu^{+}#mu^{-}","lp");
  lMrec->AddEntry(hdataUUMrec,"data2009 #mu^{+}#mu^{-}","lp");
  lMrec->Draw();

}

void show_chi2(void)
{
  TChain * mc_event   =  load_tree("../data/mc09",".root","event");
  TChain * mc_mdc     =  load_tree("../data/mc09",".root","mdc");
  TChain * mc_emc     =  load_tree("../data/mc09",".root","emc");

  TChain * data_event =  load_tree("../data/data09",".root","event");
  TChain * data_mdc   =  load_tree("../data/data09",".root","mdc");
  TChain * data_emc   =  load_tree("../data/data09",".root","emc");
	TCanvas * c = new TCanvas;
  TCut main_cut = "abs(Mrec-3.097)<0.45";
	//data_event->SetLineColor(kBlack);
	//data_event->Draw("pid_chi2>>hdpidKK(100, 0, 30)", main_cut && "KK" && "pid_chi2<30 && kin_chi2<40", "NORME");
	//mc_event->SetLineColor(kBlue);
	//mc_event->Draw("pid_chi2>>hmcpidKK(100, 0, 30)",  main_cut && "KK"  && "pid_chi2<30 && kin_chi2<40", "NORMSAME");
	//data_event->SetLineColor(kMagenta);
	//data_event->Draw("pid_chi2>>hdpiduu(100, 0, 30)", main_cut && "uu" && "pid_chi2<30 && kin_chi2<40", "NORMESAME");
	//mc_event->SetLineColor(kRed);
	//mc_event->Draw("pid_chi2>>hmcpiduu(100, 0, 30)",  main_cut && "uu" && "pid_chi2<30 && kin_chi2<40", "NORMSAME");
  //TLegend *l = new TLegend(0.8,0.8,1.0,1.0);
	//TH1 * hdpidKK = (TH1*)gDirectory->Get("hdpidKK");
	//TH1 * hmcpidKK = (TH1*)gDirectory->Get("hmcpidKK");
	//TH1 * hdpiduu = (TH1*)gDirectory->Get("hdpiduu");
	//TH1 * hmcpiduu = (TH1*)gDirectory->Get("hmcpiduu");
  //l->AddEntry(hdpidKK,"data2009 K^{+}K^{-}","lp");
  //l->AddEntry(hmcpidKK,"mc2009 K^{+}K^{-}","lp");
  //l->AddEntry(hdpiduu,"data2009 #mu^{+}#mu^{-}","lp");
  //l->AddEntry(hmcpiduu,"mc2009 #mu^{+}#mu^{-}","lp");

	new TCanvas;
	data_event->SetLineColor(kBlack);
	data_event->Draw("kin_chi2>>hdkinKK(100, 0, 50)", main_cut && "KK" && "kin_chi2<50  && pid_chi2<20", "NORME");
	TH1 * hdkinKK = (TH1*)gDirectory->Get("hdkinKK");
	mc_event->SetLineColor(kBlue);
	mc_event->Draw("kin_chi2>>hmckinKK(100, 0, 50)",  main_cut && "KK"  && "kin_chi2<50 && pid_chi2<20", "NORMSAME");
	TH1 * hmckinKK = (TH1*)gDirectory->Get("hmckinKK");
	data_event->SetLineColor(kMagenta);
	data_event->Draw("kin_chi2>>hdkinuu(100, 0, 50)", main_cut && "uu" && "kin_chi2<50  && pid_chi2<20", "NORMESAME");
	TH1 * hdkinuu = (TH1*)gDirectory->Get("hdkinuu");
	mc_event->SetLineColor(kRed);
	mc_event->Draw("kin_chi2>>hmckinuu(100, 0, 50)",  main_cut && "uu" && "kin_chi2<50  && pid_chi2<20", "NORMSAME");
	TH1 * hmckinuu = (TH1*)gDirectory->Get("hmckinuu");
  TLegend *lkin = new TLegend(0.8,0.8,1.0,1.0);
  lkin->AddEntry(hdkinKK,"data2009 K^{+}K^{-}","lp");
  lkin->AddEntry(hmckinKK,"mc2009 K^{+}K^{-}","lp");
  lkin->AddEntry(hdkinuu,"data2009 #mu^{+}#mu^{-}","lp");
  lkin->AddEntry(hmckinuu,"mc2009 #mu^{+}#mu^{-}","lp");
  lkin->Draw();
  lkin->Draw();
}


void show_radcor(const char * filename, const char * topo_name="mctopoKK", const char *event_name="eventKK", const char * cut="")
{
  TChain * mctopKK = new TChain(topo_name,topo_name);
  TChain * eventKK = new TChain(event_name,event_name);
  eventKK->AddFriend(mctopKK);
  mctopKK->AddFile(filename);
  eventKK->AddFile(filename);

  TCut main_cut(cut);

  char buf[1024];
  sprintf(buf,"%s.indexmc",topo_name);
  eventKK->SetAlias("indexmc",buf);
  sprintf(buf,"%s.pdgid",topo_name);
  eventKK->SetAlias("pdgid",buf);

  new TCanvas;
  eventKK->Draw("indexmc");
  eventKK->Draw("Mrec>>h0","mctopoKK.indexmc==6");

  gStyle->SetOptStat(kFALSE);
  TCanvas * c =  new TCanvas;
  c->SetLogy();
  eventKK->SetLineColor(kMagenta); 
  eventKK->Draw("Mrec>>h0",main_cut && "indexmc==6");
	TH1 * h0 = (TH1*)gDirectory->Get("h0");
  eventKK->GetHistogram()->GetXaxis()->SetTitle("M_{rec}(#pi^{+}#pi^{-}), GeV");
  eventKK->GetHistogram()->SetTitle("Pion recoil mass");
  eventKK->SetLineColor(kBlue); 
  eventKK->Draw("Mrec>>hpi1",main_cut && "indexmc==7&& pdgid[4]==-22","same");
	TH1 * hpi1 = (TH1*)gDirectory->Get("hpi1");
  eventKK->SetLineColor(kGreen); 
  eventKK->Draw("Mrec>>hK1",main_cut && "indexmc==7&& pdgid[6]==-22","same");
	TH1 * hK1 = (TH1*)gDirectory->Get("hK1");
  eventKK->SetLineColor(kRed); 
  eventKK->Draw("Mrec>>hK2",main_cut && "indexmc==8","same");
	TH1 * hK2 = (TH1*)gDirectory->Get("hK2");
  eventKK->SetLineColor(kBlack); 
  eventKK->Draw("Mrec>>htot",main_cut,"Esame");
	TH1 * htot = (TH1*)gDirectory->Get("htot");
  TLegend * l = new TLegend(0.8,0.8,1.0,1.0);
  l->AddEntry(h0,"#psi(2S)#rightarrow#pi^{+}#pi^{-} J/#psi #rightarrow K^{+}K^{-}","l");
  l->AddEntry(hpi1,"#psi(2S)#rightarrow#pi^{+}#pi^{-}#gamma J/#psi #rightarrow K^{+}K^{-}","l");
  l->AddEntry(hK1,"#psi(2S)#rightarrow#pi^{+}#pi^{-} J/#psi #rightarrow K^{+}K^{-} #gamma","l");
  l->AddEntry(hK2,"#psi(2S)#rightarrow#pi^{+}#pi^{-} J/#psi #rightarrow K^{+}K^{-} #gamma #gamma","l");
  l->AddEntry(htot,"#psi(2S)#rightarrow#pi^{+}#pi^{-} (#gamma)J/#psi #rightarrow K^{+}K^{-} (#gamma)","l");
  //eventKK->SetLineColor(kPink); 
  //eventKK->Draw("Minv>>hinv","","same");
  l->Draw();
}


void interf(Long64_t N=0)
{
  TChain * KK = load_tree("../../mckk09");
  TChain * bg = load_tree("../../mcK1K2009");
  TH2F * htmp = new TH2F("h","h",200,0.6,1.8,200,0.8,2);
  TCanvas * cKK = new TCanvas;
  gStyle->SetOptStat(false);
  htmp->Draw();
  TCut cut = "";
  cut = "KK && kin_chi2<40 && pid_chi2<20 && abs(Mrec-3.097)<0.006";
  KK->Draw("M012:M12>>hkk(100,0.6,1.8,100,0.8,2)",cut,"col2NORM");
	TH1 * hkk = (TH1*)gDirectory->Get("hkk");
  hkk->GetXaxis()->SetTitle("M_{#pi^{-}K^{+}}, GeV");
  hkk->GetYaxis()->SetTitle("M_{#pi^{+}#pi^{-}K^{+}}, GeV");
  TCanvas * cbg = new TCanvas;
  bg->Draw("M012:M12>>hbg(100,0.6,1.8,100,0.8,2)",cut,"col2NORM");
	TH1 * hbg = (TH1*)gDirectory->Get("hbg");
  hbg->GetXaxis()->SetTitle("M_{#pi^{-}K^{+}}, GeV");
  hbg->GetYaxis()->SetTitle("M_{#pi^{+}#pi^{-}K^{+}}, GeV");
  double sum=0;
  double sumA=0;
  double sumB=0;
  for(int i=0;i<200;i++)
    for(int j=0;j<200;j++)
    {
      double a = sqrt(hkk->GetBinContent(i,j));
      double b = sqrt(hbg->GetBinContent(i,j));
      sum+=a*b;
      sumA+=a*a;
      sumB+=b*b;
    }
  cout << sum << endl;
  cout << "Norm A = " << sumA << endl;
  cout << "Norm B = " << sumB << endl;
  new TCanvas;
  KK->SetMarkerColor(kBlack);
  KK->Draw("M012:M12",cut,"cont");
  //bg->SetMarkerColor(kRed);
  //bg->Draw("M012:M12",cut,"contsame");
}

void show_Kmu_momentum(void)
{
  TChain * mc = load_tree("../../mckkuu09");
  TChain * data = load_tree("../../2009/");
  TH2F * h = new TH2F("h","h",200,1.2,1.8,200,0,0.03);
  h->Draw();
  h->GetXaxis()->SetTitle("p, GeV/c");
  TCut cut = "pid_chi2<20 && kin_chi2<40 && abs(Mrec-3.097)<0.1";
  cut="";
  mc->SetLineColor(kBlue);
  mc->Draw("p[2]>>hmcKK" ,"KK"&& cut,"NORMSAME");
	TH1 * hmcKK = (TH1*)gDirectory->Get("hmcKK");
  mc->SetLineColor(kRed);
  mc->Draw("p[2]>>hmcuu" ,"uu"&& cut,"NORMSAME");
	TH1 * hmcuu = (TH1*)gDirectory->Get("hmcuu");
  data->SetLineColor(kBlack);
  data->SetMarkerColor(kBlack);
  data->Draw("p[2]>>hdKK","KK"&& cut,"ENORMSAME");
	TH1 * hdKK = (TH1*)gDirectory->Get("hdKK");
  data->SetLineColor(kGreen);
  data->SetMarkerColor(kGreen);
  data->Draw("p[2]>>hduu","uu"&& cut,"ENORMSAME");
	TH1 * hduu = (TH1*)gDirectory->Get("hduu");
  TLegend * l = new TLegend(0.8,0.8,1.0,1.0);
  l->AddEntry(hmcKK,"MC09 KK","lp");
  l->AddEntry(hmcuu,"MC09 #mu#mu","lp");
  l->AddEntry(hdKK,"data09 KK","lp");
  l->AddEntry(hduu,"data09 #mu#mu","lp");
  l->Draw();
}

void show_Ep(void)
{
  TChain * mc = load_tree("../../mckkuu09");
  TChain * data = load_tree("../../2009/");
  TChain * mcMdc = load_tree("../../mckkuu09",".root","mdc");
  TChain * dataMdc = load_tree("../../2009/",".root","mdc");
  mc->AddFriend(mcMdc);
  data->AddFriend(dataMdc);
  TH2F * h = new TH2F("h","h",200,1.2,1.8,200,0,0.03);
  //h->Draw();
  //h->GetXaxis()->SetTitle("p, GeV/c");
  TCut cut = "pid_chi2<20 && kin_chi2<40 && abs(Mrec-3.097)<0.1";
  mc->SetLineColor(kRed);
  mc->Draw("mdc.E[2]/mdc.p[2]>>hmcuu(200,0,0.8)" ,"uu"&& cut,"NORM");
	TH1 * hmcuu = (TH1*)gDirectory->Get("hmcuu");
  mc->SetLineColor(kBlue);
  mc->Draw("mdc.E[2]/mdc.p[2]>>hmcKK(100,0,0.8)" ,"KK"&& cut,"NORMSAME");
	TH1 * hmcKK = (TH1*)gDirectory->Get("hmcKK");
  data->SetLineColor(kGreen);
  data->SetMarkerColor(kGreen);
  data->Draw("mdc.E[2]/mdc.p[2]>>hduu(200,0,0.8)","uu"&& cut,"ENORMSAME");
	TH1 * hduu = (TH1*)gDirectory->Get("hduu");
  data->SetLineColor(kBlack);
  data->SetMarkerColor(kBlack);
  data->Draw("mdc.E[2]/mdc.p[2]>>hdKK(100,0,0.8)","KK"&& cut,"ENORMSAME");
	TH1 * hdKK = (TH1*)gDirectory->Get("hdKK");
  TLegend * l = new TLegend(0.8,0.8,1.0,1.0);
  l->AddEntry(hmcKK,"MC09 KK","lp");
  l->AddEntry(hmcuu,"MC09 #mu#mu","lp");
  l->AddEntry(hdKK,"data09 KK","lp");
  l->AddEntry(hduu,"data09 #mu#mu","lp");
  l->Draw();
}


void load(void)
{
	//gROOT->Reset();
  //gSystem->Load("mctop/libMyEvent.so");
  ////std::cout << gSystem->GetMakeSharedLib() << endl;
  ////gSystem->SetMakeSharedLib("-std=c++11");
  //gSystem->AddIncludePath("-I$HOME/work -I./mctop");
	//gSystem->CompileMacro("CrystalBall.cpp","kO","","/tmp");
	////gSystem->CompileMacro("mctopo.C","kO","","/tmp");
	//gSystem->CompileMacro("analize.C","kO","","/tmp");
}

//add n to every bin of histogram
void add_uniform_events(TH1 * h,  double n)
{
	int nbins = h->GetNbinsX();
	for(int i=0;i<nbins;i++)
	{
		h->SetBinContent(i,  h->GetBinContent(i) + n/nbins);
	}
}



//#include <RooFit.h>
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

void roofit(TTree *tree)
{
	gSystem->Load("libRooFit") ;
	//using namespace RooFit ;

  // S e t u p   m o d e l 
  // ---------------------

  // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
  RooRealVar Mrec("Mrec","Mrec",3.06,3.14) ;
  RooRealVar mean("mean","mean of gaussian",3.097,3.0,3.2) ;
  RooRealVar sigma("sigma","width of gaussian",0.001,0.0001,0.005) ;

  RooRealVar lambda("lambda","asymetry",0.001, -1, 1);
  RooRealVar rho_left("rho_left","left tail",0.001, 0, 100);
  RooRealVar rho_right("rho_right","right tail",0.001, 0, 100);

  // Build gaussian p.d.f in terms of x,mean and sigma
  //RooGaussian gauss("gauss","gaussian PDF",Mrec,mean,sigma) ;  
  RooBukinPdf gauss("bukin","bukin PDF",Mrec,mean,sigma, lambda, rho_left, rho_right) ;  

	RooRealVar poly_c1("poly_c1","coefficient of x^1 term",0,-10,10);
	RooPolynomial bkgd_poly("bkgd_poly", "linear function for background", Mrec, RooArgList(poly_c1));	

  // Construct plot frame in 'x'
  //RooPlot* xframe = Mrec.frame(Title("Gaussian p.d.f.")) ;

	RooRealVar peak_yield("peak_yield", "yield signal peak", 3000, 0, 1000000);
	RooRealVar bkgd_yield("bkgd_yield", "yield of background", 500, 0, 1000000);
	RooArgList shapes;
	RooArgList yields;
	shapes.add(bkgd_poly);      yields.add(bkgd_yield);
	shapes.add(gauss);  yields.add(peak_yield);
	RooAddPdf  totalPdf("totalPdf", "sum of signal and background PDF's", shapes, yields);


  // P l o t   m o d e l   a n d   c h a n g e   p a r a m e t e r   v a l u e s
  // ---------------------------------------------------------------------------

  // Plot gauss in frame (i.e. in x) 
  //gauss.plotOn(xframe) ;

  // Change the value of sigma to 3
  //sigma.setVal(3) ;

  // Plot gauss in frame (i.e. in x) and draw frame on canvas
  //gauss.plotOn(xframe,LineColor(kRed)) ;
  

  // G e n e r a t e   e v e n t s 
  // -----------------------------

  // Generate a dataset of 1000 events in x from gauss
  //RooDataSet* data = gauss.generate(x,10000) ;  
	RooArgSet ntupleVarSet(Mrec);
	RooDataSet * data = new RooDataSet("tree", "tree",  tree,  ntupleVarSet);
  
  // Make a second plot frame in x and draw both the 
  // data and the p.d.f in the frame
  RooPlot* xframe2 = Mrec.frame(Title("Gaussian p.d.f. with data")) ;
  data->plotOn(xframe2) ;
  totalPdf.plotOn(xframe2) ;
	totalPdf.plotOn(xframe2,Components(bkgd_poly),LineStyle(kDashed)) ;
  

  // F i t   m o d e l   t o   d a t a
  // -----------------------------

  // Fit pdf to data
  //gauss.fitTo(*data) ;
	totalPdf.fitTo(*data, Extended());

  // Print values of mean and sigma (that now reflect fitted values and errors)
  mean.Print() ;
  sigma.Print() ;

  // Draw all frames on a canvas
  TCanvas* c = new TCanvas("rf101_basics","rf101_basics",800,400) ;
  //c->Divide(2) ;
  //c->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe->GetYaxis()->SetTitleOffset(1.6) ; xframe->Draw() ;
  c->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe2->GetYaxis()->SetTitleOffset(1.6) ; xframe2->Draw() ;
	
}

void roofit(TH1F * his)
{
	gSystem->Load("libRooFit") ;
	//using namespace RooFit ;

  // S e t u p   m o d e l 
  // ---------------------
	RooRealProxy realProxy;

  // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
  RooRealVar Mrec("Mrec","Mrec",-45,+45) ;
  RooRealVar mean("mean","mean of gaussian",0,-45,45) ;
  RooRealVar sigma("sigma","width of gaussian",1.3,0.01,10) ;

  RooRealVar lambda("lambda","asymetry",0.001, -1, 1);
  RooRealVar rho_left("rho_left","left tail",0.001, 0, 100);
  RooRealVar rho_right("rho_right","right tail",0.001, 0, 100);

  // Build gaussian p.d.f in terms of x,mean and sigma
  //RooGaussian gauss("gauss","gaussian PDF",Mrec,mean,sigma) ;  
  RooBukinPdf gauss("bukin","bukin PDF",Mrec,mean,sigma, lambda, rho_left, rho_right) ;  

	RooRealVar poly_c1("poly_c1","coefficient of x^1 term",0,-10,10);
	RooPolynomial bkgd_poly("bkgd_poly", "linear function for background", Mrec, RooArgList(poly_c1));	

  // Construct plot frame in 'x'
  //RooPlot* xframe = Mrec.frame(Title("Gaussian p.d.f.")) ;

	RooRealVar peak_yield("peak_yield", "yield signal peak", 3000, 0, 1000000);
	RooRealVar bkgd_yield("bkgd_yield", "yield of background", 500, 0, 1000000);
	RooArgList shapes;
	RooArgList yields;
	shapes.add(bkgd_poly);      yields.add(bkgd_yield);
	shapes.add(gauss);  yields.add(peak_yield);
	RooAddPdf  totalPdf("totalPdf", "sum of signal and background PDF's", shapes, yields);


  // P l o t   m o d e l   a n d   c h a n g e   p a r a m e t e r   v a l u e s
  // ---------------------------------------------------------------------------

  // Plot gauss in frame (i.e. in x) 
  //gauss.plotOn(xframe) ;

  // Change the value of sigma to 3
  //sigma.setVal(3) ;

  // Plot gauss in frame (i.e. in x) and draw frame on canvas
  //gauss.plotOn(xframe,LineColor(kRed)) ;
  

	// Create category observable c that serves as index for the ROOT histograms
  RooCategory rooCategory("c","c") ;
  rooCategory.defineType("his") ;
  // G e n e r a t e   e v e n t s 
  // -----------------------------

  // Generate a dataset of 1000 events in x from gauss
  //RooDataSet* data = gauss.generate(x,10000) ;  
	RooArgSet ntupleVarSet(Mrec);
	//RooDataSet * data = new RooDataSet("tree", "tree",  tree,  ntupleVarSet);
	//RooDataSet * data = new RooDataSet("dh", "dh", Mrec, Index(rooCategory),  Import("his",*his));
  RooRealVar x("x","x",-45,+45) ;
	RooDataHist * data = new RooDataHist("dh", "dh", Mrec, Import(*his));
  
  // Make a second plot frame in x and draw both the 
  // data and the p.d.f in the frame
  RooPlot* xframe2 = Mrec.frame(Title("Gaussian p.d.f. with data")) ;
  data->plotOn(xframe2) ;
  totalPdf.plotOn(xframe2) ;
	totalPdf.plotOn(xframe2,Components(bkgd_poly),LineStyle(kDashed)) ;
  

  // F i t   m o d e l   t o   d a t a
  // -----------------------------

  // Fit pdf to data
  //gauss.fitTo(*data) ;
	//totalPdf.fitTo(*data, Extended());
	totalPdf.fitTo(*data, FitOptions("qmh"));

  // Print values of mean and sigma (that now reflect fitted values and errors)
  mean.Print() ;
  sigma.Print() ;

  // Draw all frames on a canvas
  TCanvas* c = new TCanvas("rf101_basics","rf101_basics",800,400) ;
  //c->Divide(2) ;
  //c->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe->GetYaxis()->SetTitleOffset(1.6) ; xframe->Draw() ;
  c->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe2->GetYaxis()->SetTitleOffset(1.6) ; xframe2->Draw() ;
	
}

void myroofit(TH1F * his)
{
	gSystem->Load("libRooFit") ;
	gSystem->Load("libFit") ;
	//RooRealProxy realProxy;

  // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
  RooRealVar Mrec("Mrec","Mrec",-45,+45);
  RooRealVar sigma("sigma","sigma",1.3,0.01,10, "MeV") ;

	std::vector<RooRealVar> staple(7);
	for(int i=0;i<staple.size();i++)
	{
		char buf[128];
		sprintf(buf, "staple-%d",i );
		staple[i] = RooRealVar(buf, buf, -45, +45);
	}
	std::vector<RooRealVar> N(2);
	N[0] = RooRealVar("nl", "Left power", 2,  0.1, 100);
	N[1] = RooRealVar("nr", "Right power", 2,  0.1, 100);

	ModifiedCrystalBall mcb("ModCB", "Modified CrystalBall",  Mrec,  sigma,  
			staple[0], 
			staple[1], 
			staple[2], 
			staple[3], 
			staple[4], 
			staple[5], 
			staple[6], 
			N[0], 
			N[1]);

	RooRealVar poly_c1("poly_c1","coefficient of x^1 term",0,-10,10);
	RooPolynomial bkgd_poly("bkgd_poly", "linear function for background", Mrec, RooArgList(poly_c1));	

  // Construct plot frame in 'x'
  //RooPlot* xframe = Mrec.frame(Title("Gaussian p.d.f.")) ;

	RooRealVar peak_yield("peak_yield", "yield signal peak", 3000, 0, 1000000);
	RooRealVar bkgd_yield("bkgd_yield", "yield of background", 500, 0, 1000000);
	RooArgList shapes;
	RooArgList yields;
	shapes.add(bkgd_poly);      yields.add(bkgd_yield);
	shapes.add(mcb);  yields.add(peak_yield);
	RooAddPdf  totalPdf("totalPdf", "sum of signal and background PDF's", shapes, yields);


  // P l o t   m o d e l   a n d   c h a n g e   p a r a m e t e r   v a l u e s
  // ---------------------------------------------------------------------------

  // Plot gauss in frame (i.e. in x) 
  //gauss.plotOn(xframe) ;

  // Change the value of sigma to 3
  //sigma.setVal(3) ;

  // Plot gauss in frame (i.e. in x) and draw frame on canvas
  //gauss.plotOn(xframe,LineColor(kRed)) ;
  

	// Create category observable c that serves as index for the ROOT histograms
  RooCategory rooCategory("c","c") ;
  rooCategory.defineType("his") ;
  // G e n e r a t e   e v e n t s 
  // -----------------------------

  // Generate a dataset of 1000 events in x from gauss
  //RooDataSet* data = gauss.generate(x,10000) ;  
	RooArgSet ntupleVarSet(Mrec);
	//RooDataSet * data = new RooDataSet("tree", "tree",  tree,  ntupleVarSet);
	//RooDataSet * data = new RooDataSet("dh", "dh", Mrec, Index(rooCategory),  Import("his",*his));
	RooDataHist * data = new RooDataHist("dh", "dh", Mrec, Import(*his));
  
  // Make a second plot frame in x and draw both the 
  // data and the p.d.f in the frame
  RooPlot* xframe2 = Mrec.frame(Title("Gaussian p.d.f. with data")) ;
  data->plotOn(xframe2) ;
  totalPdf.plotOn(xframe2) ;
	totalPdf.plotOn(xframe2,Components(bkgd_poly),LineStyle(kDashed)) ;
  

  // F i t   m o d e l   t o   d a t a
  // -----------------------------

  // Fit pdf to data
  //gauss.fitTo(*data) ;
	//totalPdf.fitTo(*data, Extended());
	totalPdf.fitTo(*data, FitOptions("qmh"));

  // Print values of mean and sigma (that now reflect fitted values and errors)
  //mean.Print() ;
  sigma.Print() ;
	for(auto s : staple)
	{
		s.Print();
	}
	for(auto n : N)
	{
		n.Print();
	}

  // Draw all frames on a canvas
  TCanvas* c = new TCanvas("rf101_basics","rf101_basics",800,400) ;
  //c->Divide(2) ;
  //c->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe->GetYaxis()->SetTitleOffset(1.6) ; xframe->Draw() ;
  c->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe2->GetYaxis()->SetTitleOffset(1.6) ; xframe2->Draw() ;
	
}


void gaus_fit(TH1 * h)
{
	TF1 * fun = new TF1("fun",  "gaus(0)+gaus(3)+gaus(6)+pol0(9)");
	for(int i=0;i<3;i++)
	{
		fun->SetParameter(3*i+1, 0);
	}
	fun->SetParameter(0, h->GetMaximum());
	fun->SetParameter(3, h->GetMaximum()/10);
	fun->SetParameter(6, h->GetMaximum()/100);

	fun->SetParameter(2, 1);
	fun->SetParameter(5, 5);
	fun->SetParameter(8, 25);
	fun->SetParameter(9, 0);
	h->Fit("fun");

}
