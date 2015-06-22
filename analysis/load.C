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
  TChain * mc_event   =  load_tree("../../mc09-2",".root","event");
  TChain * mc_mdc     =  load_tree("../../mc09-2",".root","mdc");
  TChain * mc_emc     =  load_tree("../../mc09-2",".root","emc");

  TChain * data_event =  load_tree("../../2009",".root","event");
  TChain * data_mdc   =  load_tree("../../2009",".root","mdc");
  TChain * data_emc   =  load_tree("../../2009",".root","emc");
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
