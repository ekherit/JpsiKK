// =====================================================================================
//
//       Filename:  select.cpp
//
//    Description:  Select the data
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
#include <list>
#include <string> 
#include <regex>
#include <unordered_map>

#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/program_options.hpp>
#include <boost/functional/hash.hpp>


#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TCut.h>
#include <TH1F.h>
#include <TH1D.h>

#include "libFit.h"

#include "mctopo/mctopo.h"

using namespace std;

list<string> get_file_list_in_dir(string dirname)
{
  list<string> lst;
  using namespace boost::filesystem;
  boost::filesystem::path dir(dirname);
  regex re(".+.root");
  smatch match;
  if(is_directory(dir))
  {
    for (directory_iterator itr(dir); itr!=directory_iterator(); ++itr)
    {
      auto file = absolute(itr->path());
      auto & name = file.string();
      if(file.extension()==".root") lst.push_back(name);
    }
  }
  else 
  {
    if(is_regular_file(dir))
    {
      auto file = absolute(dir);
      auto & name = file.string();
      if(file.extension()==".root") lst.push_back(name);
    }
  }
  lst.sort();
  return lst;
}


TTree * load_tree(const list<string> &  lst)
{
  TChain * event = new TChain("event","event");

  list<TChain*> chain_list;
  chain_list.push_back(new TChain("mdc","Main Drift Chamber"));
  chain_list.push_back(new TChain("mc","Monte Carlo information"));
  chain_list.push_back(new TChain("mctopo","Monte Carlo topology"));

  for(auto & name : lst)
  {
    event->AddFile(name.c_str());
    for(auto  chain : chain_list)
    {
      chain -> AddFile(name.c_str());
    }
  }
  for(auto  chain : chain_list)
  {
    event->AddFriend(chain);
  }
  return event;
}

TTree * load_tree(string tree_name, const list<string> &  lst)
{
  TChain * tree = new TChain(tree_name.c_str(),tree_name.c_str());
  for(auto & name : lst)
  {
    std::cout << "Adding file... " << name;
    int rc =  tree->AddFile(name.c_str());
    if(rc==1) std::cout << " OK" << std::endl;
    else std::cout << " ERROR" << std::endl;
  }
  return tree;
}

TTree * load_tree(string tree_name, string file_name)
{
  TChain * tree = new TChain(tree_name.c_str(),tree_name.c_str());
  tree->AddFile(file_name.c_str());
  return tree;
}

const double MSCALE=1000;
const double MJPSI_SHIFT=3.097;
double mshift(double m) 
{
  return MSCALE*(m-MJPSI_SHIFT);
}
struct Index_t
{
  int channel;
  int charge;
  int tracks;
  bool operator==(const Index_t & other) const
  {
    //return std::tie(channel,charge,tracks) == std::tie(other.channel,other.charge,other.tracks);
    return channel==other.channel && charge == other.charge && tracks==other.tracks;
  }

  std::size_t hash(void) const
  {
    std::size_t seed = 0;
    boost::hash_combine(seed, channel);
    boost::hash_combine(seed, charge);
    boost::hash_combine(seed, tracks);
    return seed;
  }
};

struct IndexHash_t
{
  std::size_t operator()(Index_t const& i) const
  {
    return i.hash();
  }
};

#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>
extern void InitGui();
VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
TROOT root("select","select", initfuncs);

enum {KAON, MUON};

int main(int argc, char ** argv)
{
  namespace po=boost::program_options;
  po::options_description opt_desc("Allowed options");
  std::string tree_name;
  //std::string tree_file;
  std::vector< std::string> files = {"sample.root"};
  std::string output_file;
  std::string event_tree_name;
  std::string mctopo_tree_name;
  std::string tree_suffix;
  unsigned long long NMAX;
  unsigned long MAX_EVERY;
  opt_desc.add_options()
    ("help,h","Print this help")
    //("input", po::value<std::string>(&tree_file)->default_value("sample.root"), "Root file (.root) with the data")
    ("input", po::value<std::vector< std::string> >(&files), "Root file (.root) with the data")
    ("output", po::value<std::string>(&output_file)->default_value("result.root"), "Output file")
    ("event_tree",  po::value<std::string>(&event_tree_name)->default_value("event"), "Event tree name")
    ("mctopo_tree", po::value<std::string>(&mctopo_tree_name)->default_value("mctopo"), "mctopo tree name")
    ("tree_suffix", po::value<std::string>(&tree_suffix)->default_value(""), "suffix for the tree")
    ("max_every", po::value<unsigned long>(&MAX_EVERY)->default_value(1e4), "Maximum every for printing")
    ("N", po::value<unsigned long long>(&NMAX)->default_value(std::numeric_limits<unsigned long long>::max()), "Maximum event number to proceed")
    ("rad",  "Fit by rad gaus")
    ("simple",  "Simple model gaus + power + exp")
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


  if(opt.count("help"))
  {
    std::cout << "Usage: " << argv[0] <<  " <root_file>" << std::endl;
    std::clog << opt_desc;
    return 0;
  }

  std::list < std::string > file_list;
  for( auto dir : files)
  {
    file_list.merge(get_file_list_in_dir(dir));
  };


  RootEvent  event(load_tree(event_tree_name+tree_suffix,file_list));
  RootMCTopo mctopo(load_tree(mctopo_tree_name+tree_suffix,file_list));
  //RootMC     mc(load_tree("mc"+tree_suffix,file_list));
  RootMdc    mdc(load_tree("mdc"+tree_suffix,file_list));

  //TApplication theApp("root_app", &argc, argv);

  int NBINS_KK= 200;
  int NBINS_UU= 500;
  double MRANGE=0.09;
  double MAX_RECOIL_MASS=MJPSI_SHIFT+MRANGE*0.5;
  double MIN_RECOIL_MASS=MJPSI_SHIFT-MRANGE*0.5;

  TFile file(output_file.c_str(),"RECREATE");
  auto hMrecKK = new TH1F("hMrecKK","#pi^{+}#pi^{-} recoil mass for KK channel",NBINS_KK,mshift(MIN_RECOIL_MASS),mshift(MAX_RECOIL_MASS));
  auto hMrecUU = new TH1F("hMrecUU","#pi^{+}#pi^{-} recoil mass for uu channel",NBINS_UU,mshift(MIN_RECOIL_MASS),mshift(MAX_RECOIL_MASS));


  std::unordered_map<Index_t, TH1D*, IndexHash_t > hMrec;

  
  hMrec[{KAON,0,4}] = new TH1D("hMrecKK","#pi^{+}#pi^{-} recoil mass for KK channel",NBINS_KK,mshift(MIN_RECOIL_MASS),mshift(MAX_RECOIL_MASS));
  hMrec[{MUON,0,4}] = new TH1D("hMrecUU","#pi^{+}#pi^{-} recoil mass for UU channel",NBINS_UU,mshift(MIN_RECOIL_MASS),mshift(MAX_RECOIL_MASS));
  hMrec[{KAON,-1,3}] = new TH1D("hMrecKm3","#pi^{+}#pi^{-} recoil mass for negative K channel with 3+4 tracks",NBINS_KK,mshift(MIN_RECOIL_MASS),mshift(MAX_RECOIL_MASS));
  hMrec[{MUON,-1,3}] = new TH1D("hMrecUm3","#pi^{+}#pi^{-} recoil mass for negative U channel with 3+4 tracks",NBINS_UU,mshift(MIN_RECOIL_MASS),mshift(MAX_RECOIL_MASS));

  hMrec[{KAON,1,3}] = new TH1D("hMrecKp3","#pi^{+}#pi^{-} recoil mass for positive K channel with 3+4 tracks",NBINS_KK,mshift(MIN_RECOIL_MASS),mshift(MAX_RECOIL_MASS));
  hMrec[{MUON,1,3}] = new TH1D("hMrecUp3","#pi^{+}#pi^{-} recoil mass for positive U channel with 3+4 tracks",NBINS_UU,mshift(MIN_RECOIL_MASS),mshift(MAX_RECOIL_MASS));

  hMrec[{KAON,-1,4}] = new TH1D("hMrecKm4","#pi^{+}#pi^{-} recoil mass for negative K channel with 4 tracks",NBINS_KK,mshift(MIN_RECOIL_MASS),mshift(MAX_RECOIL_MASS));
  hMrec[{MUON,-1,4}] = new TH1D("hMrecUm4","#pi^{+}#pi^{-} recoil mass for negative U channel with 4 tracks",NBINS_UU,mshift(MIN_RECOIL_MASS),mshift(MAX_RECOIL_MASS));

  hMrec[{KAON,1,4}] = new TH1D("hMrecKp4","#pi^{+}#pi^{-} recoil mass for positive K channel with 4 tracks",NBINS_KK,mshift(MIN_RECOIL_MASS),mshift(MAX_RECOIL_MASS));
  hMrec[{MUON,1,4}] = new TH1D("hMrecUp4","#pi^{+}#pi^{-} recoil mass for positive U channel with 4 tracks",NBINS_UU,mshift(MIN_RECOIL_MASS),mshift(MAX_RECOIL_MASS));


  std::cout << "Init event_treeKK" << std::endl;


  auto event_treeKK = event.fChain->CloneTree(0);
  event_treeKK->SetName("eventKK");
  event_treeKK->SetTitle("KK events");

  auto mctopo_treeKK = mctopo.fChain->CloneTree(0);
  mctopo_treeKK->SetName("mctopoKK");
  mctopo_treeKK->SetTitle("Monte Carlo topology for KK selection");


  auto event_treeK = event.fChain->CloneTree(0);
  event_treeK->SetName("eventK");
  event_treeK->SetTitle("K events");

  auto mctopo_treeK = mctopo.fChain->CloneTree(0);
  mctopo_treeK->SetName("mctopoK");
  mctopo_treeK->SetTitle("Monte Carlo topology for K events");

  auto mdc_treeK = mdc.fChain->CloneTree(0);
  mdc_treeK->SetName("mdcK");
  mdc_treeK->SetTitle("mdc info for K events");

  auto event_treeUU = event.fChain->CloneTree(0);
  event_treeUU->SetName("eventUU");
  event_treeUU->SetTitle("UU events");

  auto mctopo_treeUU = mctopo.fChain->CloneTree(0);
  mctopo_treeUU->SetName("mctopoUU");
  mctopo_treeUU->SetTitle("Monte Carlo topology for UU selection");

  auto event_treeU = event.fChain->CloneTree(0);
  event_treeU->SetName("eventU");
  event_treeU->SetTitle("U events");

  auto mctopo_treeU = mctopo.fChain->CloneTree(0);
  mctopo_treeU->SetName("mctopoU");
  mctopo_treeU->SetTitle("Monte Carlo topology for U selection");

  Long64_t nentries = event.fChain->GetEntriesFast();
  //std::cout << "Number of entries in tree (GetEntriesFast): " << nentries << std::endl;
  //std::cout << "Number of entries (GetEntries)" << event.fChain->GetEntries() << std::endl;

  Long64_t nbytes = 0, nb = 0;
  double PID_CHI2=20;
  double KIN_CHI2=40;

  Long64_t N0=0;
  Long64_t NKK=0;
  Long64_t Nuu=0;





  std::unordered_map<Index_t, long int, IndexHash_t > theCounter;

  //TTree::SetMaxTreeSize(100000000000LL);
  
  for (Long64_t jentry=0; jentry<nentries && jentry<NMAX;jentry++)
  {
    Long64_t ientry;
    ientry = event.LoadTree(jentry);
    if (ientry < 0)  break;
    //ientry = mc.LoadTree(jentry);
    //if (ientry < 0)
    //{
    //  cerr << "ERROR: mc.LoadTree() retururn " << ientry << endl;
    //}
    ientry = mctopo.LoadTree(jentry);
    if (ientry < 0) 
    {
      cerr << "ERROR: mctopo.LoadTree() retururn " << ientry << endl;
    }
    ientry = mdc.LoadTree(jentry);

    nb = event.fChain->GetEntry(jentry);   nbytes += nb;
    //nb = mc.fChain->GetEntry(jentry);
    nb = mctopo.fChain->GetEntry(jentry);
    nb = mdc.fChain->GetEntry(jentry);
    N0++; //count total number of events proceed
    if(MIN_RECOIL_MASS <= event.Mrec && event.Mrec <= MAX_RECOIL_MASS)
      if(event.pid_chi2 <= PID_CHI2)
        if(event.kin_chi2 <= KIN_CHI2)
        {
          auto hash = mctopo_hash(&mctopo);
          mctopo.hash = hash;
          Index_t index = {event.channel, event.sign, event.ngtrack};
          if(event.KK==1) 
          {
            NKK++;
            //hMrecKK->Fill(mshift(event.Mrec));
            //hpid_chi2KK->Fill(pid_chi2);
            //hkin_chi2KK->Fill(kin_chi2);
            event_treeKK->Fill();
            mctopo_treeKK->Fill();
            hMrec[index]->Fill(mshift(event.Mrec));
          }
          if(event.uu==1) 
          {
            //hMrecUU->Fill(mshift(event.Mrec));
            //hpid_chi2UU->Fill(pid_chi2);
            //hkin_chi2UU->Fill(kin_chi2);
            event_treeUU->Fill();
            mctopo_treeUU->Fill();
            Nuu++;
            hMrec[index]->Fill(mshift(event.Mrec));
          }
          //if(event.K==1 && event.ngntrack==0 && event.depth[2]<40 && event.depth[3] < 40)
          //if(event.K==1 && event.ngntrack==0 && mdc.E[2]/mdc.p[2]>0.26 && mdc.E[3]/mdc.p[3]>0.26)
          if(event.K==1)
          {
            if( fabs(cos(event.theta[2])) < 0.8 && fabs(cos(event.theta[3])) <0.8)
            {
              mctopo_treeK->Fill();
              event_treeK->Fill();
              mdc_treeK->Fill();
              hMrec[index]->Fill(mshift(event.Mrec));
            }
          }
          if(event.u==1)
          {
            event_treeU->Fill();
            mctopo_treeU->Fill();
            hMrec[index]->Fill(mshift(event.Mrec));
          }
          theCounter[index]++;
        }
    static unsigned long every = 1;
    if(jentry > every*10 && every<MAX_EVERY) every*=10;
    if(jentry % every==0 || jentry==nentries-1)
    {
      std::cout << setw(15) << jentry << setw(15) << NKK << setw(15) << Nuu;
      std::cout << "  hash = " << setw(10) << std::hex << mctopo_hash(&mctopo) << "    " << std::dec << mctopo_info(&mctopo);
      std::cout << std::endl;
    }
  }
  //hMrecKK->Write();
  //hMrecUU->Write();
  for(auto  his : hMrec)
  {
    his.second->Write();
  }

  event_treeKK->Write();
  mctopo_treeKK->Write();

  event_treeK->Write();
  mctopo_treeK->Write();
  mdc_treeK->Write();

  event_treeUU->Write();
  mctopo_treeUU->Write();
  mctopo_treeU->Write();
  event_treeU->Write();


  //make usefull cuts 
  TCut signal_cut("hash==0x3031ea55 || hash==0x846b22bf || hash==0xcdffabb3 || hash==0xecdbe915");
  signal_cut.SetName("signal_cut");
  //signal_cut.SetTitle("signal: Ψ(2S) → π+π-(J/Ψ → μ+μ-)");
  TCut bg1_cut("hash==0xcfe7a549 || hash==0x34525dce || hash==0x4d7eec07 || hash==0x59476175 || hash==0x758ebb38 || hash==0x7d41c6b4 || hash==0x8a4d7453 || hash==0xcb9192a5 || hash==0x34525dce || hash==0xcfe7a549 || hash==0xdbde283b || hash==0xffd88ffa");
  bg1_cut.SetName("bg1_cut");
  //bg1_cut.SetTitle("background: Ψ(2S) → π+π-(J/Ψ → π+(ρ(770)- → (π0 → ɣɣ)π-)) + rad + c.c.");

  bg1_cut.Write();
  signal_cut.Write();
  file.Close();

  //std::cout << theCounter[{0,0,4}] << std::endl;
  //theApp.Run();
}
