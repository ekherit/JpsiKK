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
#include <boost/format.hpp>


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
  //selection parameters
  double MRANGE;
  double PID_CHI2;
  double KIN_CHI2;
  opt_desc.add_options()
    ("help,h","Print this help")
    //("input", po::value<std::string>(&tree_file)->default_value("sample.root"), "Root file (.root) with the data")
    ("input", po::value<std::vector< std::string> >(&files), "Root file (.root) with the data")
    ("output,o", po::value<std::string>(&output_file)->default_value("result.root"), "Output file")
    ("event_tree",  po::value<std::string>(&event_tree_name)->default_value("event"), "Event tree name")
    ("mctopo_tree", po::value<std::string>(&mctopo_tree_name)->default_value("mctopo"), "mctopo tree name")
    ("tree_suffix", po::value<std::string>(&tree_suffix)->default_value(""), "suffix for the tree")
    ("max_every", po::value<unsigned long>(&MAX_EVERY)->default_value(1e4), "Maximum every for printing")
    ("N", po::value<unsigned long long>(&NMAX)->default_value(std::numeric_limits<unsigned long long>::max()), "Maximum event number to proceed")
    ("mrange", po::value<double>(&MRANGE)->default_value(0.09), "Pion recoil mass range")
    ("kin_chi2", po::value<double>(&KIN_CHI2)->default_value(40), "Kinematic chi square cut")
    ("pid_chi2", po::value<double>(&PID_CHI2)->default_value(20), "Particle id cut")
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
  double MAX_RECOIL_MASS=MJPSI_SHIFT+MRANGE*0.5;
  double MIN_RECOIL_MASS=MJPSI_SHIFT-MRANGE*0.5;

  TFile file(output_file.c_str(),"RECREATE");
  //auto hMrecKK = new TH1F("hMrecKK","#pi^{+}#pi^{-} recoil mass for KK channel",NBINS_KK,mshift(MIN_RECOIL_MASS),mshift(MAX_RECOIL_MASS));
  //auto hMrecUU = new TH1F("hMrecUU","#pi^{+}#pi^{-} recoil mass for uu channel",NBINS_UU,mshift(MIN_RECOIL_MASS),mshift(MAX_RECOIL_MASS));


  std::unordered_map<Index_t, TTree*, IndexHash_t > event_tree;
  std::unordered_map<Index_t, TTree*, IndexHash_t > mctopo_tree;
  std::unordered_map<Index_t, TTree*, IndexHash_t > mdc_tree;

  std::unordered_map<Index_t, TH1D*, IndexHash_t > hMrec;

  int POSNEG=2;
  
  auto make_tree = [&](TTree * parent_tree, string  name,  string  title)
  {
    TTree * tree  = parent_tree->CloneTree(0);
    tree->SetName(name.c_str());
    tree->SetTitle(title.c_str());
    return tree;
  };

  auto make_hMrec = [&](string name, string  title) 
  {
    return new TH1D(name.c_str(),title.c_str(),NBINS_KK,mshift(MIN_RECOIL_MASS),mshift(MAX_RECOIL_MASS));
  };

  //main four track channel
  hMrec[{KAON,0,4}] = make_hMrec("hMrecKK","#pi^{+}#pi^{-} recoil mass for KK channel");
  hMrec[{MUON,0,4}] = make_hMrec("hMrecUU","#pi^{+}#pi^{-} recoil mass for UU channel");

  hMrec[{KAON,-1,3}] = make_hMrec("hMrecKm3","#pi^{+}#pi^{-} recoil mass for negative K channel with 3 tracks");
  hMrec[{MUON,-1,3}] = make_hMrec("hMrecUm3","#pi^{+}#pi^{-} recoil mass for negative U channel with 3 tracks");

  hMrec[{KAON,1,3}] = make_hMrec("hMrecKp3","#pi^{+}#pi^{-} recoil mass for positive K channel with 3 tracks");
  hMrec[{MUON,1,3}] = make_hMrec("hMrecUp3","#pi^{+}#pi^{-} recoil mass for positive U channel with 3 tracks");

  hMrec[{KAON,-1,4}] = make_hMrec("hMrecKm4","#pi^{+}#pi^{-} recoil mass for negative K channel with 4 tracks");
  hMrec[{MUON,-1,4}] = make_hMrec("hMrecUm4","#pi^{+}#pi^{-} recoil mass for negative U channel with 4 tracks");

  hMrec[{KAON,1,4}] = make_hMrec("hMrecKp4","#pi^{+}#pi^{-} recoil mass for positive K channel with 4 tracks");
  hMrec[{MUON,1,4}] = make_hMrec("hMrecUp4","#pi^{+}#pi^{-} recoil mass for positive U channel with 4 tracks");

  hMrec[{KAON,POSNEG,3+4}] = make_hMrec("hMrecK34","#pi^{+}#pi^{-} recoil mass for negative K channel with 3 or 4 tracks");
  hMrec[{KAON,POSNEG,4}]   = make_hMrec("hMrecK4", "#pi^{+}#pi^{-} recoil mass for negative K channel with 4 tracks");

  hMrec[{MUON,POSNEG,3+4}] = make_hMrec("hMrecU34","#pi^{+}#pi^{-} recoil mass for Muon channel with 3 or 4 tracks");
  hMrec[{MUON,POSNEG,4}]   = make_hMrec("hMrecU4","#pi^{+}#pi^{-} recoil mass for Muon channel with 4 tracks");



  std::cout << "Init event_treeKK" << std::endl;


  event_tree[{KAON,0,4}] = make_tree(event.fChain,"eventKK","KK events");
  event_tree[{MUON,0,4}] = make_tree(event.fChain,"eventUU","UU events");


  event_tree[{KAON,-1,4}] = make_tree(event.fChain,"eventKm4","negative K events with 4 tracks");
  event_tree[{MUON,-1,4}] = make_tree(event.fChain,"eventUm4","negative U events with 4 tracks");

  event_tree[{KAON,1,4}] = make_tree(event.fChain,"eventKp4","positive K events with 4 tracks");
  event_tree[{MUON,1,4}] = make_tree(event.fChain,"eventUp4","positive U events with 4 tracks");

  event_tree[{KAON,-1,3}] = make_tree(event.fChain,"eventKm3","negative K events with 3 tracks");
  event_tree[{MUON,-1,3}] = make_tree(event.fChain,"eventUm3","negative U events with 3 tracks");

  event_tree[{KAON,1,3}] = make_tree(event.fChain,"eventKp3","positive K events with 3 tracks");
  event_tree[{MUON,1,3}] = make_tree(event.fChain,"eventUp3","positive U events with 3 tracks");


  event_tree[{KAON,POSNEG,3+4}] = make_tree(event.fChain,"eventK34","K events with 3 or 4 tracks");
  event_tree[{MUON,POSNEG,3+4}] = make_tree(event.fChain,"eventU34","U events with 3 or 4 tracks");

  event_tree[{KAON,POSNEG,4}] = make_tree(event.fChain,"eventK4","K events with 4 tracks");
  event_tree[{MUON,POSNEG,4}] = make_tree(event.fChain,"eventU4","U events with 4 tracks");


  mctopo_tree[{KAON,0,4}] = make_tree(mctopo.fChain,"mctopoKK","KK mctopos");
  mctopo_tree[{MUON,0,4}] = make_tree(mctopo.fChain,"mctopoUU","UU mctopos");


  mctopo_tree[{KAON,-1,4}] = make_tree(mctopo.fChain,"mctopoKm4","negative K mctopos with 4 tracks");
  mctopo_tree[{MUON,-1,4}] = make_tree(mctopo.fChain,"mctopoUm4","negative U mctopos with 4 tracks");

  mctopo_tree[{KAON,1,4}] = make_tree(mctopo.fChain,"mctopoKp4","positive K mctopos with 4 tracks");
  mctopo_tree[{MUON,1,4}] = make_tree(mctopo.fChain,"mctopoUp4","positive U mctopos with 4 tracks");

  mctopo_tree[{KAON,-1,3}] = make_tree(mctopo.fChain,"mctopoKm3","negative K mctopos with 3 tracks");
  mctopo_tree[{MUON,-1,3}] = make_tree(mctopo.fChain,"mctopoUm3","negative U mctopos with 3 tracks");

  mctopo_tree[{KAON,1,3}] = make_tree(mctopo.fChain,"mctopoKp3","positive K mctopos with 3 tracks");
  mctopo_tree[{MUON,1,3}] = make_tree(mctopo.fChain,"mctopoUp3","positive U mctopos with 3 tracks");


  mctopo_tree[{KAON,POSNEG,3+4}] = make_tree(mctopo.fChain,"mctopoK34","K mctopos with 3 or 4 tracks");
  mctopo_tree[{MUON,POSNEG,3+4}] = make_tree(mctopo.fChain,"mctopoU34","U mctopos with 3 or 4 tracks");

  mctopo_tree[{KAON,POSNEG,4}] = make_tree(mctopo.fChain,"mctopoK4","K mctopos with 4 tracks");
  mctopo_tree[{MUON,POSNEG,4}] = make_tree(mctopo.fChain,"mctopoU4","U mctopos with 4 tracks");





  auto mdc_treeK = mdc.fChain->CloneTree(0);
  mdc_treeK->SetName("mdcK");
  mdc_treeK->SetTitle("mdc info for K events");

  Long64_t nentries = event.fChain->GetEntriesFast();
  //std::cout << "Number of entries in tree (GetEntriesFast): " << nentries << std::endl;
  //std::cout << "Number of entries (GetEntries)" << event.fChain->GetEntries() << std::endl;

  Long64_t nbytes = 0, nb = 0;

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
    if(event.run < 0) 
    {
      //ientry = mc.LoadTree(jentry);
      //if (ientry < 0)
      //{
      //  cerr << "ERROR: mc.LoadTree() retururn " << ientry << endl;
      //}
      ientry = mctopo.LoadTree(jentry);
      if (ientry < 0) 
      {
        cerr << "ERROR: mctopo.LoadTree() retururn " << ientry << endl;
        continue;
      }
    }
    ientry = mdc.LoadTree(jentry);

    nb = event.fChain->GetEntry(jentry);   nbytes += nb;
    if(nb < 0 )
    {
      cout <<  ientry << "run = " << event.run << " nb=" << nb << endl;
      continue;
    }
    //nb = mc.fChain->GetEntry(jentry);
    if(event.run<0)
    {
      nb = mctopo.fChain->GetEntry(jentry);
    }
    nb = mdc.fChain->GetEntry(jentry);
    N0++; //count total number of events proceed
    unsigned long hash=0;
    if(event.run<0)  
    {
      hash = mctopo_hash(&mctopo);
      mctopo.hash = hash;
    }

    bool RecoilCut = MIN_RECOIL_MASS <= event.Mrec && event.Mrec <= MAX_RECOIL_MASS;
    bool KinematicCut = event.kin_chi2 <= KIN_CHI2;
    bool PidCut = event.pid_chi2 <= PID_CHI2;
    bool ThetaCut =  fabs(cos(event.theta[2])) < 0.8 && fabs(cos(event.theta[3])) <0.8;
    if(    
           RecoilCut 
        && KinematicCut 
        && PidCut 
        && ThetaCut
      )
    {
      if(event.KK == 1 || event.uu ==1)
      {
        Index_t index = {event.channel, event.sign, event.ngtrack};
        event_tree[index]->Fill();
        if(event.run<0) mctopo_tree[index]->Fill();
        hMrec[index]->Fill(mshift(event.Mrec));
        theCounter[index]++;
      }
      if(event.K == 1 || event.u==1)
      {
        if(event.ngntrack==0  && event.kin_chi2<10)
        {
          Index_t index = {event.channel, event.sign, event.ngtrack};
          event_tree[index]->Fill();
          if(event.run<0) mctopo_tree[index]->Fill();
          hMrec[index]->Fill(mshift(event.Mrec));
          theCounter[index]++;

          if(event.ngtrack == 3 || event.ngtrack == 4) 
          {
            index = {event.channel, POSNEG, 3+4};
            event_tree[index]->Fill();
            if(event.run<0) mctopo_tree[index]->Fill();
            hMrec[index]->Fill(mshift(event.Mrec));
            theCounter[index]++;
          }
          if(event.ngtrack == 4 )
          {
            index = {event.channel, POSNEG, 4};
            event_tree[index]->Fill();
            if(event.run<0) mctopo_tree[index]->Fill();
            hMrec[index]->Fill(mshift(event.Mrec));
            theCounter[index]++;
          }
        }
      }
    }
    static unsigned long every = 1;
    static unsigned long head = 0;
    if(jentry > every*10 && every<MAX_EVERY) every*=10;
    if(jentry % every==0 || jentry==nentries-1)
    {
      if(head++%20 == 0)
      {
        std::cout << boost::format("#%14s %8s %10s %10s %10s %10s %10s %10s")
          % "entry"
          % "run"
          % "NKK"
          % "NK4"
          % "NK3+4"
          % "NUU"
          % "NU4"
          % "NU3+4";
        std::cout << endl;
      }
      std::cout << boost::format("%15d %8d %10d %10d %10d %10d %10d %10d") % jentry % event.run 
        %  theCounter[{KAON,0,4}] 
        %  theCounter[{KAON,POSNEG,4}] 
        %  theCounter[{KAON,POSNEG,3+4}] 
        %  theCounter[{MUON,0,4}] 
        %  theCounter[{MUON,POSNEG,4}] 
        %  theCounter[{MUON,POSNEG,3+4}];
      if(event.run<0)
      {
        std::cout << "  hash = " << setw(10) << std::hex << hash << "    " << std::dec << mctopo_info(&mctopo);
      }
      std::cout << std::endl;
    }
  }
  for(auto  his : hMrec)
  {
    his.second->Write();
  }


  for(auto & h : hMrec) h.second->Write();
  for(auto & t : event_tree) t.second->Write();
  for(auto & t : mctopo_tree) t.second->Write();
  



  //make usefull cuts 
  TCut KK_cut("hash==0x3031ea55 || hash==0x846b22bf || hash==0xcdffabb3 || hash==0xecdbe915");
  KK_cut.SetName("KK");

  TCut UU_cut("hash==0x1397e7e7 || hash == 0x8ac60398 || hash == 0xaca004d7 || hash == 0xcf7de4a7 || hash == 0xdaa0af44");
  UU_cut.SetName("uu");

  TCut bg1_cut("hash==0xcfe7a549 || hash==0x34525dce || hash==0x4d7eec07 || hash==0x59476175 || hash==0x758ebb38 || hash==0x7d41c6b4 || hash==0x8a4d7453 || hash==0xcb9192a5 || hash==0x34525dce || hash==0xcfe7a549 || hash==0xdbde283b || hash==0xffd88ffa");
  bg1_cut.SetName("bg1_cut");
  //bg1_cut.SetTitle("background: Ψ(2S) → π+π-(J/Ψ → π+(ρ(770)- → (π0 → ɣɣ)π-)) + rad + c.c.");

  bg1_cut.Write();
  KK_cut.Write();
  UU_cut.Write();

  file.Close();

  //std::cout << theCounter[{0,0,4}] << std::endl;
  //theApp.Run();
}
