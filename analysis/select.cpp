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
#include <boost/program_options.hpp>
#include <boost/functional/hash.hpp>


#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TCut.h>
#include <TH1F.h>

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
  for (directory_iterator itr(dir); itr!=directory_iterator(); ++itr)
  {
    auto file = absolute(itr->path());
    auto & name = file.string();
    if(file.extension()==".root") lst.push_back(name);
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

TTree * load_tree(string tree_name, string file_name)
{
  TChain * event = new TChain(tree_name.c_str(),tree_name.c_str());
  event->AddFile(file_name.c_str());
  return event;
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

int main(int argc, char ** argv)
{
  namespace po=boost::program_options;
  po::options_description opt_desc("Allowed options");
  std::string tree_name;
  std::string tree_file;
  unsigned long long N;
  unsigned long MAX_EVERY;
  opt_desc.add_options()
    ("help,h","Print this help")
    ("input", po::value<std::string>(&tree_file)->default_value("sample.root"), "Root file (.root) with the data")
    ("max_every", po::value<unsigned long>(&MAX_EVERY)->default_value(1e4), "Maximum every for printing")
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

  //auto lst = get_file_list_in_dir(tree_file);
  //for(auto & name : lst) cout<< name << endl;
  //TTree * event = load_tree(lst);

  TApplication theApp("root_app", &argc, argv);
  //TCut cut = "kin_chi2 < 40 && pid_chi2 < 40 && Mrec > 3.04 && Mrec<3.14";
  //event->Draw("Mrec",cut && "KK");
  //new TCanvas;
  //event->Draw("Mrec",cut );

  string filename = tree_file;
  RootEvent event(load_tree("event",filename));
  RootMC    mc(load_tree("mc",filename));
  RootMdc   mdc(load_tree("mdc",filename));
  //RootMCTopo   mctopo(load_tree("mctopo",filename));
  RootMCTopo mctopo(load_tree("mctopo",filename));
  int NBINS_KK= 200;
  int NBINS_UU= 500;
  double MRANGE=0.09;
  double MAX_RECOIL_MASS=MJPSI_SHIFT+MRANGE*0.5;
  double MIN_RECOIL_MASS=MJPSI_SHIFT-MRANGE*0.5;
  TFile file("result.root","RECREATE");
  auto hMrecKK = new TH1F("hMrecKK","#pi^{+}#pi^{-} recoil mass for KK channel",NBINS_KK,mshift(MIN_RECOIL_MASS),mshift(MAX_RECOIL_MASS));
  auto hMrecUU = new TH1F("hMrecUU","#pi^{+}#pi^{-} recoil mass for uu channel",NBINS_UU,mshift(MIN_RECOIL_MASS),mshift(MAX_RECOIL_MASS));

  std::cout << "Init event_treeKK" << std::endl;

  auto event_treeKK = event.fChain->CloneTree(0);
  event_treeKK->SetName("eventKK");
  event_treeKK->SetTitle("KK events");

  auto mctopo_treeKK = mctopo.fChain->CloneTree(0);
  mctopo_treeKK->SetName("mctopoKK");
  mctopo_treeKK->SetTitle("Monte Carlo topology for KK selection");

  auto event_treeUU = event.fChain->CloneTree(0);
  event_treeUU->SetName("eventUU");
  event_treeUU->SetTitle("UU events");

  auto mctopo_treeUU = mctopo.fChain->CloneTree(0);
  mctopo_treeUU->SetName("mctopoUU");
  mctopo_treeUU->SetTitle("Monte Carlo topology for UU selection");

  Long64_t nentries = event.fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  double PID_CHI2=20;
  double KIN_CHI2=40;

  Long64_t N0=0;
  Long64_t NKK=0;
  Long64_t Nuu=0;

  const int KAON = 0;
  const int MUON = 1;
  //std::map< int, std::map< int, std::map < int, long int> > > theCounter;




  std::unordered_map<Index_t, long int, IndexHash_t > theCounter;

  
  for (Long64_t jentry=0; jentry<nentries;jentry++)
  {
    Long64_t ientry;
    ientry = event.LoadTree(jentry);
    if (ientry < 0) 
    {
      cerr << "ERROR: event.LoadTree() retururn " << ientry << endl;
      exit(1);
    }
    ientry = mc.LoadTree(jentry);
    if (ientry < 0) 
    {
      cerr << "ERROR: mc.LoadTree() retururn " << ientry << endl;
      exit(1);
    }
    ientry = mctopo.LoadTree(jentry);
    if (ientry < 0) 
    {
      cerr << "ERROR: mctopo.LoadTree() retururn " << ientry << endl;
      exit(1);
    }

    nb = event.fChain->GetEntry(jentry);   nbytes += nb;
    nb = mc.fChain->GetEntry(jentry);
    nb = mctopo.fChain->GetEntry(jentry);
    N0++; //count total number of events proceed
    if(MIN_RECOIL_MASS <= event.Mrec && event.Mrec <= MAX_RECOIL_MASS)
      if(event.pid_chi2 <= PID_CHI2)
        if(event.kin_chi2 <= KIN_CHI2)
        {
          if(event.KK==1) 
          {
            NKK++;
            hMrecKK->Fill(mshift(event.Mrec));
            //hpid_chi2KK->Fill(pid_chi2);
            //hkin_chi2KK->Fill(kin_chi2);
            event_treeKK->Fill();
            mctopo_treeKK->Fill();
          }
          if(event.uu==1) 
          {
            hMrecUU->Fill(mshift(event.Mrec));
            //hpid_chi2UU->Fill(pid_chi2);
            //hkin_chi2UU->Fill(kin_chi2);
            event_treeUU->Fill();
            mctopo_treeUU->Fill();
            Nuu++;
          }
          Index_t index = {event.channel, event.sign, event.ngtrack};
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
  hMrecKK->Write();
  hMrecUU->Write();
  event_treeKK->Write();
  mctopo_treeKK->Write();
  event_treeUU->Write();
  mctopo_treeUU->Write();
  file.Close();

  //std::cout << theCounter[{0,0,4}] << std::endl;
  theApp.Run();
}
