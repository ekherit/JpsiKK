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

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>


#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TCut.h>

#include "libFit.h"

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
  opt_desc.add_options()
    ("help,h","Print this help")
    ("input", po::value<std::string>(&tree_file), "Root file (.root) with the data")
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

  if(opt.count("help") && !opt.count("input"))
  {
    std::cout << "Usage: CombMcbFit <root_file>" << std::endl;
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

  string filename = "sample.root";
  RootEvent event(load_tree("event",filename));
  RootMC    mc(load_tree("mc",filename));
  RootMdc   mdc(load_tree("mdc",filename));

  Long64_t nentries = event.fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
  {
    Long64_t ientry;
    ientry = event.LoadTree(jentry);
    if (ientry < 0) break;
    ientry = mc.LoadTree(jentry);
    if (ientry < 0) break;
    nb = event.fChain->GetEntry(jentry);   nbytes += nb;
    nb = mc.fChain->GetEntry(jentry);
    // if (Cut(ientry) < 0) continue;
    if(event.kin_chi2<40)
    {
      cout << jentry << " " << event.Mrec << " " << event.p[0] << " " << mc.p[0] << " " << (event.p[0]-mc.p[0])/event.p[0]*100 << endl;
    }

  }

  //Selector s(0);
  //event->Process(&s);
	theApp.Run();
}
