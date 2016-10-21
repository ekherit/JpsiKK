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
#include <limits>

#include <boost/format.hpp>
#include <boost/program_options.hpp>
//#include <boost/algorithm/string.hpp>

#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>
extern void InitGui();
VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
TROOT root("kkfit","kkfit", initfuncs);


#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TTree.h>
#include <TCut.h>

#include "hashlist.h"

#include "fit.h"

int main(int argc,  char ** argv)
{
  std::ios_base::sync_with_stdio(false);
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
  std::string nbin_str;
  std::map<std::string,int> Nbin;
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
    ("nobgslope","No bg slope")
    ("nobg","No background events")
    ("nograd","No radiative gaus")
    ("gaus_rad","Gaus rad")
    ("par",po::value<std::string>(&OPT_PARAM_CONFIG_FILE), "Parameters configuration file")
 //   ("Nbin",po::value<int>(&Nbin),"Override number of  bins")
    ("nbin",po::value<std::string>(&nbin_str),"Comma separated number of bins")
    ;
  po::positional_options_description pos;
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

  OPT_NOBGSLOPE=opt.count("nobgslope");
  OPT_NOBG = opt.count("nobg");
  OPT_NOGAUSRAD = opt.count("nograd");

  std::cout << " " << suffix << " " << file_name << std::endl;
  
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
    
    int N  = h->GetNbinsX();
    auto it = Nbin.find(suffix);
    if( it != Nbin.end()) N = it->second;

    double Mmin = h->GetXaxis()->GetXmin();
    double Mmax = h->GetXaxis()->GetXmax();
    auto varexp = boost::format("(Mrec-3.097)*1000 >> Mrec%s(%d,%f,%f)") % suffix % N % Mmin % Mmax; 
    std::cout << "varexp = " << varexp << std::endl;
    tree->Draw(varexp.str().c_str(),hash_cut && cut.c_str(),"goff");
    h = (TH1F*) tree->GetHistogram();
    h->SetName(suffix.c_str());
    h->SetTitle(title.c_str());
    return h;
  };

  std::list<std::string> strBinsLst;
  if(opt.count("nbin"))
  {
    boost::split(strBinsLst,nbin_str,boost::is_any_of(", "));
  }

  if(opt.count("suffix"))
  {
    std::vector<std::string> suffix_list;
    boost::split(suffix_list,suffix,boost::is_any_of(", "));
    std::list<TH1*> his_lst;
    auto nbin = strBinsLst.begin();
    for(auto str : suffix_list)
    {
      if(nbin!=strBinsLst.end())
      {
        std::cout << *nbin << std::endl;
        Nbin[str] = stoi(*nbin);
        nbin++;
      }
      his_lst.push_back(create_histogram(str));
    }
    fit(his_lst);
    theApp.Run();
  }
	theApp.Run();
}
