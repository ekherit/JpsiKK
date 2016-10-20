// =====================================================================================
//
//       Filename:  CombMcbFit
//
//    Description:  Simultaneous fit of KK and uu events with the same signal shape
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

#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>
extern void InitGui();
VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
TROOT root("polarimeter","polarimeter", initfuncs);


#include <TFile.h>

#include "fit.h"

int main(int argc,  char ** argv)
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
  //if(argc<2) return 1;
	TFile file(tree_file.c_str());
	TH1 * hisKK = (TH1*)file.Get("hMrecKK");
	TH1 * hisUU = (TH1*)file.Get("hMrecUU");
	//TTree * treeKK = (TTree*)file.Get("eventKK");
	//TTree * treeUU = (TTree*)file.Get("eventUU");
	TApplication theApp("root_app", &argc, argv);
  std::list<TH1*> hlst;
  hlst.push_back(hisKK);
  hlst.push_back(hisUU);
  fit(hisKK,hisUU);
  //fit(hlst);
	theApp.Run();
}
