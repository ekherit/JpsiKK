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

#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>
extern void InitGui();
VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
TROOT root("kksel","kksel", initfuncs);

#include "libFit.h"

#include "mctopo/mctopo.h"

#include "hashlist.h"

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

//parameters for shift and scale the Mrec (recoil mass)
const double MSCALE=1000; //will be MeV unis
const double MJPSI_SHIFT=3.097; //shift for JPsi mass
//number of bins in histogram
const int NBINS_KK= 200; //for Kaon events
const int NBINS_UU= 500; //for Muon events

double MRANGE=0.09; //Cut range for the recoil mass
double MAX_RECOIL_MASS=MJPSI_SHIFT+MRANGE*0.5;
double MIN_RECOIL_MASS=MJPSI_SHIFT-MRANGE*0.5;

double mshift(double m) 
{
  return MSCALE*(m-MJPSI_SHIFT);
}

enum {KAON, MUON};

enum {NEUTRAL=0, NEGATIVE=-1, POSITIVE=+1, POSITIVE_AND_NEGATIVE=2};
const int POSNEG=2;

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



TTree *  make_tree(TTree * parent_tree, string  name,  string  title)
{
  TTree * tree  = parent_tree->CloneTree(0);
  tree->SetName(name.c_str());
  tree->SetTitle(title.c_str());
  return tree;
};

TH1D *  make_hMrec(string name, string  title, int Nbins) 
{
  return new TH1D(name.c_str(),title.c_str(),Nbins,mshift(MIN_RECOIL_MASS),mshift(MAX_RECOIL_MASS));
};


/*
 * =====================================================================================
 *        Class:  Result_t
 *  Description:  Contain the result of the selection: histogram, data and monte carlo
 *  event tree
 * =====================================================================================
 */
struct Result_t
{
  TTree * event_tree;  //main event tree
  TTree * mctopo_tree; //corresponding mc truth information
  TTree * mis_tree; // missed data
  TH1D   * hMrec;      //histogram with pipi recoil mass
  Index_t index;       //the navigation complex index (Id of the data)
  Long64_t count;      //number of events

  Result_t(void)
  {
    event_tree = nullptr;
    mctopo_tree = nullptr;
    mis_tree = nullptr;
    hMrec = nullptr;
    index = {0,0,0};
    count =0;
  }

  Result_t & operator=(const Result_t & r)
  {
    event_tree=r.event_tree;
    mctopo_tree=r.mctopo_tree;
    mis_tree=r.mis_tree;
    hMrec = r.hMrec;
    index = r.index;
    count = r.count;
    return *this;
  }

  Result_t(Index_t idx, TTree * event, TTree * mctopo, RootMdc * mdc)
  {
    index = idx;
    count = 0;
    std::string chan;
    std::string sign;
    std::string ntrk;
    std::string title;
    int Nbins=0;
    switch(idx.channel)
    {
      case KAON:
        chan="K";
        Nbins = NBINS_KK;
        break;
      case MUON:
        chan="U";
        Nbins = NBINS_UU;
        break;
    }
    switch(idx.tracks)
    {
      case 3:
        ntrk = "3";
        break;
      case 4:
        ntrk = "4";
        break;
      case 3+4:
        ntrk = "34";
        break;
    }
    switch(idx.charge)
    {
      case -1:
        sign = "m";
        title = "negative " + chan + "  events with " + ntrk + " tracks";
        break;
      case 0:
        sign = chan;
        ntrk = "";
        title = chan+chan+ " channel";
        break;
      case 1:
        sign = "p";
        title = "positive " + chan + "  events with " + ntrk + " tracks";
        break;
      case 2:
        sign = "";
        title = "positive and negative " + chan + "  events with " + ntrk + " tracks";
        break;
    }
    std::string suffix = chan + sign + ntrk;
    std::string his_name = "hMrec" + suffix;
    std::string tree_name = "event" + suffix;
    std::string mctopo_tree_name = "mctopo"+ suffix;
    std::string mis_tree_name = "mis"+ suffix;


    hMrec = make_hMrec(his_name,"#pi^{+}#pi^{-} recoil mass for " + title, Nbins);
    event_tree  = make_tree(event,tree_name,"events for " + title);
    mctopo_tree = make_tree(mctopo,mctopo_tree_name,"Monte Carlo events for " + title);
    mis_tree = new TTree(mis_tree_name.c_str(),("Missing data" + title).c_str());
    mis_tree->Branch("Eemc",&(mdc->E),"Eemc[4]/D");
    std::cout << "Init result item: " << boost::format("(%-1d,%2d,%2d) %-4s") % index.channel % index.charge % index.tracks % suffix << std::endl;
  }
  void Fill(double Mrec, int run =0,RootMdc * mdc=nullptr)
  {
    event_tree->Fill();
    if(run<0) mctopo_tree->Fill();
    hMrec->Fill(mshift(Mrec));
    mis_tree->Fill();
    count++;
  };

  void Write(void)
  {
    hMrec->Write();
    event_tree->Write();
    mctopo_tree->Write();
    mis_tree->Write();
  }

};


int main(int argc, char ** argv)
{
  std::ios_base::sync_with_stdio(false);
  namespace po=boost::program_options;
  po::options_description opt_desc("Allowed options");
  std::string tree_name;
  std::vector< std::string> files = {"sample.root"};
  std::string output_file;
  std::string event_tree_name;
  std::string mctopo_tree_name;
  std::string tree_suffix;
  unsigned long long NMAX;
  unsigned long MAX_EVERY;
  //selection parameters
  double PID_CHI2;
  double KIN_CHI2;
  double KIN_CHI2_1C; //kin chi2 for 1C kinematic fit
  opt_desc.add_options()
    ("help,h","Print this help")
    ("input", po::value<std::vector< std::string> >(&files), "Root file (.root) with the data")
    ("output,o", po::value<std::string>(&output_file)->default_value("result.root"), "Output file")
    ("event_tree",  po::value<std::string>(&event_tree_name)->default_value("event"), "Event tree name")
    ("mctopo_tree", po::value<std::string>(&mctopo_tree_name)->default_value("mctopo"), "mctopo tree name")
    ("tree_suffix", po::value<std::string>(&tree_suffix)->default_value(""), "suffix for the tree")
    ("max_every", po::value<unsigned long>(&MAX_EVERY)->default_value(1e4), "Maximum every for printing")
    ("N", po::value<unsigned long long>(&NMAX)->default_value(std::numeric_limits<unsigned long long>::max()), "Maximum event number to proceed")
    ("mrange", po::value<double>(&MRANGE)->default_value(0.09), "Pion recoil mass range")
    ("kin_chi2", po::value<double>(&KIN_CHI2)->default_value(40), "Kinematic chi square cut")
    ("kin_chi2_1c", po::value<double>(&KIN_CHI2_1C)->default_value(10), "Kinematic chi square cut for 3 particles fit")
    ("pid_chi2", po::value<double>(&PID_CHI2)->default_value(20), "Particle id cut")
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

  // calculate the maximum and minimum recoil mass
  MAX_RECOIL_MASS=MJPSI_SHIFT+MRANGE*0.5;
  MIN_RECOIL_MASS=MJPSI_SHIFT-MRANGE*0.5;

  //read input files
  std::list < std::string > file_list;
  for( auto dir : files)
  {
    file_list.merge(get_file_list_in_dir(dir));
  };


  //load the data
  RootEvent  event(load_tree(event_tree_name+tree_suffix,file_list));
  RootMCTopo mctopo(load_tree(mctopo_tree_name+tree_suffix,file_list));
  RootMdc    mdc(load_tree("mdc"+tree_suffix,file_list));

  TFile file(output_file.c_str(),"RECREATE");


  //here the result will be stored
  std::unordered_map<Index_t, Result_t, IndexHash_t> R; 

  //define helper function for initialization
  auto InitResultItem = [&R,&event, &mctopo, &mdc](Index_t idx)
  {
    R[idx] = Result_t(idx,event.fChain,mctopo.fChain, &mdc);
  };

  //initialize the result data
  std::cout << "Init result items map"<< std::endl;
  for(int  chan : {KAON,MUON}) //loop over channel
  {
    for(int sign : {0,-1,1,2}) //over charge
    {
      /* 
       *  0    : is for 4C kenematic fit with totaly 4 tracks
       *  1,-1 : 1C kinematic for 3 tracks, charge corresponds the sign of the system of two tracks
       *  2    : 1C kinematic fit for 3 tracks, but here both charges of system put together
       */ 
      for(int ntrk : {3,4,3+4}) //over tracks number
      {
        if(sign == 0 && (ntrk ==3 || ntrk==(4+3))) 
        {
          continue;
        }
        InitResultItem({chan,sign,ntrk});
      }
    }
  }

  Long64_t nentries = event.fChain->GetEntriesFast();

  std::cout << "Total number of events in the data files: " << nentries << std::endl;

  Long64_t nbytes = 0, nb = 0;
  
  //main loop over events
  for (Long64_t jentry=0; jentry<nentries && jentry<NMAX;jentry++)
  {
    Long64_t ientry;
    ientry = event.LoadTree(jentry);
    if (ientry < 0)  break;
    if(event.run < 0) 
    {
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
    unsigned long hash=0;
    if(event.run<0)  
    {
      hash = mctopo_hash(&mctopo);
      mctopo.hash = hash;
    }

    //selection options
    bool RecoilCut = MIN_RECOIL_MASS <= event.Mrec && event.Mrec <= MAX_RECOIL_MASS;
    bool KinematicCut = event.kin_chi2 <= KIN_CHI2;
    bool PidCut = event.pid_chi2 <= PID_CHI2;
    bool ThetaCut =  fabs(cos(event.theta[2])) < 0.8 && fabs(cos(event.theta[3])) <0.8;
    bool Common1CCut   = 
         event.ngntrack == 0  
      && (event.ngtrack == 3 || event.ngtrack == 4)
      && event.kin_chi2 < KIN_CHI2_1C;

    bool EpCut[4] = {false,false,false,false};
    for(int i=2;i<4;i++)
    {
      double ep = mdc.E[i]/mdc.p[i];
      EpCut[i] =  0.1 < ep && ep < 0.25;
    };

    bool Muon4TracksEpCut = EpCut[2] && EpCut[3]; 
    bool Muon3TracksEpCut = (event.sign < 0 && EpCut[2]) || (event.sign > 0 && EpCut[3]);

    auto fill_data = [&R,&event,&mdc] (int s, int n)
    {
      R[{event.channel, s, n}].Fill(event.Mrec, event.run, &mdc);
    };

    auto posneg_fill = [&](void) 
    {
      fill_data(event.sign, event.ngtrack);
      fill_data(POSNEG,3+4);
      fill_data(POSNEG,event.ngtrack);
    };

    if(    RecoilCut 
        && KinematicCut 
        && PidCut 
        && ThetaCut
      )
    {

      //Case of Muon channel with 4 tracks with 4C kinematic fit.
      //Additional cut for the muons will be applied
      //see EpCut above    0.1  < E/p < 0.25
      //this will suppress half of the pions
      if(event.uu == 1 && Muon4TracksEpCut)
      {
        fill_data(event.sign, event.ngtrack);
      }

      //Kaon 4C kinematic case with 4 tracks
      if(event.KK == 1) 
      {
        fill_data(event.sign, event.ngtrack);
      }

      //1C kinematic fit (using 3tracks only) but in real has 3 or 4 tracks 
      if(Common1CCut)
      {
        if (event.u == 1 && Muon3TracksEpCut)
        {
          posneg_fill();
        }
        if(event.K == 1)
        {
          posneg_fill();
        }
      }
    }

    //print progress
    //logarifmicaly suppress output
    //at first print every, then every 10, then every 100... until every 10000
    static unsigned long every = 1;
    static unsigned long head = 0;
    if(jentry > every*10 && every<MAX_EVERY) every*=10;
    if(jentry % every==0 || jentry==nentries-1)
    {
      if(head++%20 == 0)
      {
        std::cout << boost::format("#%14s %8s %10s %10s %10s %10s %10s %10s %10s %10s")
          % "entry"
          % "run"
          % "NKK"
          % "NK4"
          % "NK3"
          % "NK3+4"
          % "NUU"
          % "NU4"
          % "NU3"
          % "NU3+4";
        std::cout << endl;
      }
      std::cout << boost::format("%15d %8d %10d %10d %10d %10d %10d %10d %10d %10d") % jentry % event.run 
        %  R[{KAON,0,4}].count 
        %  R[{KAON,POSNEG,4}].count 
        %  R[{KAON,POSNEG,3}].count 
        %  R[{KAON,POSNEG,3+4}].count 
        %  R[{MUON,0,4}].count 
        %  R[{MUON,POSNEG,4}].count 
        %  R[{MUON,POSNEG,3}].count 
        %  R[{MUON,POSNEG,3+4}].count;
      if(event.run<0)
      {
        std::cout << "  hash = " << setw(10) << std::hex << hash << "    " << std::dec << mctopo_info(&mctopo);
      }
      std::cout << std::endl;
    }
  }

  //save result into root file
  for(auto & r : R)
  {
    r.second.Write();
  }
  

  //prepare usefull cuts 
  std::map<std::string, TCut> UsefulCuts;
  auto do_cut = [&UsefulCuts](std::string name)
  {
    UsefulCuts[name] = make_cut({name});
    UsefulCuts[name].SetName(name.c_str());
  };

  do_cut("KK");
  do_cut("UU");
  do_cut("jpsi");
  do_cut("nojpsi");

  //save those useful cuts into root file
  for(auto & cut : UsefulCuts)
  {
    cut.second.Write();
  }

  //close file
  file.Close();

}
