/*
 * =====================================================================================
 *
 *       Filename:  JpsiKK.cxx
 *
 *    Description:  Event selection for measurement branching fraction Jpsi->K+K-
 *    via Psi(2S)->Jpsi(->K+K-)pi+pi- decay
 *
 *        Version:  1.0
 *        Created:  2014-11-06
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physics
 *
 * =====================================================================================
 */
#include <stdexcept>

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "VertexFit/IVertexDbSvc.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/ISvcLocator.h"

#include "EventModel/EventModel.h"
#include "EventModel/Event.h"

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
//#include "DstEvent/TofHitStatus.h"
#include "EventModel/EventHeader.h"

#include "McTruth/McParticle.h"
//#include "EventNavigator/EventNavigator.h"


#include "TMath.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif
#include "JpsiKK.h"
#include "Sphericity.h"




#include <algorithm>
#include <limits>



#include "Defs.h"
#include "PhysConst.h"
#include "Utils.h"
#include "MyPid.h"
#include "SelectionHelper.h"

#include "KinematicFit.h"


const  double XMASS[5] = {KAON_MASS, MUON_MASS, ELECTRON_MASS, PION_MASS, PROTON_MASS};
double BEAM_CENTER_MASS_ENERGY = PSIP_MASS;

JpsiKK::JpsiKK(const std::string& name, ISvcLocator* pSvcLocator) :
  Algorithm(name, pSvcLocator)
{
  declareProperty("CENTER_MASS_ENERGY", cfg.CENTER_MASS_ENERGY = 0); //GeV
  declareProperty("MIN_CHARGED_TRACKS", cfg.MIN_CHARGED_TRACKS=3); //minimum number of charged tracks in selection
  declareProperty("MAX_CHARGED_TRACKS", cfg.MAX_CHARGED_TRACKS=4); //maximum number of charged tracks in selection
  declareProperty("MAX_NEUTRAL_TRACKS", cfg.MAX_NEUTRAL_TRACKS=1000); //maximum number of good charged tracks in selection

  //good charged track configuration
  declareProperty("IP_MAX_Z",      cfg.IP_MAX_Z = 10.0); //cm
  declareProperty("IP_MAX_RHO",    cfg.IP_MAX_RHO = 1.0); //cm
  declareProperty("MAX_COS_THETA", cfg.MAX_COS_THETA = 0.93);

  //good neutral track configuration
  //endcup calorimeter
  declareProperty("EMC_ENDCUP_MIN_COS_THETA", cfg.EMC_ENDCUP_MIN_COS_THETA = 0.86);
  declareProperty("EMC_ENDCUP_MAX_COS_THETA", cfg.EMC_ENDCUP_MAX_COS_THETA = 0.92);
  declareProperty("EMC_ENDCUP_MIN_ENERGY",    cfg.EMC_ENDCUP_MIN_ENERGY = 0.05);
  //barrel calorimeter
  declareProperty("EMC_BARREL_MAX_COS_THETA", cfg.EMC_BARREL_MAX_COS_THETA = 0.8);
  declareProperty("EMC_BARREL_MIN_ENERGY",    cfg.EMC_BARREL_MIN_ENERGY = 0.025);

  declareProperty("NEUTRAL_CLOSE_CHARGED_ANGLE",    cfg.NEUTRAL_CLOSE_CHARGED_ANGLE = 10);


  declareProperty("MAX_MUON_EP_RATIO", cfg.MAX_MUON_EP_RATIO = 0.26);
  declareProperty("MIN_MUON_EP_RATIO", cfg.MIN_MUON_EP_RATIO = 0);

  declareProperty("MAX_KAON_EP_RATIO", cfg.MAX_KAON_EP_RATIO = 0.8);
  declareProperty("MIN_KAON_EP_RATIO", cfg.MIN_KAON_EP_RATIO = 0);

  declareProperty("MAX_PION_MOMENTUM", cfg.MAX_PION_MOMENTUM = 0.45); //GeV
  declareProperty("MIN_PION_MOMENTUM", cfg.MIN_PION_MOMENTUM = 0); //GeV

  declareProperty("MIN_KAON_MOMENTUM", cfg.MIN_KAON_MOMENTUM = 1.0); //GeV
  declareProperty("MAX_KAON_MOMENTUM", cfg.MAX_KAON_MOMENTUM = 2.0); //GeV

  declareProperty("MIN_MUON_MOMENTUM", cfg.MIN_MUON_MOMENTUM = 1.0); //GeV
  declareProperty("MAX_MUON_MOMENTUM", cfg.MAX_MUON_MOMENTUM = 2.0); //GeV

  declareProperty("MIN_RECOIL_MASS", cfg.MIN_RECOIL_MASS = 3.0); //GeV
  declareProperty("MAX_RECOIL_MASS", cfg.MAX_RECOIL_MASS = 3.2); //GeV
  declareProperty("MIN_INVARIANT_MASS", cfg.MIN_INVARIANT_MASS = 3.0); //GeV
  declareProperty("MAX_INVARIANT_MASS", cfg.MAX_INVARIANT_MASS = 3.2); //GeV
  declareProperty("MIN_KAON_MISSING_MASS", cfg.MIN_KAON_MISSING_MASS = 0.1); //GeV^2
  declareProperty("MAX_KAON_MISSING_MASS", cfg.MAX_KAON_MISSING_MASS = 0.6); //GeV^2
  declareProperty("MIN_MUON_MISSING_MASS", cfg.MIN_MUON_MISSING_MASS = 0); //GeV^2
  declareProperty("MAX_MUON_MISSING_MASS", cfg.MAX_MUON_MISSING_MASS = 0.1); //GeV^2


  declareProperty("MIN_MISSING_MASS", cfg.MIN_MISSING_MASS = -0.1); //GeV^2
  declareProperty("MAX_MISSING_MASS", cfg.MAX_MISSING_MASS = +0.1); //GeV^2

  declareProperty("MAX_KIN_CHI2", cfg.MAX_KIN_CHI2 = 40); //GeV^2
  declareProperty("MAX_PID_CHI2", cfg.MAX_PID_CHI2 = 40); //GeV^2

  declareProperty("PASS_KIN_PID_CUT", cfg.PASS_KIN_PID_CUT = false); //GeV^2

  declareProperty("FILL_MUC", cfg.FILL_MUC = true); 
  declareProperty("FILL_MDC", cfg.FILL_MDC = true); 
  declareProperty("FILL_EMC", cfg.FILL_EMC = true); 
  declareProperty("FILL_TOF", cfg.FILL_TOF = true); 
  declareProperty("FILL_DEDX", cfg.FILL_DEDX = true); 
  declareProperty("FILL_NEUTRAL", cfg.FILL_NEUTRAL = true); 
}


StatusCode JpsiKK::initialize(void)
{
  MsgStream log(msgSvc(), name());

  log << MSG::INFO << "in initialize()" << endmsg;
  event_proceed=0;
  event_write = 0;
  theCounter.clear();
	if(cfg.CENTER_MASS_ENERGY == 0) cfg.CENTER_MASS_ENERGY = PSIP_MASS;
	BEAM_CENTER_MASS_ENERGY = cfg.CENTER_MASS_ENERGY;
  std::cout << "JpsiKK.CENTER_MASS_ENERGY=" << BEAM_CENTER_MASS_ENERGY << std::endl;
  

	try
	{
		init_tuple(fEvent,      "FILE1/event",  "Signal events pi+pi- K+K-, or pi+pi- mu+mu-");
		init_tuple(fMdc,        "FILE1/mdc",    "Main Drift Chamber");
		init_tuple(fDedx,       "FILE1/dedx",   "Dedx info for signal");
		init_tuple(fEmc,        "FILE1/emc",    "Electromagnetic Calorimeter");
		init_tuple(fTof,        "FILE1/tof",    "Time of Flight");
		init_tuple(fNeutral,    "FILE1/neutral","Good neutral tracks");
		init_tuple(fMC,         "FILE1/mc",     "Monte Carlo truth information");
		init_tuple(fMCTopo,     "FILE1/mctopo", "Monte Carlo truth information topology");
		init_tuple(fMuc,        "FILE1/muc",    "MUon chamber information");
		//init_tuple(this, fPid,        "FILE1/pid","particle id");
	}
	catch(std::runtime_error & error)
	{
		log << MSG::ERROR << error.what() << endmsg;
		return StatusCode::FAILURE;
	}
  return StatusCode::SUCCESS;
}

StatusCode JpsiKK::execute()
{
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "executing" << endreq;

  SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
  static unsigned long  every = 1;
  if(event_proceed > every*10 && every<1e3) every*=10;
  if(event_proceed % every == 0 )
  {
		static unsigned long nprints  = 0;
    printSelectionDigest(nprints % 20 == 0);
		nprints ++;
  }
  event_proceed++;

  //  Get information about reconstructed events
  SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
  SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol);


  //fill initial value of the selected event

  good_charged_tracks=createGoodChargedTrackList(cfg, evtRecEvent, evtRecTrkCol);
  good_neutral_tracks=createGoodNeutralTrackList2(cfg, evtRecEvent, evtRecTrkCol);

  //SELECTION CODE
  if( cfg.MAX_NEUTRAL_TRACKS < good_neutral_tracks.size()) return StatusCode::SUCCESS;
  //SELECTION CODE
  if( good_charged_tracks.size() < cfg.MIN_CHARGED_TRACKS || cfg.MAX_CHARGED_TRACKS < good_charged_tracks.size()) return StatusCode::SUCCESS;

  
  TrackList_t charged_tracks; //selected tracks for specific cut
  TrackList_t positive_charged_tracks; //selected tracks for specific cut
  TrackList_t negative_charged_tracks; //selected tracks for specific cut
  TrackList_t positive_pion_tracks; //selected pion tracks for specific cut
  TrackList_t negative_pion_tracks; //selected pion tracks for specific cut
  TrackList_t other_positive_tracks; //other positive tracks for specific cut
  TrackList_t other_negative_tracks; //other positive tracks for specific cut
  for(TrackList_t::iterator track=good_charged_tracks.begin(); track!=good_charged_tracks.end(); track++)
  {
    //EvtRecTrackIterator & itTrk = *track;
    EvtRecTrack *Trk = *track;
    if(!Trk->isMdcTrackValid()) continue; 
    RecMdcTrack *mdcTrk = Trk->mdcTrack();
    double c = fabs(cos(mdcTrk->theta()));
    double p = mdcTrk->p();
    double q = mdcTrk->charge();
    bool barrel = c < cfg.EMC_BARREL_MAX_COS_THETA;
    if(barrel) 
    {
			//fill charged tracks list. It will not be used in selection
      charged_tracks.push_back(Trk);
			//fill positive and negative chargned track list not used in selection
			if(q>0) positive_charged_tracks.push_back(Trk);
			if(q<0) negative_charged_tracks.push_back(Trk);

			//preselect pion candidates
			if(in(p, cfg.MIN_PION_MOMENTUM, cfg.MAX_PION_MOMENTUM)) 
			{
				if(q>0) positive_pion_tracks.push_back(Trk);
				if(q<0) negative_pion_tracks.push_back(Trk);
			}
			//preselect muon and kaon candidates
			if(in(p, std::min(cfg.MIN_KAON_MOMENTUM, cfg.MIN_MUON_MOMENTUM), std::max(cfg.MAX_KAON_MOMENTUM, cfg.MAX_MUON_MOMENTUM)))
			{
				if(Trk->isEmcShowerValid())
				{
					if(q>0) other_positive_tracks.push_back(Trk);
					if(q<0) other_negative_tracks.push_back(Trk);
				}
			}
    }
  }


  //create pion pairs
  TrackPairList_t pion_pairs;
  for(TrackList_t::iterator i=negative_pion_tracks.begin(); i!=negative_pion_tracks.end(); ++i)
    for(TrackList_t::iterator j=positive_pion_tracks.begin(); j!=positive_pion_tracks.end(); ++j)
    {
      TrackPair_t pair(*i,*j);
      if( in(getPionRecoilMass(*i, *j),  cfg.MIN_RECOIL_MASS, cfg.MAX_RECOIL_MASS)) 
      {
        pion_pairs.push_back(pair);
      }
    }
  //SELECTION CODE we must have at list one pion pair
  if(pion_pairs.empty()) return StatusCode::SUCCESS; //we must have at list one pion pair
  TrackPair_t pion_pair = pion_pairs.front();
	if(pion_pairs.size()>1)
	{
		//find the best pion pair using closest value to JPSI_MASS
		for(TrackPairList_t::iterator p=pion_pairs.begin();p!=pion_pairs.end();p++)
		{
			if(fabs(getPionRecoilMass(p->first,  p->second)- JPSI_MASS) <  fabs(getPionRecoilMass(pion_pair.first,  pion_pair.second) - JPSI_MASS)) pion_pair = *p;
		}
	}
	//now we have one pion pair candidate

  //SELECTION CODE
  //if(positive_charged_tracks.size()!=2 || negative_charged_tracks.size()!=2) return StatusCode::SUCCESS;

  //if no other particles
	if(other_negative_tracks.empty() && other_positive_tracks.empty()) return StatusCode::SUCCESS;

  bool pos = !other_positive_tracks.empty();
  bool neg = !other_negative_tracks.empty();



  std::list<SelectionHelper_t> sh_list;
  if(neg)
  {
    sh_list.push_back(SelectionHelper_t(cfg));
    SelectionHelper_t & sh = sh_list.back();
		sh.kinfit(pion_pair, other_negative_tracks.front());
		sh.totalPass(cfg.PASS_KIN_PID_CUT);

  }
  if(pos)
  {
    sh_list.push_back(SelectionHelper_t(cfg));
    SelectionHelper_t & sh = sh_list.back();
		sh.kinfit(pion_pair, other_positive_tracks.front());
		sh.totalPass(cfg.PASS_KIN_PID_CUT);
  }
  if(pos && neg)
  {
    sh_list.push_back(SelectionHelper_t(cfg));
    SelectionHelper_t & sh = sh_list.back();
    sh.kinfit(pion_pair, other_negative_tracks.front(), other_positive_tracks.front());
		sh.totalPass(cfg.PASS_KIN_PID_CUT);
  }

	std::list<int> pid_list;
	pid_list.push_back(ID_KAON);
	pid_list.push_back(ID_MUON);
  for(std::list<SelectionHelper_t>::iterator it = sh_list.begin(); it!=sh_list.end(); it++)
  {
    SelectionHelper_t & sh = *it;
    for(std::list<int>::iterator chan = pid_list.begin(); chan!=pid_list.end() ; chan++)
    {
      bool pass  = sh.pass && sh.channel == *chan;
      if(pass)
      {
        std::vector<HepLorentzVector> Pkf = sh.getMomentum(*chan);
        TrackVector_t Tracks = sh.tracks; 
        fEvent.nkinbg=2;
        //kinematic fit for suppress Ψ(2S) → π+π-(J/Ψ → π-(ρ(770)+ → (π0 → ɣɣ)π+)) = ɣɣπ+π+π-π- (final state)
        fEvent.kin_chi2_bg[0] = sh.kinfit_4pi_2g(good_neutral_tracks);
        fEvent.kin_chi2_bg[1] = sh.kinfit_4pi_2g_missed_gamma(good_neutral_tracks);
        int charge=0;
        if(Tracks.size()==3) 
        {
          charge = getMdcTrack(Tracks[2])->charge();
          Tracks.push_back(0); //workout 3trk case
        }
        if(charge<0 && pos)
        {
          Tracks[3] = other_positive_tracks.front();
        }
        if(charge>0)
        {
          std::swap(Pkf[2],  Pkf[3]);
          std::swap(Tracks[2],  Tracks[3]);
          if(neg)
          {
            Tracks[2] = other_negative_tracks.front();
          }
        }

        fEvent.channel=*chan;
        fEvent.sign=charge;
        fEvent.K=0;
        fEvent.KK=0;
        fEvent.u=0;
        fEvent.uu=0;
        switch(*chan)
        {
          case ID_KAON:
            if(charge==0) 
            {
              fEvent.KK=1;
            }
            else 
            {
              fEvent.K = 1;
            }
            break;
          case ID_MUON:
            if(charge==0) 
            {
              fEvent.uu=1;
            }
            else 
            {
              fEvent.u = 1;
            }
            break;
        }

        theCounter[*chan][charge][good_charged_tracks.size()]++;

        fEvent.kin_chi2 = sh.getKinChi2(*chan);
        fEvent.pid_chi2 = sh.getPidChi2(*chan);
        fEvent.pid_prob = sh.prob[*chan];
        fEvent.mypid_chi2 = sh.mypid_chi2[*chan];

        for(int pid=0;pid<5;pid++)
        {
          fEvent.kchi[pid] = sh.getKinChi2(pid);
          fEvent.pchi[pid] = sh.getPidChi2(pid);
          fEvent.prob[pid] = sh.prob[pid];
          fEvent.mypchi[pid] = sh.mypid_chi2[pid];
          HepLorentzVector p[4];
          for(int n = 2; n<4; n++)
          {
            if(Tracks[n]!=0) 
            {
              RecMdcKalTrack *mdcTrk = Tracks[n]->mdcKalTrack();
              p[n] =  mdcTrk->p4(XMASS[pid]);
            }
            else
            {
              p[n] = sh.KF[pid].P[n];
            }
          }
          fEvent.kM23[pid] = (p[2]+p[3]).m();
        }

        fEvent.run=eventHeader->runNumber();
        fEvent.event=eventHeader->eventNumber();
        fEvent.time=eventHeader->time();
        fEvent.npid = 5;

        fEvent.ngood_charged_track = good_charged_tracks.size();
        fEvent.ngood_neutral_track = good_neutral_tracks.size();
        fEvent.npositive_track = positive_charged_tracks.size();
        fEvent.nnegative_track = negative_charged_tracks.size();
        fEvent.npositive_pions = positive_pion_tracks.size();
        fEvent.nnegative_pions = negative_pion_tracks.size();
        fEvent.npion_pairs = pion_pairs.size();
        try
        {
          fillTuples(Pkf, Tracks, evtRecEvent, evtRecTrkCol);
          writeTuples();
        }
        catch(std::runtime_error & error)
        {
          log << MSG::ERROR  << error.what() << endmsg;
          return StatusCode::FAILURE;
        }

      }
    }
  }
  return StatusCode::SUCCESS;
}



StatusCode JpsiKK::finalize()
{
  std::cout << "============== Selection Result =====================" << std::endl;
  printSelectionDigest(true);
  std::cout.precision(6);
  std::cout.setf( std::ios::fixed, std:: ios::floatfield );
  std::cout << "main selection efficiency:" << std::endl;
  std::cout << "    ε4C(K+K-)  = " << double(theCounter[ID_KAON][0][4])/event_proceed << std::endl;
  std::cout << "    ε4C(u+u-) = " << double(theCounter[ID_MUON][0][4])/event_proceed << std::endl;
  std::cout << "    ε4C(K+K-)/ε4C(u+u-) = " << double(theCounter[ID_KAON][0][4])/double(theCounter[ID_MUON][0][4]) << std::endl;
  std::cout << "tracking efficiency: ε(+track) = N4(-track)/(N3(-track)+N4(-track))" << std::endl;
  std::map< int , std::map < int, double> > eps_trk; //tracking efficiency
  std::map< int , std::map < int, double> > eps_trk_error; //tracking efficiency

  for(int ch = 0; ch<2; ch++)
  {
    int q;
    q = -1;
    eps_trk[ch][q]  =  double(theCounter[ch][-q][4])/double(theCounter[ch][-q][3]  + theCounter[ch][-q][4]);
    eps_trk_error[ch][q] = sqrt(eps_trk[ch][q]*(1.0 - eps_trk[ch][q])/double(theCounter[ch][-q][3]  + theCounter[ch][-q][4]));

    q = +1;
    eps_trk[ch][q]  =  double(theCounter[ch][-q][4])/double(theCounter[ch][-q][3]  + theCounter[ch][-q][4]);
    eps_trk_error[ch][q] = sqrt(eps_trk[ch][q]*(1.0 - eps_trk[ch][q])/double(theCounter[ch][-q][3]  + theCounter[ch][-q][4]));

    q = 0;
    eps_trk[ch][q]  =  double(theCounter[ch][-1][4] + theCounter[ch][+1][4] )/double(theCounter[ch][-1][3]  + theCounter[ch][-1][4] + theCounter[ch][+1][3]  + theCounter[ch][+1][4]);
    eps_trk_error[ch][q] = sqrt(eps_trk[ch][q]*(1.0 - eps_trk[ch][q])/double(theCounter[ch][-1][3]  + theCounter[ch][-1][4] + theCounter[ch][+1][3]  + theCounter[ch][+1][4]));
  }

  std::cout << "  ε(K-)     = "  <<  eps_trk[ID_KAON][-1] << " ± " << eps_trk_error[ID_KAON][-1]<< std::endl;
  std::cout << "  ε(K+)     = "  <<  eps_trk[ID_KAON][+1] << " ± " << eps_trk_error[ID_KAON][+1]<< std::endl;
  std::cout << "  ε(u-)     = "  <<  eps_trk[ID_MUON][-1] << " ± "  << eps_trk_error[ID_MUON][-1]<< std::endl;
  std::cout << "  ε(u+)     = "  <<  eps_trk[ID_MUON][+1] << " ± " << eps_trk_error[ID_MUON][+1]<< std::endl;
  std::cout << "  ε(K)      = "  <<  eps_trk[ID_KAON][0] << " ± " << eps_trk_error[ID_KAON][0]<< std::endl;
  std::cout << "  ε(u)      = "  <<  eps_trk[ID_MUON][0] << " ± " << eps_trk_error[ID_MUON][0]<< std::endl;
  std::cout << "  ε(K)/ε(u) = "  <<  eps_trk[ID_KAON][0]/eps_trk[ID_MUON][0] << " ± " <<  eps_trk[ID_KAON][0]/eps_trk[ID_MUON][0] * sqrt(sq(eps_trk_error[ID_KAON][0]/eps_trk[ID_KAON][0]) + sq(eps_trk_error[ID_MUON][0]/eps_trk[ID_MUON][0]))<< std::endl;
  return StatusCode::SUCCESS;
}

template <class A>
void JpsiKK::init_tuple(A & a,  const char * dir, const char * title)
{
  NTuplePtr nt(this->ntupleSvc(), dir);
  if(nt) a.tuple = nt;
  else
  {
    a.tuple = this->ntupleSvc()->book(dir, CLID_ColumnWiseTuple, title);
    if(a.tuple)
    {
      a.init_tuple();
    }
    else
    {
			char buf[1024];
			sprintf(buf, "   Cannot book N-tuple: %d",  long(a.tuple));
			throw std::runtime_error(buf);
    }
  }
}

void  JpsiKK::fillTuples(const std::vector<CLHEP::HepLorentzVector> & Pkf,  TrackVector_t & Tracks,
		SmartDataPtr<EvtRecEvent>    & evtRecEvent, 
		SmartDataPtr<EvtRecTrackCol> & evtRecTrkCol
    )
{
  SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(),  EventModel::MC::McParticleCol);
	fEvent.fill(Pkf);
  //fEvent.Mpi0 = getPi0Mass(good_neutral_tracks);
	//fPid.ntrack=4;
	fMdc.T.ntrack=4;
	fDedx.ntrack=4;
	fEmc.ntrack=4;
	fTof.ntrack=4;
  fMuc.ntrack=4;
	//fMdc.M.Mrec = get_recoil__mass(Tracks[0], Tracks[1], PION_MASS,  cfg.CENTER_MASS_ENERGY);
	for(int i=0;i<4;i++)
	{
		if(Tracks[i]==0) continue;
    fEvent.fill(i,Tracks[i]);
    if(cfg.FILL_EMC)  fEmc.fill (i,  Tracks[i], evtRecEvent, evtRecTrkCol);
		if(cfg.FILL_MDC)  fMdc.fill (i,  Tracks[i]);
		if(cfg.FILL_DEDX) fDedx.fill(i,  Tracks[i]);
		if(cfg.FILL_TOF)  fTof.fill (i,  Tracks[i]);
    if(cfg.FILL_MUC)  fMuc.fill (i,  Tracks[i]);
	}
	fMdc.fill_mass(Tracks, Pkf);
	//Monte Carlo information
	if(fEvent.run<0)
	{
		
		if(!mcParticleCol)
		{
			throw std::runtime_error("Could not retrieve McParticelCol");
		}
		fMCTopo.fill(mcParticleCol);
		fMC.fill(Pkf, mcParticleCol);
	}
	fNeutral.fill(good_neutral_tracks,evtRecEvent, evtRecTrkCol);
}

void JpsiKK::writeTuples(void)
{
	if(fEvent.run<0)
	{
		fMC.write();
		fMCTopo.write();
	}
  fEvent.write();
  //fPid.write();
  if(cfg.FILL_MDC)  fMdc.write();
  if(cfg.FILL_EMC)  fEmc.write();
  if(cfg.FILL_DEDX) fDedx.write();
  if(cfg.FILL_TOF)  fTof.write();
  if(cfg.FILL_MUC)  fMuc.write();
  if(cfg.FILL_NEUTRAL) fNeutral.tuple->write();
  event_write++;
}

void JpsiKK::printSelectionDigest(bool head)
{
  int big_width=20;
  int width=15;
  if(head)
  {
    std::cout << setw(big_width) << "# event proceed";
    std::cout << " " << setw(big_width)    << "event written";
    std::cout << " " << setw(width) << "N4(K+K-)";
    std::cout << " " << setw(width) << "N3(K+)";
    std::cout << " " << setw(width) << "N4(K+)";
    std::cout << " " << setw(width) << "N3(K-)";
    std::cout << " " << setw(width) << "N4(K-)";
    std::cout << " " << setw(width) << "N4(u+u-)";
    std::cout << " " << setw(width) << "N3(u+)";
    std::cout << " " << setw(width) << "N4(u+)";
    std::cout << " " << setw(width) << "N3(u-)";
    std::cout << " " << setw(width) << "N4(u-)";
    std::cout << " " << setw(width) << "chi2(3pi)";
    std::cout << endl;
  }
  std::cout << setw(big_width)    << event_proceed;
  std::cout << " " << setw(big_width) << event_write;
  int Pids[2] = {ID_KAON,ID_MUON};
  for(int pid = 0; pid < 2; pid++)
  {
    std::cout << " " << setw(width) << theCounter[pid][0][4];
    std::cout << " " << setw(width) << theCounter[pid][+1][3];
    std::cout << " " << setw(width) << theCounter[pid][+1][4];
    std::cout << " " << setw(width) << theCounter[pid][-1][3];
    std::cout << " " << setw(width) << theCounter[pid][-1][4];
  }
  std::cout << " " << setw(width) << kin_chi2_3pi;
  std::cout << endl;
}


