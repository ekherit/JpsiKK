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


#include "ParticleID/ParticleID.h"


#include <algorithm>
#include <limits>

int JpsiKK::RootEmc::ARRAY_SIZE = 100;


#include "Defs.h"
#include "PhysConst.h"
#include "Utils.h"
#include "MyPid.h"
#include "SelectionHelper.h"

#include "KinematikFit.h"

enum
{
	OTHER_NO_TRACK = 0, 
	OTHER_NEGATIVE_TRACK=0x1, 
	OTHER_POSITIVE_TRACK=0x2, 
	OTHER_TWO_TRACKS=0x3
};

enum
{
	CHAN_KAONS = ID_KAON, 
	CHAN_MUONS = ID_MUON, 
	CHAN_KAON_MUON  = 10, 
	CHAN_MUON_KAON  = 11
};


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

}

//this is service function for fast book ntuple
template <class A>
StatusCode init_tuple(JpsiKK * alg, A & a,  const char * dir, const char * title, MsgStream & log)
{
  StatusCode status;
  NTuplePtr nt(alg->ntupleSvc(), dir);
  if(nt) a.tuple = nt;
  else
  {
    a.tuple = alg->ntupleSvc()->book(dir, CLID_ColumnWiseTuple, title);
    if(a.tuple)
    {
      return a.init_tuple();
    }
    else
    {
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(a.tuple) << endmsg;
      return StatusCode::FAILURE;
    }
  }
  return status;
}

StatusCode JpsiKK::initialize(void)
{
  MsgStream log(msgSvc(), name());

  log << MSG::INFO << "in initialize()" << endmsg;
  event_proceed=0;
  event_write = 0;
  event_with_kaons=0;
  event_with_muons=0;
  event_with_pions=0;
  event_with_electrons=0;
  event_with_protons=0;
  good_kinematic_fit=0;
	if(cfg.CENTER_MASS_ENERGY == 0) cfg.CENTER_MASS_ENERGY = PSIP_MASS;

  StatusCode status;
  status = init_tuple(this, fEvent,  "FILE1/event","Signal events pi+pi- K+K-, or pi+pi- mu+mu-",log);
//  status = init_tuple(this, fPid,    "FILE1/pid","particle id",log);
  status = init_tuple(this, fMdc,    "FILE1/mdc","Mdc info for signal",log);
  status = init_tuple(this, fDedx,   "FILE1/dedx","Dedx info for signal",log);
  status = init_tuple(this, fEmc,    "FILE1/emc","Emc info for signal",log);
  status = init_tuple(this, fTof,    "FILE1/tof","Tof info for signal",log);
  status = init_tuple(this, fNeutral,"FILE1/neutral","Good neutral tracks",log);
  status = init_tuple(this, fMC,     "FILE1/mc","Monte Carlo truth information",log);
  status = init_tuple(this, fMCTopo,     "FILE1/mctopo","Monte Carlo truth information topology",log);

  return status;
}


StatusCode JpsiKK::RootEvent::init_tuple(void)
{
  StatusCode status;
  status = tuple->addItem ("run", run); //run number
  status = tuple->addItem ("event", event); //event number
  status = tuple->addItem ("time", time); //event time
  status = tuple->addItem ("ngtrack", ngood_charged_track); //good charged track in event
  status = tuple->addItem ("ngntrack", ngood_neutral_track); //good neutral track in event
  status = tuple->addItem ("nptrack", npositive_track); //good positive charged track in event
  status = tuple->addItem ("nntrack", nnegative_track); //good negative charged track in event
  status = tuple->addItem ("nppions", npositive_pions); //good poitive pion tracks in event
  status = tuple->addItem ("nnpions", nnegative_pions); //good negative pion track in event
  status = tuple->addItem ("npion_pairs", npion_pairs); //number of pions paris in event

  status = tuple->addItem ("sign", sign); //number of pions paris in event
  status = tuple->addItem ("channel", channel); //decay channel of the J/psi
  status = tuple->addItem ("KK", KK); //KK decay channel of the J/psi
  status = tuple->addItem ("uu", uu); //mu-mu decay channel of the J/psi
  status = tuple->addItem ("Ku", Ku); //Kmu or muK events

	status = M.add_to_tuple(tuple);
  status = tuple->addItem ("Mee",   kM[ID_ELECTRON]);
  status = tuple->addItem ("MKK",   kM[ID_KAON]);
  status = tuple->addItem ("Muu",   kM[ID_MUON]);
  status = tuple->addItem ("Mpp",   kM[ID_PROTON]);
  status = tuple->addItem ("Mpipi", kM[ID_PION]);

	status = T.add_to_tuple(tuple);

  status = tuple->addItem ("kin_chi2", kin_chi2); 
  status = tuple->addItem ("pid_chi2", pid_chi2); 

  status = tuple->addIndexedItem ("kchi",  T.ntrack, kchi);
  status = tuple->addIndexedItem ("pchi",  T.ntrack, pchi);

  status = tuple->addIndexedItem ("probe",  T.ntrack,  prob[ID_ELECTRON]);
  status = tuple->addIndexedItem ("probmu", T.ntrack,  prob[ID_MUON]);
  status = tuple->addIndexedItem ("probpi", T.ntrack,  prob[ID_PION]);
  status = tuple->addIndexedItem ("probk",  T.ntrack,  prob[ID_KAON]);
  status = tuple->addIndexedItem ("probp",  T.ntrack,  prob[ID_PROTON]);


  return status;
}

void JpsiKK::RootEvent::init(void)
{
  T.ntrack=4;
}

StatusCode JpsiKK::RootPid::init_tuple(void)
{
  StatusCode status;
  status = tuple->addItem ("Mee", M[ID_ELECTRON]);
  status = tuple->addItem ("MKK", M[ID_KAON]);
  status = tuple->addItem ("Muu", M[ID_MUON]);
  status = tuple->addItem ("Mpp", M[ID_PROTON]);
  status = tuple->addItem ("Mpipi", M[ID_PION]);

  status = tuple->addItem ("kMee", kM[ID_ELECTRON]);
  status = tuple->addItem ("kMKK", kM[ID_KAON]);
  status = tuple->addItem ("kMuu", kM[ID_MUON]);
  status = tuple->addItem ("kMpp", kM[ID_PROTON]);
  status = tuple->addItem ("kMpipi", kM[ID_PION]);
  //info prof ParticleID package
  status = tuple->addItem ("ntrack", ntrack,0,4); //array size must be = 4
  status = tuple->addIndexedItem ("probe",  ntrack, prob[ID_ELECTRON]);
  status = tuple->addIndexedItem ("probmu",  ntrack, prob[ID_MUON]);
  status = tuple->addIndexedItem ("probpi",  ntrack, prob[ID_PION]);
  status = tuple->addIndexedItem ("probk",  ntrack, prob[ID_KAON]);
  status = tuple->addIndexedItem ("probp",  ntrack, prob[ID_PROTON]);
  //my particle id information
  status = tuple->addIndexedItem ("chi2e",  ntrack, chi2[ID_ELECTRON]);
  status = tuple->addIndexedItem ("chi2mu",  ntrack, chi2[ID_MUON]);
  status = tuple->addIndexedItem ("chi2pi",  ntrack, chi2[ID_PION]);
  status = tuple->addIndexedItem ("chi2k",  ntrack, chi2[ID_KAON]);
  status = tuple->addIndexedItem ("chi2p",  ntrack, chi2[ID_PROTON]);
  return status;
}

void JpsiKK::RootPid::init(void)
{
  ntrack=4;
}

StatusCode JpsiKK::RootMdc::init_tuple(void)
{
  StatusCode status;
	status = M.add_to_tuple(tuple);
	status = T.add_to_tuple(tuple);
  return status;
}


void JpsiKK::RootMdc::init(void)
{
  T.ntrack=4;
}


StatusCode JpsiKK::RootEmc::init_tuple(void)
{
  StatusCode status;
  status = tuple->addItem ("ntrack",       ntrack,0,ARRAY_SIZE); //good nuetral track in event
  status = tuple->addIndexedItem ("E",     ntrack, E);
  status = tuple->addIndexedItem ("theta", ntrack, theta);
  status = tuple->addIndexedItem ("phi",   ntrack, phi);
  status = tuple->addIndexedItem ("time",  ntrack, time);
  return status;
}

void JpsiKK::RootEmc::init(void)
{
  ntrack=4;
  for(int i=0;i<ntrack;i++)
  {
    E[i] = 0;
    theta[i] = -1000;;
    phi[i] = -1000;
    time[i] = -1000;
  }
}

StatusCode JpsiKK::RootDedx::init_tuple(void)
{
  StatusCode status;
  status = tuple->addItem ("ntrack", ntrack,0,4); 
  status = tuple->addIndexedItem ("chie",  ntrack, chie);
  status = tuple->addIndexedItem ("chimu",  ntrack, chimu);
  status = tuple->addIndexedItem ("chipi",  ntrack, chipi);
  status = tuple->addIndexedItem ("chik",  ntrack, chik);
  status = tuple->addIndexedItem ("chip",  ntrack, chip);
  status = tuple->addIndexedItem ("probPH",  ntrack, probPH);
  status = tuple->addIndexedItem ("normPH",  ntrack, normPH);
  //status = tuple->addIndexedItem ("probe",  ntrack, probe);
  //status = tuple->addIndexedItem ("probmu",  ntrack, probmu);
  //status = tuple->addIndexedItem ("probpi",  ntrack, probpi);
  //status = tuple->addIndexedItem ("probk",  ntrack, probk);
  //status = tuple->addIndexedItem ("probp",  ntrack, probp);
  return status;
}

void JpsiKK::RootDedx::init(void)
{
  ntrack=4;
  for(int i=0;i<ntrack;i++)
  {
    chie[i]=-1000;
    chimu[i]=-1000;
    chipi[i]=-1000;
    chik[i]=-1000;
    chip[i]=-1000;
    probPH[i]=0;
    normPH[i]=0;
  }
}

StatusCode JpsiKK::RootTof::init_tuple(void)
{
  StatusCode status;
  status = tuple->addItem ("ntrack", ntrack,0,4); 
  status = tuple->addIndexedItem ("ID",  ntrack, tofID);
  status = tuple->addIndexedItem ("t",  ntrack, t);
  status = tuple->addIndexedItem ("dt",  ntrack, dt);
  status = tuple->addIndexedItem ("t0",  ntrack, t0);
  status = tuple->addIndexedItem ("chie",  ntrack, chie);
  status = tuple->addIndexedItem ("chimu",  ntrack, chimu);
  status = tuple->addIndexedItem ("chipi",  ntrack, chipi);
  status = tuple->addIndexedItem ("chik",  ntrack, chik);
  status = tuple->addIndexedItem ("chip",  ntrack, chip);
  status = tuple->addIndexedItem ("beta",  ntrack, beta);
  status = tuple->addIndexedItem ("te",  ntrack, te);
  status = tuple->addIndexedItem ("tmu",  ntrack, tmu);
  status = tuple->addIndexedItem ("tpi",  ntrack, tpi);
  status = tuple->addIndexedItem ("tk",  ntrack, tk);
  status = tuple->addIndexedItem ("tp",  ntrack, tp);
  return status;
}

void JpsiKK::RootTof::init(void)
{
  ntrack=4;
  for(int i=0;i<ntrack;i++)
  {
    tofID[i]=-1000;
    t[i]=-1000;
    dt[i]=-1000;
    t0[i]=-1000;
    chie[i]=-1000;
    chimu[i]=-1000;
    chipi[i]=-1000;
    chik[i]=-1000;
    chip[i]=-1000;
    beta[i]=-1000;
    te[i]=-1000;
    tmu[i]=-1000;
    tpi[i]=-1000;
    tk[i]=-1000;
    tp[i]=-1000;
  }
}

StatusCode JpsiKK::RootMC::init_tuple(void)
{
  StatusCode status;
  status = tuple->addItem ("psip_decay", psip_decay); //flag for psip decay
  status = tuple->addItem ("jpsi_decay", jpsi_decay); //flag for jpsi decay 
  status = tuple->addItem ("KK", KK);               //KK event
  status = tuple->addItem ("uu", uu);               //mu mu event
  status = tuple->addItem ("oo", oo);               //other event
  status = tuple->addItem ("ntrack", ntrack,0,4); 
  status = tuple->addIndexedItem ("id",    ntrack, pid);
  status = tuple->addIndexedItem ("q",     ntrack, q);
  status = tuple->addIndexedItem ("E",     ntrack, E);
  status = tuple->addIndexedItem ("p",     ntrack, p);
  status = tuple->addIndexedItem ("px",    ntrack, px);
  status = tuple->addIndexedItem ("py",    ntrack, py);
  status = tuple->addIndexedItem ("pz",    ntrack, pz);
  status = tuple->addIndexedItem ("pt",    ntrack, pt);
  status = tuple->addIndexedItem ("theta", ntrack, theta);
  status = tuple->addIndexedItem ("phi",   ntrack, phi);
  return status;
}

void JpsiKK::RootMC::init(void)
{
  ntrack=4;
}

StatusCode JpsiKK::RootMCTopo::init_tuple(void)
{
  StatusCode status;
  status = tuple->addItem("indexmc", m_idxmc, 0, 100);
  status = tuple->addIndexedItem("pdgid", m_idxmc, m_pdgid);
  status = tuple->addIndexedItem("motheridx", m_idxmc, m_motheridx);
  status = tuple->addIndexedItem("idx", m_idxmc, m_idx);
  status = tuple->addItem("hash", m_hash);
  return status;
}

void JpsiKK::RootMCTopo::init(void)
{
  m_idxmc=0;
}



StatusCode JpsiKK::execute()
{
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "executing" << endreq;

  SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
  fEvent.run=eventHeader->runNumber();
  fEvent.event=eventHeader->eventNumber();
  fEvent.time=eventHeader->time();
  bool isprint=false;
  if(event_proceed<10) isprint=true;
  if(10 <= event_proceed && event_proceed < 100 && event_proceed % 10 ==0) isprint=true;
  if(100 <= event_proceed && event_proceed < 1000 && event_proceed % 100 ==0) isprint = true;
  if(1000 <= event_proceed && event_proceed < 10000 && event_proceed % 1000 ==0) isprint = true;
  if(10000 <= event_proceed && event_proceed % 10000 ==0) isprint = true;
  if(isprint)
  {
    std::cout << "proceed event: " << setw(15) << event_proceed<<", ";
    std::cout << "event write: "   << setw(15) << event_write<< ",  ";
    std::cout << "kaons: "   << setw(15) << event_with_kaons << ",  ";
    std::cout << "muons "   << setw(15) << event_with_muons << ",  ";
    std::cout << "good knm fit "   << setw(15) << good_kinematic_fit;
    std::cout << std::endl;
  }
  event_proceed++;

  //  Get information about reconstructed events
  SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
  SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol);
  SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(),  EventModel::MC::McParticleCol);


  //fill initial value of the selected event
  fEvent.init();

  std::list<EvtRecTrackIterator> good_charged_tracks=createGoodChargedTrackList(cfg, evtRecEvent, evtRecTrkCol);
  std::list<EvtRecTrackIterator> good_neutral_tracks=createGoodNeutralTrackList(cfg, evtRecEvent, evtRecTrkCol);

  //SELECTION CODE
  if( cfg.MAX_NEUTRAL_TRACKS < good_neutral_tracks.size()) return StatusCode::SUCCESS;
  //SELECTION CODE
  if( good_charged_tracks.size() < cfg.MIN_CHARGED_TRACKS || cfg.MAX_CHARGED_TRACKS < good_charged_tracks.size()) return StatusCode::SUCCESS;

	cout <<  "Number of good charged tracks: " << good_charged_tracks.size() << endl;
  
  TrackList_t charged_tracks; //selected tracks for specific cut
  TrackList_t positive_charged_tracks; //selected tracks for specific cut
  TrackList_t negative_charged_tracks; //selected tracks for specific cut
  TrackList_t positive_pion_tracks; //selected pion tracks for specific cut
  TrackList_t negative_pion_tracks; //selected pion tracks for specific cut
  TrackList_t other_positive_tracks; //other positive tracks for specific cut
  TrackList_t other_negative_tracks; //other positive tracks for specific cut
  for(TrackList_t::iterator track=good_charged_tracks.begin(); track!=good_charged_tracks.end(); track++)
  {
    EvtRecTrackIterator & itTrk = *track;
    if(!(*itTrk)->isMdcTrackValid()) continue; 
    RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
    double c = fabs(cos(mdcTrk->theta()));
    double p = mdcTrk->p();
    double q = mdcTrk->charge();
    bool barrel = c < cfg.EMC_BARREL_MAX_COS_THETA;
    if(barrel) 
    {
      if(q>0) 
      {
        positive_charged_tracks.push_back(itTrk);
        if(cfg.MIN_PION_MOMENTUM < p && p<cfg.MAX_PION_MOMENTUM) 
        {
          positive_pion_tracks.push_back(itTrk);
        }
        if(p>std::min(cfg.MIN_KAON_MOMENTUM, cfg.MIN_MUON_MOMENTUM))
        {
          if((*itTrk)->isEmcShowerValid())
          {
            other_positive_tracks.push_back(itTrk);
          }
        }
      }
      if(q<0) 
      {
        negative_charged_tracks.push_back(itTrk);
        if(cfg.MIN_PION_MOMENTUM < p &&  p<cfg.MAX_PION_MOMENTUM) 
        {
          negative_pion_tracks.push_back(itTrk);
        }
        if(p>std::min(cfg.MIN_KAON_MOMENTUM, cfg.MIN_MUON_MOMENTUM))
        {
          if((*itTrk)->isEmcShowerValid())
          {
            other_negative_tracks.push_back(itTrk);
          }
        }
      }
      charged_tracks.push_back(itTrk);
    }
  }

  //create pion pairs
  TrackPairList_t pion_pairs;
  for(TrackList_t::iterator i=negative_pion_tracks.begin(); i!=negative_pion_tracks.end(); ++i)
    for(TrackList_t::iterator j=positive_pion_tracks.begin(); j!=positive_pion_tracks.end(); ++j)
    {
      TrackPair_t pair(*i,*j);
      double M_recoil = get_recoil__mass(pair, PION_MASS,  cfg.CENTER_MASS_ENERGY);
      if(cfg.MIN_RECOIL_MASS < M_recoil && M_recoil < cfg.MAX_RECOIL_MASS) 
      {
        pion_pairs.push_back(pair);
      }
    }
	cout <<  "Number of pion pairs: " << pion_pairs.size() << endl;
  //SELECTION CODE we must have at list one pion pair
  if(pion_pairs.empty()) return StatusCode::SUCCESS; //we must have at list one pion pair
	//if(pion_pairs.size()!=1) return StatusCode::SUCCESS; //exacly one pion pair
  TrackPair_t pion_pair = pion_pairs.front();
	if(pion_pairs.size()>1)
	{
		//find the best pion pair using closest value to JPSI_MASS
		for(TrackPairList_t::iterator p=pion_pairs.begin();p!=pion_pairs.end();p++)
		{
			if(fabs(get_recoil__mass(*p,PION_MASS, cfg.CENTER_MASS_ENERGY) - JPSI_MASS) <  fabs(get_recoil__mass(pion_pair,PION_MASS, cfg.CENTER_MASS_ENERGY) - JPSI_MASS)) pion_pair = *p;
		}
	}
	//now we have one pion pair candidate

  //SELECTION CODE
  //keep only specific signature
  //if(positive_charged_tracks.size()!=2 || negative_charged_tracks.size()!=2) return StatusCode::SUCCESS;

  log << MSG::ERROR << "good charged tracks: " << charged_tracks.size() << " (" << negative_charged_tracks.size() << ", " << positive_charged_tracks.size() << endmsg;
  log << MSG::ERROR << "pions: " << negative_pion_tracks.size()  << ", " << positive_pion_tracks.size() << endmsg;
  log << MSG::ERROR << "other: " << other_negative_tracks.size()  << ", " << other_positive_tracks.size() << endmsg;
  log << MSG::ERROR << "pion pairs: " << pion_pairs.size() << endmsg;


  //if no other particles
	if(other_negative_tracks.empty() && other_positive_tracks.empty()) return StatusCode::SUCCESS;



	//SelectionHelper_t KFP(CENTER_MASS_ENERGY);
	SelectionHelper_t negative_sh(cfg.CENTER_MASS_ENERGY, evtRecTrkCol->end());
	SelectionHelper_t positive_sh(cfg.CENTER_MASS_ENERGY, evtRecTrkCol->end());

	if(!other_negative_tracks.empty()) 
	{
		kinfit(pion_pairs,  other_negative_tracks,  negative_sh);
		negative_sh.totalPass(cfg);
	}

	if(!other_positive_tracks.empty()) 
	{
		kinfit(pion_pairs,  other_positive_tracks,  positive_sh);
		positive_sh.totalPass(cfg);
	}

	TrackVector_t Tracks;
	std::vector<HepLorentzVector> Pkf;
	SelectionHelper_t * sh;

	fEvent.sign = (int(negative_sh()) << 1 ) + int(positive_sh());
	fEvent.KK = 0;
	fEvent.uu = 0;
	fEvent.Ku = 0;
	switch(fEvent.sign)
	{
		case OTHER_NO_TRACK:
			return StatusCode::SUCCESS;
			break;

		case OTHER_NEGATIVE_TRACK: //one negative track
			fEvent.channel = negative_sh.channel;
			Tracks = negative_sh.tracks;
			Pkf = negative_sh.P;
			//add missing positive tracks
			Tracks.push_back(negative_sh.end);
			sh = & negative_sh;
			break;

		case OTHER_POSITIVE_TRACK: //one positive track
			fEvent.channel = positive_sh.channel;
			Tracks = positive_sh.tracks;
			Pkf = positive_sh.P;
			//add missing negative tracks
			Tracks.push_back(positive_sh.end);
			//negative tracks go first
			std::swap(Tracks[2], Tracks[3]);
			std::swap(Pkf[2],  Pkf[3]);
			sh = & positive_sh;
			break;

		case OTHER_TWO_TRACKS:
			if(negative_sh.channel == ID_KAON && positive_sh.channel == ID_KAON)
			{
				fEvent.channel = CHAN_KAONS;
			}

			if(negative_sh.channel == ID_MUON && positive_sh.channel == ID_MUON)
			{
				fEvent.channel = CHAN_MUONS;
			}

			if(negative_sh.channel == ID_KAON && positive_sh.channel == ID_MUON)
			{
				fEvent.channel = CHAN_KAON_MUON;
			}

			if(negative_sh.channel == 1 && positive_sh.channel == 0)
			{
				fEvent.channel = CHAN_MUON_KAON;
			}
			Pkf = negative_sh.P;
			Tracks = negative_sh.tracks;
			Tracks.push_back(positive_sh.Tracks[2]);
			sh = & negative_sh;
			break;
		default:
			return StatusCode::SUCCESS;
			break;
	}

	switch(fEvent.channel)
	{
		case CHAN_KAONS:
			fEvent.KK = 1;
			break;
		case CHAN_MUONS:
			fEvent.uu = 1;
			break;
		case CHAN_KAON_MUON:
		case CHAN_MUON_KAON:
			fEvent.Ku = 1;
			break;
		default:
			break;
	}

  //now fill the tuples

  //some statistics information
	std::cerr << "DEBUG: BEFORE header fillin" << std::endl;
  fEvent.ngood_charged_track = good_charged_tracks.size();
  fEvent.ngood_neutral_track = good_neutral_tracks.size();
  fEvent.npositive_track = positive_charged_tracks.size();
  fEvent.nnegative_track = negative_charged_tracks.size();
  fEvent.npositive_pions = positive_pion_tracks.size();
  fEvent.nnegative_pions = negative_pion_tracks.size();
  fEvent.npion_pairs = pion_pairs.size();
  // fill the decay channel of the J/psi 0 - kaons, 1 --muons
  //fEvent.channel = channel; 
  fEvent.kin_chi2 = sh.kin_chi2; //kinematic_chi2;
  fEvent.pid_chi2 = sh.mypid_chi2[fEvent.channel]; //pchi2[channel];
	std::cerr << "DEBUG: BEFORE Pkf inv masses" << std::endl;
  fEvent.Minv = (Pkf[2]+Pkf[3]).m();
  fEvent.M012 = (Pkf[0]+Pkf[1]+Pkf[2]).m();
  fEvent.M013 = (Pkf[0]+Pkf[1]+Pkf[3]).m();
  fEvent.M03 =  (Pkf[0]+Pkf[3]).m();
  fEvent.M12 =  (Pkf[1]+Pkf[2]).m();
  fEvent.M01 =  (Pkf[0]+Pkf[1]).m();
  HepLorentzVector P_psip(cfg.CENTER_MASS_ENERGY*sin(0.011),0,0,cfg.CENTER_MASS_ENERGY); //initial vector of psip
	std::cerr << "DEBUG: BEFORE Mrecoil" << std::endl;
  fEvent.Mrecoil = (P_psip - Pkf[0] - Pkf[1]).m();

  fEvent.ntrack = 4;
	std::cerr << "DEBUG: BEFORE fEvent filling" << std::endl;
  for ( int i=0;i<4;i++)
  {
    fEvent.T.q[i]  = i%2 == 0 ? -1 : +1;
    fEvent.T.E[i]  = Pkf[i].e();
    fEvent.T.px[i] = Pkf[i].px();
    fEvent.T.py[i] = Pkf[i].py();
    fEvent.T.pz[i] = Pkf[i].pz();
    fEvent.T.p[i]  = sqrt(sq(Pkf[i].px())+sq(Pkf[i].py())+sq(Pkf[i].pz()));
    fEvent.T.pt[i] = sqrt(sq(Pkf[i].px())+sq(Pkf[i].py()));
    fEvent.T.theta[i]= Pkf[i].theta();
    fEvent.T.phi[i] = Pkf[i].phi();
    fEvent.T.x[i]=0;
    fEvent.T.y[i]=0;
    fEvent.T.z[i]=0;
    fEvent.T.r[i]=0;
    fEvent.T.vxy[i]=0;
    fEvent.T.vz[i]=0;
    fEvent.T.vphi[i]=0;
  }
	std::cerr << "DEBUG: BEFORE particle id" << std::endl;

  ParticleID * PID = ParticleID::instance();
  PID->init();
  PID->setMethod(PID->methodProbability());
  PID->setChiMinCut(4);
  PID->usePidSys(PID->useDedx() || PID->useTof());
  PID->identify(PID->all()); 

  fEvent.ntrack=4;
  //fPid.ntrack=4;
  fMdc.ntrack=4;
  fDedx.ntrack=4;
  fEmc.ntrack=4;
  fTof.ntrack=4;
  fMdc.Mrecoil = get_recoil__mass(pion_pair, PION_MASS,  cfg.CENTER_MASS_ENERGY);
  //fMdc.Minv    = sqrt(get_invariant_mass2(result_pair,XMASS[channel]));
  //EvtRecTrackIterator itTrk[4] = {pion_pair.first, pion_pair.second, result_pair.first, result_pair.second};
	std::cerr << "DEBUG: BEFORE loop" << std::endl;
  for(int i=0;i<4;i++)
  {
		std::cerr << "DEBUG: track " << i << std::endl;
		if(Tracks[i]==evtRecTrkCol->end()) continue;
		std::cerr << "DEBUG: Before emcTrk:" << std::endl;
    if(i>1)
    {
      RecEmcShower *emcTrk = (*Tracks[i])->emcShower();
      fEmc.E[i] = emcTrk->energy();
      fEmc.theta[i] = emcTrk->theta();
      fEmc.phi[i] = emcTrk->phi();
      fEmc.time[i] = emcTrk->time();
      fMdc.E[i] = fEmc.E[i];
    }
		std::cerr << "DEBUG: Before mdcTrk:" << std::endl;
    RecMdcTrack  *mdcTrk = (*Tracks[i])->mdcTrack();
    fMdc.T.trackId[i] = mdcTrk->trackId();
    fMdc.T.q[i] = mdcTrk->charge(); 
    fMdc.T.p[i] = mdcTrk->p();
    fMdc.T.px[i]= mdcTrk->px();
    fMdc.T.py[i]= mdcTrk->py();
    fMdc.T.pz[i]= mdcTrk->pz();
    fMdc.T.theta[i]= mdcTrk->theta();
    fMdc.T.phi[i] = mdcTrk->phi();
    fMdc.T.x[i]  = mdcTrk->x();
    fMdc.T.y[i]  = mdcTrk->y();
    fMdc.T.z[i]  = mdcTrk->z();
    fMdc.T.x[i]  = mdcTrk->x();
    fMdc.T.y[i]  = mdcTrk->y();
    fMdc.T.z[i]  = mdcTrk->z();
    double rvxy,rvz,rvphi;
		std::cerr << "DEBUG: calculate vertex:" << std::endl;
    calculate_vertex(mdcTrk,rvxy,rvz,rvphi); 
    fMdc.T.vxy[i] = rvxy;
    fMdc.T.vz[i]  = rvz; 
    fMdc.T.vphi[i] = rvphi; 


    //dedx information
		std::cerr << "DEBUG: DedX:" << std::endl;
    if((*Tracks[i])->isMdcDedxValid())
    {
      RecMdcDedx* dedxTrk = (*Tracks[i])->mdcDedx();
      fDedx.chie[i] = dedxTrk->chiE();
      fDedx.chimu[i] = dedxTrk->chiMu();
      fDedx.chipi[i] = dedxTrk->chiPi();
      fDedx.chik[i] = dedxTrk->chiK();
      fDedx.chip[i] = dedxTrk->chiP();
      //fDedx.ghit[i] = dedxTrk->numGoodHits();
      //fDedx.thit[i] = dedxTrk->numTotalHits();
      fDedx.probPH[i] = dedxTrk->probPH();
      fDedx.normPH[i] = dedxTrk->normPH();
      //fDedx.e[i] = dedxTrk->getDedxExpect(0);
      //fDedx.mu[i] = dedxTrk->getDedxExpect(1);
      //fDedx.pi[i] = dedxTrk->getDedxExpect(2);
      //fDedx.K[i] = dedxTrk->getDedxExpect(3);
      //fDedx.p[i] = dedxTrk->getDedxExpect(4);
      //fDedx.pid[i]=dedxTrk->particleId();
    }
		std::cerr << "DEBUG: Before TOF:" << std::endl;
    if((*Tracks[i])->isTofTrackValid())
    {
      SmartRefVector<RecTofTrack> tofTrkCol = (*Tracks[i])->tofTrack();
      SmartRefVector<RecTofTrack>::iterator tofTrk = tofTrkCol.begin();
      TofHitStatus *hitst = new TofHitStatus;
      std::vector<int> tofecount;
      int goodtofetrk=0;
      for(tofTrk = tofTrkCol.begin(); tofTrk!=tofTrkCol.end(); tofTrk++,goodtofetrk++)
      {
        unsigned int st = (*tofTrk)->status();
        hitst->setStatus(st);
        //if(  (hitst->is_barrel()) ) continue;
        if( !(hitst->is_counter()) ) continue;
        tofecount.push_back(goodtofetrk);
      }
      delete hitst;
      if(!tofecount.empty())
      {
        tofTrk = tofTrkCol.begin()+tofecount[0];
        fTof.tofID[i] = (*tofTrk)->tofID();
        fTof.t0[i] = (*tofTrk)->t0();
        fTof.t[i] = (*tofTrk)->tof();
        fTof.dt[i] = (*tofTrk)->errtof();
        fTof.beta[i] = (*tofTrk)->beta();
        fTof.te[i] = (*tofTrk)->texpElectron();
        fTof.tmu[i]= (*tofTrk)->texpMuon();
        fTof.tpi[i]= (*tofTrk)->texpPion();
        fTof.tk[i] = (*tofTrk)->texpKaon();
        fTof.tp[i] = (*tofTrk)->texpProton();
        if(fTof.dt[i]>0)
        {
          fTof.chie[i]  = (fTof.t[i] - fTof.te[i])  /  fTof.dt[i];
          fTof.chimu[i] = (fTof.t[i] - fTof.tmu[i]) /  fTof.dt[i];
          fTof.chipi[i] = (fTof.t[i] - fTof.tpi[i]) /  fTof.dt[i];
          fTof.chik[i]  = (fTof.t[i] - fTof.tk[i])  /  fTof.dt[i];
          fTof.chip[i]  = (fTof.t[i] - fTof.tp[i])  /  fTof.dt[i];
        }
      }
    }

    //fill particle id

		std::cerr << "DEBUG: Before PID:" << std::endl;
    PID->setRecTrack((*Tracks[i]));
    PID->calculate();
    if(PID->IsPidInfoValid())
    {
      //fPid.prob[ID_ELECTRON][i] = PID->probElectron();
      //fPid.prob[ID_MUON][i]     = PID->probMuon();
      //fPid.prob[ID_PION][i]     = PID->probPion();
      //fPid.prob[ID_KAON][i]     = PID->probKaon();
      //fPid.prob[ID_PROTON][i]   = PID->probProton();
    }
		std::cerr << "DEBUG: Before chi2:" << std::endl;
    vector<double> chi2 = get_chi2(Tracks[i]);
    for(int pid=0;pid<5;pid++)
    {
      //fPid.chi2[pid][i]   = chi2[pid];
    }

  }
	std::cerr << "DEBUG: Before fPID:" << std::endl;
  for(int i=0;i<5;i++)
  {
    //fPid.M[i]    = sqrt(get_invariant_mass2(result_pair,XMASS[i]));
    HepLorentzVector p1(Pkf[2].vect(), XMASS[i]);
    HepLorentzVector p2(Pkf[3].vect(), XMASS[i]);
    //fPid.kM[i] = (p1+p2).m();
  };


  if(fEvent.run<0)
  {
    //check the MC information
    if(!mcParticleCol)
    {
      log << MSG::ERROR << "Could not retrieve McParticelCol" << endreq;
      return StatusCode::FAILURE;
    }
    //Fill MC TOPO INFORMATION
    bool psipDecay = false;
    int rootIndex = -1;
    Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
    int m_numParticle = 0;
    for (; iter_mc != mcParticleCol->end(); iter_mc++)
    {
      if ((*iter_mc)->primaryParticle()) continue;
      if (!(*iter_mc)->decayFromGenerator()) continue;
      if ((*iter_mc)->particleProperty()==100443)
      {
        psipDecay = true;
        rootIndex = (*iter_mc)->trackIndex();
      }
      if (!psipDecay) continue;
      int pdgid = (*iter_mc)->particleProperty();
      int mcidx = ((*iter_mc)->mother()).trackIndex() - rootIndex;
      fMCTopo.m_pdgid[m_numParticle] = pdgid;
      fMCTopo.m_motheridx[m_numParticle] = mcidx;
      fMCTopo.m_idx[m_numParticle] = (*iter_mc)->trackIndex()-rootIndex;
      fMCTopo.m_hash=0; //no hash calculation now
      m_numParticle += 1;
    }
    fMCTopo.m_idxmc = m_numParticle;

    //Fill my mc truth information
    fMC.ntrack=4;
    HepLorentzVector MCPpion[2];
    HepLorentzVector MCPkaon_or_muon[2];
    bool pi_minus(false);
    bool pi_plus(false);
    bool K_minus(false);
    bool K_plus(false);
    bool mu_minus(false);
    bool mu_plus(false);
    int mytrack=0;
    fMC.psip_decay = 0;
    fMC.jpsi_decay = 0;
    fMC.KK = 0;
    fMC.uu = 0;
    fMC.oo = 0;
    psipDecay=false;
    rootIndex=-1;
    for (iter_mc = mcParticleCol->begin(); iter_mc != mcParticleCol->end(); iter_mc++)
    {
      if ((*iter_mc)->primaryParticle()) continue;
      if (!(*iter_mc)->decayFromGenerator()) continue;
      //if ( ((*iter_mc)->mother()).trackIndex()<3 ) continue;
      if ((*iter_mc)->particleProperty()==100443)
      {
        psipDecay = true;
        fMC.psip_decay=1;
        rootIndex = (*iter_mc)->trackIndex();
      }
      if (!psipDecay) continue;
      if ((*iter_mc)->particleProperty()==443)
      {
        fMC.jpsi_decay=1;
      }
      if (fMC.jpsi_decay!=1) continue;
      if((*iter_mc)->particleProperty() == +211) 
      {
        pi_plus=true;
        MCPpion[1] = (*iter_mc)->initialFourMomentum();
        fMC.pid[1]=211;
        mytrack++;
      }
      if((*iter_mc)->particleProperty() == -211) 
      {
        MCPpion[0] = (*iter_mc)->initialFourMomentum();
        pi_minus=true;
        fMC.pid[0]=-211;
        mytrack++;
      }
      if( ! pi_plus && !pi_minus) continue; //keep only psip to Jpsi pi pi decay
      switch((*iter_mc)->particleProperty())
      {
        case -13:
          MCPkaon_or_muon[0] = (*iter_mc)->initialFourMomentum();
          mu_minus = true;
          mytrack++;
          fMC.pid[2]=-13;
          break;
        case +13:
          MCPkaon_or_muon[1] = (*iter_mc)->initialFourMomentum();
          mu_plus = true;
          mytrack++;
          fMC.pid[3]=13;
          break;
        case -321:
          MCPkaon_or_muon[0] = (*iter_mc)->initialFourMomentum();
          K_minus=true;
          mytrack++;
          fMC.pid[2]=-321;
          break;
        case +321:
          MCPkaon_or_muon[1] = (*iter_mc)->initialFourMomentum();
          mytrack++;
          K_plus = true;
          fMC.pid[3]=321;
          break;
      };
    }
    if(K_plus && K_minus) 
    {
      fMC.KK=1;
      fMC.oo=0;
    }
    if(mu_plus && mu_minus) 
    {
      fMC.uu=1;
      fMC.oo=0;
    }
    if(fMC.KK==1 && fMC.uu==1) fMC.oo=1;
    if(mytrack!=4) fMC.oo=1;
    if(fMC.KK==1 || fMC.uu==1)
    {
      vector<HepLorentzVector> P(4);
      P[0]=MCPpion[0];
      P[1]=MCPpion[1];
      P[2]=MCPkaon_or_muon[0];
      P[3]=MCPkaon_or_muon[1];
      for(int i=0;i<4;i++)
      {
        fMC.q[i] = 0; 
        fMC.E[i] = P[i].e();
        fMC.p[i] = P[i].rho();
        fMC.px[i]= P[i].px();
        fMC.py[i]= P[i].py();
        fMC.pz[i]= P[i].pz();
        fMC.pt[i]= P[i].perp();
        fMC.theta[i]= P[i].theta();
        fMC.phi[i] = P[i].phi();
      }
    }
    else
    {
      Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
      for (; iter_mc != mcParticleCol->end(); iter_mc++)
      {
        if ((*iter_mc)->primaryParticle()) continue;
        if (!(*iter_mc)->decayFromGenerator()) continue;
        HepLorentzVector P = (*iter_mc)->initialFourMomentum();
        int pid = (*iter_mc)->particleProperty();
        Hep3Vector p_mc = P.vect();
        for(int i=0;i<4;i++)
        {
          Hep3Vector p_rec = Pkf[i].vect();
          Hep3Vector dp = p_rec - p_mc;
          if(dp.mag()/std::min(p_rec.mag(), p_mc.mag()) < 0.05)
          {
            fMC.pid[i] = pid;
            //fMC.q[i] = 0; 
            fMC.E[i] = P.e();
            fMC.p[i] = P.rho();
            fMC.px[i]= P.px();
            fMC.py[i]= P.py();
            fMC.pz[i]= P.pz();
            fMC.pt[i]= P.perp();
            fMC.theta[i]= P.theta();
            fMC.phi[i] = P.phi();
          }
        }
      }
    }
  }


  fNeutral.ntrack=std::min(good_neutral_tracks.size(), size_t(RootEmc::ARRAY_SIZE));
  int idx=0;
  for(list<EvtRecTrackIterator>::iterator track=good_neutral_tracks.begin(); track!=good_neutral_tracks.end(); track++)
  {
    EvtRecTrackIterator  itTrk = *track;
    RecEmcShower *emcTrk = (*itTrk)->emcShower();
    fNeutral.E[idx]  =  emcTrk->energy();
    fNeutral.theta[idx] =  emcTrk->theta();
    fNeutral.phi[idx] =  emcTrk->phi();
    fNeutral.time[idx] = emcTrk->time();
    idx++;
  }

  if(fEvent.run<0) 
  {
    fMC.tuple->write();
    fMCTopo.tuple->write();
  }
  fEvent.tuple->write();
  //fPid.tuple->write();
  fMdc.tuple->write();
  fEmc.tuple->write();
  fDedx.tuple->write();
  fTof.tuple->write();
  fNeutral.tuple->write();
  event_write++;
  return StatusCode::SUCCESS;
}

StatusCode JpsiKK::finalize()
{
  std::cout << "Event proceed: " << event_proceed << std::endl;
  std::cout << "Event selected: " << event_write << std::endl;
  std::cout << "Event with kaons: " << event_with_kaons << std::endl;
  std::cout << "Event with muons: " << event_with_muons << std::endl;
  std::cout << "Event with pions: " << event_with_pions << std::endl;
  std::cout << "Event with electron: " << event_with_electrons << std::endl;
  std::cout << "Event with proton: " << event_with_protons << std::endl;
  std::cout << "Good kinematic fits: " << good_kinematic_fit << std::endl;
  return StatusCode::SUCCESS;
}
