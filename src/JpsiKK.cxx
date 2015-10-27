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




#include <algorithm>
#include <limits>



#include "Defs.h"
#include "PhysConst.h"
#include "Utils.h"
#include "MyPid.h"
#include "SelectionHelper.h"

#include "KinematicFit.h"

enum
{
	OTHER_NEGATIVE_TRACK=-1, 
	OTHER_POSITIVE_TRACK=+1, 
	OTHER_TWO_TRACKS=0, 
	OTHER_NO_TRACK=99
};

enum
{
	CHAN_KAONS = ID_KAON, 
	CHAN_MUONS = ID_MUON, 
	CHAN_KAON_MUON  = 10, 
	CHAN_MUON_KAON  = 11
};

double XMASS[5] = {KAON_MASS, MUON_MASS, ELECTRON_MASS, PION_MASS, PROTON_MASS};

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


StatusCode JpsiKK::initialize(void)
{
  MsgStream log(msgSvc(), name());

  log << MSG::INFO << "in initialize()" << endmsg;
  event_proceed=0;
  event_write = 0;
  event_with_kaons=0;
  event_with_muons=0;
  event_with_kaons_and_muons=0;
  event_with_pions=0;
  event_with_electrons=0;
  event_with_protons=0;
  good_kinematic_fit=0;
	if(cfg.CENTER_MASS_ENERGY == 0) cfg.CENTER_MASS_ENERGY = PSIP_MASS;

  StatusCode status;
  status = init_tuple(this, fEvent,      "FILE1/event","Signal events pi+pi- K+K-, or pi+pi- mu+mu-",log);
//  status = init_tuple(this, fPid,        "FILE1/pid","particle id",log);
  status = init_tuple(this, fMdc,        "FILE1/mdc","Mdc info for signal",log);
  status = init_tuple(this, fDedx,       "FILE1/dedx","Dedx info for signal",log);
  status = init_tuple(this, fEmc,        "FILE1/emc","Emc info for signal",log);
  status = init_tuple(this, fTof,        "FILE1/tof","Tof info for signal",log);
  status = init_tuple(this, fNeutral,    "FILE1/neutral","Good neutral tracks",log);
  status = init_tuple(this, fMC,         "FILE1/mc","Monte Carlo truth information",log);
  status = init_tuple(this, fMCTopo,     "FILE1/mctopo","Monte Carlo truth information topology",log);

  return status;
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
    std::cout << "Kmu "   << setw(15) << event_with_kaons_and_muons << ",  ";
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

  //log << MSG::ERROR << "good charged tracks: " << charged_tracks.size() << " (" << negative_charged_tracks.size() << ", " << positive_charged_tracks.size() << endmsg;
  //log << MSG::ERROR << "pions: " << negative_pion_tracks.size()  << ", " << positive_pion_tracks.size() << endmsg;
  //log << MSG::ERROR << "other: " << other_negative_tracks.size()  << ", " << other_positive_tracks.size() << endmsg;
  //log << MSG::ERROR << "pion pairs: " << pion_pairs.size() << endmsg;


  //if no other particles
	if(other_negative_tracks.empty() && other_positive_tracks.empty()) return StatusCode::SUCCESS;


	SelectionHelper_t negative_sh(cfg);
	SelectionHelper_t positive_sh(cfg);


	if(!other_negative_tracks.empty()) 
	{
		//kinfit(pion_pair,  other_negative_tracks,  negative_sh);
		negative_sh.kinfit(pion_pair,  other_negative_tracks.front());
		negative_sh.totalPass();
	}

	if(!other_positive_tracks.empty()) 
	{
		//kinfit(pion_pair,  other_positive_tracks,  positive_sh);
		positive_sh.kinfit(pion_pair,  other_positive_tracks.front());
		positive_sh.totalPass();
	}

	TrackVector_t Tracks;
	std::vector<HepLorentzVector> Pkf;
	SelectionHelper_t * sh;

	fEvent.sign = int(positive_sh.pass) - int(negative_sh.pass);
	fEvent.KK = 0;
	fEvent.uu = 0;
	fEvent.Ku = 0;
	switch(fEvent.sign)
	{
		case OTHER_NEGATIVE_TRACK: //one negative track
			fEvent.channel = negative_sh.channel;
			Tracks = negative_sh.tracks;
			Pkf = negative_sh.KF[negative_sh.channel].P;
			//add missing positive tracks
			Tracks.push_back(evtRecTrkCol->end());
			sh = & negative_sh;
			break;

		case OTHER_POSITIVE_TRACK: //one positive track
			fEvent.channel = positive_sh.channel;
			Tracks = positive_sh.tracks;
			Pkf = positive_sh.KF[positive_sh.channel].P;
			//add missing negative tracks
			Tracks.push_back(evtRecTrkCol->end());
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
			Pkf = negative_sh.KF[negative_sh.channel].P;
			Tracks = negative_sh.tracks;
			Tracks.push_back(positive_sh.tracks[2]);
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
			event_with_kaons++;
			break;
		case CHAN_MUONS:
			fEvent.uu = 1;
			event_with_muons++;
			break;
		case CHAN_KAON_MUON:
		case CHAN_MUON_KAON:
			fEvent.Ku = 1;
			event_with_kaons_and_muons++;
			break;
		default:
			return StatusCode::SUCCESS;
			break;
	}



  //now fill the tuples

  //some statistics information
  fEvent.ngood_charged_track = good_charged_tracks.size();
  fEvent.ngood_neutral_track = good_neutral_tracks.size();
  fEvent.npositive_track = positive_charged_tracks.size();
  fEvent.nnegative_track = negative_charged_tracks.size();
  fEvent.npositive_pions = positive_pion_tracks.size();
  fEvent.nnegative_pions = negative_pion_tracks.size();
  fEvent.npion_pairs = pion_pairs.size();
  // fill the decay channel of the J/psi 0 - kaons, 1 --muons
  //fEvent.channel = channel; 
  fEvent.kin_chi2 = sh -> KF[sh->channel].chi2; //kinematic_chi2;

	fEvent.npid = 5;
	for(int pid=0;pid<5;pid++)
	{
		fEvent.kchi[pid] = sh->KF[pid].chi2;
		fEvent.pchi[pid] = sh->pid_chi2[pid];
		fEvent.kM23[pid] = (sh->KF[pid].P[2] + sh->KF[pid].P[3]).m();
	}

	switch(fEvent.channel)
	{
		case CHAN_KAONS:
		case CHAN_MUONS:
			fEvent.pid_chi2 = sh -> pid_chi2[fEvent.channel]; //pchi2[channel];
			break;
		case CHAN_KAON_MUON:
		case CHAN_MUON_KAON:
			fEvent.pid_chi2 = 0.5*(positive_sh.pid_chi2[positive_sh.channel] + negative_sh.pid_chi2[negative_sh.channel]) ; //pchi2[channel];
			break;
	}

	//define initial four-momentum
  HepLorentzVector Pcm(cfg.CENTER_MASS_ENERGY*sin(0.011),0,0,cfg.CENTER_MASS_ENERGY); 

  fEvent.M.Mrec = (Pcm - Pkf[0] - Pkf[1]).m();

  fEvent.M.M012 = (Pkf[0]+Pkf[1]+Pkf[2]).m();
  fEvent.M.M013 = (Pkf[0]+Pkf[1]+Pkf[3]).m();
  fEvent.M.M023 = (Pkf[0]+Pkf[2]+Pkf[3]).m();
  fEvent.M.M123 = (Pkf[1]+Pkf[2]+Pkf[3]).m();

  fEvent.M.M03 =  (Pkf[0]+Pkf[3]).m();
  fEvent.M.M12 =  (Pkf[1]+Pkf[2]).m();
  fEvent.M.M01 =  (Pkf[0]+Pkf[1]).m();
  fEvent.M.M23 =  (Pkf[2]+Pkf[3]).m();


  fEvent.T.ntrack = 4;
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


  fEvent.T.ntrack=4;
  //fPid.ntrack=4;
  fMdc.T.ntrack=4;
  fDedx.ntrack=4;
  fEmc.ntrack=4;
  fTof.ntrack=4;
  fMdc.M.Mrec = get_recoil__mass(pion_pair, PION_MASS,  cfg.CENTER_MASS_ENERGY);
  //fMdc.Minv    = sqrt(get_invariant_mass2(result_pair,XMASS[channel]));
  //EvtRecTrackIterator itTrk[4] = {pion_pair.first, pion_pair.second, result_pair.first, result_pair.second};
  for(int i=0;i<4;i++)
  {
		if(Tracks[i]==evtRecTrkCol->end()) continue;
    if(i>1)
    {
      RecEmcShower *emcTrk = (*Tracks[i])->emcShower();
      fEmc.E[i] = emcTrk->energy();
      fEmc.theta[i] = emcTrk->theta();
      fEmc.phi[i] = emcTrk->phi();
      fEmc.time[i] = emcTrk->time();
      fMdc.T.E[i] = fEmc.E[i];
    }
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
    calculate_vertex(mdcTrk,rvxy,rvz,rvphi); 
    fMdc.T.vxy[i] = rvxy;
    fMdc.T.vz[i]  = rvz; 
    fMdc.T.vphi[i] = rvphi; 


    //dedx information
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

    //PID->setRecTrack((*Tracks[i]));
    //PID->calculate();
    //if(PID->IsPidInfoValid())
    //{
    //  //fPid.prob[ID_ELECTRON][i] = PID->probElectron();
    //  //fPid.prob[ID_MUON][i]     = PID->probMuon();
    //  //fPid.prob[ID_PION][i]     = PID->probPion();
    //  //fPid.prob[ID_KAON][i]     = PID->probKaon();
    //  //fPid.prob[ID_PROTON][i]   = PID->probProton();
    //}
    //vector<double> chi2 = get_chi2(Tracks[i]);
    //for(int pid=0;pid<5;pid++)
    //{
    //  //fPid.chi2[pid][i]   = chi2[pid];
    //}

  }
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
  std::cout << "Event with kaons and muons: " << event_with_kaons_and_muons << std::endl;
  std::cout << "Event with pions: " << event_with_pions << std::endl;
  std::cout << "Event with electron: " << event_with_electrons << std::endl;
  std::cout << "Event with proton: " << event_with_protons << std::endl;
  std::cout << "Good kinematic fits: " << good_kinematic_fit << std::endl;
  return StatusCode::SUCCESS;
}

template <class A>
inline StatusCode JpsiKK::init_tuple(JpsiKK * alg, A & a,  const char * dir, const char * title, MsgStream & log)
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
