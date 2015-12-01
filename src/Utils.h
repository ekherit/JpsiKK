// =====================================================================================
//
//       Filename:  Utils.h
//
//    Description:  Some usefull functions
//
//        Version:  1.0
//        Created:  19.10.2015 21:12:17
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================


#pragma once

#include <exception>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;

//#include "GaudiKernel/IDataProviderSvc.h"
//#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/Bootstrap.h"

#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/Helix.h"

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"

#include "utils.h"
#include "Defs.h"
#include "SelectionConfig.h"
#include "PhysConst.h"


inline HepLorentzVector getTotalMomentum(double Wcm = BEAM_CENTER_MASS_ENERGY)
{
	return HepLorentzVector(Wcm*tan(0.5*BEPC_CROSSING_ANGLE),0,0,Wcm/cos(0.5*BEPC_CROSSING_ANGLE));
	//I think formula below is wrong 
	//return HepLorentzVector(Wcm*sin(0.5*BEPC_CROSSING_ANGLE),0,0,Wcm);
}

inline void calculate_vertex(RecMdcTrack *mdcTrk, double & ro, double  & z, double phi)
{
  ro = -9999;
  z = -9999;
  phi = -9999;
  /*  Reconstruct the vertex */
  Hep3Vector xorigin(0,0,0);
  IVertexDbSvc*  vtxsvc;
  Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
  if(vtxsvc->isVertexValid())
  {
    double* dbv = vtxsvc->PrimaryVertex(); 
    double*  vv = vtxsvc->SigmaPrimaryVertex();  
    xorigin.setX(dbv[0]);
    xorigin.setY(dbv[1]);
    xorigin.setZ(dbv[2]);
  }
  /* Vertex game. copy from rhophi analysis */
  double phi0=mdcTrk->helix(1);
  double xv=xorigin.x();
  double yv=xorigin.y();
  //double Rxy=(mdc.x[i]-xv)*cos(phi0)+(mdc.y[i]-yv)*sin(phi0);
  //mdc.r[i]=Rxy;
  HepVector a = mdcTrk->helix();
  HepSymMatrix Ea = mdcTrk->err();
  HepPoint3D point0(0.,0.,0.);   // the initial point for MDC recosntruction
  HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]); 
  VFHelix helixip(point0,a,Ea); 
  helixip.pivot(IP);
  HepVector vecipa = helixip.a();
  double  Rvxy0=fabs(vecipa[0]);  //the nearest distance to IP in xy plane
  double  Rvz0=vecipa[3];         //the nearest distance to IP in z direction
  double  Rvphi0=vecipa[1];
  ro=Rvxy0;
  z=Rvz0;
  phi=Rvphi0;
}


inline std::list<EvtRecTrackIterator> createGoodChargedTrackList(
		SelectionConfig & cfg, 
		SmartDataPtr<EvtRecEvent>    & evtRecEvent, 
		SmartDataPtr<EvtRecTrackCol> & evtRecTrkCol
		)
{
  std::list<EvtRecTrackIterator> good_charged_tracks;
  for(unsigned i = 0; i < evtRecEvent->totalCharged(); i++)
  {
    EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
    if(!(*itTrk)->isMdcTrackValid()) continue;  //use only valid charged tracks
    RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();  //main drift chambe
    //calculate interaction point distance
    double rvxy,rvz,rvphi;
    calculate_vertex(mdcTrk,rvxy,rvz,rvphi); //find distance to interaction point
    bool IP_track = fabs(rvz)< cfg.IP_MAX_Z && fabs(rvxy)<cfg.IP_MAX_RHO;  //tracks begin near interaction point
    bool good_track = IP_track && fabs(cos(mdcTrk->theta()))<cfg.MAX_COS_THETA; //track is good
    if(good_track) good_charged_tracks.push_back(itTrk);
  }
	return good_charged_tracks;
}

inline std::list<EvtRecTrackIterator> createGoodNeutralTrackList(
		SelectionConfig & cfg, 
		SmartDataPtr<EvtRecEvent>    & evtRecEvent, 
		SmartDataPtr<EvtRecTrackCol> & evtRecTrkCol
		)
{
	std::list<EvtRecTrackIterator> good_neutral_tracks;
	//collect good neutral track
	for(int i = evtRecEvent->totalCharged(); i<evtRecEvent->totalTracks(); i++)
	{
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		if(!(*itTrk)->isEmcShowerValid()) continue; //keep only valid neutral tracks
		RecEmcShower *emcTrk = (*itTrk)->emcShower();
		double c =  fabs(cos(emcTrk->theta())); //abs cos theta
		double E  =  emcTrk->energy();
		bool hit_barrel = (c <= cfg.EMC_BARREL_MAX_COS_THETA);
		bool hit_endcup = (cfg.EMC_ENDCUP_MIN_COS_THETA <=c) && (c <= cfg.EMC_ENDCUP_MAX_COS_THETA);
		//barrel and endcup calorimeters have different energy threshold
		bool barrel_good_track = hit_barrel && (E > cfg.EMC_BARREL_MIN_ENERGY);
		bool endcup_good_track = hit_endcup && (E > cfg.EMC_ENDCUP_MIN_ENERGY);
		if(barrel_good_track  || endcup_good_track) 
		{
			//cout << "Energy of good neutral track: " << E << endl;
			good_neutral_tracks.push_back(itTrk);
		}
	}
	return good_neutral_tracks;
}




inline double getMissingMass2(double Wcm,  TrackVector_t & T, std::vector<int> & pid)
{
	HepLorentzVector Ptotal = getTotalMomentum(Wcm);
	std::vector<HepLorentzVector> P(T.size());
	HepLorentzVector Psum;
	for(int i=0; i<T.size(); i++)
	{
		if(!(*T[i])->isMdcTrackValid()) throw std::runtime_error("Bad track at calculating missing mass (getMissingMass2)");
    //RecMdcTrack *mdcTrk = (*T[i])->mdcTrack();
    RecMdcKalTrack *mdcTrk = (*T[i])->mdcKalTrack();
		P[i] =  mdcTrk->p4(XMASS[pid[i]]);
		Psum+=P[i];
	}
  HepLorentzVector Pmis = Ptotal - Psum;
  return Pmis.m2();
}

inline double getMissingMass2(TrackVector_t & T, std::vector<int> & pid)
{
	return getMissingMass2(BEAM_CENTER_MASS_ENERGY, T, pid);
}

inline double getInvariantMass2(TrackVector_t & T, std::vector<int> & pid)
{
	HepLorentzVector Psum;
	for(int i=0; i<T.size(); i++)
	{
		if(!(*T[i])->isMdcTrackValid()) throw std::runtime_error("Bad track at calculating invariant mass (getInvariantMass2)");
    //RecMdcTrack *mdcTrk = (*T[i])->mdcTrack();
    RecMdcKalTrack *mdcTrk = (*T[i])->mdcKalTrack();
		Psum += (mdcTrk->p4(XMASS[pid[i]]));;
	}
  return Psum.m2();
}

inline double getInvariantMass2(int pid1,  EvtRecTrackIterator & t1,  int pid2, EvtRecTrackIterator & t2)
{
	std::vector<int> pids(2);
	pids[0] = pid1;
	pids[1] = pid2;
	TrackVector_t T(2);
	T[0] = t1;
	T[1] = t2;
	return getInvariantMass2(T, pids);
}

inline double getInvariantMass2(int pid, EvtRecTrackIterator & track, const HepLorentzVector & v)
{
  if(!(*track)->isMdcTrackValid()) throw std::runtime_error("Bad track at calculating invariant mass (getInvariantMass2)");
  //RecMdcTrack *mdcTrk = (*track)->mdcTrack();
  RecMdcKalTrack *mdcTrk = (*track)->mdcKalTrack();
  return (v + mdcTrk->p4(XMASS[pid])).m2();
}





inline double getPionRecoilMass(double Wcm,  EvtRecTrackIterator & t1,  EvtRecTrackIterator & t2)
{
	std::vector<int> pid(2, ID_PION);
	TrackVector_t T(2);
	T[0] = t1;
	T[1] = t2;
	return sqrt(getMissingMass2(Wcm, T,  pid));
}

inline double getPionRecoilMass(EvtRecTrackIterator & t1,  EvtRecTrackIterator & t2)
{
	std::vector<int> pid(2, ID_PION);
	TrackVector_t T(2);
	T[0] = t1;
	T[1] = t2;
	return sqrt(getMissingMass2(T,  pid));
}

RecMdcTrack * getMdcTrack(EvtRecTrackIterator & track)
{
  EvtRecTrackIterator & itTrk = *track;
  if(!(*itTrk)->isMdcTrackValid()) throw std::runtime_error("No MDC info"); 
  return (*itTrk)->mdcTrack();
}
