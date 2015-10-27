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

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;

#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/Helix.h"

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"

#include "utils.h"
#include "Defs.h"
#include "SelectionConfig.h"

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



inline double get_invariant_mass2(TrackPair_t & pair, double mass)
{
  EvtRecTrackIterator  itTrk[2] = {pair.first, pair.second};
  HepLorentzVector  P[2];
  for(int k=0;k<2;k++)
  {
    if(!(*itTrk[k])->isMdcTrackValid()) break; 
    RecMdcTrack *mdcTrk = (*itTrk[k])->mdcTrack();
    P[k] = mdcTrk->p4(mass);
  }
  HepLorentzVector P_sum = P[0]+P[1];
  return P_sum.m2();
}

inline double get_recoil__mass(EvtRecTrackIterator & trk1, EvtRecTrackIterator & trk2, double mass,  double W /*  total energy */)
{
  EvtRecTrackIterator  itTrk[2] = {trk1, trk2};
  HepLorentzVector  P[2];
  for(int k=0;k<2;k++)
  {
    if(!(*itTrk[k])->isMdcTrackValid()) break; 
    RecMdcTrack *mdcTrk = (*itTrk[k])->mdcTrack();
    P[k] = mdcTrk->p4(mass);
  }
  HepLorentzVector P_psip(W*sin(0.011),0,0,W); //initial vector of psip
  HepLorentzVector P_sum = P[0]+P[1];
  HepLorentzVector P_recoil = P_psip - P_sum;
  return P_recoil.m();
}

inline double get_recoil__mass(std::pair<EvtRecTrackIterator, EvtRecTrackIterator> p, double mass,  double W)
{
  return get_recoil__mass(p.first, p.second, mass,  W);
}


inline double get_missing_mass(std::pair<EvtRecTrackIterator, EvtRecTrackIterator> pions, std::pair<EvtRecTrackIterator, EvtRecTrackIterator> kaons,  double W)
{
  EvtRecTrackIterator  PionTrk[2] = {pions.first, pions.second};
  EvtRecTrackIterator  KaonTrk[2] = {kaons.first, kaons.second};
  HepLorentzVector P_psip(W*sin(0.011),0,0,W); //initial vector of psip
  HepLorentzVector  pionP[2];
  HepLorentzVector  kaonP[2];
  for(int k=0;k<2;k++)
  {
    pionP[k] = HepLorentzVector(0,0,0,0);
    if(!(*PionTrk[k])->isMdcTrackValid()) break; 
    RecMdcTrack *mdcTrk = (*PionTrk[k])->mdcTrack();
    pionP[k] = mdcTrk->p4(PION_MASS);
  }
  for(int k=0;k<2;k++)
  {
    kaonP[k] = HepLorentzVector(0,0,0,0);
    if(!(*KaonTrk[k])->isMdcTrackValid()) break; 
    RecMdcTrack *mdcTrk = (*KaonTrk[k])->mdcTrack();
    kaonP[k] = mdcTrk->p4(KAON_MASS);
  }
  HepLorentzVector Pmis = P_psip - pionP[0] - pionP[1] - kaonP[0] - kaonP[1];
  return Pmis.m2();
}
