// =====================================================================================
//
//       Filename:  KinematicFit.h
//
//    Description:  
//
//        Version:  1.0
//        Created:  20.10.2015 11:22:00
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#pragma once
#include <iostream>

#include "VertexFit/KinematicFit.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"

#include "PhysConst.h"
#include "Utils.h"
#include "Defs.h"

#include <limits>

struct KinematicFit_t
{
	bool success;
	double chi2;
	std::vector<HepLorentzVector> P;
	KinematicFit_t(void)
	{
		chi2=std::numeric_limits<double>::max();
		success = false;
	}
};

bool vertex_fit(const std::vector<WTrackParameter> & input_tracks,  std::vector<WTrackParameter> & output_tracks)
{
  //vertex fit - уточним вершины
  //initial vertex
  HepPoint3D vx(0., 0., 0.);
  HepSymMatrix Evx(3, 0);
  double bx = 1E+6;
  double by = 1E+6;
  double bz = 1E+6;
  Evx[0][0] = bx*bx;
  Evx[1][1] = by*by;
  Evx[2][2] = bz*bz;
  VertexParameter vxpar;
  vxpar.setVx(vx);
  vxpar.setEvx(Evx);

  //Vetex fitter
  VertexFit* vtxfit = VertexFit::instance();
  vtxfit->init();
	std::vector<int> index_list(input_tracks.size());
	for(int i=0; i<input_tracks.size(); i++)
	{
		vtxfit->AddTrack(i, input_tracks[i]);
		index_list[i] = i;
	}
  vtxfit->AddVertex(0, vxpar,index_list);
  if(!vtxfit->Fit()) return false;
  vtxfit->Swim(0);

	output_tracks.resize(input_tracks.size());
	for(int i=0; i<output_tracks.size(); i++)
	{
		output_tracks[i] = vtxfit->wtrk(i);
	}
	return true;
}




bool kinematic_fit(
		std::vector<int> & pids,  //  particle id for all supposed tracks in event 
		const std::vector<EvtRecTrack*> & Tracks,
		KinematicFit_t & kft
    )
{
  kft.P.resize(pids.size());
  kft.success = false;
	std::vector<RecMdcKalTrack*> KalTrk(Tracks.size());
	std::vector<WTrackParameter> WTrk(Tracks.size());
	for(int i=0;i<Tracks.size();i++)
	{
		KalTrk[i] = Tracks[i]->mdcKalTrack();
		switch(pids[i])
		{
			case ID_KAON:
				WTrk[i] = WTrackParameter(XMASS[pids[i]], KalTrk[i]->getZHelixK(),  KalTrk[i]->getZErrorK());
				break;
			case ID_MUON:
				WTrk[i] = WTrackParameter(XMASS[pids[i]], KalTrk[i]->getZHelixMu(), KalTrk[i]->getZErrorMu());
				break;
			case ID_ELECTRON:
				WTrk[i] = WTrackParameter(XMASS[pids[i]], KalTrk[i]->getZHelixE(),  KalTrk[i]->getZErrorE());
				break;
			case ID_PION:
				WTrk[i] = WTrackParameter(XMASS[pids[i]], KalTrk[i]->getZHelix(),   KalTrk[i]->getZError());
				break;
			case ID_PROTON:
				WTrk[i] = WTrackParameter(XMASS[pids[i]], KalTrk[i]->getZHelixP(),  KalTrk[i]->getZErrorP());
				break;
		}
	}

  //vertex fit - уточним вершины
	std::vector<WTrackParameter> VertexWTrk;
	if(!vertex_fit(WTrk, VertexWTrk)) return false;

  KalmanKinematicFit * kmfit = KalmanKinematicFit::instance();
  kmfit->init();
  //kmfit->setChisqCut(1000);
  for(int i=0;i<WTrk.size();i++)
  {
    kmfit->AddTrack(i,VertexWTrk[i]);
  }
	//add missed tracks
	for(int i=WTrk.size(); i<pids.size(); i++)
	{
		kmfit->AddMissTrack(i,XMASS[pids[i]]);
	}
	//total momentum of tracks
  kmfit->AddFourMomentum(0,  getTotalMomentum());
  if(!kmfit->Fit(0)) 
	{
		return false;
	}
  bool oksq = kmfit->Fit();
  if(oksq) 
  {
    kft.chi2  = kmfit->chisq();
		//kft.wtracks = kmfit->wTrackInfit();
    for(int i=0;i<kft.P.size();i++)
    {
      kft.P[i] = kmfit->pfit(i);
    }
		kft.success = true;
  }
  return oksq;
}


inline std::vector<KinematicFit_t> kinfit(const std::vector<EvtRecTrack*> & Tracks)
{
		/* Tracks[0] - pi-
		 * Tracks[1] - pi+
		 * Tracks[2] - K/mu-
		 * Tracks[3] - K/mu+
		 */
	std::vector<KinematicFit_t> theKFV(5); //5 hypotesis
	std::vector<int> pids(4);
	pids[0] = ID_PION;
	pids[1] = ID_PION;
	for(int pid=0;pid<theKFV.size(); pid++)
	{
		pids[2] = pid;
		pids[3] = pid;
		kinematic_fit(pids, Tracks, theKFV[pid]);
	}
	return theKFV;
}



inline double kinfit_3pi(
    TrackVector_t & Tq,  //charged tracks
    TrackList_t &  T0,   //neutral tracks
    double & chi2)       //list of all neutral tracks
{
  chi2=300;
  double Mpi0=-10; //best pi0 mass
	std::vector<RecMdcKalTrack*> KalTrk(Tq.size());
	std::vector<WTrackParameter> WTrk(Tq.size());
	for(int i=0;i<Tq.size();i++)
	{
		KalTrk[i] = Tq[i]->mdcKalTrack();
    WTrk[i] = WTrackParameter(XMASS[ID_KAON], KalTrk[i]->getZHelix(),   KalTrk[i]->getZError());
	}
	std::vector<WTrackParameter> VertexWTrk;
	if(!vertex_fit(WTrk, VertexWTrk)) return Mpi0;

  KalmanKinematicFit * kmfit = KalmanKinematicFit::instance();
  if(T0.empty()) return Mpi0;

  //now loop over  neutral tracks and find best
  for(TrackList_t::iterator it1 = T0.begin() ; it1 != T0.end() ; it1++)
  {
    TrackList_t::iterator it2=it1;
    it2++;
    for(it2; it2 != T0.end() ; it2++)
    {
      kmfit->init();
      for(int i=0;i<Tq.size();i++)
      {
        kmfit->AddTrack(i,VertexWTrk[i]);
      }
      if(Tq.size()==3)
      {
        kmfit->AddMissTrack(3,XMASS[ID_PION]);
      }

      HepLorentzVector Pg[2]; //photon four-momentum
      RecEmcShower * emcTrk1=(*it1)->emcShower();
      RecEmcShower * emcTrk2=(*it2)->emcShower();

      kmfit->AddTrack(4,0,emcTrk1);
      kmfit->AddTrack(5,0,emcTrk2);

      kmfit->AddResonance(0,0.1349766, 4,5); //pi0 particle
      //kmfit->AddResonance(1,JPSI_MASS, 2,3,4,5); //jpsi particle
      kmfit->AddFourMomentum(1,  getTotalMomentum()); //total momeunum
      //if(!kmfit->Fit(0)) continue;
      //if(!kmfit->Fit(1)) continue;
      //if(!kmfit->Fit(2)) continue;
      bool oksq = kmfit->Fit();
      std::cout << "E1 = " << emcTrk1->energy() << " E2=" << emcTrk2->energy() << " oksq=" << oksq << " " << kmfit->shisq() << " " ;
      if(oksq)
      {
        if(kmfit->chisq() < chi2)
        {
          chi2 =  kmfit->chisq();
          Pg[0] = kmfit->pfit(4);
          Pg[1] = kmfit->pfit(5);
          Mpi0 = (Pg[0]+Pg[1]).m();
          std::cout << " Mpi0=" << Mpi0 << std::endl;
        }
      }
    }
  }
  return Mpi0;
}



