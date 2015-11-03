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
	std::vector<WTrackParameter>  wtracks;
	std::vector<HepLorentzVector> P;
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
		const std::vector<EvtRecTrackIterator> & Tracks,
		KinematicFit_t & kft
    )
//bool kinematic_fit(
//		std::vector<int> & pids,  //  particle id for all supposed tracks in event 
//		const std::vector<EvtRecTrackIterator> & Tracks,  // list of registered tracks it supposed to be less then pids size
//    std::vector<HepLorentzVector> & P,     //4momentum fit result
//    double & chi2,  //chi2 result of the fit, 
//		const double W //center of mass energy
//    )
{
	std::vector<RecMdcKalTrack*> KalTrk(Tracks.size());
	std::vector<WTrackParameter> WTrk(Tracks.size());
	for(int i=0;i<Tracks.size();i++)
	{
		KalTrk[i] = (*Tracks[i])->mdcKalTrack();
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
		kft.success = false;
		return false;
	}
  bool oksq = kmfit->Fit();
  if(oksq) 
  {
    kft.chi2  = kmfit->chisq();
		kft.wtracks = kmfit->wTrackInfit();
		kft.wtracks.resize(pids.size());
		kft.P.resize(pids.size());
    for(int i=0;i<P.size();i++)
    {
      kft.P[i] = kmfit->pfit(i);
			cout << "pfit[" << i << "]" << kft.P[i].px() << " " << kft.P[i].py() << " " << kft.P[i].pz() << " " << kft.P[i].E() << endl; 
			cout << "wtrp[" << i << "]" << kmfit->wTrackInfit()[i].p().px() << " " << kmfit->wTrackInfit()[i].p().py() << " " << kmfit->wTrackInfit()[i].p().pz() << " " << kmfit->wTrackInfit()[i].p().E() << endl; 
    }
		kft.success = true;
  }
  return oksq;
}


inline std::vector<KinematicFit_t> kinfit(const std::vector<EvtRecTrackIterator> & Tracks)
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
		KinematicFit_t & k = theKFV[pid];
		k.chi2=std::numeric_limits<double>::max();
		k.success = false;
		kinematic_fit(pids, Tracks, k);
	}
	return theKFV;
}



