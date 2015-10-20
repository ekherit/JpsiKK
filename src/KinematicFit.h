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

#include "VertexFit/KinematicFit.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"

#include "PhysConst.h"
#include "Utils.h"
#include "Defs.h"

bool vertex_fit(const std::vector<WTrackParameter> & input_tracks,  std::vector<WTrackParameter> output_tracks)
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
  //vtxfit->Fit();
  vtxfit->Swim(0);

	output_tracks.resize(input_tracks.size());
	for(int i=0; i<output_tracks.size(); i++)
	{
		output_tracks[i] = vtxfit->wtrk(i);
	}
	return true;
}

/* Kinematic fit for specific pairs */
bool kinematic_fit(
    int PID,  // hypotezies 0 - kaon, 1 - muon
    TrackPair_t  &  pion_pair,  //pion pair
    TrackPair_t  & other_pair,  // other pair (kaon or muon)
    std::vector<HepLorentzVector> & P,  //4momentum fit result
    double & chi2,  //chi2 result of the fit, 
		double W
    )
{
  P.resize(4);
  EvtRecTrackIterator  PionTrk[2] = {pion_pair.first, pion_pair.second};
  EvtRecTrackIterator  OtherTrk[2] = {other_pair.first, other_pair.second};
  RecMdcKalTrack * PionKalTrk[2];
  RecMdcKalTrack * OtherKalTrk[2];
  WTrackParameter PionWTrk[2];
  WTrackParameter OtherWTrk[2];

  //For positive and negative charged pair create track parameters (WTrackParameter)
  //in different hyptoizies scpecified by PID
  for(int i=0;i<2;i++)
  {
    PionKalTrk[i] = (*PionTrk[i])->mdcKalTrack();
    OtherKalTrk[i] = (*OtherTrk[i])->mdcKalTrack();
    //pions have to be only pions
    PionWTrk[i] = WTrackParameter(PION_MASS, PionKalTrk[i]->getZHelix(), PionKalTrk[i]->getZError());
    //create other
    switch(PID)
    {
      case ID_KAON:
        OtherWTrk[i] = WTrackParameter(XMASS[PID], OtherKalTrk[i]->getZHelixK(), OtherKalTrk[i]->getZErrorK());
        break;
      case ID_MUON:
        OtherWTrk[i] = WTrackParameter(XMASS[PID], OtherKalTrk[i]->getZHelixMu(), OtherKalTrk[i]->getZErrorMu());
        break;
      case ID_ELECTRON:
        OtherWTrk[i] = WTrackParameter(XMASS[PID], OtherKalTrk[i]->getZHelixE(), OtherKalTrk[i]->getZErrorE());
        break;
      case ID_PION:
        OtherWTrk[i] = WTrackParameter(XMASS[PID], OtherKalTrk[i]->getZHelix(), OtherKalTrk[i]->getZError());
        break;
      case ID_PROTON:
        OtherWTrk[i] = WTrackParameter(XMASS[PID], OtherKalTrk[i]->getZHelixP(), OtherKalTrk[i]->getZErrorP());
        break;
    }
  }
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
  //add tracks. I know the first two one must be pions
  vtxfit->AddTrack(0,  PionWTrk[0]);
  vtxfit->AddTrack(1,  PionWTrk[1]);
  vtxfit->AddTrack(2,  OtherWTrk[0]);
  vtxfit->AddTrack(3,  OtherWTrk[1]);
  vtxfit->AddVertex(0, vxpar,0, 1, 2,3);
  if(!vtxfit->Fit(0)) return false;
  vtxfit->Fit();
  vtxfit->Swim(0);

  KalmanKinematicFit * kmfit = KalmanKinematicFit::instance();
  //kmfit->setIterNumber(10000);
  //kmfit->setChisqCut(10000);

  kmfit->init();
  for(int i=0;i<4;i++)
  {
    kmfit->AddTrack(i,vtxfit->wtrk(i));
  }
  HepLorentzVector Pcmf(W*sin(0.011) /* 40.546 MeV*/,0,0,W); //initial vector of center of mass frame
  kmfit->AddFourMomentum(0,  Pcmf);
  //kmfit->AddTotalEnergy(0,PSIP_MASS,0,1,2,3);
  //kmfit->AddResonance(1, JPSI_MASS, 2, 3);
  //kmfit->AddResonance(1, PSIP_MASS, 0, 1, 2,3);
  if(!kmfit->Fit(0)) return false;
  bool oksq = kmfit->Fit();
  if(oksq) 
  {
    chi2  = kmfit->chisq();
    for(int i=0;i<4;i++)
    {
      P[i] = kmfit->pfit(i);
    }
  }
  return oksq;
}




bool kinematic_fit(
    int PID, //hypotizes 0 -- KAONS, 1 - MUONS
    TrackPairList_t  & pion_pairs,  //list of pion pairs  (signle pair in simple case)
    TrackPairList_t &  other_pairs, //list of other pairs  (signle pair in simple case)
    std::vector<HepLorentzVector> & P, //?
    double & chi2,  //result of the fit
    TrackPair_t & result_pion_pair,  //best pion pair
    TrackPair_t & result_other_pair,   //best other pair
		double W
    )
{
  chi2=std::numeric_limits<double>::max();
  if(pion_pairs.empty() || other_pairs.empty()) return false;
  bool GoodKinematikFit=false;
  //loop over all pion pairs and all other pairs
  for(TrackPairList_t::iterator pion_pair=pion_pairs.begin(); pion_pair!=pion_pairs.end();pion_pair++)
    for(TrackPairList_t::iterator other_pair=other_pairs.begin(); other_pair!=other_pairs.end();other_pair++)
    {
      std::vector<HepLorentzVector> P_tmp;
      double chi2_tmp=std::numeric_limits<double>::max();
      bool oksq=kinematic_fit(PID, *pion_pair, *other_pair, P_tmp, chi2_tmp, W);
      if(oksq) 
      {
        if(chi2_tmp < chi2)
        {
          GoodKinematikFit = true;
          chi2  = chi2_tmp;
          P = P_tmp;
          result_pion_pair = *pion_pair;
          result_other_pair = * other_pair;
        }
      }
    } 
  return GoodKinematikFit;
}



bool kinematic_fit(
    int PID,  // hypotezies 0 - kaon, 1 - muon,  ...
		const std::vector<EvtRecTrackIterator> & Tracks, 
		/* Tracks[0] - pi-
		 * Tracks[1] - pi+
		 * Tracks[2] - K/mu-
		 * Tracks[3] - K/mu+
		 */
    std::vector<HepLorentzVector> & P,     //4momentum fit result
    double & chi2,  //chi2 result of the fit, 
		double W
    )
{
	std::vector<RecMdcKalTrack*> KalTrk(Tracks.size());
	std::vector<WTrackParameter> WTrk(Tracks.size());
	for(int i=0;i<Tracks.size();i++)
	{
		KalTrk[i] = (*Tracks[i])->mdcKalTrack();
		if(i<2) WTrk[i] = WTrackParameter(PION_MASS, KalTrk[i]->getZHelix(), KalTrk[i]->getZError());
		else
		{
			switch(PID)
			{
				case ID_KAON:
					WTrk[i] = WTrackParameter(XMASS[PID], KalTrk[i]->getZHelixK(),  KalTrk[i]->getZErrorK());
					break;
				case ID_MUON:
					WTrk[i] = WTrackParameter(XMASS[PID], KalTrk[i]->getZHelixMu(), KalTrk[i]->getZErrorMu());
					break;
				case ID_ELECTRON:
					WTrk[i] = WTrackParameter(XMASS[PID], KalTrk[i]->getZHelixE(),  KalTrk[i]->getZErrorE());
					break;
				case ID_PION:
					WTrk[i] = WTrackParameter(XMASS[PID], KalTrk[i]->getZHelix(),   KalTrk[i]->getZError());
					break;
				case ID_PROTON:
					WTrk[i] = WTrackParameter(XMASS[PID], KalTrk[i]->getZHelixP(),  KalTrk[i]->getZErrorP());
					break;
			}
		}
	}

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
	for(int i=0;i<WTrk.size();i++)
	{
		vtxfit->AddTrack(i, WTrk[i]);
	}
	if(WTrk.size()<4) vtxfit->AddMissTrack(3,XMASS[PID]);
  //add tracks. I know the first two one must be pions
  vtxfit->AddVertex(0, vxpar,0, 1, 2, 3);
  if(!vtxfit->Fit(0)) return false;
  vtxfit->Fit();
  vtxfit->Swim(0);

  KalmanKinematicFit * kmfit = KalmanKinematicFit::instance();
  //kmfit->setIterNumber(10000);
  //kmfit->setChisqCut(10000);

  kmfit->init();
  for(int i=0;i<4;i++)
  {
    kmfit->AddTrack(i,vtxfit->wtrk(i));
  }
  HepLorentzVector Pcmf(W*sin(0.011) /* 40.546 MeV*/,0,0,W); //initial vector of center of mass frame
  kmfit->AddFourMomentum(0,  Pcmf);
  if(!kmfit->Fit(0)) return false;
  bool oksq = kmfit->Fit();
  if(oksq) 
  {
    chi2  = kmfit->chisq();
		P.resize(4);
    for(int i=0;i<P.size();i++)
    {
      P[i] = kmfit->pfit(i);
    }
  }
  return oksq;
}


bool kinfit(
		const std::vector<EvtRecTrackIterator> & Tracks,  
		int & channel,  
		double & chi2,  
		std::vector<HepLorentzVector> & P,  
		const double CENTER_MASS_ENERGY
		)
{
	bool goodfit=false;
	for(int pid=0;pid<2;pid++)
	{
		double chi2_tmp=1e100;
		std::vector<HepLorentzVector> P_tmp;
		std::cout << "Before kinematic_fit" << endl;
		bool fit_result = kinematic_fit(pid, Tracks, P_tmp, chi2_tmp, CENTER_MASS_ENERGY);
		if(fit_result)
		{
			goodfit=true;
			if(chi2_tmp<chi2)
			{
				channel = pid;
				chi2 = chi2_tmp;
				P = P_tmp;
			}
		}
	}
	return goodfit;
}

bool kinfit(SelectionHelper_t & kfp)
{
	kfp.good_kinematic_fit = kinfit(kfp.tracks,  kfp.channel,  kfp.kin_chi2,  kfp.P,  kfp.W);
	return kfp.good_kinematic_fit;
}


bool kinfit(
		TrackPair_t & pion_pair,
		TrackList_t & other_tracks, 
		SelectionHelper_t & kfp
		)
{
	SelectionHelper_t tmp_kfp(kfp);
	tmp_kfp.tracks.resize(3);
	tmp_kfp.tracks[0]=pion_pair.first;
	tmp_kfp.tracks[1]=pion_pair.second;
	for(TrackList_t::iterator i=kfp.tracks.begin(); i!=kfp.tracks.end(); ++i)
	{
		EvtRecTrackIterator track = *i;
		tmp_kfp.tracks[2] = track;
		if(kinfit(tmp_kfp))
		{
			kfp.good_kinematic_fit = true;
			if(tmp_kfp.kin_chi2 < kfp.kin_chi2)
			{
				kfp = tmp_kfp;
			}
		}
	}
	//kfp.tracks.push_back(kfp.tracks[2]);
	return  kfp.good_kinematic_fit;
}
