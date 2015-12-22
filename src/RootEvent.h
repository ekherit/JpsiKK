// =====================================================================================
//
//       Filename:  RootMCTopo.h
//
//    Description:  
//
//        Version:  1.0
//        Created:  27.10.2015 16:52:08
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#pragma once

#include <vector>

#include "RootTuple.h"
#include "RootMass.h"
#include "RootTrack.h"

#include "CLHEP/Vector/LorentzVector.h"

// =====================================================================================
//        Class:  RootEvent
//  Description:  Main event information supposed to used for selection
// =====================================================================================
struct RootEvent : public RootTuple
{
	NTuple::Item<long>    run; //run number
	NTuple::Item<long>    event; //event number 
	NTuple::Item<long>    time; //time of the event
	NTuple::Item<long>    ngood_charged_track;     //number of good charged tracks in event
	NTuple::Item<long>    ngood_neutral_track;     //number of good neutral tracks in event
	NTuple::Item<long>    npositive_track; //number of positive charged tracks
	NTuple::Item<long>    nnegative_track; //number of positive charged tracks
	NTuple::Item<long>    npositive_pions; //number of positive pions
	NTuple::Item<long>    nnegative_pions; //number of negative pions
	NTuple::Item<long>    npion_pairs; //total number of found pion pairs

	NTuple::Item<long>    sign;   //signature of the event shows missed particle  01(K- or mu-), 10 (K+ or mu+),  11 (KK or uu or Ku or uK)
	NTuple::Item<long>    channel;     //J/psi decay channel 0 -- kaons, 1 -- muons,  10 -- muK,  11 - Kmu
	NTuple::Item<long>    K;     //KK Jpsi decay event
	NTuple::Item<long>    u;     //KK Jpsi decay event
	NTuple::Item<long>    KK;     //KK Jpsi decay event
	NTuple::Item<long>    uu;     //MuMu event


	Mass_t  M;

	//here will be result of the kinematic fit and particle id
	NTuple::Item<double>  kin_chi2;    //kinematic chi2
	NTuple::Item<double>  pid_chi2;    //ParticleID chi2
	NTuple::Item<double>  pid_prob;    //ParticleID probability
	NTuple::Item<double>  mypid_chi2; //My particle id chi2

	Track_t T;

	NTuple::Item<long>    npid;      //number of particle id's
	NTuple::Array<double> kchi;      //kinematik chi2
	NTuple::Array<double> pchi;      //particle id chi square
	NTuple::Array<double> prob;      //probability from ParticleID
	NTuple::Array<double> mypchi;      //probability from ParticleID

	NTuple::Array<double> kM23;      //Invariant mass for track 23 for different hyptotesis

  NTuple::Item<double> Mpi0; //invariant mass of two gamma most close to pi0 mezon
  //NTuple::Item<double> kin_chi2_3pi; //invariant mass of two gamma most close to pi0 mezon

	NTuple::Item<long>    nkinbg;    //number background kinematic fits
  NTuple::Array<double> kin_chi2_bg; //kinematic fit for the background
  /*
   * 0 -  Ψ(2S) → π+π-(J/Ψ → π-(ρ(770)+ → (π0 → ɣɣ)π+)) = ɣɣπ+π+π-π- (final state)
   * 1 -  same and missed one gamma
   * 
   * */

	virtual void init(void);
	virtual void init_tuple(void);

	virtual void fill(const std::vector<HepLorentzVector> & Pkf);
	virtual void fill(int i,  EvtRecTrack * track);
};
