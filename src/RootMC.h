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

#include "RootTuple.h"
#include "McParticle.h"

#include "CLHEP/Vector/LorentzVector.h"

struct RootMC : public RootTuple
{
	NTuple::Item<long>    psip_decay; //is event from psip decay
	NTuple::Item<long>    jpsi_decay; //is event from jpsi decay
	NTuple::Item<long>    KK; //is event KK
	NTuple::Item<long>    uu; //is event mumu
	NTuple::Item<long>    oo; //other other particle
	NTuple::Item<long>    ntrack;  //size of the array = 4: [pi-,pi+,K-,K+]
	NTuple::Array<double> pid; //particle id
	NTuple::Array<double> q; //charge of the track
	NTuple::Array<double> E,p; //energy and momentum
	NTuple::Array<double> px,py,pz; //momentum
	NTuple::Array<double> pt; //transvese momentum
	NTuple::Array<double> theta,phi;
	virtual void init(void);
	virtual void init_tuple(void);
	virtual void fill(EvtRecTrackIterator * track){};
	virtual void fill(const std::vector<CLHEP::HepLorentzVector> & Pkf,  Event::McParticleCol *);
};
