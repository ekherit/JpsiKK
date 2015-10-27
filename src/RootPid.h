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

struct RootPid : public RootTuple
{
	NTuple::Item<double>  M[5]; //invariant mass of highmomentum track  based on Mdc
	NTuple::Item<double>  kM[5]; //invariant mass of highmomentum track  based on kinematic recons
	NTuple::Item<long>    ntrack;  //size of the array = 4: [pi-,pi+,K-,K+]
	NTuple::Array<double> prob[5]; //probability of track to be e,mu,pi,k or p
	NTuple::Array<double> chi2[5];  
	virtual void init(void);
	virtual StatusCode init_tuple(void);
};
