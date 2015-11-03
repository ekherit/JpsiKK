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
#include "RootTrack.h"
#include "RootMass.h"

#include "Defs.h"

struct RootMdc : public RootTuple
{
	//Mass_t M;
	Track_t T;

	NTuple::Item<double> Mrec;
	NTuple::Item<long> npid;
	NTuple::Array<double>  M23;
	NTuple::Array<double>  M12;
	NTuple::Array<double>  M03;
	NTuple::Array<double>  Mmis;

	virtual void init(void);
	virtual void init_tuple(void);
	virtual void fill(int i,  EvtRecTrackIterator & track);
	virtual void fill_mass(TrackVector_t&, EvtRecTrackIterator&);
};
