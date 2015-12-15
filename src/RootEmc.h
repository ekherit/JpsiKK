// =====================================================================================
//
//       Filename:  RootEmc.h
//
//    Description:  Root adaptor for Emc
//
//        Version:  1.0
//        Created:  27.10.2015 15:23:24
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#pragma once

#include "RootTuple.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"

struct RootEmc : public RootTuple
{
	static int ARRAY_SIZE; 
	NTuple::Item<long> ntrack;
	NTuple::Array<double> E;
	NTuple::Array<double> theta;
	NTuple::Array<double> phi;
	NTuple::Array<double> time;
  NTuple::Array<double> dangle; //angle to close charged
	virtual void init(void);
	virtual void init_tuple(void);
	virtual void fill(EvtRecTrack * track) {};

	virtual void fill(int i,  EvtRecTrack * track,
		SmartDataPtr<EvtRecEvent>    & evtRecEvent, 
		SmartDataPtr<EvtRecTrackCol> & evtRecTrkCol);
      
	virtual void fill(list<EvtRecTrack*> & tracks, 
		SmartDataPtr<EvtRecEvent>    & evtRecEvent, 
		SmartDataPtr<EvtRecTrackCol> & evtRecTrkCol);
};
