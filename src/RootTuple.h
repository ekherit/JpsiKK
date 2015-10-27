// =====================================================================================
//
//       Filename:  RootTuple.h
//
//    Description:  Interface class for root tuples
//
//        Version:  1.0
//        Created:  27.10.2015 15:24:06
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#pragma once 

#include "GaudiKernel/NTuple.h"
#include "EvtRecEvent/EvtRecTrack.h"

struct RootTuple
{
	public:
		NTuple::Tuple * tuple; //tuple
		virtual ~RootTuple(void){};
		virtual void init(void)=0;
		virtual StatusCode init_tuple(void)=0;
		virtual void fill(EvtRecTrackIterator & track) = 0;
};
