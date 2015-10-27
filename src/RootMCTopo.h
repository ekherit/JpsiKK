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
#include "EventModel/Event.h"

struct RootMCTopo : public RootTuple
{
	NTuple::Item  <int>  m_idxmc;
	NTuple::Array <int>  m_pdgid;
	NTuple::Array <int>  m_motheridx;
	NTuple::Array <int>  m_idx;
	NTuple::Item <unsigned long>  m_hash;
	virtual void init(void);
	virtual StatusCode init_tuple(void);
	virtual void fill(Event::McParticleCol *);
};
