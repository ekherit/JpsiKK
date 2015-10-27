// =====================================================================================
//
//       Filename:  RootTof.h
//
//    Description:  TOF tuple
//
//        Version:  1.0
//        Created:  27.10.2015 16:49:00
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#pragma once

#include "RootTuple.h"

struct RootTof : public RootTuple
{
	NTuple::Item<long> ntrack;
	NTuple::Array<double> tofID;
	NTuple::Array<double> t0;
	NTuple::Array<double> t; //tof time
	NTuple::Array<double> dt; //error of tof time
	NTuple::Array<double> beta;  
	NTuple::Array<double> te;  //electron expected time
	NTuple::Array<double> tmu; //muon
	NTuple::Array<double> tpi; //pion
	NTuple::Array<double> tk;  //kaon
	NTuple::Array<double> tp;  //proton

	NTuple::Array<double> chie;  
	NTuple::Array<double> chimu; 
	NTuple::Array<double> chipi; 
	NTuple::Array<double> chik;  
	NTuple::Array<double> chip;  
	virtual void init(void);
	virtual StatusCode init_tuple(void);
	virtual void fill(EvtRecTrackIterator & track);
};


