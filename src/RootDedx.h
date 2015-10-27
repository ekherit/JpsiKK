// =====================================================================================
//
//       Filename:  RootDedx.h
//
//    Description:  Dedx Root tuple
//
//        Version:  1.0
//        Created:  27.10.2015 16:46:54
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#pragma once

#include "RootTuple.h"

struct RootDedx : public RootTuple
{
	NTuple::Item<long>    ntrack;  //size of the array = 4: [pi-,pi+,K-,K+]
	NTuple::Array<double> chie;  //chi e
	NTuple::Array<double> chimu; //chi e
	NTuple::Array<double> chipi; //chi e
	NTuple::Array<double> chik;  //chi e
	NTuple::Array<double> chip;  //chi e
	NTuple::Array<double> probPH;  //хрень какая-то
	NTuple::Array<double> normPH;
	//NTuple::Array<double> probe;  //prob e
	//NTuple::Array<double> probmu; //prob e
	//NTuple::Array<double> probpi; //prob e
	//NTuple::Array<double> probk;  //prob e
	//NTuple::Array<double> probp;  //prob e
	virtual void init(void);
	virtual void init_tuple(void);
	void fill(int i,  EvtRecTrackIterator & track);
};
