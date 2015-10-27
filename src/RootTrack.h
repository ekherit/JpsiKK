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


struct Track_t
{
	NTuple::Item<long>    ntrack;  //size of the array = 4: [pi-,pi+,K-,K+]
	NTuple::Array<long>   trackId; //id of the track
	NTuple::Array<double> q; //charge of the track
	NTuple::Array<double> E;
	NTuple::Array<double> p;
	NTuple::Array<double> px;
	NTuple::Array<double> py;
	NTuple::Array<double> pz;
	NTuple::Array<double> pt; //transvese momentum
	NTuple::Array<double> theta,phi;
	NTuple::Array<double> x, y, z, r; //poca coordinate of track
	NTuple::Array<double> vxy, vz, vphi; //poca coordinate of track
	virtual StatusCode add_to_tuple(NTuple::Tuple * tuple)
	{
		StatusCode status;
		status = tuple->addItem ("ntrack", ntrack,0,4); //array size must be = 4
		status = tuple->addIndexedItem ("trackId",   ntrack, trackId);
		status = tuple->addIndexedItem ("q",     ntrack, q);
		status = tuple->addIndexedItem ("E",     ntrack, E);
		status = tuple->addIndexedItem ("p",     ntrack, p);
		status = tuple->addIndexedItem ("px",    ntrack, px);
		status = tuple->addIndexedItem ("py",    ntrack, py);
		status = tuple->addIndexedItem ("pz",    ntrack, pz);
		status = tuple->addIndexedItem ("pt",    ntrack, pt);
		status = tuple->addIndexedItem ("theta", ntrack, theta);
		status = tuple->addIndexedItem ("phi",   ntrack, phi);
		status = tuple->addIndexedItem ("x",     ntrack, x);
		status = tuple->addIndexedItem ("y",     ntrack, y);
		status = tuple->addIndexedItem ("z",     ntrack, z);
		status = tuple->addIndexedItem ("r",     ntrack, r);
		status = tuple->addIndexedItem ("vxy",   ntrack, vxy);
		status = tuple->addIndexedItem ("vz",    ntrack, vz);
		status = tuple->addIndexedItem ("vphi",  ntrack, vphi);
		return status;
	}

};
