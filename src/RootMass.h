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

#include "GaudiKernel/NTuple.h"

struct Mass_t
{
	NTuple::Item<double>  Mrec;  //pion recoil mass
	NTuple::Item<double>  M012; //M(pi pi K/mu-)
	NTuple::Item<double>  M013; //M(pi pi K/mu+)
	NTuple::Item<double>  M023; //M(pi- KK/uu)
	NTuple::Item<double>  M123; //M(pi+ KK/uu)

	NTuple::Item<double>  M03; //invariant mass of Kaon and pion
	NTuple::Item<double>  M12; //invariant mass of kaon and pion
	NTuple::Item<double>  M01; //invariant mass of pion
	NTuple::Item<double>  M23; //invariant mass of kaons or muons

	NTuple::Array<double> Mmis;    //missing invariant mass

	virtual void add_to_tuple(NTuple::Tuple * tuple)
	{
		tuple->addItem ("Mrec", Mrec); 
		tuple->addItem ("M012", M012); 
		tuple->addItem ("M013", M013); 
		tuple->addItem ("M023", M023); 
		tuple->addItem ("M123", M123); 

		tuple->addItem ("M03", M03); 
		tuple->addItem ("M12", M12); 
		tuple->addItem ("M01", M01); 
		tuple->addItem ("M23", M23); 
	}
};
