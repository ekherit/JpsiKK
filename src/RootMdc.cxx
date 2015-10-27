// =====================================================================================
//
//       Filename:  RootMdc.cxx
//
//    Description:  
//
//        Version:  1.0
//        Created:  27.10.2015 17:02:36
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#include "RootMdc.h"
#include "Utils.h"

void RootMdc::init_tuple(void)
{
	M.add_to_tuple(tuple);
	T.add_to_tuple(tuple);
}


void RootMdc::init(void)
{
  T.ntrack=4;
}


void RootMdc::fill(int i, EvtRecTrackIterator & track)
{
	if(!(*track)->isMdcTrackValid()) return; 
	RecMdcTrack  *mdcTrk = (*track)->mdcTrack();
	T.trackId[i] = mdcTrk->trackId();
	T.q[i] = mdcTrk->charge(); 
	T.p[i] = mdcTrk->p();
	T.px[i]= mdcTrk->px();
	T.py[i]= mdcTrk->py();
	T.pz[i]= mdcTrk->pz();
	T.theta[i]= mdcTrk->theta();
	T.phi[i] = mdcTrk->phi();
	T.x[i]  = mdcTrk->x();
	T.y[i]  = mdcTrk->y();
	T.z[i]  = mdcTrk->z();
	T.x[i]  = mdcTrk->x();
	T.y[i]  = mdcTrk->y();
	T.z[i]  = mdcTrk->z();
	double rvxy,rvz,rvphi;
	calculate_vertex(mdcTrk,rvxy,rvz,rvphi); 
	T.vxy[i] = rvxy;
	T.vz[i]  = rvz; 
	T.vphi[i] = rvphi; 

	if((*track)->isEmcShowerValid())
	{
		RecEmcShower *emcTrk = (*track)->emcShower();
		T.E[i] = emcTrk->energy();
	}
	else
	{
		T.E[i] = 0;
	}
}
