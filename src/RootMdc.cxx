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

StatusCode RootMdc::init_tuple(void)
{
  StatusCode status;
	status = M.add_to_tuple(tuple);
	status = T.add_to_tuple(tuple);
  return status;
}


void RootMdc::init(void)
{
  T.ntrack=4;
}


void RootMdc::fill(EvtRecTrackIterator & track)
{
	RecMdcTrack  *mdcTrk = (*track)->mdcTrack();
	fMdc.T.trackId[i] = mdcTrk->trackId();
	fMdc.T.q[i] = mdcTrk->charge(); 
	fMdc.T.p[i] = mdcTrk->p();
	fMdc.T.px[i]= mdcTrk->px();
	fMdc.T.py[i]= mdcTrk->py();
	fMdc.T.pz[i]= mdcTrk->pz();
	fMdc.T.theta[i]= mdcTrk->theta();
	fMdc.T.phi[i] = mdcTrk->phi();
	fMdc.T.x[i]  = mdcTrk->x();
	fMdc.T.y[i]  = mdcTrk->y();
	fMdc.T.z[i]  = mdcTrk->z();
	fMdc.T.x[i]  = mdcTrk->x();
	fMdc.T.y[i]  = mdcTrk->y();
	fMdc.T.z[i]  = mdcTrk->z();
	double rvxy,rvz,rvphi;
	calculate_vertex(mdcTrk,rvxy,rvz,rvphi); 
	fMdc.T.vxy[i] = rvxy;
	fMdc.T.vz[i]  = rvz; 
	fMdc.T.vphi[i] = rvphi; 
}
