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
	//M.add_to_tuple(tuple);
	T.add_to_tuple(tuple);
	tuple->addItem ("Mrec", Mrec); //array size must be = 4
	tuple->addItem ("npid", npid,0,5); 
	tuple->addIndexedItem ("M23", npid, M23);
	tuple->addIndexedItem ("M12", npid, M12);
	tuple->addIndexedItem ("M03", npid, M03);
	tuple->addIndexedItem ("Mmis", npid, Mmis);
}


void RootMdc::init(void)
{
  T.ntrack=4;
	npid =5;
}


void RootMdc::fill(int i, EvtRecTrackIterator & track)
{
	if(!(*track)->isMdcTrackValid()) return; 
	//RecMdcTrack  *mdcTrk = (*track)->mdcTrack();
  RecMdcKalTrack * mdcTrk = (*track)->mdcKalTrack();
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
	calculate_vertex((*track)->mdcTrack(),rvxy,rvz,rvphi); 
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

void RootMdc::fill_mass(TrackVector_t & tracks,   EvtRecTrackIterator & end, const std::vector<HepLorentzVector> & P)
{
	npid=5;

	for(int pid =0; pid <5;pid++)
	{
		if(tracks[2] != end && tracks[3] != end)
		{
			std::vector<int> pids(4);
			pids[0] = ID_PION;
			pids[1] = ID_PION;
			pids[2] = pid;
			pids[3] = pid;
			M23[pid] = sqrt(getInvariantMass2(pids[2], tracks[2], pids[3], tracks[3]));
			M12[pid] = sqrt(getInvariantMass2(pids[1], tracks[1], pids[2], tracks[2]));
			M03[pid] = sqrt(getInvariantMass2(pids[0], tracks[0], pids[3], tracks[3]));
			Mmis[pid] = sqrt(getMissingMass2(tracks, pids));
		}
		else
		{
			std::vector<int> pids(3);
			pids[0] = ID_PION;
			pids[1] = ID_PION;
			pids[2] = pid;
			TrackVector_t T(3);
			T[0] = tracks[0];
			T[1] = tracks[1];
			if(tracks[2] != end)
			{
				T[2] = tracks[2];
				M12[pid]  = sqrt(getInvariantMass2(pids[1], tracks[1], pid, tracks[2]));
				M03[pid]  = 0;
        M23[pid]  = sqrt(getInvariantMass2(pid,tracks[2], P[3]));
			}
			if(tracks[3] != end)
			{
				T[2] = tracks[3];
				M12[pid]  = 0;
				M03[pid]  = sqrt(getInvariantMass2(pids[0], tracks[0], pid, tracks[3]));
        M23[pid]  = sqrt(getInvariantMass2(pid,tracks[3], P[2]));
			}
			//M23[pid]  = 0;
			Mmis[pid]  = sqrt(getMissingMass2(T, pids));
		}
	}
	Mrec = getPionRecoilMass(tracks[0],  tracks[1]);
}
