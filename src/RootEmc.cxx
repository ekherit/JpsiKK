// =====================================================================================
//
//       Filename:  RootEmc.cxx
//
//    Description:  
//
//        Version:  1.0
//        Created:  27.10.2015 15:32:04
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#include <algorithm>

#include "RootEmc.h"

int RootEmc::ARRAY_SIZE = 100;

void RootEmc::init_tuple(void)
{
  tuple->addItem ("ntrack",       ntrack, 0, ARRAY_SIZE); //good nuetral track in event
  tuple->addIndexedItem ("E",     ntrack, E);
  tuple->addIndexedItem ("theta", ntrack, theta);
  tuple->addIndexedItem ("phi",   ntrack, phi);
  tuple->addIndexedItem ("time",  ntrack, time);
}

void RootEmc::init(void)
{
  ntrack=4;
  for(int i=0;i<ntrack;i++)
  {
    E[i] = 0;
    theta[i] = -1000;;
    phi[i] = -1000;
    time[i] = -1000;
  }
}

void RootEmc::fill(int i,  EvtRecTrack * track)
{
	if(track->isEmcShowerValid())
	{
		RecEmcShower *emcTrk = track->emcShower();
		E[i] = emcTrk->energy();
		theta[i] = emcTrk->theta();
		phi[i] = emcTrk->phi();
		time[i] = emcTrk->time();
	}
}


void RootEmc::fill(list<EvtRecTrack*> & good_neutral_tracks)
{
	ntrack=std::min(good_neutral_tracks.size(), size_t(ARRAY_SIZE));
	int idx=0;
	for(list<EvtRecTrack*>::iterator track=good_neutral_tracks.begin(); track!=good_neutral_tracks.end(); track++)
	{
		EvtRecTrack   *Trk = *track;
		RecEmcShower *emcTrk = Trk->emcShower();
		E[idx]  =  emcTrk->energy();
		theta[idx] =  emcTrk->theta();
		phi[idx] =  emcTrk->phi();
		time[idx] = emcTrk->time();
		idx++;
	}
}
