// =====================================================================================
//
//       Filename:  RootTof.cxx
//
//    Description:  
//
//        Version:  1.0
//        Created:  27.10.2015 16:50:38
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#include "RootTof.h"
#include "DstEvent/TofHitStatus.h"

void RootTof::init_tuple(void)
{
  tuple->addItem ("ntrack", ntrack,0,4); 
  tuple->addIndexedItem ("ID",  ntrack, tofID);
  tuple->addIndexedItem ("t",  ntrack, t);
  tuple->addIndexedItem ("dt",  ntrack, dt);
  tuple->addIndexedItem ("t0",  ntrack, t0);
  tuple->addIndexedItem ("chie",  ntrack, chie);
  tuple->addIndexedItem ("chimu",  ntrack, chimu);
  tuple->addIndexedItem ("chipi",  ntrack, chipi);
  tuple->addIndexedItem ("chik",  ntrack, chik);
  tuple->addIndexedItem ("chip",  ntrack, chip);
  tuple->addIndexedItem ("beta",  ntrack, beta);
  tuple->addIndexedItem ("te",  ntrack, te);
  tuple->addIndexedItem ("tmu",  ntrack, tmu);
  tuple->addIndexedItem ("tpi",  ntrack, tpi);
  tuple->addIndexedItem ("tk",  ntrack, tk);
  tuple->addIndexedItem ("tp",  ntrack, tp);
}

void RootTof::init(void)
{
  ntrack=4;
  for(int i=0;i<ntrack;i++)
  {
    tofID[i]=-1000;
    t[i]=-1000;
    dt[i]=-1000;
    t0[i]=-1000;
    chie[i]=-1000;
    chimu[i]=-1000;
    chipi[i]=-1000;
    chik[i]=-1000;
    chip[i]=-1000;
    beta[i]=-1000;
    te[i]=-1000;
    tmu[i]=-1000;
    tpi[i]=-1000;
    tk[i]=-1000;
    tp[i]=-1000;
  }
}


void RootTof::fill(int i,  EvtRecTrack * track)
{
	if(track->isTofTrackValid())
	{
		SmartRefVector<RecTofTrack> tofTrkCol = track->tofTrack();
		SmartRefVector<RecTofTrack>::iterator tofTrk = tofTrkCol.begin();
		TofHitStatus *hitst = new TofHitStatus;
		std::vector<int> tofecount;
		int goodtofetrk=0;
		for(tofTrk = tofTrkCol.begin(); tofTrk!=tofTrkCol.end(); tofTrk++,goodtofetrk++)
		{
			unsigned int st = (*tofTrk)->status();
			hitst->setStatus(st);
			//if(  (hitst->is_barrel()) ) continue;
			if( !(hitst->is_counter()) ) continue;
			tofecount.push_back(goodtofetrk);
		}
		delete hitst;
		if(!tofecount.empty())
		{
			tofTrk = tofTrkCol.begin()+tofecount[0];
			tofID[i] = (*tofTrk)->tofID();
			t0[i] = (*tofTrk)->t0();
			t[i] = (*tofTrk)->tof();
			dt[i] = (*tofTrk)->errtof();
			beta[i] = (*tofTrk)->beta();
			te[i] = (*tofTrk)->texpElectron();
			tmu[i]= (*tofTrk)->texpMuon();
			tpi[i]= (*tofTrk)->texpPion();
			tk[i] = (*tofTrk)->texpKaon();
			tp[i] = (*tofTrk)->texpProton();
			if(dt[i]>0)
			{
				chie[i]  = (t[i] - te[i])  /  fabs(dt[i]);
				chimu[i] = (t[i] - tmu[i]) /  fabs(dt[i]);
				chipi[i] = (t[i] - tpi[i]) /  fabs(dt[i]);
				chik[i]  = (t[i] - tk[i])  /  fabs(dt[i]);
				chip[i]  = (t[i] - tp[i])  /  fabs(dt[i]);
			}
		}
	}

}
