// =====================================================================================
//
//       Filename:  RootDedx.cxx
//
//    Description:  
//
//        Version:  1.0
//        Created:  27.10.2015 16:47:35
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#include "RootDedx.h"

void RootDedx::init_tuple(void)
{
  tuple->addItem ("ntrack", ntrack,0,4); 
  tuple->addIndexedItem ("chie",  ntrack, chie);
  tuple->addIndexedItem ("chimu",  ntrack, chimu);
  tuple->addIndexedItem ("chipi",  ntrack, chipi);
  tuple->addIndexedItem ("chik",  ntrack, chik);
  tuple->addIndexedItem ("chip",  ntrack, chip);
  tuple->addIndexedItem ("probPH",  ntrack, probPH);
  tuple->addIndexedItem ("normPH",  ntrack, normPH);
  //tuple->addIndexedItem ("probe",  ntrack, probe);
  //tuple->addIndexedItem ("probmu",  ntrack, probmu);
  //tuple->addIndexedItem ("probpi",  ntrack, probpi);
  //tuple->addIndexedItem ("probk",  ntrack, probk);
  //tuple->addIndexedItem ("probp",  ntrack, probp);
}

void RootDedx::init(void)
{
  ntrack=4;
  for(int i=0;i<ntrack;i++)
  {
    chie[i]=-1000;
    chimu[i]=-1000;
    chipi[i]=-1000;
    chik[i]=-1000;
    chip[i]=-1000;
    probPH[i]=0;
    normPH[i]=0;
  }
}

void RootDedx::fill(int i, EvtRecTrack * track)
{
	if(track->isMdcDedxValid())
	{
		RecMdcDedx* dedxTrk = track->mdcDedx();
		chie[i] = dedxTrk->chiE();
		chimu[i] = dedxTrk->chiMu();
		chipi[i] = dedxTrk->chiPi();
		chik[i] = dedxTrk->chiK();
		chip[i] = dedxTrk->chiP();
		//fDedx.ghit[i] = dedxTrk->numGoodHits();
		//fDedx.thit[i] = dedxTrk->numTotalHits();
		probPH[i] = dedxTrk->probPH();
		normPH[i] = dedxTrk->normPH();
		//fDedx.e[i] = dedxTrk->getDedxExpect(0);
		//fDedx.mu[i] = dedxTrk->getDedxExpect(1);
		//fDedx.pi[i] = dedxTrk->getDedxExpect(2);
		//fDedx.K[i] = dedxTrk->getDedxExpect(3);
		//fDedx.p[i] = dedxTrk->getDedxExpect(4);
		//fDedx.pid[i]=dedxTrk->particleId();
	}
}
