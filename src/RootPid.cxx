// =====================================================================================
//
//       Filename:  RootPid.cxx
//
//    Description:  
//
//        Version:  1.0
//        Created:  27.10.2015 17:01:58
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#include "RootPid.h"
#include "PhysConst.h"

void RootPid::init_tuple(void)
{
  tuple->addItem ("Mee", M[ID_ELECTRON]);
  tuple->addItem ("MKK", M[ID_KAON]);
  tuple->addItem ("Muu", M[ID_MUON]);
  tuple->addItem ("Mpp", M[ID_PROTON]);
  tuple->addItem ("Mpipi", M[ID_PION]);

  tuple->addItem ("kMee", kM[ID_ELECTRON]);
  tuple->addItem ("kMKK", kM[ID_KAON]);
  tuple->addItem ("kMuu", kM[ID_MUON]);
  tuple->addItem ("kMpp", kM[ID_PROTON]);
  tuple->addItem ("kMpipi", kM[ID_PION]);
  //info prof ParticleID package
  tuple->addItem ("ntrack", ntrack,0,4); //array size must be = 4
  tuple->addIndexedItem ("probe",  ntrack, prob[ID_ELECTRON]);
  tuple->addIndexedItem ("probmu",  ntrack, prob[ID_MUON]);
  tuple->addIndexedItem ("probpi",  ntrack, prob[ID_PION]);
  tuple->addIndexedItem ("probk",  ntrack, prob[ID_KAON]);
  tuple->addIndexedItem ("probp",  ntrack, prob[ID_PROTON]);
  //my particle id information
  tuple->addIndexedItem ("chi2e",  ntrack, chi2[ID_ELECTRON]);
  tuple->addIndexedItem ("chi2mu",  ntrack, chi2[ID_MUON]);
  tuple->addIndexedItem ("chi2pi",  ntrack, chi2[ID_PION]);
  tuple->addIndexedItem ("chi2k",  ntrack, chi2[ID_KAON]);
  tuple->addIndexedItem ("chi2p",  ntrack, chi2[ID_PROTON]);
}


void RootPid::init(void)
{
  ntrack=4;
}

void RootPid::fill(int i,  EvtRecTrackIterator & track)
{
	//PID->setRecTrack((*Tracks[i]));
	//PID->calculate();
	//if(PID->IsPidInfoValid())
	//{
	//  //fPid.prob[ID_ELECTRON][i] = PID->probElectron();
	//  //fPid.prob[ID_MUON][i]     = PID->probMuon();
	//  //fPid.prob[ID_PION][i]     = PID->probPion();
	//  //fPid.prob[ID_KAON][i]     = PID->probKaon();
	//  //fPid.prob[ID_PROTON][i]   = PID->probProton();
	//}
	//vector<double> chi2 = get_chi2(Tracks[i]);
	//for(int pid=0;pid<5;pid++)
	//{
	//  //fPid.chi2[pid][i]   = chi2[pid];
	//}
  //for(int i=0;i<5;i++)
  //{
  //  //fPid.M[i]    = sqrt(get_invariant_mass2(result_pair,XMASS[i]));
  //  HepLorentzVector p1(Pkf[2].vect(), XMASS[i]);
  //  HepLorentzVector p2(Pkf[3].vect(), XMASS[i]);
  //  //fPid.kM[i] = (p1+p2).m();
  //};
}

