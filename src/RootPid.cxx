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

StatusCode RootPid::init_tuple(void)
{
  StatusCode status;
  status = tuple->addItem ("Mee", M[ID_ELECTRON]);
  status = tuple->addItem ("MKK", M[ID_KAON]);
  status = tuple->addItem ("Muu", M[ID_MUON]);
  status = tuple->addItem ("Mpp", M[ID_PROTON]);
  status = tuple->addItem ("Mpipi", M[ID_PION]);

  status = tuple->addItem ("kMee", kM[ID_ELECTRON]);
  status = tuple->addItem ("kMKK", kM[ID_KAON]);
  status = tuple->addItem ("kMuu", kM[ID_MUON]);
  status = tuple->addItem ("kMpp", kM[ID_PROTON]);
  status = tuple->addItem ("kMpipi", kM[ID_PION]);
  //info prof ParticleID package
  status = tuple->addItem ("ntrack", ntrack,0,4); //array size must be = 4
  status = tuple->addIndexedItem ("probe",  ntrack, prob[ID_ELECTRON]);
  status = tuple->addIndexedItem ("probmu",  ntrack, prob[ID_MUON]);
  status = tuple->addIndexedItem ("probpi",  ntrack, prob[ID_PION]);
  status = tuple->addIndexedItem ("probk",  ntrack, prob[ID_KAON]);
  status = tuple->addIndexedItem ("probp",  ntrack, prob[ID_PROTON]);
  //my particle id information
  status = tuple->addIndexedItem ("chi2e",  ntrack, chi2[ID_ELECTRON]);
  status = tuple->addIndexedItem ("chi2mu",  ntrack, chi2[ID_MUON]);
  status = tuple->addIndexedItem ("chi2pi",  ntrack, chi2[ID_PION]);
  status = tuple->addIndexedItem ("chi2k",  ntrack, chi2[ID_KAON]);
  status = tuple->addIndexedItem ("chi2p",  ntrack, chi2[ID_PROTON]);
  return status;
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

