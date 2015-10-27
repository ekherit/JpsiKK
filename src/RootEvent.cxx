// =====================================================================================
//
//       Filename:  RootEvent.cxx
//
//    Description:  
//
//        Version:  1.0
//        Created:  27.10.2015 16:59:57
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#include "RootEvent.h"
#include "PhysConst.h"
#include "utils.h"

void RootEvent::init_tuple(void)
{
  tuple->addItem ("run", run); //run number
  tuple->addItem ("event", event); //event number
  tuple->addItem ("time", time); //event time
  tuple->addItem ("ngtrack", ngood_charged_track); //good charged track in event
  tuple->addItem ("ngntrack", ngood_neutral_track); //good neutral track in event
  tuple->addItem ("nptrack", npositive_track); //good positive charged track in event
  tuple->addItem ("nntrack", nnegative_track); //good negative charged track in event
  tuple->addItem ("nppions", npositive_pions); //good poitive pion tracks in event
  tuple->addItem ("nnpions", nnegative_pions); //good negative pion track in event
  tuple->addItem ("npion_pairs", npion_pairs); //number of pions paris in event
  tuple->addItem ("ntrack", T.ntrack,0,4);     

  tuple->addItem ("sign", sign); //number of pions paris in event
  tuple->addItem ("channel", channel); //decay channel of the J/psi
  tuple->addItem ("KK", KK); //KK decay channel of the J/psi
  tuple->addItem ("uu", uu); //mu-mu decay channel of the J/psi
  tuple->addItem ("Ku", Ku); //Kmu or muK events

	M.add_to_tuple(tuple);
  //tuple->addItem ("Mee",   kM[ID_ELECTRON]);
  //tuple->addItem ("MKK",   kM[ID_KAON]);
  //tuple->addItem ("Muu",   kM[ID_MUON]);
  //tuple->addItem ("Mpp",   kM[ID_PROTON]);
  //tuple->addItem ("Mpipi", kM[ID_PION]);

	T.add_to_tuple(tuple);

  tuple->addItem ("kin_chi2", kin_chi2); 
  tuple->addItem ("pid_chi2", pid_chi2); 

  tuple->addItem ("npid", npid,0,5);     
  tuple->addIndexedItem ("kchi",  npid, kchi);
  tuple->addIndexedItem ("pchi",  npid, pchi);
  tuple->addIndexedItem ("kM",  npid,     kM23);

  //tuple->addIndexedItem ("probe",  T.ntrack,  prob[ID_ELECTRON]);
  //tuple->addIndexedItem ("probmu", T.ntrack,  prob[ID_MUON]);
  //tuple->addIndexedItem ("probpi", T.ntrack,  prob[ID_PION]);
  //tuple->addIndexedItem ("probk",  T.ntrack,  prob[ID_KAON]);
  //tuple->addIndexedItem ("probp",  T.ntrack,  prob[ID_PROTON]);
}

void RootEvent::init(void)
{
  T.ntrack=4;
	npid = 5; //number of particle id hypotesis
}

void RootEvent::fill(const std::vector<HepLorentzVector> & Pkf,  double CENTER_MASS_ENERGY)
{
  HepLorentzVector Pcm(CENTER_MASS_ENERGY*sin(0.011),0,0,CENTER_MASS_ENERGY); 

  M.Mrec = (Pcm - Pkf[0] - Pkf[1]).m();

  M.M012 = (Pkf[0]+Pkf[1]+Pkf[2]).m();
  M.M013 = (Pkf[0]+Pkf[1]+Pkf[3]).m();
  M.M023 = (Pkf[0]+Pkf[2]+Pkf[3]).m();
  M.M123 = (Pkf[1]+Pkf[2]+Pkf[3]).m();

  M.M03 =  (Pkf[0]+Pkf[3]).m();
  M.M12 =  (Pkf[1]+Pkf[2]).m();
  M.M01 =  (Pkf[0]+Pkf[1]).m();
  M.M23 =  (Pkf[2]+Pkf[3]).m();

	T.ntrack = Pkf.size();
	for ( int i=0;i<4;i++)
	{
		T.q[i]  = i%2 == 0 ? -1 : +1;
		T.E[i]  = Pkf[i].e();
		T.px[i] = Pkf[i].px();
		T.py[i] = Pkf[i].py();
		T.pz[i] = Pkf[i].pz();
		T.p[i]  = sqrt(sq(Pkf[i].px())+sq(Pkf[i].py())+sq(Pkf[i].pz()));
		T.pt[i] = sqrt(sq(Pkf[i].px())+sq(Pkf[i].py()));
		T.theta[i]= Pkf[i].theta();
		T.phi[i] = Pkf[i].phi();
		T.x[i]=0;
		T.y[i]=0;
		T.z[i]=0;
		T.r[i]=0;
		T.vxy[i]=0;
		T.vz[i]=0;
		T.vphi[i]=0;
	}
}
