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

StatusCode RootEvent::init_tuple(void)
{
  StatusCode status;
  status = tuple->addItem ("run", run); //run number
  status = tuple->addItem ("event", event); //event number
  status = tuple->addItem ("time", time); //event time
  status = tuple->addItem ("ngtrack", ngood_charged_track); //good charged track in event
  status = tuple->addItem ("ngntrack", ngood_neutral_track); //good neutral track in event
  status = tuple->addItem ("nptrack", npositive_track); //good positive charged track in event
  status = tuple->addItem ("nntrack", nnegative_track); //good negative charged track in event
  status = tuple->addItem ("nppions", npositive_pions); //good poitive pion tracks in event
  status = tuple->addItem ("nnpions", nnegative_pions); //good negative pion track in event
  status = tuple->addItem ("npion_pairs", npion_pairs); //number of pions paris in event
  //status = tuple->addItem ("ntrack", T.ntrack,0,4);     

  status = tuple->addItem ("sign", sign); //number of pions paris in event
  status = tuple->addItem ("channel", channel); //decay channel of the J/psi
  status = tuple->addItem ("KK", KK); //KK decay channel of the J/psi
  status = tuple->addItem ("uu", uu); //mu-mu decay channel of the J/psi
  status = tuple->addItem ("Ku", Ku); //Kmu or muK events

	status = M.add_to_tuple(tuple);
  //status = tuple->addItem ("Mee",   kM[ID_ELECTRON]);
  //status = tuple->addItem ("MKK",   kM[ID_KAON]);
  //status = tuple->addItem ("Muu",   kM[ID_MUON]);
  //status = tuple->addItem ("Mpp",   kM[ID_PROTON]);
  //status = tuple->addItem ("Mpipi", kM[ID_PION]);

	status = T.add_to_tuple(tuple);

  status = tuple->addItem ("kin_chi2", kin_chi2); 
  status = tuple->addItem ("pid_chi2", pid_chi2); 

  status = tuple->addItem ("npid", npid,0,5);     
  status = tuple->addIndexedItem ("kchi",  npid, kchi);
  status = tuple->addIndexedItem ("pchi",  npid, pchi);
  status = tuple->addIndexedItem ("kM",  npid,     kM23);

  status = tuple->addIndexedItem ("probe",  T.ntrack,  prob[ID_ELECTRON]);
  status = tuple->addIndexedItem ("probmu", T.ntrack,  prob[ID_MUON]);
  status = tuple->addIndexedItem ("probpi", T.ntrack,  prob[ID_PION]);
  status = tuple->addIndexedItem ("probk",  T.ntrack,  prob[ID_KAON]);
  status = tuple->addIndexedItem ("probp",  T.ntrack,  prob[ID_PROTON]);


  return status;
}

void RootEvent::init(void)
{
  T.ntrack=4;
	npid = 5; //number of particle id hypotesis
}
