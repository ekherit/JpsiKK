/*
 * =====================================================================================
 *
 *       Filename:  JpsiKK.h
 *
 *    Description:  Event selection for measurement branching fraction Jpsi->K+K-
 *    via Psi(2S)->Jpsi(->K+K-)pi+pi- decay *
 *    
 *        Version:  1.0
 *        Created:  2014-11-06
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physics
 *
 * =====================================================================================
 */

#ifndef IBN_BES3_JPSIKK_H
#define IBN_BES3_JPSIKK_H
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"

#include "EventModel/EventHeader.h"

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"

#include "SelectionConfig.h"
#include "RootAdaptors.h"
#include "Defs.h"

class JpsiKK : public Algorithm 
{
	public:
  JpsiKK(const std::string& name, ISvcLocator* pSvcLocator);
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();  

	SelectionConfig cfg;

  private:

	//some online counters
	long int event_proceed;
	long int event_write;
  //long int event_with_kaons_and_muons;
  //long int good_kinematic_fit;
  //long int event_with_pions;
  //long int event_with_protons;
  //long int event_with_electrons;

  long int event_with_kaons;  //count 4C kinematic fit for kaons
  long int event_with_muons;  //count 4C kinematic fit for muons

  long int kaons_event_number;
  long int muons_event_number;
  std::map< int, std::map< int
  struct 
  {
    long int number
  }
  long int positive_kaon_event_number; //count only positive kinematic fit
  long int negative_kaon_event_number; //count only negative kinematic fit
  long int positive_muon_event_number;
  long int negative_muon_event_number;

	//long int nprints;

  RootEvent  fEvent;   //signal event essential information
  //RootPid    fPid;     //Paritlce id information
  RootMdc    fMdc;     //Mdc information
  RootDedx   fDedx;    //DeDx for the event
  RootEmc    fEmc;    //Emc infromation for the event
  RootTof    fTof;    //TOF infromation for the event
  RootMC     fMC;     //Monte Carlo truth
  RootMCTopo fMCTopo; //Monte Carlo topology
  RootEmc    fNeutral; //neutral tracks

	template <class A> 
		void init_tuple(A & a,  const char * dir, const char * title);


	void fillTuples(const std::vector<CLHEP::HepLorentzVector> & Pkf,  TrackVector_t & Tracks);
	void writeTuples(void);


	std::list<EvtRecTrackIterator> good_neutral_tracks;
	std::list<EvtRecTrackIterator> good_charged_tracks;
	EvtRecTrackIterator tracks_end; //a flag of the end of the track
};

#endif
