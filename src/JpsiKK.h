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
  //counter of the selected events theCounter[channel][charge][number of tracks]
  std::map< int, std::map< int, std::map < int, long int> > > theCounter;

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
  RootMuc    fMuc; //Muon chamber

	template <class A> 
		void init_tuple(A & a,  const char * dir, const char * title);


	void fillTuples(const std::vector<CLHEP::HepLorentzVector> & Pkf,  TrackVector_t & Tracks);
	void writeTuples(void);


	std::list<EvtRecTrack*> good_neutral_tracks;
	std::list<EvtRecTrack*> good_charged_tracks;

  void printSelectionDigest(bool head=false);
};

#endif
