/*
 * =====================================================================================
 *
 *       Filename:  JPsi.h
 *
 *    Description:  Multihadron event selection for j/psi and psi prime resonance.
 *
 *        Version:  1.0
 *        Created:  04/27/2010 02:47:32 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan Nikolaev (), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physics
 *
 * =====================================================================================
 */

#ifndef IBN_TAUEMU_H
#define IBN_TAUEMU_H
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"

#include <TMatrixD.h>
#include <vector>
#include <algorithm>

#include "EventModel/EventHeader.h"

#include "ibn/averager.h"

//#include "EvtRecEvent/EvtRecTrack.h"
class JpsiKK : public Algorithm 
{
	public:
  JpsiKK(const std::string& name, ISvcLocator* pSvcLocator);
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();  

	private:
	int MIN_CHARGED_TRACKS; //minimum good charged tracks in selection
	int MAX_CHARGED_TRACKS; //maximum good charged tracks in selection
	int MAX_NEUTRAL_TRACKS; //maximum good neutral tracks in selection
	double IP_MAX_RHO, IP_MAX_Z; //interection point cut
  double MAX_COS_THETA; //maximum  cos(theta) for good charged track

  double EMC_ENDCUP_MIN_COS_THETA;
  double EMC_ENDCUP_MAX_COS_THETA;
  double EMC_ENDCUP_MIN_ENERGY;
  double EMC_BARREL_MAX_COS_THETA;
  double EMC_BARREL_MIN_ENERGY;

  double MAX_MUON_EP_RATIO;
  double MAX_KAON_EP_RATIO;

  double MAX_PION_MOMENTUM; //maximum pion momentum

  double MIN_RECOIL_MASS; //minimum recoil mass cut
  double MAX_RECOIL_MASS; //minimum recoil mass cut

  double MIN_KAON_MOMENTUM; //minimum kaon momentum
  double MAX_KAON_MOMENTUM; //maximum pion momentum

  double MIN_MUON_MOMENTUM; //minimum kaon momentum
  double MAX_MUON_MOMENTUM; //maximum pion momentum

  double MIN_INVARIANT_MASS; //minimum invariant  mass cut
  double MAX_INVARIANT_MASS; //manimum invariant  mass cut

  double MIN_KAON_MISSING_MASS;   //minimum kaon missing mass
  double MAX_KAON_MISSING_MASS;   //minimum kaon missing mass
  double MIN_MUON_MISSING_MASS;   //minimum muon missing mass
  double MAX_MUON_MISSING_MASS;   //minimum muon missing mass

  double MIN_MISSING_MASS; 
  double MAX_MISSING_MASS; 

	long int event_proceed;
	long int event_write;
  long int event_with_kaons;
  long int event_with_muons;

  struct RootEvent
  {
    NTuple::Item<long>    run; //run number
    NTuple::Item<long>    event; //event number 
    NTuple::Item<long>    time; //time of the event
    NTuple::Item<long>    ngood_charged_track;     //number of good charged tracks in event
    NTuple::Item<long>    ngood_neutral_track;     //number of good neutral tracks in event
    NTuple::Item<long>    npositive_track; //number of positive charged tracks
    NTuple::Item<long>    nnegative_track; //number of positive charged tracks
    NTuple::Item<long>    npositive_pions; //number of positive pions
    NTuple::Item<long>    nnegative_pions; //number of negative pions
    NTuple::Item<long>    npion_pairs; //total number of found pion pairs
    NTuple::Item<long>    channel;     //J/psi decay channel 0 -- K+K-, 1 -- mu+mu-
    NTuple::Item<double>  Mrecoil;  //pion recoil mass
    NTuple::Item<double>  M2missing; //missing square invariant mass
    NTuple::Item<double>  Minv; //invariant mass two charged particles

    NTuple::Item<long>    npid;     //number of particle hypo 0-K,1-mu,2-e,3-pi,4-p
    NTuple::Array<double> M;     

    NTuple::Item<long>    ntrack;  //size of the array = 4: [pi-,pi+,K-,K+]
    NTuple::Array<long> index; //index of track
    NTuple::Array<double> q; //charge of the track
    NTuple::Array<double> E;
    NTuple::Array<double> p;
    NTuple::Array<double> px;
    NTuple::Array<double> py;
    NTuple::Array<double> pz;
    NTuple::Array<double> pt; //transvese momentum
    NTuple::Array<double> theta,phi;
    NTuple::Array<double> x, y, z, r; //poca coordinate of track
    NTuple::Array<double> vxy, vz, vphi; //poca coordinate of track

    NTuple::Tuple * tuple; //tuple
    void init(void);
    StatusCode init_tuple(void);
  };

  struct RootNeutralTrack
  {
    NTuple::Item<long> ntrack; //number of good neutral track
    NTuple::Array<double> E;
    NTuple::Array<double> theta;
    NTuple::Array<double> phi;
    NTuple::Tuple * tuple; //tuple
    void init(void);
    StatusCode init_tuple(void);
  };

  RootEvent fEvent;
  RootNeutralTrack fNeutral;
};

#endif
