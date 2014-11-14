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

  public:
  struct RootTuple
  {
    public:
    NTuple::Tuple * tuple; //tuple
    virtual ~RootTuple(void){};
    virtual void init(void)=0;
    virtual StatusCode init_tuple(void)=0;
  };

  struct RootEvent : public RootTuple
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
    NTuple::Array<double> prob;     

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
    NTuple::Array<double> probe,probmu,probpi,probk,probp; //probability of track to be e,mu,pi,k or p
    virtual void init(void);
    virtual StatusCode init_tuple(void);
  };

  struct RootDedx : public RootTuple
  {
    NTuple::Item<long>    ntrack;  //size of the array = 4: [pi-,pi+,K-,K+]
    NTuple::Array<double> chie;  //chi e
    NTuple::Array<double> chimu; //chi e
    NTuple::Array<double> chipi; //chi e
    NTuple::Array<double> chik;  //chi e
    NTuple::Array<double> chip;  //chi e
    NTuple::Array<double> probPH;  //хрень какая-то
    NTuple::Array<double> normPH;
    //NTuple::Array<double> probe;  //prob e
    //NTuple::Array<double> probmu; //prob e
    //NTuple::Array<double> probpi; //prob e
    //NTuple::Array<double> probk;  //prob e
    //NTuple::Array<double> probp;  //prob e
    virtual void init(void);
    virtual StatusCode init_tuple(void);
  };

  struct RootEmc : public RootTuple
  {
    NTuple::Item<long> ntrack;
    NTuple::Array<double> E;
    NTuple::Array<double> theta;
    NTuple::Array<double> phi;
    NTuple::Array<double> time;
    virtual void init(void);
    virtual StatusCode init_tuple(void);
  };

  struct RootTof : public RootTuple
  {
    NTuple::Item<long> ntrack;
    NTuple::Array<double> tofID;
    NTuple::Array<double> t0;
    NTuple::Array<double> tof;
    NTuple::Array<double> errtof;
    NTuple::Array<double> chie;  
    NTuple::Array<double> chimu; 
    NTuple::Array<double> chipi; 
    NTuple::Array<double> chik;  
    NTuple::Array<double> chip;  
    virtual void init(void);
    virtual StatusCode init_tuple(void);
  };


  private:

  RootEvent fEvent; //signal event essential information
  RootDedx  fDedx; //dedx for the event
  RootEmc   fEmc;  //emc infromation for the event
  RootEmc   fTof;  //emc infromation for the event
  RootEmc   fNeutral; //neutral tracks
};

#endif
