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

class JpsiKK : public Algorithm 
{
	public:
  JpsiKK(const std::string& name, ISvcLocator* pSvcLocator);
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();  

	SelectionConfig cfg;

	private:

	long int event_proceed;
	long int event_write;
  long int event_with_kaons;
  long int event_with_muons;
  long int good_kinematic_fit;
  long int event_with_pions;
  long int event_with_protons;
  long int event_with_electrons;
	protected:

  public:
  struct RootTuple
  {
    public:
    NTuple::Tuple * tuple; //tuple
    virtual ~RootTuple(void){};
    virtual void init(void)=0;
    virtual StatusCode init_tuple(void)=0;
  };


	struct Mass_t
	{
    NTuple::Item<double>  Mrec;  //pion recoil mass
    NTuple::Item<double>  M012; //M(pi pi K/mu-)
    NTuple::Item<double>  M013; //M(pi pi K/mu+)
    NTuple::Item<double>  M023; //M(pi- KK/uu)
    NTuple::Item<double>  M123; //M(pi+ KK/uu)

    NTuple::Item<double>  M03; //invariant mass of Kaon and pion
    NTuple::Item<double>  M12; //invariant mass of kaon and pion
    NTuple::Item<double>  M01; //invariant mass of pion
    NTuple::Item<double>  M23; //invariant mass of kaons or muons

    NTuple::Array<double> Mmis;    //missing invariant mass

    virtual StatusCode add_to_tuple(NTuple::Tuple * tuple)
		{
			StatusCode status;
			status = tuple->addItem ("Mrec", Mrec); 
			status = tuple->addItem ("M012", M012); 
			status = tuple->addItem ("M013", M013); 
			status = tuple->addItem ("M023", M023); 
			status = tuple->addItem ("M123", M123); 

			status = tuple->addItem ("M03", M03); 
			status = tuple->addItem ("M12", M12); 
			status = tuple->addItem ("M01", M01); 
			status = tuple->addItem ("M23", M23); 
			return status;
		}
	};

	struct Track_t
	{
    NTuple::Item<long>    ntrack;  //size of the array = 4: [pi-,pi+,K-,K+]
    NTuple::Array<long>   trackId; //id of the track
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
    virtual StatusCode add_to_tuple(NTuple::Tuple * tuple)
		{
			StatusCode status;
			status = tuple->addItem ("ntrack", ntrack,0,4); //array size must be = 4
			status = tuple->addIndexedItem ("trackId",   ntrack, trackId);
			status = tuple->addIndexedItem ("q",     ntrack, q);
			status = tuple->addIndexedItem ("E",     ntrack, E);
			status = tuple->addIndexedItem ("p",     ntrack, p);
			status = tuple->addIndexedItem ("px",    ntrack, px);
			status = tuple->addIndexedItem ("py",    ntrack, py);
			status = tuple->addIndexedItem ("pz",    ntrack, pz);
			status = tuple->addIndexedItem ("pt",    ntrack, pt);
			status = tuple->addIndexedItem ("theta", ntrack, theta);
			status = tuple->addIndexedItem ("phi",   ntrack, phi);
			status = tuple->addIndexedItem ("x",     ntrack, x);
			status = tuple->addIndexedItem ("y",     ntrack, y);
			status = tuple->addIndexedItem ("z",     ntrack, z);
			status = tuple->addIndexedItem ("r",     ntrack, r);
			status = tuple->addIndexedItem ("vxy",   ntrack, vxy);
			status = tuple->addIndexedItem ("vz",    ntrack, vz);
			status = tuple->addIndexedItem ("vphi",  ntrack, vphi);
			return status;
		}
	};


	// =====================================================================================
	//        Class:  RootEvent
	//  Description:  Main event information supposed to used for selection
	// =====================================================================================
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

		NTuple::Item<long>    sign;   //signature of the event shows missed particle  01(K- or mu-), 10 (K+ or mu+),  11 (KK or uu or Ku or uK)
    NTuple::Item<long>    channel;     //J/psi decay channel 0 -- kaons, 1 -- muons,  10 -- muK,  11 - Kmu
    NTuple::Item<long>    KK;     //KK Jpsi decay event
    NTuple::Item<long>    uu;     //MuMu event
    NTuple::Item<long>    Ku;     //K mu or mu K events 


		Mass_t  M;
    NTuple::Item<double> kM[5];    //Invariant mass for track 23 for different hyptotesis

		//here will be result of the kinematic fit and particle id
    NTuple::Item<double>  kin_chi2; //kinematic chi2
    NTuple::Item<double>  pid_chi2; //my prob  chi2

		Track_t T;

    NTuple::Array<double> kchi;    //kinematik chi2
    NTuple::Array<double> pchi;    //my particle id chi square

    NTuple::Array<double> prob[5]; //probability from ParticleID
    NTuple::Array<double> pchi2[5];//my particle id chi square different hypo

    virtual void init(void);
    virtual StatusCode init_tuple(void);
  };


	// =====================================================================================
	//        Class:  RootMdc
	//  Description:  Discribes track information from Mdc before Vertex and 
	//  Kinematik fit.
	// =====================================================================================
  struct RootMdc : public RootTuple
  {
		Mass_t M;
		Track_t T;
    virtual void init(void);
    virtual StatusCode init_tuple(void);
  };


  struct RootPid : public RootTuple
  {
    NTuple::Item<double>  M[5]; //invariant mass of highmomentum track  based on Mdc
    NTuple::Item<double>  kM[5]; //invariant mass of highmomentum track  based on kinematic recons
    NTuple::Item<long>    ntrack;  //size of the array = 4: [pi-,pi+,K-,K+]
    NTuple::Array<double> prob[5]; //probability of track to be e,mu,pi,k or p
    NTuple::Array<double> chi2[5];  
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
    static int ARRAY_SIZE; 
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
    NTuple::Array<double> t; //tof time
    NTuple::Array<double> dt; //error of tof time
    NTuple::Array<double> beta;  
    NTuple::Array<double> te;  //electron expected time
    NTuple::Array<double> tmu; //muon
    NTuple::Array<double> tpi; //pion
    NTuple::Array<double> tk;  //kaon
    NTuple::Array<double> tp;  //proton

    NTuple::Array<double> chie;  
    NTuple::Array<double> chimu; 
    NTuple::Array<double> chipi; 
    NTuple::Array<double> chik;  
    NTuple::Array<double> chip;  
    virtual void init(void);
    virtual StatusCode init_tuple(void);
  };


  struct RootMC : public RootTuple
  {
    NTuple::Item<long>    psip_decay; //is event from psip decay
    NTuple::Item<long>    jpsi_decay; //is event from jpsi decay
    NTuple::Item<long>    KK; //is event KK
    NTuple::Item<long>    uu; //is event mumu
    NTuple::Item<long>    oo; //other other particle
    NTuple::Item<long>    ntrack;  //size of the array = 4: [pi-,pi+,K-,K+]
    NTuple::Array<double> pid; //particle id
    NTuple::Array<double> q; //charge of the track
    NTuple::Array<double> E,p; //energy and momentum
    NTuple::Array<double> px,py,pz; //momentum
    NTuple::Array<double> pt; //transvese momentum
    NTuple::Array<double> theta,phi;
    virtual void init(void);
    virtual StatusCode init_tuple(void);
  };


  struct RootMCTopo : public RootTuple
  {
    NTuple::Item  <int>  m_idxmc;
    NTuple::Array <int>  m_pdgid;
    NTuple::Array <int>  m_motheridx;
    NTuple::Array <int>  m_idx;
    NTuple::Item <unsigned long>  m_hash;
    virtual void init(void);
    virtual StatusCode init_tuple(void);
  };



  private:

  RootEvent  fEvent;   //signal event essential information
  //RootPid    fPid;     //Paritlce id information
  RootMdc    fMdc;     //Mdc information
  RootDedx   fDedx;    //DeDx for the event
  RootEmc    fEmc;    //Emc infromation for the event
  RootTof    fTof;    //TOF infromation for the event
  RootMC     fMC;     //Monte Carlo truth
  RootMCTopo fMCTopo; //Monte Carlo topology
  RootEmc    fNeutral; //neutral tracks
};

#endif
