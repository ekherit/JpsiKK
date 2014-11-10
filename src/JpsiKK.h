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

class JpsiKK : public Algorithm 
{
	public:
  JpsiKK(const std::string& name, ISvcLocator* pSvcLocator);
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();  

	private:
	int MIN_CHARGED_TRACKS; //minimum charged tracks in selection
	int MAX_CHARGED_TRACKS; //minimum charged tracks in selection
	double IP_MAX_RHO, IP_MAX_Z; //interection point cut
  double MAX_COS_THETA; //maximum  cos(theta) for good charged track

  double MAX_PION_MOMENTUM; //maximum pion momentum
  double MIN_KAON_MOMENTUM; //minimum kaon momentum
  double MAX_KAON_MOMENTUM; //maximum pion momentum
  double MIN_MUON_MOMENTUM; //minimum kaon momentum
  double MAX_MUON_MOMENTUM; //maximum pion momentum

  double MIN_RECOIL_MASS; //minimum recoil mass cut
  double MAX_RECOIL_MASS; //minimum recoil mass cut
  double MIN_KAON_MISSING_MASS;   //minimum kaon missing mass
  double MAX_KAON_MISSING_MASS;   //minimum kaon missing mass
  double MIN_MUON_MISSING_MASS;   //minimum muon missing mass
  double MAX_MUON_MISSING_MASS;   //minimum muon missing mass


  int STRICT_TAU_CUT;
	int USE_IPCUT; //use interection point cut
	int CHECK_TOF; //use toff.
  int CHECK_DEDX;
  int CHECK_MUC;
  int CHECK_MC;
  double MC_DP; //momentum difference to identify particles
	double IPR;
	int IPTRACKS; //tracks number from interection point
	int MAX_TRACK_NUMBER; //minimum charged tracks in selection
	long int event_proceed;
	long int event_write;
  long int tau_events; //number of mu ivents for track #1
  long int bhabha_events; 
  long int gg_events; 

	NTuple::Tuple * main_tuple;//main tuple 
	NTuple::Item<long> m_time; //time when events is writed (unixtime)
	NTuple::Item<long> m_ntrack; //total number of tracks.
	NTuple::Item<long> m_nchtr; //number of charged tracks
	NTuple::Item<long> m_nneutr;//number of neutral tracks
  NTuple::Item<double> m_Etotal; //total energy deposition
  NTuple::Item<double> m_Eemc; //total energy deposition EMC only

  /*  Main Drift Chamber Information */
  struct MDC_t
  {
    NTuple::Item<long>    ntrack; //number of charged tracks.
    NTuple::Item<long>    ngood_track; //number of good tracks.
    NTuple::Array<double> p; //Momentum
    NTuple::Array<double> px,py,pz; //Componets of momentum
    NTuple::Array<double> pt; //transvese momentum
    NTuple::Array<double> x, y, z, r; //poca coordinate of track
    NTuple::Array<double> rvxy, rvz, rvphi; //poca coordinate of track
    NTuple::Array<double> theta,phi;
    NTuple::Array<double> q; //charge of the track
    NTuple::Array<long>   isemc; //has emc information
    /*  EMC section for this charged track */
    NTuple::Array<double> ncrstl;
    NTuple::Array<double> cellId;
    NTuple::Array<double> status;
    NTuple::Array<double> module;
    NTuple::Array<double> E,dE;
    NTuple::Array<double> temc; // Time of central crystal in the shower.
    NTuple::Array<double> M;
    NTuple::Array<double> ismu; //has muon track information
    NTuple::Array<double> istof; //has tof track information
    /*  Additional section */
    NTuple::Array<double> X, Y, Z; //pivot 
    NTuple::Item<long>   nip; //Track number from interaction point (same condition as for highest track)
    NTuple::Item<double> Mrec; //Recoil mass of the pions
    NTuple::Item<double> MKK; //Invariant mass of the Kaons
    NTuple::Item<double> Mmumu; //Invariant mass of the muons
    NTuple::Item<double> Mmiss; //Miss invariant mass
    NTuple::Item<double> Eemc; //total energy using emc
    NTuple::Item<double> Emdc; //total energy using only mdc
    NTuple::Item<double> S; //Sphericity
    NTuple::Item<double> ccos; //cos between two high energy tracks.
    NTuple::Item<double> atheta, aphi; //acolinearity 
    NTuple::Item<double> acompl; //acomplanarity
    /*  particle identificaiton */
    NTuple::Item<int> jpsi_decay_channel;  //0 - KK channel, 1 -- mumu channel
    NTuple::Array<double>   probe; //probability of electron
    NTuple::Array<double>   probmu; //probability for track to be muon
    NTuple::Array<double>   probK; //probability for track to be Kaon
    NTuple::Array<double>   probpi; //probability for track to be pion
    NTuple::Array<double>   probp; //probability for track to be proton
  };

  /* ElecroMagnetic Calorimeter Information */
  struct EMC_t
  {
    NTuple::Item<long>    ntrack; //number of neutral tracks.
    NTuple::Item<long>    ngood_track; //number of good neutral tracks.
    NTuple::Item<long>    ngood_charged_track; //number of good charged tracks.

    NTuple::Item<double> S; //sphericity
    NTuple::Item<double> ccos; //aclolinearity
    NTuple::Item<double> atheta; //theta acolinearity theta_0+theta_1 - pi
    NTuple::Item<double> aphi; //phi aclolinearity  abs(phi_0-phi_1)-pi
    //arrays
    NTuple::Array<double>   status; //status status=1: single seed cluster; status=2: splitted from multi-seeds cluster.
    NTuple::Array<double>   ncrstl; //Number of crystals in the shower
    NTuple::Array<double>   cellId; //Central crystal’s identifier
    NTuple::Array<double>   module; //module=0: east endcap;  module=1: barrel;  module=2: west endcap.
    NTuple::Array<double> x, y, z; //coordinates of claster
    NTuple::Array<double> theta, phi;  //angles
    NTuple::Array<double> E,dE; // energy deposition and error
    NTuple::Array<double> t; //Time of central crystal in the shower.
    NTuple::Item<double>  Etotal;
    void init(void);
    StatusCode init_tuple(NTuple::Tuple * tuple);
    int MAX_TRACK_NUMBER;
  };
	
  NTuple::Tuple * mdc_tuple;
  MDC_t mdc;

  NTuple::Tuple * emc_tuple;
  EMC_t emc;

  struct MUC_t
  {
    NTuple::Item<long>    ntrack; //number of selected tracks. Some of them dont have muc information
    NTuple::Array<double>   status; //status status=1: single seed cluster; status=2: splitted from multi-seeds cluster.
    NTuple::Array<double>   type; //seed mode. 0: ext, 1: emc, 2: muc
    NTuple::Array<double>   depth; //depth of the track in iron
    NTuple::Array<double>   chi2; //chi2 of the fit
    NTuple::Array<double>   ndf; //degree of freedom
    NTuple::Array<double>   distance; //distance between ext track and fired strip in 1st layer of MUC
    NTuple::Array<double>   phi; //delta phi between mdc momentum and direction of MUC track
    NTuple::Array<double>   nhit; //number of hits
    NTuple::Array<double>   nlayer; //number of hits
    NTuple::Array<double>   nhitmax; //max number of hits in layer
    NTuple::Array<double>   brlast; //last layer in barrel
    NTuple::Array<double>   eclast; //last layer in end cup
    void init(void);
    StatusCode init_tuple(NTuple::Tuple * tuple);
    int MAX_TRACK_NUMBER;
  };

  MUC_t muc;

  NTuple::Tuple * muc_tuple;






  struct DEDX_t
  {
	  NTuple::Item<long> ntrack;
	  NTuple::Array<double>   pid; //particle id (dedx)
	  NTuple::Array<double> chie; //chi2_dEdx for electron
	  NTuple::Array<double> chimu;//chi2_dEdx for muon
	  NTuple::Array<double> chipi;//chi2_dEdx pion
	  NTuple::Array<double> chik; //chi2_dEdx kaon
	  NTuple::Array<double> chip; //chi2_dEdx proton
	  NTuple::Array<double> ghit; //number of good de/dx hits (excluding overflow)
	  NTuple::Array<double> thit; //nuber of total de/dx including overflow
	  NTuple::Array<double> probPH; //most probable pulse height from trucated mean 
	  NTuple::Array<double> normPH; //normalized pulse height
	  NTuple::Array<double> e; //dedx for electron
	  NTuple::Array<double> mu; //muon
	  NTuple::Array<double> pi; //pion
	  NTuple::Array<double> K; //Kaon
	  NTuple::Array<double> p; //proton
  };
  NTuple::Tuple * dedx_tuple;
  DEDX_t dedx;


  struct TOF_t
  {
	  NTuple::Item<long> ntrack; //number of tracks
	  NTuple::Array<double>  trackID; //track id 
	  NTuple::Array<double>  tofID; //TOF counter ID
	  NTuple::Array<double>  tofTrackID; //Cluster (Group of TOF Hits) ID
	  NTuple::Array<double>  status; //Raw/Readout/Counter/TrackEnd/Layer/Barrel
	  NTuple::Array<double>  path;   //distance of flight
	  NTuple::Array<double>  zrhit;  //hit position , for barrel, it is z, for endcap, it is r
	  NTuple::Array<double>  ph;     //Pulse height (adc) unit: channel
	  NTuple::Array<double>  tof;    //Time of flight, unit: to be determined
	  NTuple::Array<double>  errtof; //Error of time of flight
	  NTuple::Array<double>  beta;   //Beta value of the particle
	  NTuple::Array<double>  texpe;  //Expected time of flight. electron 2 pion, 3 kaon, 4 proton
	  NTuple::Array<double>  texpmu; //Expected time of flight muon
	  NTuple::Array<double>  texppi; //Expected time of fligh pion
	  NTuple::Array<double>  texpK;  //Expected time of fligh Kaon
	  NTuple::Array<double>  texpp;  //Expected time of fligh proton
	  NTuple::Array<double>  toffsete; //Time offset of electron
	  NTuple::Array<double>  toffsetmu;//Time offset of muon
	  NTuple::Array<double>  toffsetpi;//Time offset of pion
	  NTuple::Array<double>  toffsetK; //Time offset of kaon
	  NTuple::Array<double>  toffsetp; //Time offset of proton
	  NTuple::Array<double>  toffsetap;//Time offset of anti proton
	  NTuple::Array<double>  sigmae;  //Time resolution(sigma) of electron
	  NTuple::Array<double>  sigmamu; //Time resolution(sigma) of muon
	  NTuple::Array<double>  sigmapi; //Time resolution(sigma) of pion
	  NTuple::Array<double>  sigmaK;  //Time resolution(sigma) of kaon
	  NTuple::Array<double>  sigmap;  //Time resolution(sigma) of proton
	  NTuple::Array<double>  sigmaap; //Time resolution(sigma) of antiproton
	  NTuple::Array<long>    quality; 
	  /* Data quality of reconstruction.
        0: ZT-ZTDC didnot match
        1: good charged track
        2: neutral track with good hit
        3: no hit in counter
        4: two hits in counter
        5: more than two hits in counter
        6: only single end output of one layer
        7: two hits in counter with bad match with ZTDC
        10: initialize
	   */
	  NTuple::Array<double> t0;   //Event start time
	  NTuple::Array<double> errt0;//Error of event start time
	  NTuple::Array<double> errz; //Error of hit position (for neutral track)
	  NTuple::Array<double> phi;  //Hit position, phi angle (barrel) (for neutral track)
	  NTuple::Array<double> errphi; //Error of hit position phi
	  NTuple::Array<double> E;      //Energy deposit in TOF (for neutral track)
	  NTuple::Array<double> errE;   //error of Energy deposit in TOF (for neutral track)
  };

  NTuple::Tuple * tof_tuple;
  TOF_t tof;

  /*  Monte Carlo Information */
  struct MC_t
  {
	  NTuple::Item<long> ntrack; //number of tracks
	  NTuple::Array<double>  id; //particle pdg number
    NTuple::Array<double> px; //momentum
    NTuple::Array<double> py; //momentum
    NTuple::Array<double> pz; //momentum
    NTuple::Array<double> p; //momentum
    NTuple::Array<double> E; //Energy
    NTuple::Tuple * tuple;
    void init(void);
    StatusCode init_tuple(NTuple::Tuple *, int max_track_number);
    int MAX_TRACK_NUMBER;
  };

  MC_t mc;

  void InitData(long number_charged_track, long number_neutral_track);

	//gamma-gamma annihilation selection
	long gg_event_writed;
	NTuple::Tuple * gg_tuple;
  EMC_t gg;


  /*  Averagin number of tracks */
  NTuple::Tuple * head_tuple;
  NTuple::Item<long> head_event_number;
  NTuple::Item<long> head_event_selected;
  NTuple::Item<long> head_run;
  NTuple::Item<double> head_ncharged_tracks;
  NTuple::Item<double> head_ncharged_tracks_rms;
  NTuple::Item<double> head_nneutral_tracks;
  NTuple::Item<double> head_nneutral_tracks_rms;
  NTuple::Item<double> head_ntotal_tracks;
  NTuple::Item<double> head_ntotal_tracks_rms;
  NTuple::Item<long>   ngood_pion_pairs;
  
  ibn::averager <double> nchtr_a; //averager for number of charged tracks
  ibn::averager <double> nntr_a; //averager for number of neutral tracks
  ibn::averager <double> nttr_a; //averager for number of total tracks
};

#endif
