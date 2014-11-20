/*
 * =====================================================================================
 *
 *       Filename:  JpsiKK.cxx
 *
 *    Description:  Event selection for measurement branching fraction Jpsi->K+K-
 *    via Psi(2S)->Jpsi(->K+K-)pi+pi- decay
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

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "VertexFit/IVertexDbSvc.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/ISvcLocator.h"

#include "EventModel/EventModel.h"
#include "EventModel/Event.h"

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "DstEvent/TofHitStatus.h"
#include "EventModel/EventHeader.h"

#include "McTruth/McParticle.h"
//#include "EventNavigator/EventNavigator.h"


#include "TMath.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif
#include "JpsiKK.h"
#include "Sphericity.h"


#include "VertexFit/KinematicFit.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"
#include "ParticleID/ParticleID.h"


#include <algorithm>
#include <limits>

int JpsiKK::RootEmc::ARRAY_SIZE = 100;

const double PI_MESON_MASS = 0.13957018; //GeV
const double PION_MASS     = 0.13957018; //GeV
const double MUON_MASS     = 0.105658389; //GeV
const double KAON_MASS     = 0.493677; //GeV
const double ELECTRON_MASS = 0.000510999;//GeV
const double PROTON_MASS   = 0.93827231;//GeV

const double JPSI_MASS = 3.096916; //GeV
const double PSIP_MASS = 3.686093; //GeV

enum              {ID_KAON=0, ID_MUON=1, ID_ELECTRON=2, ID_PION=3, ID_PROTON=4};
double XMASS[5] = {KAON_MASS, MUON_MASS, ELECTRON_MASS, PION_MASS, PROTON_MASS};

inline double sq(double x) { return x*x; }

typedef std::pair<EvtRecTrackIterator, EvtRecTrackIterator> TrackPair_t;
typedef std::list<TrackPair_t> TrackPairList_t;
typedef std::list<EvtRecTrackIterator> TrackList_t;
typedef std::vector<EvtRecTrackIterator> TrackVector_t;


TrackVector_t make_track_vector(TrackPair_t & pair1, TrackPair_t & pair2)
{
  TrackVector_t V(4);
  V[0] = pair1.first;
  V[1] = pair1.second;
  V[2] = pair2.first;
  V[3] = pair2.second;
  return V;
}

TrackVector_t make_track_vector(TrackPair_t & pair1)
{
  TrackVector_t V(2);
  V[0] = pair1.first;
  V[1] = pair1.second;
  return V;
}

JpsiKK::JpsiKK(const std::string& name, ISvcLocator* pSvcLocator) :
  Algorithm(name, pSvcLocator)
{
  declareProperty("MAX_CHARGED_TRACKS", MAX_CHARGED_TRACKS=4); //maximum number of charged tracks in selection
  declareProperty("MIN_CHARGED_TRACKS", MIN_CHARGED_TRACKS=4); //minimum number of charged tracks in selection
  declareProperty("MAX_NEUTRAL_TRACKS", MAX_NEUTRAL_TRACKS=0); //maximum number of good charged tracks in selection

  //good charged track configuration
  declareProperty("IP_MAX_Z", IP_MAX_Z = 10.0); //cm
  declareProperty("IP_MAX_RHO", IP_MAX_RHO = 1.0); //cm
  declareProperty("MAX_COS_THETA", MAX_COS_THETA = 0.93);

  //good neutral track configuration
  //endcup calorimeter
  declareProperty("EMC_ENDCUP_MIN_COS_THETA", EMC_ENDCUP_MIN_COS_THETA = 0.86);
  declareProperty("EMC_ENDCUP_MAX_COS_THETA", EMC_ENDCUP_MAX_COS_THETA = 0.92);
  declareProperty("EMC_ENDCUP_MIN_ENERGY", EMC_ENDCUP_MIN_ENERGY = 0.05);
  //barrel calorimeter
  declareProperty("EMC_BARREL_MAX_COS_THETA", EMC_BARREL_MAX_COS_THETA = 0.8);
  declareProperty("EMC_BARREL_MIN_ENERGY", EMC_BARREL_MIN_ENERGY = 0.025);

  declareProperty("MAX_MUON_EP_RATIO", MAX_MUON_EP_RATIO = 0.26);
  declareProperty("MAX_KAON_EP_RATIO", MAX_KAON_EP_RATIO = 0.8);

  declareProperty("MAX_PION_MOMENTUM", MAX_PION_MOMENTUM = 0.45); //GeV
  declareProperty("MIN_KAON_MOMENTUM", MIN_KAON_MOMENTUM = 1.0); //GeV
  declareProperty("MAX_KAON_MOMENTUM", MAX_KAON_MOMENTUM = 2.0); //GeV
  declareProperty("MIN_MUON_MOMENTUM", MIN_MUON_MOMENTUM = 1.0); //GeV
  declareProperty("MAX_MUON_MOMENTUM", MAX_MUON_MOMENTUM = 2.0); //GeV

  declareProperty("MIN_RECOIL_MASS", MIN_RECOIL_MASS = 3.0); //GeV
  declareProperty("MAX_RECOIL_MASS", MAX_RECOIL_MASS = 3.2); //GeV
  declareProperty("MIN_INVARIANT_MASS", MIN_INVARIANT_MASS = 3.0); //GeV
  declareProperty("MAX_INVARIANT_MASS", MAX_INVARIANT_MASS = 3.2); //GeV
  declareProperty("MIN_KAON_MISSING_MASS", MIN_KAON_MISSING_MASS = 0.1); //GeV^2
  declareProperty("MAX_KAON_MISSING_MASS", MAX_KAON_MISSING_MASS = 0.6); //GeV^2
  declareProperty("MIN_MUON_MISSING_MASS", MIN_MUON_MISSING_MASS = 0); //GeV^2
  declareProperty("MAX_MUON_MISSING_MASS", MAX_MUON_MISSING_MASS = 0.1); //GeV^2

  declareProperty("MIN_MISSING_MASS", MIN_MISSING_MASS = -0.1); //GeV^2
  declareProperty("MAX_MISSING_MASS", MAX_MISSING_MASS = +0.1); //GeV^2
}

//this is service function for fast book ntuple
template <class A>
StatusCode init_tuple(JpsiKK * alg, A & a,  const char * dir, const char * title, MsgStream & log)
{
  StatusCode status;
  NTuplePtr nt(alg->ntupleSvc(), dir);
  if(nt) a.tuple = nt;
  else
  {
    a.tuple = alg->ntupleSvc()->book(dir, CLID_ColumnWiseTuple, title);
    if(a.tuple)
    {
      return a.init_tuple();
    }
    else
    {
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(a.tuple) << endmsg;
      return StatusCode::FAILURE;
    }
  }
  return status;
}

StatusCode JpsiKK::initialize(void)
{
  MsgStream log(msgSvc(), name());

  log << MSG::INFO << "in initialize()" << endmsg;
  event_proceed=0;
  event_write = 0;
  event_with_kaons=0;
  event_with_muons=0;

  StatusCode status;
  status = init_tuple(this, fEvent,  "FILE1/event","Signal events pi+pi- K+K-, or pi+pi- mu+mu-",log);
  status = init_tuple(this, fPid,  "FILE1/pid","particle id",log);
  status = init_tuple(this, fMdc,    "FILE1/mdc","Mdc info for signal",log);
  status = init_tuple(this, fDedx,   "FILE1/dedx","Dedx info for signal",log);
  status = init_tuple(this, fEmc,    "FILE1/emc","Emc info for signal",log);
  status = init_tuple(this, fTof,    "FILE1/tof","Tof info for signal",log);
  status = init_tuple(this, fNeutral,"FILE1/neutral","Good neutral tracks",log);

  return status;
}


StatusCode JpsiKK::RootEvent::init_tuple(void)
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
  status = tuple->addItem ("channel", channel); //decay channel of the J/psi
  status = tuple->addItem ("KK", KK); //KK decay channel of the J/psi
  status = tuple->addItem ("uu", uu); //mu-mu decay channel of the J/psi
  status = tuple->addItem ("Mrec", Mrecoil); 
  status = tuple->addItem ("Minv", Minv); 
  status = tuple->addItem ("kin_chi2", kin_chi2); 
  status = tuple->addItem ("pid_chi2", pid_chi2); 

  status = tuple->addItem ("ntrack", ntrack,0,4); //array size must be = 4
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

void JpsiKK::RootEvent::init(void)
{
  ntrack=4;
}

StatusCode JpsiKK::RootPid::init_tuple(void)
{
  StatusCode status;
  status = tuple->addItem ("Mee", M[ID_ELECTRON]);
  status = tuple->addItem ("MKK", M[ID_KAON]);
  status = tuple->addItem ("Muu", M[ID_MUON]);
  status = tuple->addItem ("Mpp", M[ID_PROTON]);
  status = tuple->addItem ("Mpipi", M[ID_PION]);
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

void JpsiKK::RootPid::init(void)
{
  ntrack=4;
}

StatusCode JpsiKK::RootMdc::init_tuple(void)
{
  StatusCode status;
  status = tuple->addItem ("Mrec", Mrecoil); 
  status = tuple->addItem ("Minv", Minv); 
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


void JpsiKK::RootMdc::init(void)
{
  ntrack=4;
}


StatusCode JpsiKK::RootEmc::init_tuple(void)
{
  StatusCode status;
  status = tuple->addItem ("ntrack",       ntrack,0,ARRAY_SIZE); //good nuetral track in event
  status = tuple->addIndexedItem ("E",     ntrack, E);
  status = tuple->addIndexedItem ("theta", ntrack, theta);
  status = tuple->addIndexedItem ("phi",   ntrack, phi);
  status = tuple->addIndexedItem ("time",  ntrack, time);
  return status;
}

void JpsiKK::RootEmc::init(void)
{
  ntrack=4;
  for(int i=0;i<ntrack;i++)
  {
    E[i] = 0;
    theta[i] = -1000;;
    phi[i] = -1000;
    time[i] = -1000;
  }
}

StatusCode JpsiKK::RootDedx::init_tuple(void)
{
  StatusCode status;
  status = tuple->addItem ("ntrack", ntrack,0,4); 
  status = tuple->addIndexedItem ("chie",  ntrack, chie);
  status = tuple->addIndexedItem ("chimu",  ntrack, chimu);
  status = tuple->addIndexedItem ("chipi",  ntrack, chipi);
  status = tuple->addIndexedItem ("chik",  ntrack, chik);
  status = tuple->addIndexedItem ("chip",  ntrack, chip);
  status = tuple->addIndexedItem ("probPH",  ntrack, probPH);
  status = tuple->addIndexedItem ("normPH",  ntrack, normPH);
  //status = tuple->addIndexedItem ("probe",  ntrack, probe);
  //status = tuple->addIndexedItem ("probmu",  ntrack, probmu);
  //status = tuple->addIndexedItem ("probpi",  ntrack, probpi);
  //status = tuple->addIndexedItem ("probk",  ntrack, probk);
  //status = tuple->addIndexedItem ("probp",  ntrack, probp);
  return status;
}

void JpsiKK::RootDedx::init(void)
{
  ntrack=4;
  for(int i=0;i<ntrack;i++)
  {
    chie[i]=-1000;
    chimu[i]=-1000;
    chipi[i]=-1000;
    chik[i]=-1000;
    chip[i]=-1000;
    probPH[i]=0;
    normPH[i]=0;
  }
}

StatusCode JpsiKK::RootTof::init_tuple(void)
{
  StatusCode status;
  status = tuple->addItem ("ntrack", ntrack,0,4); 
  status = tuple->addIndexedItem ("ID",  ntrack, tofID);
  status = tuple->addIndexedItem ("t",  ntrack, t);
  status = tuple->addIndexedItem ("dt",  ntrack, dt);
  status = tuple->addIndexedItem ("t0",  ntrack, t0);
  status = tuple->addIndexedItem ("chie",  ntrack, chie);
  status = tuple->addIndexedItem ("chimu",  ntrack, chimu);
  status = tuple->addIndexedItem ("chipi",  ntrack, chipi);
  status = tuple->addIndexedItem ("chik",  ntrack, chik);
  status = tuple->addIndexedItem ("chip",  ntrack, chip);
  status = tuple->addIndexedItem ("beta",  ntrack, beta);
  status = tuple->addIndexedItem ("te",  ntrack, te);
  status = tuple->addIndexedItem ("tmu",  ntrack, tmu);
  status = tuple->addIndexedItem ("tpi",  ntrack, tpi);
  status = tuple->addIndexedItem ("tk",  ntrack, tk);
  status = tuple->addIndexedItem ("tp",  ntrack, tp);
  return status;
}

void JpsiKK::RootTof::init(void)
{
  ntrack=4;
  for(int i=0;i<ntrack;i++)
  {
    tofID[i]=-1000;
    t[i]=-1000;
    dt[i]=-1000;
    t0[i]=-1000;
    chie[i]=-1000;
    chimu[i]=-1000;
    chipi[i]=-1000;
    chik[i]=-1000;
    chip[i]=-1000;
    beta[i]=-1000;
    te[i]=-1000;
    tmu[i]=-1000;
    tpi[i]=-1000;
    tk[i]=-1000;
    tp[i]=-1000;
  }
}


void calculate_vertex(RecMdcTrack *mdcTrk, double & ro, double  & z, double phi)
{
  ro = -9999;
  z = -9999;
  phi = -9999;
  /*  Reconstruct the vertex */
  Hep3Vector xorigin(0,0,0);
  IVertexDbSvc*  vtxsvc;
  Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
  if(vtxsvc->isVertexValid())
  {
    double* dbv = vtxsvc->PrimaryVertex(); 
    double*  vv = vtxsvc->SigmaPrimaryVertex();  
    xorigin.setX(dbv[0]);
    xorigin.setY(dbv[1]);
    xorigin.setZ(dbv[2]);
  }
  /* Vertex game. copy from rhophi analysis */
  double phi0=mdcTrk->helix(1);
  double xv=xorigin.x();
  double yv=xorigin.y();
  //double Rxy=(mdc.x[i]-xv)*cos(phi0)+(mdc.y[i]-yv)*sin(phi0);
  //mdc.r[i]=Rxy;
  HepVector a = mdcTrk->helix();
  HepSymMatrix Ea = mdcTrk->err();
  HepPoint3D point0(0.,0.,0.);   // the initial point for MDC recosntruction
  HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]); 
  VFHelix helixip(point0,a,Ea); 
  helixip.pivot(IP);
  HepVector vecipa = helixip.a();
  double  Rvxy0=fabs(vecipa[0]);  //the nearest distance to IP in xy plane
  double  Rvz0=vecipa[3];         //the nearest distance to IP in z direction
  double  Rvphi0=vecipa[1];
  ro=Rvxy0;
  z=Rvz0;
  phi=Rvphi0;
}


double get_invariant_mass2(TrackPair_t & pair, double mass)
{
  EvtRecTrackIterator  itTrk[2] = {pair.first, pair.second};
  HepLorentzVector  P[2];
  for(int k=0;k<2;k++)
  {
    if(!(*itTrk[k])->isMdcTrackValid()) break; 
    RecMdcTrack *mdcTrk = (*itTrk[k])->mdcTrack();
    P[k] = mdcTrk->p4(mass);
  }
  HepLorentzVector P_sum = P[0]+P[1];
  return P_sum.m2();
}

double get_recoil__mass(EvtRecTrackIterator & trk1, EvtRecTrackIterator & trk2, double mass)
{
  EvtRecTrackIterator  itTrk[2] = {trk1, trk2};
  HepLorentzVector  P[2];
  for(int k=0;k<2;k++)
  {
    if(!(*itTrk[k])->isMdcTrackValid()) break; 
    RecMdcTrack *mdcTrk = (*itTrk[k])->mdcTrack();
    P[k] = mdcTrk->p4(mass);
  }
  HepLorentzVector P_psip(0.040546,0,0,PSIP_MASS); //initial vector of psip
  HepLorentzVector P_sum = P[0]+P[1];
  HepLorentzVector P_recoil = P_psip - P_sum;
  return P_recoil.m();
}

double get_recoil__mass(std::pair<EvtRecTrackIterator, EvtRecTrackIterator> p, double mass)
{
  return get_recoil__mass(p.first, p.second, mass);
}


double get_missing_mass(std::pair<EvtRecTrackIterator, EvtRecTrackIterator> pions, std::pair<EvtRecTrackIterator, EvtRecTrackIterator> kaons)
{
  EvtRecTrackIterator  PionTrk[2] = {pions.first, pions.second};
  EvtRecTrackIterator  KaonTrk[2] = {kaons.first, kaons.second};
  HepLorentzVector P_psip(0.040546,0,0,PSIP_MASS); //initial vector of psip
  HepLorentzVector  pionP[2];
  HepLorentzVector  kaonP[2];
  for(int k=0;k<2;k++)
  {
    pionP[k] = HepLorentzVector(0,0,0,0);
    if(!(*PionTrk[k])->isMdcTrackValid()) break; 
    RecMdcTrack *mdcTrk = (*PionTrk[k])->mdcTrack();
    pionP[k] = mdcTrk->p4(PION_MASS);
  }
  for(int k=0;k<2;k++)
  {
    kaonP[k] = HepLorentzVector(0,0,0,0);
    if(!(*KaonTrk[k])->isMdcTrackValid()) break; 
    RecMdcTrack *mdcTrk = (*KaonTrk[k])->mdcTrack();
    kaonP[k] = mdcTrk->p4(KAON_MASS);
  }
  HepLorentzVector Pmis = P_psip - pionP[0] - pionP[1] - kaonP[0] - kaonP[1];
  return Pmis.m2();
}

SmartRefVector<RecTofTrack>::iterator  getTofTrk(EvtRecTrackIterator itTrk, bool & isTofValid)
{
  SmartRefVector<RecTofTrack> tofTrkCol = (*itTrk)->tofTrack();
  SmartRefVector<RecTofTrack>::iterator tofTrk = tofTrkCol.begin();
  TofHitStatus *hitst = new TofHitStatus;
  std::vector<int> tofecount;
  int goodtofetrk=0;
  for(tofTrk = tofTrkCol.begin(); tofTrk!=tofTrkCol.end(); tofTrk++,goodtofetrk++)
  {
    unsigned int st = (*tofTrk)->status();
    hitst->setStatus(st);
    //if( !hitst->is_claster() ) continue;
    //if(  (hitst->is_barrel()) ) continue;
    if( !(hitst->is_counter()) ) continue;
    tofecount.push_back(goodtofetrk);
  }
  delete hitst;
  if(!tofecount.empty()) 
  {
    tofTrk = tofTrkCol.begin()+tofecount[0];
    isTofValid = true;
  }
  return tofTrk;
}

//vector<SmartRefVector<RecTofTrack>::iterator>  getTofTrk(EvtRecTrackIterator itTrk, bool & isTofValid)
//{
//  SmartRefVector<RecTofTrack> tofTrkCol = (*itTrk)->tofTrack();
//  SmartRefVector<RecTofTrack>::iterator tofTrk = tofTrkCol.begin();
//  vector<SmartRefVector<RecTofTrack>::iterator> result_tof;
//  TofHitStatus *hitst = new TofHitStatus;
//  //std::vector<int> tofecount;
//  //int goodtofetrk=0;
//  for(tofTrk = tofTrkCol.begin(); tofTrk!=tofTrkCol.end(); tofTrk++,goodtofetrk++)
//  {
//    unsigned int st = (*tofTrk)->status();
//    hitst->setStatus(st);
//    //if(  (hitst->is_barrel()) ) continue;
//    if( !(hitst->is_counter()) ) continue;
//    tofecount.push_back(goodtofetrk);
//  }
//  delete hitst;
//  if(!tofecount.empty()) 
//  {
//    tofTrk = tofTrkCol.begin()+tofecount[0];
//    isTofValid = true;
//  }
//  return tofTrk;
//}


vector<double> get_chi2(EvtRecTrackIterator & itTrk)
{
  vector<double> chi2(5,99999);
  if(!(*itTrk)->isMdcTrackValid()) return chi2;
  if(!(*itTrk)->isMdcDedxValid())  return chi2;
  for(int i=0;i<5;i++) chi2[i]=0;
  RecMdcTrack * mdcTrk  = (*itTrk)->mdcTrack();

  //dedx information
  RecMdcDedx  * dedxTrk = (*itTrk)->mdcDedx();
  chi2[ID_KAON]     +=   sq(dedxTrk->chiK());
  chi2[ID_MUON]     +=   sq(dedxTrk->chiMu());
  chi2[ID_ELECTRON] +=   sq(dedxTrk->chiE());
  chi2[ID_PION]     +=   sq(dedxTrk->chiPi());
  chi2[ID_PROTON]   +=   sq(dedxTrk->chiP());

  //tof information
  if(!(*itTrk)->isTofTrackValid()) return chi2;
  bool isTofValid=false;
  SmartRefVector<RecTofTrack>::iterator tofTrk = getTofTrk(itTrk, isTofValid);
  if(isTofValid)
  {
    double t = (*tofTrk)->tof();  //flight time
    double dt = (*tofTrk)->errtof(); //error of flight time
    chi2[ID_KAON]     +=   sq(((*tofTrk)->texpKaon()-t)/dt);
    chi2[ID_MUON]     +=   sq(((*tofTrk)->texpMuon()-t)/dt);
    chi2[ID_ELECTRON] +=   sq(((*tofTrk)->texpElectron()-t)/dt);
    chi2[ID_PION]     +=   sq(((*tofTrk)->texpPion()-t)/dt);
    chi2[ID_PROTON]   +=   sq(((*tofTrk)->texpProton()-t)/dt);
  }
  return chi2;
}


vector< vector<double> > get_chi2(TrackPair_t & pion_pair, TrackPair_t & kaon_pair)
{
  EvtRecTrackIterator  itTrk[4] = {pion_pair.first, pion_pair.second, kaon_pair.first, kaon_pair.second};
  vector< vector<double> > chi2(4);
  for(int track=0;track<chi2.size();track++)
  {
    chi2[track] = get_chi2(itTrk[track]);
  }
  return chi2;
}

vector<double> get_chi2(std::pair<EvtRecTrackIterator, EvtRecTrackIterator> & pair)
{
  EvtRecTrackIterator  itTrk[2] = {pair.first, pair.second};
  vector<double> chi2(5,0);
  for(int track=0;track<2;track++)
  {
    vector<double> tmp_chi2=get_chi2(itTrk[track]);
    for(int i=0;i<5;i++)
    {
      chi2[i]  += tmp_chi2[i]; 
    }
  }
  return chi2;
}

bool kinematic_fit(int PID, TrackPair_t  & pion_pair, TrackPair_t &  other_pair, std::vector<HepLorentzVector> & P, double & chi2)
{
  P.resize(4);
  EvtRecTrackIterator  PionTrk[2] = {pion_pair.first, pion_pair.second};
  EvtRecTrackIterator  OtherTrk[2] = {other_pair.first, other_pair.second};
  RecMdcKalTrack * PionKalTrk[2];
  RecMdcKalTrack * OtherKalTrk[2];
  WTrackParameter PionWTrk[2];
  WTrackParameter OtherWTrk[2];

  for(int i=0;i<2;i++)
  {
    PionKalTrk[i] = (*PionTrk[i])->mdcKalTrack();
    OtherKalTrk[i] = (*OtherTrk[i])->mdcKalTrack();
    PionWTrk[i] = WTrackParameter(PION_MASS, PionKalTrk[i]->getZHelix(), PionKalTrk[i]->getZError());
    switch(PID)
    {
      case ID_KAON:
        OtherWTrk[i] = WTrackParameter(XMASS[PID], OtherKalTrk[i]->getZHelixK(), OtherKalTrk[i]->getZErrorK());
        break;
      case ID_MUON:
        OtherWTrk[i] = WTrackParameter(XMASS[PID], OtherKalTrk[i]->getZHelixMu(), OtherKalTrk[i]->getZErrorMu());
        break;
      case ID_ELECTRON:
        OtherWTrk[i] = WTrackParameter(XMASS[PID], OtherKalTrk[i]->getZHelixE(), OtherKalTrk[i]->getZErrorE());
        break;
      case ID_PION:
        OtherWTrk[i] = WTrackParameter(XMASS[PID], OtherKalTrk[i]->getZHelix(), OtherKalTrk[i]->getZError());
        break;
      case ID_PROTON:
        OtherWTrk[i] = WTrackParameter(XMASS[PID], OtherKalTrk[i]->getZHelixP(), OtherKalTrk[i]->getZErrorP());
        break;
    }
  }
  //now try the kinematic fit
  //initial vertex
  HepPoint3D vx(0., 0., 0.);
  cout << "x: " << PionWTrk[0].x().x() << " " << PionWTrk[0].x().y() << " " << PionWTrk[0].x().z() << endl;
  //error matrix inital valu
  HepSymMatrix Evx(3, 0);
  double bx = 1E+6;
  double by = 1E+6;
  double bz = 1E+6;
  Evx[0][0] = bx*bx;
  Evx[1][1] = by*by;
  Evx[2][2] = bz*bz;
  VertexParameter vxpar;
  vxpar.setVx(vx);
  vxpar.setEvx(Evx);


  //Vetex fitter
  VertexFit* vtxfit = VertexFit::instance();
  vtxfit->init();
  //add tracks. I know the first two one must be pions
  vtxfit->AddTrack(0,  PionWTrk[0]);
  vtxfit->AddTrack(1,  PionWTrk[1]);
  vtxfit->AddTrack(2,  OtherWTrk[0]);
  vtxfit->AddTrack(3,  OtherWTrk[1]);
  vtxfit->AddVertex(0, vxpar,0, 1, 2,3);
  //if(!vtxfit->Fit(0)) return false;
  //if(!vtxfit->Fit(1)) return false;
  vtxfit->Fit();
  //vtxfit->Swim(0);
  cout << "After vertex fit: x: " << vtxfit->wtrk(0).x().x() << " " << vtxfit->wtrk(0).x().y() << " " << vtxfit->wtrk(0).x().z() << endl;

  //KinematicFit * kmfit = KinematicFit::instance();
  KalmanKinematicFit * kmfit = KalmanKinematicFit::instance();
  kmfit->init();
  for(int i=0;i<4;i++)
  {
    kmfit->AddTrack(i,vtxfit->wtrk(i));
  }
  HepLorentzVector Pcmf(0.040546,0,0,PSIP_MASS); //initial vector of center of mass frame
  kmfit->AddFourMomentum(0,  Pcmf);
  //kmfit->AddResonance(0, PSIP_MASS, 0, 1, 2,3);
  //kmfit->AddResonance(1, JPSI_MASS, 2, 3);
  //if(!kmfit->Fit(0)) return false;
  //kmfit->Fit(0);
  bool oksq = kmfit->Fit();
  if(oksq) 
  {
    chi2  = kmfit->chisq();
    for(int i=0;i<4;i++)
    {
      P[i] = kmfit->pfit(i);
    }
  }
  return oksq;
}

bool kinematic_fit(int PID, TrackPairList_t  & pion_pairs, TrackPairList_t &  other_pairs, std::vector<HepLorentzVector> & P, double & chi2, TrackPair_t & result_pion_pair, TrackPair_t & result_other_pair)
{
  if(pion_pairs.empty() || other_pairs.empty()) return false;
  chi2=1e100;
  bool GoodKinematikFit=false;

  for(TrackPairList_t::iterator pion_pair=pion_pairs.begin(); pion_pair!=pion_pairs.end();pion_pair++)
    for(TrackPairList_t::iterator other_pair=other_pairs.begin(); other_pair!=other_pairs.end();other_pair++)
    {
      std::vector<HepLorentzVector> P_tmp;
      double chi2_tmp=3e100;
      bool oksq=kinematic_fit(PID, *pion_pair, *other_pair, P_tmp, chi2_tmp);
      if(oksq) 
      {
        if(chi2_tmp < chi2)
        {
          GoodKinematikFit = true;
          chi2  = chi2_tmp;
          P = P_tmp;
          result_pion_pair = *pion_pair;
          result_other_pair = * other_pair;
        }
      }
    } 
  return GoodKinematikFit;
}

StatusCode JpsiKK::execute()
{
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "executing" << endreq;

  SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
  fEvent.run=eventHeader->runNumber();
  fEvent.event=eventHeader->eventNumber();
  fEvent.time=eventHeader->time();
  bool isprint=false;
  if(event_proceed<10) isprint=true;
  if(10 <= event_proceed && event_proceed < 100 && event_proceed % 10 ==0) isprint=true;
  if(100 <= event_proceed && event_proceed < 1000 && event_proceed % 100 ==0) isprint = true;
  if(1000 <= event_proceed && event_proceed < 10000 && event_proceed % 1000 ==0) isprint = true;
  if(10000 <= event_proceed && event_proceed % 10000 ==0) isprint = true;
  if(isprint)
  {
    std::cout << "proceed event: " << setw(15) << event_proceed<<", ";
    std::cout << "event write: "   << setw(15) << event_write<< ",  ";
    std::cout << "kaons: "   << setw(15) << event_with_kaons << ",  ";
    std::cout << "muons "   << setw(15) << event_with_muons;
    std::cout << std::endl;
  }
  event_proceed++;

  //  Get information about reconstructed events
  SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
  SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol);
  SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(),  EventModel::MC::McParticleCol);

//  SmartDataPtr<EventNavigator> navigator (eventSvc(),"/Event/Navigator");
//  if( ! navigator )
//    {
//      log << MSG::ERROR << " Unable to retrieve EventNavigator" << endreq;
//      return StatusCode::FAILURE;
//    }
//
  //log << MSG::INFO << "EventNavigator object" << endreq;
  //navigator->Print();
  

  //fill initial value of the selected event
  fEvent.init();

  //list of good charged tracks
  std::list<EvtRecTrackIterator> good_charged_tracks;
  //list of good netutral tracks
  std::list<EvtRecTrackIterator> good_neutral_tracks;

  //collect  good charged track
  for(unsigned i = 0; i < evtRecEvent->totalCharged(); i++)
  {
    EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
    if(!(*itTrk)->isMdcTrackValid()) continue;  //use only valid charged tracks
    RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();  //main drift chambe
    //calculate interaction point distance
    double rvxy,rvz,rvphi;
    calculate_vertex(mdcTrk,rvxy,rvz,rvphi); //find distance to interaction point
    bool IP_track = fabs(rvz)< IP_MAX_Z && fabs(rvxy)<IP_MAX_RHO;  //tracks begin near interaction point
    bool good_track = IP_track && fabs(cos(mdcTrk->theta()))<MAX_COS_THETA; //track is good
    if(good_track) good_charged_tracks.push_back(itTrk);
  }

  //collect good neutral track
  for(int i = evtRecEvent->totalCharged(); i<evtRecEvent->totalTracks(); i++)
  {
    EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
    if(!(*itTrk)->isEmcShowerValid()) continue; //keep only valid neutral tracks
    RecEmcShower *emcTrk = (*itTrk)->emcShower();
    double c =  fabs(cos(emcTrk->theta())); //abs cos theta
    double E  =  emcTrk->energy();
    bool hit_barrel = (c <= EMC_BARREL_MAX_COS_THETA);
    bool hit_endcup = (EMC_ENDCUP_MIN_COS_THETA <=c) && (c <= EMC_ENDCUP_MAX_COS_THETA);
    //barrel and endcup calorimeters have different energy threshold
    bool barrel_good_track = hit_barrel && (E > EMC_BARREL_MIN_ENERGY);
    bool endcup_good_track = hit_endcup && (E > EMC_ENDCUP_MIN_ENERGY);
    if(barrel_good_track  || endcup_good_track) 
    {
      //cout << "Energy of good neutral track: " << E << endl;
      good_neutral_tracks.push_back(itTrk);
    }
  }

  //print good charged track index
  //cout << "Good charged track: ";
  //for(list<EvtRecTrackIterator>::iterator i=good_charged_tracks.begin();i!=good_charged_tracks.end();i++)
  //{
  //  cout << *i - evtRecTrkCol->begin() << " ";
  //}
  //cout << endl;

  
  //log << MSG::ERROR << "MAX_NEUTRAL_TRACKS  = " << MAX_NEUTRAL_TRACKS << endmsg;
  //log << MSG::ERROR << "MIN_CHARGED_TRACKS  = " << MIN_CHARGED_TRACKS << endmsg;
  //log << MSG::ERROR << "MAX_CHARGED_TRACKS  = " << MAX_CHARGED_TRACKS << endmsg;
  //log << MSG::ERROR << "good charged tracks: " << good_charged_tracks.size() <<  ",  neutral tracks: " << good_neutral_tracks.size() << endmsg;

  //SELECTION CODE
  if( MAX_NEUTRAL_TRACKS < good_neutral_tracks.size()) return StatusCode::SUCCESS;
  //SELECTION CODE
  if( good_charged_tracks.size() < MIN_CHARGED_TRACKS || MAX_CHARGED_TRACKS < good_charged_tracks.size()) return StatusCode::SUCCESS;


  
  TrackList_t charged_tracks; //selected tracks for specific cut
  TrackList_t positive_charged_tracks; //selected tracks for specific cut
  TrackList_t negative_charged_tracks; //selected tracks for specific cut
  TrackList_t positive_pion_tracks; //selected pion tracks for specific cut
  TrackList_t negative_pion_tracks; //selected pion tracks for specific cut
  TrackList_t other_positive_tracks; //other positive tracks for specific cut
  TrackList_t other_negative_tracks; //other positive tracks for specific cut
  for(TrackList_t::iterator track=good_charged_tracks.begin(); track!=good_charged_tracks.end(); track++)
  {
    EvtRecTrackIterator & itTrk = *track;
    if(!(*itTrk)->isMdcTrackValid()) continue; 
    RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
    double c = fabs(cos(mdcTrk->theta()));
    double p = mdcTrk->p();
    double q = mdcTrk->charge();
    bool barrel = c < EMC_BARREL_MAX_COS_THETA;
    if(barrel) 
    {
      if(q>0) 
      {
        positive_charged_tracks.push_back(itTrk);
        if(p<MAX_PION_MOMENTUM) 
        {
          positive_pion_tracks.push_back(itTrk);
        }
        if(p>std::min(MIN_KAON_MOMENTUM, MIN_MUON_MOMENTUM))
        {
          if((*itTrk)->isEmcShowerValid())
          {
            other_positive_tracks.push_back(itTrk);
          }
        }
      }
      if(q<0) 
      {
        negative_charged_tracks.push_back(itTrk);
        if(p<MAX_PION_MOMENTUM) 
        {
          negative_pion_tracks.push_back(itTrk);
        }
        if(p>std::min(MIN_KAON_MOMENTUM, MIN_MUON_MOMENTUM))
        {
          if((*itTrk)->isEmcShowerValid())
          {
            other_negative_tracks.push_back(itTrk);
          }
        }
      }
      charged_tracks.push_back(itTrk);
    }
  }


  //SELECTION CODE
  //keep only specific signature
  if(positive_charged_tracks.size()!=2 || negative_charged_tracks.size()!=2) return StatusCode::SUCCESS;
  //SELECTION CODE
  if(negative_pion_tracks.empty() || positive_pion_tracks.empty()) return StatusCode::SUCCESS;

  //typedef std::list< std::pair<EvtRecTrackIterator, EvtRecTrackIterator> > PairList_t;
  //create pion pairs
  TrackPairList_t pion_pairs;
  for(TrackList_t::iterator i=negative_pion_tracks.begin(); i!=negative_pion_tracks.end(); ++i)
    for(TrackList_t::iterator j=positive_pion_tracks.begin(); j!=positive_pion_tracks.end(); ++j)
    {
      TrackPair_t pair(*i,*j);
      double M_recoil = get_recoil__mass(pair, PION_MASS);
      if(MIN_RECOIL_MASS < M_recoil && M_recoil < MAX_RECOIL_MASS) 
      {
        pion_pairs.push_back(pair);
      }
    }
  //SELECTION CODE
  if(pion_pairs.empty()) return StatusCode::SUCCESS;
  //log << MSG::ERROR << "good charged tracks: " << charged_tracks.size() << " (" << negative_charged_tracks.size() << ", " << positive_charged_tracks.size() << endmsg;
  //log << MSG::ERROR << "pions: " << negative_pion_tracks.size()  << ", " << positive_pion_tracks.size() << endmsg;
  //log << MSG::ERROR << "other: " << other_negative_tracks.size()  << ", " << other_positive_tracks.size() << endmsg;
  //log << MSG::ERROR << "pion pairs: " << pion_pairs.size() << endmsg;

  //find the best pion pair using closest value to JPSI_MASS
  TrackPair_t pion_pair = pion_pairs.front();
  for(TrackPairList_t::iterator p=pion_pairs.begin();p!=pion_pairs.end();p++)
  {
    if(fabs(get_recoil__mass(*p,PION_MASS) - JPSI_MASS) <  fabs(get_recoil__mass(pion_pair,PION_MASS) - JPSI_MASS)) pion_pair = *p;
  }

  //make kaon or muon pairs
  TrackPairList_t muon_pairs;
  TrackPairList_t kaon_pairs;
  TrackPairList_t other_pairs;
  for(TrackList_t::iterator i=other_negative_tracks.begin(); i!=other_negative_tracks.end(); ++i)
    for(TrackList_t::iterator j=other_positive_tracks.begin(); j!=other_positive_tracks.end(); ++j)
    {
      TrackPair_t pair(*i,*j);
      EvtRecTrackIterator  itTrk[2] = {pair.first, pair.second};
      double Ep[2]; // E/p ratio
      double E[2];
      double p[2];
      for(int k=0;k<2;k++)
      {
        if(!(*itTrk[k])->isMdcTrackValid() || ! (*itTrk[k])->isEmcShowerValid()) 
        {
          log << MSG::ERROR << "Invalid mdc info for track.Exiting" << endmsg;
          return StatusCode::FAILURE;
        }
        //SELECTION CODE:
        RecMdcTrack *mdcTrk = (*itTrk[k])->mdcTrack();
        RecEmcShower *emcTrk = (*itTrk[k])->emcShower();
        E[k] = emcTrk->energy();
        p[k] = mdcTrk->p();
        if(p[k]!=0) Ep[k] = E[k]/p[k];
        else Ep[k]=999999;
      }
      double M[5]={0,0,0,0,0};
      for(int pid=0;pid<5;pid++)
      {
        M[pid]=get_invariant_mass2(pair,XMASS[pid]);
        if(M[pid]>0) M[pid] = sqrt(M[pid]);
        else M[pid] = 0;
      }
      //SELECTION CODE KAON CASE
      if(MIN_INVARIANT_MASS <  M[0]   && M[0]  < MAX_INVARIANT_MASS)
      {
        if(Ep[0] < MAX_KAON_EP_RATIO && Ep[1] < MAX_KAON_EP_RATIO)
        {
          if(MIN_KAON_MOMENTUM < p[0] && p[0] < MAX_KAON_MOMENTUM)
          {
            if(MIN_KAON_MOMENTUM < p[1] && p[1] < MAX_KAON_MOMENTUM)
            {
              kaon_pairs.push_back(pair);
            }
          }
        }
      }
      //SELECTION CODE
      if(MIN_INVARIANT_MASS <  M[1]   && M[1]  < MAX_INVARIANT_MASS)
      {
        if(Ep[0] < MAX_MUON_EP_RATIO && Ep[1] < MAX_MUON_EP_RATIO)
        {
          if(MIN_MUON_MOMENTUM < p[0] && p[0] < MAX_MUON_MOMENTUM)
          {
            if(MIN_MUON_MOMENTUM < p[1] && p[1] < MAX_MUON_MOMENTUM)
            {
              muon_pairs.push_back(pair);
            }
          }
        }
      }
      //SELECT ANY OTHER PAIRS
      //suppress electron case
      if(Ep[0] < MAX_KAON_EP_RATIO && Ep[1] < MAX_KAON_EP_RATIO)
      {
        if(MIN_MUON_MOMENTUM < p[0] && p[0] < MAX_MUON_MOMENTUM)
        {
          if(MIN_MUON_MOMENTUM < p[1] && p[1] < MAX_MUON_MOMENTUM)
          {
            other_pairs.push_back(pair);
          }
        }
      }
    }

  int channel=-1; //default no channel
  TrackPair_t result_pair;
  //SELECTION CODE
  if(!muon_pairs.empty() || !kaon_pairs.empty()) 
  {
    //the best pair which is closer to JPSI
    if(!kaon_pairs.empty()) result_pair = kaon_pairs.front();
    if(!muon_pairs.empty()) result_pair = muon_pairs.front();
    for(TrackPairList_t::iterator p=kaon_pairs.begin();p!=kaon_pairs.end();p++)
    {
      if(fabs(sqrt(get_invariant_mass2(*p,KAON_MASS)) - JPSI_MASS) 
          <=  fabs(sqrt(get_invariant_mass2(result_pair,KAON_MASS)) - JPSI_MASS)) 
      {
        result_pair = *p;
        channel=ID_KAON; //setup kaon channel
      }
    }
    for(TrackPairList_t::iterator p=muon_pairs.begin();p!=muon_pairs.end();p++)
    {
      if(fabs(sqrt(get_invariant_mass2(*p,MUON_MASS)) - JPSI_MASS) 
          <=  fabs(sqrt(get_invariant_mass2(result_pair,MUON_MASS)) - JPSI_MASS)) 
      {
        result_pair = *p;
        channel=ID_MUON; //setup muon channel
      }
    }
    if(channel<0) 
    {
      log << MSG::ERROR << "Must be some channel but it's not" << endmsg;
      return StatusCode::FAILURE; 
    }
    vector<double> chi2 = get_chi2(result_pair);
    bool muon_probable=true;
    bool kaon_probable=true;
    for(int i=0;i<5;i++)
    {
      if(i!=ID_KAON) kaon_probable = kaon_probable && chi2[ID_KAON] < chi2[i];
      if(i!=ID_MUON) muon_probable = muon_probable && chi2[ID_MUON] < chi2[i];
    }
    //switch(channel)
    //{
    //  case ID_KAON:
    //    if(kaon_probable)
    //    {
    //      fEvent.KK=1;
    //      fEvent.uu=0;
    //      event_with_kaons++;
    //    }
    //    break;
    //  case ID_MUON:
    //    if(muon_probable)
    //    {
    //      fEvent.KK=0;
    //      fEvent.uu=1;
    //      event_with_muons++;
    //    }
    //    break;
    //}
  }
  else
  {
    return StatusCode::SUCCESS;
  }


  if(other_pairs.empty()) return StatusCode::SUCCESS;
  bool GoodKinematikFit=false;
  double kinematic_chi2=2e100;
  std::vector<HepLorentzVector> Pkf(4);
  for(int pid = 0;pid<5;pid++)
  {
    std::vector<HepLorentzVector> P_tmp(4);
    TrackPair_t pion_pr;
    TrackPair_t other_pr;
    double chi2_tmp;
    bool fit_result = kinematic_fit(pid, pion_pairs, other_pairs, P_tmp, chi2_tmp, pion_pr, other_pr);
    cout << fit_result << " " << pid << " " << P_tmp[0].rho() << " " << P_tmp[1].rho() << " " << P_tmp[2].rho() << " " << P_tmp[3].rho() << endl;
    if(fit_result)
    {
      GoodKinematikFit = true;
      if(chi2_tmp<kinematic_chi2)
      {
        //save this better fit result
        pion_pair=pion_pr;
        result_pair = other_pr;
        channel = pid;
        kinematic_chi2 = chi2_tmp;
        Pkf=P_tmp;
      }
    }
  }
  if(!GoodKinematikFit) return StatusCode::SUCCESS;
  switch(channel)
  {
    case ID_KAON:
        fEvent.KK=1;
        fEvent.uu=0;
        event_with_kaons++;
      break;
    case ID_MUON:
        fEvent.KK=0;
        fEvent.uu=1;
        event_with_muons++;
      break;
    default:
      return StatusCode::SUCCESS;
      break;
  }
  //now we have best pion_pair, best kaon/muon pair (result_pair), four-momentum
  //of all particles after kinematik fit


  //now fill the tuples

  //some statistics information
  fEvent.ngood_charged_track = good_charged_tracks.size();
  fEvent.ngood_neutral_track = good_neutral_tracks.size();
  fEvent.npositive_track = positive_charged_tracks.size();
  fEvent.nnegative_track = negative_charged_tracks.size();
  fEvent.npositive_pions = positive_pion_tracks.size();
  fEvent.nnegative_pions = negative_pion_tracks.size();
  fEvent.npion_pairs = pion_pairs.size();
  // fill the decay channel of the J/psi 0 - kaons, 1 --muons
  fEvent.channel = channel; 
  fEvent.kin_chi2 = kinematic_chi2;
  fEvent.pid_chi2 = get_chi2(result_pair)[channel];
  fEvent.Minv = (Pkf[2]+Pkf[3]).m();
  HepLorentzVector P_psip(0.040546,0,0,PSIP_MASS); //initial vector of psip
  fEvent.Mrecoil = (P_psip - Pkf[0] - Pkf[1]).m();

  fEvent.ntrack = 4;
  for ( int i=0;i<4;i++)
  {
    fEvent.q[i]  = i%2 == 0 ? -1 : +1;
    fEvent.E[i]  = Pkf[i].e();
    fEvent.px[i] = Pkf[i].px();
    fEvent.py[i] = Pkf[i].py();
    fEvent.pz[i] = Pkf[i].pz();
    fEvent.p[i]  = sqrt(sq(Pkf[i].px())+sq(Pkf[i].py())+sq(Pkf[i].pz()));
    fEvent.pt[i] = sqrt(sq(Pkf[i].px())+sq(Pkf[i].py()));
    fEvent.theta[i]= Pkf[i].theta();
    fEvent.phi[i] = Pkf[i].phi();
    fEvent.x[i]=0;
    fEvent.y[i]=0;
    fEvent.z[i]=0;
    fEvent.r[i]=0;
    fEvent.vxy[i]=0;
    fEvent.vz[i]=0;
    fEvent.vphi[i]=0;
  }


  ParticleID * PID = ParticleID::instance();
  PID->init();
  PID->setMethod(PID->methodProbability());
  PID->setChiMinCut(4);
  PID->usePidSys(PID->useDedx() || PID->useTof());
  PID->identify(PID->all()); 

  fEvent.ntrack=4;
  fPid.ntrack=4;
  fMdc.ntrack=4;
  fDedx.ntrack=4;
  fEmc.ntrack=4;
  fTof.ntrack=4;
  fMdc.Mrecoil = get_recoil__mass(pion_pair, PION_MASS);
  fMdc.Minv    = sqrt(get_invariant_mass2(result_pair,XMASS[channel]));
  EvtRecTrackIterator itTrk[4] = {pion_pair.first, pion_pair.second, result_pair.first, result_pair.second};
  for(int i=0;i<4;i++)
  {
    if(!(*itTrk[i])->isMdcTrackValid() || !(*itTrk[i])->isEmcShowerValid() ) 
    {
      //pions could not have EMC information
      if(i>1)
      {
        log << MSG::ERROR << "Somthing wrong in selection. This track must have Mdc and Emc information but it is not. exiting." << endmsg;
        return StatusCode::FAILURE;
      }
    }
    if((*itTrk[i])->isEmcShowerValid())
    {
      RecEmcShower *emcTrk = (*itTrk[i])->emcShower();
      fEmc.E[i] = emcTrk->energy();
      fEmc.theta[i] = emcTrk->theta();
      fEmc.phi[i] = emcTrk->phi();
      fEmc.time[i] = emcTrk->time();
      fMdc.E[i] = fEmc.E[i];
    }
    RecMdcTrack  *mdcTrk = (*itTrk[i])->mdcTrack();
    fMdc.trackId[i] = mdcTrk->trackId();
    fMdc.q[i] = mdcTrk->charge(); 
    fMdc.p[i] = mdcTrk->p();
    fMdc.px[i]= mdcTrk->px();
    fMdc.py[i]= mdcTrk->py();
    fMdc.pz[i]= mdcTrk->pz();
    fMdc.theta[i]= mdcTrk->theta();
    fMdc.phi[i] = mdcTrk->phi();
    fMdc.x[i]  = mdcTrk->x();
    fMdc.y[i]  = mdcTrk->y();
    fMdc.z[i]  = mdcTrk->z();
    fMdc.x[i]  = mdcTrk->x();
    fMdc.y[i]  = mdcTrk->y();
    fMdc.z[i]  = mdcTrk->z();
    double rvxy,rvz,rvphi;
    calculate_vertex(mdcTrk,rvxy,rvz,rvphi); 
    fMdc.vxy[i] = rvxy;
    fMdc.vz[i]  = rvz; 
    fMdc.vphi[i] = rvphi; 

    //McParticleVector particles = navigator->getMcParticles((*itTrk[i])->mdcKalTrack());
    //if(!particles.empty())
    //{
    //  cout <<"REC: " << fEvent.channel << " MC: ";
    //  for(int k = 0;k<particles.size();k++)
    //  {
    //     cout << particles[k]->particleProperty() << " "; 
    //  }
    //  cout << endl;
    //  //cout << "Retrieved " << particles.size() << " McParticles for for MdcKalTrack # " 
    //  //  << mdcTrk->trackId() << " of reconstructed momentum " << mdcTrk->p() << " GeV/c (PID=" 
    //  //  << endl;
    //}
    //  //<< particles.front()->particleProperty() << endl;
    

    

    //dedx information
    if((*itTrk[i])->isMdcDedxValid())
    {
      RecMdcDedx* dedxTrk = (*itTrk[i])->mdcDedx();
      fDedx.chie[i] = dedxTrk->chiE();
      fDedx.chimu[i] = dedxTrk->chiMu();
      fDedx.chipi[i] = dedxTrk->chiPi();
      fDedx.chik[i] = dedxTrk->chiK();
      fDedx.chip[i] = dedxTrk->chiP();
      //fDedx.ghit[i] = dedxTrk->numGoodHits();
      //fDedx.thit[i] = dedxTrk->numTotalHits();
      fDedx.probPH[i] = dedxTrk->probPH();
      fDedx.normPH[i] = dedxTrk->normPH();
      //fDedx.e[i] = dedxTrk->getDedxExpect(0);
      //fDedx.mu[i] = dedxTrk->getDedxExpect(1);
      //fDedx.pi[i] = dedxTrk->getDedxExpect(2);
      //fDedx.K[i] = dedxTrk->getDedxExpect(3);
      //fDedx.p[i] = dedxTrk->getDedxExpect(4);
      //fDedx.pid[i]=dedxTrk->particleId();
    }
    if((*itTrk[i])->isTofTrackValid())
    {
      SmartRefVector<RecTofTrack> tofTrkCol = (*itTrk[i])->tofTrack();
      SmartRefVector<RecTofTrack>::iterator tofTrk = tofTrkCol.begin();
      TofHitStatus *hitst = new TofHitStatus;
      std::vector<int> tofecount;
      int goodtofetrk=0;
      for(tofTrk = tofTrkCol.begin(); tofTrk!=tofTrkCol.end(); tofTrk++,goodtofetrk++)
      {
        unsigned int st = (*tofTrk)->status();
        hitst->setStatus(st);
        //if(  (hitst->is_barrel()) ) continue;
        if( !(hitst->is_counter()) ) continue;
        tofecount.push_back(goodtofetrk);
      }
      delete hitst;
      if(!tofecount.empty())
      {
        tofTrk = tofTrkCol.begin()+tofecount[0];
        fTof.tofID[i] = (*tofTrk)->tofID();
        fTof.t0[i] = (*tofTrk)->t0();
        fTof.t[i] = (*tofTrk)->tof();
        fTof.dt[i] = (*tofTrk)->errtof();
        fTof.beta[i] = (*tofTrk)->beta();
        fTof.te[i] = (*tofTrk)->texpElectron();
        fTof.tmu[i]= (*tofTrk)->texpMuon();
        fTof.tpi[i]= (*tofTrk)->texpPion();
        fTof.tk[i] = (*tofTrk)->texpKaon();
        fTof.tp[i] = (*tofTrk)->texpProton();
        if(fTof.dt[i]>0)
        {
          fTof.chie[i]  = (fTof.t[i] - fTof.te[i])  /  fTof.dt[i];
          fTof.chimu[i] = (fTof.t[i] - fTof.tmu[i]) /  fTof.dt[i];
          fTof.chipi[i] = (fTof.t[i] - fTof.tpi[i]) /  fTof.dt[i];
          fTof.chik[i]  = (fTof.t[i] - fTof.tk[i])  /  fTof.dt[i];
          fTof.chip[i]  = (fTof.t[i] - fTof.tp[i])  /  fTof.dt[i];
        }
      }
    }
    if(fEvent.run<0 && false)
    {
      int m_numParticle(0), m_true_pid(0);
      if(!mcParticleCol)
      {
        log << MSG::ERROR << "Could not retrieve McParticelCol" << endreq;
        return StatusCode::FAILURE;
      }
      else
      {
        bool psipDecay(false);
        int rootIndex(-1);
        Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
        for (; iter_mc != mcParticleCol->end(); iter_mc++)
        {
          if ((*iter_mc)->primaryParticle()) continue;
          if (!(*iter_mc)->decayFromGenerator()) continue;
          //if ( ((*iter_mc)->mother()).trackIndex()<3 ) continue;
          if ((*iter_mc)->particleProperty()==100443)
          {
            psipDecay = true;
            rootIndex = (*iter_mc)->trackIndex();
          }
          if (!psipDecay) continue;
          int mcidx = ((*iter_mc)->mother()).trackIndex() - rootIndex;
          int pdgid = (*iter_mc)->particleProperty();
          cout << pdgid << " ";
          //m_pdgid[m_numParticle] = pdgid;
          //m_motheridx[m_numParticle] = mcidx;
          //m_numParticle ++;    

          //if(!(*iter_mc)->leafParticle()) continue;
          //if((*iter_mc)->particleProperty() == 211) m_true_pionp = (*iter_mc)->initialFourMomentum().vect().mag();
          //if((*iter_mc)->particleProperty() == -211) m_true_pionm = (*iter_mc)->initialFourMomentum().vect().mag();
        }
        //m_idxmc = m_numParticle;
      }
      cout << endl;
    }

    //fill particle id

    PID->setRecTrack(*itTrk[i]);
    PID->calculate();
    if(PID->IsPidInfoValid())
    {
      fPid.prob[ID_ELECTRON][i] = PID->probElectron();
      fPid.prob[ID_MUON][i]     = PID->probMuon();
      fPid.prob[ID_PION][i]     = PID->probPion();
      fPid.prob[ID_KAON][i]     = PID->probKaon();
      fPid.prob[ID_PROTON][i]   = PID->probProton();
    }
    vector<double> chi2 = get_chi2(itTrk[i]);
    for(int pid=0;pid<5;pid++)
    {
      fPid.chi2[pid][i]   = chi2[pid];
    }

  }
  for(int i=0;i<5;i++)
  {
    fPid.M[i]    = sqrt(get_invariant_mass2(result_pair,XMASS[i]));
  }



  fNeutral.ntrack=std::min(good_neutral_tracks.size(), size_t(RootEmc::ARRAY_SIZE));
  int idx=0;
  for(list<EvtRecTrackIterator>::iterator track=good_neutral_tracks.begin(); track!=good_neutral_tracks.end(); track++)
  {
    EvtRecTrackIterator  itTrk = *track;
    RecEmcShower *emcTrk = (*itTrk)->emcShower();
    fNeutral.E[idx]  =  emcTrk->energy();
    fNeutral.theta[idx] =  emcTrk->theta();
    fNeutral.phi[idx] =  emcTrk->phi();
    fNeutral.time[idx] = emcTrk->time();
    idx++;
  }

  fEvent.tuple->write();
  fPid.tuple->write();
  fMdc.tuple->write();
  fEmc.tuple->write();
  fDedx.tuple->write();
  fTof.tuple->write();
  fNeutral.tuple->write();
  event_write++;
  return StatusCode::SUCCESS;
}

StatusCode JpsiKK::finalize()
{
  std::cout << "Event proceed: " << event_proceed << std::endl;
  std::cout << "Event selected: " << event_write << std::endl;
  std::cout << "Event with kaons: " << event_with_kaons << std::endl;
  std::cout << "Event with muons: " << event_with_muons << std::endl;
  return StatusCode::SUCCESS;
}
