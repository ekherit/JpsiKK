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
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"
#include "ParticleID/ParticleID.h"


const double PI_MESON_MASS =0.13957018; //GeV
const double PION_MASS =0.13957018; //GeV
const double MUON_MASS = 0.105658389; //GeV
const double KAON_MASS = 0.493677; //GeV
const double ELECTRON_MASS = 0.000510999;//GeV
const double PROTON_MASS = 0.93827231;//GeV

const double JPSI_MASS = 3.096916; //GeV
const double PSIP_MASS = 3.686093; //GeV

inline double sq(double x) { return x*x; }

double XMASS[5] = {KAON_MASS, MUON_MASS, ELECTRON_MASS, PION_MASS, PROTON_MASS};

JpsiKK::JpsiKK(const std::string& name, ISvcLocator* pSvcLocator) :
  Algorithm(name, pSvcLocator)
{
  declareProperty("MAX_CHARGED_TRACKS", MAX_CHARGED_TRACKS=4); //maximum number of charged tracks in selection
  declareProperty("MIN_CHARGED_TRACKS", MIN_CHARGED_TRACKS=4); //minimum number of charged tracks in selection
  declareProperty("MAX_NEUTRAL_TRACKS", MAX_NEUTRAL_TRACKS=0); //maximum number of good charged tracks in selection

  //good charged track configuration
  declareProperty("IP_MAX_Z", IP_MAX_Z = 10.0); //cm?
  declareProperty("IP_MAX_RHO", IP_MAX_RHO = 1.0); //cm?
  declareProperty("MAX_COS_THETA", MAX_COS_THETA = 0.93);

  //good neutral track configuration
  //endcup calorimeter
  declareProperty("EMC_ENDCUP_MIN_COS_THETA", EMC_ENDCUP_MIN_COS_THETA = 0.86);
  declareProperty("EMC_ENDCUP_MAX_COS_THETA", EMC_ENDCUP_MAX_COS_THETA = 0.92);
  declareProperty("EMC_ENDCUP_MIN_ENERGY", EMC_ENDCUP_MIN_ENERGY = 0.05);
  //barrel calorimeter
  declareProperty("EMC_BARREL_MAX_COS_THETA", EMC_BARREL_MAX_COS_THETA = 0.8);
  declareProperty("EMC_BARREL_MIN_ENERGY", EMC_BARREL_MIN_ENERGY = 0.025);

  //not electron
  declareProperty("MAX_EP_RATIO", MAX_EP_RATIO = 0.26);

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


StatusCode JpsiKK::initialize(void)
{
  MsgStream log(msgSvc(), name());

  log << MSG::INFO << "in initialize()" << endmsg;
  event_proceed=0;
  event_write = 0;
  event_with_kaons=0;
  event_with_muons=0;

  StatusCode status;

  NTuplePtr nt_event(ntupleSvc(), "FILE1/event");
  if(nt_event) fEvent.tuple = nt_event;
  else
  {
    fEvent.tuple = ntupleSvc()->book("FILE1/event", CLID_ColumnWiseTuple, "Signal events pi+pi- K+K-, or pi+pi- mu+mu-");
    log << MSG::INFO << "After book" << endmsg;
    if(fEvent.tuple)
    {
      status = fEvent.init_tuple();
    }
    else
    {
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(fEvent.tuple) << endmsg;
      return StatusCode::FAILURE;
    }
  }
  return StatusCode::SUCCESS;
}


StatusCode JpsiKK::RootEvent::init_tuple(void)
{
  StatusCode status;
  status = tuple->addItem ("run", run); //run number
  status = tuple->addItem ("event", event); //event number
  status = tuple->addItem ("time", time); //event time
  status = tuple->addItem ("ngoodtrack", ngood_track); //good charged track in event
  status = tuple->addItem ("nptrack", npositive_track); //good positive charged track in event
  status = tuple->addItem ("nntrack", nnegative_track); //good negative charged track in event
  status = tuple->addItem ("nppions", npositive_pions); //good poitive pion tracks in event
  status = tuple->addItem ("nnpions", nnegative_pions); //good negative pion track in event
  status = tuple->addItem ("npion_pairs", npion_pairs); //number of pions paris in event
  status = tuple->addItem ("channel", channel); //decay channel of the J/psi
  status = tuple->addItem ("Mrec", Mrecoil); 
  status = tuple->addItem ("Minv", Minv); 
  status = tuple->addItem ("M2mis", M2missing); 
  status = tuple->addItem ("npid", npid,0,5); 
  status = tuple->addIndexedItem ("M", npid, M); 

  status = tuple->addItem ("ntrack", ntrack,0,2); //array size must be =2
  //pions information
  status = tuple->addIndexedItem ("idx",   ntrack, index);
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


double get_invariant_mass2(std::pair<EvtRecTrackIterator,EvtRecTrackIterator> & pair, double mass)
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
    std::cout << "proceed event: " << setw(15) << event_proceed;
    std::cout << std::endl;
  }
  event_proceed++;

  /*  Get information about reconstructed events */
  SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
  SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol);
  SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(),  EventModel::MC::McParticleCol);

  fEvent.init();

  std::list<EvtRecTrackIterator> good_charged_tracks;
  std::list<EvtRecTrackIterator> good_neutral_tracks;

  //now count good charged track
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
  //now count good neutral track
  for(int i = evtRecEvent->totalCharged(); i<evtRecEvent->totalTracks(); i++)
  {
    EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
    if(!(*itTrk)->isEmcShowerValid()) continue; //keep only valid neutral tracks
    RecEmcShower *emcTrk = (*itTrk)->emcShower();
    double c =  fabs(cos(emcTrk->theta())); //abs cos theta
    double E  =  emcTrk->energy();
    bool barrel = (c <= EMC_BARREL_MAX_COS_THETA);
    bool endcup = (EMC_ENDCUP_MIN_COS_THETA <=c) && (c <= EMC_ENDCUP_MAX_COS_THETA);
    bool good_track = ( barrel && E > EMC_BARREL_MIN_ENERGY) || (endcup && EMC_ENDCUP_MIN_ENERGY);
    if(good_track) good_neutral_tracks.push_back(itTrk);
  }

  //print good charged track index
  //cout << "Good charged track: ";
  //for(list<EvtRecTrackIterator>::iterator i=good_charged_tracks.begin();i!=good_charged_tracks.end();i++)
  //{
  //  cout << *i - evtRecTrkCol->begin() << " ";
  //}
  //cout << endl;

  log << MSG::INFO << "MAX_NEUTRAL_TRACKS  = " << MAX_NEUTRAL_TRACKS << endmsg;
  log << MSG::INFO << "MIN_CHARGED_TRACKS  = " << MIN_CHARGED_TRACKS << endmsg;
  log << MSG::INFO << "MAX_CHARGED_TRACKS  = " << MAX_CHARGED_TRACKS << endmsg;
  log << MSG::INFO << "good charged tracks: " << charged_tracks.size() <<  ",  neutral tracks: " << good_neutral_tracks.size() << endmsg;

  //SELECTION CODE
  if( MAX_NEUTRAL_TRACKS < good_neutral_tracks.size()) return StatusCode::SUCCESS;
  //SELECTION CODE
  if( good_charged_tracks.size() < MIN_CHARGED_TRACKS || MAX_CHARGED_TRACKS < good_charged_tracks.size()) return StatusCode::SUCCESS;


  
  std::list<EvtRecTrackIterator> charged_tracks; //selected tracks for specific cut
  std::list<EvtRecTrackIterator> positive_charged_tracks; //selected tracks for specific cut
  std::list<EvtRecTrackIterator> negative_charged_tracks; //selected tracks for specific cut
  std::list<EvtRecTrackIterator> positive_pion_tracks; //selected pion tracks for specific cut
  std::list<EvtRecTrackIterator> negative_pion_tracks; //selected pion tracks for specific cut
  std::list<EvtRecTrackIterator> other_positive_tracks; //other positive tracks for specific cut
  std::list<EvtRecTrackIterator> other_negative_tracks; //other positive tracks for specific cut
  for(std::list<EvtRecTrackIterator>::iterator track=good_charged_tracks.begin(); track!=good_charged_tracks.end(); track++)
  {
    EvtRecTrackIterator & itTrk = *track;
    if(!(*itTrk)->isMdcTrackValid()) continue; 
    if(!(*itTrk)->isEmcShowerValid()) continue; 
    RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
    RecEmcShower *emcTrk = (*itTrk)->emcShower();
    double c = fabs(cos(mdcTrk->theta()));
    double E = emcTrk->energy();
    double p = mdcTrk->p();
    double q = mdcTrk->charge();
    bool barrel = c < EMC_BARREL_MAX_COS_THETA;
    bool not_electron_pion = E/p < MAX_EP_RATIO; 
    if(barrel) 
    {
      if(q>0) 
      {
        positive_charged_tracks.push_back(itTrk);
        if(p<MAX_PION_MOMENTUM) 
        {
          positive_pion_tracks.push_back(itTrk);
        }
        if(p>MIN_KAON_MOMENTUM && not_electron_pion)
        {
          other_positive_tracks.push_back(itTrk);
        }
      }
      if(q<0) 
      {
        negative_charged_tracks.push_back(itTrk);
        if(p<MAX_PION_MOMENTUM) 
        {
          negative_pion_tracks.push_back(itTrk);
        }
        if(p>MIN_KAON_MOMENTUM && not_electron_pion)
        {
          other_negative_tracks.push_back(itTrk);
        }
      }
      charged_tracks.push_back(itTrk);
    }
  }

  
  log << MSG::INFO << "good charged tracks: " << charged_tracks.size() << " (" << negative_charged_tracks.size() << ", " << positive_charged_tracks.size() << endmsg;
  log << MSG::INFO << "pions: " << negative_pion_tracks.size()  << ", " << positive_pion_tracks.size() << endmsg;
  log << MSG::INFO << "other: " << other_negative_tracks.size()  << ", " << other_positive_tracks.size() << endmsg;


  //SELECTION CODE
  //keep only specific signature
  if(positive_charged_tracks.size()!=2 || negative_charged_tracks.size()!=2) return StatusCode::SUCCESS;
  //SELECTION CODE
  if(negative_pion_tracks.empty() || positive_pion_tracks.empty()) return StatusCode::SUCCESS;

  typedef std::list< std::pair<EvtRecTrackIterator, EvtRecTrackIterator> > PairList_t;
  //create pion pairs
  PairList_t pion_pairs;
  for(list<EvtRecTrackIterator>::iterator i=negative_pion_tracks.begin(); i!=negative_pion_tracks.end(); ++i)
    for(list<EvtRecTrackIterator>::iterator j=positive_pion_tracks.begin(); j!=positive_pion_tracks.end(); ++j)
    {
      std::pair<EvtRecTrackIterator,EvtRecTrackIterator> pair(*i,*j);
      double M_recoil = get_recoil__mass(pair, PION_MASS);
      if(MIN_RECOIL_MASS < M_recoil && M_recoil < MAX_RECOIL_MASS) 
      {
        pion_pairs.push_back(pair);
      }
    }
  log << MSG::INFO << "pion pairs: " << pion_pairs.size() << endmsg;
  //SELECTION CODE
  if(pion_pairs.empty()) return StatusCode::SUCCESS;

  //find the best pion pair using closest value to JPSI_MASS
  std::pair<EvtRecTrackIterator,EvtRecTrackIterator> pion_pair = pion_pairs.front();
  for(PairList_t::iterator p=pion_pairs.begin();p!=pion_pairs.end();p++)
  {
    if(fabs(get_recoil__mass(*p,PION_MASS) - JPSI_MASS) <  fabs(get_recoil__mass(pion_pair,PION_MASS) - JPSI_MASS)) pion_pair = *p;
  }

  //make kaon or muon pairs
  PairList_t muon_pairs;
  PairList_t kaon_pairs;
  for(list<EvtRecTrackIterator>::iterator i=other_negative_tracks.begin(); i!=other_negative_tracks.end(); ++i)
    for(list<EvtRecTrackIterator>::iterator j=other_positive_tracks.begin(); j!=other_positive_tracks.end(); ++j)
    {
      std::pair<EvtRecTrackIterator,EvtRecTrackIterator> pair(*i,*j);
      double M[5]={0,0,0,0,0};
      for(int pid=0;pid<5;pid++)
      {
        M[pid]=get_invariant_mass2(pair,XMASS[pid]);
        if(M[pid]>0) M[pid] = sqrt(M[pid]);
        else M[pid] = 0;
      }
      if(MIN_INVARIANT_MASS <  M[0] - JPSI_MASS  && M[0] - JPSI_MASS < MAX_INVARIANT_MASS)   kaon_pairs.push_back(pair);
      if(MIN_INVARIANT_MASS <  M[1] - JPSI_MASS  && M[1] - JPSI_MASS < MAX_INVARIANT_MASS)   muon_pairs.push_back(pair);
    }

  log << MSG::INFO << "kaon pairs: " << kaon_pairs.size() << ",  muon pairs: " << muon_pairs.size() << endmsg;
  //SELECTION CODE
  if(muon_pairs.empty() && kaon_pairs.empty()) return StatusCode::SUCCESS;


  //the best pair which is closer to JPSI
  std::pair<EvtRecTrackIterator,EvtRecTrackIterator> result_pair;
  if(!kaon_pairs.empty()) result_pair = kaon_pairs.front();
  if(!muon_pairs.empty()) result_pair = muon_pairs.front();
  int channel=-1; //default no channel
  for(PairList_t::iterator p=kaon_pairs.begin();p!=kaon_pairs.end();p++)
  {
    if(fabs(sqrt(get_invariant_mass2(*p,KAON_MASS)) - JPSI_MASS) 
        <=  fabs(sqrt(get_invariant_mass2(result_pair,KAON_MASS)) - JPSI_MASS)) 
    {
      result_pair = *p;
      channel=0; //setup kaon channel
    }
  }
  for(PairList_t::iterator p=muon_pairs.begin();p!=muon_pairs.end();p++)
  {
    if(fabs(sqrt(get_invariant_mass2(*p,MUON_MASS)) - JPSI_MASS) 
        <=  fabs(sqrt(get_invariant_mass2(result_pair,MUON_MASS)) - JPSI_MASS)) 
    {
      result_pair = *p;
      channel=1; //setup muon channel
    }
  }
  if(channel<0) 
  {
    clog << MSG::WARNING << "Must be some channel but it's not" << endmsg;
    return StatusCode::FAILURE; 
  }
  switch(channel)
  {
    case 0:
      event_with_kaons++;
      break;
    case 1:
      event_with_muons++;
      break;
  }


  //now fill the tuples

  //some statistics information
  fEvent.ngood_track = good_charged_tracks.size();
  fEvent.npositive_track = positive_charged_tracks.size();
  fEvent.nnegative_track = negative_charged_tracks.size();
  fEvent.npositive_pions = positive_pion_tracks.size();
  fEvent.nnegative_pions = negative_pion_tracks.size();
  fEvent.npion_pairs = pion_pairs.size();
  // fill the decay channel of the J/psi 0 - kaons, 1 --muons
  fEvent.channel = channel; 

  fEvent.Mrecoil = get_recoil__mass(pion_pair, PION_MASS);
  fEvent.Minv    = sqrt(get_invariant_mass2(result_pair,XMASS[channel]));
  fEvent.M2missing = get_missing_mass(pion_pair,result_pair);

  fEvent.ntrack=4;
  EvtRecTrackIterator itTrk[4] = {pion_pair.first, pion_pair.second, result_pair.first, result_pair.second};
  for(int i=0;i<4;i++)
  {
    if(!(*itTrk[i])->isMdcTrackValid() || !(*itTrk[i])->isEmcShowerValid()) 
    {
      clog << MSG::ERROR << "Somthing wrong in selection. This track must have Mdc and Emc information but it is not. exiting." << endmsg;
      return StatusCode::FAILURE;
    }
    RecMdcTrack  *mdcTrk = (*itTrk[i])->mdcTrack();
    RecEmcShower *emcTrk = (*itTrk[i])->emcShower();
    fEvent.index[i] = itTrk[i]-evtRecTrkCol->begin(); 
    fEvent.q[i] = mdcTrk->charge(); 
    fEvent.E[i] = emcTrk->energy();
    fEvent.p[i] = mdcTrk->p();
    fEvent.px[i]= mdcTrk->px();
    fEvent.py[i]= mdcTrk->py();
    fEvent.pz[i]= mdcTrk->pz();
    fEvent.theta[i]= mdcTrk->theta();
    fEvent.phi[i] = mdcTrk->phi();
    fEvent.x[i]  = mdcTrk->x();
    fEvent.y[i]  = mdcTrk->y();
    fEvent.z[i]  = mdcTrk->z();
    double rvxy,rvz,rvphi;
    calculate_vertex(mdcTrk,rvxy,rvz,rvphi); 
    fEvent.vxy[i] = rvxy;
    fEvent.vz[i]  = rvz; 
    fEvent.vphi[i] = rvphi; 
  }
  fEvent.npid=5;
  for(int i=0;i<5;i++)
  {
    fEvent.M[i] = sqrt(get_invariant_mass2(result_pair,XMASS[i]));
  }

  fEvent.tuple->write();
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
