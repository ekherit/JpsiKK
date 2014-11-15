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
  status = init_tuple(this, fEvent,"FILE1/event","Signal events pi+pi- K+K-, or pi+pi- mu+mu-",log);
  status = init_tuple(this, fDedx,"FILE1/dedx","Dedx info for signal",log);
  status = init_tuple(this, fEmc,"FILE1/emc","Emc info for signal",log);
  status = init_tuple(this, fTof,"FILE1/tof","Tof info for signal",log);
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
  status = tuple->addItem ("Mrec", Mrecoil); 
  status = tuple->addItem ("Minv", Minv); 
  status = tuple->addItem ("M2mis", M2missing); 
  status = tuple->addItem ("npid", npid,0,5); 
  status = tuple->addIndexedItem ("M", npid, M); 
  status = tuple->addIndexedItem ("prob", npid, prob); 

  status = tuple->addItem ("ntrack", ntrack,0,4); //array size must be = 4
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

  //particle id information
  status = tuple->addIndexedItem ("probe",  ntrack, probe);
  status = tuple->addIndexedItem ("probmu",  ntrack, probmu);
  status = tuple->addIndexedItem ("probpi",  ntrack, probpi);
  status = tuple->addIndexedItem ("probk",  ntrack, probk);
  status = tuple->addIndexedItem ("probp",  ntrack, probp);
  return status;
}

void JpsiKK::RootEvent::init(void)
{
}


StatusCode JpsiKK::RootEmc::init_tuple(void)
{
  StatusCode status;
  status = tuple->addItem ("ntrack", ntrack,0,100); //good nuetral track in event
  status = tuple->addIndexedItem ("E",     ntrack, E);
  status = tuple->addIndexedItem ("theta", ntrack, theta);
  status = tuple->addIndexedItem ("phi",   ntrack, phi);
  status = tuple->addIndexedItem ("time",     ntrack, time);
  return status;
}

void JpsiKK::RootEmc::init(void)
{
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
}
StatusCode JpsiKK::RootTof::init_tuple(void)
{
  StatusCode status;
  status = tuple->addItem ("ntrack", ntrack,0,4); 
  status = tuple->addIndexedItem ("ID",  ntrack, tofID);
  status = tuple->addIndexedItem ("t",  ntrack, tof);
  status = tuple->addIndexedItem ("dt",  ntrack, errtof);
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

  //  Get information about reconstructed events
  SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
  SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol);
  SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(),  EventModel::MC::McParticleCol);

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
  //SELECTION CODE
  if(pion_pairs.empty()) return StatusCode::SUCCESS;
  //log << MSG::ERROR << "good charged tracks: " << charged_tracks.size() << " (" << negative_charged_tracks.size() << ", " << positive_charged_tracks.size() << endmsg;
  //log << MSG::ERROR << "pions: " << negative_pion_tracks.size()  << ", " << positive_pion_tracks.size() << endmsg;
  //log << MSG::ERROR << "other: " << other_negative_tracks.size()  << ", " << other_positive_tracks.size() << endmsg;
  //log << MSG::ERROR << "pion pairs: " << pion_pairs.size() << endmsg;

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
      EvtRecTrackIterator  itTrk[2] = {pair.first, pair.second};
      double Ep[2];
      for(int k=0;k<2;k++)
      {
        if(!(*itTrk[k])->isMdcTrackValid() || ! (*itTrk[k])->isEmcShowerValid()) 
        {
          log << MSG::ERROR << "Invalid mdc info for track.Exiting" << endmsg;
          return StatusCode::FAILURE;
        }
        //SELECTION CODE: no EMC information
        //if(! (*itTrk[k])->isEmcShowerValid()) goto SKIP_THIS_PAIR;
        RecMdcTrack *mdcTrk = (*itTrk[k])->mdcTrack();
        RecEmcShower *emcTrk = (*itTrk[k])->emcShower();
        double E = emcTrk->energy();
        double p = mdcTrk->p();
        Ep[k] = E/p;
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
        //SELECTION CODE
        if(Ep[0] < MAX_KAON_EP_RATIO && Ep[1] < MAX_KAON_EP_RATIO)
        {
          kaon_pairs.push_back(pair);
        }
      }
      //SELECTION CODE
      if(MIN_INVARIANT_MASS <  M[1]   && M[1]  < MAX_INVARIANT_MASS)
      {
        //SELECTION CODE
        if(Ep[0] < MAX_MUON_EP_RATIO && Ep[1] < MAX_MUON_EP_RATIO)
        {
          muon_pairs.push_back(pair);
        }
      }
    }

  int channel=-1; //default no channel
  std::pair<EvtRecTrackIterator,EvtRecTrackIterator> result_pair;
  //log << MSG::ERROR << "kaon pairs: " << kaon_pairs.size() << ",  muon pairs: " << muon_pairs.size() << endmsg;
  //SELECTION CODE
  if(!muon_pairs.empty() || !kaon_pairs.empty()) 
  {
    //the best pair which is closer to JPSI
    if(!kaon_pairs.empty()) result_pair = kaon_pairs.front();
    if(!muon_pairs.empty()) result_pair = muon_pairs.front();
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
      log << MSG::ERROR << "Must be some channel but it's not" << endmsg;
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
  }
  else
  {
    return StatusCode::SUCCESS;
  }


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

  fEvent.Mrecoil = get_recoil__mass(pion_pair, PION_MASS);
  fEvent.Minv    = sqrt(get_invariant_mass2(result_pair,XMASS[channel]));
  fEvent.M2missing = get_missing_mass(pion_pair,result_pair);

  ParticleID * PID = ParticleID::instance();
  PID->init();
  PID->setMethod(PID->methodProbability());
  PID->setChiMinCut(4);
  PID->usePidSys(PID->useDedx());
  PID->identify(PID->all()); 

  fEvent.ntrack=4;
  fDedx.ntrack=4;
  fEmc.ntrack=4;
  fTof.ntrack=4;
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
      fEvent.E[i] = emcTrk->energy();
      fEmc.E[i] = emcTrk->energy();
      fEmc.theta[i] = emcTrk->theta();
      fEmc.phi[i] = emcTrk->phi();
      fEmc.time[i] = emcTrk->time();
    }
    else
    {
      fEvent.E[i] = 0;
      fEmc.E[i] = 0;
      fEmc.theta[i] = -1000;;
      fEmc.phi[i] = -1000;
      fEmc.time[i] = -1000;
    }
    RecMdcTrack  *mdcTrk = (*itTrk[i])->mdcTrack();
    fEvent.index[i] = itTrk[i]-evtRecTrkCol->begin(); 
    fEvent.q[i] = mdcTrk->charge(); 
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

    PID->setRecTrack(*itTrk[i]);
    PID->calculate();
    if(PID->IsPidInfoValid())
    {
      fEvent.probe[i] =  PID->probElectron();
      fEvent.probmu[i] = PID->probMuon();
      fEvent.probpi[i] = PID->probPion();
      fEvent.probk[i] =  PID->probKaon();
      fEvent.probp[i] =  PID->probProton();
    }
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
    else
    {
      fDedx.chie[i] = -1000;
      fDedx.chimu[i] = -1000;
      fDedx.chipi[i] = -1000;
      fDedx.chik[i] = -1000; 
      fDedx.chip[i] = -1000; 
      fDedx.probPH[i] = -1000;
      fDedx.normPH[i] = -1000; 
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
        fTof.tof[i] = (*tofTrk)->tof();
        fTof.errtof[i] = (*tofTrk)->errtof();
        fTof.beta[i] = (*tofTrk)->beta();
        fTof.te[i] = (*tofTrk)->texpElectron();
        fTof.tmu[i]= (*tofTrk)->texpMuon();
        fTof.tpi[i]= (*tofTrk)->texpPion();
        fTof.tk[i] = (*tofTrk)->texpKaon();
        fTof.tp[i] = (*tofTrk)->texpProton();

        if(fTof.errtof[i]>0)
        {
          fTof.chie[i] = (fTof.tof[i]-fTof.te[i])/fTof.errtof[i];
          fTof.chimu[i]= (fTof.tof[i] -fTof.tmu[i])/fTof.errtof[i];
          fTof.chipi[i]= (fTof.tof[i] -fTof.tpi[i])/fTof.errtof[i];
          fTof.chik[i] = (fTof.tof[i] -fTof.tk[i])/fTof.errtof[i];
          fTof.chip[i] = (fTof.tof[i] -fTof.tp[i])/fTof.errtof[i];
        }
        else
        {
          fTof.chie[i] =0; 
          fTof.chimu[i]=0; 
          fTof.chipi[i]=0; 
          fTof.chik[i] =0; 
          fTof.chip[i] =0; 
        }
        cout << i << " " << fTof.chik[i] << " " << fTof.chimu[i] << endl;
      }
    }
  }
  fEvent.npid=5;
  for(int i=0;i<5;i++)
  {
    fEvent.M[i] = sqrt(get_invariant_mass2(result_pair,XMASS[i]));
  }
  fEvent.prob[0] =  fEvent.probk[2]*fEvent.probk[3];
  fEvent.prob[1] =  fEvent.probmu[2]*fEvent.probmu[3];
  fEvent.prob[2] =  fEvent.probe[2]*fEvent.probe[3];
  fEvent.prob[3] =  fEvent.probpi[2]*fEvent.probpi[3];
  fEvent.prob[4] =  fEvent.probp[2]*fEvent.probp[3];

  fNeutral.ntrack=good_neutral_tracks.size();
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
