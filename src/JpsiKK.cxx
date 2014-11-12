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
  declareProperty("MAX_CHARGED_TRACKS", MAX_CHARGED_TRACKS=10); //maximum number of charged tracks in selection
  declareProperty("MIN_CHARGED_TRACKS", MIN_CHARGED_TRACKS=3); //minimum number of charged tracks in selection

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
  declareProperty("MIN_KAON_MOMENTUM", MAX_KAON_MOMENTUM = 1.0); //GeV
  declareProperty("MAX_KAON_MOMENTUM", MAX_KAON_MOMENTUM = 2.0); //GeV
  declareProperty("MIN_MUON_MOMENTUM", MAX_MUON_MOMENTUM = 1.0); //GeV
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
  good_pion_pairs_number=0;
  good_high_mom_pairs_number=0;

  StatusCode status;

  NTuplePtr nt_event(ntupleSvc(), "FILE1/event");
  if(nt_event) fEvent.tuple = nt_event;
  else
  {
    fEvent.tuple = ntupleSvc()->book("FILE1/event", CLID_ColumnWiseTuple, "Signal events pi+pi- K+K-, or pi+pi- mu+mu-");
    if(fEvent.tuple)
    {
      status = fEvent.init_tuple(muc_tuple);
    }
    else
    {
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(fEvent.tuple) << endmsg;
      return StatusCode::FAILURE;
    }
  }
  return StatusCode::SUCCESS;
}


StatusCode JpsiKK::RootEvent::init_tuple(NTuple::Tuple * tuple)
{
  StatusCode status;
  status = tuple->addItem ("ntrack", ntrack); //good charged track in event
  status = tuple->addItem ("nptrack", npositive_track); //good positive charged track in event
  status = tuple->addItem ("nntrack", nnegative_track); //good negative charged track in event
  status = tuple->addItem ("nppions", npositive_pions); //good poitive pion tracks in event
  status = tuple->addItem ("nnpions", nnegative_pions); //good negative pion track in event
  status = tuple->addItem ("npion_pairs", npion_pairs); //number of pions paris in event
  status = tuple->addItem ("channel", channel); //decay channel of the J/psi
  status = tuple->addItem ("ngood_pions", ngood_pions); //decay channel of the J/psi
  //pions information
  status = tuple->addIndexedItem ("pidx", ngood_pions, pions.index);
  status = tuple->addIndexedItem ("pq", ngood_pions, pions.q);
  status = tuple->addIndexedItem ("pE", ngood_pions, pions.E);
  status = tuple->addIndexedItem ("pp", ngood_pions, pions.p);
  status = tuple->addIndexedItem ("ppx", ngood_pions, pions.px);
  status = tuple->addIndexedItem ("ppy", ngood_pions, pions.py);
  status = tuple->addIndexedItem ("ppz", ngood_pions, pions.pz);
  status = tuple->addIndexedItem ("ppt", ngood_pions, pions.pt);
  status = tuple->addIndexedItem ("ptheta", ngood_pions, pions.theta);
  status = tuple->addIndexedItem ("pphi", ngood_pions, pions.phi);
  status = tuple->addIndexedItem ("px", ngood_pions, pions.x);
  status = tuple->addIndexedItem ("py", ngood_pions, pions.y);
  status = tuple->addIndexedItem ("pz", ngood_pions, pions.z);
  status = tuple->addIndexedItem ("pr", ngood_pions, pions.r);
  status = tuple->addIndexedItem ("pvxy", ngood_pions, pions.vxy);
  status = tuple->addIndexedItem ("pvz", ngood_pions, pions.vz);
  status = tuple->addIndexedItem ("pvphi", ngood_pions, pions.vphi);
  //kaons or muons information
  status = tuple->addIndexedItem ("kidx", ngood_pions, kmuons.index);
  status = tuple->addIndexedItem ("kq", ngood_pions, kmuons.q);
  status = tuple->addIndexedItem ("kE", ngood_pions, kmuons.E);
  status = tuple->addIndexedItem ("kp", ngood_pions, kmuons.p);
  status = tuple->addIndexedItem ("kpx", ngood_pions, kmuons.px);
  status = tuple->addIndexedItem ("kpy", ngood_pions, kmuons.py);
  status = tuple->addIndexedItem ("kpz", ngood_pions, kmuons.pz);
  status = tuple->addIndexedItem ("kpt", ngood_pions, kmuons.pt);
  status = tuple->addIndexedItem ("ktheta", ngood_pions, kmuons.theta);
  status = tuple->addIndexedItem ("kphi", ngood_pions, kmuons.phi);
  status = tuple->addIndexedItem ("kx", ngood_pions, kmuons.x);
  status = tuple->addIndexedItem ("ky", ngood_pions, kmuons.y);
  status = tuple->addIndexedItem ("kz", ngood_pions, kmuons.z);
  status = tuple->addIndexedItem ("kr", ngood_pions, kmuons.r);
  status = tuple->addIndexedItem ("kvxy", ngood_pions, kmuons.vxy);
  status = tuple->addIndexedItem ("kvz", ngood_pions, kmuons.vz);
  status = tuple->addIndexedItem ("kvphi", ngood_pions, kmuons.vphi);
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

StatusCode JpsiKK::execute()
{
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in execute()" << endreq;

  SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
  int runNo=eventHeader->runNumber();
  //if(runNo<0) CHECK_MC = true;
  int event=eventHeader->eventNumber();
  head_event_number=event;
  head_run=runNo;
  time_t t=eventHeader->time();
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


  InitData(evtRecEvent->totalCharged(), evtRecEvent->totalNeutral());

  //typedef std::multimap <double, unsigned> mmap_t;
  //typedef std::pair <double, unsigned> pair_t;
  //mmap_t pmap;

  std::list<EvtRecTrackIterator> good_charged_tracks;
  std::list<EvtRecTrackIterator> good_neutral_tracks;

  //now count good charged track
  for(unsigned i = 0; i < evtRecEvent->totalCharged(); ++i)
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
  for(int i = evtRecEvent->totalCharged(); i<evtRecEvent->totalTracks(); ++i)
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


  //number of good neutral tracks must be 0
  if(!good_neutral_tracks.empty()) return StatusCode::SUCCESS;

  
  std::list<EvtRecTrackIterator> charged_tracks; //selected tracks for specific cut
  std::list<EvtRecTrackIterator> positive_charged_tracks; //selected tracks for specific cut
  std::list<EvtRecTrackIterator> negative_charged_tracks; //selected tracks for specific cut
  std::list<EvtRecTrackIterator> positive_pion_tracks; //selected pion tracks for specific cut
  std::list<EvtRecTrackIterator> negative_pion_tracks; //selected pion tracks for specific cut
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
    bool not_electron = E/p < MAX_EP_RATIO; 
    if(barrel & not_electron) 
    {
      if(q>0) 
      {
        positive_charged_tracks.push_back(itTrk);
        if(p<MAX_PION_MOMENTUM) 
        {
          positive_pion_tracks.push_back(itTrk);
        }
      }
      if(q<0) 
      {
        if(p<MAX_PION_MOMENTUM) 
        {
          negative_pion_tracks.push_back(itTrk);
        }
        negative_charged_tracks.push_back(itTrk);
      }
      charged_tracks.push_back(itTrk);
    }
  }

  //keep only specific signature
  if(positive_charged_tracks.size()!=2 || negative_charged_tracks.size()!=2) return StatusCode::SUCCESS;

  if(negative_pion_tracks.empty() || positive_pion_tracks.empty()) return StatusCode::SUCCESS;
  std::list< std::pair<EvtRecTrackIterator, EvtRecTrackIterator> > pion_pairs;
  //create pion pairs
  for(list<EvtRecTrackIterator>::iterator i=negative_pion_tracks.begin(); i!=positive_pion_tracks.end(); ++i)
    for(list<EvtRecTrackIterator>::iterator j=negative_pion_tracks.begin(); j!=positive_pion_tracks.end(); ++j)
    {
      std::pair<EvtRecTrackIterator,EvtRecTrackIterator> pair(*i,*j);
      double M_recoil = get_recoil__mass(pair, PION_MASS);
      if(MIN_RECOIL_MASS < M_recoil && M_recoil < MAX_RECOIL_MASS) 
      {
        pion_pairs.push_back(pair);
      }
    }

  if(pion_pairs.empty()) return StatusCode::SUCCESS;


  //now fill the pion information
  fEvent.ntrack = good_charged_tracks.size();
  fEvent.npositive_track = positive_charged_tracks.size();
  fEvent.nnegative_track = negative_charged_tracks.size();
  fEvent.npositive_pions = positive_pion_tracks.size();
  fEvent.nnegative_track = negative_pion_tracks.size();
  fEvent.npion_pairs = pion_pairs.size();
  fEvent.ngood_pions = 2;
  fEvent.channel = -1; //yet not identify other particles

  fEvent.Mrec = get_recoil__mass(pion_pairs.front(), PION_MASS);

  fEvent.tuple->write();
  event_write++;
  return StatusCode::SUCCESS;
}

StatusCode JpsiKK::finalize()
{
  std::cout << "Event proceed: " << event_proceed << std::endl;
  std::cout << "Event selected: " << event_write << std::endl;
  return StatusCode::SUCCESS;
}
// for particle id look /ihepbatch/bes/alex/workarea/Analysis/Physics/PsiPrime/G2MuMuAlg-00-00-01/PipiJpsiAlg/src
