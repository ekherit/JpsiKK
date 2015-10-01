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
//const double PSIP_MASS = 3.686093; //GeV 
const double PSIP_MASS = 3.686109; //GeV PDG-2014

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
  declareProperty("CENTER_MASS_ENERGY", CENTER_MASS_ENERGY = 0); //GeV
  declareProperty("MIN_CHARGED_TRACKS", MIN_CHARGED_TRACKS=3); //minimum number of charged tracks in selection
  declareProperty("MAX_CHARGED_TRACKS", MAX_CHARGED_TRACKS=4); //maximum number of charged tracks in selection
  declareProperty("MAX_NEUTRAL_TRACKS", MAX_NEUTRAL_TRACKS=1000); //maximum number of good charged tracks in selection

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
  declareProperty("MIN_MUON_EP_RATIO", MIN_MUON_EP_RATIO = 0);

  declareProperty("MAX_KAON_EP_RATIO", MAX_KAON_EP_RATIO = 0.8);
  declareProperty("MIN_KAON_EP_RATIO", MIN_KAON_EP_RATIO = 0);

  declareProperty("MAX_PION_MOMENTUM", MAX_PION_MOMENTUM = 0.45); //GeV
  declareProperty("MIN_PION_MOMENTUM", MIN_PION_MOMENTUM = 0); //GeV

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
  event_with_pions=0;
  event_with_electrons=0;
  event_with_protons=0;
  good_kinematic_fit=0;
	if(CENTER_MASS_ENERGY == 0) CENTER_MASS_ENERGY = PSIP_MASS;

  StatusCode status;
  status = init_tuple(this, fEvent,  "FILE1/event","Signal events pi+pi- K+K-, or pi+pi- mu+mu-",log);
  status = init_tuple(this, fPid,    "FILE1/pid","particle id",log);
  status = init_tuple(this, fMdc,    "FILE1/mdc","Mdc info for signal",log);
  status = init_tuple(this, fDedx,   "FILE1/dedx","Dedx info for signal",log);
  status = init_tuple(this, fEmc,    "FILE1/emc","Emc info for signal",log);
  status = init_tuple(this, fTof,    "FILE1/tof","Tof info for signal",log);
  status = init_tuple(this, fNeutral,"FILE1/neutral","Good neutral tracks",log);
  status = init_tuple(this, fMC,     "FILE1/mc","Monte Carlo truth information",log);
  status = init_tuple(this, fMCTopo,     "FILE1/mctopo","Monte Carlo truth information topology",log);

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
  status = tuple->addItem ("sign", sign); //number of pions paris in event
  status = tuple->addItem ("channel", channel); //decay channel of the J/psi
  status = tuple->addItem ("KK", KK); //KK decay channel of the J/psi
  status = tuple->addItem ("uu", uu); //mu-mu decay channel of the J/psi
  status = tuple->addItem ("Mrec", Mrecoil); 
  status = tuple->addItem ("Minv", Minv); 
  status = tuple->addItem ("M012", M012); 
  status = tuple->addItem ("M013", M013); 
  status = tuple->addItem ("M03", M03); 
  status = tuple->addItem ("M12", M12); 
  status = tuple->addItem ("M01", M01); 
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

  status = tuple->addItem ("kMee", kM[ID_ELECTRON]);
  status = tuple->addItem ("kMKK", kM[ID_KAON]);
  status = tuple->addItem ("kMuu", kM[ID_MUON]);
  status = tuple->addItem ("kMpp", kM[ID_PROTON]);
  status = tuple->addItem ("kMpipi", kM[ID_PION]);
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

StatusCode JpsiKK::RootMC::init_tuple(void)
{
  StatusCode status;
  status = tuple->addItem ("psip_decay", psip_decay); //flag for psip decay
  status = tuple->addItem ("jpsi_decay", jpsi_decay); //flag for jpsi decay 
  status = tuple->addItem ("KK", KK);               //KK event
  status = tuple->addItem ("uu", uu);               //mu mu event
  status = tuple->addItem ("oo", oo);               //other event
  status = tuple->addItem ("ntrack", ntrack,0,4); 
  status = tuple->addIndexedItem ("id",    ntrack, pid);
  status = tuple->addIndexedItem ("q",     ntrack, q);
  status = tuple->addIndexedItem ("E",     ntrack, E);
  status = tuple->addIndexedItem ("p",     ntrack, p);
  status = tuple->addIndexedItem ("px",    ntrack, px);
  status = tuple->addIndexedItem ("py",    ntrack, py);
  status = tuple->addIndexedItem ("pz",    ntrack, pz);
  status = tuple->addIndexedItem ("pt",    ntrack, pt);
  status = tuple->addIndexedItem ("theta", ntrack, theta);
  status = tuple->addIndexedItem ("phi",   ntrack, phi);
  return status;
}

void JpsiKK::RootMC::init(void)
{
  ntrack=4;
}

StatusCode JpsiKK::RootMCTopo::init_tuple(void)
{
  StatusCode status;
  status = tuple->addItem("indexmc", m_idxmc, 0, 100);
  status = tuple->addIndexedItem("pdgid", m_idxmc, m_pdgid);
  status = tuple->addIndexedItem("motheridx", m_idxmc, m_motheridx);
  status = tuple->addIndexedItem("idx", m_idxmc, m_idx);
  status = tuple->addItem("hash", m_hash);
  return status;
}

void JpsiKK::RootMCTopo::init(void)
{
  m_idxmc=0;
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

std::list<EvtRecTrackIterator> JpsiKK::createGoodChargedTrackList(
		SmartDataPtr<EvtRecEvent>    & evtRecEvent, 
		SmartDataPtr<EvtRecTrackCol> & evtRecTrkCol
		)
{
  std::list<EvtRecTrackIterator> good_charged_tracks;
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
	return good_charged_tracks;
}

std::list<EvtRecTrackIterator> JpsiKK::createGoodNeutralTrackList(
		SmartDataPtr<EvtRecEvent>    & evtRecEvent, 
		SmartDataPtr<EvtRecTrackCol> & evtRecTrkCol
		)
{
	std::list<EvtRecTrackIterator> good_neutral_tracks;
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
	return good_neutral_tracks;
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

double get_recoil__mass(EvtRecTrackIterator & trk1, EvtRecTrackIterator & trk2, double mass,  double W /*  total energy */)
{
  EvtRecTrackIterator  itTrk[2] = {trk1, trk2};
  HepLorentzVector  P[2];
  for(int k=0;k<2;k++)
  {
    if(!(*itTrk[k])->isMdcTrackValid()) break; 
    RecMdcTrack *mdcTrk = (*itTrk[k])->mdcTrack();
    P[k] = mdcTrk->p4(mass);
  }
  HepLorentzVector P_psip(W*sin(0.011),0,0,W); //initial vector of psip
  HepLorentzVector P_sum = P[0]+P[1];
  HepLorentzVector P_recoil = P_psip - P_sum;
  return P_recoil.m();
}

double get_recoil__mass(std::pair<EvtRecTrackIterator, EvtRecTrackIterator> p, double mass,  double W)
{
  return get_recoil__mass(p.first, p.second, mass,  W);
}


double get_missing_mass(std::pair<EvtRecTrackIterator, EvtRecTrackIterator> pions, std::pair<EvtRecTrackIterator, EvtRecTrackIterator> kaons,  double W)
{
  EvtRecTrackIterator  PionTrk[2] = {pions.first, pions.second};
  EvtRecTrackIterator  KaonTrk[2] = {kaons.first, kaons.second};
  HepLorentzVector P_psip(W*sin(0.011),0,0,W); //initial vector of psip
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
    if( !hitst->is_cluster() ) continue;
    //if(  (hitst->is_barrel()) ) continue;
    //if( !(hitst->is_counter()) ) continue;
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


bool vertex_fit(const std::vector<WTrackParameter> & input_tracks,  std::vector<WTrackParameter> output_tracks)
{
  //vertex fit - уточним вершины
  //initial vertex
  HepPoint3D vx(0., 0., 0.);
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
	std::vector<int> index_list(input_tracks.size());
	for(int i=0; i<input_tracks.size(); i++)
	{
		vtxfit->AddTrack(i, input_tracks[i]);
		index_list[i] = i;
	}
  vtxfit->AddVertex(0, vxpar,index_list);
  if(!vtxfit->Fit()) return false;
  //vtxfit->Fit();
  vtxfit->Swim(0);

	output_tracks.resize(input_tracks.size());
	for(int i=0; i<output_tracks.size(); i++)
	{
		output_tracks[i] = vtxfit->wtrk(i);
	}
	return true;
}

/* Kinematic fit for specific pairs */
bool kinematic_fit(
    int PID,  // hypotezies 0 - kaon, 1 - muon
    TrackPair_t  &  pion_pair,  //pion pair
    TrackPair_t  & other_pair,  // other pair (kaon or muon)
    std::vector<HepLorentzVector> & P,  //4momentum fit result
    double & chi2,  //chi2 result of the fit, 
		double W
    )
{
  P.resize(4);
  EvtRecTrackIterator  PionTrk[2] = {pion_pair.first, pion_pair.second};
  EvtRecTrackIterator  OtherTrk[2] = {other_pair.first, other_pair.second};
  RecMdcKalTrack * PionKalTrk[2];
  RecMdcKalTrack * OtherKalTrk[2];
  WTrackParameter PionWTrk[2];
  WTrackParameter OtherWTrk[2];

  //For positive and negative charged pair create track parameters (WTrackParameter)
  //in different hyptoizies scpecified by PID
  for(int i=0;i<2;i++)
  {
    PionKalTrk[i] = (*PionTrk[i])->mdcKalTrack();
    OtherKalTrk[i] = (*OtherTrk[i])->mdcKalTrack();
    //pions have to be only pions
    PionWTrk[i] = WTrackParameter(PION_MASS, PionKalTrk[i]->getZHelix(), PionKalTrk[i]->getZError());
    //create other
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
  //vertex fit - уточним вершины
  //initial vertex
  HepPoint3D vx(0., 0., 0.);
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
  if(!vtxfit->Fit(0)) return false;
  vtxfit->Fit();
  vtxfit->Swim(0);

  KalmanKinematicFit * kmfit = KalmanKinematicFit::instance();
  //kmfit->setIterNumber(10000);
  //kmfit->setChisqCut(10000);

  kmfit->init();
  for(int i=0;i<4;i++)
  {
    kmfit->AddTrack(i,vtxfit->wtrk(i));
  }
  HepLorentzVector Pcmf(W*sin(0.011) /* 40.546 MeV*/,0,0,W); //initial vector of center of mass frame
  kmfit->AddFourMomentum(0,  Pcmf);
  //kmfit->AddTotalEnergy(0,PSIP_MASS,0,1,2,3);
  //kmfit->AddResonance(1, JPSI_MASS, 2, 3);
  //kmfit->AddResonance(1, PSIP_MASS, 0, 1, 2,3);
  if(!kmfit->Fit(0)) return false;
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



bool kinematic_fit(
    int PID,  // hypotezies 0 - kaon, 1 - muon,  ...
		const std::vector<EvtRecTrackIterator> & Tracks, 
		/* Tracks[0] - pi-
		 * Tracks[1] - pi+
		 * Tracks[2] - K/mu-
		 * Tracks[3] - K/mu+
		 */
    std::vector<HepLorentzVector> & P,     //4momentum fit result
    double & chi2,  //chi2 result of the fit, 
		double W
    )
{
	std::vector<RecMdcKalTrack*> KalTrk(Tracks.size());
	std::vector<WTrackParameter> WTrk(Tracks.size());
	for(int i=0;i<Tracks.size();i++)
	{
		KalTrk[i] = (*Tracks[i])->mdcKalTrack();
		if(i<2) WTrk[i] = WTrackParameter(PION_MASS, KalTrk[i]->getZHelix(), KalTrk[i]->getZError());
		else
		{
			switch(PID)
			{
				case ID_KAON:
					WTrk[i] = WTrackParameter(XMASS[PID], KalTrk[i]->getZHelixK(),  KalTrk[i]->getZErrorK());
					break;
				case ID_MUON:
					WTrk[i] = WTrackParameter(XMASS[PID], KalTrk[i]->getZHelixMu(), KalTrk[i]->getZErrorMu());
					break;
				case ID_ELECTRON:
					WTrk[i] = WTrackParameter(XMASS[PID], KalTrk[i]->getZHelixE(),  KalTrk[i]->getZErrorE());
					break;
				case ID_PION:
					WTrk[i] = WTrackParameter(XMASS[PID], KalTrk[i]->getZHelix(),   KalTrk[i]->getZError());
					break;
				case ID_PROTON:
					WTrk[i] = WTrackParameter(XMASS[PID], KalTrk[i]->getZHelixP(),  KalTrk[i]->getZErrorP());
					break;
			}
		}
	}

  //vertex fit - уточним вершины
  //initial vertex
  HepPoint3D vx(0., 0., 0.);
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
	for(int i=0;i<WTrk.size();i++)
	{
		vtxfit->AddTrack(i, WTrk[i]);
	}
	if(WTrk.size()<4) vtxfit->AddMissTrack(3,XMASS[PID]);
  //add tracks. I know the first two one must be pions
  vtxfit->AddVertex(0, vxpar,0, 1, 2, 3);
  if(!vtxfit->Fit(0)) return false;
  vtxfit->Fit();
  vtxfit->Swim(0);

  KalmanKinematicFit * kmfit = KalmanKinematicFit::instance();
  //kmfit->setIterNumber(10000);
  //kmfit->setChisqCut(10000);

  kmfit->init();
  for(int i=0;i<4;i++)
  {
    kmfit->AddTrack(i,vtxfit->wtrk(i));
  }
  HepLorentzVector Pcmf(W*sin(0.011) /* 40.546 MeV*/,0,0,W); //initial vector of center of mass frame
  kmfit->AddFourMomentum(0,  Pcmf);
  if(!kmfit->Fit(0)) return false;
  bool oksq = kmfit->Fit();
  if(oksq) 
  {
    chi2  = kmfit->chisq();
		P.resize(4);
    for(int i=0;i<P.size();i++)
    {
      P[i] = kmfit->pfit(i);
    }
  }
  return oksq;
}

bool kinfit(
		const std::vector<EvtRecTrackIterator> & Tracks,  
		int & channel,  
		double & chi2,  
		HepLorentzVector & P,  
		const double CENTER_MASS_ENERGY
		)
{
	bool goodfit=false;
	for(int pid=0;pid<2;pid++)
	{
		double chi2_tmp=1e100;
		HepLorentzVector P_tmp;
		bool fit_result = kinematic_fit(pid, Tracks, P_tmp, chi2_tmp, CENTER_MASS_ENERGY);
		if(fit_result)
		{
			goodfit=true;
			if(chi2_tmp<chi2)
			{
				channel = pid;
				chi2 = chi2_tmp;
				P = P_tmp;
			}
		}
	}
	return goodfit;
}

bool kinematic_fit(
    int PID, //hypotizes 0 -- KAONS, 1 - MUONS
    TrackPairList_t  & pion_pairs,  //list of pion pairs  (signle pair in simple case)
    TrackPairList_t &  other_pairs, //list of other pairs  (signle pair in simple case)
    std::vector<HepLorentzVector> & P, //?
    double & chi2,  //result of the fit
    TrackPair_t & result_pion_pair,  //best pion pair
    TrackPair_t & result_other_pair,   //best other pair
		double W
    )
{
  chi2=std::numeric_limits<double>::max();
  if(pion_pairs.empty() || other_pairs.empty()) return false;
  bool GoodKinematikFit=false;
  //loop over all pion pairs and all other pairs
  for(TrackPairList_t::iterator pion_pair=pion_pairs.begin(); pion_pair!=pion_pairs.end();pion_pair++)
    for(TrackPairList_t::iterator other_pair=other_pairs.begin(); other_pair!=other_pairs.end();other_pair++)
    {
      std::vector<HepLorentzVector> P_tmp;
      double chi2_tmp=std::numeric_limits<double>::max();
      bool oksq=kinematic_fit(PID, *pion_pair, *other_pair, P_tmp, chi2_tmp, W);
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
    std::cout << "muons "   << setw(15) << event_with_muons << ",  ";
    std::cout << "good knm fit "   << setw(15) << good_kinematic_fit;
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

  std::list<EvtRecTrackIterator> good_charged_tracks=createGoodChargedTrackList(evtRecEvent, evtRecTrkCol);
  std::list<EvtRecTrackIterator> good_neutral_tracks=createGoodNeutralTrackList(evtRecEvent, evtRecTrkCol);

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
        if(MIN_PION_MOMENTUM < p && p<MAX_PION_MOMENTUM) 
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
        if(MIN_PION_MOMENTUM < p &&  p<MAX_PION_MOMENTUM) 
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

  //create pion pairs
  TrackPairList_t pion_pairs;
  for(TrackList_t::iterator i=negative_pion_tracks.begin(); i!=negative_pion_tracks.end(); ++i)
    for(TrackList_t::iterator j=positive_pion_tracks.begin(); j!=positive_pion_tracks.end(); ++j)
    {
      TrackPair_t pair(*i,*j);
      double M_recoil = get_recoil__mass(pair, PION_MASS,  CENTER_MASS_ENERGY);
      if(MIN_RECOIL_MASS < M_recoil && M_recoil < MAX_RECOIL_MASS) 
      {
        pion_pairs.push_back(pair);
      }
    }
  //SELECTION CODE we must have at list one pion pair
  if(pion_pairs.empty()) return StatusCode::SUCCESS; //we must have at list one pion pair
	//if(pion_pairs.size()!=1) return StatusCode::SUCCESS; //exacly one pion pair
  TrackPair_t pion_pair = pion_pairs.front();
	if(pion_pairs.size()>1)
	{
		//find the best pion pair using closest value to JPSI_MASS
		for(TrackPairList_t::iterator p=pion_pairs.begin();p!=pion_pairs.end();p++)
		{
			if(fabs(get_recoil__mass(*p,PION_MASS, CENTER_MASS_ENERGY) - JPSI_MASS) <  fabs(get_recoil__mass(pion_pair,PION_MASS, CENTER_MASS_ENERGY) - JPSI_MASS)) pion_pair = *p;
		}
	}
	//now we have one pion pair candidate

  //SELECTION CODE
  //keep only specific signature
  //if(positive_charged_tracks.size()!=2 || negative_charged_tracks.size()!=2) return StatusCode::SUCCESS;

  //log << MSG::ERROR << "good charged tracks: " << charged_tracks.size() << " (" << negative_charged_tracks.size() << ", " << positive_charged_tracks.size() << endmsg;
  //log << MSG::ERROR << "pions: " << negative_pion_tracks.size()  << ", " << positive_pion_tracks.size() << endmsg;
  //log << MSG::ERROR << "other: " << other_negative_tracks.size()  << ", " << other_positive_tracks.size() << endmsg;
  //log << MSG::ERROR << "pion pairs: " << pion_pairs.size() << endmsg;


  //if no other particles
	if(other_negative_tracks.empty() && other_positive_tracks.empty()) return StatusCode::SUCCESS;

	int sign = (int(!other_positive_tracks.empty()) << 1 ) + int(!other_negative_tracks.empty());

	bool GoodKinematikFit=false;
	double kinematic_chi2=2e100;
	int channel = -1;
	TrackPair_t result_pair;
	std::vector<HepLorentzVector> Pkf(4); //four momentum after kinematic fit
	std::vector<EvtRecTrackIterator> Tracks;
	Tracks.reserve(4);
	
  //one charged particle is missing
	if(other_positive_tracks.empty() || other_negative_tracks.empty())
	{
		TrackList_t * tracks=0;
		if(!other_positive_tracks.empty()) tracks = &other_positive_tracks;
		if(!other_negative_tracks.empty()) tracks = &other_negative_tracks;
		int tmp_chan=-1;
		double tmp_chi2=2e100;
		std::vector<EvtRecTrackIterator> tmp_Tracks(3);
		tmp_Tracks.resize(3);
		tmp_Tracks[0]=pion_pair.first;
		tmp_Tracks[1]=pion_pair.second;
		std::vector<HepLorentzVector> tmp_P;
		for(TrackList_t::iterator i=tracks->begin(); i!=tracks->end(); ++i)
		{
			EvtRecTrackIterator track = *i;
			tmp_Tracks[2] = track;
			bool fit_result=kinfit(tmp_Tracks, tmp_chan, tmp_chi2,tmp_P,  CENTER_MASS_ENERGY);
			if(fit_result)
			{
				GoodKinematikFit=true;
				if(tmp_chi2 < kinematic_chi2)
				{
					kinematic_chi2 = tmp_chi2;
					channel = tmp_chan;
					Pkf = tmp_P;
					Tracks = tmp_Tracks;
				}
			}
		}
		//convert resulting 4-momentum P to Pkf with corresponding sign
		if(sign == 0x4) 
		{
			std::swap(Pkf[2], Pkf[3]);
			if(Tracks.size()==3) 
			{
				Tracks.push_back(evtRecTrkCol->end());
				std::swap(Tracks[2], Tracks[3]);
			}
		}
	}
	else //positive and negartive particles exists then find best pair
	{
		std::vector<EvtRecTrackIterator> tmp_Tracks(3);
		tmp_Tracks.resize(4);
		tmp_Tracks[0]=pion_pair.first;
		tmp_Tracks[1]=pion_pair.second;
		int tmp_chan=-1;
		double tmp_chi2=2e100;
		std::vector<HepLorentzVector> tmp_P;
		//make other pairs
		//TrackPairList_t other_pairs;
		for(TrackList_t::iterator i=other_negative_tracks.begin(); i!=other_negative_tracks.end(); ++i)
		{
			EvtRecTrackIterator first_track = *i;
			for(TrackList_t::iterator j=other_positive_tracks.begin(); j!=other_positive_tracks.end(); ++j)
			{
				EvtRecTrackIterator second_track = *j;
				tmp_Tracks[2] = *i;
				tmp_Tracks[3] = *j;
				bool fit_result=kinfit(tmp_Tracks, tmp_chan, tmp_chi2,tmp_P,  CENTER_MASS_ENERGY);
				if(fit_result)
				{
					GoodKinematikFit=true;
					if(tmp_chi2 < kinematic_chi2)
					{
						kinematic_chi2 = tmp_chi2;
						channel = tmp_chan;
						Pkf = tmp_P;
						Tracks = tmp_Tracks;
					}
				}
			} //end of loop for positive tracks
		} //end of loop for negative tracks 
	}

	//SELECTION CODE
	if(!GoodKinematikFit) return StatusCode::SUCCESS;
	if(kinematic_chi2>200) return StatusCode::SUCCESS;
	good_kinematic_fit++;

	double MIN_MOMENTUM[5] = { MIN_KAON_MOMENTUM,  MIN_MUON_MOMENTUM,  0, 0, 0}; 
	double MAX_MOMENTUM[5] = { MAX_KAON_MOMENTUM,  MAX_MUON_MOMENTUM,  0, 0, 0}; 
	double MIN_EP_RATIO[5] = { MIN_KAON_EP_RATIO,  MIN_MUON_EP_RATIO,  0, 0, 0}; 
	double MAX_EP_RATIO[5] = { MAX_KAON_EP_RATIO,  MAX_MUON_EP_RATIO,  0, 0, 0}; 
	//SELECTION CODE SUPRESS ELECTRONS
	for(int i=2;i<4;i++)
	{
		if(Tracks[i]==evtRecTrkCol->end()) continue;
		RecMdcTrack  * mdcTrk = (*Tracks[i])->mdcTrack();
		RecEmcShower * emcTrk = (*Tracks[i])->emcShower();
		double EpRatio = emcTrk->energy()/mdcTrk->p();
		if( EpRatio < MIN_EP_RATIO[channel]      || MAX_EP_RATIO[channel] < EpRatio )     return StatusCode::SUCCESS;
		if( mdcTrk->p() < MIN_MOMENTUM[channel]  || MAX_MOMENTUM[channel] < mdcTrk->p() ) return StatusCode::SUCCESS;
	}

	vector<double> pchi2(5, 0);
	int scale=0;
	for(int i=2;i<4;i++)
	{
		if(Tracks[i]==evtRecTrkCol->end()) continue;
		vector<double> chi2 = get_chi2(Tracks[i]);
		scale++;
		for(int pid =0;pid<5;pid++)
		{
			pchi2[pid]+=chi2[pid];
		}
	}
	//normalization of chi2
	for(int pid=0;pid<5;pid++)
	{
		pchi2[pid]/scale*2;
	}

	//SELECTION CODE
	if( pchi2[channel] > 200 ) return StatusCode::SUCCESS;
	if( pchi2[channel] > pchi2[ID_KAON] )     return StatusCode::SUCCESS;
	if( pchi2[channel] > pchi2[ID_MUON] )     return StatusCode::SUCCESS;
	//if( pchi2[channel] > pchi2[ID_PION] )     return StatusCode::SUCCESS;
	//if( pchi2[channel] > pchi2[ID_PROTON] )   return StatusCode::SUCCESS;
	//if( pchi2[channel] > pchi2[ID_ELECTRON] ) return StatusCode::SUCCESS;
	
	fEvent.sign = sign;

	switch(channel)
	{
		case ID_KAON:
			if( pchi2[channel] > pchi2[ID_PION] )     return StatusCode::SUCCESS;
			fEvent.KK=1;
			fEvent.uu=0;
			event_with_kaons++;
			break;
		case ID_MUON:
			fEvent.KK=0;
			fEvent.uu=1;
			event_with_muons++;
			break;
		case ID_ELECTRON:
			event_with_electrons++;
			return StatusCode::SUCCESS;
			break;
		case ID_PION:
			event_with_pions++;
			return StatusCode::SUCCESS;
			break;
		case ID_PROTON:
			event_with_protons++;
			return StatusCode::SUCCESS;
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
  fEvent.M012 = (Pkf[0]+Pkf[1]+Pkf[2]).m();
  fEvent.M013 = (Pkf[0]+Pkf[1]+Pkf[3]).m();
  fEvent.M03 =  (Pkf[0]+Pkf[3]).m();
  fEvent.M12 =  (Pkf[1]+Pkf[2]).m();
  fEvent.M01 =  (Pkf[0]+Pkf[1]).m();
  HepLorentzVector P_psip(CENTER_MASS_ENERGY*sin(0.011),0,0,CENTER_MASS_ENERGY); //initial vector of psip
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
  fMdc.Mrecoil = get_recoil__mass(pion_pair, PION_MASS,  CENTER_MASS_ENERGY);
  fMdc.Minv    = sqrt(get_invariant_mass2(result_pair,XMASS[channel]));
  //EvtRecTrackIterator itTrk[4] = {pion_pair.first, pion_pair.second, result_pair.first, result_pair.second};
  for(int i=0;i<4;i++)
  {
		if(Tracks[i]==evtRecTrkCol->end()) continue;
    if(i>1)
    {
      RecEmcShower *emcTrk = (*Tracks[i])->emcShower();
      fEmc.E[i] = emcTrk->energy();
      fEmc.theta[i] = emcTrk->theta();
      fEmc.phi[i] = emcTrk->phi();
      fEmc.time[i] = emcTrk->time();
      fMdc.E[i] = fEmc.E[i];
    }
    RecMdcTrack  *mdcTrk = (*Tracks[i])->mdcTrack();
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


    //dedx information
    if((*Tracks[i])->isMdcDedxValid())
    {
      RecMdcDedx* dedxTrk = (*Tracks[i])->mdcDedx();
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
    if((*Tracks[i])->isTofTrackValid())
    {
      SmartRefVector<RecTofTrack> tofTrkCol = (*Tracks[i])->tofTrack();
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

    //fill particle id

    PID->setRecTrack(Tracks[i]);
    PID->calculate();
    if(PID->IsPidInfoValid())
    {
      fPid.prob[ID_ELECTRON][i] = PID->probElectron();
      fPid.prob[ID_MUON][i]     = PID->probMuon();
      fPid.prob[ID_PION][i]     = PID->probPion();
      fPid.prob[ID_KAON][i]     = PID->probKaon();
      fPid.prob[ID_PROTON][i]   = PID->probProton();
    }
    vector<double> chi2 = get_chi2(Tracks[i]);
    for(int pid=0;pid<5;pid++)
    {
      fPid.chi2[pid][i]   = chi2[pid];
    }

  }
  for(int i=0;i<5;i++)
  {
    fPid.M[i]    = sqrt(get_invariant_mass2(result_pair,XMASS[i]));
    HepLorentzVector p1(Pkf[2].vect(), XMASS[i]);
    HepLorentzVector p2(Pkf[3].vect(), XMASS[i]);
    fPid.kM[i] = (p1+p2).m();
  };


  if(fEvent.run<0)
  {
    //check the MC information
    if(!mcParticleCol)
    {
      log << MSG::ERROR << "Could not retrieve McParticelCol" << endreq;
      return StatusCode::FAILURE;
    }
    //Fill MC TOPO INFORMATION
    bool psipDecay = false;
    int rootIndex = -1;
    Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
    int m_numParticle = 0;
    for (; iter_mc != mcParticleCol->end(); iter_mc++)
    {
      if ((*iter_mc)->primaryParticle()) continue;
      if (!(*iter_mc)->decayFromGenerator()) continue;
      if ((*iter_mc)->particleProperty()==100443)
      {
        psipDecay = true;
        rootIndex = (*iter_mc)->trackIndex();
      }
      if (!psipDecay) continue;
      int pdgid = (*iter_mc)->particleProperty();
      int mcidx = ((*iter_mc)->mother()).trackIndex() - rootIndex;
      fMCTopo.m_pdgid[m_numParticle] = pdgid;
      fMCTopo.m_motheridx[m_numParticle] = mcidx;
      fMCTopo.m_idx[m_numParticle] = (*iter_mc)->trackIndex()-rootIndex;
      fMCTopo.m_hash=0; //no hash calculation now
      m_numParticle += 1;
    }
    fMCTopo.m_idxmc = m_numParticle;

    //Fill my mc truth information
    fMC.ntrack=4;
    HepLorentzVector MCPpion[2];
    HepLorentzVector MCPkaon_or_muon[2];
    bool pi_minus(false);
    bool pi_plus(false);
    bool K_minus(false);
    bool K_plus(false);
    bool mu_minus(false);
    bool mu_plus(false);
    int mytrack=0;
    fMC.psip_decay = 0;
    fMC.jpsi_decay = 0;
    fMC.KK = 0;
    fMC.uu = 0;
    fMC.oo = 0;
    psipDecay=false;
    rootIndex=-1;
    for (iter_mc = mcParticleCol->begin(); iter_mc != mcParticleCol->end(); iter_mc++)
    {
      if ((*iter_mc)->primaryParticle()) continue;
      if (!(*iter_mc)->decayFromGenerator()) continue;
      //if ( ((*iter_mc)->mother()).trackIndex()<3 ) continue;
      if ((*iter_mc)->particleProperty()==100443)
      {
        psipDecay = true;
        fMC.psip_decay=1;
        rootIndex = (*iter_mc)->trackIndex();
      }
      if (!psipDecay) continue;
      if ((*iter_mc)->particleProperty()==443)
      {
        fMC.jpsi_decay=1;
      }
      if (fMC.jpsi_decay!=1) continue;
      if((*iter_mc)->particleProperty() == +211) 
      {
        pi_plus=true;
        MCPpion[1] = (*iter_mc)->initialFourMomentum();
        fMC.pid[1]=211;
        mytrack++;
      }
      if((*iter_mc)->particleProperty() == -211) 
      {
        MCPpion[0] = (*iter_mc)->initialFourMomentum();
        pi_minus=true;
        fMC.pid[0]=-211;
        mytrack++;
      }
      if( ! pi_plus && !pi_minus) continue; //keep only psip to Jpsi pi pi decay
      switch((*iter_mc)->particleProperty())
      {
        case -13:
          MCPkaon_or_muon[0] = (*iter_mc)->initialFourMomentum();
          mu_minus = true;
          mytrack++;
          fMC.pid[2]=-13;
          break;
        case +13:
          MCPkaon_or_muon[1] = (*iter_mc)->initialFourMomentum();
          mu_plus = true;
          mytrack++;
          fMC.pid[3]=13;
          break;
        case -321:
          MCPkaon_or_muon[0] = (*iter_mc)->initialFourMomentum();
          K_minus=true;
          mytrack++;
          fMC.pid[2]=-321;
          break;
        case +321:
          MCPkaon_or_muon[1] = (*iter_mc)->initialFourMomentum();
          mytrack++;
          K_plus = true;
          fMC.pid[3]=321;
          break;
      };
    }
    if(K_plus && K_minus) 
    {
      fMC.KK=1;
      fMC.oo=0;
    }
    if(mu_plus && mu_minus) 
    {
      fMC.uu=1;
      fMC.oo=0;
    }
    if(fMC.KK==1 && fMC.uu==1) fMC.oo=1;
    if(mytrack!=4) fMC.oo=1;
    if(fMC.KK==1 || fMC.uu==1)
    {
      vector<HepLorentzVector> P(4);
      P[0]=MCPpion[0];
      P[1]=MCPpion[1];
      P[2]=MCPkaon_or_muon[0];
      P[3]=MCPkaon_or_muon[1];
      for(int i=0;i<4;i++)
      {
        fMC.q[i] = 0; 
        fMC.E[i] = P[i].e();
        fMC.p[i] = P[i].rho();
        fMC.px[i]= P[i].px();
        fMC.py[i]= P[i].py();
        fMC.pz[i]= P[i].pz();
        fMC.pt[i]= P[i].perp();
        fMC.theta[i]= P[i].theta();
        fMC.phi[i] = P[i].phi();
      }
    }
    else
    {
      Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
      for (; iter_mc != mcParticleCol->end(); iter_mc++)
      {
        if ((*iter_mc)->primaryParticle()) continue;
        if (!(*iter_mc)->decayFromGenerator()) continue;
        HepLorentzVector P = (*iter_mc)->initialFourMomentum();
        int pid = (*iter_mc)->particleProperty();
        Hep3Vector p_mc = P.vect();
        for(int i=0;i<4;i++)
        {
          Hep3Vector p_rec = Pkf[i].vect();
          Hep3Vector dp = p_rec - p_mc;
          if(dp.mag()/std::min(p_rec.mag(), p_mc.mag()) < 0.05)
          {
            fMC.pid[i] = pid;
            //fMC.q[i] = 0; 
            fMC.E[i] = P.e();
            fMC.p[i] = P.rho();
            fMC.px[i]= P.px();
            fMC.py[i]= P.py();
            fMC.pz[i]= P.pz();
            fMC.pt[i]= P.perp();
            fMC.theta[i]= P.theta();
            fMC.phi[i] = P.phi();
          }
        }
      }
    }
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

  if(fEvent.run<0) 
  {
    fMC.tuple->write();
    fMCTopo.tuple->write();
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
  std::cout << "Event with pions: " << event_with_pions << std::endl;
  std::cout << "Event with electron: " << event_with_electrons << std::endl;
  std::cout << "Event with proton: " << event_with_protons << std::endl;
  std::cout << "Good kinematic fits: " << good_kinematic_fit << std::endl;
  return StatusCode::SUCCESS;
}
