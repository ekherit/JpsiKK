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


const double PI_MESON_MASS=0.13957018; //GeV


const double EMS_THRESHOLD = 0.05; //GeV
const double MAX_MOMENTUM  = 2.5; //GeV

const double EMC_ENDCUP_THRESHOLD=0.05;
const double EMC_BARREL_THRESHOLD=0.025;

inline double sq(double x) { return x*x; }

JpsiKK::JpsiKK(const std::string& name, ISvcLocator* pSvcLocator) :
  Algorithm(name, pSvcLocator)
{
  declareProperty("MAX_CHARGED_TRACKS", MAX_CHARGED_TRACKS=10); //maximum number of charged tracks in selection
  declareProperty("MIN_CHARGED_TRACKS", MIN_CHARGED_TRACKS=3); //minimum number of charged tracks in selection
  declareProperty("IP_MAX_Z", IP_MAX_Z = 10.0); //cm?
  declareProperty("IP_MAX_RHO", IP_MAX_RHO = 1.0); //cm?
  declareProperty("MAX_COS_THETA", MAX_COS_THETA = 0.93); //cm?
  declareProperty("CHECK_TOF", CHECK_TOF=1);
  declareProperty("CHECK_DEDX", CHECK_DEDX = 1);
  declareProperty("CHECK_MUC", CHECK_MUC = 1);
  declareProperty("CHECK_MC", CHECK_MC = 1);
  declareProperty("MC_DP", MC_DP = 0.01);
  declareProperty("USE_IPCUT", USE_IPCUT=1); //to use interection point cut.
  declareProperty("IPR", IPR=1); //Interaction point cut distance.
  declareProperty("IPTRACKS", IPTRACKS=4); //number of tracks from interection point
  declareProperty("MAX_TRACK_NUMBER", MAX_TRACK_NUMBER=4); //maximum number of charged tracks
  declareProperty("STRICT_TAU_CUT", STRICT_TAU_CUT=1); //maximum number of charged tracks

  emc.MAX_TRACK_NUMBER = MAX_TRACK_NUMBER;
  gg.MAX_TRACK_NUMBER = MAX_TRACK_NUMBER;
  muc.MAX_TRACK_NUMBER = MAX_TRACK_NUMBER;
}


StatusCode JpsiKK::initialize(void)
{
  MsgStream log(msgSvc(), name());

  log << MSG::INFO << "in initialize()" << endmsg;
  event_proceed=0;
  event_write = 0;
  tau_events=0;
  bhabha_events=0;
  gg_event_writed=0;

  StatusCode status;
  NTuplePtr my_nt(ntupleSvc(), "FILE1/mhadr");
  if(my_nt) main_tuple=my_nt;
  else
  {
    main_tuple = ntupleSvc()->book("FILE1/mhadr", CLID_ColumnWiseTuple, "Multihadron tree plus bhabha");
    if(main_tuple)
    {
      //common
      status=main_tuple->addItem("t", m_time);
      status=main_tuple->addItem("ntrack", m_ntrack);
      status=main_tuple->addItem("nchtrk", m_nchtr);
      status=main_tuple->addItem("nneutrk", m_nneutr);
      status=main_tuple->addItem("Etotal", m_Etotal);
      status=main_tuple->addItem("Eemc", m_Eemc);
    }
    else
    {
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(main_tuple) << endmsg;
      return StatusCode::FAILURE;
    }
  }

  /*  MDC ntuple intialization */
  NTuplePtr nt_mdc(ntupleSvc(), "FILE1/mdc");
  if(nt_mdc) mdc_tuple=nt_mdc;
  else
  {
    mdc_tuple = ntupleSvc()->book("FILE1/mdc", CLID_ColumnWiseTuple, "Charged  tracks");
    if(mdc_tuple)
    {
      status = mdc_tuple->addItem ("ntrack", mdc.ntrack, 0, MAX_TRACK_NUMBER);
      status = mdc_tuple->addItem ("ngt", mdc.ngood_track, 0, MAX_TRACK_NUMBER); //number of good tracks
      status = mdc_tuple->addItem ("Emdc", mdc.Emdc);
      status = mdc_tuple->addItem ("Eemc", mdc.Eemc);
      status = mdc_tuple->addItem ("nip", mdc.nip);
      status = mdc_tuple->addIndexedItem ("p", mdc.ntrack, mdc.p);
      status = mdc_tuple->addIndexedItem ("pt", mdc.ntrack, mdc.pt);
      status = mdc_tuple->addIndexedItem ("px", mdc.ntrack, mdc.px);
      status = mdc_tuple->addIndexedItem ("py", mdc.ntrack, mdc.py);
      status = mdc_tuple->addIndexedItem ("pz", mdc.ntrack, mdc.pz);
      status = mdc_tuple->addIndexedItem ("theta", mdc.ntrack, mdc.theta);
      status = mdc_tuple->addIndexedItem ("phi", mdc.ntrack, mdc.phi);
      status = mdc_tuple->addIndexedItem ("x", mdc.ntrack, mdc.x);
      status = mdc_tuple->addIndexedItem ("y", mdc.ntrack, mdc.y);
      status = mdc_tuple->addIndexedItem ("z", mdc.ntrack, mdc.z);
      status = mdc_tuple->addIndexedItem ("r", mdc.ntrack, mdc.r);
      status = mdc_tuple->addIndexedItem ("rvxy", mdc.ntrack, mdc.rvxy);
      status = mdc_tuple->addIndexedItem ("rvz", mdc.ntrack, mdc.rvz);
      status = mdc_tuple->addIndexedItem ("rvphi", mdc.ntrack, mdc.rvphi);
      status = mdc_tuple->addIndexedItem ("ismu", mdc.ntrack, mdc.ismu);
      status = mdc_tuple->addIndexedItem ("q", mdc.ntrack, mdc.q);
      // EMC information for charged tracks
      status = mdc_tuple->addIndexedItem ("isemc", mdc.ntrack, mdc.isemc);
      status = mdc_tuple->addIndexedItem ("E", mdc.ntrack, mdc.E);
      status = mdc_tuple->addIndexedItem ("dE", mdc.ntrack, mdc.dE);
      status = mdc_tuple->addIndexedItem ("temc", mdc.ntrack, mdc.temc);
      status = mdc_tuple->addIndexedItem ("ncrstl", mdc.ntrack, mdc.ncrstl);
      status = mdc_tuple->addIndexedItem ("cellId", mdc.ntrack, mdc.cellId);
      status = mdc_tuple->addIndexedItem ("status", mdc.ntrack, mdc.status);
      status = mdc_tuple->addIndexedItem ("module", mdc.ntrack, mdc.module);
      status = mdc_tuple->addIndexedItem ("M", mdc.ntrack, mdc.M);
      status = mdc_tuple->addIndexedItem ("istof", mdc.ntrack, mdc.istof);
      status = mdc_tuple->addIndexedItem ("X", mdc.ntrack, mdc.X);
      status = mdc_tuple->addIndexedItem ("Y", mdc.ntrack, mdc.Y);
      status = mdc_tuple->addIndexedItem ("Z", mdc.ntrack, mdc.Z);

      /*  sphericity part */
      status = mdc_tuple->addItem("Mrec", mdc.Mrec);
      status = mdc_tuple->addItem("MKK", mdc.MKK);
      status = mdc_tuple->addItem("Mmumu", mdc.Mmumu);
      status = mdc_tuple->addItem("Mmiss", mdc.Mmiss);
      status = mdc_tuple->addItem("S", mdc.S);
      status = mdc_tuple->addItem("ccos", mdc.ccos);
      status = mdc_tuple->addItem("atheta", mdc.atheta);
      status = mdc_tuple->addItem("aphi", mdc.aphi);
      status = mdc_tuple->addItem("acompl", mdc.acompl);

      /*  particle identification part */
      status = mdc_tuple->addItem("chan", mdc.jpsi_decay_channel);

      status = mdc_tuple->addIndexedItem ("probe", mdc.ntrack, mdc.probe);
      status = mdc_tuple->addIndexedItem ("probmu", mdc.ntrack, mdc.probmu);
      status = mdc_tuple->addIndexedItem ("probpi", mdc.ntrack, mdc.probpi);
      status = mdc_tuple->addIndexedItem ("probK", mdc.ntrack, mdc.probK);
      status = mdc_tuple->addIndexedItem ("probp", mdc.ntrack, mdc.probp);
    }
    else
    {
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(mdc_tuple) << endmsg;
      return StatusCode::FAILURE;
    }
  }

  NTuplePtr nt_emc(ntupleSvc(), "FILE1/emc");
  if(nt_emc) emc_tuple=nt_emc;
  else
  {
    emc_tuple = ntupleSvc()->book("FILE1/emc", CLID_ColumnWiseTuple, "Netutral track");
    if(emc_tuple)
    {
      status = emc.init_tuple(emc_tuple);
    }
    else
    {
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(emc_tuple) << endmsg;
      return StatusCode::FAILURE;
    }
  }

  NTuplePtr nt_muc(ntupleSvc(), "FILE1/muc");
  if(nt_muc) muc_tuple=nt_muc;
  else
  {
    muc_tuple = ntupleSvc()->book("FILE1/muc", CLID_ColumnWiseTuple, "Muon system");
    if(muc_tuple)
    {
      status = muc.init_tuple(muc_tuple);
    }
    else
    {
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(muc_tuple) << endmsg;
      return StatusCode::FAILURE;
    }
  }

  NTuplePtr nt2(ntupleSvc(), "FILE1/dedx");
  if(nt2) dedx_tuple=nt2;
  else
  {
    dedx_tuple = ntupleSvc()->book("FILE1/dedx", CLID_ColumnWiseTuple, "dedx information");
    if(dedx_tuple)
    {
      status = dedx_tuple->addItem ("dedx.ntrack", dedx.ntrack, 0, MAX_TRACK_NUMBER);
      status = dedx_tuple->addIndexedItem ("pid", dedx.ntrack, dedx.pid );
      status = dedx_tuple->addIndexedItem ("chie", dedx.ntrack, dedx.chie );
      status = dedx_tuple->addIndexedItem ("chimu",dedx.ntrack, dedx.chimu );
      status = dedx_tuple->addIndexedItem ("chipi", dedx.ntrack, dedx.chipi );
      status = dedx_tuple->addIndexedItem ("chik", dedx.ntrack, dedx.chik );
      status = dedx_tuple->addIndexedItem ("chip", dedx.ntrack, dedx.chip );
      status = dedx_tuple->addIndexedItem ("ghit", dedx.ntrack, dedx.ghit );
      status = dedx_tuple->addIndexedItem ("thit", dedx.ntrack, dedx.thit );
      status = dedx_tuple->addIndexedItem ("probPH", dedx.ntrack, dedx.probPH );
      status = dedx_tuple->addIndexedItem ("normPH", dedx.ntrack, dedx.normPH );
      status = dedx_tuple->addIndexedItem ("dedx_e", dedx.ntrack, dedx.e );
      status = dedx_tuple->addIndexedItem ("dedx_mu", dedx.ntrack, dedx.mu );
      status = dedx_tuple->addIndexedItem ("dedx_pi", dedx.ntrack, dedx.pi );
      status = dedx_tuple->addIndexedItem ("dedx_K", dedx.ntrack, dedx.K );
      status = dedx_tuple->addIndexedItem ("dedx_p", dedx.ntrack, dedx.p );
    }
    else
    {
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(dedx_tuple) << endmsg;
      return StatusCode::FAILURE;
    }
  }

  if(CHECK_MC)
  {
    NTuplePtr nt(ntupleSvc(), "FILE1/mc");
    if(nt) mc.tuple=nt;
    else
    {
      mc.tuple = ntupleSvc()->book("FILE1/mc", CLID_ColumnWiseTuple, "Monte Carlo information");
      if(mc.tuple)
      {
        status = mc.init_tuple(mc.tuple, MAX_TRACK_NUMBER);
      }
      else
      {
        log << MSG::ERROR << "    Cannot book N-tuple:" << long(mc.tuple) << endmsg;
        return StatusCode::FAILURE;
      }
    }
  }

  if(CHECK_TOF)
  {

    NTuplePtr nt_tof(ntupleSvc(), "FILE1/tof");
    if(nt_tof) tof_tuple=nt2;
    else
    {
      tof_tuple = ntupleSvc()->book("FILE1/tof", CLID_ColumnWiseTuple, "tof information");
      if(tof_tuple)
      {
        status = tof_tuple->addItem ("ntrack", tof.ntrack, 0, MAX_TRACK_NUMBER);
        status = tof_tuple->addIndexedItem ("trackID", tof.ntrack, tof.trackID );
        status = tof_tuple->addIndexedItem ("tofID", tof.ntrack, tof.tofID);
        status = tof_tuple->addIndexedItem ("tofTrackID", tof.ntrack, tof.tofTrackID);
        status = tof_tuple->addIndexedItem ("status", tof.ntrack, tof.status);
        status = tof_tuple->addIndexedItem ("path", tof.ntrack, tof.path);
        status = tof_tuple->addIndexedItem ("zrhit", tof.ntrack, tof.zrhit);
        status = tof_tuple->addIndexedItem ("ph", tof.ntrack, tof.ph);
        status = tof_tuple->addIndexedItem ("tof", tof.ntrack, tof.tof);
        status = tof_tuple->addIndexedItem ("errtof", tof.ntrack, tof.errtof);
        status = tof_tuple->addIndexedItem ("beta", tof.ntrack, tof.beta);
        status = tof_tuple->addIndexedItem ("texpe", tof.ntrack, tof.texpe);
        status = tof_tuple->addIndexedItem ("texpmu", tof.ntrack, tof.texpmu);
        status = tof_tuple->addIndexedItem ("texppi", tof.ntrack, tof.texppi);
        status = tof_tuple->addIndexedItem ("texpK", tof.ntrack, tof.texpK);
        status = tof_tuple->addIndexedItem ("texpp", tof.ntrack, tof.texpp);
        status = tof_tuple->addIndexedItem ("toffsete", tof.ntrack, tof.toffsete);
        status = tof_tuple->addIndexedItem ("toffsetmu", tof.ntrack, tof.toffsetmu);
        status = tof_tuple->addIndexedItem ("toffsetpi", tof.ntrack, tof.toffsetpi);
        status = tof_tuple->addIndexedItem ("toffsetK", tof.ntrack, tof.toffsetK);
        status = tof_tuple->addIndexedItem ("toffsetp", tof.ntrack, tof.toffsetp);
        status = tof_tuple->addIndexedItem ("toffsetap", tof.ntrack, tof.toffsetap);
        status = tof_tuple->addIndexedItem ("sigmae", tof.ntrack, tof.sigmae);
        status = tof_tuple->addIndexedItem ("sigmamu", tof.ntrack, tof.sigmamu);
        status = tof_tuple->addIndexedItem ("sigmapi", tof.ntrack, tof.sigmapi);
        status = tof_tuple->addIndexedItem ("sigmaK", tof.ntrack, tof.sigmaK);
        status = tof_tuple->addIndexedItem ("sigmap", tof.ntrack, tof.sigmap);
        status = tof_tuple->addIndexedItem ("sigmaap", tof.ntrack, tof.sigmaap);
        status = tof_tuple->addIndexedItem ("quality", tof.ntrack, tof.quality);
        status = tof_tuple->addIndexedItem ("t0", tof.ntrack, tof.t0);
        status = tof_tuple->addIndexedItem ("errt0", tof.ntrack, tof.errt0);
        status = tof_tuple->addIndexedItem ("errz", tof.ntrack, tof.errz);
        status = tof_tuple->addIndexedItem ("phi", tof.ntrack, tof.phi);
        status = tof_tuple->addIndexedItem ("errphi", tof.ntrack, tof.errphi);
        status = tof_tuple->addIndexedItem ("E", tof.ntrack, tof.E);
        status = tof_tuple->addIndexedItem ("errE", tof.ntrack, tof.errE);
      }
      else
      {
        log << MSG::ERROR << "    Cannot book N-tuple:" << long(tof_tuple) << endmsg;
        return StatusCode::FAILURE;
      }
    }
  }

  NTuplePtr nt_gg(ntupleSvc(), "FILE1/gg");
  if(nt_gg) gg_tuple=nt_gg;
  else
  {
    gg_tuple = ntupleSvc()->book("FILE1/gg", CLID_ColumnWiseTuple, "gamma-gamma annihilation");
    if(gg_tuple)
    {
      status = gg.init_tuple(gg_tuple);
    }
    else
    {
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(gg_tuple) << endmsg;
      return StatusCode::FAILURE;
    }
  }

  NTuplePtr nt_head(ntupleSvc(), "FILE1/head");
  if(nt_head) head_tuple=nt_head;
  else
  {
    head_tuple = ntupleSvc()->book("FILE1/head", CLID_ColumnWiseTuple, "Head");
    if(head_tuple)
    {
      status = head_tuple->addItem ("run", head_run);
      status = head_tuple->addItem ("n", head_event_number);
      status = head_tuple->addItem ("nsel", head_event_selected);
      status = head_tuple->addItem ("nchtr", head_ncharged_tracks);
      status = head_tuple->addItem ("nchtr_rms", head_ncharged_tracks_rms);
      status = head_tuple->addItem ("nntr", head_nneutral_tracks);
      status = head_tuple->addItem ("nntr_rms", head_nneutral_tracks_rms);
      status = head_tuple->addItem ("nttr", head_ntotal_tracks);
      status = head_tuple->addItem ("nttr_rms", head_ntotal_tracks_rms);
    }
    else
    {
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(head_tuple) << endmsg;
      return StatusCode::FAILURE;
    }
  }

  nchtr_a.reset();
  nntr_a.reset();
  nttr_a.reset();

  return StatusCode::SUCCESS;
}

void JpsiKK::EMC_t::init(void)
{
  // emc information init.
  ntrack=0;
  ngood_track=0;
  ngood_charged_track=0; //number of good charged tracks
  Etotal=0;
  S=-1000;
  ccos=-1000;
  atheta=-1000;
  aphi=-1000;
  for(int i=0; i<MAX_TRACK_NUMBER; i++)
  {
    status[i]=-1000;
    ncrstl[i]=-1000;
    cellId[i]=-1000;
    module[i]=-1000;
    E[i]=-1000;
    dE[i]=-1000;
    x[i]=-1000;
    y[i]=-1000;
    z[i]=-1000;
    theta[i]=-1000;
    phi[i]=-1000;
    t[i]=-10000;
  }
}

StatusCode JpsiKK::EMC_t::init_tuple(NTuple::Tuple * tuple)
{
  StatusCode status;
  status = tuple->addItem ("ntrack", ntrack, 0, MAX_TRACK_NUMBER);
  status = tuple->addItem ("ngt", ngood_track, 0, MAX_TRACK_NUMBER); //number of good charged tracks.
  status = tuple->addItem ("ngct",    ngood_charged_track, 0, MAX_TRACK_NUMBER); //number of good charged tracks
  status = tuple->addItem ("Etotal", Etotal);
  status = tuple->addItem ("ccos", ccos);
  status = tuple->addItem ("S", S); //sphericity
  status = tuple->addItem ("atheta", atheta);
  status = tuple->addItem ("aphi", aphi);
  //arrays
  status = tuple->addIndexedItem ("module",  ntrack, module );
  status = tuple->addIndexedItem ("status", ntrack, EMC_t::status ); //name interference of status.
  status = tuple->addIndexedItem ("ncrstl", ntrack, ncrstl );
  status = tuple->addIndexedItem ("cellId", ntrack, cellId );
  status = tuple->addIndexedItem ("x", ntrack, x );
  status = tuple->addIndexedItem ("y", ntrack, y );
  status = tuple->addIndexedItem ("z", ntrack, z );
  status = tuple->addIndexedItem ("E", ntrack, E );
  status = tuple->addIndexedItem ("dE",ntrack, dE );
  status = tuple->addIndexedItem ("theta", ntrack, theta );
  status = tuple->addIndexedItem ("phi", ntrack, phi);
  status = tuple->addIndexedItem ("t", ntrack, t);
  return status;
}

StatusCode JpsiKK::MUC_t::init_tuple(NTuple::Tuple * tuple)
{
  StatusCode status;
  status = tuple->addItem ("ntrack", ntrack, 0, MAX_TRACK_NUMBER);
  //arrays
  status = tuple->addIndexedItem ("status",  ntrack, MUC_t::status );
  status = tuple->addIndexedItem ("type", ntrack, type);
  status = tuple->addIndexedItem ("depth", ntrack, depth );
  status = tuple->addIndexedItem ("chi2", ntrack,  chi2 );
  status = tuple->addIndexedItem ("ndf", ntrack, ndf );
  status = tuple->addIndexedItem ("distance", ntrack, distance );
  status = tuple->addIndexedItem ("phi", ntrack, phi);
  status = tuple->addIndexedItem ("nhit", ntrack, nhit);
  status = tuple->addIndexedItem ("nlayer", ntrack, nlayer);
  status = tuple->addIndexedItem ("nhitmax", ntrack, nhitmax);
  status = tuple->addIndexedItem ("brlast", ntrack, brlast);
  status = tuple->addIndexedItem ("eclast", ntrack, eclast);
  return status;
}

void JpsiKK::MUC_t::init(void)
{
  // emc information init.
  ntrack=0;
  for(int i=0; i<MAX_TRACK_NUMBER; i++)
  {
    status[i]=-1000;
    type[i]=-1000;
    depth[i]=-1000;
    chi2[i]=-1000;
    ndf[i]=-1000;
    distance[i]=-1000;
    phi[i] = -1000;
    nhit[i] = -1000;
    nlayer[i] = -1000;
    nhitmax[i] = -1000;
    brlast[i] = -1000;
    eclast[i] = -1000;
  }
}

void JpsiKK::MC_t::init(void)
{
  // emc information init.
  ntrack=0;
  for(int i=0; i<MAX_TRACK_NUMBER; i++)
  {
    E[i]=0;
    px[i]=-1000;
    py[i]=-1000;
    pz[i]=-1000;
    p[i]=-1000;
    id[i]=0;
  }
}

StatusCode JpsiKK::MC_t::init_tuple(NTuple::Tuple * tpl, const int NMAX)
{
  MAX_TRACK_NUMBER = NMAX;
  tuple = tpl;
  StatusCode status;
  status = tuple->addItem ("ntrack", ntrack, 0, NMAX);
  status = tuple->addIndexedItem ("E",  ntrack, E);
  status = tuple->addIndexedItem ("p",  ntrack, p);
  status = tuple->addIndexedItem ("px",  ntrack, px);
  status = tuple->addIndexedItem ("py",  ntrack, py);
  status = tuple->addIndexedItem ("pz",  ntrack, pz);
  status = tuple->addIndexedItem ("id",  ntrack, id);
  return status;
}

void JpsiKK::InitData(long nchtrack, long nneutrack)
{
  m_ntrack=nchtrack+nneutrack;
  m_nchtr=nchtrack;
  m_nneutr=nneutrack;
  m_Etotal=0;
  m_Eemc=0;
  //mdc track informaion init
  mdc.nip=0;
  mdc.Mrec=-1000;
  mdc.MKK=-1000;
  mdc.Mmumu=-1000;
  mdc.Mmiss=-1000;
  mdc.Eemc=0;
  mdc.Emdc=0;
  mdc.S=0;
  mdc.ccos=-1000;
  mdc.atheta=-1000;
  mdc.aphi=-1000;
  mdc.acompl=-1000;
  mdc.jpsi_decay_channel = -1000;
  mdc.ntrack=0;
  mdc.ngood_track=0; //number of good tracks
  for(int i=0;i<MAX_TRACK_NUMBER; i++)
  {
    mdc.p[i]=-1000;
    mdc.px[i]=-1000;
    mdc.py[i]=-1000;
    mdc.pz[i]=-1000;
    mdc.pt[i]=-1000;
    mdc.x[i]=-1000;
    mdc.y[i]=-1000;
    mdc.z[i]=-1000;
    mdc.r[i]=-1000;
    mdc.rvxy[i]=-1000;
    mdc.rvz[i]=-1000;
    mdc.rvphi[i]=-1000;
    mdc.theta[i]=-1000;
    mdc.phi[i]=-1000;
    mdc.q[i]=-1000;
    mdc.isemc[i]=-1000;
    mdc.temc[i]=-1000;
    mdc.ncrstl[i]=-1000;
    mdc.cellId[i]=-1000;
    mdc.status[i]=-1000;
    mdc.module[i]=-1000;
    mdc.E[i]=-1000;
    mdc.dE[i]=-1000;
    mdc.M[i]=-1000;
    mdc.ismu[i]=-1000;
    mdc.istof[i]=-1000;
    mdc.X[i]=-1000;
    mdc.Y[i]=-1000;
    mdc.Z[i]=-1000;

    //dedx information
    dedx.pid[i]=-1000;
    dedx.chie[i] = -1000;
    dedx.chimu[i] = -1000;
    dedx.chipi[i] = -1000;
    dedx.chik[i] = -1000;
    dedx.chip[i] = -1000;
    dedx.ghit[i] = -1000;
    dedx.thit[i] = -1000;
    dedx.probPH[i] = -1000;
    dedx.normPH[i] = -1000;
    dedx.e[i]=-1000;
    dedx.mu[i]=-1000;
    dedx.pi[i]=-1000;
    dedx.K[i]=-1000;
    dedx.p[i]=-1000;

    /* probability part */
    mdc.probe[i]=-1000;
    mdc.probmu[i]=-1000;
    mdc.probpi[i]=-1000;
    mdc.probK[i]=-1000;
    mdc.probp[i]=-1000;
  }

  emc.init();
  gg.init();
  for(int i=0;i<MAX_TRACK_NUMBER;i++)
  {
    if(CHECK_TOF)
    {
      tof.trackID[i]=-1000;
      tof.tofID[i]=-1000;
      tof.tofTrackID[i]=-1000;
      tof.status[i] =-1000;
      tof.path[i]  =-1000;
      tof.zrhit[i]  =-1000;
      tof.ph[i]  =-1000;
      tof.tof[i]  =-1000;
      tof.errtof[i]  =-1000;
      tof.beta[i]  =-1000;
      tof.texpe[i]  =-1000;
      tof.texpmu[i]  =-1000;
      tof.texppi[i]  =-1000;
      tof.texpK[i]  =-1000;
      tof.texpp[i]  =-1000;
      tof.toffsete[i]  =-1000;
      tof.toffsetmu[i]  =-1000;
      tof.toffsetpi[i]  =-1000;
      tof.toffsetK[i]  =-1000;
      tof.toffsetp[i]  =-1000;
      tof.toffsetap[i]  =-1000;
      tof.sigmae[i]  =-1000;
      tof.sigmamu[i]  =-1000;
      tof.sigmapi[i]  =-1000;
      tof.sigmaK[i]  =-1000;
      tof.sigmap[i]  =-1000;
      tof.sigmaap[i]  =-1000;
      tof.quality[i] = -1000;
      tof.t0[i]  =-1000;
      tof.errt0[i]  =-1000;
      tof.errz[i]  =-1000;
      tof.phi[i]  =-1000;
      tof.E[i]  =-1000;
      tof.errE[i]  =-1000;
    }
  }
  if(CHECK_MC) mc.init();
}


void calculate_vertex(RecMdcTrack *mdcTrk, double & ro, double  & z, double phi)
{
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
    std::cout << "  tau-mu:" << setw(15) << tau_events;
    std::cout << "  ee:" << setw(15) << bhabha_events;
    std::cout << "  gg:" << setw(15) << gg_event_writed;
    std::cout << std::endl;
  }
  event_proceed++;

  /*  Get information about reconstructed events */
  SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
  SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol);
  SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(),  EventModel::MC::McParticleCol);


  InitData(evtRecEvent->totalCharged(), evtRecEvent->totalNeutral());

  nchtr_a.add(evtRecEvent->totalCharged());
  nntr_a.add(evtRecEvent->totalNeutral());
  nttr_a.add(evtRecEvent->totalTracks());


  typedef std::multimap <double, unsigned> mmap_t;
  typedef std::pair <double, unsigned> pair_t;
  mmap_t pmap;
  unsigned good_charged_tracks = 0;


  if(evtRecEvent->totalCharged()  >= MIN_CHARGED_TRACKS)
  {
    /*  loop over charged track */
    mdc.ntrack=0;
    //look thru the charged tracks and sort them on energy
    //count good charged tracks wich come from interaction point
    //and has cos(theta) < 0.93
    for(unsigned idx = 0; idx < evtRecEvent->totalCharged(); idx++)
    {
      EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + idx;
      if(!(*itTrk)->isMdcTrackValid()) continue;  //use only valid charged tracks
      RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();  //main drift chambe
      //calculate interaction point distance
      double rvxy=-9999,rvz=-9999,rvphi=-9999;
      calculate_vertex(mdcTrk,rvxy,rvz,rvphi); //find distance to interaction point
      bool is_fromIP = fabs(rvz)< IP_MAX_Z && fabs(rvxy)<IP_MAX_RHO;  //tracks begin near interaction point
      bool is_good_track = is_fromIP && fabs(cos(mdcTrk->theta()))<MAX_COS_THETA; //track is good
      if(!is_good_track) continue;
      //if(!(*itTrk)->isEmcShowerValid()) continue; //charged track must have energy deposition in EMC
      if((*itTrk)->isEmcShowerValid()) 
      {
        RecEmcShower *emcTrk = (*itTrk)->emcShower();
        double E = emcTrk->energy();
        //Emap.insert(pair_t(E,idx));
      }
      pmap.insert(pair_t(mdcTrk->p(),idx));
    }
    std::cout << "Number of good charged tracks: " << pmap.size() << " MIN_CHARGED_TRACKS = " << MIN_CHARGED_TRACKS << ", MAX_CHARGED_TRACKS = " << MAX_CHARGED_TRACKS << std::endl;
    /* Two or more charged tracks witch signal in EMC */
    good_charged_tracks=pmap.size();
    //if no valid charged tracks
    //if(MIN_CHARGED_TRACKS <  Emap.size() || Emap.size() < MAX_CHARGED_TRACKS) goto SKIP_CHARGED;
    if(good_charged_tracks < MIN_CHARGED_TRACKS   || MAX_CHARGED_TRACKS < good_charged_tracks) goto SKIP_CHARGED;
    std::cout << "Number of good charged tracks: " << MIN_CHARGED_TRACKS << " < " << good_charged_tracks << " < " << MAX_CHARGED_TRACKS  << std::endl;

    //now fill the arrayes using indexes sorted by energy
    mdc.ntrack =pmap.size(); //save number of good charged tracks
    muc.ntrack =pmap.size();
    dedx.ntrack=pmap.size();
    tof.ntrack =pmap.size();
    Sphericity S;

    //particle id 
    ParticleID *pid = ParticleID::instance();
    //loop over tracks oredered by energy
    int gidx=0; //good charged track idx
    int npip=0; //number of positive pions
    int npin=0; //number of negative pions
    int nKp=0; //number of found positive kaons
    int nKm=0; //number of found negative kaons
    int nmup=0; //number of found positive muons
    int nmum=0; //number of found negative muons
    int nK=0; //total number of kaons
    int nmu=0; //total number of muons
    int npi=0; //total number of pions
    int pip_idx=-999; //pion index
    int pin_idx=-999; //pion index
    int Kmup_idx=-999;
    int Kmum_idx=-999;
    for(mmap_t::reverse_iterator ri=pmap.rbegin(); ri!=pmap.rend(); ++ri)
    {
      EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + ri->second;
      RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();  //main drift chambe
      //RecEmcShower *emcTrk = (*itTrk)->emcShower(); //Electro Magnet Calorimeer
      double rvxy=-9999,rvz=-9999,rvphi=-9999;
      calculate_vertex(mdcTrk,rvxy,rvz,rvphi); //find distance to interaction point
      //select good tracks before
      //bool is_fromIP = fabs(rvz)<10 && fabs(rvxy)<1.0;  //tracks begin near interaction point
      //bool is_good_track = is_fromIP && fabs(cos(mdcTrk->theta()))<0.93; //track is good
      //if(!is_good_track) continue;
      int i = gidx; //now fill

      //fill vertex information
      mdc.rvxy[i]=rvxy;
      mdc.rvz[i]=rvz;
      mdc.rvphi[i]=rvphi;

      //fil track information
      mdc.p[i]     =  mdcTrk->p();
      mdc.pt[i]    =  mdcTrk->p()*sin(mdcTrk->theta());
      mdc.px[i]    =  mdcTrk->px();
      mdc.py[i]    =  mdcTrk->py();
      mdc.pz[i]    =  mdcTrk->pz();
      mdc.theta[i] =  mdcTrk->theta();
      mdc.phi[i]   =  mdcTrk->phi();
      mdc.q[i]     =  mdcTrk->charge();
      mdc.x[i]     =  mdcTrk->x();
      mdc.y[i]     =  mdcTrk->y();
      mdc.z[i]     =  mdcTrk->z();

      /*  Particle identification game */
      pid->init();
      pid->setMethod(pid->methodProbability());
      pid->setChiMinCut(4);

      pid->setRecTrack(*itTrk);
      //pid->usePidSys((pid->useMuc() | pid->useEmc()) | pid->useDedx()); // use PID sub-system
      //pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2() | pid->useTofE() | pid->useTofQ() | pid->useEmc() | pid->useMuc());
      //pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2() | pid->useTofE() | pid->useTofQ() | pid->useMuc());
      pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2());
      //pid->identify(pid->onlyMuon() | pid->onlyElectron()); 
      pid->identify(pid->all()); 
      pid->calculate();
      if(pid->IsPidInfoValid())
      {
        mdc.probe[i] =  pid->probElectron();
        mdc.probmu[i] = pid->probMuon();
        mdc.probpi[i] = pid->probPion();
        mdc.probK[i] =  pid->probKaon();
        mdc.probp[i] =  pid->probProton();
      }

      //if momentum below 0.5 GeV it could be pions
      if(mdc.p[i] < 0.5)
      {
        npi++;
        mdc.M[i]=0.13956995;
        if(mdc.q[i]>0) 
        {
          npip++;
          pip_idx = i;
        }
        else
        { 
          pin_idx = i;
          npin++;
        }
      }

      //check for pions or muons
      if(mdc.p[i]>1.0)
      {
        if (
            mdc.probmu[i] > mdc.probe[i] &&
            mdc.probmu[i] > mdc.probpi[i] &&
            mdc.probmu[i] > mdc.probK[i] &&
            mdc.probmu[i] > mdc.probp[i]
           )
        {
          nmu++;
          mdc.M[i]=0.105658389;
          if(mdc.q[i]>0) 
          {
            nmup++;
            Kmup_idx = i;
          }
          else
          {
            nmum++;
            Kmum_idx = i;
          }
        }
        if (
            mdc.probK[i] > mdc.probe[i] &&
            mdc.probK[i] > mdc.probpi[i] &&
            mdc.probK[i] > mdc.probmu[i] &&
            mdc.probK[i] > mdc.probp[i]
           )
        {
          nK++;
          mdc.M[i]=0.493677;
          if(mdc.q[i]>0) 
          {
            nKp++;
            Kmup_idx = i;
          }
          else
          {
            nKm++;
            Kmum_idx = i;
          }
        }
      }


      mdc.Emdc+=sqrt(mdc.p[i]*mdc.p[i]+PI_MESON_MASS*PI_MESON_MASS);


      /* Calculate sphericity tensor */
      S.add(mdcTrk->p3());

      if((*itTrk)->isEmcShowerValid())
      {
        mdc.isemc[i] = 1;
        RecEmcShower *emcTrk = (*itTrk)->emcShower(); //Electro Magnet Calorimeer
        // Add EMC information
        mdc.E[i]     =  emcTrk->energy();
        mdc.dE[i]    =  emcTrk->dE();
        mdc.ncrstl[i] = emcTrk->numHits();
        mdc.status[i] = emcTrk->status();
        mdc.cellId[i] = emcTrk->cellId();
        mdc.module[i] = emcTrk->module();
        mdc.temc[i] = emcTrk->time();
        mdc.Eemc+=mdc.E[i]; //Accumulate energy deposition
      }
      else
      {
        mdc.isemc[i]=0;
      }

      HepLorentzVector P(mdc.px[i], mdc.py[i], mdc.pz[i], mdc.E[i]);
      //mdc.M[i]=P.m();

      /* Check muon system information for this track */
      mdc.ismu[i]=(*itTrk)->isMucTrackValid();
      if((*itTrk)->isMucTrackValid() && CHECK_MUC==1)
      {
        RecMucTrack *mucTrk = (*itTrk)->mucTrack();  //main drift chambe
        muc.status[i]= mucTrk->status();
        muc.type[i]= mucTrk->type();
        muc.depth[i]= mucTrk->depth();
        muc.chi2[i]= mucTrk->chi2();
        muc.ndf[i]= mucTrk->dof();
        muc.distance[i]= mucTrk->distance();
        muc.phi[i]= mucTrk->deltaPhi();
        muc.nhit[i] = mucTrk->numHits();
        muc.nlayer[i] = mucTrk->numLayers();
        muc.nhitmax[i] = mucTrk->maxHitsInLayer();
        muc.brlast[i] = mucTrk->brLastLayer();
        muc.eclast[i] = mucTrk->ecLastLayer();
      }


      /* dEdx information */
      if(CHECK_DEDX == 1 && (*itTrk)->isMdcDedxValid())
      {
        RecMdcDedx* dedxTrk = (*itTrk)->mdcDedx();
        dedx.chie[i] = dedxTrk->chiE();
        dedx.chimu[i] = dedxTrk->chiMu();
        dedx.chipi[i] = dedxTrk->chiPi();
        dedx.chik[i] = dedxTrk->chiK();
        dedx.chip[i] = dedxTrk->chiP();
        dedx.ghit[i] = dedxTrk->numGoodHits();
        dedx.thit[i] = dedxTrk->numTotalHits();
        dedx.probPH[i] = dedxTrk->probPH();
        dedx.normPH[i] = dedxTrk->normPH();
        dedx.e[i] = dedxTrk->getDedxExpect(0);
        dedx.mu[i] = dedxTrk->getDedxExpect(1);
        dedx.pi[i] = dedxTrk->getDedxExpect(2);
        dedx.K[i] = dedxTrk->getDedxExpect(3);
        dedx.p[i] = dedxTrk->getDedxExpect(4);
        dedx.pid[i]=dedxTrk->particleId();
      }

      /* check TOF information */
      mdc.istof[i]=(*itTrk)->isTofTrackValid();
      if(CHECK_TOF && mdc.istof[i])
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
          //if(  (hitst->is_barrel()) ) continue;
          //if( !(hitst->is_counter()) ) continue;
          tofecount.push_back(goodtofetrk);
        }
        delete hitst;
        if(tofecount.size()>0) //not tof2 track or more than 1 tracks
        {
          tofTrk = tofTrkCol.begin()+tofecount[0];

          tof.trackID[i]=(*tofTrk)->trackID();
          tof.tofID[i]=(*tofTrk)->tofID();
          tof.tofTrackID[i]=(*tofTrk)->tofTrackID();
          tof.status[i] = (*tofTrk)->status();
          tof.path[i]  = (*tofTrk)->path();
          tof.zrhit[i]  = (*tofTrk)->zrhit();
          tof.ph[i]  = (*tofTrk)->ph();
          tof.tof[i]  = (*tofTrk)->tof();
          tof.errtof[i]  = (*tofTrk)->errtof();
          tof.beta[i]  = (*tofTrk)->beta();
          tof.texpe[i]  = (*tofTrk)->texpElectron();
          tof.texpmu[i]  = (*tofTrk)->texpMuon();
          tof.texppi[i]  = (*tofTrk)->texpPion();
          tof.texpK[i]  = (*tofTrk)->texpKaon();
          tof.texpp[i]  = (*tofTrk)->texpProton();
          tof.toffsete[i]  = (*tofTrk)->toffsetElectron();
          tof.toffsetmu[i]  = (*tofTrk)->toffsetMuon();
          tof.toffsetpi[i]  = (*tofTrk)->toffsetPion();
          tof.toffsetK[i]  = (*tofTrk)->toffsetKaon();
          tof.toffsetp[i]  = (*tofTrk)->toffsetProton();
          tof.toffsetap[i]  = (*tofTrk)->toffsetAntiProton();
          tof.sigmae[i]  = (*tofTrk)->sigmaElectron();
          tof.sigmamu[i]  = (*tofTrk)->sigmaMuon();
          tof.sigmapi[i]  = (*tofTrk)->sigmaPion();
          tof.sigmaK[i]  = (*tofTrk)->sigmaKaon();
          tof.sigmap[i]  = (*tofTrk)->sigmaProton();
          tof.sigmaap[i]  = (*tofTrk)->sigmaAntiProton();
          tof.t0[i]  = (*tofTrk)->t0();
          tof.errt0[i]  = (*tofTrk)->errt0();
          tof.errz[i]  = (*tofTrk)->errz();
          tof.phi[i]  = (*tofTrk)->phi();
          tof.E[i]  = (*tofTrk)->energy();
          tof.errE[i]  = (*tofTrk)->errenergy();
        }
      }

      if(CHECK_MC)
      {
        for(Event::McParticleCol::iterator ip=mcParticleCol->begin(); ip!=mcParticleCol->end(); ++ip)
        {
          Event::McParticle * p = *ip;
          int mc_track_id = p->trackIndex();
          int pid = p->particleProperty();
          HepLorentzVector P4 = p->initialFourMomentum();
          Hep3Vector P3 = P4.vect();
          double angle = P3.angle(mdcTrk->p3());
//          Hep3Vector dP = P3 - mdcTrk->p3();
          if(angle<MC_DP)
          {
            mc.px[i] = P3.x();
            mc.py[i] = P3.y();
            mc.pz[i] = P3.z();
            mc.p[i] = P3.mag();
            mc.E[i] = P4.e();
            mc.id[i] = pid;
            //cout << "mc track=" << mdcTrk->trackId() <<   " mc track=" << mc_track_id << " pid=" << pid <<endl;
          }
        }
      }
      gidx++;
    }
    mdc.ntrack=gidx;
    mdc.ngood_track = gidx;
    if(CHECK_MC) mc.ntrack=gidx;

    /* Calculate acolinearity  for two tracks with big enrgies */
    //Hep3Vector p0(mdc.px[0], mdc.py[0],mdc.pz[0]);
    //Hep3Vector p1(mdc.px[1], mdc.py[1],mdc.pz[1]);
    //mdc.ccos = p0.dot(p1)/(p0.mag()*p1.mag());
    //mdc.atheta = mdc.theta[0]+mdc.theta[1] - M_PI;
    //mdc.aphi =  fabs(mdc.phi[0]-mdc.phi[1]) - M_PI;
    //mdc.acompl = (mdc.px[0]*mdc.py[1]-mdc.py[0]*mdc.px[1])/(mdc.p[0]*mdc.p[1]);
    ////normalize sphericity tensor
    //S.norm();
    ///* fill sphericity */
    //mdc.S = S();
    
    //we always should have two pions.This is taging of J/psi particle
    if(npip!=1 || npin!=1)  goto SKIP_CHARGED;
    //calculate pion energy 
    double Epin = sqrt(mdc.p[pin_idx]*mdc.p[pin_idx] + mdc.M[pin_idx]*mdc.M[pin_idx]);
    double Epip = sqrt(mdc.p[pip_idx]*mdc.p[pip_idx] + mdc.M[pip_idx]*mdc.M[pip_idx]);
    HepLorentzVector P_psip(0.040546,0,0,3.686); //initial vector of psip
    HepLorentzVector P_pip(mdc.px[pip_idx],mdc.py[pip_idx],mdc.pz[pip_idx], Epip); //pion vector
    HepLorentzVector P_pin(mdc.px[pin_idx],mdc.py[pin_idx],mdc.pz[pin_idx], Epin); //pion vector
    HepLorentzVector P_recoil = P_psip  - P_pip - P_pin;
    mdc.Mrec = P_recoil.m(); //recoil mass of two pions
    if(mdc.Mrec < 3.0 || 3.2 < mdc.Mrec) goto SKIP_CHARGED;

    //and other must be Kaons or muons
    //if( (nKp!=1 || nKm!=1) && (nmup!=1 || nmum!=1)) goto SKIP_CHARGED;
    HepLorentzVector Pp(0,0,0,0);
    HepLorentzVector Pm(0,0,0,0);
    if(nKp>0 || nmup>0)
    {
      double Ep = sqrt(mdc.p[Kmup_idx]*mdc.p[Kmup_idx] + mdc.M[Kmup_idx]*mdc.M[Kmup_idx]);
      Pp=HepLorentzVector(mdc.px[Kmup_idx],mdc.py[Kmup_idx],mdc.pz[Kmup_idx], Ep);
    }
    if(nKm>0 || nmum>0)
    {
      double Em = sqrt(mdc.p[Kmum_idx]*mdc.p[Kmum_idx] + mdc.M[Kmum_idx]*mdc.M[Kmum_idx]);
      Pm=HepLorentzVector(mdc.px[Kmum_idx],mdc.py[Kmum_idx],mdc.pz[Kmum_idx], Em);
    }
    HepLorentzVector P = Pp + Pm;
    HepLorentzVector Pmis = P_psip - P_pip - P_pin - Pp - Pm;
    mdc.Mmiss = Pmis.m2();
    if( (nKp!=1 && nKm!=1) || (nmup!=1 && nmum!=1))
    {
    //calculate the for momentum
      if(nKp==1) mdc.MKK = P.m();
      if(nmup==1) mdc.Mmumu = P.m();
    }

    //tag KK decay channel. We could register only 1 kaon
    if( nKp < 2 && nKm < 2 && 0 < (nKp + nKm) && (nKp + nKm) < 3 && (nmup + nmum) == 0 ) mdc.jpsi_decay_channel = 0;
    //tag mumu decay channel. We could register only one muon
    if( nmup <2 && nmum < 2 && 0 < (nmup + nmum) && (nmup + nmum) < 3 && (nKp + nKm) == 0 ) mdc.jpsi_decay_channel = 1;
    //could not find required configuration
    if(mdc.jpsi_decay_channel < 0) goto SKIP_CHARGED;
    

    /* ================================================================================= */
    /*  fill data for neutral tracks */
    int track=0; //index for neutral tracks
    emc.Etotal=0;
    emc.ngood_charged_track=good_charged_tracks;
    mmap_t Emap; //multi map to sort good tracks by energy order
    Emap.clear();
    pmap.clear();
    //calculate good tracks
    for(int idx = evtRecEvent->totalCharged(); idx<evtRecEvent->totalTracks(); idx++)
    {
      EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + idx;
      if(!(*itTrk)->isEmcShowerValid()) continue;
      RecEmcShower *emcTrk = (*itTrk)->emcShower();
      double c =  fabs(cos(emcTrk->theta())); //abs cos theta
      double E  =  emcTrk->energy();
      bool barrel = c <= 0.8;
      bool endcup = (0.86 <=c) && (c <=0.92);
      if( (E<EMC_BARREL_THRESHOLD && barrel) || (E<EMC_ENDCUP_THRESHOLD && endcup) ) continue; 
      Emap.insert(pair_t(E,idx));
    }
    if(Emap.size()>0 && MAX_TRACK_NUMBER < Emap.size()) goto SKIP_CHARGED;  //no good neutral tracks
    int gnidx=0; //good neutrack track idx
    for(mmap_t::reverse_iterator ri=Emap.rbegin(); ri!=Emap.rend(); ++ri, ++gnidx)
    {
      EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + ri->second;
      RecEmcShower *emcTrk = (*itTrk)->emcShower();
      emc.status[gnidx] = emcTrk->status();
      emc.ncrstl[gnidx] = emcTrk->numHits();
      emc.cellId[gnidx] = emcTrk->cellId();
      emc.module[gnidx] = emcTrk->module();
      emc.x[gnidx] = emcTrk->x();
      emc.y[gnidx] = emcTrk->y();
      emc.z[gnidx] = emcTrk->z();
      emc.phi[gnidx] = emcTrk->phi();
      emc.theta[gnidx] = emcTrk->theta();
      emc.E[gnidx]  =  emcTrk->energy();
      emc.dE[gnidx] =  emcTrk->dE();
      emc.Etotal+=emcTrk->energy();
      emc.t[gnidx] = emcTrk->time();
    }
    emc.ntrack=gnidx;


    m_nchtr=evtRecEvent->totalCharged();
    m_nneutr=evtRecEvent->totalNeutral();
    m_ntrack=evtRecEvent->totalCharged()+evtRecEvent->totalNeutral();
    m_Etotal = emc.Etotal+mdc.Emdc;
    m_Eemc = emc.Etotal+mdc.Eemc;
    m_time = eventHeader->time();


    //оставим два трека только для электронов и мюонов
    //if(emc.ntrack>0) goto SKIP_CHARGED;
    //if(mdc.ntrack<2 || mdc.ntrack>2) goto SKIP_CHARGED;

    vector <bool> ise(2);
    vector <bool> ismu(2);
    for(int i=0;i<ise.size();++i)
    { 
      double evp = mdc.E[i] /mdc.p[i];
      //ise[i] = mdc.ismu[i] == 0 && evp > 0.9 && evp < 1.05;
      ise[i] =  evp > 0.8 && evp < 1.2;
      //ismu[i] = mdc.ismu[i] == 1 && mdc.E[i] > 0.1 && mdc.E[i]<0.4;
      ismu[i] =  mdc.E[i] > 0.1 && mdc.E[i]<0.4 && evp < 0.8;
    }
    bool charge = (mdc.q[0]*mdc.q[1]) <  0;
    bool kinem = mdc.p[0] < 1.5 && mdc.p[1] < 1.5;
    bool tau_sig  = ( (ise[0] && ismu[1]) || (ise[1] && ismu[0]) ) && charge && kinem;
    bool bhabha_sig = ( (ise[0] && ise[1]) ) && charge;
    if(bhabha_sig) bhabha_events++;
    if(tau_sig) tau_events++;

    //if(!charge) goto SKIP_CHARGED;

    //if(!tau_sig && !bhabha_sig) goto SKIP_CHARGED;  
    //if(!tau_sig && STRICT_TAU_CUT == 1) goto SKIP_CHARGED;  
    //if(! (ise[0] || ise[1] || ismu[0] || ismu[1]) && STRICT_TAU_CUT == 1) goto SKIP_CHARGED;  

    /* now fill the data */
    main_tuple->write();
    mdc_tuple->write();
    emc_tuple->write();
    if(CHECK_DEDX) dedx_tuple->write();
    if(CHECK_MUC) muc_tuple->write();
    if(CHECK_TOF) tof_tuple->write();
    if(CHECK_MC) mc.tuple->write();
    event_write++;
  }
//selection of gamma-gamma events
SKIP_CHARGED:
  mmap_t Emap;
  gg.ngood_charged_track = good_charged_tracks;
  if(good_charged_tracks==0)
  {
    //select and sort only good neutral tracks.
    Emap.clear();
    for(int track = evtRecEvent->totalCharged(); track < evtRecEvent->totalTracks(); track++)
    {
      EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + track;
      if(!(*itTrk)->isEmcShowerValid()) continue;
      RecEmcShower *emcTrk = (*itTrk)->emcShower();
      double E = emcTrk->energy();
      double c =  fabs(cos(emcTrk->theta())); //abs cos theta
      bool barrel = c <= 0.8;
      bool endcup = (0.86 <=c) && (c <=0.92);
      //save only good photons
      if( (E<EMC_BARREL_THRESHOLD && barrel) || (E<EMC_ENDCUP_THRESHOLD && endcup) ) continue; 
      Emap.insert(pair_t(E,track));
    }
    //Select exactly two good photons
    if(Emap.size() < 2) goto SKIP_GG;
    if(MAX_TRACK_NUMBER < Emap.size()) goto SKIP_GG;
    gg.ngood_track=Emap.size();
    gg.ntrack= Emap.size();
    Sphericity SS;
    vector < Hep3Vector> R(Emap.size());
    int idx=0;
    for(mmap_t::reverse_iterator ri=Emap.rbegin(); ri!=Emap.rend(); ++ri)
    {
      EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + ri->second;
      assert((*itTrk)->isEmcShowerValid()); //check that EMS data is present
      RecEmcShower *emcTrk = (*itTrk)->emcShower();
      gg.x[idx] = emcTrk->x();
      gg.y[idx] = emcTrk->y();
      gg.z[idx] = emcTrk->z();
      gg.theta[idx] = emcTrk->theta();
      gg.phi[idx] = emcTrk->phi();
      gg.E[idx]  =  emcTrk->energy();
      gg.dE[idx] =  emcTrk->dE();
      gg.status[idx] = emcTrk->status();
      gg.ncrstl[idx] = emcTrk->numHits();
      gg.module[idx] = emcTrk->module();
      gg.cellId[idx] = emcTrk->cellId();
      gg.Etotal+=gg.E[idx];
      /* Calculate sphericity tensor */
      R[idx] = Hep3Vector(emcTrk->x(),emcTrk->y(),emcTrk->z());
      SS.add(R[idx]);
      idx++;
    }
    SS.norm();
    gg.S = SS();


    //calculate colliniarity of two high energy tracks
    gg.ccos = R[0].dot(R[1])/(R[0].mag()*R[1].mag());
    gg.atheta = gg.theta[0]+gg.theta[1] - M_PI;
    gg.aphi = fabs(gg.phi[0]-gg.phi[1]) - M_PI;

    bool acol = fabs(gg.atheta) < 0.05 && -0.06 < gg.aphi && gg.aphi<0.02;
    bool barrel = fabs(cos(gg.theta[0])) < 0.8 && fabs(cos(gg.theta[1])) < 0.8;
    bool highE = gg.E[0] > 0.8*1.4 && gg.E[1] > 0.8*1.4  && gg.E[0] < 1.2*1.9 && gg.E[1] < 1.2*1.9;

    bool isgg =  gg.ngood_track==2 && acol && barrel && highE;
    if(!isgg) goto SKIP_GG;

    gg_tuple->write();
    gg_event_writed++;
  }
SKIP_GG:
  return StatusCode::SUCCESS;
}

StatusCode JpsiKK::finalize()
{
  std::cout << "Event proceed: " << event_proceed << std::endl;
  std::cout << "Event selected: " << event_write << std::endl;
  std::cout << "Tau candidates: " << tau_events << endl;
  std::cout << "Bhabha candidates: " << bhabha_events << endl;
  std::cout << "Gamma-Gamma candidates: " << gg_event_writed << endl;
  std::cout << "Selection efficiency: " << event_write/double(event_proceed) << std::endl;
  std::cout << "Average number of total tracks: " << nttr_a.average() << ", rms=" << nttr_a.rms() << endl;
  std::cout << "Average number of charged tracks: " << nchtr_a.average() << ", rms=" << nchtr_a.rms() << endl;
  std::cout << "Average number of neutral tracks: " << nntr_a.average() << ", rms=" << nntr_a.rms() << endl;
  head_event_selected=event_write;
  head_ncharged_tracks=nchtr_a.average();
  head_ncharged_tracks_rms=nchtr_a.rms();
  head_nneutral_tracks=nntr_a.average();
  head_nneutral_tracks_rms=nntr_a.rms();
  head_ntotal_tracks=nttr_a.average();
  head_ntotal_tracks_rms=nttr_a.rms();
  head_tuple->write();
  return StatusCode::SUCCESS;
}


// for particle id look /ihepbatch/bes/alex/workarea/Analysis/Physics/PsiPrime/G2MuMuAlg-00-00-01/PipiJpsiAlg/src
