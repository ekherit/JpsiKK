// =====================================================================================
//
//       Filename:  MyPid.h
//
//    Description:  
//
//        Version:  1.0
//        Created:  19.10.2015 21:16:19
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#pragma once 

#include <vector>

#include "EvtRecEvent/EvtRecTrack.h"
#include "DstEvent/TofHitStatus.h"

SmartRefVector<RecTofTrack>::iterator  getTofTrk(EvtRecTrack * track, bool & isTofValid)
{
  SmartRefVector<RecTofTrack> tofTrkCol = track->tofTrack();
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



vector<double> get_chi2(EvtRecTrack * Trk)
{
  vector<double> chi2(5,99999);
  if(!Trk->isMdcTrackValid()) return chi2;
  if(!Trk->isMdcDedxValid())  return chi2;
  for(int i=0;i<5;i++) chi2[i]=0;
  RecMdcTrack * mdcTrk  = Trk->mdcTrack();

  //dedx information
  RecMdcDedx  * dedxTrk = Trk->mdcDedx();
  chi2[ID_KAON]     +=   sq(dedxTrk->chiK());
  chi2[ID_MUON]     +=   sq(dedxTrk->chiMu());
  chi2[ID_ELECTRON] +=   sq(dedxTrk->chiE());
  chi2[ID_PION]     +=   sq(dedxTrk->chiPi());
  chi2[ID_PROTON]   +=   sq(dedxTrk->chiP());

  //tof information
  if(!Trk->isTofTrackValid()) return chi2;
  bool isTofValid=false;
  SmartRefVector<RecTofTrack>::iterator tofTrk = getTofTrk(Trk, isTofValid);
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
  EvtRecTrack *  Trk[4] = {pion_pair.first, pion_pair.second, kaon_pair.first, kaon_pair.second};
  vector< vector<double> > chi2(4);
  for(int track=0;track<chi2.size();track++)
  {
    chi2[track] = get_chi2(Trk[track]);
  }
  return chi2;
}

vector<double> get_chi2(TrackPair_t & pair)
{
  EvtRecTrack*  Trk[2] = {pair.first, pair.second};
  vector<double> chi2(5,0);
  for(int track=0;track<2;track++)
  {
    vector<double> tmp_chi2=get_chi2(Trk[track]);
    for(int i=0;i<5;i++)
    {
      chi2[i]  += tmp_chi2[i]; 
    }
  }
  return chi2;
}
