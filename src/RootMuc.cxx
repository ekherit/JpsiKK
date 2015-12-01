#include "RootMuc.h"

void RootMuc::init_tuple(void)
{
  tuple->addItem ("ntrack",       ntrack, 0, 4); //good nuetral track in event
  tuple->addIndexedItem ("status",   ntrack, status);
  tuple->addIndexedItem ("type",     ntrack, type);
  tuple->addIndexedItem ("depth",    ntrack, depth);
  tuple->addIndexedItem ("chi2",     ntrack, chi2);
  tuple->addIndexedItem ("ndf",      ntrack, ndf);
  tuple->addIndexedItem ("distance", ntrack, distance);
  tuple->addIndexedItem ("phi",      ntrack, phi);
  tuple->addIndexedItem ("nhit",     ntrack, nhit);
  tuple->addIndexedItem ("nlayer",   ntrack, nlayer);
  tuple->addIndexedItem ("nhitmax",  ntrack, nhitmax);
  tuple->addIndexedItem ("brlast",   ntrack, brlast);
  tuple->addIndexedItem ("eclast",   ntrack, eclast);
}


void RootMuc::init(void)
{
  ntrack=4;
}

void RootMuc::fill(int i,  EvtRecTrackIterator & itTrk)
{
  if((*itTrk)->isMucTrackValid())
  {
    RecMucTrack *mucTrk = (*itTrk)->mucTrack();
    status[i]= mucTrk->status();
    type[i]= mucTrk->type();
    depth[i]= mucTrk->depth();
    chi2[i]= mucTrk->chi2();
    ndf[i]= mucTrk->dof();
    distance[i]= mucTrk->distance();
    phi[i]= mucTrk->deltaPhi();
    nhit[i] = mucTrk->numHits();
    nlayer[i] = mucTrk->numLayers();
    nhitmax[i] = mucTrk->maxHitsInLayer();
    brlast[i] = mucTrk->brLastLayer();
    eclast[i] = mucTrk->ecLastLayer();
  }
}

