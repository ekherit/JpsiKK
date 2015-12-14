#pragma once

#include "RootTuple.h"
#include "RootMass.h"
#include "RootTrack.h"

struct RootMuc : public RootTuple
{
  NTuple::Item<long>    valid; //number of selected tracks. Some of them dont have muc information
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
  virtual void init(void);
  virtual void init_tuple(void);
  virtual void fill(int, EvtRecTrack *);
};
