// =====================================================================================
//
//       Filename:  Defs.h
//
//    Description:  
//
//        Version:  1.0
//        Created:  20.10.2015 11:29:50
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================
#pragma once

#include <list>
#include <vector>
#include <map>

#include "EvtRecEvent/EvtRecTrack.h"

typedef std::pair<EvtRecTrackIterator, EvtRecTrackIterator> TrackPair_t;
typedef std::list<TrackPair_t> TrackPairList_t;
typedef std::list<EvtRecTrackIterator> TrackList_t;
typedef std::vector<EvtRecTrackIterator> TrackVector_t;

inline TrackVector_t make_track_vector(TrackPair_t & pair1, TrackPair_t & pair2)
{
  TrackVector_t V(4);
  V[0] = pair1.first;
  V[1] = pair1.second;
  V[2] = pair2.first;
  V[3] = pair2.second;
  return V;
}

inline TrackVector_t make_track_vector(TrackPair_t & pair1)
{
  TrackVector_t V(2);
  V[0] = pair1.first;
  V[1] = pair1.second;
  return V;
}
