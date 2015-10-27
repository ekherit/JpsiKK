// =====================================================================================
//
//       Filename:  RootAdaptors.h
//
//    Description:  list of root adaptors
//
//        Version:  1.0
//        Created:  27.10.2015 17:14:43
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#pragma once 

#include "RootEvent.h"
#include "RootMdc.h"
#include "RootDedx.h"
#include "RootTof.h"
#include "RootMCTopo.h"
#include "RootMC.h"
#include "RootPid.h"

template <class A>
inline StatusCode init_tuple(JpsiKK * alg, A & a,  const char * dir, const char * title, MsgStream & log)
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
