// =====================================================================================
//
//       Filename:  RootMCTopo.cxx
//
//    Description:  
//
//        Version:  1.0
//        Created:  27.10.2015 17:03:52
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#include "RootMCTopo.h"
StatusCode RootMCTopo::init_tuple(void)
{
  StatusCode status;
  status = tuple->addItem("indexmc", m_idxmc, 0, 100);
  status = tuple->addIndexedItem("pdgid", m_idxmc, m_pdgid);
  status = tuple->addIndexedItem("motheridx", m_idxmc, m_motheridx);
  status = tuple->addIndexedItem("idx", m_idxmc, m_idx);
  status = tuple->addItem("hash", m_hash);
  return status;
}

void RootMCTopo::init(void)
{
  m_idxmc=0;
}
