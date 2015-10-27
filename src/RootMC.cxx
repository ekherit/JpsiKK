// =====================================================================================
//
//       Filename:  RootMC.cxx
//
//    Description:  
//
//        Version:  1.0
//        Created:  27.10.2015 17:03:06
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#include "RootMC.h"

StatusCode RootMC::init_tuple(void)
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

void RootMC::init(void)
{
  ntrack=4;
}
