// =====================================================================================
//
//       Filename:  RootEmc.cxx
//
//    Description:  
//
//        Version:  1.0
//        Created:  27.10.2015 15:32:04
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#include "RootEmc.h"

StatusCode RootEmc::init_tuple(void)
{
  StatusCode status;
  status = tuple->addItem ("ntrack",       ntrack, 0, ARRAY_SIZE); //good nuetral track in event
  status = tuple->addIndexedItem ("E",     ntrack, E);
  status = tuple->addIndexedItem ("theta", ntrack, theta);
  status = tuple->addIndexedItem ("phi",   ntrack, phi);
  status = tuple->addIndexedItem ("time",  ntrack, time);
  return status;
}

void RootEmc::init(void)
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
