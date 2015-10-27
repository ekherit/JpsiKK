// =====================================================================================
//
//       Filename:  RootTof.cxx
//
//    Description:  
//
//        Version:  1.0
//        Created:  27.10.2015 16:50:38
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#include "RootTof.h"

StatusCode RootTof::init_tuple(void)
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

void RootTof::init(void)
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
