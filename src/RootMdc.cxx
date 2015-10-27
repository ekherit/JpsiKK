// =====================================================================================
//
//       Filename:  RootMdc.cxx
//
//    Description:  
//
//        Version:  1.0
//        Created:  27.10.2015 17:02:36
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#include "RootMdc.h"

StatusCode RootMdc::init_tuple(void)
{
  StatusCode status;
	status = M.add_to_tuple(tuple);
	status = T.add_to_tuple(tuple);
  return status;
}


void RootMdc::init(void)
{
  T.ntrack=4;
}
