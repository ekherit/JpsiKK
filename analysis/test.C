// =====================================================================================
//
//       Filename:  test.C
//
//    Description:  
//
//        Version:  1.0
//        Created:  21.06.2015 23:17:22
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#include "mctopo/libMcTopo.h"
#include "CrystalBall.h"

#include "analize.h"
void test(void)
{
	gROOT->Reset();
  gSystem->Load("mctopo/libMcTopo.so");
  gSystem->AddIncludePath("-I$HOME/work -I./mctop");
	gSystem->CompileMacro("CrystalBall.cpp","kO","","/tmp");
	gSystem->CompileMacro("analize.C","kO","","/tmp");
	gSystem->CompileMacro("load.C","kO","","/tmp");
}

