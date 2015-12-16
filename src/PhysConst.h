// =====================================================================================
//
//       Filename:  PhysConst.h
//
//    Description:  Physics Constat used during selection
//
//        Version:  1.0
//        Created:  19.10.2015 21:06:25
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================
#pragma once

const double PION_MASS     = 0.13957018; //GeV
const double MUON_MASS     = 0.105658389; //GeV
const double KAON_MASS     = 0.493677; //GeV
const double ELECTRON_MASS = 0.000510999;//GeV
const double PROTON_MASS   = 0.93827231;//GeV

const double PI0_MASS      = 0.1349766; //GeV

const double JPSI_MASS     = 3.096916; //GeV
const double PSIP_MASS     = 3.686109; //GeV PDG-2014

enum  {ID_KAON=0, ID_MUON=1, ID_ELECTRON=2, ID_PION=3, ID_PROTON=4};

extern const double XMASS[5];


const double BEPC_CROSSING_ANGLE=0.022;
extern double BEAM_CENTER_MASS_ENERGY;
