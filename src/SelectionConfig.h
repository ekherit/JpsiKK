// =====================================================================================
//
//       Filename:  SelectionConfig.h
//
//    Description:  
//
//        Version:  1.0
//        Created:  20.10.2015 12:26:35
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#pragma once

struct SelectionConfig
{
	double CENTER_MASS_ENERGY;      //center mass energy

	int MIN_CHARGED_TRACKS; //minimum good charged tracks in selection
	int MAX_CHARGED_TRACKS; //maximum good charged tracks in selection
	int MAX_NEUTRAL_TRACKS; //maximum good neutral tracks in selection
	double IP_MAX_RHO, IP_MAX_Z; //interection point cut
  double MAX_COS_THETA; //maximum  cos(theta) for good charged track

  double EMC_ENDCUP_MIN_COS_THETA;
  double EMC_ENDCUP_MAX_COS_THETA;
  double EMC_ENDCUP_MIN_ENERGY;
  double EMC_BARREL_MAX_COS_THETA;
  double EMC_BARREL_MIN_ENERGY;

  double MAX_MUON_EP_RATIO;
  double MAX_KAON_EP_RATIO;

  double MIN_MUON_EP_RATIO;
  double MIN_KAON_EP_RATIO;

  double MAX_PION_MOMENTUM; //maximum pion momentum
  double MIN_PION_MOMENTUM; //maximum pion momentum

  double MIN_RECOIL_MASS; //minimum recoil mass cut
  double MAX_RECOIL_MASS; //minimum recoil mass cut

  double MIN_KAON_MOMENTUM; //minimum kaon momentum
  double MAX_KAON_MOMENTUM; //maximum pion momentum

  double MIN_MUON_MOMENTUM; //minimum kaon momentum
  double MAX_MUON_MOMENTUM; //maximum pion momentum

  double MIN_INVARIANT_MASS; //minimum invariant  mass cut
  double MAX_INVARIANT_MASS; //manimum invariant  mass cut

  double MIN_KAON_MISSING_MASS;   //minimum kaon missing mass
  double MAX_KAON_MISSING_MASS;   //minimum kaon missing mass
  double MIN_MUON_MISSING_MASS;   //minimum muon missing mass
  double MAX_MUON_MISSING_MASS;   //minimum muon missing mass

  double MIN_MISSING_MASS; 
  double MAX_MISSING_MASS; 

	double MAX_KIN_CHI2; //maximum chi2 for kinematic fit
	double MAX_PID_CHI2; //maximum chi2 for my PID
};
