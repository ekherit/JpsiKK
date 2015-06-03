// =====================================================================================
//
//       Filename:  CrystalBall.h
//
//    Description:   Modified Double Crystal Ball function
//
//        Version:  1.0
//        Created:  10.03.2015 16:07:29
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#ifndef IBN_CRYSTAL_BALL_H
#define IBN_CRYSTAL_BALL_H

#include <vector>

// ===  FUNCTION  ======================================================================
//         
//         Name:  ModifiedDoubleCrystalBall
//
//  Description:   Parameters:
//  P0  - number of signal event
//  P1  - mean of the gauss
//  P2  - dispersion of the gaus
//  P3  - left shift of the exp tait
//  P4  - left shift of the power tail relative to exp shift (P3)
//  P5  - power of power tail
//  P6  - right shift of the exp tait
//  P7  - right shift of the power tail relative to exp shift (P6)
//  P8  - power of right power tail
//  P9  - number of background event
//  P10 - slope of the background event
//  P11 - minimum fit range
//  P12 - maxinum fit range
//  P13 - number of bins in histogram
// =====================================================================================

#include <TH1F.h>

extern double ModifiedDoubleCrystalBall(const double* X, const double* P);

class CrystalBallFitter2;


extern const double * Fit(TH1F * his);
extern std::vector<double> Fit2(TH1F * his);

#endif
