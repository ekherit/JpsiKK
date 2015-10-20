// =====================================================================================
//
//       Filename:  utils.h
//
//    Description:  math utils
//
//        Version:  1.0
//        Created:  20.10.2015 11:35:28
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================
#pragma once

inline double sq(double x)
{
	return x*x;
}

inline bool in(double x,  double a,  double b)
{
	return a < x && x < b;
}

inline bool out(double x,  double a,  double b)
{
	return !in(x, a, b);
}
