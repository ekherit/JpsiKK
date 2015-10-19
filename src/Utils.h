// =====================================================================================
//
//       Filename:  Utils.h
//
//    Description:  Some usefull functions
//
//        Version:  1.0
//        Created:  19.10.2015 21:12:17
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

inline bool in(double x,  double a,  double b)
{
	return a < x && x < b;
}

inline bool out(double x,  double a,  double b)
{
	return !in(x, a, b);
}
