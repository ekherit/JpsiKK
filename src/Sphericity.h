/*
 * =====================================================================================
 *
 *       Filename:  Sphericity.h
 *
 *    Description:  Calculate Sphericity
 *
 *        Version:  1.0
 *        Created:  05.03.2012 17:27:58
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physics, Novosibirsk, Russia
 *
 * =====================================================================================
 */

#ifndef IBN_SPHERICITY_H
#define IBN_SPHERICITY_H

#include <TMatrixD.h>
#include <TMatrixDEigen.h>
/*  this class is used for calculating sphericity */

class Sphericity
{
  double sum2; //sum of squared
  public:
  TMatrixD S; //sphericity tensor
    Sphericity(void):  S(TMatrixD(3,3))
    {
      //clear 
      for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
          S[i][j]=0;
      sum2=0;
    };
    void add(double x, double y, double z)
    {
      std::vector <double> p(3);
      p[0] = x;
      p[1] = y;
      p[2] = z;
      for(int k=0;k<3;k++)
        for(int m=0;m<3;m++)
          S[k][m]+=p[k]*p[m];
      sum2+=x*x+y*y+z*z;
    }

    void add(const Hep3Vector & p)
    {
      for(int k=0;k<3;k++)
        for(int m=0;m<3;m++)
        {
          S[k][m]+=p[k]*p[m];
        }
      sum2+=p.mag2();
    }

    void norm(void)
    {
      for(int k=0;k<3;k++)
        for(int m=0;m<3;m++)
          S[k][m]/=sum2;
    }

    double operator()(void)
    {
      TMatrixDEigen Stmp(S);
      const TVectorD & eval = Stmp.GetEigenValuesRe();
      std::vector<double> v(3);
      for(int i=0;i<3;i++) v[i]=eval[i];
      std::sort(v.begin(), v.end());
      if(!(v[0]<=v[1] && v[1]<=v[2])) //test the order of eigenvalues
      {
        cerr << "Bad sphericity" << endl;
        exit(1);
      }
      double sphericity = 1.5*(v[0]+v[1]);
      return sphericity;
    }
};
double Sphericity2(TMatrixD & S)
{
  TMatrixDEigen Stmp(S);
  const TVectorD & eval = Stmp.GetEigenValuesRe();
  std::vector<double> v(3);
  for(int i=0;i<3;i++) v[i]=eval[i];
  std::sort(v.begin(), v.end());
  if(!(v[0]<=v[1] && v[1]<=v[2])) //test the order of eigenvalues
  {
    cerr << "Bad sphericity" << endl;
    exit(1);
  }
  double sphericity = 1.5*(v[0]+v[1]);
  return sphericity;
}
#endif
