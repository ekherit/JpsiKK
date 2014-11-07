#ifndef IBN_INTEGRAL_H
#define IBN_INTEGRAL_H

#include <complex>
#include <algorithm>
#include <cmath>
using namespace std;
/* 
   Здесь собраны различные алгоритмы для вычисления интеграла 
   */


#include <iostream>

namespace ibn 
{


  /* интегралл прямоугольниками */
  template <class F> double integrate0(F  f, double min, double max, double step)
  {
    double sum=0;
    for(double x = min; x < max; x+=step)	{
      sum+= f(x);
    }
    return sum*step;
  }

  // интегралл трапециями 
  template <class F> double integrate1(F  f, double min, double max, double step)
  {
    double sum=0;
    for(double x = min; x < max; x+=step)
    {
      sum+= f(x) + f(x+step);
    }
    return sum*step/2.;
  }

  // интеграл параболами
  template <class F> double integrate2(F  f, double min, double max, double step)
  {
    double sum=0;
    for(double x = min; x < max; x+=step)
    {
      sum+= f(x) + f(x+step) + 4.*f(x+step/2.) ;
    }
    return sum*step/6.;
  }


  // интеграл с помощью полинома третьей степени
  template <class F> double integrate3(F  f, double min, double max, double step)
  {
    double sum=0;
    for(double x = min; x < max; x+=step)
    {
      sum+= f(x) + f(x+step) + 3.* ( f(x+step/3.)  + f(x+2./3.*step) ) ;
    }
    return sum*step/8.;
  }


  //алгоритм из cernlib
  //интеграл 16-ти и 8 точеченой гауссовой квадратурой.
  template <class F> double dgaus(F  func, double a, double b, double epsilon)
  {
    const double Z1 = 1;
    const double HF = Z1/2;
    const double CST = 5*Z1/1000;

    double x[12] = { 0.96028985649753623,  0.79666647741362674,
      0.52553240991632899,  0.18343464249564980,
      0.98940093499164993,  0.94457502307323258,
      0.86563120238783174,  0.75540440835500303,
      0.61787624440264375,  0.45801677765722739,
      0.28160355077925891,  0.09501250983763744
    };

    double w[12] = { 0.10122853629037626,  0.22238103445337447,
      0.31370664587788729,  0.36268378337836198,
      0.02715245941175409,  0.06225352393864789,
      0.09515851168249278,  0.12462897125553387,
      0.14959598881657673,  0.16915651939500254,
      0.18260341504492359,  0.18945061045506850
    };

    register double h, aconst, bb, aa, c1, c2, u, s8, s16, f1, f2;
    double xx[1];
    int i;

    //InitArgs(xx,params);

    h = 0;
    if (b == a) return h;
    aconst = CST/fabs(b-a);
    bb = a;
CASE1:
    aa = bb;
    bb = b;
CASE2:
    c1 = HF*(bb+aa);
    c2 = HF*(bb-aa);
    s8 = 0;
    for (i=0;i<4;i++) {
      u     = c2*x[i];
      xx[0] = c1+u;
      f1    = func(xx[0]);
      xx[0] = c1-u;
      f2    = func(xx[0]);
      s8   += w[i]*(f1 + f2);
    }
    s16 = 0;
    for (i=4;i<12;i++) {
      u     = c2*x[i];
      xx[0] = c1+u;
      f1    = func(xx[0]);
      xx[0] = c1-u;
      f2    = func(xx[0]);
      s16  += w[i]*(f1 + f2);
    }
    s16 = c2*s16;
    //	if (fabs(s16-c2*s8) <= epsilon*(fabs(s16))) {
    if (fabs(s16-c2*s8) <= epsilon*(1. + fabs(s16))) {
      h += s16;
      if(bb != b) goto CASE1;
    } else {
      bb = c1;
      if(1. + aconst*fabs(c2) != 1) goto CASE2;
      h = s8;  //this is a crude approximation (cernlib function returned 0 !)
    }
    return h;
  }


  //Многомерный интеграл
  template <class F> double dgaus(F  func, double *a, double *b, double eps, double &relerr)
  {
    double ctr[15], wth[15], wthl[15], z[15];
    int n = func.dimension();
    const double xl2 = 0.358568582800318073;
    const double xl4 = 0.948683298050513796;
    const double xl5 = 0.688247201611685289;
    const double w2  = 980./6561;
    const double w4  = 200./19683;
    const double wp2 = 245./486;
    const double wp4 = 25./729;

    double wn1[14] = {     -0.193872885230909911, -0.555606360818980835,
      -0.876695625666819078, -1.15714067977442459,  -1.39694152314179743,
      -1.59609815576893754,  -1.75461057765584494,  -1.87247878880251983,
      -1.94970278920896201,  -1.98628257887517146,  -1.98221815780114818,
      -1.93750952598689219,  -1.85215668343240347,  -1.72615963013768225
    };

    double wn3[14] = {     0.0518213686937966768,  0.0314992633236803330,
      0.0111771579535639891,-0.00914494741655235473,-0.0294670527866686986,
      -0.0497891581567850424,-0.0701112635269013768, -0.0904333688970177241,
      -0.110755474267134071, -0.131077579637250419,  -0.151399685007366752,
      -0.171721790377483099, -0.192043895747599447,  -0.212366001117715794
    };

    double wn5[14] = {         0.871183254585174982e-01,  0.435591627292587508e-01,
      0.217795813646293754e-01,  0.108897906823146873e-01,  0.544489534115734364e-02,
      0.272244767057867193e-02,  0.136122383528933596e-02,  0.680611917644667955e-03,
      0.340305958822333977e-03,  0.170152979411166995e-03,  0.850764897055834977e-04,
      0.425382448527917472e-04,  0.212691224263958736e-04,  0.106345612131979372e-04
    };

    double wpn1[14] = {   -1.33196159122085045, -2.29218106995884763,
      -3.11522633744855959, -3.80109739368998611, -4.34979423868312742,
      -4.76131687242798352, -5.03566529492455417, -5.17283950617283939,
      -5.17283950617283939, -5.03566529492455417, -4.76131687242798352,
      -4.34979423868312742, -3.80109739368998611, -3.11522633744855959
    };

    double wpn3[14] = {     0.0445816186556927292, -0.0240054869684499309,
      -0.0925925925925925875, -0.161179698216735251,  -0.229766803840877915,
      -0.298353909465020564,  -0.366941015089163228,  -0.435528120713305891,
      -0.504115226337448555,  -0.572702331961591218,  -0.641289437585733882,
      -0.709876543209876532,  -0.778463648834019195,  -0.847050754458161859
    };

    double result = 0;
    double abserr = 0;
    int ifail = 3;
    if (n < 2 || n > 15) return 0;

    double twondm = std::pow(2.0,n);
    int ifncls = 0;
    bool ldv   = false;
    int irgnst = 2*n+3;
    int irlcls = int(twondm) +2*n*(n+1)+1;
    int isbrgn = irgnst;
    int isbrgs = irgnst;

    // The original algorithm expected a parameter MAXPTS
    //   where MAXPTS = Maximum number of function evaluations to be allowed.
    //   Here we set MAXPTS to 1000*(the lowest possible value)
    int maxpts = 1000*irlcls;
    int minpts = 1;

    // The original agorithm expected a working space array WK of length IWK
    // with IWK Length ( >= (2N + 3) * (1 + MAXPTS/(2**N + 2N(N + 1) + 1))/2).
    // Here, this array is allocated dynamically

    int iwk = irgnst*(1 +maxpts/irlcls)/2;
    double *wk = new double[iwk+10];
    int j;
    for (j=0;j<n;j++) {
      ctr[j] = (b[j] + a[j])*0.5;
      wth[j] = (b[j] - a[j])*0.5;
    }

    double rgnvol, sum1, sum2, sum3, sum4, sum5, difmax, f2, f3, dif;
    double rgncmp=0, rgnval, rgnerr;
    int j1, k, l, m, idvaxn=0, idvax0=0, isbtmp, isbtpp;

    //InitArgs(z,fParams);

L20:
    rgnvol = twondm;
    for (j=0;j<n;j++) {
      rgnvol *= wth[j];
      z[j]    = ctr[j];
    }
    sum1 = func(z); //evaluate function

    difmax = 0;
    sum2   = 0;
    sum3   = 0;
    for (j=0;j<n;j++) {
      z[j]    = ctr[j] - xl2*wth[j];
      f2      = func(z);
      z[j]    = ctr[j] + xl2*wth[j];
      f2     += func(z);
      wthl[j] = xl4*wth[j];
      z[j]    = ctr[j] - wthl[j];
      f3      = func(z);
      z[j]    = ctr[j] + wthl[j];
      f3     += func(z);
      sum2   += f2;
      sum3   += f3;
      dif     = fabs(7*f2-f3-12*sum1);
      difmax  = max(dif, difmax);
      if (difmax == dif) idvaxn = j+1;
      z[j]    = ctr[j];
    }

    sum4 = 0;
    for (j=1;j<n;j++) {
      j1 = j-1;
      for (k=j;k<n;k++) {
        for (l=0;l<2;l++) {
          wthl[j1] = -wthl[j1];
          z[j1]    = ctr[j1] + wthl[j1];
          for (m=0;m<2;m++) {
            wthl[k] = -wthl[k];
            z[k]    = ctr[k] + wthl[k];
            sum4 += func(z);
          }
        }
        z[k] = ctr[k];
      }
      z[j1] = ctr[j1];
    }

    sum5 = 0;
    for (j=0;j<n;j++) {
      wthl[j] = -xl5*wth[j];
      z[j] = ctr[j] + wthl[j];
    }
L90:
    sum5 += func(z);
    for (j=0;j<n;j++) {
      wthl[j] = -wthl[j];
      z[j] = ctr[j] + wthl[j];
      if (wthl[j] > 0) goto L90;
    }

    rgncmp  = rgnvol*(wpn1[n-2]*sum1+wp2*sum2+wpn3[n-2]*sum3+wp4*sum4);
    rgnval  = wn1[n-2]*sum1+w2*sum2+wn3[n-2]*sum3+w4*sum4+wn5[n-2]*sum5;
    rgnval *= rgnvol;
    rgnerr  = fabs(rgnval-rgncmp);
    result += rgnval;
    abserr += rgnerr;
    ifncls += irlcls;

    if (ldv)
    {
L110:
      isbtmp = 2*isbrgn;
      if (isbtmp > isbrgs) goto L160;
      if (isbtmp < isbrgs)
      {
        isbtpp = isbtmp + irgnst;
        if (wk[isbtmp-1] < wk[isbtpp-1]) isbtmp = isbtpp;
      }
      if (rgnerr >= wk[isbtmp-1]) goto L160;
      for (k=0;k<irgnst;k++) {
        wk[isbrgn-k-1] = wk[isbtmp-k-1];
      }
      isbrgn = isbtmp;
      goto L110;
    }
L140:
    isbtmp = (isbrgn/(2*irgnst))*irgnst;
    if (isbtmp >= irgnst && rgnerr > wk[isbtmp-1])
    {
      for (k=0;k<irgnst;k++)
      {
        wk[isbrgn-k-1] = wk[isbtmp-k-1];
      }
      isbrgn = isbtmp;
      goto L140;
    }

L160:
    wk[isbrgn-1] = rgnerr;
    wk[isbrgn-2] = rgnval;
    wk[isbrgn-3] = double(idvaxn);
    for (j=0;j<n;j++)
    {
      isbtmp = isbrgn-2*j-4;
      wk[isbtmp]   = ctr[j];
      wk[isbtmp-1] = wth[j];
    }
    if (ldv)
    {
      ldv = false;
      ctr[idvax0-1] += 2*wth[idvax0-1];
      isbrgs += irgnst;
      isbrgn  = isbrgs;
      goto L20;
    }
    relerr = abserr/fabs(result);
    if (isbrgs+irgnst > iwk) ifail = 2;
    if (ifncls+2*irlcls > maxpts) ifail = 1;
    if (relerr < eps && ifncls >= minpts) ifail = 0;
    if (ifail == 3)
    {
      ldv = true;
      isbrgn  = irgnst;
      abserr -= wk[isbrgn-1];
      result -= wk[isbrgn-2];
      idvax0  = int(wk[isbrgn-3]);
      for (j=0;j<n;j++) {
        isbtmp = isbrgn-2*j-4;
        ctr[j] = wk[isbtmp];
        wth[j] = wk[isbtmp-1];
      }
      wth[idvax0-1]  = 0.5*wth[idvax0-1];
      ctr[idvax0-1] -= wth[idvax0-1];
      goto L20;
    }
    // IFAIL On exit:
    //     0 Normal exit.  . At most MAXPTS calls to the function F were performed.
    //     1 MAXPTS is too small for the specified accuracy EPS. RESULT and RELERR
    //              contain the values obtainable for the specified value of MAXPTS.
    //
    delete [] wk;
    //   int nfnevl = ifncls; //number of function evaluations performed.
    //is_nan(result);
    return result;         //an approximate value of the integral
  }









  template <class F1, class F2> double svertka(F1 &f1, F2 &f2, double S, double min, double max, double epsilon)
  {
    class F
    {
      F1 & ff1;
      F2 & ff2;
      const double SS;
      public:
      F(F1 & ff1_, F2 & ff2_, double SS_) : ff1(ff1_), ff2(ff2_), SS(SS_) {}
      double operator()(double s)
      {
        return ff1(s)*ff2(SS - s);
      }
    };
    F fun(f1,f2, S);
    return dgaus(fun, min, max, epsilon);
  }

  template <class F1, class F2, class Z> double svertka_reg(F1 &f1, F2 &f2, Z &z, double S, double min, double max, double epsilon)	{
    class F 	{
      F1 & ff1;
      F2 & ff2;
      Z & zz;
      const double SS;
      public:
      F(F1 & ff1_, F2 & ff2_, Z & zz_ , double SS_) : ff1(ff1_), ff2(ff2_), zz(zz_), SS(SS_) {}
      double operator()(double s)
      {
        return ff1(zz(s))*ff2(SS - zz(s));
      }
    };
    F fun(f1,f2, S);
    return dgaus(fun, min, max, epsilon);
  }



  template <class F> complex<double> dgaus_comp(F  func, double a, double b, double epsilon)
  {
    const double Z1 = 1;
    const double HF = Z1/2;
    const double CST = 5*Z1/1000;

    double x[12] = { 0.96028985649753623,  0.79666647741362674,
      0.52553240991632899,  0.18343464249564980,
      0.98940093499164993,  0.94457502307323258,
      0.86563120238783174,  0.75540440835500303,
      0.61787624440264375,  0.45801677765722739,
      0.28160355077925891,  0.09501250983763744
    };

    double w[12] = { 0.10122853629037626,  0.22238103445337447,
      0.31370664587788729,  0.36268378337836198,
      0.02715245941175409,  0.06225352393864789,
      0.09515851168249278,  0.12462897125553387,
      0.14959598881657673,  0.16915651939500254,
      0.18260341504492359,  0.18945061045506850
    };

    double  aconst, bb, aa, c1, c2, u;
    complex<double> h,s8,s16,f1,f2;
    double xx[1];
    int i;

    //InitArgs(xx,params);

    h = 0;
    if (b == a) return h;
    aconst = CST/fabs(b-a);
    bb = a;
CASE1:
    aa = bb;
    bb = b;
CASE2:
    c1 = HF*(bb+aa);
    c2 = HF*(bb-aa);
    s8 = 0;
    for (i=0;i<4;i++)
    {
      u     = c2*x[i];
      xx[0] = c1+u;
      f1    = func(xx[0]);
      xx[0] = c1-u;
      f2    = func(xx[0]);
      s8   += w[i]*(f1 + f2);
    }
    s16 = 0;
    for (i=4;i<12;i++)
    {
      u     = c2*x[i];
      xx[0] = c1+u;
      f1    = func(xx[0]);
      xx[0] = c1-u;
      f2    = func(xx[0]);
      s16  += w[i]*(f1 + f2);
    }
    s16 = c2*s16;
    if (abs(s16-c2*s8) <= epsilon*(1. + abs(s16)))
    {
      h += s16;
      if(bb != b) goto CASE1;
    } else
    {
      bb = c1;
      if(1. + aconst*fabs(c2) != 1) goto CASE2;
      h = s8;  //this is a crude approximation (cernlib function returned 0 !)
    }
    return h;
  }

  template <class F> double derivative(F f, double x, double epsilon)
  {
    double df1,df2;
    if(epsilon == 0) return 0;
    do
    {
      epsilon/=2.;
      df1 = f(x+2.*epsilon) - f(x);
      df2 = f(x+epsilon) - f(x);
    } while ( fabs(df2 - df1) > epsilon*(1.+fabs(df2)));	
  }

  template <class F> 
    class Integrate
    {
      const int method;
      public:
      Integrate(int n): method(n) {}
      double operator()(F f,double min, double max, double step)
      {
        switch(method)
        {
          case 0:
            return integrate0(f,min, max,step);
            break;
          case 1:
            return integrate1(f,min,max,step);
            break;
          case 2:
            return integrate2(f,min,max,step);
            break;
          case 3:
            return integrate3(f,min,max,step);	
            break;
        }
      }
    };
  }
#endif
