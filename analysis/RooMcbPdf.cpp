#include <iostream>
#include <math.h> 

#include <boost/format.hpp>

using namespace std;

#include "Riostream.h" 

#include "RooMcbPdf.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include "TMath.h" 

//ClassImp(RooMcbPdf)

RooMcbPdf::RooMcbPdf(const char *name, const char *title,
	      RooAbsReal & _X, //invariant mass 
				RooAbsReal & _Sigma,   //Common sigma
				RooAbsReal & _s1,   //staple
				RooAbsReal & _s2,   //staple 
				RooAbsReal & _s3,   //staple 
				RooAbsReal & _s4,   //staple 
				RooAbsReal & _s5,   //staple 
				RooAbsReal & _s6,   //staple 
				RooAbsReal & _s7,   //staple 
				RooAbsReal & _n1,   //staple 
				RooAbsReal & _n2   //staple 
				) : RooAbsPdf(name,title), 
   			fX("x","x",this,_X),
   			fSigma("sigma","sigma",this,_Sigma), 
   			fStaple1("staple1","staple1",this,_s1), 
   			fStaple2("staple2","staple2",this,_s2), 
   			fStaple3("staple3","staple3",this,_s3), 
   			fStaple4("staple4","staple4",this,_s4), 
   			fStaple5("staple5","staple5",this,_s5), 
   			fStaple6("staple6","staple6",this,_s6), 
   			fStaple7("staple7","staple7",this,_s7), 
   			fN1("nl","nl",this,_n1), 
   			fN2("nr","nr",this,_n2) 
{
	fType = TYPE_GEPE;
}

RooMcbPdf::RooMcbPdf(const char *name, const char *title,
	      RooAbsReal & _X, //invariant mass 
				RooAbsReal & _Sigma,   //Common sigma
				RooAbsReal & _s1,   //staple
				RooAbsReal & _s2,   //staple 
				RooAbsReal & _s3,   //staple 
				RooAbsReal & _s4,   //staple 
				RooAbsReal & _s5,   //staple 
				RooAbsReal & _n1,   //staple 
				RooAbsReal & _n2   //staple 
				) : RooAbsPdf(name,title), 
   			fX("x","x",this,_X),
   			fSigma("sigma","sigma",this,_Sigma), 
   			fStaple1("staple1","staple1",this,_s1), 
   			fStaple2("staple2","staple2",this,_s2), 
   			fStaple3("staple3","staple3",this,_s3), 
   			fStaple4("staple4","staple4",this,_s4), 
   			fStaple5("staple5","staple5",this,_s5), 
   			fN1("nl","nl",this,_n1), 
   			fN2("nr","nr",this,_n2) 
{
	fType = TYPE_GPE;
}

	


RooMcbPdf::RooMcbPdf(const RooMcbPdf& other, const char* name) :
  RooAbsPdf(other,name), 
   fX("x",this,other.fX), 
	 fSigma("sigma",this,other.fSigma), 
	 fStaple1("staple1",this,other.fStaple1), 
	 fStaple2("staple2",this,other.fStaple2), 
	 fStaple3("staple3",this,other.fStaple3), 
	 fStaple4("staple4",this,other.fStaple4), 
	 fStaple5("staple5",this,other.fStaple5), 
	 fStaple6("staple6",this,other.fStaple6), 
	 fStaple7("staple7",this,other.fStaple7), 
	 fN1("nl",this,other.fN1), 
	 fN2("nr",this,other.fN2) 
{
	fType = other.fType;
}



Double_t RooMcbPdf::evaluate() const 
{
	initArgs();
	double x = (fX - mean)/fSigma;
	int i = x < 0 ? 0 : 1; //tail index
	double y = TMath::Abs(x); //abs of x
	if(y>c[i]) return C[i]*exp(-(y-c[i])*kc[i]);
	if(y>b[i]) return B[i]*pow(1.0 + kb[i]*(y-b[i]),  -n[i]);
	if(y>a[i]) return A[i]*exp(-(y-a[i])*ka[i]);
	return exp(-0.5*y*y);
}


void RooMcbPdf::initArgs(void) const
{
	switch(fType)
	{
		case TYPE_GPE:
			mean_index = 2;
			staple.resize(5);
			break;
		case TYPE_GEPE:
			mean_index = 3;
			staple.resize(7);
			break;
	};
	//init staple points
	staple[0] = fStaple1;
	staple[1] = fStaple2;
	staple[2] = fStaple3;
	staple[3] = fStaple4;
	staple[4] = fStaple5;
	if(fType == TYPE_GEPE)
	{
		staple[5] = fStaple6;
		staple[6] = fStaple7;
	}
	//sort them
	std::sort(std::begin(staple), std::end(staple));
	mean = staple[mean_index];
	for(auto & s : staple)
	{
		s = (s-mean)/fSigma;
	}
	n[0] = fN1;
	n[1] = fN2;


	for(int i=0;i<2;i++)
	{
		int sign = i ==0 ? -1 : +1;
		a[i] = std::abs(staple[mean_index+1*sign]); 
		switch(fType)
		{
			case TYPE_GPE:
				b[i]=a[i];
				c[i] = std::abs(staple[mean_index+2*sign]);
				break;
			case TYPE_GEPE:
				b[i] = std::abs(staple[mean_index+2*sign]);
				c[i] = std::abs(staple[mean_index+3*sign]);
				break;
		}
	}

	for(int i=0;i<2;i++)
	{
		ka[i] =  a[i];
		kb[i] =  ka[i]/n[i];
		kc[i] =  n[i]*kb[i]/(1.0 + kb[i]*(c[i]-b[i]));
		A[i]  =  TMath::Exp(-0.5*a[i]*a[i]);
		B[i]  =  A[i]*TMath::Exp(- ka[i]*(b[i]-a[i]));
		C[i]  =  B[i]*pow( 1.0  + kb[i]*(c[i]-b[i]),  - n[i]);
		//cout << "A = " << A[i] << endl;
		//cout << "B = " << B[i] << endl;
		//cout << "C = " << C[i] << endl;
		//cout << i << "xb = " << xb[i] << endl;
	}
}


//_____________________________________________________________________________
Int_t RooMcbPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* ) const
{
  if (matchArgs(allVars,analVars,fX)) return 1 ;
  return 0;
}



//_____________________________________________________________________________
Double_t RooMcbPdf::analyticalIntegral(Int_t code, const char* rangeName) const
{
	initArgs();
	//cout << fX.max(rangeName) << " " << fX.min(rangeName) << endl;
	double xmax = (fX.max(rangeName) - mean)/fSigma;
	double xmin = (fX.min(rangeName) - mean)/fSigma;
	double IG[2]; //gauss part of integral
	double IA[2]; //exp part of integral
	double IB[2]; //power part of integral
	double IC[2]; //exp tail part of integral
	double I=0; //sum of previouse one
	double xrange[2]  = { fabs(xmin), fabs(xmax)};
	//std::cout << "xmin = " << xmin << " xmax = " << xmax << endl;
	for(int i=0;i<2;i++)
	{
		//cout << "a["<<i<<"]=" << a[i] << endl;
		//cout << "b["<<i<<"]=" << b[i] << endl;
		//cout << "c["<<i<<"]=" << c[i] << endl;
		//cout << "n["<<i<<"]=" << n[i] << endl;
		if(false) 
		{
			I+=xrange[i];
		}
		else 
		{
			//cout << " A["<< i<<"] = " << A[i] << endl;
			//cout << " B["<<i<<"] = " << B[i] << endl;
			//cout << " C["<<i<<"] = " << C[i] << endl;
			//cout << " ka["<<i<<"] = " << ka[i] << endl;
			//cout << " kc["<<i<<"] = " << kc[i] << endl;
			//cout << " kb["<<i<<"] = " << kb[i] << endl;
			if(xrange[i] < a[i])
			{
				IG[i]  =  sqrt(M_PI*0.5)*TMath::Erf(xrange[i]/sqrt(2.0));
				IA[i]  =  0;
				IB[i]  =  0;
				IC[i]  =  0;
			}
			if(a[i] <= xrange[i] < b[i])
			{
				IG[i]  =  sqrt(M_PI*0.5)*TMath::Erf(a[i]/sqrt(2.0));
				IA[i]  =  A[i]/ka[i]*( 1.0 -  TMath::Exp(- ka[i]*(xrange[i]-a[i])));
				IB[i]  =  0;
				IC[i]  =  0;
			}
			if(b[i] <= xrange[i] < c[i])
			{
				IG[i]  =  sqrt(M_PI*0.5)*TMath::Erf(a[i]/sqrt(2.0));
				//cout << IG[i] << endl;
				IA[i]  =  A[i]/ka[i]*( 1.0 -  TMath::Exp(- ka[i]*(b[i]-a[i])));
				if(n[i]==1.0) IB[i] = B[i]/kb[i]*log(1.0 +kb[i]*(xrange[i]-b[i]));
				else IB[i] = B[i]/kb[i]/(n[i]-1.0)*(1.0  - pow( 1.0 + kb[i]*(xrange[i]-b[i]), - n[i] + 1.));
			}
			if(c[i] <= xrange[i])
			{
				if(ka[i]>0)
				{
					IG[i]  =  sqrt(M_PI*0.5)*TMath::Erf(a[i]/sqrt(2.0));
					IA[i]  =  A[i]/ka[i]*( 1.0 -  TMath::Exp(- ka[i]*(b[i]-a[i])));
					if(n[i]>1.0) IB[i] = B[i]/kb[i]/(n[i]-1.0)*(1.0  - pow( 1.0 + kb[i]*(c[i]-b[i]), - n[i] + 1.));
					else IB[i] = B[i]/kb[i]*log(1.0 +kb[i]*(c[i]-b[i]));
					IC[i]  =  C[i]/kc[i]*( 1.0  - TMath::Exp(- kc[i]*(xrange[i] - c[i])));
				}
				else
				{
					IG[i] = 0;
					IA[i] = A[i]*(b[i]-a[i]);
					IB[i] = B[i]*(c[i]-b[i]);
					IC[i] = C[i]*(xrange[i]-c[i]);
				}
				//cout << " hit max range" << endl;
			}
			else
			{
				cout << "ERROR: wrong integral hit range" << endl;
			}
			//cout << "IG[" << i << "]=" << IG[i] << endl;
			//cout << "IA[" << i << "]=" << IA[i] << endl;
			//cout << "IB[" << i << "]=" << IB[i] << endl;
			//cout << "IC[" << i << "]=" << IC[i] << endl;
			I+=IG[i] + IA[i] + IB[i] + IC[i];
		}
	}
	//cout << "Itotal = " << I << endl;
	return I*fSigma;
}

RooBgPdf::RooBgPdf(const char *name, const char *title,
	      RooAbsReal & _X, //invariant mass 
				RooAbsReal & _B,   //Common sigma
        double x1,
        double x2
				) : RooAbsPdf(name,title), 
   			fX("x","x",this,_X),
   			fB("B","B",this,_B), 
        fXmin(x1),
        fXmax(x2)
{
}

RooBgPdf::RooBgPdf(const RooBgPdf& other, const char* name) :
  RooAbsPdf(other,name), 
   fX("x",this,other.fX), 
	 fB("B",this,other.fB),
   fXmin(other.fXmin),
   fXmax(other.fXmax)
{
}

Int_t RooBgPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* ) const
{
  if (matchArgs(allVars,analVars,fX)) return 1 ;
  return 0;
}



Double_t RooBgPdf::analyticalIntegral(Int_t code, const char* rangeName) const
{
	double xmax = fX.max(rangeName);
	double xmin = fX.min(rangeName);
  return fabs( log ( fabs(xmin - fB) / fabs(xmax - fB) ) );
}

Double_t RooBgPdf::evaluate() const 
{
  return pow ( fabs(fX - fB), -1);
}


RooMcb2Pdf::RooMcb2Pdf(const char *name, const char *title,
	      RooAbsReal & _X, //invariant mass 
				RooAbsReal & _Mean,   //mean
				RooAbsReal & _Sigma,   //Common sigma
				RooAbsReal & _s1,   //staple 
				RooAbsReal & _s2,   //staple 
				RooAbsReal & _s3,   //staple 
				RooAbsReal & _s4,   //staple 
				RooAbsReal & _s5,   //staple 
				RooAbsReal & _s6,   //staple 
				RooAbsReal & _s7,   //staple 
				RooAbsReal & _s8,   //staple 
				RooAbsReal & _n1,   //staple 
				RooAbsReal & _n2,   //staple 
        double xmin,
        double xmax
				) : RooAbsPdf(name,title), 
   			fX("x","x",this,_X),
   			fMean("mean","mean",this,_Mean), 
   			fSigma("sigma","sigma",this,_Sigma), 
   			fStaple1("staple1","staple1",this,_s1), 
   			fStaple2("staple2","staple2",this,_s2), 
   			fStaple3("staple3","staple3",this,_s3), 
   			fStaple4("staple4","staple4",this,_s4), 
   			fStaple5("staple5","staple5",this,_s5), 
   			fStaple6("staple6","staple6",this,_s6), 
   			fStaple7("staple7","staple7",this,_s7), 
   			fStaple8("staple8","staple8",this,_s8), 
   			fN1("nl","nl",this,_n1), 
   			fN2("nr","nr",this,_n2),
        fXmin(xmin),
        fXmax(xmax)
{
}

RooMcb2Pdf::RooMcb2Pdf(const char *name, const char *title,
	      RooAbsReal & _X, //invariant mass 
				RooAbsReal & _Mean,   //mean
				RooAbsReal & _Sigma,   //Common sigma
				std::vector<RooRealVar> & _s,   //staple 
				RooAbsReal & _n1,   
				RooAbsReal & _n2,   
        double xmin,
        double xmax
				) :  
  // 			fX("x","x",this,_X),
  // 			fMean("mean","mean",this,_Mean), 
  // 			fSigma("sigma","sigma",this,_Sigma), 
  // 			fStaple1("staple1","staple1",this,_s[0]), 
  // 			fStaple2("staple2","staple2",this,_s[1]), 
  // 			fStaple3("staple3","staple3",this,_s[2]), 
  // 			fStaple4("staple4","staple4",this,_s[3]), 
  // 			fStaple5("staple5","staple5",this,_s[4]), 
  // 			fStaple6("staple6","staple6",this,_s[5]), 
  // 			fStaple7("staple7","staple7",this,_s[6]), 
  // 			fStaple8("staple8","staple8",this,_s[7]), 
  // 			fN1("nl","nl",this,_n1), 
  // 			fN2("nr","nr",this,_n2),
  //      fXmin(xmin),
  //      fXmax(xmax)
  RooMcb2Pdf(name,title,
          _X,
          _Mean, 
          _Sigma, 
          _s[0],
          _s[1],
          _s[2],
          _s[3],
          _s[4],
          _s[5],
          _s[6],
          _s[7], 
          _n1, 
          _n2, 
          xmin, 
          xmax)
{
}

	
RooMcb2Pdf::RooMcb2Pdf(const RooMcb2Pdf& other, const char* name) :
  RooAbsPdf(other,name), 
   fX("x",this,other.fX), 
	 fMean("mean",this,other.fMean), 
	 fSigma("sigma",this,other.fSigma), 
	 fStaple1("staple1",this,other.fStaple1), 
	 fStaple2("staple2",this,other.fStaple2), 
	 fStaple3("staple3",this,other.fStaple3), 
	 fStaple4("staple4",this,other.fStaple4), 
	 fStaple5("staple5",this,other.fStaple5), 
	 fStaple6("staple6",this,other.fStaple6), 
	 fStaple7("staple7",this,other.fStaple7), 
	 fStaple8("staple8",this,other.fStaple8), 
	 fN1("nl",this,other.fN1), 
	 fN2("nr",this,other.fN2),
   fXmin(other.fXmin),
   fXmax(other.fXmax)
{
}



Double_t RooMcb2Pdf::evaluate() const 
{
	initArgs();
	double x = (fX - fMean)/fSigma;
	int i = x < 0 ? 0 : 1; //tail index 0 for left tail, 1 for right tail
	double y = TMath::Abs(x); //abs of x
	if(y>c[i]) return C[i]*exp(-(y-c[i])*kc[i]);
	if(y>b[i]) return B[i]*pow(1.0 + kb[i]*(y-b[i]),  -n[i]);
	if(y>a[i]) return A[i]*exp(-(y-a[i])*ka[i]);
	return exp(-0.5*y*y);
}


void RooMcb2Pdf::initArgs(void) const
{
	//Init left staple points
  staple[0].resize(3);
  staple[1].resize(3);

  double SL[2] = {fStaple1 + fStaple2 + fStaple3 + fStaple4, fStaple5 + fStaple6 + fStaple7 + fStaple8};

  double Norm[2] = { (fMean - fXmin)/SL[0]/fSigma, (fXmax -fMean)/SL[1]/fSigma  };

  //calculate left staple
  staple[0][0] =  Norm[0]*(fStaple1);
  staple[0][1] =  Norm[0]*(fStaple1 + fStaple2);
  staple[0][2] =  Norm[0]*(fStaple1 + fStaple2 + fStaple3);


  //calculate right staple
  staple[1][0] =  Norm[1]*(fStaple5);
  staple[1][1] =  Norm[1]*(fStaple5 + fStaple6);
  staple[1][2] =  Norm[1]*(fStaple5 + fStaple6 + fStaple7);

	n[0] = fN1;
	n[1] = fN2;

  // [0, a] is gaussian part
  // [a, b] exponential part
  // [b, c] polinomial  part
  // [c, max] exponential part
  for(int i=0;i<2;i++)
  {
    a[i] = staple[i][0];
    b[i] = staple[i][1];
    c[i] = staple[i][2];
  }
  /* 
  std::cout << a[0] << " " << a[1] << std::endl;
  std::cout << b[0] << " " << b[1] << std::endl;
  std::cout << c[0] << " " << c[1] << std::endl;
  */

	for(int i=0;i<2;i++)
	{
		ka[i] =  a[i];
		kb[i] =  ka[i]/n[i];
		kc[i] =  n[i]*kb[i]/(1.0 + kb[i]*(c[i]-b[i]));
		A[i]  =  TMath::Exp(-0.5*a[i]*a[i]);
		B[i]  =  A[i]*TMath::Exp(- ka[i]*(b[i]-a[i]));
		C[i]  =  B[i]*pow( 1.0  + kb[i]*(c[i]-b[i]),  - n[i]);
	}
}


Int_t RooMcb2Pdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* ) const
{
  //std::cout << "in getAnayticalIntegral" << std::endl;
  if (matchArgs(allVars,analVars,fX)) return 1 ;
  return 0;
}


Double_t RooMcb2Pdf::analyticalIntegral(Int_t code, const char* rangeName) const
{
	initArgs();
	//cout << fX.max(rangeName) << " " << fX.min(rangeName) << endl;
	double xmax = (fX.max(rangeName) - fMean)/fSigma;
	double xmin = (fX.min(rangeName) - fMean)/fSigma;
	double IG[2]; //gauss part of integral
	double IA[2]; //exp part of integral
	double IB[2]; //power part of integral
	double IC[2]; //exp tail part of integral
	double I=0; //sum of previouse one
	double xrange[2]  = { fabs(xmin), fabs(xmax)};
	//std::cout << "xmin = " << xmin << " xmax = " << xmax << endl;
	for(int i=0;i<2;i++)
	{
		//cout << "a["<<i<<"]=" << a[i] << endl;
		//cout << "b["<<i<<"]=" << b[i] << endl;
		//cout << "c["<<i<<"]=" << c[i] << endl;
		//cout << "n["<<i<<"]=" << n[i] << endl;
		if(false) 
		{
			I+=xrange[i];
		}
		else 
		{
			//cout << " A["<< i<<"] = " << A[i] << endl;
			//cout << " B["<<i<<"] = " << B[i] << endl;
			//cout << " C["<<i<<"] = " << C[i] << endl;
			//cout << " ka["<<i<<"] = " << ka[i] << endl;
			//cout << " kc["<<i<<"] = " << kc[i] << endl;
			//cout << " kb["<<i<<"] = " << kb[i] << endl;
			if(xrange[i] < a[i])
			{
				IG[i]  =  sqrt(M_PI*0.5)*TMath::Erf(xrange[i]/sqrt(2.0));
				IA[i]  =  0;
				IB[i]  =  0;
				IC[i]  =  0;
			}
			if(a[i] <= xrange[i] < b[i])
			{
				IG[i]  =  sqrt(M_PI*0.5)*TMath::Erf(a[i]/sqrt(2.0));
				IA[i]  =  A[i]/ka[i]*( 1.0 -  TMath::Exp(- ka[i]*(xrange[i]-a[i])));
				IB[i]  =  0;
				IC[i]  =  0;
			}
			if(b[i] <= xrange[i] < c[i])
			{
				IG[i]  =  sqrt(M_PI*0.5)*TMath::Erf(a[i]/sqrt(2.0));
				//cout << IG[i] << endl;
				IA[i]  =  A[i]/ka[i]*( 1.0 -  TMath::Exp(- ka[i]*(b[i]-a[i])));
				if(n[i]==1.0) IB[i] = B[i]/kb[i]*log(1.0 +kb[i]*(xrange[i]-b[i]));
				else IB[i] = B[i]/kb[i]/(n[i]-1.0)*(1.0  - pow( 1.0 + kb[i]*(xrange[i]-b[i]), - n[i] + 1.));
			}
			if(c[i] <= xrange[i])
			{
				if(ka[i]>0)
				{
					IG[i]  =  sqrt(M_PI*0.5)*TMath::Erf(a[i]/sqrt(2.0));
					IA[i]  =  A[i]/ka[i]*( 1.0 -  TMath::Exp(- ka[i]*(b[i]-a[i])));
					if(n[i]>1.0) IB[i] = B[i]/kb[i]/(n[i]-1.0)*(1.0  - pow( 1.0 + kb[i]*(c[i]-b[i]), - n[i] + 1.));
					else IB[i] = B[i]/kb[i]*log(1.0 +kb[i]*(c[i]-b[i]));
					IC[i]  =  C[i]/kc[i]*( 1.0  - TMath::Exp(- kc[i]*(xrange[i] - c[i])));
				}
				else
				{
					IG[i] = 0;
					IA[i] = A[i]*(b[i]-a[i]);
					IB[i] = B[i]*(c[i]-b[i]);
					IC[i] = C[i]*(xrange[i]-c[i]);
				}
				//cout << " hit max range" << endl;
			}
			else
			{
        if(fabs(c[i]-xrange[i]) > 1e-14)
        {

          std::cerr << "WARNING: wrong integral hit range: "<< boost::format("c[%d] (%20.15f) > xrange[%d](%20.15f)") % i % c[i] % i % xrange[i] << " c-xrange = " <<  c[i]-xrange[i] << std::endl;
        }
			}
			//cout << "IG[" << i << "]=" << IG[i] << endl;
			//cout << "IA[" << i << "]=" << IA[i] << endl;
			//cout << "IB[" << i << "]=" << IB[i] << endl;
			//cout << "IC[" << i << "]=" << IC[i] << endl;
			I+=IG[i] + IA[i] + IB[i] + IC[i];
		}
	}
//  std::clog << "rangeName = " << rangeName << "  code = " << code << std::endl;
//  std::clog << " xmin = " << xmin*fSigma+fMean << " xmax = " << xmax*fSigma + fMean   << "  Itotal = " << I*fSigma << std::endl;
	return I*fSigma;
}

