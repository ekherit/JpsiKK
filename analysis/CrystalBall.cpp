// =====================================================================================
//
//       Filename:  CrystalBall.cpp
//
//    Description:   Implementation Modified Double Crystal Ball function
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


#include <TMath.h>
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include "CrystalBall.h"


#include <memory>
#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>
#include <list>
#include <algorithm>

#include <ibn/integral.h>


#include <cmath>
using namespace std;


#include <TVirtualFitter.h>

inline double sq(double x) {return x*x;}
inline double cb(double x) {return x*x*x;}
inline double qd(double x) {return sq(sq(x));}

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
//

double ModifiedDoubleCrystalBall(const double* X, const double* P)
{
	double x    = X[0];
	double Nsig  = P[0];
	double mu    = P[1];//mean of the gaus
	double sigma = P[2]; //sigma of the gaus


	//tails 0 - left tail,  1 - right tail
	//double beta[2] = {TMath::Abs(P[3]),           TMath::Abs(P[6])};
	//double alfa[2] = {TMath::Abs(P[4]) + beta[0], TMath::Abs(P[7]) + beta[1]};
	//double n[2] =    {TMath::Abs(P[5]),           TMath::Abs(P[8])};
	double beta[2] = {P[3],           P[6]};
	double alfa[2] = {P[4] + beta[0], P[7] + beta[1]};
	double n[2] =    {P[5],P[8]};

	//background
	double NBG = P[9] - Nsig;
	double pbg = P[10]; //background slope in sigma units

	double xmin = P[11];
	double xmax = P[12];
	double Nbins = P[13];

	//some constants
	double C[2]; //correspond exp tail
	double A[2]; //correspond stepennoi tail
	double B[2]; //correspon stepennoi tail handy coef
	for(int i=0;i<2;i++)
	{
		C[i] = TMath::Exp(0.5*beta[i]*beta[i]);
		A[i] = C[i]*pow(n[i]/beta[i], n[i]) * TMath::Exp(- beta[i]*alfa[i]);
		B[i] = n[i]/beta[i] - alfa[i];
	}

	//scale and shift the x
	x = (x-mu)/sigma;

	double xrange[2] = {-(xmin-mu)/sigma, (xmax-mu)/sigma };

	int t = x < 0 ? 0 : 1; //tail index
	double y = TMath::Abs(x);
	double cb=0; //crystal ball function result gaus is 1

	if(y >= 0  && y <= beta[t]) cb = TMath::Exp( - 0.5*y*y);
	if(y > beta[t] && y <= alfa[t]) cb = C[t]*TMath::Exp(- beta[t] * y);
	if(y > alfa[t]) cb = A[t]*pow(B[t] + y,  - n[t]);

	double Igaus[2]; //gauss part of integral
	double Iexp[2]; //exp part of integral
	double Ipow[2]; //power part of integral
	double I=0; //sum of previouse one

	for(int i=0;i<2;i++)
	{
		//calculate partial integral
		if( xrange[i] <= beta[i] ) //only gause in range
		{
			Igaus[i] = sqrt(TMath::Pi()*0.5)*TMath::Erf(xrange[i]/sqrt(2.0));
			Iexp[i] = 0;
			Ipow[i] = 0;
		}
		if( beta[i] <  xrange[i] && xrange[i] <= alfa[i]) //gaus and part of exp tail in range
		{
			Igaus[i] = sqrt(TMath::Pi()*0.5)*TMath::Erf(beta[i]/sqrt(2.0));
			Iexp[i] =  C[i]/beta[i]*( TMath::Exp(-beta[i]*beta[i]) -  TMath::Exp(-beta[i]*xrange[i]) );
			Ipow[i] =0;
		}
		if(alfa[i] < xrange[i]) //all part of function inside the range
		{
			Igaus[i] = sqrt(TMath::Pi()*0.5)*TMath::Erf(beta[i]/sqrt(2.0));
			Iexp[i] =  C[i]/beta[i]*( TMath::Exp(-beta[i]*beta[i]) -  TMath::Exp(-beta[i]*alfa[i]));
			if(n[i]==1.0) Ipow[i] = A[i]*log((B[i] + xrange[i])/(B[i] + alfa[i]));
			else Ipow[i] = A[i]/(n[i]-1.0)*(pow( B[i] + alfa[i], - n[i] + 1.)  - pow( B[i] + xrange[i], - n[i] + 1.));
		}
		//calculate total integral
		I+=Igaus[i] + Iexp[i] +Ipow[i];
	}
	//calculate the scale
	double dx = (xrange[0] + xrange[1])/Nbins;

	//the amplitude of signal
	double Ns = Nsig/I*dx;

	//the amplitude of background
	double Nbg =  NBG*dx/(xrange[0]+xrange[1])/(1.0 + 0.5*pbg*(xrange[1]-xrange[0]));

	//result of my double crystal ball function
	double result = Ns*cb + Nbg*(1.0 + pbg*x);

	//test the result to suppress bad calculation
	//if(TMath::IsNaN(result) || !TMath::Finite(result)) return std::numeric_limits<double>::max();
	if(TMath::IsNaN(result) || result <= 0) return 0;
	if(!TMath::Finite(result)) return 1e100;
	cout << Nsig << " result = " << result << endl;
	return result;
}


#include <Minuit2/MnUserParameters.h>
#include <Minuit2/FCNBase.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnScan.h>
#include <Minuit2/MnApplication.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/FunctionMinimum.h>

using namespace ROOT::Minuit2;

using namespace ROOT::Minuit2;

namespace ibn
{
	class regulator
	{
		double maximum=1e100;
		vector <double> par;
		bool hit=false;
		bool bad(double v) const 
		{
			return std::isnan(v) ||  std::isinf(v);
		}

		bool good(double v) const 
		{
			return !std::isnan(v)  && !std::isinf(v);
		}

		public:
		regulator(void)
		{
		}

		regulator(double fmax)
		{
			maximum = fmax;
		}

		regulator & operator=(const regulator &r )
		{
			maximum = r.maximum;
			hit = r.hit;
			par  = r.par;
			return *this;
		}

		double operator()(double x,  const vector<double> & p)
		{
			//cout << "chi2 = " << setw(10) << x;
			//for(auto pp : p) cout << setw(12) << pp;
			//cout << endl;
			if(good(x) && x<=maximum) 
			{
			//	cout << "good x " << endl;
				return x;
			}
			if(hit)
			{
				//double old = x;
				x=maximum;
				for(unsigned i=0;i<par.size();i++)
				{
					if(good(p[i])) x+= pow(p[i]-par[i], 2.0);
				}
			//	cout << "raw chi2 = " << old << " regular chi2 = "<< x << endl;
				return x;
			}
			else
			{
				if(good(x) && x>maximum)
				{
			//		cout << "save hit chi2 " << x << endl;
					hit=true;
					par=p;
					//maximum = x;
				}
			//	cout << "no hit: return maximum" << endl;
				return maximum;
			}
		}
	};
};

class CrystalBallFitter  : public ROOT::Minuit2::FCNBase
{ 
	TH1F * his;
	TF1 * fun;
	mutable double xmin; //minimum fit range
	mutable double xmax; //maximum fit range
	mutable double Nbins; //number of bins
	mutable double N0; //total number of events in his
	bool debug;
  bool opt_integrate;

	ROOT::Minuit2::MnUserParameters inipar;
	ROOT::Minuit2::MnUserParameters minpar;

	vector<double> fit_result;
	vector< std::pair<double, double> > par_error;

	list<TGraph *> graph_list;
	

	mutable ibn::regulator reg;


	public: 
	CrystalBallFitter(TH1F * h) 
	{
    opt_integrate = false;
		debug =false;
		his = h;
		xmin = his->GetXaxis()->GetXmin();
		xmax = his->GetXaxis()->GetXmax();
		Nbins = his->GetNbinsX();
		N0 = his->GetEntries();
		fun =0;
    inipar.Add("Nsig",  N0,   -N0*0.05); 
    inipar.Add("mean",  -0.368,   1); 
    inipar.Add("sigma",  1.258,   1); 
    inipar.Add("bl",     1.75,   1); 
    inipar.Add("al-bl",  0,     0.1);
    inipar.Add("nl",     3.1,   0.5);
    inipar.Add("br",     1.2,   0.1);
    inipar.Add("ar-br",  0,   0.1);
    inipar.Add("nr",     2.4,   0.1);
    //inipar.Add("kbg",    -0.003746,    1.0/(xmax-xmin));
    inipar.Add("kbg",    0,    1.0/(xmax-xmin));

		inipar.SetLimits("Nsig",  0, N0);
		inipar.SetLimits("mean", xmin, xmax);
		inipar.SetLowerLimit("nl",    1);
		inipar.SetLowerLimit("nr",    1);
		inipar.SetLowerLimit("bl",    0);
		inipar.SetLowerLimit("br",    0);
		inipar.SetLowerLimit("al-bl", 0);
		inipar.SetLowerLimit("ar-br", 0);
		inipar.SetLowerLimit("sigma", 0.1);
		inipar.SetLimits("kbg", -2.0/(xmax-xmin), 2.0/(xmax-xmin));
		inipar.Fix("kbg");
		inipar.Fix("al-bl");
		inipar.Fix("ar-br");
		//inipar.Fix("nr");
		//inipar.Fix("nl");
		//inipar.SetLimits("nl", 1, 10);
		//inipar.SetLimits("nr", 1, 10);
		//double sigma_min = 1;
		//inipar.SetLimits("sigma", sigma_min,  xmax - xmin);
		//inipar.SetLimits("bl", 0.01, 3);
		//inipar.SetLimits("br", 0.01, 3);
		//inipar.SetLimits("al-bl", 0, 3);
		//inipar.SetLimits("ar-br", 0, 3);
	}


	~CrystalBallFitter(void)
	{
		//delete min;
		//delete fun;
		//for(auto g : graph_list) delete g;
	}
  
  void SetIntegrate(bool b) { opt_integrate = b; }

	void Minos(FunctionMinimum & minimum)
	{
		minpar = minimum.UserParameters();
		MnMinos minos(*this, minimum);
		par_error.resize(inipar.Params().size());
		fit_result.resize(inipar.Params().size());
		for (unsigned i=0;i<par_error.size();i++)
		{
			MinuitParameter p  = minpar.Parameter(i);
			fit_result[i] = minpar.Parameter(i).Value();
			if(!p.IsFixed())
			{
				par_error[i] = minos(i, (unsigned)1e9);
				if(p.HasLowerLimit()  &&  p.Value() + par_error[i].first < p.LowerLimit()) par_error[i].first = p.LowerLimit() -p.Value();  
				if(p.HasUpperLimit()  &&  p.Value() + par_error[i].second > p.UpperLimit()) par_error[i].second = p.UpperLimit() - p.Value();  
			}
		}
	}

	void Print(FunctionMinimum & minimum)
	{
		cout << "fcn = " << minimum.Fval() << endl;
		for(unsigned i=0;i<par_error.size();i++)
		{
			cout << setw(5) << i << setw(10) << inipar.Name(i) << setw(15) << fit_result[i] << setw(15) << par_error[i].first << setw(15) << par_error[i].second << endl;
		}
		cout << "Total number of events: " <<  N0 << endl;
		double psig = fit_result[0]/N0;
		cout << setw(25) << "Number of signal events: " <<  setw(15) << fit_result[0] << setw(15) << -sqrt(sq(par_error[0].first) + psig*psig*N0)   << setw(15) << sqrt(sq(par_error[0].second) + psig*psig*N0) << endl;
		cout << setw(25) << "Number of signal events(2): " <<  setw(15) << fit_result[0] << setw(15) << -sqrt(N0*psig)  << setw(15) << sqrt(N0*psig) << endl;
		cout << setw(25) << "Number of background events: " << setw(15) << N0 - fit_result[0] << setw(5) << "+-"   << setw(15) << sqrt(N0-fit_result[0]) << endl;
	}

	void Draw(void)
	{
		fun = new TF1("cbf",*this,  xmin, xmax, inipar.Params().size());
		fun->SetLineColor(kRed);
		fun->SetParameters(&minpar.Params()[0]);
    fun->SetLineWidth(1);
		fun->Draw("same");
	}


	TGraph*  DrawGraph(std::vector<std::pair<double, double>> & v, const char * name = "")
	{
		auto canvas = new TCanvas(name, name);
    canvas->cd();
		TGraph * g = new TGraph;
		sort(begin(v), end(v), [](const pair<double, double> & p1,  const pair<double, double> & p2) { return p1.first < p2.first;});
		for(unsigned i=0;i<v.size();i++) g->SetPoint(i, v[i].first,  v[i].second);
		g->Draw("al");
		g->SetName(name);
		g->SetTitle(name);
		graph_list.push_back(g);
		return g;
	}


	void Scan(const ROOT::Minuit2::MnUserParameters &  upar)
	{
		MnScan scan(*this, upar, 2);
		//vector<double> e(10, 1);
		//vector<double> p = {2663.75,  -0.365978,    1.27853,    1.57638 ,   1.87204,    1.00001,     1.3674,    4.56482,    3.14814,  0.0396767};
		//MnScan scan(*this, p, e);
		for(unsigned i=0;i<upar.Params().size();i++)
		{
			auto v = scan.Scan(i, 1000);
			DrawGraph(v, upar.Name(i));
		}
	}

	void Fit(void)
	{
		//MnStrategy strategy(2);
		//strategy.SetGradientNCycles(1000000);
		//strategy.SetHessianGradientNCycles(1000000);
		//Scan(inipar);
		//new TCanvas;
		his->Draw();
    MnMigrad migrad(*this, inipar, 2);
		auto minimum = migrad(unsigned(1e9));
    //FunctionMinimum minimum = migrad(unsigned(1e9));
		cout << "NumCall: " << migrad.NumOfCalls() << endl;
		cout << "Nfcn: " << minimum.NFcn() << endl;
		if(!minimum.IsValid()) 
		{
			cout << "Minimum invalid " << endl;
			cout << minimum << endl;
			minpar = minimum.UserParameters();
			Scan(minpar);
			Draw();
			return;
		}
		Minos(minimum);
		//MnMinos minos(*this, minimum);
		//par_error.resize(inipar.Params().size());
		//fit_result.resize(inipar.Params().size());
		//for (unsigned i=0;i<par_error.size();i++)
		//{
		//	MinuitParameter p  = minpar.Parameter(i);
		//	fit_result[i] = minpar.Parameter(i).Value();
		//	if(!p.IsFixed())
		//	{
		//		par_error[i] = minos(i, (unsigned)1e9);
		//		if(p.HasLowerLimit()  &&  p.Value() + par_error[i].first < p.LowerLimit()) par_error[i].first = p.LowerLimit() -p.Value();  
		//		if(p.HasUpperLimit()  &&  p.Value() + par_error[i].second > p.UpperLimit()) par_error[i].second = p.UpperLimit() - p.Value();  
		//	}
		//}
		Print(minimum);
		Draw();
		//Scan(minpar);
	}

	double operator()( const std::vector<double> & par) const 
	{
		double chi2=0;
		for(int i=0; i<his->GetNbinsX(); i++)
		{
			double n = his->GetBinContent(i);
      double mu;

      if(opt_integrate)
      {
        double xbin_low=his->GetBinLowEdge(i);
        double xbin_up = xbin_low + his->GetBinWidth(i);
        auto f = [&par,this](double z) { return ModifiedDoubleCrystalBall(&z, &par[0]); };
        mu = ibn::dgaus(f, xbin_low, xbin_up,1e-10);
      }
      else
      {
        double x = his->GetBinCenter(i);
        mu = ModifiedDoubleCrystalBall(&x, &par[0]);
      }

			double dchi2;
			if(n!=0)
			{
				//if(mu<=0) return 1e100;
				dchi2 = 2*(mu - n + n*log(n/mu));
			}
			else
			{
				dchi2 = 2*(mu-n);
			}
			chi2+=dchi2;
		}
		//add limits for beta and alfa
		//double mean =  par[1];
		//double sigma =  par[2];
		//double beta[2] = {par[3], par[6]};
		//double alfa[2] = {par[4] + beta[0], par[7] + beta[1]};
		//beta[0] = mean - beta[0]*sigma; 
		//beta[1] = mean + beta[1]*sigma; 
		//alfa[0] = mean - alfa[0]*sigma; 
		//alfa[1] = mean + alfa[1]*sigma; 

		///double n[2] = {par[5],  par[8]};
		//if(alfa[0] <= xmin) 
		//{
		//	chi2 *= exp(pow((alfa[0]-xmin)/sigma, 2.0)); 
		//	//chi2 += pow(n[0]-1, 2);
		//}
		//if(alfa[1] >= xmax) 
		//{
		//	chi2 *= exp(pow((alfa[1]-xmax)/sigma, 2.0)); 
		//	//chi2 += pow(n[1]-1, 2);
		//}
		//if(alfa[0] <= xmin) return 1e100;
		//if(alfa[1] >= xmax) return 1e100;


		//if(isnan(chi2) || isinf(chi2)) return 0.5*std::numeric_limits<double>::max();
		//if(std::isnan(chi2) || std::isinf(chi2)) chi2= 1e100;
		//cout << chi2 << endl;
		//for(int i=0;i<10;i++) cout << setw(15) << par[i];
		//cout << setw(15) << result << endl;
		return reg(chi2, par);
		return chi2;
	}

	double Up(void) const { return 1.0;}

	//this is parameterized function
	double operator()(const double * x ,  double * p) const 
	{ 
		return ModifiedDoubleCrystalBall(x, p);
	} 



	double ModifiedDoubleCrystalBall(const double* X, const double* P) const
	{
		double x    =  X[0];
		double Nsig  = P[0];
		double mean    = P[1];//mean of the gaus
		double sigma = P[2]; //sigma of the gaus



		//tails 0 - left tail,  1 - right tail
		double beta[2] = {P[3],           P[6]};
		double alfa[2] = {P[4] + beta[0], P[7] + beta[1]};
		double n[2] =    {P[5],P[8]};

		double pbg = P[9]; //background slope in sigma units

		//background
		double NBG = N0 - Nsig;


		//some constants
		double C[2]; //correspond exp tail
		double A[2]; //correspond stepennoi tail
		double B[2]; //correspon stepennoi tail handy coef
		for(int i=0;i<2;i++)
		{
			C[i] = TMath::Exp(0.5*beta[i]*beta[i]);
			A[i] = C[i]*pow(n[i]/beta[i], n[i]) * TMath::Exp(- beta[i]*alfa[i]);
			B[i] = n[i]/beta[i] - alfa[i];
		}

		//scale and shift the x
		x = (x-mean)/sigma;

		double xrange[2] = {-(xmin-mean)/sigma, (xmax-mean)/sigma };

		int t = x < 0 ? 0 : 1; //tail index
		double y = TMath::Abs(x); //abs of x
		double cb=0; //crystal ball function result gaus is 1

		if(y >= 0  && y <= beta[t]) cb = TMath::Exp( - 0.5*y*y);
		if(y > beta[t] && y <= alfa[t]) cb = C[t]*TMath::Exp(- beta[t] * y);
		if(y > alfa[t]) cb = A[t]*pow(B[t] + y,  - n[t]);


		double Igaus[2]; //gauss part of integral
		double Iexp[2]; //exp part of integral
		double Ipow[2]; //power part of integral
		double I=0; //sum of previouse one

		for(int i=0;i<2;i++)
		{
			//calculate partial integral
			if( xrange[i] <= beta[i] ) //only gause in range
			{
				Igaus[i] = sqrt(TMath::Pi()*0.5)*TMath::Erf(xrange[i]/sqrt(2.0));
				Iexp[i] = 0;
				Ipow[i] = 0;
			}
			if( beta[i] <  xrange[i] && xrange[i] <= alfa[i]) //gaus and part of exp tail in range
			{
				Igaus[i] = sqrt(TMath::Pi()*0.5)*TMath::Erf(beta[i]/sqrt(2.0));
				Iexp[i] =  C[i]/beta[i]*( TMath::Exp(-beta[i]*beta[i]) -  TMath::Exp(-beta[i]*xrange[i]) );
				Ipow[i] =0;
			}
			if(alfa[i] < xrange[i]) //all part of function inside the range
			{
				Igaus[i] = sqrt(TMath::Pi()*0.5)*TMath::Erf(beta[i]/sqrt(2.0));
				Iexp[i] =  C[i]/beta[i]*( TMath::Exp(-beta[i]*beta[i]) -  TMath::Exp(-beta[i]*alfa[i]));
				if(n[i]==1.0) Ipow[i] = A[i]*log((B[i] + xrange[i])/(B[i] + alfa[i]));
				else Ipow[i] = A[i]/(n[i]-1.0)*(pow( B[i] + alfa[i], - n[i] + 1.)  - pow( B[i] + xrange[i], - n[i] + 1.));
			}
			//calculate total integral
			I+=Igaus[i] + Iexp[i] +Ipow[i];
		}
		//calculate the scale
		double dx = (xrange[0] + xrange[1])/Nbins;


		//the amplitude of signal
		double Ns = Nsig/I*dx;

		//the amplitude of background
		//double Nbg =  NBG*dx/(xrange[0]+xrange[1])/(1.0 + 0.5*pbg*(xrange[1]-xrange[0]));
		double Nbg =  NBG/Nbins;
		//cout << Nbg << " " << xrange[1] - xrange[0] << endl;

		//result of my double crystal ball function
		double result = Ns*cb + Nbg*(1.0 + pbg* (X[0] - 0.5*(xmax+xmin)));

		//test the result to suppress bad calculation
		//if(TMath::IsNaN(result) || result <=0) return 0;
		//if(!TMath::Finite(result)) return 1e100;
		//for(int i=0;i<10;i++) cout << setw(15) << P[i];
		//cout << setw(15) << result << endl;
		return result;
	}


  std::vector<double> get_result(void) const 
  {
    std::vector<double> buf(3);
    buf[0] = fit_result[0];
    double psig =  fit_result[0]/N0;
    buf[1] = sqrt(sq(par_error[0].first) + psig*psig*N0);
    buf[2] = sqrt(sq(par_error[0].second) + psig*psig*N0);
    return buf;
  }
};



std::vector<double> Fit(TH1F * his)
{
	CrystalBallFitter cb(his);
	cb.Fit();
	return cb.get_result();
}
