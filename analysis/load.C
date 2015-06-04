#include <regex>
#include <math.h>

#include <list>
#include <TRandom.h>

#include "CrystalBall.h"
//TRandom RANDOM;

list<string> list_files(const char *dirname=".", const char *ext=".root")
{

	list<string> flist;
	cout << dirname << endl;
	TSystemDirectory dir(dirname, dirname);
	TList *files = dir.GetListOfFiles();
	if (files) 
	{
		TSystemFile *file;
		TString fname;
		TIter next(files);
		while ((file=(TSystemFile*)next()))
		{
			fname = file->GetName();
			if (!file->IsDirectory() && fname.EndsWith(ext))
			{
				flist.push_back(fname.Data());
			}
		}
	}
	flist.sort();
	return flist;
}
//TChain * load_tree(const char *dirname=".", const char *ext=".root")
TChain * load_tree(string dirname=".", string ext=".root")
{
	TChain * event = new TChain("event", "event");
	TString pwd(gSystem->pwd());
	TSystemDirectory dir(dirname.c_str(), dirname.c_str());
	TList *files = dir.GetListOfFiles();
	gSystem->cd(pwd.Data());
	if (files) 
	{
		TSystemFile *file;
		TString fname;
		TIter next(files);
		while ((file=(TSystemFile*)next()))
		{
			fname = file->GetName();
			if (!file->IsDirectory() && fname.EndsWith(ext.c_str()))
			{
				string f = dirname + "/" + fname.Data();
				event->AddFile(f.c_str());
				//flist.push_back(fname.Data());
			}
		}
	}
	return event;
}

double CrystalBall(double* x, double* par)
{
	//http://en.wikipedia.org/wiki/Crystal_Ball_function
	double xcur = x[0];
	double alpha = par[0];
	double n = par[1];
	double mu = par[2];
	double sigma = par[3];
	double N = par[4];
	TF1* exp = new TF1("exp","exp(x)",1e-20,1e20);
	double A; double B;
	if (alpha < 0)
	{
		A = pow((n/(-1*alpha)),n)*exp->Eval((-1)*alpha*alpha/2);
		B = n/(-1*alpha) + alpha;
	}
	else
	{
		A = pow((n/alpha),n)*exp->Eval((-1)*alpha*alpha/2);
		B = n/alpha - alpha;
	}
	double f;
	if ((xcur-mu)/sigma > (-1)*alpha)
		f = N*exp->Eval((-1)*(xcur-mu)*(xcur-mu)/ (2*sigma*sigma));
	else
		f = N*A*pow((B- (xcur-mu)/sigma),(-1*n));
	delete exp;
	return f;
}

double DoulbeCrystalBall(double* X, double* P)//left and right
{
	double x = X[0];
	double N = P[0];
	double mu = P[1];//mean of the gaus
	double sigma = P[2]; //sigma of the gaus

	double alfa_l = P[3];
	double n_l = P[4];

	double alfa_r = P[5];
	double n_r = P[6];

	double bg = P[7];
	//scale and shift the x
	x = (x-mu)/sigma;


	double A_l = pow(n_l/TMath::Abs(alfa_l), n_l) * TMath::Exp(-alfa_l*alfa_l/2.0);
	double A_r = pow(n_r/TMath::Abs(alfa_r), n_r) * TMath::Exp(-alfa_r*alfa_r/2.0);
	double B_l = n_l / TMath::Abs(alfa_l)  - TMath::Abs(alfa_l);
	double B_r = n_r / TMath::Abs(alfa_r)  - TMath::Abs(alfa_r);

	double result;
	if(x < - alfa_l) 
	{
		result =  N*A_l*pow(B_l - x,  - n_l);
	}
	else 
	{
		if(x > alfa_r) 
		{
			result =  N*A_r*pow(x + B_r,  - n_r);
		}
		else 
		{
			result = N*TMath::Exp( - x*x/2.0);
		}
	}
	if(TMath::IsNaN(result) || !TMath::Finite(result)) return 0;
	return result+N*bg;
}

double TripleCrystalBall(double* X, double* P)//left and right
{
	double x = X[0];
	double N = P[0];
	double mu = P[1];//mean of the gaus
	double sigma = P[2]; //sigma of the gaus

	double alfa_l = P[3];
	double n_l = P[4];

	double alfa_r = P[5];
	double n_r = P[6];

	double alfa_r2 = P[7]+alfa_r;
	double n_r2 = P[8];

	double bg = P[9];
	//scale and shift the x
	x = (x-mu)/sigma;


	double A_l = pow(n_l/TMath::Abs(alfa_l), n_l) * TMath::Exp(-alfa_l*alfa_l/2.0);
	double A_r = pow(n_r/TMath::Abs(alfa_r), n_r) * TMath::Exp(-alfa_r*alfa_r/2.0);
	double B_l = n_l / TMath::Abs(alfa_l)  - TMath::Abs(alfa_l);
	double B_r = n_r / TMath::Abs(alfa_r)  - TMath::Abs(alfa_r);

	double B_r2 = n_r2/n_r*(alfa_r2 + B_r) - alfa_r2;
	double A_r2 = A_r*pow(alfa_r2 + B_r2,  n_r2)*pow(alfa_r2 + B_r,  - n_r);

	double result;
	if(x < - alfa_l) 
	{
		result =  N*A_l*pow(B_l - x,  - n_l);
	}
	if( -alfa_l <= x && x <= alfa_r)
	{
		result = N*TMath::Exp( - x*x/2.0);
	}
	if( alfa_r < x && x <= alfa_r2)
	{
		result =  N*A_r*pow(x + B_r,  - n_r);
	}
	if(  x > alfa_r2)
	{
		result =  N*A_r2*pow(x + B_r2,  - n_r2);
	}
	if(TMath::IsNaN(result) || !TMath::Finite(result)) return 0;
	return result+N*bg;
}

double ExpCrystalBall(double* X, double* P)//left and right
{
	double x = X[0];
	double N = P[0];
	double mu = P[1];//mean of the gaus
	double sigma = P[2]; //sigma of the gaus

	//left tail
	double alfa_l = TMath::Abs(P[3]);
	double n_l = TMath::Abs(P[4]);

	//exp right tail
	double alfa_e = TMath::Abs(P[5]);

	//right tail
	double alfa_r = TMath::Abs(P[6])+alfa_e;
	double n_r = TMath::Abs(P[7]);

	//background
	double bg = TMath::Abs(P[8]);

	//scale and shift the x
	x = (x-mu)/sigma;


	//calculate left tail params
	double A_l = pow(n_l/alfa_l, n_l) * TMath::Exp(-0.5*alfa_l*alfa_l);
	double B_l = n_l /alfa_l   - alfa_l;


	//calculate right exp tail params
	double A_e = TMath::Exp(0.5*alfa_e*alfa_e);

	//calculate right tail params
	double B_r = n_r / alfa_e  - alfa_r;
	double A_r = A_e*pow(n_r/alfa_e, n_r) * TMath::Exp(- alfa_r*alfa_e);

	//calculate result
	double result;
	//if on the left tail
	if(x < - alfa_l) 
	{
		result =  N*A_l*pow(B_l - x,  - n_l);
	}
	//if on the gaus
	if( -alfa_l <= x && x <= alfa_e)
	{
		result = N*TMath::Exp( - x*x/2.0);
	}
	//on the right exp tail
	if( alfa_e < x && x <= alfa_r)
	{
		result =  N*A_e*TMath::Exp(-alfa_e*x);
	}
	//on the right tail
	if(  x > alfa_r)
	{
		result =  N*A_r*pow(x + B_r,  - n_r);
	}
	if(TMath::IsNaN(result) || !TMath::Finite(result)) return 0;
	return result+N*bg;
}

double DoubleExpCrystalBall(double* X, double* P)//left and right
{
	double x = X[0];
	double N = P[0];
	double mu = P[1];//mean of the gaus
	double sigma = P[2]; //sigma of the gaus


	//exp left tail
	double alfa_el = TMath::Abs(P[3]);

	//left tail
	double alfa_l = alfa_el + TMath::Abs(P[4]);
	double n_l = TMath::Abs(P[5]);

	//exp right tail
	double alfa_er = TMath::Abs(P[6]);

	//right tail
	double alfa_r = TMath::Abs(P[7])+alfa_er;
	double n_r = TMath::Abs(P[8]);

	//background
	double bg0 = TMath::Abs(P[9]);
	double bg1 = P[10];

	//scale and shift the x
	x = (x-mu)/sigma;


	//calculate right exp tail params
	double A_el = TMath::Exp(0.5*alfa_el*alfa_el);

	//calculate left tail params
	double B_l;
	if(alfa_el == 0) B_l = n_l /alfa_l   - alfa_l;
	else B_l = n_l /alfa_el   - alfa_l;

	double A_l;
	if(alfa_el == 0) A_l = pow(n_l/alfa_l, n_l) * TMath::Exp(-0.5*alfa_l*alfa_l);
	else A_l = A_el*pow(n_l/alfa_el, n_l) * TMath::Exp(- alfa_l*alfa_el);


	//calculate right exp tail params
	double A_er = TMath::Exp(0.5*alfa_er*alfa_er);

	//calculate right tail params
	double B_r = n_r / alfa_er  - alfa_r;
	double A_r = A_er*pow(n_r/alfa_er, n_r) * TMath::Exp(- alfa_r*alfa_er);

	//calculate result
	double result;
	//if on the left tail
	if(x < - alfa_l) 
	{
		result =  N*A_l*pow(B_l - x,  - n_l);
	}
	if(alfa_el  == 0)
	{
		if(-alfa_l < x && x <=0)
		{
			result = N*TMath::Exp( - x*x/2.0);
		}
	}
	else
	{
		if( -alfa_l <= x && x< -alfa_el )
		{
			result =  N*A_el*TMath::Exp(alfa_el*x);
		}
		if( -alfa_el <= x && x <= 0)
		{
			result = N*TMath::Exp( - x*x/2.0);
		}
	}
	//if on the gaus
	if( 0 < x && x <= alfa_er)
	{
		result = N*TMath::Exp( - x*x/2.0);
	}
	//on the right exp tail
	if( alfa_er < x && x <= alfa_r)
	{
		result =  N*A_er*TMath::Exp(-alfa_er*x);
	}
	//on the right tail
	if(  x > alfa_r)
	{
		result =  N*A_r*pow(B_r + x,  - n_r);
	}
	if(TMath::IsNaN(result) || !TMath::Finite(result)) return 1e500;
	return result+bg0*(1.0 + bg1*(X[0]-mu));
}



double DoubleExpCrystalBallNorm(double* X, double* P)//left and right
{
	double x = X[0];
	double Nsig = P[0];
	double mu = P[1];//mean of the gaus
	double sigma = P[2]; //sigma of the gaus


	//exp left tail
	double alfa_el = TMath::Abs(P[3]);

	//left tail
	double alfa_l = alfa_el + TMath::Abs(P[4]);
	double n_l = TMath::Abs(P[5]);

	//exp right tail
	double alfa_er = TMath::Abs(P[6]);

	//right tail
	double alfa_r = TMath::Abs(P[7])+alfa_er;
	double n_r = TMath::Abs(P[8]);

	//background
	double Nbg = TMath::Abs(P[9]);
	//double bg0 = TMath::Abs(P[9]);
	double bg1 = P[10];

	double xmin = P[11];
	double xmax = P[12];
	double Nbins = P[13];

	double bg0=Nbg/Nbins/(1.0 + 0.5*bg1*(xmax+xmin-2*mu));

	//scale and shift the x
	x = (x-mu)/sigma;


	//calculate right exp tail params
	double A_el = TMath::Exp(0.5*alfa_el*alfa_el);

	//calculate left tail params
	double B_l;
	if(alfa_el == 0) B_l = n_l /alfa_l   - alfa_l;
	else B_l = n_l /alfa_el   - alfa_l;

	double A_l;
	if(alfa_el == 0) A_l = pow(n_l/alfa_l, n_l) * TMath::Exp(-0.5*alfa_l*alfa_l);
	else A_l = A_el*pow(n_l/alfa_el, n_l) * TMath::Exp(- alfa_l*alfa_el);


	//calculate right exp tail params
	double A_er = TMath::Exp(0.5*alfa_er*alfa_er);

	//calculate right tail params
	double B_r = n_r / alfa_er  - alfa_r;
	double A_r = A_er*pow(n_r/alfa_er, n_r) * TMath::Exp(- alfa_r*alfa_er);


	double Igaus = sqrt(TMath::Pi()/2.0)*( TMath::Erf(alfa_er/sqrt(2.0)) -  TMath::Erf(-alfa_l/sqrt(2.0)));
	double Ier = TMath::Exp(-0.5*alfa_er*alfa_er)/alfa_er*(1.0 - TMath::Exp(-alfa_er*(alfa_r-alfa_er)));
	double Ir = A_r/(n_r - 1.0)*(pow(B_r+alfa_r, -n_r+1.0) - pow(B_r+(xmax-mu)/sigma, -n_r+1.0));
	double Il = A_l/(n_l - 1.0)*(pow(B_l+alfa_l, -n_l+1.0) - pow(B_l-(xmin-mu)/sigma, -n_l+1.0));


	double N = Nsig/(Igaus+Ier+Ir+Il)/Nbins*(xmax-xmin);
	//cout << Igaus << " " << Ier << " " << Ir << " " << Il << " N = " << N << endl;

	//calculate result
	double result;
	//if on the left tail
	if(x < - alfa_l) 
	{
		result =  N*A_l*pow(B_l - x,  - n_l);
	}
	if(alfa_el  == 0)
	{
		if(-alfa_l < x && x <=0)
		{
			result = N*TMath::Exp( - x*x/2.0);
		}
	}
	else
	{
		if( -alfa_l <= x && x< -alfa_el )
		{
			result =  N*A_el*TMath::Exp(alfa_el*x);
		}
		if( -alfa_el <= x && x <= 0)
		{
			result = N*TMath::Exp( - x*x/2.0);
		}
	}
	//if on the gaus
	if( 0 < x && x <= alfa_er)
	{
		result = N*TMath::Exp( - x*x/2.0);
	}
	//on the right exp tail
	if( alfa_er < x && x <= alfa_r)
	{
		result =  N*A_er*TMath::Exp(-alfa_er*x);
	}
	//on the right tail
	if(  x > alfa_r)
	{
		result =  N*A_r*pow(B_r + x,  - n_r);
	}
	result+=bg0*(1.0 + bg1*(X[0]-mu));
	if(TMath::IsNaN(result) || !TMath::Finite(result)) return 1e500;
	return result;
}





//TChain * load_tree(const char * dirname)
//{
//	TChain * event = new TChain("event", "event");
//	list<string> lst = list_files(dirname, ".root");
//	int filenumber=0;
//	for(list<string>::iterator l = lst.begin();l!=lst.end(); l++)
//	{
//		cout << "Loading " << *l << endl;
//		event->AddFile(l->c_str());
//		filenumber++;
//	}
//	cout << filenumber << " files loaded" << endl;
//	return event;
//}


void drawMinv(void)
{
	TCut cut;
	//TCut cut="Mrec>3.09 && Mrec<3.102";
	TCanvas * c = new TCanvas();
	TChain * mcKK = load_tree("mcjpkk2009");
	TChain * mcK1K = load_tree("mcK1K2009");
	mcKK->SetMarkerColor(kBlue);
	mcK1K->SetMarkerColor(kRed);
	mcKK->SetLineColor(kBlue);
	mcK1K->SetLineColor(kRed);
	mcKK->Draw("M012:Mrec","KK" && cut, "", 1e5);
	mcK1K->Draw("M012:Mrec","KK" && cut,  "same", 2e3);
	new TCanvas;
	//mcKK->Draw("M012", "KK" && cut, "", 1e5);
	//cout << "K1K events:" << mcK1K->GetEntries("KK" && cut) << endl;
	mcKK->Draw("M012", "KK" && cut, "", 100000);
	mcK1K->Draw("M012", "KK" && cut, "same", 1000);

}


TRandom * RANDOM;



//chebyshev function
inline double T0(double x)
{
	return 1;
}

inline double T1(double x)
{
	return x;
}

inline double T2(double x)
{
	return 2.0*x*x -1.0;
}

inline double T3(double x)
{
	return 4.0*x*x*x -3*x;
}

inline double T4(double x)
{
	return 8.0*pow(x, 4)-8.0*x*x+1;
}

double fitfun_fitMrecKK(double *X,  double *p)
{
	//gaus + koshi
	double AG = p[0];
	double MeanG = p[1];
	double SigmaG = p[2];
	double AK = p[3];
	double MeanK = p[4];
	double SigmaK = p[5];
	double AG2 = p[6];
	double MeanG2 = p[7];
	double SigmaG2 = p[8];
	MeanK=MeanG+p[4];
	MeanG2=MeanG+p[7];
	double x = *X;
	double gaus = AG*TMath::Gaus(x, MeanG, SigmaG);
	double koshi = AG*AK/(1.0 + ((x-MeanK)/SigmaK)**2.0);
	double gaus2 = AG*AG2*TMath::Gaus(x, MeanG2, SigmaG2);
	//return gaus + gaus2 + koshi;
	//return gaus*(1.0 + gaus2);
	//return gaus*(1.0 + koshi);
}

double tailgaus_fitMrecKK(double x,  double A,  double Mean,  double Sigma,  double C)
{
	if( TMath::Abs(x-Mean)/Sigma < C || C<=0) return A*TMath::Gaus(x,  Mean,  Sigma);
	double B = A*TMath::Gaus(C, 0, 1);
	double L = Sigma/C;
	//cout << A << " " <<B << " " << L << endl;
	return B*TMath::Exp(-(TMath::Abs(x-Mean)-C*Sigma)/L);
}

double tailgaus_fitMrecKK(double *X,  double*P)
{
	return tailgaus_fitMrecKK(X[0],  P[0], P[1], P[2], P[3]);
}

void fitMrecKK(void)
{
	TChain * c = load_tree("mcjpkk2009");
	TCut cut = "KK && Mrec > 3.09 && Mrec < 3.104 && kin_chi2<40 && pid_chi2<20";
	c->Draw("Mrec>>hfitMrecKK(2000)",cut,"");
	TF1 *fitfun_fitMrecKK = new TF1("fitfun_fitMrecKK", &fitfun_fitMrecKK, 3.09,  3.104, 9);
	fitfun_fitMrecKK->SetParameter(0, 1569);
	fitfun_fitMrecKK->SetParameter(1, 3.09661);
	fitfun_fitMrecKK->SetParameter(2, 0.001);
	fitfun_fitMrecKK->SetParameter(3, 0.637);
	fitfun_fitMrecKK->SetParameter(4, 0);
	fitfun_fitMrecKK->SetParameter(5, 0.005);
	fitfun_fitMrecKK->SetParameter(6, 0);
	fitfun_fitMrecKK->SetParameter(7, 0);
	fitfun_fitMrecKK->SetParameter(8, 0.005);
	//fitfun_fitMrecKK->FixParameter(4, 1);
	//fitfun_fitMrecKK->FixParameter(7, 1);
	//fixing koshi
	fitfun_fitMrecKK->FixParameter(3, 1);
	fitfun_fitMrecKK->FixParameter(4, 1);
	fitfun_fitMrecKK->FixParameter(5, 1);
	//fixing gaus2
	//fitfun_fitMrecKK->FixParameter(6, 0);
	//fitfun_fitMrecKK->FixParameter(7, 0);
	//fitfun_fitMrecKK->FixParameter(8, 1);
	//fixing main gaus
	fitfun_fitMrecKK->FixParameter(6, 0);
	fitfun_fitMrecKK->FixParameter(7, 0);
	fitfun_fitMrecKK->FixParameter(8, 1);
	hfitMrecKK->Fit("fitfun_fitMrecKK", "R");
}

void fitMrecKK2(void)
{
	TChain * c = load_tree("mcjpkk2009");
	TCut cut = "KK && Mrec > 3.09 && Mrec < 3.104 && kin_chi2<40 && pid_chi2<20";
	c->Draw("Mrec>>hfitMrecKK(2000)",cut,"");
	TF1 *fun = new TF1("fitfun_fitMrecKK2", &tailgaus_fitMrecKK, 3.09,  3.104, 4);
	fun->SetParameter(0, 1569);
	fun->SetParameter(1, 3.09661);
	fun->SetParameter(2, 0.001);
	fun->SetParameter(3, 10);
	fun->SetParLimits(3, 0, 100);
	hfitMrecKK->Fit("fitfun_fitMrecKK2", "R");


	TH1F * his = new TH1F("his", "his", 2000, 3.09, 3.104);
	his->Eval(fun);
	//TH1F * his = fun->GetHistogram();
	hfitMrecKK->Add(his, -1);
}


void fitCB(void)
{
	TChain * c = load_tree("mcjpkk2009");
	TCut cut = "KK && Mrec > 3.09 && Mrec < 3.104 && kin_chi2<40 && pid_chi2<20";
	c->Draw("Mrec>>hfitMrecKK(2000)",cut,"");
	TF1 *fun = new TF1("crystal_ball_fun", &CrystalBall, 3.09,  3.104, 5);
	fun->SetParameter(0, 1);
	fun->SetParameter(1, 1);
	fun->SetParameter(2, 3.09661);
	fun->SetParameter(3, 0.001);
	fun->SetParameter(4, 1500);
	hfitMrecKK->Fit("crystal_ball_fun", "R");
}

void fitCB2(void)
{
	TCanvas *canvas = new TCanvas("Double Crystabl ball fit to the KK Monte Carlo", "Mrec_KKfit_CB2");
	canvas->SetLogy();
	TChain * c = load_tree("mcjpkk2009");
	//TCut cut = "KK && Mrec > 3.09 && Mrec < 3.104 && kin_chi2<40 && pid_chi2<20";
	TCut cut = "KK && Mrec > 3.08 && Mrec < 3.114 && kin_chi2<40 && pid_chi2<20";
	c->Draw("Mrec>>hfitMrecKK(2000)",cut,"");
	//TF1 *fun = new TF1("double_crystal_ball_fun", &DoulbeCrystalBall, 3.09,  3.104, 8);
	TF1 *fun = new TF1("double_crystal_ball_fun", &DoulbeCrystalBall, 3.08,  3.114, 8);
	//gaus
	fun->SetParameter(0, 2450);
	fun->SetParameter(1, 3.09661);
	fun->SetParameter(2, 0.001256);
	//left CB
	fun->SetParameter(3, 1.623);
	fun->SetParameter(4, 5.079);
	//right CB
	fun->SetParameter(5, 1.295);
	fun->SetParameter(6, 11.38);
	//background
	fun->SetParameter(7, 0);
	//fun->FixParameter(7, 0);
	fun->SetLineColor(kRed);
	gStyle->SetOptFit();
	hfitMrecKK->Fit("double_crystal_ball_fun", "R");
	hfitMrecKK->GetXaxis()->SetTitle("M_{rec}(#pi^{+}#pi^{-}),  GeV");
}
void fitCB3(const char * dir="mcjpkk2009", int Nbins=1000)
{
	TCanvas *canvas = new TCanvas("Triple Crystabl ball fit to the KK Monte Carlo", "Mrec_KKfit_CB3");
	canvas->SetLogy();
	TChain * c = load_tree(dir);
	//TCut cut = "KK && Mrec > 3.09 && Mrec < 3.104 && kin_chi2<40 && pid_chi2<20";
	TCut cut = "KK && Mrec > 3.08 && Mrec < 3.114 && kin_chi2<40 && pid_chi2<20";
	char  selection_buf[1024];
	sprintf(selection_buf, "Mrec>>hfitMrecKK3(%d)", Nbins);
	c->Draw(selection_buf,cut,"");
	//TF1 *fun = new TF1("triple_crystal_ball_fun", &DoulbeCrystalBall, 3.09,  3.104, 8);
	TF1 *fun = new TF1("triple_crystal_ball_fun", &TripleCrystalBall, 3.08,  3.114, 10);
	//gaus
	fun->SetParName(0, "N");
	//fun->SetParameter(0, 2450);
	fun->SetParameter(0, hfitMrecKK3->GetMaximum());
	fun->SetParName(1, "mean");
	fun->SetParameter(1, 3.09661);
	fun->SetParName(2, "sigma");
	fun->SetParameter(2, 0.001256);
	//left CB
	fun->SetParName(3, "Al");
	fun->SetParameter(3, 1.623);
	fun->SetParName(4, "Nl");
	fun->SetParameter(4, 5.079);
	//right CB
	fun->SetParName(5, "Ar");
	fun->SetParameter(5, 1.295);
	fun->SetParName(6, "Nr");
	fun->SetParameter(6, 11.38);
	//right CB2
	fun->SetParName(7, "Ar2");
	fun->SetParameter(7, 0.1);
	fun->SetParLimits(7, 0, 100);
	fun->SetParName(8, "Nr2");
	fun->SetParameter(8, 11.38);
	fun->SetParLimits(8, 0.1, 100);
	//background
	fun->SetParName(9, "BG");
	fun->SetParameter(9, 0);
	//fun->FixParameter(9, 0);
	fun->SetLineColor(kRed);
	gStyle->SetOptFit();
	hfitMrecKK3->Fit("triple_crystal_ball_fun", "R");
	hfitMrecKK3->GetXaxis()->SetTitle("M_{rec}(#pi^{+}#pi^{-}),  GeV");
}

void fitECB2(const char * dir="mcjpkk2009", int Nbins=1000,  double xmin=3.09,  double xmax=3.104)
{
	const char * fun_name = "exp_crystal_ball_fun";
	const char * his_name = "hfit_ecb";
	double shift =3.0967;
	TCanvas *canvas = new TCanvas("Doulble Crystabl ball with left exp tail  fit to the Mrec", fun_name);
	canvas->SetLogy();
	TChain * c = load_tree(dir);
	//TCut cut = "KK && Mrec > 3.09 && Mrec < 3.104 && kin_chi2<40 && pid_chi2<20";
	char cut_buf[1024];
	sprintf(cut_buf, "Mrec-%f> %f && Mrec-%f < %f", shift, xmin- shift,  shift,  xmax-shift);
	TCut cut = cut_buf;
	cut = cut && "KK && kin_chi2<40 && pid_chi2<20";
	char  selection_buf[1024];
	sprintf(selection_buf, "Mrec-%f>>%s(%d)", shift, his_name,  Nbins);
	c->Draw(selection_buf,cut,"");
	TH1F * his = c->GetHistogram();
	TF1 *fun = new TF1(fun_name, &ExpCrystalBall, xmin-shift,  xmax-shift, 9);
	//gaus
	fun->SetParName(0, "N");
	//fun->SetParameter(0, 2450);
	fun->SetParameter(0, his->GetMaximum());
	fun->SetParName(1, "mean");
	//fun->SetParameter(1, 3.09661);
	fun->SetParameter(1, 3.0967-shift);
	fun->SetParName(2, "sigma");
	fun->SetParameter(2, 0.001256);
	//left CB
	fun->SetParName(3, "Al");
	fun->SetParameter(3, 1.623);
	fun->SetParName(4, "Nl");
	fun->SetParameter(4, 5.079);
	fun->SetParLimits(4, 0, 1000);
	//exp
	fun->SetParName(5, "Ae");
	fun->SetParameter(5, 1.295);

	//right CB
	fun->SetParName(6, "Ar");
	fun->SetParameter(6, 1.295);
	//fun->FixParameter(6, 100);
	fun->SetParName(7, "Nr");
	fun->SetParameter(7, 11.38);
	//background
	fun->SetParName(8, "BG");
	fun->SetParameter(8, 0);
	//fun->FixParameter(8, 0);
	fun->SetLineColor(kRed);
	gStyle->SetOptFit();
	his->Fit(fun_name, "IR");
	his->GetXaxis()->SetTitle("M_{rec}(#pi^{+}#pi^{-}),  GeV");
}

TF1 *  fitE2CB2(const char * dir="mcjpkk2009", const char * channel="KK",  int Nbins=1000,  double xmin=3.09,  double xmax=3.104,  const char * fitopt="RL",  const char * gopt="SAMES")
{
	//int r  = RANDOM->Integer(10000);
	int r = 0;
	const double SCALE=1e3;
	char fun_name[1024];
	char his_name[1024];
	sprintf(fun_name, "fun_2E2CB_%s_%d",dir,  r);
	sprintf(his_name, "his_2E2CB_%s_%d",dir,  r);
	r++;
	double shift =3.097;
	TChain * c = load_tree(dir);
	char cut_buf[1024];
	double range = xmax -xmin;
	sprintf(cut_buf, "%s && Mrec-%f> %f && Mrec-%f < %f", channel,  shift, xmin- shift,  shift,  xmax-shift);
	TCut cut = cut_buf;
	cut = cut && "kin_chi2<40 && pid_chi2<20";
	char  selection_buf[1024];

	sprintf(selection_buf, "%f*(Mrec-%f)>>%s(%d)", SCALE, shift, his_name,  Nbins);
	//c->Draw(selection_buf,cut,"goff");
	c->Draw(selection_buf,cut, "goff");
	double N0 = c->GetSelectedRows();
	cout << "Number of selected events: " << N0 << endl;
	TH1F * his = c->GetHistogram();
	TF1 *fun = new TF1(fun_name, &DoubleExpCrystalBall, SCALE*(xmin-shift),  SCALE*(xmax-shift), 11);
	//gaus
	fun->SetParName(0, "N");
	fun->SetParameter(0, his->GetMaximum());
	fun->SetParName(1, "mean");
	fun->SetParameter(1, (3.09664-shift)*SCALE);
	fun->SetParName(2, "sigma");
	fun->SetParameter(2, 0.001237*SCALE);

	//left exp tail 
	fun->SetParName(3, "Ael");
	fun->SetParameter(3, 1.396);
	fun->SetParLimits(3, 0, 100);
	fun->FixParameter(3, 0);

	//left CB
	fun->SetParName(4, "Al-Ael");
	fun->SetParameter(4, 1.29);
	fun->SetParLimits(4, 0, 100);

	fun->SetParName(5, "Nl");
	fun->SetParameter(5, 2.162);
	fun->SetParLimits(5, 0, 1000);
	
	//right exp tail
	fun->SetParName(6, "Aer");
	fun->SetParameter(6, 1.17);
	fun->SetParLimits(6, 0, 100);

	//right CB
	fun->SetParName(7, "Ar-Aer");
	fun->SetParameter(7, 2.274);
	fun->SetParLimits(7, 0, 100);

	fun->SetParName(8, "Nr");
	fun->SetParameter(8, 1.991);
	fun->SetParLimits(8, 0, 1000);

	//background
	fun->SetParName(9, "BG0");
	fun->SetParameter(9, 0);
	fun->SetParLimits(9, 0, his->GetMaximum()*0.5);
	fun->FixParameter(9, 0);

	fun->SetParName(10, "BG1");
	fun->SetParameter(10, 0);
	fun->FixParameter(10, 0);

	fun->SetLineColor(kRed);
	fun->SetLineWidth(2);
	gStyle->SetOptFit(1111);
	//char  fit_name[1024];
	//sprintf(fit_name,  "pol1(0)+%s(2)", fun_name);
	//his->Fit(fit_name, fitopt, gopt);
	//his->Fit(fit_name, fitopt, gopt);
	his->Fit(fun_name, fitopt, gopt);
	his->Fit(fun_name, fitopt, gopt);
	char energy_unit_buf[1024];
	if(SCALE==1e6) sprintf(energy_unit_buf, "keV");
	if(SCALE==1e3) sprintf(energy_unit_buf, "MeV");
	if(SCALE==1)   sprintf(energy_unit_buf, "MeV");
	char title_buf[1024];
	if(shift !=0) sprintf(title_buf, "M_{rec}(#pi^{+}#pi^{-}) - %1.3f,  %s", shift*SCALE,  energy_unit_buf);
	else sprintf(title_buf, "M_{rec}(#pi^{+}#pi^{-}),  %s",  energy_unit_buf);
	his->GetXaxis()->SetTitle(title_buf);
	//c->Fit(fun_name, selection_buf,cut,fitopt, "E");
	double N = fun->GetParameter(0);
	double dN = fun->GetParError(0);
	double bg = fun->GetParameter(9);
	double dbg = fun->GetParError(9);
	double bg1 = fun->GetParameter(10);
	double dbg1 = fun->GetParError(10);
	
	double Nbg2 = bg*Nbins*(1.0 + 0.5*bg1*(xmax-xmin)*SCALE);
	double dNbg2 = Nbins*sqrt( pow( dbg*(1.0 + 0.5*bg1*(xmax+xmin-2*fun->GetParameter(1)-2*shift)*SCALE) , 2.0)   
			+ pow( dbg1*bg*0.5*(xmax-xmin)*SCALE , 2.0));


	double Nbg = Nbins*bg;
	double dNbg = Nbins*dbg;
	double Nkk = N0 - Nbins*bg;
	double dNkk = sqrt( N0 + dNbg*dNbg);
	cout << "Number of KK events: " <<  Nkk  << " +- " << dNkk <<   " ( +- " << dNkk/Nkk*100 << " % )" << endl;
	cout << "Number of background events: " << Nbg << " +- " <<  dNbg << " ( " << Nbg/N0*100 << " +- " << dNbg/N0*100 << " % )" << endl;
	cout << "Number of background events(2): " << Nbg2 << " +- " <<  dNbg2 << " ( " << Nbg2/N0*100 << " +- " << dNbg2/N0*100 << " % )" << endl;
	return fun;
}

TF1 *  fitE2CB2Norm(const char * dir="mcjpkk2009", const char * channel="KK",  int Nbins=1000,  double xmin=3.09,  double xmax=3.104,  const char * fitopt="RL",  const char * gopt="SAMES")
{
	//int r  = RANDOM->Integer(10000);
	int r = 0;
	const double SCALE=1e3;
	char fun_name[1024];
	char his_name[1024];
	sprintf(fun_name, "fun_2E2CBN_%s_%d",dir,  r);
	sprintf(his_name, "his_2E2CBN_%s_%d",dir,  r);
	r++;
	double shift =3.097;
	TChain * c = load_tree(dir);
	char cut_buf[1024];
	double range = xmax -xmin;
	sprintf(cut_buf, "%s && Mrec-%f> %f && Mrec-%f < %f", channel,  shift, xmin- shift,  shift,  xmax-shift);
	TCut cut = cut_buf;
	cut = cut && "kin_chi2<40 && pid_chi2<20";
	char  selection_buf[1024];

	sprintf(selection_buf, "%f*(Mrec-%f)>>%s(%d)", SCALE, shift, his_name,  Nbins);
	//c->Draw(selection_buf,cut,"goff");
	c->Draw(selection_buf,cut, "goff");
	double N0 = c->GetSelectedRows();
	cout << "Number of selected events: " << N0 << endl;
	TH1F * his = c->GetHistogram();
	TF1 *fun = new TF1(fun_name, &DoubleExpCrystalBallNorm, SCALE*(xmin-shift),  SCALE*(xmax-shift), 14);
	//gaus
	fun->SetParName(0, "Nsig");
	fun->SetParameter(0, N0);
	fun->SetParName(1, "mean");
	fun->SetParameter(1, (3.09664-shift)*SCALE);
	fun->SetParName(2, "sigma");
	fun->SetParameter(2, 0.001237*SCALE);

	//left exp tail 
	fun->SetParName(3, "Ael");
	fun->SetParameter(3, 1.396);
	fun->SetParLimits(3, 0, 100);
	fun->FixParameter(3, 0);

	//left CB
	fun->SetParName(4, "Al-Ael");
	fun->SetParameter(4, 1.29);
	fun->SetParLimits(4, 0, 100);

	fun->SetParName(5, "Nl");
	fun->SetParameter(5, 2.162);
	fun->SetParLimits(5, 1, 1000);
	
	//right exp tail
	fun->SetParName(6, "Aer");
	fun->SetParameter(6, 1.17);
	fun->SetParLimits(6, 0, 100);

	//right CB
	fun->SetParName(7, "Ar-Aer");
	fun->SetParameter(7, 2.274);
	fun->SetParLimits(7, 0, 100);

	fun->SetParName(8, "Nr");
	fun->SetParameter(8, 1.991);
	fun->SetParLimits(8, 1, 1000);

	//background
	fun->SetParName(9, "Nbg");
	fun->SetParameter(9, 0);
	fun->SetParLimits(9, 0, N0);
	//fun->FixParameter(9, N0);

	fun->SetParName(10, "kbg");
	fun->SetParameter(10, 0);
	//fun->FixParameter(10, 0);
	

	fun->SetParName(11, "xmin");
	fun->SetParName(12, "xmax");
	fun->FixParameter(11, (xmin-shift)*SCALE);
	fun->FixParameter(12, (xmax-shift)*SCALE);
	fun->SetParName(13, "Nbins");
	fun->FixParameter(13, Nbins);

	fun->SetLineColor(kRed);
	fun->SetLineWidth(2);
	gStyle->SetOptFit(1111);
	//char  fit_name[1024];
	//sprintf(fit_name,  "pol1(0)+%s(2)", fun_name);
	//his->Fit(fit_name, fitopt, gopt);
	//his->Fit(fit_name, fitopt, gopt);
	his->Fit(fun_name, fitopt, gopt);
	his->Fit(fun_name, fitopt, gopt);
	char energy_unit_buf[1024];
	if(SCALE==1e6) sprintf(energy_unit_buf, "keV");
	if(SCALE==1e3) sprintf(energy_unit_buf, "MeV");
	if(SCALE==1)   sprintf(energy_unit_buf, "MeV");
	char title_buf[1024];
	if(shift !=0) sprintf(title_buf, "M_{rec}(#pi^{+}#pi^{-}) - %1.3f,  %s", shift*SCALE,  energy_unit_buf);
	else sprintf(title_buf, "M_{rec}(#pi^{+}#pi^{-}),  %s",  energy_unit_buf);
	his->GetXaxis()->SetTitle(title_buf);
	cout<< "Number of events in fun: " << fun->Integral((xmin-shift)*SCALE, (xmax-shift)*SCALE)*Nbins/(xmax-xmin)/SCALE << endl;
	//c->Fit(fun_name, selection_buf,cut,fitopt, "E");
	double N = fun->GetParameter(0);
	double dN = fun->GetParError(0);
	double bg = fun->GetParameter(9);
	double dbg = fun->GetParError(9);
	double bg1 = fun->GetParameter(10);
	double dbg1 = fun->GetParError(10);
	
	double Nbg2 = bg*Nbins*(1.0 + 0.5*bg1*(xmax-xmin)*SCALE);
	double dNbg2 = Nbins*sqrt( pow( dbg*(1.0 + 0.5*bg1*(xmax-xmin)*SCALE) , 2.0)   
			+ pow( dbg1*bg*0.5*(xmax-xmin)*SCALE , 2.0));


	double Nbg = Nbins*bg;
	double dNbg = Nbins*dbg;
	double Nkk = N0 - Nbins*bg;
	double dNkk = sqrt( N0 + dNbg*dNbg);
	cout << "Number of KK events: " <<  Nkk  << " +- " << dNkk <<   " ( +- " << dNkk/Nkk*100 << " % )" << endl;
	cout << "Number of background events: " << Nbg << " +- " <<  dNbg << " ( " << Nbg/N0*100 << " +- " << dNbg/N0*100 << " % )" << endl;
	cout << "Number of background events(2): " << Nbg2 << " +- " <<  dNbg2 << " ( " << Nbg2/N0*100 << " +- " << dNbg2/N0*100 << " % )" << endl;
	return fun;
}



//Сравнение монте карло и данных в при подгонке массы отдачи
void fitMrec_mcdatacmp(const char * dir_mc="mc09", const char * dir_data="2009",  const char * channel="KK",  int Nbins=1000,  double xmin=3.08,  double xmax=3.114,  const char * fitopt="ERLLI")
{
	TCanvas * c = new TCanvas;
	fitE2CB2(dir_data,channel, Nbins,xmin,xmax,fitopt, "");
	fitE2CB2(dir_mc,  channel, Nbins,xmin,xmax,fitopt,"sames");
}


void testE2CB2(const char * dir="2009", int Nbins=1000,  double xmin=3.08,  double xmax=3.114,  const char * fitopt="R")
{
	double range=xmax-xmin;
	double shift = 0;
	TF1 * f =  fitE2CB2("mcjpkk2009/", Nbins,  xmin+range*0.01, xmax-range*0.01, fitopt);
	TF1 * f2 = new TF1("test_fun", &DoubleExpCrystalBall, xmin,  xmax, 10);
	const char * his_name = "test_hfit_e2cb2";
	TChain * c = load_tree(dir);
	char cut_buf[1024];
	sprintf(cut_buf, "Mrec-%f> %f && Mrec-%f < %f", shift, xmin- shift,  shift,  xmax-shift);
	TCut cut = cut_buf;
	cut = cut && "KK && kin_chi2<40 && pid_chi2<20";
	char  selection_buf[1024];
	sprintf(selection_buf, "Mrec-%f>>%s(%d)", shift, his_name,  Nbins);
	c->Draw(selection_buf,cut,"E");
	TH1F * his = c->GetHistogram();
	f2->SetParameter(0,his->GetMaximum());
	for(int i=1;i<10;i++) f2->SetParameter(i, f->GetParameter(i));
	for(int i=1;i<9;i++) f2->FixParameter(i, f->GetParameter(i));
	f2->SetParLimits(9, 0, his->GetMaximum());
	his->Fit("test_fun",  fitopt);
}



TF1 *  fitE2CB2Norm2(const char * dir="mcjpkk2009", const char * channel="KK",  int Nbins=1000,  double xmin=3.09,  double xmax=3.104,  const char * fitopt="RL",  const char * gopt="")
{
	int r = 0;
	const double SCALE=1e3;
	char fun_name[1024];
	char his_name[1024];
	sprintf(fun_name, "fun_2E2CBN_%s_%d",dir,  r);
	sprintf(his_name, "his_2E2CBN_%s_%d",dir,  r);
	r++;
	double shift =3.097;
	TChain * c = load_tree(dir);
	char cut_buf[1024];
	double range = xmax -xmin;
	sprintf(cut_buf, "%s && Mrec-%f> %f && Mrec-%f < %f", channel,  shift, xmin- shift,  shift,  xmax-shift);
	TCut cut = cut_buf;
	cut = cut && "kin_chi2<40 && pid_chi2<20";
	char  selection_buf[1024];

	sprintf(selection_buf, "%f*(Mrec-%f)>>%s(%d, %f, %f)", SCALE, shift, his_name,  Nbins,  (xmin-shift)*SCALE,  (xmax-shift)*SCALE);
	//c->Draw(selection_buf,cut,"goff");
	c->Draw(selection_buf,cut, "goff");
	double N0 = c->GetSelectedRows();
	cout << "Number of selected events: " << N0 << endl;
	TH1F * his = c->GetHistogram();
	TF1 *fun = new TF1(fun_name, &ModifiedDoubleCrystalBall, SCALE*(xmin-shift),  SCALE*(xmax-shift), 14);
	//gaus
	fun->SetParName(0, "Nsig");
	fun->SetParameter(0, N0);
	fun->SetParLimits(0, 0, N0);

	fun->SetParName(1, "mean");
	fun->SetParameter(1, (3.09664-shift)*SCALE);
	fun->SetParName(2, "sigma");
	//fun->SetParLimits(2, 0, 3*(xmax-shift)*SCALE);
	fun->SetParameter(2, 0.0013*SCALE);

	//left exp tail 
	fun->SetParName(3, "bL");
	fun->SetParameter(3, 1);
	//fun->SetParLimits(3, 0, 100);

	//left CB
	fun->SetParName(4, "aL-bL");
	fun->SetParameter(4, 0);
	fun->SetParLimits(4, 0, 100);
	//fun->FixParameter(4, 0);

	fun->SetParName(5, "nL");
	fun->SetParameter(5, 2.162);
	//fun->SetParLimits(5, 0, 1000);
	
	//right exp tail
	fun->SetParName(6, "bR");
	fun->SetParameter(6, 1);
	//fun->SetParLimits(6, 0, 100);

	//right CB
	fun->SetParName(7, "aR-bR");
	fun->SetParameter(7, 2);
	//fun->SetParLimits(7, 0, 100);

	fun->SetParName(8, "nR");
	fun->SetParameter(8, 1);
	//fun->SetParLimits(8, 0, 1000);

	//total number of envets in histogram
	fun->SetParName(9, "N0");
	fun->SetParameter(9, N0);
	fun->SetParLimits(9, 0, N0);
	fun->FixParameter(9, N0);

	//background
	fun->SetParName(10, "bgslope");
	fun->SetParameter(10, 0);
	//fun->FixParameter(10, 0);
	

	fun->SetParName(11, "xmin");
	fun->SetParName(12, "xmax");
	fun->FixParameter(11, (xmin-shift)*SCALE);
	fun->FixParameter(12, (xmax-shift)*SCALE);
	fun->SetParName(13, "Nbins");
	fun->FixParameter(13, Nbins);

	fun->SetLineColor(kRed);
	fun->SetLineWidth(2);
	gStyle->SetOptFit(1111);
	//char  fit_name[1024];
	//sprintf(fit_name,  "pol1(0)+%s(2)", fun_name);
	//his->Fit(fit_name, fitopt, gopt);
	//his->Fit(fit_name, fitopt, gopt);
	his->Fit(fun_name, fitopt, gopt);
	//his->Fit(fun_name, fitopt, gopt);
	//his->Fit(fun_name, fitopt, gopt);
	char energy_unit_buf[1024];
	if(SCALE==1e6) sprintf(energy_unit_buf, "keV");
	if(SCALE==1e3) sprintf(energy_unit_buf, "MeV");
	if(SCALE==1)   sprintf(energy_unit_buf, "MeV");
	char title_buf[1024];
	if(shift !=0) sprintf(title_buf, "M_{rec}(#pi^{+}#pi^{-}) - %1.3f,  %s", shift*SCALE,  energy_unit_buf);
	else sprintf(title_buf, "M_{rec}(#pi^{+}#pi^{-}),  %s",  energy_unit_buf);
	his->GetXaxis()->SetTitle(title_buf);
	cout<< "Number of events in fun: " << fun->Integral((xmin-shift)*SCALE, (xmax-shift)*SCALE)*Nbins/(xmax-xmin)/SCALE << endl;
	//c->Fit(fun_name, selection_buf,cut,fitopt, "E");
	double Nsig = fun->GetParameter(0);
	double dNsig = fun->GetParError(0);
	double Nbg = N0 - Nsig;
	double dNbg = dNsig;
	//cout.precission(15);
	cout << "Total number of selected events: " << N0 << endl;
	cout << "Number of signal events: " <<  long(Nsig)  << " +- " << long(dNsig) <<   " ( +- " << dNsig/Nsig*100 << " % )" << endl;
	cout << "Number of background events: " << Nbg << " +- " <<  dNbg << " ( " << Nbg/N0*100 << " +- " << dNbg/N0*100 << " % )" << endl;

	TH1F * h2  = new TH1F("h2", "h2", Nbins, (xmin-shift)*SCALE, (xmax-shift)*SCALE);
	cout << "xmin = " << (xmin-shift)*SCALE << endl;
	cout << "xmax = " << (xmax-shift)*SCALE << endl;

	double difference=0;
	for(int i=0; i<Nbins;i++)
	{
		double x = his->GetBinCenter(i);
		double f = fun->Eval(x);
		double n = his->GetBinContent(i);
		difference += (f-n);
		h2->SetBinContent(i, f-n);
	}
	cout << "Number difference between function and histogram: " << difference << endl;

	new TCanvas;
	h2->Draw();
	cout << "h2 bin number: " << h2->GetNbinsX() << "  Nbins = " << Nbins << endl;
	return fun;
}

void  fitE2CB2Norm3(const char * dir="mcjpkk2009", const char * channel="KK",  int Nbins=1000,  double xmin=3.09,  double xmax=3.104,  const char * fitopt="RL",  const char * gopt="")
{
	int r = 0;
	const double SCALE=1e3;
	char fun_name[1024];
	char his_name[1024];
	sprintf(fun_name, "fun_2E2CBN_%s_%d",dir,  r);
	sprintf(his_name, "his_2E2CBN_%s_%d",dir,  r);
	r++;
	double shift =3.097;
	TChain * c = load_tree(dir);
	char cut_buf[1024];
	double range = xmax -xmin;
	sprintf(cut_buf, "%s && Mrec-%f> %f && Mrec-%f < %f", channel,  shift, xmin- shift,  shift,  xmax-shift);
	TCut cut = cut_buf;
	cut = cut && "kin_chi2<40 && pid_chi2<20";
	char  selection_buf[1024];

	sprintf(selection_buf, "%f*(Mrec-%f)>>%s(%d, %f, %f)", SCALE, shift, his_name,  Nbins,  (xmin-shift)*SCALE,  (xmax-shift)*SCALE);
	//c->Draw(selection_buf,cut,"goff");
	c->Draw(selection_buf,cut);
	double N0 = c->GetSelectedRows();
	cout << "Number of selected events: " << N0 << endl;
	TH1F * his = c->GetHistogram();
	const  double * par = Fit(his);
	//TF1 *fun = new TF1(fun_name, &ModifiedDoubleCrystalBall, SCALE*(xmin-shift),  SCALE*(xmax-shift), 14);
	//fun->SetParameters(par);
	//fun->SetLineColor(kRed);
	//fun->Draw("same");
	CrystalBall cb(his);
	//ROOT::Math::ParamFunctor pf;
	//pf.SetFunction(cb,CrystalBall::Eval);
	TF1 * fun2 = new TF1("test_fun",  &cb, CrystalBall::Eval,   SCALE*(xmin-shift),  SCALE*(xmax-shift), 10, "CrystalBall", "Eval");
	////TF1 * fun2 = new TF1("test_fun",  pf,   SCALE*(xmin-shift),  SCALE*(xmax-shift), 14);
	fun2->SetLineColor(kBlue);
	fun2->SetParameters(par);
	fun2->Draw("same");
	//return fun;
}

void  fitE2CB2Norm4(const char * dir="mcjpkk2009", const char * channel="KK",  int Nbins=1000,  double xmin=3.09,  double xmax=3.104,  const char * fitopt="RL",  const char * gopt="")
{
	int r = 0;
	const double SCALE=1e3;
	char fun_name[1024];
	char his_name[1024];
	sprintf(fun_name, "fun_2E2CBN_%s_%d",dir,  r);
	sprintf(his_name, "his_2E2CBN_%s_%d",dir,  r);
	r++;
	double shift =3.097;
	TChain * c = load_tree(dir);
	char cut_buf[1024];
	double range = xmax -xmin;
	sprintf(cut_buf, "%s && Mrec-%f> %f && Mrec-%f < %f", channel,  shift, xmin- shift,  shift,  xmax-shift);
	TCut cut = cut_buf;
	cut = cut && "kin_chi2<40 && pid_chi2<20";
	char  selection_buf[1024];

	sprintf(selection_buf, "%f*(Mrec-%f)>>%s(%d, %f, %f)", SCALE, shift, his_name,  Nbins,  (xmin-shift)*SCALE,  (xmax-shift)*SCALE);
	//c->Draw(selection_buf,cut,"goff");
	c->Draw(selection_buf,cut);
	double N0 = c->GetSelectedRows();
	cout << "Number of selected events: " << N0 << endl;
	TH1F * his = c->GetHistogram();
	CrystalBallFitter cbf(his);
	cbf.Fit();
	//const  double * par = Fit(his);
	////TF1 *fun = new TF1(fun_name, &ModifiedDoubleCrystalBall, SCALE*(xmin-shift),  SCALE*(xmax-shift), 14);
	////fun->SetParameters(par);
	////fun->SetLineColor(kRed);
	////fun->Draw("same");
	//CrystalBall cb(his);
	////ROOT::Math::ParamFunctor pf;
	////pf.SetFunction(cb,CrystalBall::Eval);
	//TF1 * fun2 = new TF1("test_fun",  &cb, CrystalBall::Eval,   SCALE*(xmin-shift),  SCALE*(xmax-shift), 10, "CrystalBall", "Eval");
	//////TF1 * fun2 = new TF1("test_fun",  pf,   SCALE*(xmin-shift),  SCALE*(xmax-shift), 14);
	//fun2->SetLineColor(kBlue);
	//fun2->SetParameters(par);
	//fun2->Draw("same");
	//return fun;
}

void  fitE2CB2Norm5(const char * dir="mcjpkk2009", const char * channel="KK",  int Nbins=1000,  double xmin=3.09,  double xmax=3.104,  const char * fitopt="RL",  const char * gopt="")
{
	int r = 0;
	const double SCALE=1e3;
	char fun_name[1024];
	char his_name[1024];
	sprintf(fun_name, "fun_2E2CBN_%s_%d",dir,  r);
	sprintf(his_name, "his_2E2CBN_%s_%d",dir,  r);
	r++;
	double shift =3.097;
	TChain * c = load_tree(dir);
	char cut_buf[1024];
	double range = xmax -xmin;
	sprintf(cut_buf, "%s && Mrec-%f> %f && Mrec-%f < %f", channel,  shift, xmin- shift,  shift,  xmax-shift);
	TCut cut = cut_buf;
	cut = cut && "kin_chi2<40 && pid_chi2<20";
	char  selection_buf[1024];

	sprintf(selection_buf, "%f*(Mrec-%f)>>%s(%d, %f, %f)", SCALE, shift, his_name,  Nbins,  (xmin-shift)*SCALE,  (xmax-shift)*SCALE);
	//c->Draw(selection_buf,cut,"goff");
	c->Draw(selection_buf,cut, "E");
	double N0 = c->GetSelectedRows();
	cout << "Number of selected events: " << N0 << endl;
	TH1F * his = c->GetHistogram();
	Fit2(his);
	//const  double * par = Fit(his);
	////TF1 *fun = new TF1(fun_name, &ModifiedDoubleCrystalBall, SCALE*(xmin-shift),  SCALE*(xmax-shift), 14);
	////fun->SetParameters(par);
	////fun->SetLineColor(kRed);
	////fun->Draw("same");
	//CrystalBall cb(his);
	////ROOT::Math::ParamFunctor pf;
	////pf.SetFunction(cb,CrystalBall::Eval);
	//TF1 * fun2 = new TF1("test_fun",  &cb, CrystalBall::Eval,   SCALE*(xmin-shift),  SCALE*(xmax-shift), 10, "CrystalBall", "Eval");
	//////TF1 * fun2 = new TF1("test_fun",  pf,   SCALE*(xmin-shift),  SCALE*(xmax-shift), 14);
	//fun2->SetLineColor(kBlue);
	//fun2->SetParameters(par);
	//fun2->Draw("same");
	//return fun;
}

void load(void)
{
	gROOT->Reset();
  gSystem->AddIncludePath("-I$HOME/work");
	gSystem->CompileMacro("CrystalBall.cpp","kO","","/tmp");
	gSystem->CompileMacro("analize.C","kO","","/tmp");
}


void analize(const char * dir, Long64_t N=5e4)
{
	TChain * c = load_tree(dir);
  analize an;
  c->Process(&an,"",N);
}

