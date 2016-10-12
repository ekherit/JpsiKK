#ifndef MODIFIEDCRYSTALBALL
#define MODIFIEDCRYSTALBALL

#include <memory>
#include <vector>

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class RooMcbPdf : public RooAbsPdf {
public:
  RooMcbPdf() {} ;
  RooMcbPdf(const char *name, const char *title, 
	      RooAbsReal& _X, //invariant mass 
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
			);

  RooMcbPdf(const char *name, const char *title, 
	      RooAbsReal& _X, //invariant mass 
				RooAbsReal & _Sigma,   //Common sigma
				RooAbsReal & _s1,   //staple
				RooAbsReal & _s2,   //staple 
				RooAbsReal & _s3,   //staple 
				RooAbsReal & _s4,   //staple 
				RooAbsReal & _s5,   //staple 
				RooAbsReal & _n1,   //staple 
				RooAbsReal & _n2   //staple 
			);


  RooMcbPdf(const RooMcbPdf& other,const char* name=0) ;

  virtual TObject* clone(const char* newname) const { return new RooMcbPdf(*this,newname); }
  inline virtual ~RooMcbPdf() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:
  RooRealProxy fX;
  RooRealProxy fSigma;
  RooRealProxy fStaple1;
  RooRealProxy fStaple2;
  RooRealProxy fStaple3;
  RooRealProxy fStaple4;
  RooRealProxy fStaple5;
  RooRealProxy fStaple6;
  RooRealProxy fStaple7;
  RooRealProxy fN1;
  RooRealProxy fN2;

  Double_t evaluate() const ;
//   void initGenerator();
	
	virtual Bool_t selfNormalized() const
	{
		return kFALSE;
	}

private:
	void initArgs(void) const;
	mutable double mean;
	mutable int mean_index;
	mutable std::vector<double> staple; 
	mutable double n[2]; //power of tail

	mutable double A[2]; //first exp tail amplitude
	mutable double B[2]; //power tail amplitude
	mutable double C[2]; //second exp tail amplitude

//	mutable double xb[2]; //handy shift for power tail
	mutable double ka[2]; //first exp tail k
	mutable double kb[2]; //koef for power
	mutable double kc[2]; //second exp tail k


	mutable double a[2]; 
	mutable double b[2]; 
	mutable double c[2]; 
	
	enum { TYPE_GPE,  TYPE_GEPE,  } fType;
//  ClassDef(RooMcbPdf,1) 
};
 
class RooBgPdf : public RooAbsPdf 
{
public:
  RooBgPdf() {} ;
  RooBgPdf(const char *name, const char *title, 
	      RooAbsReal& _X, //invariant mass 
				RooAbsReal & _B, 
        double x1, //xmin
        double x2  //xmax
			);


  RooBgPdf(const RooBgPdf& other,const char* name=0) ;

  virtual TObject* clone(const char* newname) const { return new RooBgPdf(*this,newname); }
  inline virtual ~RooBgPdf() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:
  RooRealProxy fX;
  RooRealProxy fB;
  double fXmax,fXmin;

  Double_t evaluate() const ;
	
	virtual Bool_t selfNormalized() const
	{
		return kFALSE;
	}

private:
//ClassDef(RooBgPdf,1) 
};

class RooMcb2Pdf : public RooAbsPdf
{
public:
  RooMcb2Pdf() {} ;
  RooMcb2Pdf(const char *name, const char *title, 
	      RooAbsReal & _X, //invariant mass 
				RooAbsReal & _mean,    //position of the peack
				RooAbsReal & _Sigma,   //Common sigma
				RooAbsReal & _staple1, 
				RooAbsReal & _staple2, 
				RooAbsReal & _staple3, 
				RooAbsReal & _staple4, 
				RooAbsReal & _staple5, 
				RooAbsReal & _staple6, 
				RooAbsReal & _staple7, 
				RooAbsReal & _staple8, 
				RooAbsReal & _n1, 
				RooAbsReal & _n2, 
        double xmin,
        double xmax
			);

  RooMcb2Pdf(const RooMcb2Pdf& other,const char* name=0) ;

  virtual TObject* clone(const char* newname) const { return new RooMcb2Pdf(*this,newname); }

  inline virtual ~RooMcb2Pdf() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:
  RooRealProxy fX;
  RooRealProxy fMean;
  RooRealProxy fSigma;
  RooRealProxy fStaple1;
  RooRealProxy fStaple2;
  RooRealProxy fStaple3;
  RooRealProxy fStaple4;
  RooRealProxy fStaple5;
  RooRealProxy fStaple6;
  RooRealProxy fStaple7;
  RooRealProxy fStaple8;
  RooRealProxy fN1;
  RooRealProxy fN2;
  double fXmin;
  double fXmax;

  Double_t evaluate() const ;
	
	virtual Bool_t selfNormalized() const
	{
		return kFALSE;
	}

private:
	void initArgs(void) const;
	mutable double mean;
	mutable int mean_index;
	mutable std::vector<double> staple[2]; 
	mutable double n[2]; //power of tail

	mutable double A[2]; //first exp tail amplitude
	mutable double B[2]; //power tail amplitude
	mutable double C[2]; //second exp tail amplitude

//	mutable double xb[2]; //handy shift for power tail
	mutable double ka[2]; //first exp tail k
	mutable double kb[2]; //koef for power
	mutable double kc[2]; //second exp tail k


	mutable double a[2]; 
	mutable double b[2]; 
	mutable double c[2]; 
	
	enum { TYPE_GPE,  TYPE_GEPE,  } fType;
};
#endif
