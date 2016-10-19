/*
 * =====================================================================================
 *
 *       Filename:  fit.h
 *
 *    Description:  Fit to 
 *
 *        Version:  1.0
 *        Created:  19.10.2016 23:25:39
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */


#include <RooSimultaneous.h>
#include <RooFitResult.h>
#include <RooChi2Var.h>
#include <RooNLLVar.h>
#include <RooRealVar.h>
#include <RooPolynomial.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooCategory.h>
#include <RooAddPdf.h>
#include <TCanvas.h>
#include <RooPlot.h>
#include <TAxis.h>
#include <RooBukinPdf.h>

#include "RooMcbPdf.h"

using namespace RooFit;
using namespace std;

#include "RooMcbPdf.h"
extern void fit(TH1 * hisKK, TH1 * hisUU);
