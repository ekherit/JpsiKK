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

#ifndef IBN_JPSIKK_FIT_H
#define IBN_JPSIKK_FIT_H

#include <list>
#include <string>

using namespace std;

class TH1;

#include "RooMcbPdf.h"

extern void fit(TH1 * his);
extern void fit(TH1 * hisKK, TH1 * hisUU);
extern void fit(std::list<TH1*> & );

extern bool OPT_NOBGSLOPE; //no slope for the background
extern bool OPT_NOBG; //no background
extern bool OPT_NOGAUSRAD; //no gaus rad
extern std::string OPT_PARAM_CONFIG_FILE; 

#endif
