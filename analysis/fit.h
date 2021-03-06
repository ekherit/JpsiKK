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
class TTree;

#include "RooMcbPdf.h"

extern void fit(TH1 * his);
extern void fit(TH1 * hisKK, TH1 * hisUU);
extern void fit(std::list<TH1*> & , std::list<TTree*> & , bool = false);
extern void draw(TTree *);

extern bool OPT_NOBGSLOPE; //no slope for the background
extern bool OPT_NOBG; //no background
extern bool OPT_NOGAUSRAD; //no gaus rad
extern std::string OPT_PARAM_CONFIG_FILE; 
extern bool OPT_SEPARATE_MREC; //separate Mrec for each channel
extern std::string OPT_FIT_METHOD;
extern std::string OPT_FIT_RESULT_FILE_NAME;
extern int OPT_NGRAD; //number of gauses used for radiative effects
extern bool OPT_SHOW_FIT_RESULT; 
extern bool OPT_FIT_INTEGRATE;


#endif
