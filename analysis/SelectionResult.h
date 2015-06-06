#ifndef IBN_JPSIKK_SELECTION_RESULT_H
#define IBN_JPSIKK_SELECTION_RESULT_H

#include <TH1F.h>
//#include "Event.h"
#include <string>
struct SelectionResult
{
  Long64_t NKK=0; //number of selected KK events
  Long64_t Nuu=0; //number of selected uu events
  Long64_t N0=0; //total number of events
  const double MJPSI_SHIFT=3.097;
  const double MSCALE = 1e3;
  double MRANGE;
  double PID_CHI2;
  double KIN_CHI2;
  double MAX_RECOIL_MASS;
  double MIN_RECOIL_MASS;
  int NBINS_UU=500;
  int NBINS_KK=200;
  std::string output_file_name;

  TH1F * hMrecKK=0;
  TH1F * hMrecUU=0;
  TH1F * hpid_chi2KK=0;
  TH1F * hpid_chi2UU=0;
  TH1F * hkin_chi2KK=0;
  TH1F * hkin_chi2UU=0;
  TTree * event_treeKK=0;
  TTree * mctopo_treeKK=0;
  TTree * event_treeUU=0;
  TTree * mctopo_treeUU=0;
  TFile * file;

  SelectionResult(void)
  {
    event_treeKK=0;
    mctopo_treeKK=0;
    event_treeUU=0;
    mctopo_treeUU=0;
    output_file_name="analize_result.root";
  }

  ~SelectionResult(void)
  {
    delete file;
  }

  void Init(void)
  {
    MRANGE=0.09;
    PID_CHI2=20;
    KIN_CHI2=40;
    MAX_RECOIL_MASS=MJPSI_SHIFT+MRANGE*0.5;
    MIN_RECOIL_MASS=MJPSI_SHIFT-MRANGE*0.5;

    NKK=0;
    Nuu=0;
    N0=0;
    if(hMrecKK!=0) delete hMrecKK;
    if(hMrecUU!=0) delete hMrecUU;
    hMrecKK = new TH1F("hMrecKK","#pi^{+}#pi^{-} recoil mass for KK channel",NBINS_KK,mshift(MIN_RECOIL_MASS),mshift(MAX_RECOIL_MASS));
    hMrecUU = new TH1F("hMrecUU","#pi^{+}#pi^{-} recoil mass for uu channel",NBINS_UU,mshift(MIN_RECOIL_MASS),mshift(MAX_RECOIL_MASS));
    hpid_chi2KK = new TH1F("hpid_chi2KK","pid chi2 KK channel",300,0, PID_CHI2);
    hpid_chi2UU = new TH1F("hpid_chi2UU","pid chi2 UU channel",300,0, PID_CHI2);
    hkin_chi2KK = new TH1F("hkin_chi2KK","kin chi2 KK channel",300,0, KIN_CHI2);
    hkin_chi2UU = new TH1F("hkin_chi2UU","kin chi2 UU channel",300,0, KIN_CHI2);
  }

  double mshift(double m) const
  {
    return MSCALE*(m-MJPSI_SHIFT);
  }
};
#endif
