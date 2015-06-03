#ifndef IBN_JPSIKK_SELECTION_RESULT_H
#define IBN_JPSIKK_SELECTION_RESULT_H

struct SelectionResult
{
  Long64_t NKK=0; //number of selected KK events
  Long64_t Nuu=0; //number of selected uu events
  Long64_t N0=0; //total number of events

  TH1F * hMrecKK=0;
  TH1F * hMrecUU=0;

  Init(void)
  {
    NKK=0;
    Nuu=0;
    N0=0;
    if(hMrecKK!=0) delete hMrecKK;
    if(hMrecUU!=0) delete hMrecUU;
    hMrecKK = new TH1F("hMrecKK","#pi^{+}#pi^{-} recoil mass for KK channel",200,3.08,3.105);
    hMrecUU = new TH1F("hMrecUU","#pi^{+}#pi^{-} recoil mass for uu channel",200,3.08,3.105);
  }

};
#endif
