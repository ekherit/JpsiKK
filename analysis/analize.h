//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jun  3 10:59:36 2015 by ROOT version 5.34/05
// from TTree event/Signal events pi+pi- K+K-, or pi+pi- mu+mu-
// found on file: mc09-0008348.root
//////////////////////////////////////////////////////////

#ifndef analize_h
#define analize_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

//#include "mctopo.h"

#include <mctopo/McTopo.h>

#include "SelectionResult.h"
class analize : public TSelector , public SelectionResult {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   //mctopo   mctp;
   McTopo  mctp;

   // Declaration of leaf types
   Int_t           run;
   Int_t           event;
   Int_t           time;
   Int_t           ngtrack;
   Int_t           ngntrack;
   Int_t           nptrack;
   Int_t           nntrack;
   Int_t           nppions;
   Int_t           nnpions;
   Int_t           npion_pairs;
   Int_t           channel;
   Int_t           KK;
   Int_t           uu;
   Double_t        Mrec;
   Double_t        Minv;
   Double_t        M012;
   Double_t        M013;
   Double_t        M03;
   Double_t        M12;
   Double_t        M01;
   Double_t        kin_chi2;
   Double_t        pid_chi2;
   Int_t           ntrack;
   Double_t        q[4];   //[ntrack]
   Double_t        E[4];   //[ntrack]
   Double_t        p[4];   //[ntrack]
   Double_t        px[4];   //[ntrack]
   Double_t        py[4];   //[ntrack]
   Double_t        pz[4];   //[ntrack]
   Double_t        pt[4];   //[ntrack]
   Double_t        theta[4];   //[ntrack]
   Double_t        phi[4];   //[ntrack]
   Double_t        x[4];   //[ntrack]
   Double_t        y[4];   //[ntrack]
   Double_t        z[4];   //[ntrack]
   Double_t        r[4];   //[ntrack]
   Double_t        vxy[4];   //[ntrack]
   Double_t        vz[4];   //[ntrack]
   Double_t        vphi[4];   //[ntrack]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_time;   //!
   TBranch        *b_ngtrack;   //!
   TBranch        *b_ngntrack;   //!
   TBranch        *b_nptrack;   //!
   TBranch        *b_nntrack;   //!
   TBranch        *b_nppions;   //!
   TBranch        *b_nnpions;   //!
   TBranch        *b_npion_pairs;   //!
   TBranch        *b_channel;   //!
   TBranch        *b_KK;   //!
   TBranch        *b_uu;   //!
   TBranch        *b_Mrec;   //!
   TBranch        *b_Minv;   //!
   TBranch        *b_M012;   //!
   TBranch        *b_M013;   //!
   TBranch        *b_M03;   //!
   TBranch        *b_M12;   //!
   TBranch        *b_M01;   //!
   TBranch        *b_kin_chi2;   //!
   TBranch        *b_pid_chi2;   //!
   TBranch        *b_ntrack;   //!
   TBranch        *b_q;   //!
   TBranch        *b_E;   //!
   TBranch        *b_p;   //!
   TBranch        *b_px;   //!
   TBranch        *b_py;   //!
   TBranch        *b_pz;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_z;   //!
   TBranch        *b_r;   //!
   TBranch        *b_vxy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_vphi;   //!

   analize(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~analize() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(analize,0);
};

#endif

