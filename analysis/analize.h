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

#include "../../mctop/McTopo.h"

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

#ifdef analize_cxx
void analize::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("time", &time, &b_time);
   fChain->SetBranchAddress("ngtrack", &ngtrack, &b_ngtrack);
   fChain->SetBranchAddress("ngntrack", &ngntrack, &b_ngntrack);
   fChain->SetBranchAddress("nptrack", &nptrack, &b_nptrack);
   fChain->SetBranchAddress("nntrack", &nntrack, &b_nntrack);
   fChain->SetBranchAddress("nppions", &nppions, &b_nppions);
   fChain->SetBranchAddress("nnpions", &nnpions, &b_nnpions);
   fChain->SetBranchAddress("npion_pairs", &npion_pairs, &b_npion_pairs);
   fChain->SetBranchAddress("channel", &channel, &b_channel);
   fChain->SetBranchAddress("KK", &KK, &b_KK);
   fChain->SetBranchAddress("uu", &uu, &b_uu);
   fChain->SetBranchAddress("Mrec", &Mrec, &b_Mrec);
   fChain->SetBranchAddress("Minv", &Minv, &b_Minv);
   fChain->SetBranchAddress("M012", &M012, &b_M012);
   fChain->SetBranchAddress("M013", &M013, &b_M013);
   fChain->SetBranchAddress("M03", &M03, &b_M03);
   fChain->SetBranchAddress("M12", &M12, &b_M12);
   fChain->SetBranchAddress("M01", &M01, &b_M01);
   fChain->SetBranchAddress("kin_chi2", &kin_chi2, &b_kin_chi2);
   fChain->SetBranchAddress("pid_chi2", &pid_chi2, &b_pid_chi2);
   fChain->SetBranchAddress("ntrack", &ntrack, &b_ntrack);
   fChain->SetBranchAddress("q", q, &b_q);
   fChain->SetBranchAddress("E", E, &b_E);
   fChain->SetBranchAddress("p", p, &b_p);
   fChain->SetBranchAddress("px", px, &b_px);
   fChain->SetBranchAddress("py", py, &b_py);
   fChain->SetBranchAddress("pz", pz, &b_pz);
   fChain->SetBranchAddress("pt", pt, &b_pt);
   fChain->SetBranchAddress("theta", theta, &b_theta);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("x", x, &b_x);
   fChain->SetBranchAddress("y", y, &b_y);
   fChain->SetBranchAddress("z", z, &b_z);
   fChain->SetBranchAddress("r", r, &b_r);
   fChain->SetBranchAddress("vxy", vxy, &b_vxy);
   fChain->SetBranchAddress("vz", vz, &b_vz);
   fChain->SetBranchAddress("vphi", vphi, &b_vphi);
   mctp.Init(fChain->GetFriend("mctopo"));//->GetFriend("mctop");
}

Bool_t analize::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef analize_cxx
