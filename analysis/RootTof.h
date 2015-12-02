//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec  2 20:01:56 2015 by ROOT version 6.04/10
// from TTree tof/Time of Flight
// found on file: sample.root
//////////////////////////////////////////////////////////

#ifndef RootTof_h
#define RootTof_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class RootTof {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           ntrack;
   Double_t        ID[4];   //[ntrack]
   Double_t        t[4];   //[ntrack]
   Double_t        dt[4];   //[ntrack]
   Double_t        t0[4];   //[ntrack]
   Double_t        chie[4];   //[ntrack]
   Double_t        chimu[4];   //[ntrack]
   Double_t        chipi[4];   //[ntrack]
   Double_t        chik[4];   //[ntrack]
   Double_t        chip[4];   //[ntrack]
   Double_t        beta[4];   //[ntrack]
   Double_t        te[4];   //[ntrack]
   Double_t        tmu[4];   //[ntrack]
   Double_t        tpi[4];   //[ntrack]
   Double_t        tk[4];   //[ntrack]
   Double_t        tp[4];   //[ntrack]

   // List of branches
   TBranch        *b_ntrack;   //!
   TBranch        *b_ID;   //!
   TBranch        *b_t;   //!
   TBranch        *b_dt;   //!
   TBranch        *b_t0;   //!
   TBranch        *b_chie;   //!
   TBranch        *b_chimu;   //!
   TBranch        *b_chipi;   //!
   TBranch        *b_chik;   //!
   TBranch        *b_chip;   //!
   TBranch        *b_beta;   //!
   TBranch        *b_te;   //!
   TBranch        *b_tmu;   //!
   TBranch        *b_tpi;   //!
   TBranch        *b_tk;   //!
   TBranch        *b_tp;   //!

   RootTof(TTree *tree=0);
   virtual ~RootTof();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef RootTof_cxx
RootTof::RootTof(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("sample.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("sample.root");
      }
      f->GetObject("tof",tree);

   }
   Init(tree);
}

RootTof::~RootTof()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t RootTof::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t RootTof::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void RootTof::Init(TTree *tree)
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
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ntrack", &ntrack, &b_ntrack);
   fChain->SetBranchAddress("ID", ID, &b_ID);
   fChain->SetBranchAddress("t", t, &b_t);
   fChain->SetBranchAddress("dt", dt, &b_dt);
   fChain->SetBranchAddress("t0", t0, &b_t0);
   fChain->SetBranchAddress("chie", chie, &b_chie);
   fChain->SetBranchAddress("chimu", chimu, &b_chimu);
   fChain->SetBranchAddress("chipi", chipi, &b_chipi);
   fChain->SetBranchAddress("chik", chik, &b_chik);
   fChain->SetBranchAddress("chip", chip, &b_chip);
   fChain->SetBranchAddress("beta", beta, &b_beta);
   fChain->SetBranchAddress("te", te, &b_te);
   fChain->SetBranchAddress("tmu", tmu, &b_tmu);
   fChain->SetBranchAddress("tpi", tpi, &b_tpi);
   fChain->SetBranchAddress("tk", tk, &b_tk);
   fChain->SetBranchAddress("tp", tp, &b_tp);
   Notify();
}

Bool_t RootTof::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void RootTof::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t RootTof::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef RootTof_cxx
