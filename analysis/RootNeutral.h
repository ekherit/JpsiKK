//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec 16 09:37:40 2015 by ROOT version 6.04/12
// from TTree neutral/Good neutral tracks
// found on file: sample.root
//////////////////////////////////////////////////////////

#ifndef RootNeutral_h
#define RootNeutral_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class RootNeutral {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           ntrack;
   Double_t        E[15];   //[ntrack]
   Double_t        theta[15];   //[ntrack]
   Double_t        phi[15];   //[ntrack]
   Double_t        time[15];   //[ntrack]
   Double_t        dangle[15];   //[ntrack]

   // List of branches
   TBranch        *b_ntrack;   //!
   TBranch        *b_E;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_time;   //!
   TBranch        *b_dangle;   //!

   RootNeutral(TTree *tree=0);
   virtual ~RootNeutral();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef RootNeutral_cxx
RootNeutral::RootNeutral(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("sample.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("sample.root");
      }
      f->GetObject("neutral",tree);

   }
   Init(tree);
}

RootNeutral::~RootNeutral()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t RootNeutral::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t RootNeutral::LoadTree(Long64_t entry)
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

void RootNeutral::Init(TTree *tree)
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
   fChain->SetBranchAddress("E", E, &b_E);
   fChain->SetBranchAddress("theta", theta, &b_theta);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("time", time, &b_time);
   fChain->SetBranchAddress("dangle", dangle, &b_dangle);
   Notify();
}

Bool_t RootNeutral::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void RootNeutral::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t RootNeutral::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef RootNeutral_cxx
