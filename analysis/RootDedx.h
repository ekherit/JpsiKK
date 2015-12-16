//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec 16 09:37:40 2015 by ROOT version 6.04/12
// from TTree dedx/Dedx info for signal
// found on file: sample.root
//////////////////////////////////////////////////////////

#ifndef RootDedx_h
#define RootDedx_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class RootDedx {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           ntrack;
   Double_t        chie[4];   //[ntrack]
   Double_t        chimu[4];   //[ntrack]
   Double_t        chipi[4];   //[ntrack]
   Double_t        chik[4];   //[ntrack]
   Double_t        chip[4];   //[ntrack]
   Double_t        probPH[4];   //[ntrack]
   Double_t        normPH[4];   //[ntrack]

   // List of branches
   TBranch        *b_ntrack;   //!
   TBranch        *b_chie;   //!
   TBranch        *b_chimu;   //!
   TBranch        *b_chipi;   //!
   TBranch        *b_chik;   //!
   TBranch        *b_chip;   //!
   TBranch        *b_probPH;   //!
   TBranch        *b_normPH;   //!

   RootDedx(TTree *tree=0);
   virtual ~RootDedx();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef RootDedx_cxx
RootDedx::RootDedx(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("sample.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("sample.root");
      }
      f->GetObject("dedx",tree);

   }
   Init(tree);
}

RootDedx::~RootDedx()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t RootDedx::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t RootDedx::LoadTree(Long64_t entry)
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

void RootDedx::Init(TTree *tree)
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
   fChain->SetBranchAddress("chie", chie, &b_chie);
   fChain->SetBranchAddress("chimu", chimu, &b_chimu);
   fChain->SetBranchAddress("chipi", chipi, &b_chipi);
   fChain->SetBranchAddress("chik", chik, &b_chik);
   fChain->SetBranchAddress("chip", chip, &b_chip);
   fChain->SetBranchAddress("probPH", probPH, &b_probPH);
   fChain->SetBranchAddress("normPH", normPH, &b_normPH);
   Notify();
}

Bool_t RootDedx::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void RootDedx::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t RootDedx::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef RootDedx_cxx
