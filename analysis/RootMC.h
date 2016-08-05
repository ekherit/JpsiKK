//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec 30 11:51:45 2015 by ROOT version 6.04/12
// from TTree mc/Monte Carlo truth information
// found on file: sample.root
//////////////////////////////////////////////////////////

#ifndef RootMC_h
#define RootMC_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class RootMC {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           psip_decay;
   Int_t           jpsi_decay;
   Int_t           KK;
   Int_t           uu;
   Int_t           oo;
   Int_t           ntrack;
   Double_t        id[4];   //[ntrack]
   Double_t        q[4];   //[ntrack]
   Double_t        E[4];   //[ntrack]
   Double_t        p[4];   //[ntrack]
   Double_t        px[4];   //[ntrack]
   Double_t        py[4];   //[ntrack]
   Double_t        pz[4];   //[ntrack]
   Double_t        pt[4];   //[ntrack]
   Double_t        theta[4];   //[ntrack]
   Double_t        phi[4];   //[ntrack]

   // List of branches
   TBranch        *b_psip_decay;   //!
   TBranch        *b_jpsi_decay;   //!
   TBranch        *b_KK;   //!
   TBranch        *b_uu;   //!
   TBranch        *b_oo;   //!
   TBranch        *b_ntrack;   //!
   TBranch        *b_id;   //!
   TBranch        *b_q;   //!
   TBranch        *b_E;   //!
   TBranch        *b_p;   //!
   TBranch        *b_px;   //!
   TBranch        *b_py;   //!
   TBranch        *b_pz;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_phi;   //!

   RootMC(TTree *tree=0);
   virtual ~RootMC();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef RootMC_cxx
RootMC::RootMC(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("sample.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("sample.root");
      }
      f->GetObject("mc",tree);

   }
   Init(tree);
}

RootMC::~RootMC()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t RootMC::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t RootMC::LoadTree(Long64_t entry)
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

void RootMC::Init(TTree *tree)
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

   fChain->SetBranchAddress("psip_decay", &psip_decay, &b_psip_decay);
   fChain->SetBranchAddress("jpsi_decay", &jpsi_decay, &b_jpsi_decay);
   fChain->SetBranchAddress("KK", &KK, &b_KK);
   fChain->SetBranchAddress("uu", &uu, &b_uu);
   fChain->SetBranchAddress("oo", &oo, &b_oo);
   fChain->SetBranchAddress("ntrack", &ntrack, &b_ntrack);
   fChain->SetBranchAddress("id", id, &b_id);
   fChain->SetBranchAddress("q", q, &b_q);
   fChain->SetBranchAddress("E", E, &b_E);
   fChain->SetBranchAddress("p", p, &b_p);
   fChain->SetBranchAddress("px", px, &b_px);
   fChain->SetBranchAddress("py", py, &b_py);
   fChain->SetBranchAddress("pz", pz, &b_pz);
   fChain->SetBranchAddress("pt", pt, &b_pt);
   fChain->SetBranchAddress("theta", theta, &b_theta);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   Notify();
}

Bool_t RootMC::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void RootMC::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t RootMC::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef RootMC_cxx
