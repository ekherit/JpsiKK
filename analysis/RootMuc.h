//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec 30 11:51:45 2015 by ROOT version 6.04/12
// from TTree muc/MUon chamber information
// found on file: sample.root
//////////////////////////////////////////////////////////

#ifndef RootMuc_h
#define RootMuc_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class RootMuc {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           valid;
   Int_t           ntrack;
   Double_t        status[4];   //[ntrack]
   Double_t        type[4];   //[ntrack]
   Double_t        depth[4];   //[ntrack]
   Double_t        chi2[4];   //[ntrack]
   Double_t        ndf[4];   //[ntrack]
   Double_t        distance[4];   //[ntrack]
   Double_t        phi[4];   //[ntrack]
   Double_t        nhit[4];   //[ntrack]
   Double_t        nlayer[4];   //[ntrack]
   Double_t        nhitmax[4];   //[ntrack]
   Double_t        brlast[4];   //[ntrack]
   Double_t        eclast[4];   //[ntrack]

   // List of branches
   TBranch        *b_valid;   //!
   TBranch        *b_ntrack;   //!
   TBranch        *b_status;   //!
   TBranch        *b_type;   //!
   TBranch        *b_depth;   //!
   TBranch        *b_chi2;   //!
   TBranch        *b_ndf;   //!
   TBranch        *b_distance;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_nhit;   //!
   TBranch        *b_nlayer;   //!
   TBranch        *b_nhitmax;   //!
   TBranch        *b_brlast;   //!
   TBranch        *b_eclast;   //!

   RootMuc(TTree *tree=0);
   virtual ~RootMuc();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef RootMuc_cxx
RootMuc::RootMuc(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("sample.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("sample.root");
      }
      f->GetObject("muc",tree);

   }
   Init(tree);
}

RootMuc::~RootMuc()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t RootMuc::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t RootMuc::LoadTree(Long64_t entry)
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

void RootMuc::Init(TTree *tree)
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

   fChain->SetBranchAddress("valid", &valid, &b_valid);
   fChain->SetBranchAddress("ntrack", &ntrack, &b_ntrack);
   fChain->SetBranchAddress("status", status, &b_status);
   fChain->SetBranchAddress("type", type, &b_type);
   fChain->SetBranchAddress("depth", depth, &b_depth);
   fChain->SetBranchAddress("chi2", chi2, &b_chi2);
   fChain->SetBranchAddress("ndf", ndf, &b_ndf);
   fChain->SetBranchAddress("distance", distance, &b_distance);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("nhit", nhit, &b_nhit);
   fChain->SetBranchAddress("nlayer", nlayer, &b_nlayer);
   fChain->SetBranchAddress("nhitmax", nhitmax, &b_nhitmax);
   fChain->SetBranchAddress("brlast", brlast, &b_brlast);
   fChain->SetBranchAddress("eclast", eclast, &b_eclast);
   Notify();
}

Bool_t RootMuc::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void RootMuc::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t RootMuc::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef RootMuc_cxx
