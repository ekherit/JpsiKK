//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec 30 11:51:45 2015 by ROOT version 6.04/12
// from TTree mdc/Main Drift Chamber
// found on file: sample.root
//////////////////////////////////////////////////////////

#ifndef RootMdc_h
#define RootMdc_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class RootMdc {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           ntrack;
   Int_t           trackId[4];   //[ntrack]
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
   Double_t        depth[4];   //[ntrack]
   Double_t        Mrec;
   Int_t           npid;
   Double_t        M23[5];   //[npid]
   Double_t        M12[5];   //[npid]
   Double_t        M03[5];   //[npid]
   Double_t        Mmis[5];   //[npid]

   // List of branches
   TBranch        *b_ntrack;   //!
   TBranch        *b_trackId;   //!
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
   TBranch        *b_depth;   //!
   TBranch        *b_Mrec;   //!
   TBranch        *b_npid;   //!
   TBranch        *b_M23;   //!
   TBranch        *b_M12;   //!
   TBranch        *b_M03;   //!
   TBranch        *b_Mmis;   //!

   RootMdc(TTree *tree=0);
   virtual ~RootMdc();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef RootMdc_cxx
RootMdc::RootMdc(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("sample.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("sample.root");
      }
      f->GetObject("mdc",tree);

   }
   Init(tree);
}

RootMdc::~RootMdc()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t RootMdc::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t RootMdc::LoadTree(Long64_t entry)
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

void RootMdc::Init(TTree *tree)
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
   fChain->SetBranchAddress("trackId", trackId, &b_trackId);
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
   fChain->SetBranchAddress("depth", depth, &b_depth);
   fChain->SetBranchAddress("Mrec", &Mrec, &b_Mrec);
   fChain->SetBranchAddress("npid", &npid, &b_npid);
   fChain->SetBranchAddress("M23", M23, &b_M23);
   fChain->SetBranchAddress("M12", M12, &b_M12);
   fChain->SetBranchAddress("M03", M03, &b_M03);
   fChain->SetBranchAddress("Mmis", Mmis, &b_Mmis);
   Notify();
}

Bool_t RootMdc::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void RootMdc::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t RootMdc::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef RootMdc_cxx
