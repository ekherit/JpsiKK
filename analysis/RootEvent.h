//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Dec 17 13:24:32 2015 by ROOT version 6.04/12
// from TTree event/Signal events pi+pi- K+K-, or pi+pi- mu+mu-
// found on file: sample.root
//////////////////////////////////////////////////////////

#ifndef RootEvent_h
#define RootEvent_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class RootEvent {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

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
   Int_t           ntrack;
   Int_t           sign;
   Int_t           channel;
   Int_t           KK;
   Int_t           uu;
   Int_t           K;
   Int_t           u;
   Double_t        Mrec;
   Double_t        M012;
   Double_t        M013;
   Double_t        M023;
   Double_t        M123;
   Double_t        M03;
   Double_t        M12;
   Double_t        M01;
   Double_t        M23;
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
   Double_t        kin_chi2;
   Double_t        pid_chi2;
   Int_t           npid;
   Double_t        kchi[5];   //[npid]
   Double_t        pchi[5];   //[npid]
   Double_t        kM[5];   //[npid]
   Int_t           nkinbg;
   Double_t        kin_chi2_bg[10];   //[nkinbg]
   Double_t        Mpi0;

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
   TBranch        *b_ntrack;   //!
   TBranch        *b_sign;   //!
   TBranch        *b_channel;   //!
   TBranch        *b_KK;   //!
   TBranch        *b_uu;   //!
   TBranch        *b_K;   //!
   TBranch        *b_u;   //!
   TBranch        *b_Mrec;   //!
   TBranch        *b_M012;   //!
   TBranch        *b_M013;   //!
   TBranch        *b_M023;   //!
   TBranch        *b_M123;   //!
   TBranch        *b_M03;   //!
   TBranch        *b_M12;   //!
   TBranch        *b_M01;   //!
   TBranch        *b_M23;   //!
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
   TBranch        *b_kin_chi2;   //!
   TBranch        *b_pid_chi2;   //!
   TBranch        *b_npid;   //!
   TBranch        *b_kchi;   //!
   TBranch        *b_pchi;   //!
   TBranch        *b_kM;   //!
   TBranch        *b_nkinbg;   //!
   TBranch        *b_kin_chi2_bg;   //!
   TBranch        *b_Mpi0;   //!

   RootEvent(TTree *tree=0);
   virtual ~RootEvent();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef RootEvent_cxx
RootEvent::RootEvent(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("sample.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("sample.root");
      }
      f->GetObject("event",tree);

   }
   Init(tree);
}

RootEvent::~RootEvent()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t RootEvent::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t RootEvent::LoadTree(Long64_t entry)
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

void RootEvent::Init(TTree *tree)
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
   fChain->SetBranchAddress("ntrack", &ntrack, &b_ntrack);
   fChain->SetBranchAddress("sign", &sign, &b_sign);
   fChain->SetBranchAddress("channel", &channel, &b_channel);
   fChain->SetBranchAddress("KK", &KK, &b_KK);
   fChain->SetBranchAddress("uu", &uu, &b_uu);
   fChain->SetBranchAddress("K", &K, &b_K);
   fChain->SetBranchAddress("u", &u, &b_u);
   fChain->SetBranchAddress("Mrec", &Mrec, &b_Mrec);
   fChain->SetBranchAddress("M012", &M012, &b_M012);
   fChain->SetBranchAddress("M013", &M013, &b_M013);
   fChain->SetBranchAddress("M023", &M023, &b_M023);
   fChain->SetBranchAddress("M123", &M123, &b_M123);
   fChain->SetBranchAddress("M03", &M03, &b_M03);
   fChain->SetBranchAddress("M12", &M12, &b_M12);
   fChain->SetBranchAddress("M01", &M01, &b_M01);
   fChain->SetBranchAddress("M23", &M23, &b_M23);
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
   fChain->SetBranchAddress("kin_chi2", &kin_chi2, &b_kin_chi2);
   fChain->SetBranchAddress("pid_chi2", &pid_chi2, &b_pid_chi2);
   fChain->SetBranchAddress("npid", &npid, &b_npid);
   fChain->SetBranchAddress("kchi", kchi, &b_kchi);
   fChain->SetBranchAddress("pchi", pchi, &b_pchi);
   fChain->SetBranchAddress("kM", kM, &b_kM);
   fChain->SetBranchAddress("nkinbg", &nkinbg, &b_nkinbg);
   fChain->SetBranchAddress("kin_chi2_bg", kin_chi2_bg, &b_kin_chi2_bg);
   fChain->SetBranchAddress("Mpi0", &Mpi0, &b_Mpi0);
   Notify();
}

Bool_t RootEvent::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void RootEvent::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t RootEvent::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef RootEvent_cxx
