#define analize_cxx
// The class definition in analize.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("analize.C")
// Root > T->Process("analize.C","some options")
// Root > T->Process("analize.C+")
//

#include <iostream>
#include <iomanip>

#include "analize.h"
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>


#include "CrystalBall.h"
#include <mctopo/libMcTopo.h>
//#define mctopo_cxx
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

//mctopo mct;
void analize::Begin(TTree * )
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   SelectionResult::Init();
   output_file_name = option;

}

void analize::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t analize::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either analize::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.
   analize::GetEntry(entry);
   //std::cout << "entry = " << entry << std::endl;
   if(event_treeKK==0) 
   {
     std::cout << "Init event_treeKK" << std::endl;
     event_treeKK=fChain->CloneTree(0);
     event_treeKK->SetName("eventKK");
     event_treeKK->SetTitle("KK events");
   }
   if(mctopo_treeKK==0)
   {
     std::cout << "Init mctopo_treeKK" << std::endl;
     mctopo_treeKK=fChain->GetFriend("mctopo")->CloneTree(0);
     mctopo_treeKK->SetName("mctopoKK");
     mctopo_treeKK->SetTitle("mcTopo KK events");
   }
   if(event_treeUU==0) 
   {
     std::cout << "Init event_treeUU" << std::endl;
     event_treeUU=fChain->CloneTree(0);
     event_treeUU->SetName("eventUU");
     event_treeUU->SetTitle("UU events");
   }
   if(mctopo_treeUU==0)
   {
     std::cout << "Init mctopo_treeUU" << std::endl;
     mctopo_treeUU=fChain->GetFriend("mctopo")->CloneTree(0);
     mctopo_treeUU->SetName("mctopoUU");
     mctopo_treeUU->SetTitle("mcTopo UU events");
   }
   N0++; //count total number of events proceed
   nptrack =  test_hash2(&mctp);
   if(MIN_RECOIL_MASS <= Mrec && Mrec <= MAX_RECOIL_MASS)
     if(pid_chi2 <= PID_CHI2)
       if(kin_chi2 <= KIN_CHI2)
       {
         if(KK==1 && uu==0) 
         {
             NKK++;
             hMrecKK->Fill(mshift(Mrec));
             hpid_chi2KK->Fill(pid_chi2);
             hkin_chi2KK->Fill(kin_chi2);
             event_treeKK->Fill();
             mctopo_treeKK->Fill();
         }
         if(uu==1 && KK==0) 
         {
           hMrecUU->Fill(mshift(Mrec));
           hpid_chi2UU->Fill(pid_chi2);
           hkin_chi2UU->Fill(kin_chi2);
           event_treeUU->Fill();
           mctopo_treeUU->Fill();
           Nuu++;
         }
       }
   int p = log(N0)/log(10)-1;
   if(p<1) p=1;
   int P = pow(10,p);
   if(N0 % P == 0) 
   {
     std::cout << setw(15) << N0 << setw(15) << NKK << setw(15) << Nuu;
     std::cout << "  hash = " << setw(10) << std::hex << test_hash2(&mctp) << "    " << std::dec << info(&mctp) << std::endl;
   }

   //if(N0 % 10 == 0) std::cout << setw(15) << N0 << setw(15) << NKK << setw(15) << Nuu;
   return kTRUE;
}

void analize::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

using namespace std;

void analize::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
   cout << "Total events proceed: " << N0 << endl;
   cout << "Number of KK candidates: " << NKK  << endl;
   cout << "Number of uu candidates: " << Nuu  << endl;

   TCanvas * c =new TCanvas("cMrec","pions recoil mass");
   c->Divide(2,1);
   c->cd(1);
   hMrecKK->Draw("E");
   hMrecKK->SetDrawOption("E");
   char xaxis_title[1024];
   sprintf(xaxis_title,"M_{rec}(#pi^{+}#pi^{-}) - %6.1f, MeV", MJPSI_SHIFT*MSCALE);
   hMrecKK->GetXaxis()->SetTitle(xaxis_title);
   vector<double> resKK(3), resUU(3);
   if(NKK>20)
   {
     resKK  =  Fit(hMrecKK);
     cout << "Number of selected KK events: " << resKK[0] << " " << -resKK[1] << " +" << resKK[2] << endl;
   }
   c->cd(2);
   hMrecUU->Draw("E");
   hMrecUU->SetDrawOption("E");
   hMrecUU->GetXaxis()->SetTitle(xaxis_title);
   if(Nuu>20) 
   {
     resUU  =  Fit(hMrecUU);
     cout << "Number of selected #mu#mu events: " << resUU[0] << " " << -resUU[1] << " +" << resUU[2] << endl;
     double eps = resKK[0]/resUU[0];
     cout << "epsKK/epsUU = " << eps  << "  " << -sqrt( pow(resKK[1]/resKK[0],2) +  pow(resUU[1]/resUU[0],2)) << "  " <<  sqrt( pow(resKK[2]/resKK[0],2) +  pow(resUU[2]/resUU[0],2)) << endl;
   }
   else
   {
     cout << "Too small number of mu mu events: " << Nuu << endl;
   }

  TCanvas * cchi2 = new TCanvas("cchi2","pid and kinematic chi2");
  cchi2->Divide(2,2);
  cchi2->cd(1);
  hpid_chi2KK->Draw();
  cchi2->cd(2);
  hpid_chi2UU->Draw();
  cchi2->cd(3);
  hkin_chi2KK->Draw();
  cchi2->cd(4);
  hkin_chi2UU->Draw();
  file = new TFile(output_file_name.c_str(),"RECREATE");
  file->WriteObject(c,"cMrec");
  file->WriteObject(cchi2,"cchi2");
  file->WriteObject(hMrecKK,"hMrecKK");
  file->WriteObject(hMrecUU,"hMrecUU");
  file->WriteObject(hpid_chi2KK,"hpid_chi2KK");
  file->WriteObject(hpid_chi2UU,"hpid_chi2UU");
  file->WriteObject(hkin_chi2KK,"hkin_chi2KK");
  file->WriteObject(hkin_chi2UU,"hkin_chi2UU");
  file->WriteObject(event_treeKK,"eventKK");
  file->WriteObject(mctopo_treeKK,"mctopoKK");
  file->WriteObject(event_treeUU,"eventUU");
  file->WriteObject(mctopo_treeUU,"mctopoUU");
  //event_tree->Print();
}
