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

#include "analize.h"
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>


//#include "mctopo.h"
#include "CrystalBall.h"

//mctopo mct;

void analize::Begin(TTree * )
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   SelectionResult::Init();

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
   N0++; //count total number of events proceed

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
         }
         if(uu==1 && KK==0) 
         {
           hMrecUU->Fill(mshift(Mrec));
           hpid_chi2UU->Fill(pid_chi2);
           hkin_chi2UU->Fill(kin_chi2);
           Nuu++;
         }
       }
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

   TCanvas * c =new TCanvas;
   c->Divide(2,1);
   c->cd(1);
   hMrecKK->Draw();
   auto resKK  =  Fit2(hMrecKK);
   c->cd(2);
   hMrecUU->Draw();
   auto resUU  =  Fit2(hMrecUU);


  cout << "Number of selected KK events: " << resKK[0] << " " << -resKK[1] << " +" << resKK[2] << endl;
  cout << "Number of selected UU events: " << resUU[0] << " " << -resUU[1] << " +" << resUU[2] << endl;
  double eps = resKK[0]/resUU[0];
  cout << "epsKK/epsUU = " << eps  << "  " << -sqrt( pow(resKK[1]/resKK[0],2) +  pow(resUU[1]/resUU[0],2)) << "  " <<  sqrt( pow(resKK[2]/resKK[0],2) +  pow(resUU[2]/resUU[0],2)) << endl;
	//CrystalBallFitter2 cbKK(hMrecKK);
	//CrystalBallFitter2 cbUU(hMrecUU);
	//cbKK.Fit();
  TCanvas * cchi2 = new TCanvas;
  cchi2->Divide(2,2);
  cchi2->cd(1);
  hpid_chi2KK->Draw();
  cchi2->cd(2);
  hkin_chi2KK->Draw();
  cchi2->cd(3);
  hpid_chi2UU->Draw();
  cchi2->cd(4);
  hkin_chi2UU->Draw();
}
