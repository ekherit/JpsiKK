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

#include "analize.h"
#include <TH2.h>
#include <TStyle.h>


#include "mctopo.h"

//mctopo mct;

void analize::Begin(TTree * tree)
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
   selector::GetEntry(entry);
   N0++; //count total number of events proceed
   if(3.08 <= Mrec && Mrec <= 3.105)
   {
     if(KK==1 && uu==0) 
     {
       NKK++;
       hMrecKK->Fill(Mrec);
     }
     if(uu==1 && KK==0) 
     {
       hMrecUU->Fill(Mrec);
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
   c->cd(2);
   hMrecUU->Draw();
}
