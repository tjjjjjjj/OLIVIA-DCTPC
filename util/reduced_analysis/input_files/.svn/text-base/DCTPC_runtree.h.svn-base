//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul 14 09:00:30 2014 by ROOT version 5.34/10
// from TTree dctpc_runinfo/Run info
// found on file: hadd.root
//////////////////////////////////////////////////////////

#ifndef DCTPC_runtree_h
#define DCTPC_runtree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class DCTPC_runtree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           RunNum;
   Int_t           SetNum;
   Int_t           SequenceNum;
   Int_t           Exposure_sec;
   Int_t           Totaltracks;
   Int_t           Totaltrigs;
   Int_t           Time_startofrun_sec;
   Double_t        Pressure_torr;
   Double_t        Voltage_amp_volts;
   Double_t        Voltage_drift_volts;

   // List of branches
   TBranch        *b_runnum2;   //!
   TBranch        *b_setnum;   //!
   TBranch        *b_seqnum;   //!
   TBranch        *b_exposure;   //!
   TBranch        *b_totaltrack;   //!
   TBranch        *b_totaltrig;   //!
   TBranch        *b_time_start;   //!
   TBranch        *b_pressure;   //!
   TBranch        *b_voltage_amp;   //!
   TBranch        *b_voltage_drift;   //!

   DCTPC_runtree(TTree *tree=0);
   virtual ~DCTPC_runtree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef DCTPC_runtree_cxx
DCTPC_runtree::DCTPC_runtree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("$BigDCTPC_physics_input");
      if (!f || !f->IsOpen()) {
         f = new TFile("$BigDCTPC_physics_input");
      }
      f->GetObject("dctpc_runinfo",tree);

   }
   Init(tree);
}

DCTPC_runtree::~DCTPC_runtree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t DCTPC_runtree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t DCTPC_runtree::LoadTree(Long64_t entry)
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

void DCTPC_runtree::Init(TTree *tree)
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

   fChain->SetBranchAddress("RunNum", &RunNum, &b_runnum2);
   fChain->SetBranchAddress("SetNum", &SetNum, &b_setnum);
   fChain->SetBranchAddress("SequenceNum", &SequenceNum, &b_seqnum);
   fChain->SetBranchAddress("Exposure_sec", &Exposure_sec, &b_exposure);
   fChain->SetBranchAddress("Totaltracks", &Totaltracks, &b_totaltrack);
   fChain->SetBranchAddress("Totaltrigs", &Totaltrigs, &b_totaltrig);
   fChain->SetBranchAddress("Time_startofrun_sec", &Time_startofrun_sec, &b_time_start);
   fChain->SetBranchAddress("Pressure_torr", &Pressure_torr, &b_pressure);
   fChain->SetBranchAddress("Voltage_amp_volts", &Voltage_amp_volts, &b_voltage_amp);
   fChain->SetBranchAddress("Voltage_drift_volts", &Voltage_drift_volts, &b_voltage_drift);
   Notify();
}

Bool_t DCTPC_runtree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void DCTPC_runtree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t DCTPC_runtree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef DCTPC_runtree_cxx
