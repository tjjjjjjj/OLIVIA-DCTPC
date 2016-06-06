//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Aug 16 09:33:16 2014 by ROOT version 5.34/10
// from TTree dctpc_eventinfo/Event info
// found on file: hadd_00955_14141.root
//////////////////////////////////////////////////////////

#ifndef DCTPCTree_h
#define DCTPCTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class DCTPCTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           RunNum;
   Int_t           SetNum;
   Int_t           SequenceNum;
   Int_t           ExposureInRun_sec;
   Int_t           EventNum;
   Double_t        Image_mean_ccdadu;
   Double_t        Image_rms_ccdadu;
   Int_t           Edge;
   Int_t           BurnIn;
   Int_t           Pixels_killed;
   Int_t           LastSpark;
   Int_t           NextSpark;
   Int_t           Ntrack;
   Int_t           Ntrig;
   Double_t        Etrack_kev;
   Double_t        Etrig_kev;
   Double_t        Emesh_kev;
   Double_t        Track_mean_ccdadu;
   Double_t        Track_rms_ccdadu;
   Double_t        Track_x_pix;
   Double_t        Track_y_pix;
   Double_t        Track_x_start_pix;
   Double_t        Track_y_start_pix;
   Double_t        Track_x_end_pix;
   Double_t        Track_y_end_pix;
   Double_t        Track_range_pix;
   Double_t        Track_fitwidth_pix;
   Double_t        Track_width_pix;
   Double_t        Track_maxpixel_ccdadu;
   Int_t           Track_neighbors;
   Int_t           Track_pixels;
   Double_t        Track_phi_deg;
   Double_t        Track_skewness;
   Double_t        Anode_rms_v;
   Double_t        Mesh_rms_v;
   Double_t        Mesh_peak_v;
   Double_t        Mesh_base_v;
   Double_t        Anode_base_v;
   Double_t        Mesh_max_v;
   Double_t        Anode_max_v;
   Double_t        Veto_peak_v;
   Double_t        Mesh_starttime_samp;
   Double_t        Anode_starttime_samp;
   Double_t        Mesh_peaktime_samp;
   Double_t        Mesh_totaltime_samp;
   Double_t        Anode_R0time_samp;
   Double_t        Mesh_R0time_samp;
   Double_t        Mesh_R10time_samp;
   Double_t        Mesh_R25time_samp;
   Double_t        Mesh_R50time_samp;
   Double_t        Mesh_R75time_samp;
   Double_t        Mesh_R90time_samp;
   Double_t        Mesh_F0time_samp;
   Double_t        Mesh_F10time_samp;
   Double_t        Mesh_F25time_samp;
   Double_t        Mesh_F50time_samp;
   Double_t        Mesh_F75time_samp;
   Double_t        Mesh_F90time_samp;
   Int_t           Timenow_sec;
   Int_t           Triggerindex;
   Double_t        Triggertimestamp_samp[10];   //[ntrig]

   // List of branches
   TBranch        *b_runnum;   //!
   TBranch        *b_setnum;   //!
   TBranch        *b_seqnum;   //!
   TBranch        *b_expose;   //!
   TBranch        *b_evnum;   //!
   TBranch        *b_image_mean;   //!
   TBranch        *b_image_rms;   //!
   TBranch        *b_edge;   //!
   TBranch        *b_burnin;   //!
   TBranch        *b_pixels_killed;   //!
   TBranch        *b_last_spark;   //!
   TBranch        *b_next_spark;   //!
   TBranch        *b_ntrack;   //!
   TBranch        *b_ntrig;   //!
   TBranch        *b_Etrack;   //!
   TBranch        *b_Etrig;   //!
   TBranch        *b_Emesh;   //!
   TBranch        *b_cluster_mean;   //!
   TBranch        *b_cluster_rms;   //!
   TBranch        *b_TrackX;   //!
   TBranch        *b_TrackY;   //!
   TBranch        *b_TrackXStart;   //!
   TBranch        *b_TrackYStart;   //!
   TBranch        *b_TrackXEnd;   //!
   TBranch        *b_TrackYEnd;   //!
   TBranch        *b_range_ccd;   //!
   TBranch        *b_Trackrms;   //!
   TBranch        *b_Trackwidth;   //!
   TBranch        *b_Trackmaxpixel;   //!
   TBranch        *b_neighbors;   //!
   TBranch        *b_Tracknpixel;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_Trackmskewness;   //!
   TBranch        *b_anodeRMS;   //!
   TBranch        *b_meshRMS;   //!
   TBranch        *b_mesh_peak;   //!
   TBranch        *b_mesh_base;   //!
   TBranch        *b_anode_base;   //!
   TBranch        *b_mesh_max;   //!
   TBranch        *b_anode_max;   //!
   TBranch        *b_veto_peak;   //!
   TBranch        *b_mesh_start;   //!
   TBranch        *b_anode_start;   //!
   TBranch        *b_mesh_peaktime;   //!
   TBranch        *b_mesh_width;   //!
   TBranch        *b_anode_R0;   //!
   TBranch        *b_mesh_R0;   //!
   TBranch        *b_mesh_R10;   //!
   TBranch        *b_mesh_R25;   //!
   TBranch        *b_mesh_R50;   //!
   TBranch        *b_mesh_R75;   //!
   TBranch        *b_mesh_R90;   //!
   TBranch        *b_mesh_F0;   //!
   TBranch        *b_mesh_F10;   //!
   TBranch        *b_mesh_F25;   //!
   TBranch        *b_mesh_F50;   //!
   TBranch        *b_mesh_F75;   //!
   TBranch        *b_mesh_F90;   //!
   TBranch        *b_timenow;   //!
   TBranch        *b_triggerindex;   //!
   TBranch        *b_Triggertimestamp_samp;   //!

   DCTPCTree(TTree *tree=0);
   virtual ~DCTPCTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef DCTPCTree_cxx
DCTPCTree::DCTPCTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("$BigDCTPC_physics_input");
      if (!f || !f->IsOpen()) {
         f = new TFile("$BigDCTPC_physics_input");
      }
      f->GetObject("dctpc_eventinfo",tree);

   }
   Init(tree);
}

DCTPCTree::~DCTPCTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t DCTPCTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t DCTPCTree::LoadTree(Long64_t entry)
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

void DCTPCTree::Init(TTree *tree)
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

   fChain->SetBranchAddress("RunNum", &RunNum, &b_runnum);
   fChain->SetBranchAddress("SetNum", &SetNum, &b_setnum);
   fChain->SetBranchAddress("SequenceNum", &SequenceNum, &b_seqnum);
   fChain->SetBranchAddress("ExposureInRun_sec", &ExposureInRun_sec, &b_expose);
   fChain->SetBranchAddress("EventNum", &EventNum, &b_evnum);
   fChain->SetBranchAddress("Image_mean_ccdadu", &Image_mean_ccdadu, &b_image_mean);
   fChain->SetBranchAddress("Image_rms_ccdadu", &Image_rms_ccdadu, &b_image_rms);
   fChain->SetBranchAddress("Edge", &Edge, &b_edge);
   fChain->SetBranchAddress("BurnIn", &BurnIn, &b_burnin);
   fChain->SetBranchAddress("Pixels_killed", &Pixels_killed, &b_pixels_killed);
   fChain->SetBranchAddress("LastSpark", &LastSpark, &b_last_spark);
   fChain->SetBranchAddress("NextSpark", &NextSpark, &b_next_spark);
   fChain->SetBranchAddress("Ntrack", &Ntrack, &b_ntrack);
   fChain->SetBranchAddress("Ntrig", &Ntrig, &b_ntrig);
   fChain->SetBranchAddress("Etrack_kev", &Etrack_kev, &b_Etrack);
   fChain->SetBranchAddress("Etrig_kev", &Etrig_kev, &b_Etrig);
   fChain->SetBranchAddress("Emesh_kev", &Emesh_kev, &b_Emesh);
   fChain->SetBranchAddress("Track_mean_ccdadu", &Track_mean_ccdadu, &b_cluster_mean);
   fChain->SetBranchAddress("Track_rms_ccdadu", &Track_rms_ccdadu, &b_cluster_rms);
   fChain->SetBranchAddress("Track_x_pix", &Track_x_pix, &b_TrackX);
   fChain->SetBranchAddress("Track_y_pix", &Track_y_pix, &b_TrackY);
   fChain->SetBranchAddress("Track_x_start_pix", &Track_x_start_pix, &b_TrackXStart);
   fChain->SetBranchAddress("Track_y_start_pix", &Track_y_start_pix, &b_TrackYStart);
   fChain->SetBranchAddress("Track_x_end_pix", &Track_x_end_pix, &b_TrackXEnd);
   fChain->SetBranchAddress("Track_y_end_pix", &Track_y_end_pix, &b_TrackYEnd);
   fChain->SetBranchAddress("Track_range_pix", &Track_range_pix, &b_range_ccd);
   fChain->SetBranchAddress("Track_fitwidth_pix", &Track_fitwidth_pix, &b_Trackrms);
   fChain->SetBranchAddress("Track_width_pix", &Track_width_pix, &b_Trackwidth);
   fChain->SetBranchAddress("Track_maxpixel_ccdadu", &Track_maxpixel_ccdadu, &b_Trackmaxpixel);
   fChain->SetBranchAddress("Track_neighbors", &Track_neighbors, &b_neighbors);
   fChain->SetBranchAddress("Track_pixels", &Track_pixels, &b_Tracknpixel);
   fChain->SetBranchAddress("Track_phi_deg", &Track_phi_deg, &b_phi);
   fChain->SetBranchAddress("Track_skewness", &Track_skewness, &b_Trackmskewness);
   fChain->SetBranchAddress("Anode_rms_v", &Anode_rms_v, &b_anodeRMS);
   fChain->SetBranchAddress("Mesh_rms_v", &Mesh_rms_v, &b_meshRMS);
   fChain->SetBranchAddress("Mesh_peak_v", &Mesh_peak_v, &b_mesh_peak);
   fChain->SetBranchAddress("Mesh_base_v", &Mesh_base_v, &b_mesh_base);
   fChain->SetBranchAddress("Anode_base_v", &Anode_base_v, &b_anode_base);
   fChain->SetBranchAddress("Mesh_max_v", &Mesh_max_v, &b_mesh_max);
   fChain->SetBranchAddress("Anode_max_v", &Anode_max_v, &b_anode_max);
   fChain->SetBranchAddress("Veto_peak_v", &Veto_peak_v, &b_veto_peak);
   fChain->SetBranchAddress("Mesh_starttime_samp", &Mesh_starttime_samp, &b_mesh_start);
   fChain->SetBranchAddress("Anode_starttime_samp", &Anode_starttime_samp, &b_anode_start);
   fChain->SetBranchAddress("Mesh_peaktime_samp", &Mesh_peaktime_samp, &b_mesh_peaktime);
   fChain->SetBranchAddress("Mesh_totaltime_samp", &Mesh_totaltime_samp, &b_mesh_width);
   fChain->SetBranchAddress("Anode_R0time_samp", &Anode_R0time_samp, &b_anode_R0);
   fChain->SetBranchAddress("Mesh_R0time_samp", &Mesh_R0time_samp, &b_mesh_R0);
   fChain->SetBranchAddress("Mesh_R10time_samp", &Mesh_R10time_samp, &b_mesh_R10);
   fChain->SetBranchAddress("Mesh_R25time_samp", &Mesh_R25time_samp, &b_mesh_R25);
   fChain->SetBranchAddress("Mesh_R50time_samp", &Mesh_R50time_samp, &b_mesh_R50);
   fChain->SetBranchAddress("Mesh_R75time_samp", &Mesh_R75time_samp, &b_mesh_R75);
   fChain->SetBranchAddress("Mesh_R90time_samp", &Mesh_R90time_samp, &b_mesh_R90);
   fChain->SetBranchAddress("Mesh_F0time_samp", &Mesh_F0time_samp, &b_mesh_F0);
   fChain->SetBranchAddress("Mesh_F10time_samp", &Mesh_F10time_samp, &b_mesh_F10);
   fChain->SetBranchAddress("Mesh_F25time_samp", &Mesh_F25time_samp, &b_mesh_F25);
   fChain->SetBranchAddress("Mesh_F50time_samp", &Mesh_F50time_samp, &b_mesh_F50);
   fChain->SetBranchAddress("Mesh_F75time_samp", &Mesh_F75time_samp, &b_mesh_F75);
   fChain->SetBranchAddress("Mesh_F90time_samp", &Mesh_F90time_samp, &b_mesh_F90);
   fChain->SetBranchAddress("Timenow_sec", &Timenow_sec, &b_timenow);
   fChain->SetBranchAddress("Triggerindex", &Triggerindex, &b_triggerindex);
   fChain->SetBranchAddress("Triggertimestamp_samp", Triggertimestamp_samp, &b_Triggertimestamp_samp);
   Notify();
}

Bool_t DCTPCTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void DCTPCTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t DCTPCTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef DCTPCTree_cxx
