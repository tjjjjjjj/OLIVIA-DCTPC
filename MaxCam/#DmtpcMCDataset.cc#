#include "DmtpcMCDataset.hh"

#include <iostream>
ClassImp(DmtpcMCDataset);

DmtpcMCDataset::DmtpcMCDataset()
{
  f = 0; 
  tree = 0; 
  ncam = 0; 
  phi = 0; 
  sequence = -1;
  integ = 0; 
  tree = 0; 
  clusters = 0; 
  gain = 0; 
  pressure = 0; 
}

DmtpcMCDataset::~DmtpcMCDataset()
{

  //if (tree) tree->Delete(); 
  if (clusters)
  {
    clusters->Delete(); 
  }

  if (integ) delete[] integ; 
  if (f) 
  {
    f->Close(); 
    f->Delete();
  }
}

bool DmtpcMCDataset::loadFile(const char * filename, bool load_true_clusters)
{
  if (f)
  {
    f->Close();
    f->Delete(); 
  }

//  if (tree) tree->Delete(); 

  f = new TFile(filename); 

  if (!f->IsOpen()) return false; 
  gROOT->cd(); 
  tree = (TTree*) f->Get("Simulation"); 
  
  if (!tree) return false; 

  tree->SetBranchAddress("ncamera",&ncam); 
  tree->GetEntry(0); 

  if (integ) delete[] integ; 

  integ = new float[ncam]; 


  tree->SetBranchAddress("Integral",integ); 
  tree->SetBranchAddress("trackPhi",&phi); 
  tree->SetBranchAddress("trackTheta",&theta); 
  tree->SetBranchAddress("length",&length);
  tree->SetBranchAddress("zlength",&zlength); 
  tree->SetBranchAddress("E",&E); 
  tree->SetBranchAddress("Escint",&Escint); 
  tree->SetBranchAddress("Sequence",&sequence); 
  tree->SetBranchAddress("x",&x); 
  tree->SetBranchAddress("y",&y); 
  tree->SetBranchAddress("z",&z); 
  

  if (load_true_clusters)
  {

    if (!tree->FindBranch("trueClusterArray"))
    {
      load_true_clusters = false; 
    }
    else
    {
      tree->SetBranchAddress("trueClusterArray",&clusters); 
    }
    
  }

  tree->GetEntry(0); 

  TTree * cam_tree = (TTree*) f->Get("camera"); 
  cam_tree->SetBranchAddress("pixels",pix); 
  cam_tree->SetBranchAddress("width",width); 
  cam_tree->SetBranchAddress("cameraPosition",campos); 
  cam_tree->SetBranchAddress("gain",&gain); 
  cam_tree->SetBranchAddress("readNoise",&noise); 
  cam_tree->SetBranchAddress("bias",&bias); 
  cam_tree->GetEntry(0); 

  TTree * chamb_tree = (TTree*) f->Get("chamber"); 
  chamb_tree->SetBranchAddress("pressure",&pressure); 
  chamb_tree->GetEntry(0); 
  lengthcal = width[0] / pix[0]; 

  return true;
}
 
 

