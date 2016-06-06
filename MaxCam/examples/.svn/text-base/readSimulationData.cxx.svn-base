// for Justin Easley's study of the 3d track reconstruction from 2D plus energy


#include <iostream>

#include "TString.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"

#include "../MaxCamSRIM.hh"

using namespace std;

void readSimulation() {
  // this is sample code to read data from the simulation
  // Assumes that you have run the Simulation (see MaxCam/Simulation/)
  // and that you can access files jbattat's directory on zwicky...
  TString simfile("/net/zwicky/dmtpc/jbattat/projects/DarkMatter/MaxCam/Simulations/v1/output/data/dmtpc_run091209001.root");

  TFile f(simfile);
  f.ls();
  TTree *simtree = (TTree*)f.Get("Simulation");
  simtree->Print();
  simtree->Draw("length:E");

  // other variables of interest are:
  // Escint ([ncamera]/F) (?keV?)
  //        Energy that went into scintillation.
  //        This is the only energy detectable by the ccd) 
  // phi  (phi/F) (radians) angle of recoil in the plane of the image 
  // theta (theta/F) (radians) angle of recoil out of the plane
  // length (length/F) (?mm?) 3D recoil length
  
  // so to construct the 2D recoil length in the plane of the 
  // ccd image, you need to do trig with theta and length

  // figure out how many cameras there were in the simulation
  Int_t ncamera;
  simtree->SetBranchAddress("ncamera", &ncamera);
  simtree->GetEvent(0);
  cout << "ncamera = " << ncamera << endl;

  // extract angle and other variables from simulation
  Float_t phi;
  Float_t Escint[ncamera];
  simtree->SetBranchAddress("phi", &phi);
  simtree->SetBranchAddress("Escint", Escint);
  for (int ii=0; ii<5; ii++) {
    simtree->GetEvent(ii);
    cout << "ii, phi, Escint[0] = " << ii << ", " << phi << ", " << Escint[0] << endl;
  }

}

