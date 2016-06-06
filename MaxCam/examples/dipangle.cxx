// for Justin Easley's study of the 3d track reconstruction from 2D plus energy


#include <iostream>

#include "TString.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"

#include "../MaxCamSRIM.hh"

using namespace std;


void rangeVsEnergy() {

  // This code assumes that you have checked out and compiled MaxCam
  // and that you have access to the SRIMT tables (in MaxCam/tables)

  // get SRIM prediction for Alphas in 75 torr CF4
  // The range of a particle in a gas depends on the gas type, 
  // the gas density (pressure and temperature), the particle charge and mass
  // and the initial energy of the particle

  // An alpha particle (He nucleus) in CF4
  MaxCamSRIM *srimHe = new MaxCamSRIM("SRIM_He_in_CF4_100Torr");
  // at 75 Torr
  srimHe->setPressure(75);

  // ditto, but for a Fluorine nucleus
  MaxCamSRIM *srimF = new MaxCamSRIM("SRIM_F_in_CF4_100Torr");
  srimF->setPressure(75);
  
  TCanvas *c1 = new TCanvas("c1", "c1", 200, 10, 600, 400);

  float xmin=0.; 
  float xmax=100.;
  
  TGraph *rangeEnergyHe = srimHe->getRangeVsEnergy();
  float ymax=rangeEnergyHe->Eval(xmax)*1.1;
  rangeEnergyHe->GetXaxis()->SetLimits(xmin, xmax);
  rangeEnergyHe->SetMaximum(ymax);
  rangeEnergyHe->GetXaxis()->SetTitle("Energy [keV]");
  rangeEnergyHe->GetYaxis()->SetTitle("Range [mm]");
  rangeEnergyHe->Draw("AC");

  TGraph *rangeEnergyF = srimF->getRangeVsEnergy();
  rangeEnergyF->Draw("L");

  // what is the 3D range for a Fluorine nucleus that begins
  // with 100 keV energy?
  cout << "100 keV F recoil in CF_4 gas at 75 torr will travel (mm):  " << rangeEnergyF->Eval(100.) << endl;

  c1->Update();
}

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
}

