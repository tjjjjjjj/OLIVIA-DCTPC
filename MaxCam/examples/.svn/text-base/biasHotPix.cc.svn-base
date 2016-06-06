#include "../DmtpcDataset.hh" 
#include "../MaxCamImageTools.hh" 
#include "../MaxCamConfig.hh" 
#include "TTree.h"
#include "TFile.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "TObjString.h" 



int main (int nargs, char ** args) 
{
  if (nargs < 3) 
  {
    printf("Usage: biasHotPix rawFiles.txt out.root\n"); 
    return 1; 
  }


  std::vector<std::string> files; 
  std::ifstream files_file(args[1]); 
  std::string line; 
  while (!files_file.eof())
  {
    std::getline(files_file, line); 
    if (line.length() > 0)
      files.push_back(line); 
  }
  files_file.close(); 

  TFile f(args[2],"RECREATE"); 

  TObjString * serial = new TObjString; 
  int nhotpix; 
  double rms,rms_no; 
  double mean,mean_no; 
  double median; 
  int expo; 
  int run; 

  TTree * tree = new TTree("bias","bias"); 
  tree->Branch("serial",&serial); 
  tree->Branch("nhotpix",&nhotpix); 
  tree->Branch("rms",&rms); 
  tree->Branch("mean",&mean); 
  tree->Branch("rms_no",&rms_no); 
  tree->Branch("mean_no",&mean_no); 
  tree->Branch("median",&median); 
  tree->Branch("run",&run); 
  tree->Branch("expo",expo); 

  

  for (unsigned i = 0; i < files.size(); i++)
  {
    printf("analyzing file %s\n",files[i].c_str()); 

    DmtpcDataset d; 
    d.openRootFile(files[i].c_str()); 
    d.getEvent(0); 

    int ncam = d.getNcameras(); 


    for (int c = 0; c < ncam; c++)
    {
      TH2F * bias = d.getBiasFrame(c+1); 
      serial->SetString( d.event()->ccdConfig(c)->serialNumber); 
      expo = d.event()->ccdConfig(c)->exposureTime; 
      run = d.event()->runNumber(); 
      MaxCamImageTools::meanRMSNoOutliers(bias,mean_no,rms_no); 
      mean = bias->Integral() / (bias->GetNbinsX() * bias->GetNbinsY()); 
      rms = MaxCamImageTools::getRMS(bias); 
      median = MaxCamImageTools::median(bias); 
      
      TH2 * neigh = MaxCamImageTools::neighborRatio(bias,true,true,true); 
      nhotpix = 0; 

      for (int x = 1; x <= neigh->GetNbinsX(); x++)
      {
        for (int y = 1; y <= neigh->GetNbinsY(); y++)
        {
          if (neigh->GetBinContent(x,y) > 20)
          {
            nhotpix++;
          }
        }
      }

      neigh->Delete();
        
      f.cd(); 
      tree->Fill();


    }

  }

  tree->Write();



}





