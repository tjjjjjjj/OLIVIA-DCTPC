#include "DmtpcSkimDataset.hh"
#include "DmtpcMCDataset.hh"
#include "TFile.h" 
#include <iostream>


int main (int nargs, char ** args) 
{

  if (nargs < 4) 
  {
    std::cout << "usage: MCCompare mcfile skimfile outfile" << std::endl;; 
    return 1; 
  }


  double gain = 40; 
  double lengthcal = 0.179; 

  DmtpcSkimDataset d; 
  d.openRootFile(args[2]); 

  DmtpcMCDataset mcd; 
  mcd.loadFile(args[1],true); 

  TFile out(args[3],"RECREATE"); 

  TTree * tree = new TTree("cmp","cmp"); 
  double mcE, recoE, mcRange, recoRange, mcPhi, recoPhi, mcX, recoX, mcY, recoY, mcZ, recoWidth, cluster_rms, cluster_mean, mcInteg; 
  double dist, mc2DRange, mcTheta, mcClusterX, mcClusterY; 
  int nfound,npixel, mc_npixel; 
  int ncam, which_cam; 

  tree->Branch("mcE",&mcE); 
  tree->Branch("mcInteg",&mcInteg); 
  tree->Branch("ncam",&ncam); 
  tree->Branch("which_cam",&which_cam); 
  tree->Branch("recoE",&recoE); 
  tree->Branch("mcRange",&mcRange); 
  tree->Branch("mcTheta",&mcTheta); 
  tree->Branch("mc2DRange",&mc2DRange); 
  tree->Branch("recoRange",&recoRange); 
  tree->Branch("mcPhi",&mcPhi); 
  tree->Branch("recoPhi",&recoPhi); 
  tree->Branch("mcX",&mcX); 
  tree->Branch("mcClusterX",&mcClusterX); 
  tree->Branch("mcClusterY",&mcClusterY); 
  tree->Branch("recoX",&recoX); 
  tree->Branch("mcY",&mcY); 
  tree->Branch("recoY",&recoY); 
  tree->Branch("mcZ",&mcZ); 
  tree->Branch("npixel",&npixel); 
  tree->Branch("mc_npixel",&mc_npixel); 
  tree->Branch("recoWidth",&recoWidth); 
  tree->Branch("nfound",&nfound); 
  tree->Branch("cluster_rms",&cluster_rms); 
  tree->Branch("cluster_mean",&cluster_mean); 
  tree->Branch("dist",&dist); 

  for (int i = 0; i < d.nevents(); i++)
  {

    nfound = 0; 
    which_cam = -1; 
    ncam = 0; 
    d.getEvent(i); 
    mcd.getEvent(i); 
    mcE= mcd.getParticleE() ; 
    mcRange= mcd.getRange(); 
    mc2DRange = mcd.getRange() * sin(mcd.getTheta()); 
    mcTheta = mcd.getTheta(); 
    mcPhi= mcd.getPhi(); 
    mcX= mcd.getX(true)  / lengthcal; 
    mcY= mcd.getY(true)  / lengthcal; 
    mcZ = mcd.getZ(); 
    mcInteg = 0; 
    dist = -1; 

    for (int j = 0; j < d.event()->ncamera(); j++)
    {
      if (mcd.getE(j) > mcInteg) 
      {
        mcInteg = mcd.getE(j); 
        which_cam = j; 
        mc_npixel = (int) mcd.getCluster(j)->getCluster(0).size(); 
        mcd.getCluster(j)->getXY(0,mcClusterX,mcClusterY); 
      }

      if (mcd.getE(j) > 0) 
      {
        ncam++; 
      }
    }

    for (int j = 0; j < d.event()->ncamera(); j++)
    {
      nfound += d.event()->ntracks(j); 

      if (d.event()->ntracks(j))
      {
         recoE = d.event()->E(j,0);      
         recoRange = d.event()->range(j,0);      
         recoPhi = d.event()->phi(j,0);      
         recoX = d.event()->x(j,0);      
         recoY = d.event()->y(j,0);      
         npixel = d.event()->npixel(j,0);      
         recoWidth = sqrt(d.event()->transverse_moment(j,2,0));      
         cluster_rms = d.event()->cluster_rms(j,0); 
         cluster_mean = d.event()->cluster_mean(j,0); 
         if (j == which_cam)
         {
           dist = sqrt(pow(mcClusterX-recoX,2) + pow(mcClusterY-recoY,2));  

         }
      }
    }

    tree->Fill(); 
  }


  tree->Write(); 
  out.Close(); 

}



