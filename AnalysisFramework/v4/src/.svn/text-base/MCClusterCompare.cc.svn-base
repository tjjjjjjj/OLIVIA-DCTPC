#include "../../../MaxCam/DmtpcSkimDataset.hh"
#include "../../../MaxCam/DmtpcMCDataset.hh"
#include "../../../MaxCam/DmtpcSkimEvent.hh"
#include "../../../MaxCam/MaxCamImageTools.hh"
#include "Math/ProbFuncMathCore.h"
#include "TH2F.h"
#include "TH1.h"
#include <iostream>

const double minE = 0; 
const double maxE = 5000; 
const int Ebins = 50; 
const int cam = 0; 
const double maxDist = 32; 

const int npixbins = 200; 
const int npixmin = 0; 
const int npixmax = 200; 

const int npixratiobins = 20; 
const double npixratiomin = 0; 
const double npixratiomax = 3; 

const int seedvalbins = 500; 
const double seedvalmin = 0; 
const double seedvalmax = 30; 

const int log_prob_bins = 100; 
const double log_prob_min = -500; 
const double log_prob_max = -10; 

bool do_blur = false; 

double blur_amount = 2.7; 
int blur_radius = 9; 

int main(int nargs, char ** args)
{
  if (nargs < 4) 
  {
    std::cout << "Usage: MCCompareClusters skimfile mcfile output.root" << std::endl; 
    return 0; 
  }

  if (nargs >= 5)
  {
    do_blur = true; 
    blur_amount = atof(args[4]); 
    if (nargs >= 6)
    {
      blur_radius = atoi(args[5]); 
    }
  }


  DmtpcSkimDataset sk; 
  DmtpcMCDataset mc; 

  sk.openRootFile(args[1]); 
  mc.loadFile(args[2],true); 

  TFile outf(args[3],"RECREATE"); 

  TH2F * ecompare = new TH2F("ECompare","ECompare",Ebins,minE,maxE, Ebins, minE, maxE); 
  TH2F * distcompare = new TH2F("DistCompare","DistCompare", Ebins,minE,maxE, 10,0,maxDist); 

  TH2F * npixcmp = new TH2F("NPixCompare","NPixCompare",npixbins,npixmin, npixmax,npixbins,npixmin,npixmax); 
  TH2F * npixve = new TH2F("NPixVsE","NPixVsE",Ebins, minE, maxE, npixratiobins, npixratiomin, npixratiomax); 
  TH1F * npixspurious = new TH1F("NPixSpurious","NPixSpurious",npixbins,npixmin,npixmax); 
  TH1F * npixgood = new TH1F("NPixGood","NPixGood",npixbins,npixmin,npixmax); 

  TH1F * nfound = new TH1F("NFound","Nfound", Ebins,minE,maxE);
  TH1F * ntotal = new TH1F("NTotal","NTotal", Ebins,minE,maxE);
  TH1F * nspurious = new TH1F("NSpurious","NSpurious", Ebins,minE,maxE); 

  TH2F * seedvalve = do_blur ? new TH2F("SeedValVsE","SeedValVsE",Ebins,minE,maxE,seedvalbins,seedvalmin,seedvalmax) : 0; 
  TH1F * seedvalgood = do_blur ? new TH1F("SeedValGood","SeedValGood",seedvalbins,seedvalmin,seedvalmax) : 0 ; 
  TH1F * seedvalbad = do_blur ? new TH1F("SeedValSpurious","SeedValSpurious",seedvalbins,seedvalmin,seedvalmax) : 0 ; 

  TH2F * logprobve = do_blur ? new TH2F("LogProbVsE","LogProbVsE", Ebins, minE, maxE, log_prob_bins, log_prob_min, log_prob_max) : 0 ;  
  TH1F * logprobgood = do_blur ? new TH1F("LogProbGood","LogProbGood",log_prob_bins, log_prob_min, log_prob_max) : 0; 
  TH1F * logprobbad = do_blur ? new TH1F("LogProbSpurious","LogProbSpurious",log_prob_bins, log_prob_min, log_prob_max): 0 ; 

  TTree * tree = new TTree("mceval","mceval"); 
   

  int mc_npix, npix,ev; 
  double seedval, E, mcE, dist, wrong_p,x,y,range,mcRange,phi,mcPhi;
  tree->Branch("mc_npix",&mc_npix); 
  tree->Branch("npix",&npix); 
  tree->Branch("E",&E); 
  tree->Branch("mcE",&mcE); 
  tree->Branch("dist",&dist); 
  tree->Branch("x",&x); 
  tree->Branch("y",&y); 
  tree->Branch("ev",&ev); 
  tree->Branch("range",&range);
  tree->Branch("mcRange",&mcRange);
  tree->Branch("phi",&phi);
  tree->Branch("mcPhi",&mcPhi);

  if (do_blur)
  {
    tree->Branch("seedval",&seedval); 
    tree->Branch("wrong_p",&wrong_p); 
  }


  for (int i = 0 ; i < sk.nevents(); i++)
  {
    sk.getEvent(i); 
    mc.getEvent(i);

    
    ev = i; 
    mcE = mc.getE(cam); 
    //Increment ntotal
    if (mcE)
    {
      ntotal->Fill(mcE);
    }


    mc_npix = (int) mc.getCluster(cam)->getCluster(0).size(); 
    mcRange = mc.getRange(); 
    mcPhi = mc.getPhi(); 

    TH2F * blurred = do_blur ? MaxCamImageTools::gaussianBlur(sk.event()->cluster(cam)->getImage(), blur_radius, blur_amount ) : 0; 
    double blurred_rms, blurred_mean; 

    if (do_blur) MaxCamImageTools::meanRMSNoOutliers(blurred,blurred_mean, blurred_rms); 

    for (int t = 0; t < sk.event()->ntracks(cam); t++)
    {
      //Calculate distance 
      double truex, truey; 
      mc.getCluster(cam)->getXY(0,truex,truey); 
      double recox = sk.event()->x(cam,t); 
      double recoy = sk.event()->y(cam,t); 

      x = recox; 
      y = recoy; 

      dist = TMath::Sqrt( TMath::Power(truex - recox,2) + TMath::Power(truey - recoy,2)); 

      //Calculate values
      npix = sk.event()->npixel(cam,t);  
      range = sk.event()->range(cam,t); 
      phi =  sk.event()->phi(cam,t); 
      E = sk.event()->E(cam,t); 
     
      seedval = 0; 
      double chisq = 0; 
      if (do_blur)
      {
        std::vector<int> clust = sk.event()->cluster(cam)->getCluster(t); 
        for (std::vector<int>::iterator it = clust.begin(); it != clust.end(); it++)
        {
          double this_val = blurred->GetArray()[*it]; 
          if (this_val > seedval) seedval = this_val; 
          chisq -= 2 * TMath::Log(TMath::Erfc(this_val / blurred_rms)); 
        }
      }

      wrong_p = TMath::Log10(ROOT::Math::chisquared_cdf_c(chisq, npix -1));  

      //std::cout << wrong_p << std::endl; 
     
      //If centroid distance is too large, mark as spurious
      if (dist > maxDist)
      {
        nspurious->Fill(E); 
        npixspurious->Fill(npix);
        if (do_blur)
        {
          seedvalbad->Fill(seedval); 
          logprobbad->Fill(wrong_p);
        }
      }

      else
      {
        nfound->Fill(mcE); 
        ecompare->Fill(mcE,E); 
        distcompare->Fill(mcE, dist); 
        npixcmp->Fill(mc_npix,npix); 
        npixve->Fill(mcE,((double)npix)/mc_npix); 
        npixgood->Fill(npix); 

        if (do_blur)
        {
          logprobve->Fill(mcE,wrong_p);
          seedvalve->Fill(mcE,seedval); 
          seedvalgood->Fill(seedval); 
          logprobgood->Fill(wrong_p); 
        }
      }
      tree->Fill(); 
    }

    if (do_blur) blurred->Delete();
  }

  TH1F * efficiency = (TH1F*) nfound->Clone("Efficiency"); 
  efficiency->SetName("Efficiency"); 
  efficiency->Divide(ntotal); 

  efficiency->Write(); 
  nspurious->Write(); 
  distcompare->Write(); 
  ecompare->Write(); 
  npixcmp->Write();  
  npixve->Write(); 
  npixspurious->Write(); 
  npixgood->Write(); 

  if (do_blur)
  {
    seedvalve->Write(); 
    seedvalgood->Write(); 
    seedvalbad->Write(); 
    logprobve->Write();
    logprobgood->Write();
    logprobbad->Write();
  }

  tree->Write(); 
  outf.Close(); 

  return 0; 
}
