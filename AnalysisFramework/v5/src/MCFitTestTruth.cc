#include "DmtpcSkimDataset.hh"
#include "DmtpcMCDataset.hh"
#include "MaxCamSRIM.hh"
#include "DmtpcMath.hh"
#include "MaxCamImageTools.hh"
#include "DmtpcProjection.hh"
#include "TFile.h" 
#include <iostream>
#include <set> 
#include <vector> 

double minval = 2; 

int main (int nargs, char ** args) 
{

  if (nargs < 3) 
  {
    std::cout << "usage: MCFitTestTruth mcfile outfile [minval=2]" << std::endl;; 
    return 1; 
  }


  DmtpcDataset d; 
  d.openRootFile(args[1]); 

  DmtpcMCDataset mcd; 
  mcd.loadFile(args[1],true); 

  TFile out(args[2],"RECREATE"); 

  if (nargs > 3) minval = atof(args[3]); 

//  std::cout <<"about to load srim " << std::endl; 
  MaxCamSRIM sr("SRIM_F_in_CF4_100Torr"); 
//  std::cout <<"about to set pressure to  " << mcd.getPressure()  << std::endl; 
  sr.setPressure(mcd.getPressure()); 

  std::cout <<" Starting! " << std::endl; 
  out.cd(); 
  for (int tree_i = 0; tree_i < 2; tree_i++)
  {

    TTree * tree = new TTree(tree_i == 0 ? "fit_reco_phi": "fit_truth_phi", tree_i == 0 ? "fit_reco_phi" : "fit_truth_phi"); 

    double x, y; 
    double mcE, recoE, mcRange, recoRange, mcPhi, phi2, mcTheta, mcZ; 
    double mcProjectedRange; 
    int transverseNdof;  
    int nfound, cam; 
    int proj_length; 
    int success; 
    double mcInteg; 
    float pressure = mcd.getPressure(); 
    DmtpcProjection::SRIMLineFitParams p; 


    tree->Branch("mcE",&mcE); 
    tree->Branch("mcInteg",&mcInteg); 
    tree->Branch("mcRange",&mcRange); 
    tree->Branch("mcProjectedRange",&mcProjectedRange); 
    tree->Branch("mcPhi",&mcPhi); 
    tree->Branch("mcTheta",&mcTheta); 
    tree->Branch("mcZ",&mcZ); 
    tree->Branch("nfound",&nfound); 
    tree->Branch("phi2",&phi2); 
    tree->Branch("pressure",&pressure); 
    tree->Branch("x", &x); 
    tree->Branch("y", &y); 


    tree->Branch("recoRange",&recoRange); 

    tree->Branch("fitRange",&p.range); 
    tree->Branch("fitRangeError",&p.rangeError); 
    tree->Branch("fitOffset",&p.x0); 
    tree->Branch("fitOffsetError",&p.x0Error); 
    tree->Branch("fitOffsetErrorLow",&p.x0ErrorLow); 
    tree->Branch("fitOffsetErrorHigh",&p.x0ErrorHigh); 
    tree->Branch("fitY0",&p.y0); 
    tree->Branch("fitY0Error",&p.y0Error); 
    tree->Branch("fitY0ErrorLow",&p.y0ErrorLow); 
    tree->Branch("fitY0ErrorHigh",&p.y0ErrorHigh); 
    tree->Branch("fitY1",&p.y1); 
    tree->Branch("fitY1Error",&p.y1Error); 
    tree->Branch("fitY1ErrorLow",&p.y1ErrorLow); 
    tree->Branch("fitY1ErrorHigh",&p.y1ErrorHigh); 

    tree->Branch("fitChisq",&p.chisq); 
    tree->Branch("fitNdof",&p.ndof); 
    tree->Branch("fitTransverseChisq",&p.transverseChisq); 
    tree->Branch("fitTransverseNdof",&p.transverseNdof); 
    tree->Branch("fitSinThetaSlope",&p.sinThetaSlope); 
    tree->Branch("fitSinThetaMin",&p.sinThetaMin); 
    tree->Branch("fitPhi",&p.phi); 
    tree->Branch("fitProbability",&p.prob); 
    tree->Branch("fitTotalLikely",&p.likelyIntegralAll); 
    tree->Branch("fitSameLikely",&p.likelyIntegralSame); 
    tree->Branch("sigma", &p.sigma); 
    tree->Branch("sigmaError", &p.sigmaError); 
    tree->Branch("fitE", &p.E); 
    tree->Branch("fitEError", &p.EError); 
    tree->Branch("success",&p.success); 
    tree->Branch("gain",&p.gain); 
    tree->Branch("lengthcal",&p.lengthcal); 
    tree->Branch("recoE",&recoE); 
   
    for (int i = 0; i <d.tree()->GetEntries(); i++)
    {
      if (i % 10 == 0) std::cout << i << std::endl; 

      mcd.getEvent(i); 
      d.getEvent(i); 
      mcE= mcd.getParticleE() ; 
      mcInteg = mcd.getE(); 
      mcRange= mcd.getRange(); 
      mcPhi= mcd.getPhi(); 
      mcZ = mcd.getZ(); 
      mcTheta = mcd.getTheta(); 
      mcProjectedRange= mcd.getRange() * sin(mcd.getTheta());  
      int cam = 0; 
      
      if (mcInteg <= 0) continue; 

      nfound = 1; 
      TH2 * img = (TH2*) d.event()->ccdData(0)->Clone("img"); 
      img->Add(d.getBiasFrame(1),-1); 

      MaxCamClusterImage * clust = new  MaxCamClusterImage(img, new TTimeStamp); 

      std::vector<int> pix = mcd.getCluster(0)->getCluster(0); 
      std::set<int> new_pix; 
      for (unsigned pi = 0; pi < pix.size(); pi++)
      {
        if (mcd.getCluster(0)->getImage()->GetBinContent(pix[pi]) >= minval)
        {
          new_pix.insert(pix[pi]);  
        }
      }

      MaxCamImageTools::erode(img, &new_pix); 
      MaxCamImageTools::dilate(img, &new_pix); 
      clust->addCluster(vector<int> (new_pix.begin(), new_pix.end())); 
      clust->getXY(0,x,y);

      //   phi1 = d.event()->cluster(0)->getPhi(0); 
      phi2 = clust->getPhi2(0); 
      recoE = clust->getIntegral(0);      
      recoRange = clust->getLength2(0,phi2);      

  //         ((MaxCamClusterImage*)d.event()->cluster(0))->morphologicalOperation(0,0,6); 

      DmtpcProjection prof(clust,0,tree_i == 0 ? phi2 : mcPhi,x,y,0); 
      proj_length = prof.getLongitudinalProfile()->GetNbinsX();  
      DmtpcProjection::SRIMLineFitParams * params= prof.doSRIMLineFit(&sr, mcd.getLengthcal(), mcd.getGain(),"I"); 
      p = *params; 
      tree->Fill(); 
      delete clust; 
      delete params; 
    }

    tree->Write(); 

  }
  
  out.Close(); 

}







