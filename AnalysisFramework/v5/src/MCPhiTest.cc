#include "DmtpcSkimDataset.hh"
#include "DmtpcMCDataset.hh"
#include "DmtpcMath.hh"
#include "DmtpcProjection.hh"
#include "TFile.h" 
#include <iostream>


int main (int nargs, char ** args) 
{

  if (nargs < 4) 
  {
    std::cout << "usage: MCCompare mcfile skimfile outfile" << std::endl;; 
    return 1; 
  }


  DmtpcSkimDataset d; 
  d.openRootFile(args[2]); 

  DmtpcMCDataset mcd; 
  mcd.loadFile(args[1]); 

  TFile out(args[2],"RECREATE"); 

  TTree * tree = new TTree("phicmp","phicmp"); 

  double mcE, recoE, mcRange, recoRange, mcPhi, phi1,phi2,phi3,phi4,phi5,  mcTheta, mcZ; 
  double phi_vec; 
  double phi2_vec; 
  double phi2_12; 
  double phi2_11; 
  double phi2_13; 
  double phi2_21; 
  double phi2_22; 
  double phi2_23; 
  double phi2_01; 
  double phi2_skew; 
  double phi2_greater; 
  double phi2_slope; 
  double phi2_10; 
  int nfound, cam; 
  int proj_length; 
  double mcInteg; 

  tree->Branch("mcE",&mcE); 
  tree->Branch("mcInteg",&mcInteg); 
  tree->Branch("recoE",&recoE); 
  tree->Branch("mcRange",&mcRange); 
  tree->Branch("recoRange",&recoRange); 
  tree->Branch("mcPhi",&mcPhi); 
  tree->Branch("phi11",&phi1); 
  tree->Branch("phi2",&phi2); 
  tree->Branch("phi2_skew",&phi2_skew); 
  tree->Branch("phi2_greater",&phi2_greater); 
  tree->Branch("phi2_slope",&phi2_slope); 
  tree->Branch("proj_length",&proj_length); 
  tree->Branch("phi2_11",&phi2_11); 
  tree->Branch("phi2_12",&phi2_12); 
  tree->Branch("phi2_13",&phi2_13); 
  tree->Branch("phi2_21",&phi2_21); 
  tree->Branch("phi2_22",&phi2_22); 
  tree->Branch("phi2_23",&phi2_23); 
  tree->Branch("phi2_01",&phi2_01); 
  tree->Branch("phi2_10",&phi2_10); 
  tree->Branch("phi2_vec",&phi2_vec); 
  tree->Branch("phi_vec",&phi_vec); 
  tree->Branch("phi3",&phi3); 
  tree->Branch("phi4",&phi4); 
  tree->Branch("phi5",&phi5); 
  tree->Branch("cam",&cam); 
  tree->Branch("nfound",&nfound); 
  tree->Branch("mcTheta",&mcTheta); 
  tree->Branch("mcZ",&mcZ); 
 
  for (int i = 0; i < d.nevents(); i++)
  {

    nfound = 0; 
    mcd.getEvent(i); 
    mcE= mcd.getParticleE() ; 
    mcInteg = mcd.getE(); 
    mcRange= mcd.getRange(); 
    mcPhi= mcd.getPhi(); 
    mcZ = mcd.getZ(); 
    mcTheta = mcd.getTheta(); 
    MaxCamClusterImage * clust = mcd.getCluster(0); 

    for (int j = 0; j < d.event()->ncamera(); j++)
    {
      nfound += d.event()->ntracks(j); 

      if (d.event()->ntracks(j))
      {
         recoE = d.event()->E(j,0);      
         recoRange = d.event()->range(j,0);      

      //   phi1 = d.event()->cluster(0)->getPhi(0); 
         phi2 = d.event()->cluster(0)->getPhi2(0); 
         phi2_11 = d.event()->cluster(0)->getPhi2(0,1,1); 
         phi2_12 = d.event()->cluster(0)->getPhi2(0,1,2); 
        // phi2_13 = d.event()->cluster(0)->getPhi2(0,1,3); 
        // phi2_21 = d.event()->cluster(0)->getPhi2(0,2,1); 
        // phi2_22 = d.event()->cluster(0)->getPhi2(0,2,2); 
        // phi2_23 = d.event()->cluster(0)->getPhi2(0,2,3); 
         phi2_01 = d.event()->cluster(0)->getPhi2(0,0,1); 
         phi2_10 = d.event()->cluster(0)->getPhi2(0,1,0); 
       //  phi3 = d.event()->cluster(0)->getPhi3(0); 
       //  phi4 = d.event()->cluster(0)->getPhi4(0); 

         phi2_skew = phi2; 
         phi2_greater = phi2; 
         phi2_slope = phi2; 

         double xsum = 0; 
         double ysum = 0; 
         std::vector<int> pix = d.event()->cluster(0)->getCluster(0); 
         const TH2 * img = d.event()->cluster(0)->getImage(); 

         int bx, by, bz; 
         for (unsigned p = 0; p < pix.size(); p++)
         {
            img->GetBinXYZ(pix[p],bx,by,bz); 
            xsum += img->GetXaxis()->GetBinCenter(bx); 
            ysum += img->GetYaxis()->GetBinCenter(by); 
         }

         double xc = xsum / pix.size(); 
         double yc = ysum / pix.size(); 

         double xdiff = 0, ydiff = 0; 
         for (unsigned p = 0; p < pix.size(); p++)
         {
            img->GetBinXYZ(pix[p],bx,by,bz); 
            xdiff += img->GetXaxis()->GetBinCenter(bx) - xc; 
            ydiff += img->GetYaxis()->GetBinCenter(by) - yc; 
         }

         phi_vec = atan2(-ydiff,-xdiff); 
         phi2_vec = phi2; 

         if (cos(phi_vec - phi2) < 0) 
         {
            phi2_vec = DmtpcMath::normalizeAngle(phi2+M_PI); 
         }
         

         TH1 * prof = d.event()->cluster(0)->projectClusterInterpolate(0,phi2, "bicubic"); 

         if (prof->GetSkewness() < 0)
         {
            phi2_skew = DmtpcMath::normalizeAngle(phi2+M_PI); 
         }

         if (prof->Integral(1,prof->GetNbinsX()/2) < prof->Integral(prof->GetNbinsX() - prof->GetNbinsX()/2, prof->GetNbinsX()))
         {
            phi2_greater = DmtpcMath::normalizeAngle(phi2+M_PI); 
         }

         TF1 pol1("f","pol1"); 
         prof->Fit(&pol1); 
         if (pol1.GetParameter(1) > 0)
         {
           phi2_slope =  DmtpcMath::normalizeAngle(phi2+M_PI); 
         }

         proj_length = prof->GetNbinsX(); 

         prof->Delete();  


         //phi5 = d.event()->cluster(0)->getPhi5(0,4,0,360,0.95); 
         cam = j; 
      }
    }

    tree->Fill(); 
  }

  tree->Write(); 
  out.Close(); 

}






