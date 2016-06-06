#include "DmtpcSkimDataset.hh"
#include "DmtpcMCDataset.hh"
#include "DmtpcMath.hh"
#include "DmtpcProjection.hh"
#include "MaxCamSRIM.hh"
#include "TFile.h" 
#include <iostream>

double pressure = 60;

int main (int nargs, char ** args) 
{


  if (nargs < 4) 
  {
    std::cout << "usage: MCFitTest mcfile skimfile outfile out_dir" << std::endl;; 
    return 1; 
  }



  DmtpcSkimDataset d; 
  d.openRootFile(args[2]); 

  DmtpcMCDataset mcd; 
  mcd.loadFile(args[1]); 

  const char * outdir = args[4]; 
  const char * srim_file = nargs > 5 ? args[5] : "SRIM_F_in_CF4_100Torr"; 
  MaxCamSRIM sr(srim_file); 

  sr.setPressure(pressure); 

  TFile out(args[3],"RECREATE"); 

  TTree * tree = new TTree("fit","fit"); 

  double lengthcal = 0.16; 
  map<std::string,double> gains; 
  gains["A80334"] = 10.1; 
  gains["A80333"] = 16.5; 
  gains["110121"] = 18.4; 
  gains["100534"] = 18.6; 
  double gain; 

//  TFile gmfile("4sh_gainmap.root"); 
  map<std::string,DmtpcGainMap*> gainmaps; 
  gainmaps["A80334"] = d.getGainMap("A80334"); 
  gainmaps["A80333"] = d.getGainMap("A80333"); 
  gainmaps["110121"] = d.getGainMap("110121"); 
  gainmaps["100534"] = d.getGainMap("100534"); 

  map<std::string,double> offsets; 
  offsets["A80334"] = 3.19215826044598927; 
  offsets["A80333"] = 7.15021941682858486e-02; 
  offsets["110121"] = 1.63512125006614562e+00 ; 
  offsets["100534"] = -1.52889415813968199; 




  double mcE,  mcRange, mcPhi,  mcTheta, mcZ,mcX,mcY;
  double mcInteg; 
  double x, y; 
  double  recoE, recoRange,  phi2,  fitPhiGlobal; 
  int transverseNdof;  
  int nfound, cam, run, track,event; 
  double offset,correct; 
  int proj_length; 
  int success; 
  DmtpcProjection::SRIMLineFitParams p; 


  tree->Branch("mcE",&mcE); 
  tree->Branch("mcX",&mcX); 
  tree->Branch("mcY",&mcY); 
  tree->Branch("mcInteg",&mcInteg); 
  tree->Branch("mcRange",&mcRange); 
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
  tree->Branch("fitPhiGlobal",&fitPhiGlobal); 
  tree->Branch("fitProbability",&p.prob); 
  tree->Branch("fitProbabilityError",&p.probError); 
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
  tree->Branch("gain",&gain); 
  tree->Branch("cam",&cam); 
  tree->Branch("track",&track); 
  tree->Branch("event",&event); 
  tree->Branch("run",&run); 
  tree->Branch("correct",&correct); 
  tree->Branch("offset",&offset); 
  


 
  for (int i = 0; i < d.nevents(); i++)
  {

    nfound = 0; 
    mcd.getEvent(i); 
    d.getEvent(i); 
    mcE= mcd.getParticleE() ; 
    mcInteg = mcd.getE(); 
    mcRange= mcd.getRange(); 
    mcPhi= mcd.getPhi(); 
    mcZ = mcd.getZ(); 
    mcX = mcd.getX(); 
    mcY = mcd.getY(); 
    mcTheta = mcd.getTheta(); 
    correct = mcPhi; 
    fitPhiGlobal = DmtpcMath::normalizeAngle(p.phi-correct); 

    if (mcInteg<=0) continue; 

    event = d.event()->eventNumber();
    run = d.event()->runNumber();

    for (int j = 0; j < d.event()->ncamera(); j++)
    {
      nfound += d.event()->ntracks(j); 

      for (int t = 0; t < d.event()->ntracks(j); t++)
      {
         recoE = d.event()->E(j,t);      
         x = d.event()->x(j,t);      
         y = d.event()->y(j,t);      
         recoRange = d.event()->range(j,t);      
         std::cout << "recoE: " << recoE << std::endl; 
         if (recoE == 0 )
         {
          std::cout << "huh?!? skipping track " << recoE << std::endl; 
          continue;

         }

         if (fabs(x-mcX) > 16 && fabs(y-mcY) >16)
         {
          std::cout << "huh?!? bad match " << std::endl; 
          continue;

         }

      //   phi1 = d.event()->cluster(0)->getPhi(0); 
         phi2 = d.event()->phi(j,t); 

//         ((MaxCamClusterImage*)d.event()->cluster(0))->morphologicalOperation(0,0,6); 

         cam = j; 
         track = t; 
         std::string serial = d.event()->cameraSerialNumber(cam).Data(); 
         offset = offsets[serial]; 
         gain = gains[serial]; 

         DmtpcProjection prof(d.event(),j,t,gainmaps[serial]); 
         proj_length = prof.getLongitudinalProfile()->GetNbinsX();  
         DmtpcProjection::SRIMLineFitParams * params= prof.doSRIMLineFit(&sr,lengthcal, gain,"IVD"); 
         prof.lineFitCanvas()->SaveAs(TString::Format("%s/%04d.%d.%d.png",outdir,event,cam,track)); 
         p = *params; 



         delete params; 
         out.cd(); 
         tree->Fill(); 
      }
    }

  }

  tree->Write(); 
  out.Close(); 

}







