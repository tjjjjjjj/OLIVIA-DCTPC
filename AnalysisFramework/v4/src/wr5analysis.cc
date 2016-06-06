#include "DmtpcSkimDataset.hh"
#include "DmtpcSkimEvent.hh"
#include "DmtpcEventTable.hh"
#include "MaxCamImageTools.hh"

#include "DmtpcPulse.hh"
#include "DmtpcWaveformTools.hh"
#include "FastWfVector.hh"
#include "CspWfVector.hh"
#include "FastPulse.hh"
#include "CspPulse.hh"
#include "FastWaveform.hh"
#include "CspWaveform.hh"

#include "AnalysisCut.hh"
#include "AnalysisConfig.hh"
#include "recoilEnergy.hh"
#include "MultiVariate.hh"

#include "TF1.h"
#include "TGraphErrors.h"

#include <string.h>
#include <iostream>


using namespace std; 
using namespace MultiVariate;


AnalysisConfig* cfg=0;
TGraphErrors* gGainVar=0;


// Constants for energy corrections based on 
// differing focuses for the two cameras
double cplus1_100439 = 0.772553;
double slope_100439 = 0.000459378;
double cplus1_081264 = 0.857977;
double slope_081264 = 0.000416834;

double getEv(DmtpcSkimEvent * e, int c,int t)
{
   const char* sn= e->cameraSerialNumber(c);
   return e->EGainMap(c,t)/cfg->getEnergyCal(sn); 
}


//Get recoil energy—correct for visible vs recoil
double getEr(DmtpcSkimEvent *e,int c, int t)
{
  double Er = RecoilEnergy::visibleToRecoil(getEv(e,c,t)); 
  if(e->cameraSerialNumber(c) == "100439")
     return (sqrt(cplus1_100439*cplus1_100439+4*slope_100439*Er)-cplus1_100439)/(2*slope_100439);
  else if(e->cameraSerialNumber(c) == "081264")
     return (sqrt(cplus1_081264*cplus1_081264+4*slope_081264*Er)-cplus1_081264)/(2*slope_081264);
  else
     return Er;
  
}

// get range using fn calculated from MC
double getCalibratedRange(DmtpcSkimEvent * e, int c, int t, TF1* fn)
{
   const char* sn= e->cameraSerialNumber(c);
   return fn->Eval(e->range(c,t))*cfg->getRangeCal(sn);
}


//this may not be current
bool
passChargeCuts(DmtpcSkimEvent* ev, int c, int t, double* params)
{
  //Assumes 1 pulse / trace

   FastWfVector* tMv = (FastWfVector*) ev->waveform_vectors()->At(3);
   FastWfVector* bMv = (FastWfVector*) ev->waveform_vectors()->At(2);
   CspWfVector* tAv = (CspWfVector*) ev->waveform_vectors()->At(0);
   CspWfVector* bAv = (CspWfVector*) ev->waveform_vectors()->At(1);
   
   double Eccd = ev->E(c,t);

   int nt = tMv->size();

   for (int i = 0; i<nt; i++)
   {
      
//      FastWaveform tM = tMv->at(i);
//      FastWaveform bM = bMv->at(i);
//      CspWaveform tA = tAv->at(i);
//      CspWaveform bA = bAv->at(i);
      
      
      //Basic cuts
      bool passBasic =  bAv->at(i).getBase() > params[1] && bAv->at(i).getBase()<params[2] &&
	 tAv->at(i).getBase()>params[3] && tAv->at(i).getBase() < params[4] &&
	 tMv->at(i).getBase()>params[5] && tMv->at(i).getBase()< params[6] && 
	 bMv->at(i).getBase()>params[7] && bMv->at(i).getBase()< params[8] && 
	 tMv->at(i).getRMS()<params[9] && bMv->at(i).getRMS() < params[10] &&
	 tAv->at(i).getWfMax()<params[11]&&bAv->at(i).getWfMax()<params[12] &&
	 tMv->at(i).getWfMax()<params[13] && bMv->at(i).getWfMax()<params[14];
      if (!passBasic) continue;

      bool passTopBottom;
      if (params[0]>0){
//Top cuts
	 double Escope = tAv->at(i)[0].getPeak()*1000;
	 passTopBottom =   
	    tAv->at(i)[0].getPeak() > bAv->at(i)[0].getPeak() &&
	    tMv->at(i)[0].getIntegral() > bMv->at(i)[0].getIntegral() &&
	    tAv->at(i)[0].getRise25() - tAv->at(i)[0].getRise75() < params[15] &&
	    tMv->at(i)[0].getIntegral() > params[16] &&
	    tAv->at(i)[0].getPeak() > params[17] &&
	    //Pulse shape cuts
	    tMv->at(i)[0].getFastRise25()-tMv->at(i)[0].getFastRise90()<params[18] &&
	    tMv->at(i)[0].getFastRise10()-tMv->at(i)[0].getFastRise90()<params[19] &&
	    tMv->at(i)[0].getSlowFall50()-tMv->at(i)[0].getSlowFall90()<params[20] &&
	    tMv->at(i)[0].getSlowPeakTime()-tMv->at(i)[0].getFastPeakTime() < params[21] &&
	    tMv->at(i)[0].getTroughTime()-tMv->at(i)[0].getFastPeakTime()<params[22] &&
	    tMv->at(i)[0].getFastPeak() > params[23] * tMv->at(i)[0].getSlowPeak() &&
	    //Charge/track matching
	    fabs(Escope - params[24] - params[25]*Eccd) < params[26];
	 
	 
      }else{
	 double Escope = tAv->at(i)[0].getPeak()*1000;
	 passTopBottom =   
	    bAv->at(i)[0].getPeak() > tAv->at(i)[0].getPeak() &&
	    bMv->at(i)[0].getIntegral() > tMv->at(i)[0].getIntegral() &&
	    bAv->at(i)[0].getRise25() - bAv->at(i)[0].getRise75() < params[15] &&
	    bMv->at(i)[0].getIntegral() > params[16] &&
	    bAv->at(i)[0].getPeak() > params[17] &&
	    bMv->at(i)[0].getFastRise25()-bMv->at(i)[0].getFastRise90()<params[18] &&
	    bMv->at(i)[0].getFastRise10()-bMv->at(i)[0].getFastRise90()<params[19] &&
	    bMv->at(i)[0].getSlowFall50()-bMv->at(i)[0].getSlowFall90()<params[20] &&
	    bMv->at(i)[0].getSlowPeakTime()-bMv->at(i)[0].getFastPeakTime() < params[21] &&
	    bMv->at(i)[0].getTroughTime()-bMv->at(i)[0].getFastPeakTime()<params[22] &&
	    bMv->at(i)[0].getFastPeak() > params[23] * bMv->at(i)[0].getSlowPeak() &&
	    fabs(Escope - params[24] - params[25]*Eccd) < params[26];
      }
      
      if (passTopBottom) return true;
      
   }//loop over traces
   return false;
   
}

bool passRegionCut(DmtpcSkimEvent * ev, int c, int t,double* param)
{
  return 
       ev->x(c,t)>param[0]&&ev->x(c,t)<param[1]&&
       ev->y(c,t)>param[0]&&ev->y(c,t)<param[1]; 
} 

bool passRangeCut(DmtpcSkimEvent * ev, int c, int t,double* param)
{
   const char* sn = ev->cameraSerialNumber(c);
   return ev->range(c,t) * cfg->getRangeCal(sn) < param[0]; 
}

bool passSparkCut(DmtpcSkimEvent * ev, int c, int t, double* param)
{
   if( ev->spark(c))
      return 0;
   else
   {
      if(ev->lastspark(c) < param[0])
	 return 0;
      else
	 return 1;
   }
}

bool passNonZeroCut(DmtpcSkimEvent * ev, int c, int t, double* param)
{

   return ev->E(c,t)>param[0] && 
      ev->range(c,t) > param[0] && 
      ev->cluster_rms(c,t) > param[0] && 
      ev->neighbors(c,t) >= param[0] &&  
      ev->maxpixel(c,t) > param[0] && 
      ev->cluster_mean(c,t) > param[0] && 
      ev->npixel(c,t) > param[0]; 
}

bool passWormCut(DmtpcSkimEvent* ev, int c, int t, vector<void*> param)
{				
   MultiVariateResult * mvr = (MultiVariateResult*)param[2*c];
   double mvcut = *((double*)param[2*c+1]);
   
   if(ev->range(c,t)>0)
   {	
//      cout << mvr->getClassifier(ev,c,t) << "," << mvcut << endl;


      return mvr->getClassifier(ev,c,t)>mvcut;
   }
   return 0;
}

bool passBasicWormCut(DmtpcSkimEvent* ev, int c, int t, double* param)
{
   
   return (ev->maxpixel(c,t)<param[0] && ev->cluster_rms(c,t)<param[1] 
	   && ev->maxpixel(c,t)/ev->E(c,t)<param[2] 
	   /*&& ev->cluster(c)->getClusterPerimeter(t,false)/4.0<ev->npixel_red(c,t)*/);

}

bool passEdgeCut(DmtpcSkimEvent* ev, int c, int t, double* param)
{

   return !ev->edge(c,t);

}

bool passEnergyCut(DmtpcSkimEvent* ev, int c, int t, double* param)
{
   double E_r = getEr(ev,c,t);

   return E_r>param[0] && E_r<param[1];

}

bool passRBICut(DmtpcSkimEvent* ev, int c, int t, double* param)
{
   return !ev->isRBI(c,t, (int)param[0], (int)param[1]);
}


//Takes arguments of input skim file, cfg file, path to output file

int main (int nargs, char** args)
{

  if(nargs>2)
     cfg = new AnalysisConfig(args[2]);
  else
     cfg = new AnalysisConfig("cfg/AnalysisDefault.cfg");
  //cfg->print();

  cout << cfg->getTree() << endl;

//open dataset
  DmtpcSkimDataset ds; 
  ds.openRootFile(args[1]);
//  ds.tree()->SetBranchStatus("*trigger*",0);
//  ds.tree()->SetBranchStatus("*clusters*",0);

  //set up output file
  TString cutfilename="";
  if(nargs>3)
     cutfilename=args[3];
  else
     cutfilename = "cuts.root";
  ds.getEvent(0);
  const int ncam = ds.event()->ncamera();
  int ncamera = ncam;

  MultiVariateResult* mvr[ncam];
  double mvcut[ncam];
  TF1* rangefn[ncam];
  TFile* rangefilename = new TFile(cfg->getRangeCalFileName());

  //Draw out MV for worm cut from file
  for(int i=0; i<ncam; i++)
  {
     mvr[i] = new MultiVariateResult();
     mvr[i]->setResult(cfg->getWormFileName(ds.event()->cameraSerialNumber(i))); 
     mvcut[i] = cfg->getWormCutVal(ds.event()->cameraSerialNumber(i));
     rangefn[i] = (TF1*)rangefilename->Get(ds.event()->cameraSerialNumber(i));
  }


  //Setup cuts and set parameters
  AnalysisCut* nonZero = new AnalysisCut("nonZero",passNonZeroCut,1);
  nonZero->setParameter(0,0);
  AnalysisCut* region = new AnalysisCut("region",passRegionCut,2);
  region->setParameter(0,40);
  region->setParameter(1,984);
  AnalysisCut* range = new AnalysisCut("range",passRangeCut,1);
  range->setParameter(0,6.23545);
  AnalysisCut* worm = new AnalysisCut("worm",passWormCut);
  vector<void*> wormparams;
  for(int i=0; i<ncam; i++)
  {
     wormparams.push_back(mvr[i]);
     wormparams.push_back(&mvcut[i]);
  }
  worm->setParameters(wormparams);
  AnalysisCut* edge = new AnalysisCut("edge",passEdgeCut,0);
  AnalysisCut* energy = new AnalysisCut("energy",passEnergyCut,2);
  energy->setParameter(0,50);
  energy->setParameter(1,300);
  AnalysisCut* basicWorm = new AnalysisCut("basicWorm",passBasicWormCut,3);
  basicWorm->setParameter(0,400);
  basicWorm->setParameter(1,100);
  basicWorm->setParameter(2,0.225);
  AnalysisCut* rbi = new AnalysisCut("rbi",passRBICut,2);
  rbi->setParameter(0,2);
  rbi->setParameter(1,3);
  AnalysisCut* spark = new AnalysisCut("spark",passSparkCut,1);
  spark->setOnce(true);
  spark->setParameter(0,5);//number to veto after; 5 for 1s, 3 for 2s, 2 for 4s, 1 for 5s
  AnalysisCut* charge = new AnalysisCut("charge",passChargeCuts,27);
  charge->setParameter(0,1);
  charge->setParameter(1,-0.0025);
  charge->setParameter(2,0.015);
  charge->setParameter(3,-0.005);
  charge->setParameter(4,0.005);
  charge->setParameter(5,-0.05);
  charge->setParameter(6,-0.041);
  charge->setParameter(7,-1);
  charge->setParameter(8,1);
  charge->setParameter(9,0.004);
  charge->setParameter(10,0.0023);
  charge->setParameter(11,0.098);
  charge->setParameter(12,0.098);
  charge->setParameter(13,0.39);
  charge->setParameter(14,0.39);
  charge->setParameter(15,3e-6);
  charge->setParameter(16,1.5);
  charge->setParameter(17,0.00625);
  charge->setParameter(18,40e-9);
  charge->setParameter(19,60e-9);
  charge->setParameter(20,200e-9);
  charge->setParameter(21,400e-9);
  charge->setParameter(22,100e-9);
  charge->setParameter(23,0.7);
  charge->setParameter(24,4.01);
  charge->setParameter(25,0.0206);
  charge->setParameter(26,15);

//   TFile* gainvarfile = new TFile("series.root");
//   TGraphErrors* gGainVar = (TGraphErrors*)gainvarfile->Get("gain_v_time");
//   cfg->setEnergyCal("081264",gGainVar->Eval(ds.event()->runNumber()));
//   cout << cfg->getEnergyCal("081264");
     
  double _Erecoil[ncam][15];
  double _calibratedrange[ncam][15];
  double _wormclass[ncam][15];

  TFile* cutfile = new TFile(cutfilename,"RECREATE");
  TTree* cuts = new TTree("cuts","cuts");
  cuts->Branch("nonZero","AnalysisCut",&nonZero,128000,1);
  cuts->Branch("region","AnalysisCut",&region,128000,1);
  cuts->Branch("range","AnalysisCut",&range,128000,1);
  cuts->Branch("worm","AnalysisCut",&worm,128000,1);
  cuts->Branch("edge","AnalysisCut",&edge,128000,1);
  cuts->Branch("energy","AnalysisCut",&energy,128000,1);
  cuts->Branch("basicWorm","AnalysisCut",&basicWorm,128000,1);
  cuts->Branch("rbi","AnalysisCut",&rbi,128000,1);
  cuts->Branch("spark","AnalysisCut",&spark,128000,1);
  cuts->Branch("charge","AnalysisCut",&charge,128000,1);
  cuts->Branch("ncamera",&ncamera,"ncamera/I");
  cuts->Branch("_Erecoil",&_Erecoil,"_Erecoil[ncamera][15]/D");
  cuts->Branch("_calibratedrange",&_calibratedrange,"_calibratedrange[ncamera][15]/D");
  cuts->Branch("_wormclass",&_wormclass,"_wormclass[ncamera][15]/D");

  for(int i=0; i<ds.nevents(); i++)
  //for(int i=0; i<100; i++)
  {
     if(i%100==0) cout << "Event: " << i << endl;
     ds.getEvent(i);

     nonZero->evaluateEvent(ds.event(cfg->getTree()));
     region->evaluateEvent(ds.event(cfg->getTree()));
     range->evaluateEvent(ds.event(cfg->getTree()));
     worm->evaluateEvent(ds.event(cfg->getTree()));
     edge->evaluateEvent(ds.event(cfg->getTree()));
     energy->evaluateEvent(ds.event(cfg->getTree()));
     basicWorm->evaluateEvent(ds.event(cfg->getTree()));
     rbi->evaluateEvent(ds.event(cfg->getTree()));
     spark->evaluateEvent(ds.event(cfg->getTree()));
//     charge->evaluateEvent(ds.event("skim"));

     for(int c=0; c<ncam; c++)
     {
	for(int t=0; t<15; t++)
	{
	   _Erecoil[c][t] = getEr(ds.event(cfg->getTree()),c,t);
	   _calibratedrange[c][t] = getCalibratedRange(ds.event(cfg->getTree()),c,t,rangefn[c]);
	   if(ds.event(cfg->getTree())->range(c,t)>0)
	   {
	      _wormclass[c][t] = mvr[c]->getClassifier(ds.event(cfg->getTree()),c,t);
//	      cout << i << "," << c << "," << t << ";" << _wormclass[c][t] << endl;
	   }
	   else
	      _wormclass[c][t]=-10;
	}
     }
     cuts->Fill();
  }
  
  cuts->Write();

  cutfile->Close();
  
}
