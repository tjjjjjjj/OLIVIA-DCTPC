#include "DmtpcDataset.hh"
#include "DmtpcEvent.hh"
#include "MaxCamImageTools.hh"
#include "MaxCamConfig.hh"

#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TObjString.h"
#include "TString.h"
#include <string>
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TColor.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TKey.h"

#include <vector>
#include <map>

using namespace std;

int main(int argc, char *argv[])
{

   
   TString mcdataset = argv[1];
   TString cosmicdataset = argv[2];
   
   double blankFraction=0;
   if(argc>3)
   {
      blankFraction = atof(argv[3]);
   }
   TString outfilename = mcdataset;
   outfilename.ReplaceAll(".root","_mixed.root");
   if(argc>4)
   {
      TString fname = outfilename.Replace(0,outfilename.Last('/')+1,"");
      outfilename = argv[4];
      outfilename+="/"+fname;
   }

   
   //open MC dataset, get SN of cameras
   DmtpcDataset mc;
   mc.openRootFile(mcdataset);
   DmtpcDataset cosmic;
   cosmic.openRootFile(cosmicdataset);

   if(mc.tree()->GetEntries() > cosmic.tree()->GetEntries())
   {
      cout << "Too many MC entries for this run. Quitting" << endl;
      return 0;
   }

   gROOT->cd();
   const int nMCcam = mc.getNcameras();
   mc.getEvent(0);
   map<TString,int> MCcameras;

   //find out the purported SN in MC
   for(int u = 0; u<nMCcam; u++)
   {
      if(mc.event()->ccdConfig(u)->serialNumber != "")
      {
         MCcameras[mc.event()->ccdConfig(u)->serialNumber]=u;
         cout << "Found Camera " << mc.event()->ccdConfig(u)->serialNumber
              << " in Monte Carlo at " <<u << std::endl; 
      }
      else
         cout << "Camera " << u << " does not have an associated serial number. Not processing this camera.\n";
   }
   

   if(MCcameras.size()==0)
   {
      cout << "No MC cameras to be processed. Quitting." << endl;
      return 0;
   }


   const int nCosmicCam = cosmic.getNcameras();
   cosmic.getEvent(0);
   map<TString,int> CosmicCameras;

   //find out the SN in data
   for(int u = 0; u<nCosmicCam; u++)
   {
      CosmicCameras[cosmic.event()->ccdConfig(u)->serialNumber]=u;
      cout << "Found Camera " << cosmic.event()->ccdConfig(u)->serialNumber
      << " in Data at " <<u << std::endl; 
   }


   DmtpcDataset out;
   out.createRootFile(outfilename,"recreate");

   //check if all MC cameras have a corresponding data camera; save bias frame
   map<TString,int>::iterator it;
   for(it=MCcameras.begin(); it!=MCcameras.end(); it++)
   {
      if(!CosmicCameras.count((*it).first))
      {
         cout << "No data camera match to serial Number in Monte Carlo" << endl;
         return 0;
      }
      else
      {
         std::cout << "Using bias frame from data camera " << CosmicCameras[(*it).first] << " for camera " << (*it).second << std::endl; 
         TH2F* bias = cosmic.getBiasFrame(CosmicCameras[(*it).first]+1);
         TString name = "biasFrame";
         name+= (*it).second+1; 
         bias->SetName(name); 
         bias->SetTitle(name); 
         out.file()->cd(); 
         bias->Write(name);
         if(cosmic.event()->ccdConfig(CosmicCameras[(*it).first])->digitizeOverscan)
         {
            TH2F* biasos = cosmic.getBiasFrameOverscan(CosmicCameras[(*it).first]+1);
            name = "biasFrameOverscan";
            name+=(*it).second+1;
            biasos->Write(name);
         }
      }
   }

   out.setComment(mc.comment()->GetString());
   out.setKeyword(mc.keyword()->GetString());
   out.setLocation(mc.location()->GetString());
   
   
   TIter next(mc.file()->GetListOfKeys());  
   while (TKey * k = (TKey*) next())
   {
      cout << k->GetName() << endl;
      cout << k->GetClassName() << endl;

      TString name = k->GetName();
      TString classname = k->GetClassName();
      if (name.Contains("dmtpc") || name.Contains("biasFrame") 
          || name=="keyword" || name=="location"
          || name=="comment")
      {
         cout << "Not writing: " << k->GetName() << endl;
         continue; 
      }
      out.file()->cd();
      if(classname=="TTree")
      {
         TTree* tree = (TTree*)mc.file()->Get(k->GetName());
         TTree* copy = tree->CloneTree();
         out.file()->cd(); 
         copy->Write(k->GetName());
      }
      else
      {
         TObject * copy = k->ReadObj();
         copy->Write(k->GetName()); 
      }
      mc.file()->cd();
   } 


   out.file()->cd();
   //create tree for mixing info
   TTree* mixing = new TTree("mixing","mixing");
   
   const int nCamMix = MCcameras.size();
   int nCamMixt = nCamMix;
   int runNum;
   int eventNum;
   bool mixed[nCamMix];

   mixing->Branch("nCamMix",&nCamMixt,"nCamMix/I");
   mixing->Branch("runNumMix",&runNum,"runNumMix/I");
   mixing->Branch("eventNumMix",&eventNum,"eventNumMix/I");
   mixing->Branch("mixed",&mixed,"mixed[nCamMix]/O");

   TRandom3* rnd = new TRandom3();


   DmtpcEvent * ev = out.event(); 
   ev->rawCcdData()->ExpandCreate(nCamMix); 
   ev->rawOverscan()->ExpandCreate(nCamMix); 
   ev->ccdConfig()->ExpandCreate(nCamMix); 
   ev->experimentConfig()->ExpandCreate(mc.event()->experimentConfig()->GetEntries()); 
   for (int i = 0; i < mc.event()->experimentConfig()->GetEntries(); i++)
   {
     MaxCamChannel * cfg =  (MaxCamChannel*) (*mc.event()->experimentConfig())[i]; 
     new ((*ev->experimentConfig())[i]) MaxCamChannel(*cfg); 
   }


   for(int i=0; i<mc.tree()->GetEntries(); i++)
   {
      if(i%100==0) cout << "Event: " << i << endl;
      
      cosmic.getEvent(i);
      mc.getEvent(i);
      runNum=cosmic.event()->runNumber();
      eventNum=cosmic.event()->eventNumber()-1;

      int u=0;

      out.file()->cd(); 
      ev = out.event(); 
      ev->setEventNumber(mc.event()->eventNumber()); 
      ev->setRunNumber(mc.event()->eventNumber()); 


      bool MCHasOverscan = mc.event()->overscan()->GetEntries(); 
      bool signal = rnd->Uniform() > blankFraction; 

      
       for(it=MCcameras.begin(); it!=MCcameras.end(); it++)
       {
          
          int mcindex = (*it).second;
          int dataindex = CosmicCameras[(*it).first];


          MaxCamConfig* cfg = mc.event()->ccdConfig(mcindex); 
          ev->ccdConfig(mcindex)->~MaxCamConfig(); 
          new (ev->ccdConfig(mcindex)) MaxCamConfig(*cfg); 

          TH2S * cosmicimg = (TH2S*)  cosmic.event()->rawCcdData(dataindex); 
          ev->rawCcdData(mcindex)->~TH2(); 
          new(ev->rawCcdData(mcindex)) TH2S(*cosmicimg);

          if (signal) ev->rawCcdData(mcindex)->Add(mc.event()->rawCcdData(mcindex)); 


          TH2S * cosmicover = (TH2S*) cosmic.event()->rawOverscan(dataindex); 

          ev->rawOverscan(mcindex)->~TH2(); 
          new (ev->rawOverscan(mcindex)) TH2S(*cosmicover); 

          if(MCHasOverscan && signal)
          {
              ev->rawOverscan(mcindex)->Add(mc.event()->rawOverscan(mcindex)); 
          }
       }

     
            
      out.file()->cd(); 
      out.fill();
//      std::cout << "filled" << std::endl; 
      mixing->Fill();
   }


   out.file()->cd();
   out.write(); 
   mixing->Write();
   cosmic.closeRootFile();
   mc.closeRootFile();
   return 0;

}
