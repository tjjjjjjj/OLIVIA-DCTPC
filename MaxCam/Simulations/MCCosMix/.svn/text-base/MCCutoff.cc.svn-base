#include "DmtpcDataset.hh"
#include "DmtpcEvent.hh"
#include "MaxCamImageTools.hh"
#include "MaxCamConfig.hh"
#include "MaxCamClusterImage.hh"

#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TObjString.h"
#include "TString.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TApplication.h"
#include "TPluginManager.h"
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

   TApplication theApp("App",&argc,argv);
   
   gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
					 "*", "TStreamerInfo",
					 "RIO","TStreamerInfo()");
   TString mcdataset = theApp.Argv(1);
   TString cosmicdataset = theApp.Argv(2);
   
   double blankFraction=0;
   if(theApp.Argc()>3)
   {
      blankFraction = atof(theApp.Argv(3));
   }
   TString outfilename = mcdataset;
   outfilename.ReplaceAll(".root","_cutoff.root");
   if(theApp.Argc()>4)
   {
      TString fname = outfilename.Replace(0,outfilename.Last('/')+1,"");
      outfilename = theApp.Argv(4);
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
	      << " in Monte Carlo \n";
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
   }

   gROOT->cd();


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
	 TH2F* bias = cosmic.getBiasFrame(CosmicCameras[(*it).first]+1);
	 TString name = "biasFrame";
	 name+=MCcameras[(*it).first]+1;
	 bias->SetName(name);
	 bias->SetTitle(name);
	 bias->Write(name);
	 if(cosmic.event()->ccdConfig(CosmicCameras[(*it).first])->digitizeOverscan)
	 {
	    cosmic.file()->cd();
	    TH2F* biasos = cosmic.getBiasFrameOverscan(CosmicCameras[(*it).first]+1);
	    name = "biasFrameOverscan";
	    name+=MCcameras[(*it).first]+1;
	    biasos->SetName(name+"l");
	    biasos->SetTitle(name);
	    out.file()->cd();
	    biasos->Write(name);
	 }
      }
   }

   out.setComment(mc.comment()->GetString());
   out.setKeyword(mc.keyword()->GetString());
   out.setLocation(mc.location()->GetString());
   
   //copy over the simulation data
/*   TTree* camera = (TTree*)mc.file()->Get("camera")->Clone("camout");
   TTree* chamber = (TTree*)mc.file()->Get("chamber");
   TTree* runInfo = (TTree*)mc.file()->Get("runInfo");
*/TTree* Simulation = (TTree*)mc.file()->Get("Simulation");/*



   out.file()->cd();
   camera->Write("camera");
   chamber->Write();
   runInfo->Write(); 
   Simulation->Write();*/


   
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
   int bincut[nCamMix];

   mixing->Branch("nCamMix",&nCamMixt,"nCamMix/I");
   mixing->Branch("runNumMix",&runNum,"runNumMix/I");
   mixing->Branch("eventNumMix",&eventNum,"eventNumMix/I");
   mixing->Branch("mixed",&mixed,"mixed[nCamMix]/O");
   mixing->Branch("bincut",&bincut,"bincut[nCamMix]/I");

   DmtpcEvent* ev;
   TH2S* himage;
   out.tree()->SetBranchAddress("event",&ev);
   TRandom3* rnd = new TRandom3();
   TObjArray* trueClust = 0;
   Simulation->SetBranchAddress("trueClusterArray",&trueClust);
   
   for(int i=0; i<mc.tree()->GetEntries(); i++)
   {
      if(i%100==0) cout << "Event: " << i << endl;
      
      mc.getEvent(i);
      cosmic.getEvent(i);
      Simulation->GetEntry(i);
      ev=mc.event();
      runNum=cosmic.event()->runNumber();
      eventNum=cosmic.event()->eventNumber()-1;
      int u=0;
      
      bool overscans = ev->overscan()->GetEntries();
      
      TH2S* hoverscan;
      
      if(rnd->Uniform()>blankFraction)
      {
	 for(it=MCcameras.begin(); it!=MCcameras.end(); it++)
	 {
	    
	    int mcindex = (*it).second;
	    int dataindex = CosmicCameras[(*it).first];


	    gROOT->cd();
	    himage = (TH2S*)ev->rawCcdData(mcindex);

	    int bin=rnd->Integer(himage->GetNbinsX())+1;
	    bincut[u]=bin;
	    int side = rnd->Integer(2);

	    if(side==0)
	    {
	       for(int xb=1; xb<bin; xb++)
	       {
		  for(int yb=1; yb<=himage->GetNbinsY(); yb++)
		  {
		     himage->SetBinContent(xb,yb,0);
		  }
	       }
	    }
	    else
	    {
	       for(int xb=bin+1; xb<himage->GetNbinsX(); xb++)
	       {
		  for(int yb=1; yb<=himage->GetNbinsY(); yb++)
		  {
		     himage->SetBinContent(xb,yb,0);
		  }
	       }
	    }

	    himage->Add(cosmic.event()->rawCcdData(dataindex));
	  	    
	    if(overscans)
	    {
	       hoverscan = (TH2S*)ev->rawOverscan(mcindex);
	       hoverscan->Add(cosmic.event()->rawOverscan(dataindex));
	    }
	    else
	    {								
	       hoverscan = (TH2S*)cosmic.event()->rawOverscan(dataindex);
	       new( (*ev->rawOverscan())[mcindex] ) TH2S(*hoverscan);
	    }
	    mixed[u]=true;
	    u++;
	    out.file()->cd();
	    
	 }
      }
      else
      {
	 for(it=MCcameras.begin(); it!=MCcameras.end(); it++)
	 {
	    
	    int mcindex = (*it).second;
	    int dataindex = CosmicCameras[(*it).first];
	    gROOT->cd();
	    himage = (TH2S*)ev->rawCcdData(mcindex);
	    himage->Clear();
	    himage->Add(cosmic.event()->rawCcdData(dataindex));

	    if(ev->overscan()->GetEntries())
	    {
	       hoverscan = (TH2S*)ev->rawOverscan(mcindex);
	       hoverscan->Clear();
	       hoverscan->Add(cosmic.event()->rawOverscan(dataindex));
	    }
	    else
	    {
	       hoverscan = (TH2S*)cosmic.event()->rawOverscan(dataindex);
	       new( (*ev->rawOverscan())[mcindex] ) TH2S(*hoverscan);
	    }
	    mixed[u]=false; u++;
	    out.file()->cd();
	 }

      }
      
      out.fill();
      mixing->Fill();
      ev->overscan()->Delete();
      ev->ccdData()->Delete();
      ev->ccdConfig()->Delete();
      ev->scopeData()->Delete();
      ev->scopeDataInfo()->Delete();
      ev->experimentConfig()->Delete();
      ev->mcTrack()->Delete();
      ev->mcCcdDigi()->Delete();
      
   }

   cosmic.closeRootFile();
   mc.closeRootFile();

   out.write();
   mixing->Write();
   out.file()->cd();
   
   return 0;
}
