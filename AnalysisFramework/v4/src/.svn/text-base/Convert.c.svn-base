#include <iostream>
#include <stdlib.h>
#include "TSystem.h"
#include "TFile.h"
#include "TH2.h"
#include "stdio.h"
#include "TTree.h"

#include "/home/raytheon/work/code/MaxCam/DmtpcDataset.hh"
#include "/home/raytheon/work/code/MaxCam/MaxCamConfig.hh"

using namespace std;

int Convert(TString,int,int);

int main(int argc, char** argv)
{

  TString infolder="/short/";
  //void *dir=gSystem->OpenDirectory(infolder);

  const char* file;
  while(1){
    cout << "\nCycling...\n" << endl;
    void *dir=gSystem->OpenDirectory(infolder);
    while(file=gSystem->GetDirEntry(dir)){
      if(strcmp(file,".")!=0 && strcmp(file,"..")!=0){
	cout << "Processing " << file << endl;
	gSystem->Sleep(100);
	TString filename(file);
	int i=Convert(filename,0,0);
	TString rmcomm("rm ");
	rmcomm+=infolder;
	rmcomm+=filename;
	gSystem->Exec(rmcomm);
	cout << "Remove command " << rmcomm << endl;
	void *dir = (int*)0;
	delete dir;
      } else gSystem->Sleep(100);
    }
  }
}

int Convert(TString filename, int runno, int eventno)
{
  //
  // Input file
  std::cout << "filename " << filename << std::endl;
  TString infile="/short/";
  infile+=filename;
  TFile inFile(infile);
  cout << "After TFile" << endl;
  while(inFile.IsZombie()){
    ///Return to top of function???
    cout << "Error opening file" << endl;
    return -1;
  }


  /*  while(inFile.IsOpen()==false){
     std::cout << "File not open" << std::endl;

     return -1;
     }*/

  inFile.ls();
  //
  // Load 
  TH2S* ccd=(TH2S*) inFile.Get("ccd_0;1");
  ccd->SetDirectory(0);
  TH2F* biasFrame1=(TH2F*) inFile.Get("biasFrame1;1");
  biasFrame1->SetDirectory(0);
  TFile out("junk.root","recreate");
  ccd->Write();
  biasFrame1->Write();
  inFile.Close();
  out.ls();
  //
  // Setup output file
  TString outfile="/temp/";
  outfile+=filename;
  //outfile+=".root";
  //
  // Setup dataset
  DmtpcDataset* a=new DmtpcDataset();
  a->createRootFile(outfile,"recreate");
  a->setComment("Single event");
  a->setLocation("44-013");
  a->comment()->Write();
  a->keyword()->Write();
  a->location()->Write();
  //
  // Setup configuration
  MaxCamConfig* ccdconfig=new MaxCamConfig("ccdconfig","Andor camera configuration");
  ccdconfig->cameraID=0;
  ccdconfig->row_width=1024;
  ccdconfig->img_rows=1024;
  ccdconfig->hbin=1024;
  ccdconfig->vbin=1024;
  ccdconfig->ul_x=0;
  ccdconfig->ul_y=0;
  ccdconfig->lr_x=1024;
  ccdconfig->lr_y=1024;
  ccdconfig->CCDTemp=-20;
  ccdconfig->CCDTempSet=-20;
  ccdconfig->exposureTime=1000;
  ccdconfig->frameType=0;
  ccdconfig->nFlushes=-1;
  ccdconfig->bitDepth=65535;
  ccdconfig->daqTime=-1;
  ccdconfig->serialNumber="0";
  //
  // Write bias frame
  biasFrame1->Write("biasFrame1");
  //
  // Write event
  DmtpcEvent* ev=a->event();
  ev->setEventNumber(eventno);
  ev->setRunNumber(runno);
  new ((*a->event()->rawCcdData())[0]) TH2S(*ccd);
  new((*a->event()->ccdConfig())[0]) MaxCamConfig(*ccdconfig);
  a->fill();
  a->write();
  a->closeRootFile();
  //
  // Test section
  TFile in("test.root");
  in.ls();
  TTree* t=(TTree*) in.Get("dmtpc");
  //t->Print("all");
  std::cout << "Entries " << t->GetEntries() << endl;
  in.Close();
}
