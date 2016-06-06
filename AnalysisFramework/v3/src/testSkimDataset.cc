#include "../../../MaxCam/DmtpcSkimDataset.hh"
#include "../../../MaxCam/DmtpcSkimEvent.hh"
#include "../../../MaxCam/DmtpcDataset.hh"
#include "TTree.h"
#include "TString.h"
#include <iostream>

int main(int argn, char ** args)
{
  bool load = true;
  if (argn==3 && args[2][0]=='f') load = false; 
	DmtpcSkimDataset ds;
	ds.openRootFile(args[1]);

  ds.loadDmtpcEvent(load); 
	for (int i = 0; i < ds.nevents(); i++)
	{
		ds.getEvent(i);
    std::cout << i << std::endl; 
    if (load)
    {
      std::cout << "Pointer to orig_event: " << ds.orig_event() << std::endl; 
      std::cout << "Pointer to ccdConfig: " << ds.orig_event()->ccdConfig() << std::endl; 
      std::cout << "Pointer to ccdData: " << ds.orig_event()->ccdData() << std::endl; 
      std::cout << "Pointer to scopeData: " << ds.orig_event()->scopeData() << std::endl; 
      std::cout << "Pointer to scopeDataInfo: " << ds.orig_event()->scopeDataInfo() << std::endl; 
      std::cout << "Pointer to mcTrack: " << ds.orig_event()->mcTrack() << std::endl; 
      std::cout << "Pointer to mcCcdDigi: " << ds.orig_event()->mcCcdDigi() << std::endl; 
    }
	}
}
