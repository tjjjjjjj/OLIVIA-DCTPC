#include "../../../MaxCam/DmtpcSkimDataset.hh"
#include "../../../MaxCam/DmtpcSkimEvent.hh"
#include "../../../MaxCam/DmtpcDataset.hh"
#include "TTree.h"
#include "TString.h"
#include <iostream>

int main(int argn, char ** args)
{
	DmtpcDataset ds;
	ds.openRootFile(args[1]);

  for (int i = 0; i < ds.tree()->GetEntries();i++)
	{
		ds.getEvent(i);
    std::cout << i << std::endl; 
	}
}
