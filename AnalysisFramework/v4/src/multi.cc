#include "MultiVariate.hh"
#include <vector> 
#include <string> 
#include <iostream>
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"
#include <stdlib.h>

#include <sys/stat.h>
#include <sys/types.h>

/* Example program for using multivariate framework.  */ 


const time_t min_m_time= 1313769600;

void print_usage()
{

  std::cout << "./multi cam_id worm_file_list mc_file_list output_name [tree=skim] [mccam=0]" << std::endl;  
}

int main (int nargs, char ** args)
{


  char * tree = "skim"; 
  char * mccam = 0; 
  if (nargs < 5)
  {
    print_usage(); 
    return 0; 
  }

  if (nargs  > 5)
  {
    tree = args[5]; 
  }

  if (nargs > 6)
  {
    mccam = args[6]; 
  }

  char * cam_id = args[1]; 

  std::vector<std::string> * signal = new std::vector<std::string>; 
  std::vector<std::string> * worms = new std::vector<std::string>; 


  std::string line; 
  std::ifstream worms_file(args[2]); 
  //TChain worms_test(tree); 
  int counter = 0; 
  while (!worms_file.eof())
  {
    std::getline(worms_file,line);
    if (line.length() > 0)
    {

      struct stat statbuf; 
      if (stat(line.c_str(), &statbuf) == -1)
      {
        std::cout << line << " not found" << std::endl;
        continue; 
      }

      if (statbuf.st_mtime < min_m_time) 
      {
        std::cout << line << " too old" << std::endl;
        continue; 
      }

        worms->push_back(line); 
    }
  }
  worms_file.close();

  std::ifstream signal_file(args[3]); 
 // TChain sig_test(tree); 
 

  while(!signal_file.eof())
  {
    std::getline(signal_file,line);
    if (line.length() > 0)
    {
        signal->push_back(line); 
    }
  }

  signal_file.close(); 


  TString name(args[4]); 
  TString tmvafile = TString("TMVA_")+name+TString(".root");

  TFile * f = new TFile(tmvafile,"RECREATE"); 

  TTree * signal_ch = MultiVariate::buildMVTree("signal",f,signal,-1,mccam,tree,32); 
  TTree * worms_ch = MultiVariate::buildMVTree("worms",f,worms,-1,cam_id,tree);

  std::cout << "signal tree built" << std::endl; 
  std::cout << signal_ch->GetEntries() << std::endl; 
  MultiVariate::runMVSimple(strdup(name.Data()),f,signal_ch,worms_ch,"E<2000","E<2000");
  f->Close(); 
}
