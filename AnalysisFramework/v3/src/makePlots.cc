#include "TH2F.h"
#include "TCanvas.h"
#include "../../../MaxCam/DmtpcSkimDataset.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include "recoilEnergy.hh"
#include "TApplication.h"
#include "style.hh"
#include "TLine.h"
#include <strings.h>

/* This file is confusing... I'm using these macros 
 * to make it a little easier to add plots.  */

/* Declare all the plots here using the format 
 * TH2*F plot_var_name[num_datasets] */
#define DECL_PLOTS                                             \
  TH2F * range_v_energy[num_datasets];                         \
  TH2F * clusterrms_v_energy[num_datasets];                    \
  TH1F * cygang[num_datasets]; 

/* Here you can initialize each plot with its title and stuff 
 *  f is the counter for each dataset
 *  infile is the name of the dataset inputfile
 * */
#define INIT_PLOTS                                             \
    range_v_energy[f]= new TH2F(mkname("range_v_energy",infile), \
                                       mktitle("Range vs. Recoil Energy",infile),\
                                       100, 0,200,\
                                       100,0,20);\
    range_v_energy[f]->SetXTitle("Recoil Energy");\
    range_v_energy[f]->SetYTitle("Range"); \
    range_v_energy[f]->SetTitle(mktitle("Range vs. Recoil Energy",infile));   \
                                                                        \
    clusterrms_v_energy[f]= new TH2F(mkname("clusterrms_v_energy",infile), \
                                       mktitle("ClusterRMS vs. Recoil Energy",infile),\
                                       100, 0,200,\
                                       100,0,100);\
    clusterrms_v_energy[f]->SetXTitle("Recoil Energy");\
    clusterrms_v_energy[f]->SetYTitle("ClusterRMS"); \
    clusterrms_v_energy[f]->SetTitle(mktitle("ClusterRMS vs. Recoil Energy",infile));   \
    \
    cygang[f] = new TH1F(mkname("cygnus_angle",infile), \
                         mktitle("Cygnus Angle",infile), \
                         10,-1,1); \
    cygang[f]->SetTitle(mktitle("Cygnus Angle",infile));  

    


/* This is done inside loops to fill the plots */
#define FILL_PLOTS                                             \
 double e=RecoilEnergy::getRecoilEnergy(ev,c,t);               \
 double range = ev->range(c,t);                                \
 range = c == 0 ? range*0.143 : range*0.156;                   \
 range_v_energy[f]->Fill(e,range);                             \
 clusterrms_v_energy[f]->Fill(e,ev->cluster_rms(c,t));         \
 cygang[f]->Fill(ev->cygnus_angle(c,t)); 

/*Define your canvases here */
#define MAKE_CANVAS       \
  TCanvas canvas1("canvas1", "plot"); \
  canvas1.Divide(num_datasets,1); \
  TCanvas canvas2("canvas2", "plot"); \
  canvas2.Divide(num_datasets,1); \
  TCanvas canvas3("canvas3", "plot"); \
  canvas3.Divide(num_datasets,1); \
  for (int i = 0; i < num_datasets; i++)\
  {\
    canvas1.cd(i); \
    range_v_energy[i]->DrawCopy("colz"); \
    canvas2.cd(i); \
    clusterrms_v_energy[i]->DrawCopy("colz"); \
    canvas3.cd(i); \
    cygang[i]->DrawCopy(); \
  }


static char * mkname(char * name, char * infile)
{
   char * ret = (char*)malloc(strlen(name) + strlen(infile) + 2); 
   strcpy(ret,name);
   strcat(ret,"_"); 
   strcat(ret,infile); 
   return ret;  
}

static char * mktitle(char * title, char * infile)
{
  char * ret = (char*)malloc(strlen(title) + strlen(infile) + 4); 
  strcpy(ret,title); 
  strcat(ret," ("); 
  strcat(ret,infile); 
  strcat(ret,")"); 
  return ret; 
}

int main (int nargs, char ** args)
{
  if (nargs < 2) 
  {
    std::cout << "Please specify at least one skim file list" << std::endl; 
    return 1; 
  }
  
  TApplication * app = new TApplication("test",0,0); 
  
  int num_datasets = nargs-1; 
  
  DECL_PLOTS
 
  for (int f = 0; f < num_datasets; f++)
  {
    string skimfile; 
    char* infile = args[f+1]; 
    ifstream input_stream(infile);

    INIT_PLOTS 

    if(strcmp(infile,"pass")!=0)
    {
        while (!input_stream.eof())
        {
          std::getline(input_stream,skimfile);  
          DmtpcSkimDataset d; 
          d.openRootFile(skimfile.c_str());  
          d.loadDmtpcEvent(false);
          d.loadClusters(false); 
          for (int i = 0; i < d.nevents(); i++)
          {
            d.getEvent(i); 
            for (int c = 0; c < d.event()->ncamera(); c++)
            {
              for (int t = 0; t < d.event()->ntracks(c); t++)
              {
                 DmtpcSkimEvent * ev = d.event(); 
                 FILL_PLOTS 
              }
            }
          }
        }
    }
    else
    {
      std::getline(input_stream,skimfile);
      TFile passfile(skimfile.c_str());  
      TTree * pass = (TTree*) passfile.Get("pass"); 
      DmtpcSkimEvent * ev = new DmtpcSkimEvent; 
      int c, t; 
      pass->SetBranchAddress("event", &ev); 
      pass->SetBranchAddress("cam",&c);
      pass->SetBranchAddress("track",&t);
      for (int i = 0; i < pass->GetEntries(); i++)
      {
        pass->GetEntry(i); 
        FILL_PLOTS 
      }
    }
  }

  AnalysisStyle::setStyle();  
  MAKE_CANVAS 
  app->Run(); 
}
