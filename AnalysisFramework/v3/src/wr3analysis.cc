#include "../../../MaxCam/DmtpcSkimDataset.hh"
#include "../../../MaxCam/DmtpcSkimEvent.hh"
#include <iostream>
#include "../../../MaxCam/DmtpcEventTable.hh"
#include "analysis_constants.hh"
#include "AnalysisResult.hh"
#include "recoilEnergy.hh"
#include "MultiVariate.hh"
#include <string.h>


using namespace std; 
using namespace MultiVariate;

static char * passfilename = 0; 
static char * ratesfilename = 0; 
MultiVariateResult* mv = new MultiVariateResult();


char *col[6] = 
{ 
  "Cam0",
  "Cam1",
  "Cam0 E0-200",
  "Cam1 E0-200",
  "Cam0 E80-200",
  "Cam1 E80-200"
  
}; 

char * colname(int c)
{
  return col[c];   
}

double gain[2] = { 7.6,7.0}; 
double rangecalloc[2] =  {0.171,0.137};

double getEv(DmtpcSkimEvent * e, int c,int t)
{
  return e->E(c,t)/gain[c]; 
}

double getEr(DmtpcSkimEvent *e,int c, int t)
{
  return RecoilEnergy::visibleToRecoil(getEv(e,c,t)); 
}


bool passRegionCut(DmtpcSkimEvent * ev, int c, int t)
{
  return 
       ev->x(c,t)>40&&ev->x(c,t)<984&&
       ev->y(c,t)>40&&ev->y(c,t)<984; 
} 

bool passAlphaCut(DmtpcSkimEvent * ev, int c, int t)
{
  return ev->edge(c,t) == 0 && ev->range(c,t) * rangecalloc[c]< 5; 
}

bool passNonZeroCut(DmtpcSkimEvent * ev, int c, int t, double E_r)
{
  return E_r > 0 && ev->range(c,t) > 0 && ev->cluster_rms(c,t) > 0 && ev->neighbors(c,t) >= 0 && 
        ev->maxpixel(c,t) > 0 && ev->cluster_mean(c,t) > 0 && ev->npixel(c,t) > 0; 
}

bool passWormCut(DmtpcSkimEvent* ev, int c, int t)
{								
   double classifier = mv->getClassifier(ev,c,t);
//   return (ev->maxpixel(c,t)<250 && ev->cluster_rms(c,t)<30);
   return classifier>15;
}

bool passCharge(DmtpcSkimEvent* ev, int c)
{
   return ev->ntriggers()>0 && ev->ntriggers()<=10;
}

int main (int nargs, char** args)
{

  mv->setResult("wr3_worm"); 

  DmtpcSkimDataset ds; 
  ds.openRootFile(args[1]);
  int run = atoi(strstr(args[1],"0"));  
  if (nargs>2) passfilename = args[2];
  if (nargs>3) ratesfilename = args[3];
  TFile * ratesfile; 
  /* Set up output tables */
  DmtpcEventTable * table = new DmtpcEventTable("table","table"); 
  DmtpcEventTable * livetimetable = new DmtpcEventTable("livetime","livetime"); 
  ds.loadClusters(false); 
  ds.getEvent(0); 
  int ncam = ds.event()->ncamera(); 


 vector<int * > passing;   //camera, event, track
 UInt_t pass_info[ds.nevents()][ncam][15]; 

 memset(pass_info,0,sizeof(pass_info)); 
 
  for (int i = 0; i < ds.nevents(); i++)
  {
    ds.getEvent(i); 
    for (int c = 0; c < ncam; c++)
    {
      if (ds.event()->spark(c))
      {
        livetimetable->increment("spark",colname(c)); 
        continue;  
      }
      
      table->add("tracks",colname(c),ds.event()->ntracks(c)); 

      for (int t = 0; t < ds.event()->ntracks(c); t++)
      {
        bool pass= true; 
        DmtpcSkimEvent * e = ds.event(); 
        double E_r = getEr(e,c,t); 
        if (e->spark(c)) continue; 

        if (!passNonZeroCut(e,c,t,E_r)) pass = false; 
        else 
        {
          table->increment("valid_params",colname(c)); 
          pass_info[i][c][t] |= CONSIDERED; 
        }

        if(e->isRBI(c,t)) pass = false; 
        else
        { 
          table->increment("rbi",colname(c)); 
          pass_info[i][c][t] |= PASS_RBI; 
        }

        if (e->edge(c,t)) pass = false;
        else table->increment("!edge",colname(c)); 

        if (!passRegionCut(e,c,t)) pass = false; 
        else 
        {
          table->increment("region",colname(c)); 
          pass_info[i][c][t] |= PASS_REGION ; 
        }

        double range = e->range(c,t) * rangecalloc[c]; 
        if (range >= 5 || range <= 0 ) pass = false; 
        else 
        {
          table->increment("range",colname(c)); 
          if (!e->edge(c,t))
          {
            pass_info[i][c][t] |= PASS_ALPHA; 
          }
        }
	if(!passWormCut(e,c,t)) pass = false;
	else
	{
	   table->increment("worm",colname(c));
	   pass_info[i][c][t] |= PASS_WORM;
	}
	// if(!passCharge(e,c)) pass = false;
// 	else
// 	{
// 	   table->increment("charge",colname(c));
// 	   pass_info[i][c][t] |= PASS_CHARGE;
// 	}

        if (E_r > 0 && E_r < 200)
        {
          pass_info[i][c][t] |= E0_200; 
          if (E_r> 90)
          {
            pass_info[i][c][t] |= E90_200; 
          }
        }

        if (pass)
        {
          int *fp = new int[3]; 
          fp[0] = c; 
          fp[1] = i; 
          fp[2] = t; 
          passing.push_back(fp); 
        }
      }
    }

  }
  if (passfilename)
  {
    DmtpcSkimEvent * pass_evt = ds.event(); 
    int pass_cam; 
    int pass_track;
    int pass_run =run;
     //Check if the pass file already exists. If it doesn't, create it
    TFile * passfile = TFile::Open(passfilename,"RECREATE"); 
    passfile->cd(); 

    TTree * passTree = new TTree("skim","pass DMTPC events"); 
    passTree->Branch("event", "DmtpcSkimEvent", &pass_evt,128000,1); 
    passTree->Branch("cam", &pass_cam, "pass_cam/I"); 
    passTree->Branch("track", &pass_track, "pass_track/I"); 
    passTree->Branch("run", &pass_run, "pass_run/I"); 

    AnalysisResult * result = new AnalysisResult(ncam, ds.nevents(),(UInt_t*)pass_info); 
    cout << result->Write("result") << endl; 

    for (int p = 0; p < passing.size(); p++)
    {
      pass_cam = passing[p][0];
      pass_track = passing[p][2];
      ds.getEvent(passing[p][1]);
      pass_evt = ds.event(); 
     // cout << pass_evt->E(pass_cam,pass_track) << endl; 
      passTree->Fill(); 
      passTree->Write(); 
    }  
    passfile->Close();  
  }

  table->print(); 

  // If a ratesfilename is specified, we want to write out event rates to that file
  if (ratesfilename)
  {
    //If a lockfile is specified, then we need to obtain the lock before writing 
    
    TFile * tablefile = new TFile(ratesfilename,"RECREATE"); 
    tablefile->cd(); 
    DmtpcEventTable * writeme = new DmtpcEventTable("tab1","tab1");
    writeme->combineWith(table); 
    writeme->Write();
    DmtpcEventTable * writeme3 = new DmtpcEventTable("tab2","tab2");
    writeme3->combineWith(livetimetable); 
    writeme3->Write();
    tablefile->Write();
    tablefile->Close();   
  }


}
