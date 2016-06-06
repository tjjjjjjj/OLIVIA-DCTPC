#include "../../../MaxCam/DmtpcSkimDataset.hh"
#include "../../../MaxCam/DmtpcDataset.hh"
#include "../../../MaxCam/DmtpcSkimEvent.hh"
#include "../../../MaxCam/DmtpcEvent.hh"
#include "../../../MaxCam/DmtpcEventTable.hh"
#include "../../../MaxCam/MaxCamClusterImage.hh"
#include "TString.h"
#include <iostream>
#include <vector>
#include "TFile.h"
#include "TH2F.h"
#include "TTree.h"
#include "TROOT.h"
#include "recoilEnergy.hh"
#include "TCanvas.h"
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include "analysis_constants.hh"
#include "AnalysisResult.hh"
#include "TObjArray.h"
#include "TVector2.h"
#include "TGraph.h"
using namespace std; 


static DmtpcEventTable * table, * table2, * livetimetable; 



char *col[6] = 
{ 
  "Cam0",
  "Cam1",
  "Cam0 E0-200",
  "Cam1 E0-200",
  "Cam0 E80-200",
  "Cam1 E80-200"
  
}; 

int mincam = 0; 

static char * colname(int c)
{
  return col[c];  
}

TH2F * crms_hist[2];
TH2F * max_hist[2];

double * rms_contours_x[2];
double * rms_contours_y[2];
double * max_contours_x[2];
double * max_contours_y[2];
int rms_contours_n[2]; 
int max_contours_n[2]; 

void loadWormCutHists()
{
  TFile *  w = new TFile("wormCutHists.root");
  
  crms_hist[0] = (TH2F*) w->Get("cam0rms");
  crms_hist[1] = (TH2F*) w->Get("cam1rms");
  max_hist[0] = (TH2F*) w->Get("cam0max");
  max_hist[1] = (TH2F*) w->Get("cam1max");

  TGraph * rms0 = (TGraph*) w->Get("rms0cntr"); 
  rms_contours_n[0] = rms0->GetN(); 
  rms_contours_x[0] = (double*) malloc(sizeof(double) * rms0->GetN());
  rms_contours_y[0] = (double*) malloc(sizeof(double) * rms0->GetN());
  for (int i = 0; i < rms0->GetN(); i++) rms0->GetPoint(i, rms_contours_x[0][i], rms_contours_y[0][i]); 

  TGraph * rms1 = (TGraph*) w->Get("rms1cntr"); 
  rms_contours_n[1] = rms1->GetN(); 
  rms_contours_x[1] = (double*) malloc(sizeof(double) * rms1->GetN());
  rms_contours_y[1] = (double*) malloc(sizeof(double) * rms1->GetN());
  for (int i = 0; i < rms1->GetN(); i++) rms1->GetPoint(i, rms_contours_x[1][i], rms_contours_y[1][i]); 

  TGraph * max0 = (TGraph*) w->Get("max0cntr"); 
  max_contours_n[0] = max0->GetN(); 
  max_contours_x[0] = (double*) malloc(sizeof(double) * max0->GetN());
  max_contours_y[0] = (double*) malloc(sizeof(double) * max0->GetN());
  for (int i = 0; i < max0->GetN(); i++) max0->GetPoint(i, max_contours_x[0][i], max_contours_y[0][i]); 

  TGraph * max1 = (TGraph*) w->Get("max1cntr"); 
  max_contours_n[1] = max1->GetN(); 
  max_contours_x[1] = (double*) malloc(sizeof(double) * max1->GetN());
  max_contours_y[1] = (double*) malloc(sizeof(double) * max1->GetN());
  for (int i = 0; i < max1->GetN(); i++) max1->GetPoint(i, max_contours_x[1][i], max_contours_y[1][i]); 

}


bool passContourWormCut(DmtpcSkimEvent *ev, int c, int t, double E_r)
{
  double rms = ev->cluster_rms(c,t); 
  double max = ev->maxpixel(c,t); 

  bool rms_inside = TMath::IsInside(E_r,rms, rms_contours_n[c+mincam], rms_contours_x[c+mincam], rms_contours_y[c+mincam]); 
  bool max_inside = TMath::IsInside(E_r,max, max_contours_n[c+mincam], max_contours_x[c+mincam], max_contours_y[c+mincam]); 

  return rms_inside && max_inside; 
}


bool passNewWormCut(DmtpcSkimEvent *ev, int c, int t, double E_r)
{
  double rms = ev->cluster_rms(c,t); 
  double max = ev->maxpixel(c,t); 

  
  /* Check out of range */
  if (
      E_r < 0 || rms < 0 || max < 0 || 
      E_r > crms_hist[c+mincam]->GetXaxis()->GetXmax() || 
      E_r > max_hist[c+mincam]->GetXaxis()->GetXmax() || 
      rms > crms_hist[c+mincam]->GetYaxis()->GetXmax() || 
      max > max_hist[c+mincam]->GetYaxis()->GetXmax())
  {
      return false; 
  }


  /* Get contour level */ 
  double rms_pct = crms_hist[c+mincam]->GetBinContent(crms_hist[c+mincam]->FindBin(E_r,rms)); 
  double max_pct = max_hist[c+mincam]->GetBinContent(max_hist[c+mincam]->FindBin(E_r,max));  
//  cout << rms_pct << " , " << max_pct << endl;
  return rms_pct > 0 && rms_pct < 0.68 && max_pct > 0 && max_pct < 0.68; 

}

bool passMyWormCuts(DmtpcSkimEvent *ev, int c, int t,  double E_v)
{
  if (c+mincam == 0)
  {
    return (ev->cluster_rms(c,t) > 2.64+0.046*E_v 
        && ev->maxpixel(c,t) < 55+0.287*E_v
        && ev->cluster_rms(c,t) < 10.51 + 0.05*E_v);
  }

  if (c + mincam ==1)
  {
    return (ev->cluster_rms(c,t) > 4.45 + 0.068 * E_v 
            && ev->cluster_rms(c,t) < 22.35 + 0.038*E_v
            &&ev->maxpixel(c,t) < 80+0.35*E_v);

  }
}

bool passWormCuts(DmtpcSkimEvent * ev, int c, int t, double E_r = -1)
{
  if (E_r < 0) E_r = RecoilEnergy::getRecoilEnergy(ev->E(c,t),c+mincam); 
  if (c+mincam==0)
  { 
    return (ev->cluster_rms(c,t) < 11.8 + 0.122*E_r && 
        ev->cluster_rms(c,t) > 5.5 + 0.07*E_r &&
        ev->maxpixel(c,t) < 55+0.287*E_r);
  }
  
  else if (c+mincam==1)
  {

    return (ev->cluster_rms(c,t) < 25 + 0.07*E_r && 
        ev->cluster_rms(c,t) > 10 + 0.05*E_r &&
        ev->maxpixel(c,t) < 80+0.35*E_r);
  }
  
}


bool passRegionCut(DmtpcSkimEvent * ev, int c, int t)
{
  return 
       ev->x(c,t)>40&&ev->x(c,t)<984&&
       ev->y(c,t)>40&&ev->y(c,t)<984; 
} 

bool passAlphaCut(DmtpcSkimEvent * ev, int c, int t)
{
  return ev->edge(c,t) == 0 && ev->range(c,t) * rangecal[c+mincam]< 5; 
}

bool passNonZeroCut(DmtpcSkimEvent * ev, int c, int t, double E_r)
{
  return E_r > 0 && ev->range(c,t) > 0 && ev->cluster_rms(c,t) > 0 && ev->neighbors(c,t) >= 0 && 
        ev->maxpixel(c,t) > 0 && ev->cluster_mean(c,t) > 0 && ev->npixel(c,t) > 0; 
}

bool passCuts(DmtpcSkimEvent * ev, int c, int t, double E_r){
  double range = ev->range(c,t) * rangecal[c+mincam];
  if (!ev->spark(c) && !ev->edge(c,t)) table->increment("!edge",colname(c)); 
  if (!ev->spark(c) && !ev->edge(c,t) && passRegionCut(ev,c,t)) table->increment("region",colname(c));  
  if (!ev->spark(c) && !ev->edge(c,t) && passRegionCut(ev,c,t) && range < 5 && range > 0) table->increment("range",colname(c));  
//  if (!ev->spark(c) && !ev->edge(c,t) && passRegionCut(ev,c,t) && range <  5 && range> 0 && E_r < 200 && E_r> 0 ) table->increment("E < 200",colname(c)); 

  return(
       E_r>0&&//E_r<200&&
       range>0&&range<5&&
       ev->x(c,t)>40&&ev->x(c,t)<984&&
       ev->y(c,t)>40&&ev->y(c,t)<984&&
       ev->edge(c,t)==0&&ev->spark(c)==0&&
       ev->cluster_rms(c,t) > 0&& 
       ev->neighbors(c,t) >= 0&&
       ev->maxpixel(c,t)>0&&
       ev->cluster_mean(c,t)>0&&
       ev->npixel(c,t)>0
       );

}


static TFile * sparkreffile; 
static TTree * sparkreftree; 
static TObjArray * sparkrefarr; 
static int sparkrefevent; 
static int sparkrefcam; 

typedef enum 
{
SR_READ,
SR_WRITE,
SR_NONE
}sparkrefstatus_t; 

sparkrefstatus_t init_spark_ref(char * file)
{
  if (file ==0) return SR_NONE; 
  
  sparkreffile = new TFile(file,"READ"); 
  sparkrefarr = 0;
  if(sparkreffile->IsOpen())
  {

    cout << "Reading file " << file << endl; 
    sparkreftree = (TTree *) sparkreffile->Get("sparkref"); 
    if (sparkreftree == NULL)
    {
      std::cerr << "No sparkref tree in file " << file << std::endl; 
      sparkreffile->Close(); 
      return SR_NONE; 
    }
    
    sparkreftree->SetBranchAddress("pixels", &sparkrefarr); 
    sparkreftree->SetBranchAddress("event", &sparkrefevent); 
    sparkreftree->SetBranchAddress("cam", &sparkrefcam); 
    return SR_READ; 
  }

  else
  {
    cout << "Creating file " << file << endl; 
    sparkreffile = new TFile(file,"CREATE"); 
    sparkreftree = new TTree("sparkref", "sparkref");
    sparkreftree->Branch("pixels", &sparkrefarr); 
    sparkreftree->Branch("event", &sparkrefevent,"event/I"); 
    sparkreftree->Branch("cam", &sparkrefcam,"cam/I"); 
    
    return SR_WRITE;
  } 
} 

void write_spark_ref(vector<int *> * sathist, int event, int cam)
{
  sparkreffile->cd(); 
  sparkrefarr = new TObjArray(sathist->size());  

  for (int i = 0; i < sathist->size(); i++)
  {
    int * coords = sathist->at(i); 

    TVector2 * vec = new TVector2(coords[0],coords[1]); 
    //cout << "Writing " << vec->X() << " , " << vec->Y() << endl; 
    sparkrefarr->Add(vec); 
  }

  sparkrefevent = event; 
  sparkrefcam = cam; 
  sparkreftree->Fill(); 
}

int sparkrefindex = 0; 
void load_spark_ref(vector<int *> * sathist, int i,int  c)
{
  sparkreftree->GetEntry(sparkrefindex);

  if (sparkrefevent!=i || sparkrefcam !=c) return;

  sathist->clear(); 


  for (int i = 0; i < sparkrefarr->GetEntries(); i++)
  {
    int * p = (int*) malloc(sizeof(int) *2); 
    TVector2 * vec = (TVector2*) sparkrefarr->At(i); 
    p[0] = (int) vec->X(); 
    p[1] = (int) vec->Y(); 
   // std::cout << "adding pixel " << p[0]<< "," <<p[1]  << std::endl;
    sathist->push_back(p); 
  }

  sparkrefindex+=1; 
}

void spark_ref_save()
{
  sparkreffile->cd();
  sparkreftree->Write(); 
  sparkreffile->Close(); 
}

static int satThresh = 55000;
static int nsatThresh = 4; 
static int binThresh = 8; 

static int populateSatHist(TH2F * image, TH2F * bias, vector<int *> * sathist)
{ 
  image->Add(bias,-1);
  TH2F * s = new TH2F("sh","sh",256,0,1024,256,0,1024);
  int nsat = 0;

  int ret = 0;
  for (int m = 1; m <257; m++)
  {
    for (int n = 1; n <257; n++)
    {
      if(image->GetBinContent(m,n) > satThresh)
      {
        nsat++;   
        s->SetBinContent(m,n,1);
      }
      else
      {
        s->SetBinContent(m,n,0);  
      }
    }
  } 
  //cout << endl <<"nsat: " << nsat << endl;  
   //If enough saturated pixels, add these pixels to the bad ones
  if (nsat > nsatThresh)
  {
    for (int m = 1; m < 257; m++)
    {
      for (int n = 1; n< 257; n++)
      {
        if (s->GetBinContent(m,n) > 0)
        {
           //std::cout << "adding pixel " << m<< "," << n << std::endl;
           ret++; 
           int * px = (int *) malloc(2*sizeof(int)); 
           px[0] = m;
           px[1] = n; 
           sathist->push_back(px);  
        }
      }
    }

  }
 
  s->Delete();
  return ret; 

}


bool passSatHistCut(double x, double y, vector<int*> * sathist)
{
  int x0 = TMath::Nint(x/4); 
  int y0 = TMath::Nint(y/4); 

  if (x0 <= 0 || y0 <=0) return false; 

  for (int i = 0; i < sathist->size(); i++)
  {
    int * px = (*sathist)[i]; 
    if (abs(x0-px[0])<=binThresh && abs(y0-px[1])<=binThresh) return false; 
  } 
  return true; 
}

static char * passfilename = 0; 
static char * ratesfilename = 0; 
static char * sparkreffilename = 0; 
int main(int argn, char ** args)
{
  DmtpcSkimDataset ds;
  ds.openRootFile(args[1]);
  int run = atoi(strstr(args[1],"0"));  
  if (argn>2) passfilename = args[2]; 
  if (argn>3) ratesfilename = args[3]; 
  if (argn>4) sparkreffilename = args[4]; 
  if (argn>5) mincam = atoi(args[5]); 
 

  loadWormCutHists(); 
   

  TFile * ratesfile; 

  /* Set up output tables */
  table = new DmtpcEventTable("mytab1","mytab1"); 
  table2 = new DmtpcEventTable("tbl1","tbl1"); 
  livetimetable = new DmtpcEventTable("livetime","livetime"); 


  /* Check if sparref file exists */
  sparkrefstatus_t srstat = init_spark_ref(sparkreffilename); 
  if (srstat!=SR_READ)
  {
    ds.loadDmtpcEvent(true);
    ds.loadBiasFrames(true);  
  }
    
  ds.loadClusters(false); 
  ds.getEvent(0); 
  int ncam = ds.event()->ncamera(); 
  //This array keeps track of which tracks pass the pre_rbi cut 
  bool pass_pre_rbi_cut[ds.nevents()][ncam][15];
  
  //This array keeps track of the cumulative number of tracks passing the pre_rbi cut
  int cum_pass_pre_rbi_cut[ncam][ds.nevents()][15];  
  cout << sizeof(cum_pass_pre_rbi_cut) << endl; 
  UInt_t pass_info[ds.nevents()][ncam][15]; 

/* Initialize all prerbi cuts to not passing */
  for (int i = 0; i < ds.nevents(); i++)
  { 
    for (int j = 0; j < ncam; j++) 
    { 
      for(int k = 0; k < 15; k++)
      {
        pass_pre_rbi_cut[i][j][k]=false;
        cum_pass_pre_rbi_cut[j][i][k] = 0; 
        pass_info[i][j][k] = 0; 
      }
    }
  }

  TH2F * bias, * image; 
  vector<int *> * sathist[2];
  
  bias = NULL;
  image = NULL;
  sathist[0] = new vector<int*>;
  sathist[1] = new vector<int*>;   
  
  /* Go through and make pre_rbi_cuts */
  for (int i = 0; i < ds.nevents(); i++)
  {
   ds.getEvent(i); 
   for (int c = 0; c < ncam; c++)
   {
       
     //Do sparkref thing...
     if (ds.event()->spark(c))
     {
        //Account for double counting dead pixels from previous frames
        livetimetable->increment("spark",colname(c)); 
        if (srstat == SR_READ)
        {
          load_spark_ref(sathist[c],i,c); 
        }

        else 
        {
          bias = ds.getBiasFrame(c); 
          image = ds.orig_event()->ccdData(c); 
          int newpx = populateSatHist(image,(TH2F*)bias->Clone(),sathist[c]); 
          if (srstat == SR_WRITE && newpx > 0) write_spark_ref(sathist[c],i,c); 
        }
     } 
     else 
     {
       //table->increment("!spark",colname(c)); 
     }

     //if (ds.event()->ntracks(c) > 0) table->increment("clstrfnd",colname(c)); 
     table->add("tracks",colname(c),ds.event()->ntracks(c)); 

     for (int t = 0; t < 15; t++)
     {    
       
       //If no more tracks,set remaining tracks in cumulative pass array to last value
       //and go on
       if (t == ds.event()->ntracks(c))
       {
         for (int tt = t; tt< 15; tt++)
         {
            if (i==0 && t==0) break; 
            cum_pass_pre_rbi_cut[c][i][tt] = *(&(cum_pass_pre_rbi_cut[c][i][t]) - 1); 
         }
         break; 
       }
      
       pass_info[i][c][t] |= CONSIDERED; 

       double E_r = RecoilEnergy::getRecoilEnergy(ds.event()->E(c,t),c+mincam);

       bool passBurnin = true; 
       bool passWorm = true; 

       if (sathist[c]!=NULL)
       {
         passBurnin = passSatHistCut(ds.event()->x(c,t),
                                     ds.event()->y(c,t),
                                     sathist[c]); 
       }



       if (passBurnin)
       {
         pass_info[i][c][t] |= PASS_RBI; 
         table->increment("sparkref",colname(c));
         pass_pre_rbi_cut[i][c][t] = true; 
         if (i==0 && t==0) cum_pass_pre_rbi_cut[c][0][0] = 1; //Special case if first track passes
         else cum_pass_pre_rbi_cut[c][i][t] = *(&(cum_pass_pre_rbi_cut[c][i][t])- 1) + 1;  // Add one to the value of the last 
       }

       else
       {
         livetimetable->add("rbi_killedpix",colname(c), ds.event()->npixel(c,t));
          if(i!=0 || t!=0) cum_pass_pre_rbi_cut[c][i][t] = *(&(cum_pass_pre_rbi_cut[c][i][t]) - 1);  // same as last 
       }

       
       if (passNewWormCut(ds.event(),c,t,E_r)) pass_info[i][c][t] |= PASS_WORM;
       if (passAlphaCut(ds.event(),c,t)) pass_info[i][c][t] |= PASS_ALPHA;
       if (passRegionCut(ds.event(),c,t)) pass_info[i][c][t] |= PASS_REGION;
       if (E_r > 0 && E_r < 200) pass_info[i][c][t] |= E0_200;
       if (E_r > 80 && E_r < 200) pass_info[i][c][t] |= E90_200;

     }
   }
  }
  ds.loadDmtpcEvent(false); 
  vector<int * > pass_events;   //camera, event, track
  //cout << "done with prerbi\n";

  ///rbi cut
  for (int i = 0; i < ds.nevents(); i++)
  {
    ds.getEvent(i); 
    for (int c =0; c < ncam; c++)
    {
      for (int t=0; t<15; t++)
      {
        if ( t == ds.event()->ntracks(c)) break; 
        if (!pass_pre_rbi_cut[i][c][t]) continue;         
        int * pass_event = new int[3];
        
        pass_event[0] = c;
        pass_event[1] = i;
        pass_event[2] = t; 
        
        if (ds.event()->nburnin(c,t) == 0) 
        {
          pass_events.push_back(pass_event);
          continue;   
        }
        
        
        //Loop through last and next 20 PASSING events to see 

        int this_idx = ds.event()->burnin_this_index(c,t); 
        int howmany = 20; 
        bool pass_pos = true; 
      //  cout << endl << "====" << ds.event()->eventNumber()<< "===="  << endl; 
       // cout << "this_idx: " << this_idx << endl; 
        //First go back
        for (int b = this_idx-1; b>=0; b--)
        {
          //condition for 20 passing events before burnin
          
          int diff =  (cum_pass_pre_rbi_cut[c][ds.event()->burnin_event(c,t,b)][ds.event()->burnin_track(c,t,b)] - cum_pass_pre_rbi_cut[c][i][t]);
          if (diff < -howmany) break; 
          //If any of these passed, then fail this event. 
          if (pass_pre_rbi_cut[ds.event()->burnin_event(c,t,b)][c][ds.event()->burnin_track(c,t,b)])
          {
           //cout << "FAIL" << endl; 
           pass_pos = false;
           pass_info[i][c][t] &= ~PASS_RBI; 
           livetimetable->add("rbi_killedpix",colname(c), ds.event()->npixel(c,t));
           break; 
          }
        }
        
        if (!pass_pos)
        {
            delete pass_event; 
            continue; 
        }
 
        //Now, go forth
        for (int b = this_idx; b<ds.event()->nburnin(c,t); b++)
        {
          if (b >= ds.event()->nburnin(c,t)) break; 
          //condition for 20 passing events after burnin
          //If any of these passed, then fail this event. 
          int diff =  (cum_pass_pre_rbi_cut[c][ds.event()->burnin_event(c,t,b)][ds.event()->burnin_track(c,t,b)] - cum_pass_pre_rbi_cut[c][i][t]);
          if (diff > howmany) break; 
          if (pass_pre_rbi_cut[ds.event()->burnin_event(c,t,b)][c][ds.event()->burnin_track(c,t,b)])
          {
           pass_pos = false;
           pass_info[i][c][t] &= ~PASS_RBI; 
           livetimetable->add("rbi_killedpix",colname(c), ds.event()->npixel(c,t));
           //cout << "FAIL" << endl; 
           break; 
          }
        }
 
        if (pass_pos)
          pass_events.push_back(pass_event);
        else
          delete pass_event; 
      }
    }
  }
  
  TH2F * rve_pass = new TH2F ("rve_pass", "rve_pass", 100,0,800,40,0,20); 

  vector<int*> final_pass;  
  for (int p = 0; p < pass_events.size(); p++)
  {
    int c = pass_events[p][0];
    int i = pass_events[p][1];
    int t = pass_events[p][2];
     
    table->increment("poscut",colname(c));  
    ds.getEvent(i); 
    //cout << c << " , " << i << " , " << t<< "           " << ds.event()->eventNumber()<<endl; 
     double E = ds.event()->E(c,t); 

     double Er = RecoilEnergy::getRecoilEnergy(E,c+mincam); 
  //   double Er = RecoilEnergy::getRecoilEnergy(E,c+mincam,ds.event()->x(c,t), ds.event()->y(c,t)); 
     double E_v = mincam + c == 0 ? E / 9.5 : E/12.9; 
  //  if (passMyWormCuts(ds.event(), c, t, E_v))
     if (passNewWormCut(ds.event(),c,t,Er))
     {
       table->increment("worm",colname(c));  
       if (passCuts(ds.event(),c,t,Er))
       {
         //cout << " CAMERA: " << c << " EVENT: "<<ds.event()->eventNumber() <<" TRACK: " << t << " PASSED.    E_R:" << Er << endl; 
         table2->increment("candidate",colname(2+c)); 
         double Er_star = RecoilEnergy::getRecoilEnergy(E,c+mincam,ds.event()->x(c,t), ds.event()->y(c,t)); 
         if (Er > 80 && Er < 200) 
         { 
           table->increment(" E 80-200",colname(c)); 
           table2->increment("candidate",colname(4+c)); 
           rve_pass->Fill(RecoilEnergy::getVisibleEnergy(ds.event(),c+mincam,t), 
           ds.event()->range(c,t)*rangecal[c]); 
         } 

         if (Er_star > 80 && Er_star < 200)
         {
           table->increment(" E* 80-200", colname(c)); 
         }

         int * fp = new int[3];
         fp[0] = c; 
         fp[1] = i; 
         fp[2] = t; 
         final_pass.push_back(fp); 
       }
     }
  } 

  TH2F * rve_rbi = new TH2F ("rve_rbi", "rve_rbi", 50,0,800,40,0,20); 
  TH2F * rve_alpha = new TH2F ("rve_alpha", "rve_alpha", 50,0,6000,40,0,100); 
  TH2F * rve_alpha2 = new TH2F ("rve_alpha2", "rve_alpha2", 50,0,800,40,0,20); 
  TH2F * rve_worm = new TH2F ("rve_worm", "rve_worm", 50,0,800,40,0,20); 
  //Figure out exclusive things
  for (int i = 0; i < ds.nevents(); i++)
  {
    ds.getEvent(i); 
    for (int c = 0; c < ncam ; c++)
    {
      for (int t = 0; t < 15; t++)
      {
        UInt_t e = pass_info[i][c][t]; 

        if (e & CONSIDERED == 0) continue; 
        if (e & PASS_REGION == 0) continue;  
        

        double range = ds.event()->range(c,t) * rangecal[c]; 
        double E_v = RecoilEnergy::getVisibleEnergy(ds.event(),c,t); 
        //0-200 
        if (e & E0_200)
        {
          table2->increment("tracks",colname(2+c)); 
        } 
        if (e & E90_200)
        {
          table2->increment("tracks",colname(4+c)); 
        }
        //alpha 
        if ( (e & (PASS_RBI )) && !(e & PASS_ALPHA)) 
        {
          rve_alpha->Fill(E_v, range); 
          if (e & E0_200) table2->increment("alpha",colname(2+c));
          if (e & E90_200) table2->increment("alpha",colname(4+c));

          if (range > 1 && E_v > 50)
          {

            rve_alpha2->Fill(E_v, range); 
            if (e & E0_200) table2->increment("alpha2",colname(2+c));
            if (e & E90_200) table2->increment("alpha2",colname(4+c));

          }
           
        }
        //worm
        else if ( (e & (PASS_RBI )) && !(e & PASS_WORM)) 
        {
          rve_worm->Fill(E_v, range); 
          if (e & E0_200) table2->increment("worm",colname(2+c));
          if (e & E90_200) table2->increment("worm",colname(4+c));
        }
        //rbi
        else if ( !(e & PASS_RBI)) 
        {
          rve_rbi->Fill(E_v, range); 
          if (e & E0_200) table2->increment("rbi",colname(2+c));
          if (e & E90_200) table2->increment("rbi",colname(4+c));
        }

      }
    }
  }

  if (passfilename)
  {
    ds.loadDmtpcEvent(true); 
    DmtpcSkimEvent * pass_evt = ds.event(); 
    DmtpcEvent * orig_evt = ds.orig_event(); 
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
    passTree->Branch("orig_event", "DmtpcEvent", &orig_evt, 128000,1); 

    TH2F * rve_pass_c = (TH2F*) rve_pass->Clone(); 
    TH2F * rve_rbi_c = (TH2F*) rve_rbi->Clone(); 
    TH2F * rve_worm_c = (TH2F*) rve_worm->Clone(); 
    TH2F * rve_alpha_c = (TH2F*) rve_alpha->Clone(); 
    TH2F * rve_alpha2_c = (TH2F*) rve_alpha2->Clone(); 
    
    AnalysisResult * result = new AnalysisResult(ncam, ds.nevents(),(UInt_t*)pass_info); 
    cout << result->Write("result") << endl; 

    rve_pass_c->Write();
    rve_rbi_c->Write();
    rve_worm_c->Write();
    rve_alpha_c->Write();
    rve_alpha2_c->Write();

    for (int p = 0; p < final_pass.size(); p++)
    {
      pass_cam = final_pass[p][0];
      pass_track = final_pass[p][2];
      ds.getEvent(final_pass[p][1]);
      pass_evt = ds.event(); 
      orig_evt = ds.orig_event(); 
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
    DmtpcEventTable * writeme2 = new DmtpcEventTable("tab2","tab2");
    writeme2->combineWith(table2); 
    writeme2->Write();
    DmtpcEventTable * writeme3 = new DmtpcEventTable("tab3","tab3");
    writeme3->combineWith(livetimetable); 
    writeme3->Write();
    tablefile->Write();
    tablefile->Close();   
  }

  if (srstat==SR_WRITE)
  {
   spark_ref_save();  
  }
}

 
