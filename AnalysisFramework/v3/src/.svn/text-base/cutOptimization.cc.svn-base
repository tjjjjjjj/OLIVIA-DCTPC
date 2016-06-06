#include "TApplication.h"
#include "contourPicker.hh"
#include "AnalysisResult.hh"
#include "../../../MaxCam/DmtpcSkimEvent.hh"
#include "TChain.h"
#include "recoilEnergy.hh"
#include "TCanvas.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TList.h"
#include "TGraph.h"
#include "style.hh"
#include "TH2F.h"
#include "TF1.h"
#include "TMath.h"

static const int num_sets = 5; 

static const int contour_bins = 50; 



char * data_names[] = 
{
 "Surface Wimp Run Data",                //0
 "Surface CF252 Data",                   //1
 "Surface Worm Run",                     //2
 "CF252 MC",                             //3
 "Wimp MC",                              //4
 "RBIs from surface wimp run",           //5
 "worms from surface wimp run (npixel)", //6
 "worms from surface wimp run (cuts)"    //7
};

int width = 1200;
int height = 600;
int ncanv = 8; 
char * canvas_names[] =
{
  "cmax","crms","cneigh","ckurt","crange","c_m2_v_kurt","crange_v_kurt","cnpix"
}; 

int xmax[] = {200,200,200,200,200,20,20,20};
int xmin = 0;
int xbins = 50; 

int ymax[] = {400,100,100,5,20,5,5,100};
int ymin[] = {0,0,0,-10,0,-10,-10,0};
int ybins[] = {50,50,50,100,50,100,100,50};
  


char * merged_dir = "/net/zwicky/dmtpc/analysis/v3/ntuples/";
char * pass_prefix = "/net/zwicky/dmtpc/analysis/v3/pass/pass";


bool passPartialWormCuts(DmtpcSkimEvent * ev, int c, int t, int mincam, double E_v)
{

  if (c+mincam==0)
  { 
    return (//ev->cluster_rms(c,t) < 11.8 + 0.122*E_v && 
        ev->maxpixel(c,t) < 55+0.287*E_v);
  }
  
  else if (c+mincam==1)
  {
    return (//ev->cluster_rms(c,t) < 25 + 0.07*E_v && 
        ev->maxpixel(c,t) < 80+0.35*E_v);
  }
}

bool passWormCuts(DmtpcSkimEvent * ev, int c, int t, int mincam, double E_v)
{

  if (c+mincam==0)
  { 
    return (ev->cluster_rms(c,t) < 11.8 + 0.122*E_v && 
        ev->cluster_rms(c,t) > 5.5 + 0.07*E_v &&
        ev->maxpixel(c,t) < 55+0.287*E_v);
  }
  
  else if (c+mincam==1)
  {
    return (ev->cluster_rms(c,t) < 25 + 0.07*E_v && 
        ev->cluster_rms(c,t) > 10 + 0.05*E_v &&
        ev->maxpixel(c,t) < 80+0.35*E_v);
  }
}

bool passMyWormCuts(DmtpcSkimEvent *ev, int c, int t, int mincam, double E_v)
{
  if (c+mincam == 0)
  {
    return (ev->cluster_rms(c,t) > 1.895+0.0288*E_v 
        && ev->cluster_rms(c,t) < 12.229 + 0.0523*E_v);
  }

  if (c + mincam ==1)
  {
    return (ev->cluster_rms(c,t) > 5.723 + 0.0268 * E_v 
            && ev->cluster_rms(c,t) < 24.55 + 0.05*E_v);

  }
}


char * get_full_name(char * base)
{
  char * ret = (char*)malloc(strlen(merged_dir) + strlen(base) + 1); 
  strcpy(ret,merged_dir);
  strcpy(ret+strlen(merged_dir),base);
  return ret; 
}

char * concat(char * a, int b) 
{

  TString str(a);
  str+=b;
  char * ret = (char *) malloc(strlen(str.Data())+1); 
  strcpy(ret,str.Data()); 
  return ret;
}

static int last_run;
static AnalysisResult * result; 
static TFile * pass_handle; 
/* Loads the analysis result for the given run */
static void  loadAnalysisResult(int run)
{
  if (last_run !=-1) 
  {
    result->Delete(); 
  }
 
  TString fn(pass_prefix);
  fn+=run;
  fn+=".root";  
  pass_handle = new TFile(fn.Data());
  result = (AnalysisResult*) ((AnalysisResult *) pass_handle->Get("result"))->Clone(); 
  pass_handle->Close();
  
} 


/* Returns the minimum camera for a run. This is important
   for runs with only 1 camera */
int get_min_cam(int run)
{

  if (run > 100421005 && run <= 100421010) return 1;
  if (run > 100421015 && run <= 100421020) return 1;
  if (run > 100429005 && run <= 100429010) return 1; 
  if (run > 100502000 && run <= 100502004) return 1; 
  return 0; 


}

TGraph * get68Contour(char * name)
{

  TObjArray * cntrs = (TObjArray *) gROOT->GetListOfSpecials()->FindObject("contours"); 
  TList * lvl68 = (TList*) cntrs->At(0); 
  int index = 0; 
  int max = 0; 
  //Get contour with max number of points 
  for (int i = 0; i < lvl68->GetEntries(); i++)
  {
    TGraph * trial = (TGraph*) lvl68->At(i); 
    if (trial->GetN() > max)
    {
      max = trial->GetN(); 
      index = i; 
    }
  }
  
  TGraph * cntr68 = (TGraph*) (lvl68->At(index)); 
  return (TGraph*) cntr68->Clone(name); 

}


int main (int nargs, char **args)
{
  TApplication * app = new TApplication("cut_optimization",0,0);
  AnalysisStyle::setStyle(); 
  last_run = -1; 
  TChain * data[num_sets];

  //Set up the surface run data
  data[0] = new TChain("skim"); 
  data[0]->Add(get_full_name("merged801_850.root"));
/*  data[0]->Add(get_full_name("merged851_900.root"));
  data[0]->Add(get_full_name("merged901_950.root"));
  data[0]->Add(get_full_name("merged951_1000.root"));
  data[0]->Add(get_full_name("merged1001_1031.root"));
*/
  //Use the same chain for things derived from surface run 
  data[5] = data[0];
  data[6] = data[0];
  data[7] = data[0];
  
  //Set up the surface CF252 data
  data[1] = new TChain("skim");
  //data[1] = data[0]; 
  data[1]->Add(get_full_name("merged1032_1073.root"));

  //Set up the Surface Worm Candidate Data
  data[2] = new TChain("skim");
  data[2]->Add(get_full_name("merged779_784.root"));

  //Set up the CF252 MC Data
  data[3] = new TChain("skim");
  data[3]->Add(get_full_name("merged100429002_100429010.root"));
  data[3]->Add(get_full_name("merged100502001_100502008.root"));
  //data[3]->Add(get_full_name("merged100421001_100421010.root"));

  //Set up the Wimp MC Data
  data[4] = new TChain("skim");
  data[4]->Add(get_full_name("merged100421011_100421020.root"));
  
  //Set up the canvases..
  TCanvas * cv[2][ncanv]; 
  for (int v = 0; v < ncanv; v++) 
  {
     for (int c = 0; c < 2; c++)
     {
      cv[c][v] = new TCanvas(concat(canvas_names[v],c),concat(canvas_names[v],c),width,height); 
      cv[c][v]->Divide(4,2); 
     }
  }

  TCanvas * cam1 = new TCanvas("cam1","cam1",500,500); 
  TCanvas * cam0 = new TCanvas("cam0","cam0",500,500); 

  /*Histogram for fake contour stuff */

  TH2F * h0 = new TH2F("h0","h0",contour_bins,0,400,
                                 contour_bins,0,80);
  TH2F * h1 = new TH2F("h1","h1",contour_bins,0,400,
                                 contour_bins,0,80);
  TH2F * h2 = new TH2F("h2","h2",contour_bins,0,400,
                                 contour_bins,0,180);
  TH2F * h3 = new TH2F("h3","h3",contour_bins,0,400,
                                 contour_bins,0,180);

  TH2F * h[2][ncanv];

  for (int set = 0; set < num_sets; set++)
  {
    std::cout << "Now processing set " << data_names[set] << std::endl; 

    for (int v = 0; v < ncanv; v++)
    {
      for (int c = 0; c < 2; c++)
      {
        h[c][v] = new TH2F(concat(concat(canvas_names[v],set),c),data_names[set], 
                      xbins, xmin, xmax[v], ybins[v], ymin[v], ymax[v]);  
      }
    }
    std::cout << "Created Histograms " << std::endl; 
    DmtpcSkimEvent * ev = new DmtpcSkimEvent; 
    data[set]->SetBranchAddress("event",&ev); 

    int j = 0; //This index keeps track of our index inside each run
    for (int i = 0; i < data[set]->GetEntries(); i++)
    {
      data[set]->GetEntry(i);
      int run = ev->runNumber(); 

      int mincam = get_min_cam(run); 
      //Reset j and load the different analysis result if run number has changed
      if (run!=last_run) {
        std::cout << "Run: "<< run << std::endl; 
        j = 0;
        //Don't need analysis result for MC data..
        if (set != 3 && set!=4) 
          loadAnalysisResult(run);      
        last_run = run; 
        std::cout <<" ncam: " << ev->ncamera() << std::endl; 
        std::cout <<" min_cam: " << mincam <<std::endl; 
      } 
      
      int ncam = ev->ncamera(); 
      for (int c = 0; c < ncam ; c++)
      {
        for (int t = 0; t < 15; t++)
        {
          /** CUTS for all events!! **/
          if (ev->spark(c)) break;  //This should be covered by below... but whatever
          if (ev->ntracks(c) == t) break; 
          if (ev->edge(c,t) == 1) continue; 
          double E = ev->E(c,t); 
          double E_v = c+mincam == 0 ? E / 9.5 : E/12.9; 
          double E_r = RecoilEnergy::getRecoilEnergy(E,c+mincam); 

          if (set!=3 && E_v > 200) continue; 
        //  if (!passPartialWormCuts(ev,c,t,mincam,E_v)) continue;  
          /** RBI Cut **/
          if (result->pass_rbi(j,c,t) && set == 5) continue;
          if (!result->pass_rbi(j,c,t) && (set == 0 || set == 1 || set == 2 || set ==6 || set==7)) continue;
    //      if (set !=7 && !passWormCuts(ev,c,t,mincam)) continue;  

          /** npixel cut **/ 
//          if ((set==2 || set == 6) && ev->npixel(c,t) > 64) continue;
      
          /** cluster_rms cut */  
          if ( set==7 && (passWormCuts(ev,c,t,mincam,E_v))) continue; 

          if (set == 3)
          {
              c+mincam == 0 ? h0->Fill(E_r,ev->cluster_rms(c,t)) : h1->Fill(E_r,ev->cluster_rms(c,t)) ; 
              c+mincam== 0 ? h2->Fill(E_r,ev->maxpixel(c,t)) : h3->Fill(E_r,ev->maxpixel(c,t)) ; 
          }


          double range = ev->range(c,t) * rangecal[c+mincam]; 
          double mom2= ev->moment(c,t,2);

          double kurt = TMath::Log(mom2) - TMath::Log(ev->cluster_rms(c,t)); 

          h[c+mincam][0]->Fill(E_v,ev->maxpixel(c,t)); 
          h[c+mincam][1]->Fill(E_v,ev->cluster_rms(c,t)); 
          h[c+mincam][2]->Fill(E_v,ev->neighbors(c,t)); 
          h[c+mincam][3]->Fill(E_v,kurt); 
          h[c+mincam][4]->Fill(E_v,range);
          h[c+mincam][5]->Fill(mom2,kurt); 
          h[c+mincam][6]->Fill(range,kurt); 
          h[c+mincam][7]->Fill(range,ev->npixel(c,t)); 
        }
      }
      j++;
    }  
    
    for (int v = 0; v < ncanv; v++) 
    {
      for (int c = 0; c < 2; c++)
      {
        cv[c][v]->cd(set+1); 
        h[c][v]->DrawCopy("colz"); 
        cv[c][v]->Update(); 
      }
    }
   
    
   
  }
  
 TF1 * cam0_maxpixel = new TF1("f1","55+0.287*x",0,200);  
 TF1 * cam1_maxpixel = new TF1("f2","80+0.35*x",0,200);  

 TF1 * cam0_crms1 = new TF1("f3","11.8+0.122*x",0,200);  
 TF1 * cam1_crms1 = new TF1("f4","25+0.07*x",0,200);  
  
 TF1 * cam0_crms2 = new TF1("f5","5.5+0.07*x",0,200);  
 TF1 * cam1_crms2 = new TF1("f6","10+0.05*x",0,200);  

 TF1 * new_cam0_crms1 = new TF1("nf5","10.51+0.05*x",0,200);  
 TF1 * new_cam0_crms2 = new TF1("nf4","2.64+0.046*x",0,200);  
  
 TF1 * new_cam1_crms1 = new TF1("nf5","22.35+0.038*x",0,200);  
 TF1 * new_cam1_crms2 = new TF1("nf6","4.45+0.068*x",0,200);  

 for (int i = 1; i<=num_sets; i++)
 { 
   continue; 
   cv[0][0]->cd(i); 
   cam0_maxpixel->Draw("same"); 
   cv[1][0]->cd(i); 
   cam1_maxpixel->Draw("same"); 
   cv[0][1]->cd(i); 
   cam0_crms1->Draw("same"); 
   cam0_crms2->Draw("same"); 
   
   new_cam0_crms1->SetLineColor(kRed); 
   new_cam0_crms2->SetLineColor(kRed); 
   new_cam0_crms1->Draw("same"); 
   new_cam0_crms2->Draw("same"); 


   cv[1][1]->cd(i); 
   new_cam1_crms1->SetLineColor(kRed); 
   new_cam1_crms2->SetLineColor(kRed); 
   new_cam1_crms1->Draw("same"); 
   new_cam1_crms2->Draw("same"); 
   cam1_crms1->Draw("same"); 
   cam1_crms2->Draw("same"); 
 }
  
 for (int i = 0; i < ncanv; i++)
 {
    for (int c = 0; c< 2; c++)
    {
      cv[c][i]->Update(); 
    }
 }

/*
TCanvas * lego0 = new TCanvas("lego0", "lego0", 500,500); 
TCanvas * lego1 = new TCanvas("lego1", "lego1", 500,500); 
 lego0->cd(); 
 h0->DrawCopy("LEGO4"); 
 lego1->cd(); 
 h1->DrawCopy("LEGO4"); 
*/
 double contours[] = { 0, 0.68, 1} ;

 TCanvas * rms0 = new TCanvas("rms0","rms0",500,500); 
 TCanvas * rms1 = new TCanvas("rms1","rms1",500,500); 
 TCanvas * max0 = new TCanvas("max0","max0",500,500); 
 TCanvas * max1 = new TCanvas("max1","max1",500,500); 

 TFile savedContours("wormCutHists.root","RECREATE"); 
 ContourPicker picker(3,contours); 

 TH2F * cam0rms = picker.makeContourHist(h0,"cam0rms");
 cam0rms->Clone()->Write();

 rms0->cd();
 picker.drawHistogramContours(h0,false); 
 rms0->Update(); 
 TGraph * rms0cntr = get68Contour("rms0cntr"); 
 rms0cntr->Write(); 

 TH2F * cam1rms = picker.makeContourHist(h1,"cam1rms");
 cam1rms->Clone()->Write(); 

 rms1->cd();
 picker.drawHistogramContours(h1,false); 
 rms1->Update(); 
 TGraph * rms1cntr = get68Contour("rms1cntr"); 
 rms1cntr->Write(); 

 TH2F * cam0mp = picker.makeContourHist(h2,"cam0max");
 cam0mp->Clone()->Write(); 
 max0->cd();
 picker.drawHistogramContours(h2,false); 
 max0->Update(); 
 TGraph * max0cntr = get68Contour("max0cntr"); 
 max0cntr->Write(); 

 TH2F * cam1mp = picker.makeContourHist(h3,"cam1max");
 cam1mp->Clone()->Write(); 
 max1->cd();
 picker.drawHistogramContours(h3,false); 
 max1->Update(); 
 TGraph * max1cntr = get68Contour("max1cntr"); 
 max1cntr->Write(); 

 savedContours.Close(); 


 app->Run(); 

}
