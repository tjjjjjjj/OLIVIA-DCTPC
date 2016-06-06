#include "TApplication.h"
#include "AnalysisResult.hh"
#include "../../../MaxCam/DmtpcSkimEvent.hh"
#include "TChain.h"
#include "recoilEnergy.hh"
#include "TCanvas.h"
#include "TFile.h"
#include "style.hh"
#include "TH2F.h"
#include "TF1.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLine.h"
#include "../../../MaxCam/MaxCamSRIM.hh"
static const int num_sets = 1;


int width = 500;
int height = 500;
int ncanv = 8; 
char * canvas_names[] =
{
  "passdata","worms","rbis","alphas","cfdata","energy", "paper1", "paper2"
}; 

int nhist = 7; 

int xmax[] = {300,800,800,6000,300,300,300,200}; 
int ymax[] = {10,20,20,100,10,10,10,20}; 
  

char * merged_dir = "merged/";
char * pass_prefix = "/net/zwicky/dmtpc/analysis/v3/pass/pass";

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

int main (int nargs, char **args)
{
  TApplication * app = new TApplication("cut_optimization",0,0);
  AnalysisStyle::setStyle(); 
  last_run = -1; 
  TChain * data[num_sets];

  //Set up the surface run data
  data[0] = new TChain("skim"); 
  data[0]->Add(get_full_name("merged801_850.root"));
  //data[0]->Add(get_full_name("merged851_900.root"));
  //data[0]->Add(get_full_name("merged901_950.root"));
 // data[0]->Add(get_full_name("merged951_1000.root"));
 // data[0]->Add(get_full_name("merged1001_1031.root"));

  //Set up the surface CF252 data
  data[1] = new TChain("skim");
  //data[1] = data[0]; 
  data[1]->Add(get_full_name("merged1032_1073.root"));

  //Set up the CF252 MC Data
  data[2] = new TChain("skim");
  data[2]->Add(get_full_name("merged100429002_100429010.root"));
  data[2]->Add(get_full_name("merged100502001_100502008.root"));

  //Set up the Wimp MC Data
  data[3] = new TChain("skim");
  data[3]->Add(get_full_name("merged100421011_100421020.root"));
 
  //Set up the canvases..
  TCanvas * cv[ncanv]; 
  TH2F * h[nhist];
  for (int v = 0; v < ncanv; v++) 
  {
      cv[v] = new TCanvas(canvas_names[v],canvas_names[v],width,height); 
  }

  for (int v =0; v < nhist; v++)
  {

      h[v] = new TH2F(concat("h",v),concat("h",v),20,0, xmax[v], 20, 0,ymax[v]); 
  }
 

  for (int set = 0; set < num_sets; set++)
  {

    DmtpcSkimEvent * ev = new DmtpcSkimEvent; 
    data[set]->SetBranchAddress("event",&ev); 

    int j = 0; //This index keeps track of our index inside each run
    for (int i = 0; i < data[set]->GetEntries(); i++)
    {
      data[set]->GetEntry(i);
      int run = ev->runNumber(); 

      //Reset j and load the different analysis result if run number has changed
      if (run!=last_run) {
        std::cout << "Run: "<< run << std::endl; 
        j = 0;
        //Don't need analysis result for MC data
        if (set != 2 && set!=3) 
          loadAnalysisResult(run);      
        last_run = run; 
        std::cout <<" ncam: " << ev->ncamera() << std::endl; 
      } 
      
      int ncam = ev->ncamera(); 
      int mincam = get_min_cam(run); 
      for (int c = 0; c < ncam ; c++)
      {
        for (int t = 0; t < 15; t++)
        {
          /** CUTS for all events!! **/
          if (ev->spark(c)) break;  //This should be covered by below... but whatever
          if (ev->ntracks(c) == t) break; 
          if (ev->range(c,t) == 0) continue; 
          if (ev->x(c,t) < 40 || ev->x(c,t) > 984 || ev->y(c,t)<40 || ev->y(c,t) > 984) continue; 

          double E_r = RecoilEnergy::getRecoilEnergy(ev->E(c,t), c+mincam); 
          double range = ev->range(c,t) * rangecal[c+mincam]; 
          if (E_r == 0) continue; 
       //   if (!h[3] && E_r > 800) continue; 
          
          bool edge = (ev->edge(c,t)==1); 
          bool rbi = !result->pass_rbi(j,c,t); 
          //bool worm = !passWormCuts(ev,c,t,mincam); 
          bool alpha = range > 5; 

          if (set == 0)
          {

            if (rbi) 
            {
              h[2]->Fill(E_r,range); 
            }

            if (edge && ev->npixel(c,t) > 64 && range > 1 && E_r > 50 && !rbi) 
            {
              h[3]->Fill(E_r,range); 
            }
      
            if (ev->npixel(c,t) < 64 && !rbi)
            {
              h[1]->Fill(E_r,range); 
            }

          } 
          
    /*      if (!edge && set==2 && E_r >20)
          {
            h[5]->Fill(E_r,range); 
          }
     */       

     /*
          if (!edge && set==3 && E_r> 20)
          {
            h[6]->Fill(E_r,range); 
          }
*/
        }
      }
      j++;
    }  
  }
  TFile pass("pass0-200.root"); 
  TFile neutronpass("neutronspass0-800.root"); 
  
  TChain * wimpskim = new TChain("skim"); 
  wimpskim->Add("pass/00100000.root"); 
  wimpskim->Add("pass/00100002.root"); 
//  wimpskim->Add("pass/00100004.root"); 
//  wimpskim->Add("pass/00100005.root"); 
//  wimpskim->Add("pass/00100006.root"); 
  cout << wimpskim->GetEntries() << endl;


  TChain * cfskim = new TChain("skim"); 
  cfskim->Add("pass/100502001.root");
  cfskim->Add("pass/100502002.root");
//  cfskim->Add("pass/100502003.root");
//  cfskim->Add("pass/100502004.root");
//  cfskim->Add("pass/100502006.root");
//  cfskim->Add("pass/100502007.root");
//  cfskim->Add("pass/100502008.root");
//  cfskim->Add("pass/100429002.root");
//  cfskim->Add("pass/100429003.root");
//  cfskim->Add("pass/100429004.root");
 // cfskim->Add("pass/100429005.root");
  //cfskim->Add("pass/100429006.root");
 // cfskim->Add("pass/100429007.root");
 // cfskim->Add("pass/100429008.root");

  TTree * skim = (TTree*) pass.Get("skim");  
  TTree * nskim = (TTree*) neutronpass.Get("skim");  
  int cam; 
  int track; 
  DmtpcSkimEvent * ev1 = new DmtpcSkimEvent; 
  skim->SetBranchAddress("event",&ev1); 
  skim->SetBranchStatus("*clusters*",0);
  skim->SetBranchStatus("*orig_event*",0);
  skim->SetBranchAddress("track",&track); 
  skim->SetBranchAddress("cam",&cam); 

  DmtpcSkimEvent * ev2 = new DmtpcSkimEvent; 
  DmtpcSkimEvent * ev3 = new DmtpcSkimEvent; 
  DmtpcSkimEvent * ev4 = new DmtpcSkimEvent; 
  nskim->SetBranchStatus("*clusters*",0);
  nskim->SetBranchStatus("*orig_event*",0);
  nskim->SetBranchAddress("event",&ev2); 
  nskim->SetBranchAddress("track",&track); 
  nskim->SetBranchAddress("cam",&cam); 
  
  wimpskim->SetBranchAddress("event",&ev4); 
  wimpskim->SetBranchAddress("track",&track); 
  wimpskim->SetBranchStatus("*clusters*",0);
  wimpskim->SetBranchStatus("*orig_event*",0);
  wimpskim->SetBranchAddress("cam",&cam); 

  cfskim->SetBranchAddress("event",&ev3); 
  cfskim->SetBranchAddress("track",&track); 
  cfskim->SetBranchStatus("*clusters*",0);
  cfskim->SetBranchStatus("*orig_event*",0);
  cfskim->SetBranchAddress("cam",&cam); 

  TH1F * er_dists[2];
  TH1F * er_dist = new TH1F("er_dist","er_dist",40,0,200); 
  er_dists[0] = new TH1F("er_dist0","er_dist0",40,0,200); 
  er_dists[1] = new TH1F("er_dist1","er_dist1",40,0,200); 

  for (int i = 0; i < skim->GetEntries(); i++)
  {
     
    skim->GetEntry(i); 
    double range = ev1->range(cam,track) * rangecal[cam]; 
    double E_r = RecoilEnergy::getRecoilEnergy(ev1->E(cam,track), cam);  
    h[0]->Fill(E_r,range); 
    er_dist->Fill(E_r); 
    er_dists[cam]->Fill(E_r); 
  }

  for (int i = 0; i < nskim->GetEntries(); i++)
  {
    nskim->GetEntry(i); 
    double range = ev2->range(cam,track) * rangecal[cam]; 
    double E_r = RecoilEnergy::getRecoilEnergy(ev2->E(cam,track), cam);  
    h[4]->Fill(E_r,range); 
  }
  
  for (int i = 0; i < cfskim->GetEntries(); i++)
  {
    cfskim->GetEntry(i); 
    double range = ev3->range(cam,track) * rangecal[cam]; 
    double E_r = RecoilEnergy::getRecoilEnergy(ev3->E(cam,track), cam);  
    h[5]->Fill(E_r,range); 
  }

  for (int i = 0; i < wimpskim->GetEntries(); i++)
  {
    wimpskim->GetEntry(i); 
    double range = ev4->range(cam,track) * rangecal[cam]; 
    double E_r = RecoilEnergy::getRecoilEnergy(ev4->E(cam,track), cam);  
    h[6]->Fill(E_r,range); 
  }
  
  MaxCamSRIM he("SRIM_He_in_CF4_100Torr");
  MaxCamSRIM fl("SRIM_F_in_CF4_100Torr");
  MaxCamSRIM c("SRIM_C_in_CF4_100Torr");
  he.setPressure(75);
  fl.setPressure(75);
  c.setPressure(75);

  TGraph * he_g = he.getRangeVsEnergy(1); 
  he_g->SetLineStyle(1);
  TGraph * fl_g = fl.getRangeVsEnergy(1); 
  fl_g->SetLineStyle(2);
  TGraph * c_g = c.getRangeVsEnergy(1); 
  c_g->SetLineStyle(3);

  /* Paper Plot 1 */
  cv[7]->cd();
  h[2]->SetXTitle("Recoil Energy (keV)"); 
  h[2]->SetYTitle("Projected Range (mm)"); 
  h[1]->SetMarkerSize(0.1); 
  h[1]->SetLineColor(28); 
  h[1]->SetLineWidth(2); 
  h[2]->SetMarkerColor(6); 
  h[2]->SetLineColor(6); 
  h[2]->SetMarkerSize(0.4); 
  h[2]->SetMarkerStyle(22);  
  gStyle->SetPalette(9); 
  h[2]->Draw("colz"); 
  h[3]->SetMarkerSize(0.4); 
  h[3]->SetMarkerColor(4); 
  h[3]->Draw("same p"); 
  h[1]->Draw("same box"); 

 /*Paper Plot 2*/
  cv[6]->cd();
  h[5]->SetXTitle("Recoil Energy (keV)"); 
  h[5]->SetYTitle("Projected Range (mm)"); 
  h[0]->GetXaxis()->SetRange(0,200); 
  h[0]->SetMarkerSize(1); 
  h[0]->SetMarkerColor(kRed); 
  h[6]->SetMarkerColor(kBlue); 
  h[6]->SetMarkerSize(0.4); 
  h[5]->Draw("colz"); 
  h[6]->Draw("same p"); 
  h[0]->Draw("same p"); 

  /* wimp mc/data*/
  cv[0]->cd(); 
  h[0]->SetXTitle("Recoil Energy (keV)"); 
  h[0]->SetYTitle("Projected Range (mm)"); 
  h[0]->SetMarkerSize(0.5); 
  h[0]->SetMarkerColor(1); 
  h[6]->SetMarkerSize(0.2); 
  h[6]->SetLineColor(3); 
  h[6]->Draw("colz");  
  h[0]->Draw("same p");  


  /*Worms */
  cv[1]->cd(); 
  h[1]->SetXTitle("Recoil Energy (keV)"); 
  h[1]->SetYTitle("Projected Range (mm)"); 
  h[1]->Draw("box");  
 

  /*rbi*/ 
  cv[2]->cd(); 
  h[2]->SetXTitle("Recoil Energy (keV)"); 
  h[2]->SetYTitle("Projected Range (mm)"); 
  h[2]->Draw("box");  


  /*alphas*/ 
  cv[3]->cd(); 
  h[3]->SetXTitle("Recoil Energy (keV)"); 
  h[3]->SetYTitle("Projected Range (mm)"); 
  //h[3]->SetMarkerSize(0.1); 
  h[3]->SetMarkerColor(4); 
  h[3]->Draw("p");  



  //Cf252 mc / data
  cv[4]->cd(); 
  h[4]->SetMarkerSize(0.6);
  h[4]->SetMarkerStyle(20);
  h[4]->SetMarkerColor(2); 
  h[5]->Draw("colz");  
  h[4]->Draw("p same");  
  h[5]->SetLineColor(7); 
  h[4]->SetXTitle("Recoil Energy (keV)"); 
  h[4]->SetYTitle("Projected Range (mm)"); 

  //Stuff for histogram plot
  TFile spectrum_f("surfaceNeutronSpectrum_wr2.root"); 
  TH1F *spectrum = (TH1F*) spectrum_f.Get("recoilDistribution"); 

  TLine *l = new TLine; 
  l->SetLineStyle(1); 
  l->SetLineWidth(3); 
  l->SetLineColor(kRed); 
  cv[5]->cd();
  er_dist->SetXTitle("Recoil Energy (keV)"); 
  er_dist->SetYTitle("counts / 10 keV "); 
  spectrum->SetLineColor(kMagenta); 
 // spectrum->SetFillColor(kGreen); 

  er_dist->Draw("E1");
 // spectrum->Draw("E4 SAME");
  spectrum->Draw("SAME");
  er_dist->Draw("E1 SAME");
//  l->DrawLine(80,0,80,20);
//  l->DrawLine(200,0,200,20);
 
  TCanvas * erdists_c = new TCanvas("erdists","erdists", 800,400); 
  erdists_c->SaveAs("erdists.pdf");
  erdists_c->Divide(2,1); 
  erdists_c->cd(1);
  er_dists[0]->Draw("E1");   
  erdists_c->cd(2);
  er_dists[1]->Draw("E1");   


  for (int i = 0; i <8; i++)
  {
    cv[i]->cd();
//    he_g->Draw(); 
//    fl_g->Draw(); 
//    c_g->Draw(); 
  } 
  erdists_c->Update(); 

  for(int i = 0; i < 8; i++)
  {
    cv[i]->Update();  
    cv[i]->SaveAs(".plot.pdf");
  }

  
  app->Run(); 
}
