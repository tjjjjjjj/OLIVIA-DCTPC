#define DCTPCTree_cxx
#include "../input_files/DCTPCTree.h"
#include "../input_files/DCTPC_runtree.h"
#include "../input_files/Calibration_constants_LittleDCTPC.h";
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
using namespace std;


void DCTPCTree::Loop()
{
  TStopwatch timer;
  gROOT->SetStyle("Plain");
  gStyle->SetEndErrorSize(3);
  gStyle->SetPalette(1,0);
  gStyle->SetLineWidth(2);
  gStyle->SetHistLineWidth(2);
  gStyle->SetOptStat(kTRUE);
  gStyle->SetOptFit(kTRUE);
  TH1::AddDirectory(false);
  
  TFile *outtree = new TFile("$LittleDCTPC_physics_input");//$LittleDCTPC_physics_input
  TTree *dctreepc = (TTree*)outtree->Get("dctpc_eventinfo");
  DCTPCTree aStep(dctreepc);

  Int_t Exposure_sec;
  TTree *dctreepc2 = (TTree*)outtree->Get("dctpc_runinfo");       
  dctreepc2->SetBranchAddress("Exposure_sec", &Exposure_sec);
  Int_t nentries2 = (Int_t)dctreepc2->GetEntries();  
  double totalexposure=0.;
  
  for (int runcounter = 0; runcounter<nentries2; runcounter++)
    {
      dctreepc2->GetEntry(runcounter);
      totalexposure+=Exposure_sec;
    }
 
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TH1D *hist_rangeccd=new TH1D("CCD Range", "CCD Range", 200, 0, 200); //do I want to fit a Gaussian to this
  TH1D *hist_rangeccdlog=new TH1D("CCD Range Log Scale", "CCD Range Log Scale", 200, 0, 200);
  TH1D *hist_tracklength= new TH1D("Track Length", "Track Length", 150, 0, 300);
  
  //Length and energy plots
  TH2D *hist_range_vs_ccd=new TH2D ("Range vs CCD Energy", "Range vs CCD Energy", 200, 0, 10000, 150, 0, 150);
  TH2D *hist_range_vs_mesh=new TH2D ("Range vs Mesh Energy", "Range vs Mesh Energy", 200, 0, 10000, 150, 0, 150);
  TH2D *hist_range_vs_anode = new TH2D ("Range vs Anode Energy", "Range vs Anode Energy", 200, 0, 2000, 150, 0, 150);
  TH2D *hist_length_vs_ccd=new TH2D ("Length vs CCD Energy", "Length vs CCD Energy", 200, 0, 10000, 100, 0, 100);
  TH2D *hist_length_vs_mesh=new TH2D ("Length vs Mesh Energy", "Length vs Mesh Energy", 100, 0, 5000, 100, 0, 100);
  TH2D *hist_length_vs_anode=new TH2D ("Length vs Anode Energy", "Length vs Anode Energy", 200, 0, 2000, 100, 0, 100);
  TH2D *hist_x_length_vs_y_length=new TH2D ("Horizontal vs Vertical Length", "Horizontal vs Vertical Length", 400, 0, 400, 400, 0, 400);
  TH2D *hist_projected_length=new TH2D ("Projected Length", "Projected Length", 100, 0, 5000, 300, 0, 300);

   //////////////////////////////////////////////////////////////

  double seqcalib [19] = {0.955114 , 1.06241 , 1.07105 , 1.07669 , 1.02827 , 1.15574 , 1.18828 , 1.08779 , 0.998239 , 1.05118 , 1.03128 , 1.02639 , 0.992746 , 0.925486 , 1.1664 , 1.05791, 0.989169};

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//create event loop
  Long64_t nentries = dctreepc->GetEntries();
  cout << "Number of Entries: " << nentries << endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  int cutnum=-1;
  int neutronnumber=0;
  int passall=0;
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for (int event = 0; event<nentries; event++)
    {
      
      aStep.GetEntry(event); 
      if(event%1000000==0)
	cout<<((double)event/(double)nentries)*100.<<"%"<<endl;

      r2 = (aStep.Track_x_pix*aStep.Track_x_pix+aStep.Track_y_pix*aStep.Track_y_pix)*pow(MMPERPIXEL/10.,2);
      weight=1./pow(RADIAL_PARAM0_CALIB/(RADIAL_PARAM1_CALIB+r2),RADIAL_PARAM2_CALIB);//radial calibration 
      rho_start=(sqrt(((aStep.Track_x_start_pix*MMPERPIXEL)*(aStep.Track_x_start_pix*MMPERPIXEL))+((aStep.Track_y_start_pix*MMPERPIXEL)*(aStep.Track_y_start_pix*MMPERPIXEL))));
      rho_end=(sqrt(((aStep.Track_x_end_pix*MMPERPIXEL)*(aStep.Track_x_end_pix*MMPERPIXEL))+((aStep.Track_y_end_pix*MMPERPIXEL)*(aStep.Track_y_end_pix*MMPERPIXEL))));
      
      meshcalib=MESH_CALIB*seqcalib[aStep.SequenceNum-2];
      ccdcalib=CCD_CALIB*seqcalib[aStep.SequenceNum-2];
      anodecalib=ANODE_CALIB*seqcalib[aStep.SequenceNum-2];

      verticallength=((aStep.Mesh_R0time_samp+aStep.Mesh_F10time_samp)*driftspeed);
      horizontallength=aStep.Track_range_pix*MMPERPIXEL;
      tracklength=(sqrt((verticallength*verticallength)+(horizontallength*horizontallength)));
      lengthcalib=LENGTH_CALIB*seqcalib[aStep.SequenceNum-2];
      lengthenergy=(El_0 + (El_1*tracklength) + (El_2*tracklength*tracklength))*lengthcalib;
      ccdenergy=((Eccd_0)+(Eccd_1*aStep.Etrack_kev*weight)+(Eccd_2*aStep.Etrack_kev*weight*aStep.Etrack_kev*weight))*ccdcalib;

      delta_anodeccd=((aStep.Etrig_kev*anodecalib)-(ccdenergy))/((aStep.Etrig_kev*anodecalib)+(ccdenergy));
      delta_meshccd=((aStep.Emesh_kev*meshcalib)-(ccdenergy))/((aStep.Emesh_kev*meshcalib)+(ccdenergy));
      delta_anodemesh=((aStep.Etrig_kev*anodecalib)-(aStep.Emesh_kev*meshcalib))/((aStep.Etrig_kev*anodecalib)+(aStep.Emesh_kev*meshcalib));
      delta_meshlength=((lengthenergy)-(aStep.Emesh_kev*meshcalib))/((lengthenergy)+(aStep.Emesh_kev*meshcalib));
      delta_anodelength=((aStep.Etrig_kev*anodecalib)-(lengthenergy))/((aStep.Etrig_kev*anodecalib)+(lengthenergy));
      delta_ccdlength=((ccdenergy)-(lengthenergy))/((ccdenergy)+(lengthenergy));
      mesh_over_ccd=((aStep.Emesh_kev*meshcalib)/(ccdenergy));
           
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////
       
       if(aStep.Edge==1)
	{
	  cutnum=1;
	  continue;
	}

      if(aStep.EventNum-aStep.LastSpark<7)
	{
	  cutnum=2;
	  continue;
	}
      
      if(aStep.Track_range_pix*MMPERPIXEL<=0.)
	{
	  cutnum=3;
	  continue;
	}
      
      if(aStep.Etrig_kev<400)
	{
	  cutnum=4;
	  continue;
	}

      if(aStep.Emesh_kev<500)
	{
	  cutnum=5;
	  continue;
	}
      
      if(aStep.Etrack_kev*weight<300)
	{
	  cutnum=6;
	  continue;
	}
      
      if(aStep.Track_maxpixel_ccdadu<72 ||
				     aStep.Track_maxpixel_ccdadu>200)
	{
	  cutnum=7;
	  continue;
	}
      
      if(rho_start>100 ||
	 rho_end>100)
	{
	  cutnum=8;
	  continue;
	}
      
      /* if(tracklength>72)//change to greater than 72 for neutrons
	{continue;}*/
      
     
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	 	 
      hist_rangeccd->Fill(aStep.Track_range_pix*MMPERPIXEL);
      hist_rangeccdlog->Fill(aStep.Track_range_pix*MMPERPIXEL);
      hist_tracklength->Fill(tracklength); 
      hist_range_vs_ccd->Fill(ccdenergy, aStep.Track_range_pix*MMPERPIXEL);
      hist_range_vs_mesh->Fill(aStep.Emesh_kev*meshcalib, aStep.Track_range_pix*MMPERPIXEL);
      hist_range_vs_anode->Fill(aStep.Etrig_kev*anodecalib, aStep.Track_range_pix*MMPERPIXEL);
      hist_length_vs_ccd->Fill(ccdenergy, tracklength);
      hist_length_vs_mesh->Fill(aStep.Emesh_kev*meshcalib, tracklength);
      hist_length_vs_anode->Fill(aStep.Etrig_kev*anodecalib, tracklength);
      hist_x_length_vs_y_length->Fill(horizontallength, verticallength);
      
      //above I am filling the histograms with data that has passed all the cuts - hopefully neutrons!
         
      cutnum=0;
      neutronnumber++;
        
    }

  hist_projected_length->Fill(10, 0.512);
  hist_projected_length->Fill(15, 0.753);
  hist_projected_length->Fill(20, 0.977);
  hist_projected_length->Fill(25, 1.19);
  hist_projected_length->Fill(30, 1.38);
  hist_projected_length->Fill(35, 1.57);
  hist_projected_length->Fill(40, 1.74);
  hist_projected_length->Fill(45, 1.91);
  hist_projected_length->Fill(50, 2.07);
  hist_projected_length->Fill(55, 2.22);
  hist_projected_length->Fill(60, 2.37);
  hist_projected_length->Fill(65, 2.51);
  hist_projected_length->Fill(70, 2.65);
  hist_projected_length->Fill(80, 2.91);
  hist_projected_length->Fill(90, 3.16);
  hist_projected_length->Fill(100, 3.39);
  hist_projected_length->Fill(110, 3.62);
  hist_projected_length->Fill(120, 3.84);
  hist_projected_length->Fill(130, 4.05);
  hist_projected_length->Fill(140, 4.26);
  hist_projected_length->Fill(150, 4.46);
  hist_projected_length->Fill(160, 4.65);
  hist_projected_length->Fill(170, 4.84);
  hist_projected_length->Fill(180, 5.02);
  hist_projected_length->Fill(200, 5.37);
  hist_projected_length->Fill(250, 6.20);
  hist_projected_length->Fill(300, 6.95);
  hist_projected_length->Fill(350, 7.66);
  hist_projected_length->Fill(400, 8.34);
  hist_projected_length->Fill(450, 8.99);
  hist_projected_length->Fill(500, 9.62);
  hist_projected_length->Fill(550, 10.23);
  hist_projected_length->Fill(600, 10.84);
  hist_projected_length->Fill(650, 11.44);
  hist_projected_length->Fill(700, 12.04);
  hist_projected_length->Fill(800, 13.22);
  hist_projected_length->Fill(900, 14.42);
  hist_projected_length->Fill(1000, 15.62);
  hist_projected_length->Fill(1100, 16.85);
  hist_projected_length->Fill(1200, 18.10);
  hist_projected_length->Fill(1300, 19.38);
  hist_projected_length->Fill(1400, 20.69);
  hist_projected_length->Fill(1500, 22.02);
  hist_projected_length->Fill(1600, 23.39);
  hist_projected_length->Fill(1700, 24.80);
  hist_projected_length->Fill(1800, 26.23);
  hist_projected_length->Fill(2000, 29.22);
  hist_projected_length->Fill(2250, 33.14);
  hist_projected_length->Fill(2500, 37.30);
  hist_projected_length->Fill(2750, 41.70);
  hist_projected_length->Fill(3000, 46.33);
  hist_projected_length->Fill(3250, 51.21);
  hist_projected_length->Fill(3500, 56.32);
  hist_projected_length->Fill(3750, 61.68);
  hist_projected_length->Fill(4000, 67.28);
  hist_projected_length->Fill(4500, 79.19);
  hist_projected_length->Fill(5000, 92.05);
  hist_projected_length->Fill(5500, 105.87);
  hist_projected_length->Fill(6000, 120.62);
  hist_projected_length->Fill(6500, 136.30);
  hist_projected_length->Fill(7000, 152.89);
  hist_projected_length->Fill(8000, 188.77);
  hist_projected_length->Fill(9000, 228.07);
  hist_projected_length->Fill(10000, 270.78);
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  Double_t fitf(Double_t *x, Double_t *par)
  {
    Double_t fitval = par[0]*x[0];
    return fitval;
  }   
  
  TF1 *func1 = new TF1("pol1", "pol1", 200, 1200);//anode vs ccd
  
  TF1 *func2 = new TF1("pol1", "pol1", 0, 3000);//mesh vs ccd

  TF1 *func3 = new TF1("pol1", "pol1", 250, 1200);//anode vs mesh

  TF1 *func4 = new TF1("pol2", "pol2", 0, 5000);//projected length

  //above func 1/2/3 is where I define line of best fit for 2D energy plots

  TF1 *gaus1 = new TF1("gaus", "gaus", -1, 1);

  TF1 *gaus2 = new TF1("gaus", "gaus", -1, 1);

  TF1 *gaus3 = new TF1("gaus", "gaus", -1, 1);

  //above gaus 1/2/3 is where I define the fits for the delta distributions

   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  new TCanvas;
  hist_range_vs_ccd->SetTitle("Range vs CCD Energy");
  hist_range_vs_ccd->SetXTitle("Calibrated CCD Energy/keV");
  hist_range_vs_ccd->SetYTitle("Range/mm");
  hist_range_vs_ccd->Draw("COLZ");

  new TCanvas;
  hist_range_vs_mesh->SetTitle("Range vs Mesh Energy");
  hist_range_vs_mesh->SetXTitle("Mesh Energy/keV");
  hist_range_vs_mesh->SetYTitle("Range/mm");
  hist_range_vs_mesh->Draw("COLZ");
  hist_projected_length->SetMarkerStyle(6);
  hist_projected_length->Fit(func4, "r n");
  hist_projected_length->Draw("SAME");
  func4->Draw("SAME");

  new TCanvas;
  hist_range_vs_anode->SetTitle("Range vs Anode Enegry");
  hist_range_vs_anode->SetXTitle("Calibrated Anode Energy/keV");
  hist_range_vs_anode->SetYTitle("Range/mm");
  hist_range_vs_anode->Draw("COLZ");

  new TCanvas;
  hist_rangeccd->SetTitle("CCD Range");
  hist_rangeccd->SetXTitle("CCD Range/mm");
  hist_rangeccd->SetYTitle("Counts");
  hist_rangeccd->Draw();

  new TCanvas;
  c1_n5->SetLogy();
  hist_rangeccdlog->SetTitle("CCD Range Log Scale");
  hist_rangeccdlog->SetXTitle("CCD Range/mm");
  hist_rangeccdlog->SetYTitle("Log Counts");
  hist_rangeccdlog->Draw();
  
  new TCanvas;
  hist_tracklength->SetTitle("Track Length");
  hist_tracklength->SetXTitle("Track Length/mm");
  hist_tracklength->SetYTitle("Counts");
  hist_tracklength->Draw();

  new TCanvas;
  hist_length_vs_ccd->SetTitle("Track Length vs CCD Energy");
  hist_length_vs_ccd->SetXTitle("CCD Energy/keV");
  hist_length_vs_ccd->SetYTitle("Track Length/mm");
  hist_length_vs_ccd->Draw("COLZ");
  hist_projected_length->SetMarkerStyle(6);
  hist_projected_length->Fit(func4, "r n");
  hist_projected_length->Draw("SAME");
  func4->Draw("SAME");
  
  new TCanvas;
  hist_length_vs_mesh->SetTitle("Track Length vs Mesh Energy");
  hist_length_vs_mesh->SetXTitle("Mesh Energy/keV");
  hist_length_vs_mesh->SetYTitle("Track Length/mm");
  hist_length_vs_mesh->Draw("COLZ");
  hist_projected_length->SetMarkerStyle(6);
  hist_projected_length->Fit(func4, "r n");
  hist_projected_length->Draw("SAME");
  func4->Draw("SAME");
  
  new TCanvas;
  hist_length_vs_anode->SetTitle("Track Length vs Anode Energy");
  hist_length_vs_anode->SetXTitle("Anode Energy/keV");
  hist_length_vs_anode->SetYTitle("Track Length/mm");
  hist_length_vs_anode->Draw("COLZ");
  
  /* new TCanvas;
  hist_x_length_vs_y_length->SetTitle("Horizontal Length vs Vertical Length");
  hist_x_length_vs_y_length->SetXTitle("Horizontal Length / mm");
  hist_x_length_vs_y_length->SetYTitle("Vertical Length / mm");
  hist_x_length_vs_y_length->Draw("COLZ");
  */

   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  cout<<"Drift Speed= "<<driftspeed<<" mm/sec"<<endl;
  cout <<"Mean CCD Range= "<< hist_rangeccd->GetMean()<<" mm"<<endl;
  cout<<"Mean Track Length= "<<hist_tracklength->GetMean()<<" mm"<<endl;
 
  cout<<""<<endl;

  cout<<"There were "<<neutronnumber<<" neutrons"<<endl;
  cout<<"Total Exposure time now "<<totalexposure/(60*60*24)<<" days"<<endl;
  cout<<"There were ~  "<<neutronnumber/(totalexposure/(60*60*24))<<" neutrons per day"<<endl;

  cout<<" "<<endl;
  cout<<"The fit on the Range vs Mesh and Length vs Mesh plots is a plot of what the projected length of an alpha of a certain energy is in the DCTPC gas based on a simulation produced by Adrien"<<endl;
  
}

double fitPoly12Percent(double E)
{
  return sqrt(0.9741*0.9741*E*E+113.835*113.835)-62.675;
  //return -11.50+39.02*E+0.03459*E*E;                                                    
  //return 469 + 10*E + 0.548 *E*E -3.39e-3*E*E*E+7.5e-6*E*E*E*E;          
  
}               

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
