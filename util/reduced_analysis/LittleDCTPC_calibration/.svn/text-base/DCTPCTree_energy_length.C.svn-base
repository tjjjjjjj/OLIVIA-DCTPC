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
  gStyle->SetOptStat(kFALSE);
  //   gStyle->SetOptFit(kFALSE);
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

  //Length and energy plots
   TH2D *hist_length_vs_mesh=new TH2D ("Length vs Mesh Energy", "Length vs Mesh Energy", 100, 0, 5000, 100, 0, 100);
  TH2D *hist_projected_length=new TH2D ("Projected Length", "Projected Length", 100, 0, 5000, 300, 0, 300);
  TH2D *hist_energy_vs_length=new TH2D ("Mesh Energy vs Length", "Mesh Energy vs Length", 100, 0, 100, 100, 0, 5000);
  TH2D *hist_projected_length2=new TH2D("Projected Energy", "Projected Energy", 300, 0, 300, 100, 0, 5000);

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

      double r2 = (aStep.Track_x_pix*aStep.Track_x_pix+aStep.Track_y_pix*aStep.Track_y_pix)*pow(MMPERPIXEL/10.,2);
      double weight=1./pow(RADIAL_PARAM0_CALIB/(RADIAL_PARAM1_CALIB+r2),RADIAL_PARAM2_CALIB);//radial calibration 
      double rho_start=(sqrt(((aStep.Track_x_start_pix*MMPERPIXEL)*(aStep.Track_x_start_pix*MMPERPIXEL))+((aStep.Track_y_start_pix*MMPERPIXEL)*(aStep.Track_y_start_pix*MMPERPIXEL))));
      double rho_end=(sqrt(((aStep.Track_x_end_pix*MMPERPIXEL)*(aStep.Track_x_end_pix*MMPERPIXEL))+((aStep.Track_y_end_pix*MMPERPIXEL)*(aStep.Track_y_end_pix*MMPERPIXEL))));
      double verticallength=((aStep.Mesh_totaltime_samp)*driftspeed);
      double horizontallength=aStep.Track_range_pix*MMPERPIXEL;
      double tracklength=(sqrt((verticallength*verticallength)+(horizontallength*horizontallength)));
      double meshcalib=MESH_CALIB*seqcalib[aStep.SequenceNum-2];
      double ccdcalib=CCD_CALIB*seqcalib[aStep.SequenceNum-2];
      double anodecalib=ANODE_CALIB*seqcalib[aStep.SequenceNum-2];
      double delta_anodeccd=((aStep.Etrig_kev*anodecalib)-(aStep.Etrack_kev*weight*ccdcalib))/((aStep.Etrig_kev*anodecalib)+(aStep.Etrack_kev*weight*ccdcalib));
      double delta_meshccd=((aStep.Emesh_kev*meshcalib)-(aStep.Etrack_kev*weight*ccdcalib))/((aStep.Emesh_kev*meshcalib)+(aStep.Etrack_kev*weight*ccdcalib));
      double delta_anodemesh=((aStep.Etrig_kev*anodecalib)-(aStep.Emesh_kev*meshcalib))/((aStep.Etrig_kev*anodecalib)+(aStep.Emesh_kev*meshcalib));
      double mesh_over_ccd=((aStep.Emesh_kev*meshcalib)/(aStep.Etrack_kev*weight*ccdcalib));
           
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
      
      if(fabs(aStep.Mesh_R0time_samp)>0.5 ||
	 fabs(aStep.Mesh_R10time_samp)>0.5 ||
	 fabs(aStep.Mesh_R25time_samp)>0.5 ||
	 fabs(aStep.Mesh_R50time_samp)>0.5 ||
	 fabs(aStep.Mesh_R75time_samp)>0.5 ||
	 fabs(aStep.Mesh_R90time_samp)>0.5)
	{
	  cutnum=8;
	  continue;
	}
      
      
      if(fabs(aStep.Mesh_F0time_samp)>0.5 ||
	 fabs(aStep.Mesh_F10time_samp)>0.5 ||
	 fabs(aStep.Mesh_F25time_samp)>0.5 ||
	 fabs(aStep.Mesh_F50time_samp)>0.5 ||
	 fabs(aStep.Mesh_F75time_samp)>0.5 ||
	 fabs(aStep.Mesh_F90time_samp)>0.5)
	{
	  cutnum=9;
	  continue;
	}
      
      if(aStep.Etrig_kev<1250 &&
			 TMath::Abs(delta_anodeccd)>0.4)
	{
	  cutnum=10;
	  continue;
	}
      
      if(TMath::Abs(delta_meshccd)>0.4)
	{
	  cutnum=11;
	  continue;
	}
      
      if(aStep.Etrig_kev<1250 &&
			 TMath::Abs(delta_anodemesh)>0.2)
	{
      cutnum=12;
      continue;
	}

      if(rho_start>100 ||
	 rho_end>100)
	{
	  cutnum=13;
	  continue;
	}
      
      if(aStep.Track_range_pix*MMPERPIXEL>70)
	{
	  cutnum=14;
	  continue;
	}
      
      if(aStep.Emesh_kev*meshcalib>4000 ||
	 aStep.Emesh_kev>5900)
	{
	  cutnum=15;
	  continue;
	}
     
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	 	 
     
      hist_length_vs_mesh->Fill(aStep.Emesh_kev*meshcalib, tracklength);
      hist_energy_vs_length->Fill(tracklength, aStep.Emesh_kev*meshcalib);
     
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

  hist_projected_length2->Fill( 0.512, 10);
  hist_projected_length2->Fill( 0.753, 15);
  hist_projected_length2->Fill( 0.977, 20);
  hist_projected_length2->Fill( 1.19, 25);
  hist_projected_length2->Fill( 1.38, 30);
  hist_projected_length2->Fill( 1.57, 35);
  hist_projected_length2->Fill( 1.74, 40);
  hist_projected_length2->Fill( 1.91, 45);
  hist_projected_length2->Fill( 2.07, 50);
  hist_projected_length2->Fill( 2.22, 55);
  hist_projected_length2->Fill( 2.37, 60);
  hist_projected_length2->Fill( 2.51, 65);
  hist_projected_length2->Fill( 2.65, 70);
  hist_projected_length2->Fill( 2.91, 80);
  hist_projected_length2->Fill( 3.16, 90);
  hist_projected_length2->Fill( 3.39, 100);
  hist_projected_length2->Fill( 3.62, 110);
  hist_projected_length2->Fill( 3.84, 120);
  hist_projected_length2->Fill( 4.05, 130);
  hist_projected_length2->Fill( 4.26, 140);
  hist_projected_length2->Fill(4.46, 150);
  hist_projected_length2->Fill( 4.65, 160);
  hist_projected_length2->Fill( 4.84, 170);
  hist_projected_length2->Fill( 5.02, 180);
  hist_projected_length2->Fill( 5.37, 200);
  hist_projected_length2->Fill( 6.20, 250);
  hist_projected_length2->Fill( 6.95, 300);
  hist_projected_length2->Fill( 7.66, 350);
  hist_projected_length2->Fill( 8.34, 400);
  hist_projected_length2->Fill( 8.99, 450);
  hist_projected_length2->Fill( 9.62, 500);
  hist_projected_length2->Fill( 10.23, 550);
  hist_projected_length2->Fill( 10.84, 600);
  hist_projected_length2->Fill( 11.44, 650);
  hist_projected_length2->Fill( 12.04, 700);
  hist_projected_length2->Fill( 13.22, 800);
  hist_projected_length2->Fill( 14.42, 900);
  hist_projected_length2->Fill( 15.62, 1000);
  hist_projected_length2->Fill( 16.85, 1100);
  hist_projected_length2->Fill( 18.10, 1200);
  hist_projected_length2->Fill( 19.38, 1300);
  hist_projected_length2->Fill( 20.69, 1400);
  hist_projected_length2->Fill( 22.02, 1500);
  hist_projected_length2->Fill( 23.39, 1600);
  hist_projected_length2->Fill( 24.80, 1700);
  hist_projected_length2->Fill( 26.23, 1800);
  hist_projected_length2->Fill( 29.22, 2000);
  hist_projected_length2->Fill( 33.14, 2250);
  hist_projected_length2->Fill( 37.30, 2500);
  hist_projected_length2->Fill( 41.70, 2750);
  hist_projected_length2->Fill( 46.33, 3000);
  hist_projected_length2->Fill( 51.21, 3250);
  hist_projected_length2->Fill( 56.32, 3500);
  hist_projected_length2->Fill( 61.68, 3750);
  hist_projected_length2->Fill( 67.28, 4000);
  hist_projected_length2->Fill( 79.19, 4500);
  hist_projected_length2->Fill( 92.05, 5000);
  hist_projected_length2->Fill( 105.87, 5500);
  hist_projected_length2->Fill( 120.62, 6000);
  hist_projected_length2->Fill( 136.30, 6500);
  hist_projected_length2->Fill( 152.89, 7000);
  hist_projected_length2->Fill( 188.77, 8000);
  hist_projected_length2->Fill(228.07, 9000);
  hist_projected_length2->Fill( 270.78, 10000);
  
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

  TF1 *func5 = new TF1("pol2", "pol2", 0, 5000);//polynomial length vs energy fit

  TF1 *func6 = new TF1("pol2", "pol2", 0, 100);

  TF1 *func7 = new TF1("pol2", "pol2", 0, 100);

  //above func 1/2/3 is where I define line of best fit for 2D energy plots

  TF1 *gaus1 = new TF1("gaus", "gaus", -1, 1);

  TF1 *gaus2 = new TF1("gaus", "gaus", -1, 1);

  TF1 *gaus3 = new TF1("gaus", "gaus", -1, 1);

  //above gaus 1/2/3 is where I define the fits for the delta distributions

  TArc *innerveto = new TArc(512*MMPERPIXEL,512*MMPERPIXEL,124);
  TArc *outerveto = new TArc(512*MMPERPIXEL,512*MMPERPIXEL,134);

  //above is where I define the veto ring shape for track position plot

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   
  new TCanvas;
  hist_length_vs_mesh->SetTitle("Track Length vs Mesh Energy");
  hist_length_vs_mesh->SetXTitle("Mesh Energy/keV");
  hist_length_vs_mesh->SetYTitle("Track Length/mm");
  hist_length_vs_mesh->Fit(func5, "r n");
  func5->SetLineColor(kRed);
  hist_length_vs_mesh->Draw("COLZ");
  func5->Draw("SAME");
  hist_projected_length->SetMarkerStyle(6);
  hist_projected_length->Fit(func4, "r n");
  hist_projected_length->Draw("SAME");
  func4->Draw("SAME");

  new TCanvas;
  hist_energy_vs_length->SetTitle("Mesh Energy vs Track Length");
  hist_energy_vs_length->SetXTitle("Track Length/mm");
  hist_energy_vs_length->SetYTitle("Mesh Energy/keV");
  hist_energy_vs_length->Fit(func6, "r n");
  func6->SetLineColor(kRed);
  hist_energy_vs_length->Draw("COLZ");
  func6->Draw("SAME");
  hist_projected_length2->SetMarkerStyle(6);
  hist_projected_length2->Fit(func7, "r n");
  hist_projected_length2->Draw("SAME");
  func7->Draw("SAME");

   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   cout<<""<<endl;

  cout<<"There were "<<neutronnumber<<" neutrons"<<endl;
  cout<<"Total Exposure time now "<<totalexposure/(60*60*24)<<" days"<<endl;
  cout<<"There were ~  "<<neutronnumber/(totalexposure/(60*60*24))<<" neutrons per day"<<endl;

  cout<<"Projected length fit parameter 0 : "<<func4->GetParameter(0)<<endl;
  cout<<"Projected length fit parameter 1 : "<<func4->GetParameter(1)<<endl;
  cout<<"Projected length fit parameter 2 : "<<func4->GetParameter(2)<<endl;
  cout<<"Length polynomial fit parameter 0 : "<<func5->GetParameter(0)<<endl;
  cout<<"Length polynomial fit parameter 1 : "<<func5->GetParameter(1)<<endl;
  cout<<"Length polynomial fit parameter 2 : "<<func5->GetParameter(2)<<endl;

  cout<<" "<<endl;


  cout<<"Projected energy fit parameter 0 : "<<func7->GetParameter(0)<<endl;
  cout<<"Projected energy fit parameter 1 : "<<func7->GetParameter(1)<<endl;
  cout<<"Projected energy fit parameter 2 : "<<func7->GetParameter(2)<<endl;
  cout<<"Reconstruct energy parameter 0 : "<<func6->GetParameter(0)<<endl;
  cout<<"Reconstruct energy parameter 1 : "<<func6->GetParameter(1)<<endl;
  cout<<"Reconstruct energy parameter 2 : "<<func6->GetParameter(2)<<endl;

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
  
