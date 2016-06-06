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
  
  ofstream neutrons;
  neutrons.open("Neutrons.txt");
  ofstream cuts;
  cuts.open("Cuts.txt");


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //Energy plots
  TH1D *hist_energy_trig=new TH1D("Anode Energy","Anode Energy",2000,0,2000); 
  TH1D *hist_energy_track=new TH1D("CCD Energy", "CCD Energy", 50, 0, 10000);//this is going so high because i still need to cut the double alpha events which only the CCD picks up
  TH1D *hist_energy_mesh=new TH1D("Mesh Energy", "Mesh Energy", 50, 0, 10000);
  TH1D *hist_energy_trig_log=new TH1D ("Anode Energy Log Scale", "Anode Energy Log Scale", 2000, 0, 2000);
  TH1D *hist_energy_track_log=new TH1D("CCD Energy Log Scale", "CCD Energy Log Scale", 50, 0, 10000);
  TH1D *hist_energy_mesh_log=new TH1D("Mesh Energy Log Scale", "Mesh Energy Log Scale", 50, 0, 10000);

  //Comparing energy plots
  TH1D *hist_delta_anodeccd=new TH1D("Anode Energy Minus CCD Energy","Anode Energy Minus CCD Energy",50,-1.,1.);
  TH1D *hist_delta_meshccd=new TH1D("Mesh Energy Minus CCD Energy", "Mesh Energy Minus CCD Energy", 10, -1., 1.);
  TH1D *hist_delta_anodemesh=new TH1D("Anode Energy Minus Mesh Energy", "Anode Energy Minus Mesh Energy", 10, -1., 1.);
  TH2D *hist_mesh_over_ccd=new TH2D("Mesh Energy/CCD Energy vs CCD Energy", "Mesh Energy/CCD Energy vs CCD Energy", 300, 0, 3, 2000, 0, 10000);
  TH2D *hist_energy_anodeccd=new TH2D("Anode Energy vs CCD Energy", "Anode Energy vs CCD Energy", 100, 0, 2000., 100, 0, 2000.);
  TH2D *hist_energy_meshccd = new TH2D ("Mesh Energy vs CCD Energy", "Mesh Energy vs CCD Energy", 100, 0, 10000., 100, 0, 10000.);
  TH2D *hist_energy_anodemesh = new TH2D ("Mesh Energy vs Anode Energy", "Mesh Energy vs Anode Energy", 100, 0, 5000, 100, 0, 2000);

  //Rise and Fall Time Plots
   TH1D *hist_meshrise_0=new TH1D("R0", "R0", 100, 0, 0.00001);
  TH1D *hist_meshrise_10=new TH1D("R10", "R10", 100, 0, 0.00001);
  TH1D *hist_meshrise_25=new TH1D("R25", "R25", 100, 0, 0.00001);
  TH1D *hist_meshrise_50=new TH1D("R50", "R50", 100, 0, 0.00001);
  TH1D *hist_meshrise_75=new TH1D("R75", "R75", 100, 0, 0.00001);
  TH1D *hist_meshrise_90=new TH1D("R90", "R90", 100, 0, 0.00001);
  TH1D *hist_meshfall_0=new TH1D("F0", "F0", 100, 0, 0.00001);
  TH1D *hist_meshfall_10=new TH1D("F10", "F10", 100, 0, 0.00001);
  TH1D *hist_meshfall_25=new TH1D("F25", "F25", 100, 0, 0.00001);
  TH1D *hist_meshfall_50=new TH1D("F50", "F50", 100, 0, 0.00001);
  TH1D *hist_meshfall_75=new TH1D("F75", "F75", 100, 0, 0.00001);
  TH1D *hist_meshfall_90=new TH1D("F90", "F90", 100, 0, 0.00001);
  
  //Miscellaneous Plots
  TH1D *hist_maxpixel=new TH1D("Max Pixel", "Max Pixel", 100, 0, 1000);
  TH1D *hist_trackangle = new TH1D ("Track Angle", "Track Angle", 360, -180., 180.);
  TH1D *hist_cutnum=new TH1D("Cut Number","Cut Number",20,0,20);
  TH2D *hist_pos=new TH2D("Track Position","Track Position",100,-53.4,220,100,-53.4,220);
  TH2D *hist_startpos = new TH2D("Track Start Position", "Track Start Position", 100, -53.4., 220., 100, -53.4., 220.);
  TH2D *hist_endpos = new TH2D("Track End Position", "Track End Position", 100, -53.4., 220., 100, -53.4., 220.);
  TH1D *hist_rangeccd=new TH1D("CCD Range", "CCD Range", 200, 0, 200); //do I want to fit a Gaussian to this
  TH1D *hist_rangeccdlog=new TH1D("CCD Range Log Scale", "CCD Range Log Scale", 200, 0, 200);
  TH1D *hist_rho_start=new TH1D("Rho Start", "Rho Start", 50, 0, 150);
  TH1D *hist_rho_end=new TH1D("Rho End", "Rho End", 50, 0, 150);
  TH2D *hist_rhostart_vs_rhoend=new TH2D("Rho Start vs Rho End", "Rho Start vs Rho End", 150, 0, 150, 150, 0, 150);
  TH1D *hist_rho_start_sq= new TH1D("Rho Start Squared", "Rho Start Squared", 250, 0, 22500);
  TH1D *hist_rho_end_sq= new TH1D("Rho End Squared", "Rho End Squared", 250, 0, 22500);
  TH1D *hist_skewness= new TH1D("Track Skewness", "Track Skewness", 100, 0, 1.5);
  TH1D *hist_tracklength= new TH1D("Track Length", "Track Length", 150, 0, 300);
  TH1D *hist_meshpeak= new TH1D("Mesh Peak Value/Emesh", "Mesh Peak Value/Emesh", 200, 0, 10e-04);

  //Length and energy plots
  TH2D *hist_range_vs_ccd=new TH2D ("Range vs CCD Energy", "Range vs CCD Energy", 200, 0, 10000, 150, 0, 150);
  TH2D *hist_range_vs_mesh=new TH2D ("Range vs Mesh Energy", "Range vs Mesh Energy", 200, 0, 10000, 150, 0, 150);
  TH2D *hist_range_vs_anode = new TH2D ("Range vs Anode Energy", "Range vs Anode Energy", 200, 0, 2000, 150, 0, 150);
  TH2D *hist_length_vs_ccd=new TH2D ("Length vs CCD Energy", "Length vs CCD Energy", 200, 0, 10000, 150, 0, 300);
  TH2D *hist_length_vs_mesh=new TH2D ("Length vs Mesh Energy", "Length vs Mesh Energy", 200, 0, 10000, 150, 0, 300);
  TH2D *hist_length_vs_anode=new TH2D ("Length vs Anode Energy", "Length vs Anode Energy", 200, 0, 2000, 150, 0, 300);
  TH2D *hist_x_length_vs_y_length=new TH2D ("Horizontal vs Vertical Length", "Horizontal vs Vertical Length", 400, 0, 400, 400, 0, 400);

  //Neutron number plots
  TH2D *hist_neutronsperrun=new TH2D ("Number of Neutrons per Run", "Neutrons per Run", 8000, 9000, 17000, 1400, 0, 1400);
  
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
      if(event%10000==0)
      cout<<((double)event/(double)nentries)*100.<<"%"<<endl;
    
      hist_cutnum->Fill(cutnum);
      cuts<<"Cut: "<<aStep.RunNum<<" "<<aStep.EventNum<<endl;
  
      double r2 = (aStep.Track_x_pix*aStep.Track_x_pix+aStep.Track_y_pix*aStep.Track_y_pix)*pow(MMPERPIXEL/10.,2);
      double weight=1./pow(RADIAL_PARAM0_CALIB/(RADIAL_PARAM1_CALIB+r2),RADIAL_PARAM2_CALIB);//radial calibration 
      double rho_start=(sqrt(((aStep.Track_x_start_pix*MMPERPIXEL)*(aStep.Track_x_start_pix*MMPERPIXEL))+((aStep.Track_y_start_pix*MMPERPIXEL)*(aStep.Track_y_start_pix*MMPERPIXEL))));
      double rho_end=(sqrt(((aStep.Track_x_end_pix*MMPERPIXEL)*(aStep.Track_x_end_pix*MMPERPIXEL))+((aStep.Track_y_end_pix*MMPERPIXEL)*(aStep.Track_y_end_pix*MMPERPIXEL))));
      double verticallength=((aStep.Mesh_R0time_samp-aStep.Mesh_R90time_samp)*driftspeed);
      double horizontallength=aStep.Track_range_pix*MMPERPIXEL;
      double tracklength=(sqrt((verticallength*verticallength)+(horizontallength*horizontallength)));
      double meshcalib=MESH_CALIB*seqcalib[aStep.SequenceNum-2];
      double ccdcalib=CCD_CALIB*seqcalib[aStep.SequenceNum-2];
      double anodecalib=ANODE_CALIB*seqcalib[aStep.SequenceNum-2];
      double delta_anodeccd=((aStep.Etrig_kev*anodecalib)-(aStep.Etrack_kev*weight*ccdcalib))/((aStep.Etrig_kev*anodecalib)+(aStep.Etrack_kev*weight*ccdcalib));
      double delta_meshccd=((aStep.Emesh_kev*meshcalib)-(aStep.Etrack_kev*weight*ccdcalib))/((aStep.Emesh_kev*meshcalib)+(aStep.Etrack_kev*weight*ccdcalib));
      double delta_anodemesh=((aStep.Etrig_kev*anodecalib)-(aStep.Emesh_kev*meshcalib))/((aStep.Etrig_kev*anodecalib)+(aStep.Emesh_kev*meshcalib));
      double mesh_over_ccd=((aStep.Emesh_kev*meshcalib)/(aStep.Etrack_kev*weight*ccdcalib));
      
      
      //When nudsk0001 is back up implement this cout to check through to see if we can cut events with anode less than 140 keV: if(aStep.Etrig_kev*0.958675<140){cout<<"Check I am just noise: "<<aStep.RunNum<<" "<<aStep.EventNum<<endl;
      //When nudsk0001 is back up also implement this cout to check the match of mesh and CCD events if((aStep.Emesh_kev*0.645707-(aStep.Etrack_kev*weight*1.0528)/(aStep.Emesh_kev*0.645707+(aStep.Etrack_kev*weight*1.0528)<0.4{cout<<"Check my mesh and CCD match! <<"aStep.RunNum<<" "<<aStep.EventNum<<endl;
      //When nudsk0001 is back up also implement this cout to check the match of anode and CCD events if(((aStep.Etirg_kev*0.958675)-(aStep.Etrack_kev*weight*1.0528)/((aStep.Etrig_kev*0.958675)+(aStep.Etrack_kev*weight*1.0528)<0.4{cout<<"Check my anode and CCD match! <<"aStep.RunNum<<" "<<aStep.EventNum<<endl;
      
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
	}//check this cut!!

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
      
      if(aStep.Track_maxpixel_ccdadu>200)
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
      
      // if(aStep.Emesh_kev*meshcalib>2000/3000/4000/5000)
      //{
      //cutnum=14;
      //continue;
      //} //apply this cut if you only want to look at neutrons

      // if(aStep.Ntrig!=1)
      //{
      //  cutnum=15;
      //  continue;
      //}
     	 
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	 

	  hist_delta_anodeccd->Fill((aStep.Etrig_kev*(anodecalib)-(aStep.Etrack_kev*weight*ccdcalib))/(aStep.Etrig_kev*(anodecalib)+(aStep.Etrack_kev*weight*ccdcalib)));
	  hist_delta_anodemesh->Fill((aStep.Etrig_kev*anodecalib)-(aStep.Emesh_kev*meshcalib))/((aStep.Etrig_kev*anodecalib)+(aStep.Emesh_kev*meshcalib));
	  hist_delta_meshccd->Fill((aStep.Emesh_kev*meshcalib)-(aStep.Etrack_kev*weight*ccdcalib))/((aStep.Emesh_kev*meshcalib)+(aStep.Etrack_kev*weight*ccdcalib));
	  hist_mesh_over_ccd->Fill(mesh_over_ccd, aStep.Etrack_kev*ccdcalib*weight);

	  hist_energy_anodeccd->Fill(aStep.Etrack_kev*weight*ccdcalib ,aStep.Etrig_kev*anodecalib);
	  hist_energy_anodemesh->Fill(aStep.Emesh_kev*meshcalib, aStep.Etrig_kev*anodecalib);
	  hist_energy_meshccd->Fill(aStep.Etrack_kev*weight*ccdcalib, aStep.Emesh_kev*meshcalib); 
	    	
	  hist_energy_trig->Fill(aStep.Etrig_kev*anodecalib);
	  hist_energy_track->Fill(aStep.Etrack_kev*weight*(ccdcalib));
	  hist_energy_mesh->Fill(aStep.Emesh_kev*meshcalib);
	  hist_energy_trig_log->Fill(aStep.Etrig_kev*anodecalib);
	  hist_energy_track_log->Fill(aStep.Etrack_kev*weight*(ccdcalib));
	  hist_energy_mesh_log->Fill(aStep.Emesh_kev*(meshcalib));
	  
	  hist_pos->Fill((aStep.Track_x_pix+512)*MMPERPIXEL,(aStep.Track_y_pix+512)*MMPERPIXEL);
	  hist_startpos->Fill((aStep.Track_x_start_pix+512)*MMPERPIXEL, (aStep.Track_y_start_pix+512)*MMPERPIXEL);
	  hist_endpos->Fill((aStep.Track_x_end_pix+512)*MMPERPIXEL, (aStep.Track_y_end_pix+512)*MMPERPIXEL);

	  hist_meshrise_0->Fill(aStep.Mesh_R0time_samp*secondspersample);
	  hist_meshrise_10->Fill(aStep.Mesh_R10time_samp*secondspersample);
	  hist_meshrise_25->Fill(aStep.Mesh_R25time_samp*secondspersample);
	  hist_meshrise_50->Fill(aStep.Mesh_R50time_samp*secondspersample);
	  hist_meshrise_75->Fill(aStep.Mesh_R75time_samp*secondspersample);
	  hist_meshrise_90->Fill(aStep.Mesh_R90time_samp*secondspersample);
	  hist_meshfall_0->Fill(aStep.Mesh_F0time_samp*secondspersample);
	  hist_meshfall_10->Fill(aStep.Mesh_F10time_samp*secondspersample);
	  hist_meshfall_25->Fill(aStep.Mesh_F25time_samp*secondspersample);
	  hist_meshfall_50->Fill(aStep.Mesh_F50time_samp*secondspersample);
	  hist_meshfall_75->Fill(aStep.Mesh_F75time_samp*secondspersample);
	  hist_meshfall_90->Fill(aStep.Mesh_F90time_samp*secondspersample);	  
	  
	  hist_rangeccd->Fill(aStep.Track_range_pix*MMPERPIXEL);
	  hist_rangeccdlog->Fill(aStep.Track_range_pix*MMPERPIXEL);
	  hist_maxpixel->Fill(aStep.Track_maxpixel_ccdadu);
	  hist_trackangle->Fill(aStep.Track_phi_deg);
	  hist_rho_start->Fill(rho_start);
	  hist_rho_end->Fill(rho_end);
	  hist_rhostart_vs_rhoend->Fill(rho_start, rho_end);
	  hist_rho_start_sq->Fill(rho_start*rho_start);
	  hist_rho_end_sq->Fill(rho_end*rho_end);
	  hist_skewness->Fill(TMath::Abs(aStep.Track_skewness));
      	  hist_tracklength->Fill(tracklength);
	  hist_meshpeak->Fill(aStep.Mesh_peak_v/aStep.Emesh_kev*meshcalib);
	  
	  hist_range_vs_ccd->Fill(aStep.Etrack_kev*weight*(ccdcalib), aStep.Track_range_pix*MMPERPIXEL);
	  hist_range_vs_mesh->Fill(aStep.Emesh_kev*meshcalib, aStep.Track_range_pix*MMPERPIXEL);
	  hist_range_vs_anode->Fill(aStep.Etrig_kev*anodecalib, aStep.Track_range_pix*MMPERPIXEL);
       	  hist_length_vs_ccd->Fill(aStep.Etrack_kev*weight*ccdcalib, tracklength);
       	  hist_length_vs_mesh->Fill(aStep.Emesh_kev*meshcalib, tracklength);
	  hist_length_vs_anode->Fill(aStep.Etrig_kev*anodecalib, tracklength);
	  hist_x_length_vs_y_length->Fill(horizontallength, verticallength);
	  
	  //above I am filling the histograms with data that has passed all the cuts - hopefully neutrons!
	  
	  neutrons<<"Neutron: "<<aStep.RunNum<<", "<<aStep.EventNum<<endl;
	  
	  cutnum=0;
	  neutronnumber++;

	  hist_neutronsperrun->Fill(aStep.RunNum, neutronnumber);
	 // cout<<"run number "<<aStep.RunNum<<" neutron number "<<neutronnumber<<endl;
	  
	 // if(aStep.Track_range_pix*MMPERPIXEL>70)
	  // {cout<<aStep.RunNum<<" "<<aStep.EventNum<<endl;}

    }
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  Double_t fitf(Double_t *x, Double_t *par)
  {
    Double_t fitval = par[0]*x[0];
    return fitval;
  }   
  
  TF1 *func1 = new TF1("pol1", "pol1", 200, 1200);//anode vs ccd
  
  TF1 *func2 = new TF1("pol1", "pol1", 0, 3000);//mesh vs ccd

  TF1 *func3 = new TF1("pol1", "pol1", 250, 1200);//anode vs mesh

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
  hist_neutronsperrun->SetTitle("Number of Neutrons per Run");
  hist_neutronsperrun->SetXTitle("Run Number");
  hist_neutronsperrun->SetYTitle("Total Number of Neutrons");
  hist_neutronsperrun->SetMarkerStyle(6);
  hist_neutronsperrun->Draw();//it looks weird because the for loop runs over events from 10000 first then goes back to 9000s

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

  new TCanvas;
  hist_range_vs_anode->SetTitle("Range vs Anode Enegry");
  hist_range_vs_anode->SetXTitle("Calibrated Anode Energy/keV");
  hist_range_vs_anode->SetYTitle("Range/mm");
  hist_range_vs_anode->Draw("COLZ");

  new TCanvas;
  hist_delta_anodeccd->SetTitle("Anode Energy Minus CCD Energy");
  hist_delta_anodeccd->SetXTitle("Calibrated Anode Energy/keV Minus Calibrated CCD Energy/keV");
  hist_delta_anodeccd->SetYTitle("Counts");
  hist_delta_anodeccd->Draw();
  hist_delta_anodeccd->Fit(gaus1, "r n");
  gaus1->SetLineColor(kBlue);
  gaus1->Draw("SAME");

  new TCanvas;
  hist_delta_meshccd->SetTitle("Mesh Energy Minus CCD Energy");
  hist_delta_meshccd->SetXTitle("Mesh Energy/keV Minus Calibrated CCD Energy/keV");
  hist_delta_meshccd->SetYTitle("Counts");
  hist_delta_meshccd->Draw();
  hist_delta_meshccd->Fit(gaus2, "r n");
  gaus2->SetLineColor(kBlue);
  gaus2->Draw("SAME");

  new TCanvas;
  hist_delta_anodemesh->SetTitle("Anode Energy Minus Mesh Energy");
  hist_delta_anodemesh->SetXTitle("Calibrated Anode Energy/keV Minus Mesh Energy/keV");
  hist_delta_anodemesh->SetYTitle("Counts");
  hist_delta_anodemesh->Draw();
  hist_delta_anodemesh->Fit(gaus3, "r n");
  gaus3->SetLineColor(kBlue);
  gaus3->Draw("SAME");

  new TCanvas;
  hist_energy_anodeccd->SetTitle("Anode Energy vs CCD Energy");
  hist_energy_anodeccd->SetXTitle("Calibrated CCD Energy/keV");
  hist_energy_anodeccd->SetYTitle("Calibrated Anode Energy/keV");
  hist_energy_anodeccd->SetMarkerStyle(6);
  hist_energy_anodeccd->Draw("COLZ");
  // hist_energy_anodeccd->Fit(func1, "r n");
  //func1->SetLineColor(kRed);
  //func1->Draw("SAME");

  new TCanvas;
  hist_energy_meshccd->SetTitle("Mesh Energy vs CCD Energy");
  hist_energy_meshccd->SetXTitle("Calibrated CCD Energy/keV");
  hist_energy_meshccd->SetYTitle("Mesh Energy/keV");
  hist_energy_meshccd->SetMarkerStyle(6);
  hist_energy_meshccd->Draw("COLZ");
  //hist_energy_meshccd->Fit(func2, "r n");
  //func2->SetLineColor(kRed);
  //func2->Draw("SAME");

  new TCanvas;
  hist_energy_anodemesh->SetTitle("Anode Energy vs Mesh Energy");
  hist_energy_anodemesh->SetXTitle("Mesh Energy/keV");
  hist_energy_anodemesh->SetYTitle("Calibrated Anode Energy/keV");
  hist_energy_anodemesh->SetMarkerStyle(6);
  hist_energy_anodemesh->Draw("COLZ");
  //hist_energy_anodemesh->Fit(func3, "r n");
  //func3->SetLineColor(kRed);
  //func3->Draw("SAME");

  new TCanvas;
  c1_n11->SetLogy();
  hist_energy_trig_log->SetXTitle("Anode Energy/keV");
  hist_energy_trig_log->SetYTitle("Log Counts");
  hist_energy_trig_log->SetTitle("Calibrated Anode Energy Log Scale");
  hist_energy_trig_log->Draw();

  new TCanvas;
  hist_energy_trig->SetXTitle("Calibrated Anode Energy/keV");
  hist_energy_trig->SetYTitle("Counts");
  hist_energy_trig->SetTitle("Calibrated Anode Energy");
  hist_energy_trig->Draw();
  
  new TCanvas;
  c1_n13->SetLogy();
  hist_energy_track_log->SetXTitle("CCD Energy/keV");
  hist_energy_track_log->SetYTitle("Log Counts");
  hist_energy_track_log->SetTitle("Calibrated CCD Energy Log Scale");
  hist_energy_track_log->Draw();

  new TCanvas;
  hist_energy_track->SetTitle("Calibrated CCD Energy");
  hist_energy_track->SetXTitle("Calibrated CCD Energy/keV");
  hist_energy_track->SetYTitle("Counts");
  hist_energy_track->Draw();

  new TCanvas;
  c1_n15->SetLogy();
  hist_energy_mesh_log->SetXTitle("Mesh Energy/keV");
  hist_energy_mesh_log->SetYTitle("Log Counts");
  hist_energy_mesh_log->SetTitle("Calibrated Mesh Energy Log Scale");
  hist_energy_mesh_log->Draw();

  new TCanvas;
  hist_energy_mesh->SetTitle("Calibrated Mesh Energy");
  hist_energy_mesh->SetXTitle("Calibrated Mesh Energy/keV");
  hist_energy_mesh->SetYTitle("Counts");
  hist_energy_mesh->Draw();

  new TCanvas;
  hist_pos->SetTitle("Track Centre Position");
  hist_pos->SetXTitle("Track X/mm");
  hist_pos->SetYTitle("Track Y/mm");
  hist_pos->SetMarkerStyle(8);
  hist_pos->Draw("COLZ");
  innerveto->SetFillStyle(0);
  innerveto->SetLineColor(2);
  innerveto->Draw("SAME");
  outerveto->SetFillStyle(0);
  outerveto->SetLineColor(2);
  outerveto->Draw("SAME");

   new TCanvas;
  hist_startpos->SetTitle("Track Start Position");
  hist_startpos->SetXTitle("Track X/mm");
  hist_startpos->SetYTitle("Track Y/mm");
  hist_startpos->SetMarkerStyle(8);
  hist_startpos->Draw("COLZ");
  innerveto->SetFillStyle(0);
  innerveto->SetLineColor(2);
  innerveto->Draw("SAME");
  outerveto->SetFillStyle(0);
  outerveto->SetLineColor(2);
  outerveto->Draw("SAME");

   new TCanvas;
  hist_endpos->SetTitle("Track End  Position");
  hist_endpos->SetXTitle("Track X/mm");
  hist_endpos->SetYTitle("Track Y/mm");
  hist_endpos->SetMarkerStyle(8);
  hist_endpos->Draw("COLZ");
  innerveto->SetFillStyle(0);
  innerveto->SetLineColor(2);
  innerveto->Draw("SAME");
  outerveto->SetFillStyle(0);
  outerveto->SetLineColor(2);
  outerveto->Draw("SAME");
  
  new TCanvas;
  hist_rangeccd->SetTitle("CCD Range");
  hist_rangeccd->SetXTitle("CCD Range/mm");
  hist_rangeccd->SetYTitle("Counts");
  hist_rangeccd->Draw();

  new TCanvas;
  c1_n21->SetLogy();
  hist_rangeccdlog->SetTitle("CCD Range Log Scale");
  hist_rangeccdlog->SetXTitle("CCD Range/mm");
  hist_rangeccdlog->SetYTitle("Log Counts");
  hist_rangeccdlog->Draw();
  
  new TCanvas;
  hist_trackangle->SetTitle("Track Phi");
  hist_trackangle->SetXTitle("Track Phi");
  hist_trackangle->SetYTitle("Counts");
  hist_trackangle->Draw();

  new TCanvas;
  hist_maxpixel->SetTitle("Maximum Pixel/adu");
  hist_maxpixel->SetXTitle("Maximum Pixel/adu");
  hist_maxpixel->SetYTitle("Counts");
  hist_maxpixel->Draw();

  new TCanvas;
  hist_rho_start->SetTitle("Rho Start/mm");
  hist_rho_start->SetXTitle("Rho Start/mm");
  hist_rho_start->SetYTitle("Counts");
  hist_rho_start->Draw();

  new TCanvas;
  hist_rho_start_sq->SetTitle("Rho Start Squared/mm");
  hist_rho_start_sq->SetXTitle("Rho Start Squared/mm");
  hist_rho_start_sq->SetYTitle("Counts");
  hist_rho_start_sq->Draw();

  new TCanvas;
  hist_rho_end->SetTitle("Rho End/mm");
  hist_rho_end->SetXTitle("Rho End/mm");
  hist_rho_end->SetYTitle("Counts");
  hist_rho_end->Draw();

  new TCanvas;
  hist_rho_end_sq->SetTitle("Rho End Squared/mm");
  hist_rho_end_sq->SetXTitle("Rho End Squared/mm");
  hist_rho_end_sq->SetYTitle("Counts");
  hist_rho_end_sq->Draw();

  new TCanvas;
  hist_rhostart_vs_rhoend->SetTitle("Rho Start vs Rho End /mm");
  hist_rhostart_vs_rhoend->SetXTitle("Rho Start/mm");
  hist_rhostart_vs_rhoend->SetYTitle("Rho End/mm");
  hist_rhostart_vs_rhoend->Draw("COLZ");

  new TCanvas;
  hist_skewness->SetTitle("Track Skewness");
  hist_skewness->SetXTitle("Track Skewness");
  hist_skewness->SetYTitle("Counts");
  hist_skewness->Draw();
 
  new TCanvas;
  hist_mesh_over_ccd->SetTitle("Mesh Energy keV/CCD Energy keV vs CCD Energy keV");
  hist_mesh_over_ccd->SetXTitle("Mesh Energy/keV / CCD Energy/keV");
  hist_mesh_over_ccd->SetYTitle("CCD Energy/keV");
  hist_mesh_over_ccd->Draw("COLZ");

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
  
  new TCanvas;
  hist_length_vs_mesh->SetTitle("Track Length vs Mesh Energy");
  hist_length_vs_mesh->SetXTitle("Mesh Energy/keV");
  hist_length_vs_mesh->SetYTitle("Track Length/mm");
  hist_length_vs_mesh->Draw("COLZ");
  
  new TCanvas;
  hist_length_vs_anode->SetTitle("Track Length vs Anode Energy");
  hist_length_vs_anode->SetXTitle("Anode Energy/keV");
  hist_length_vs_anode->SetYTitle("Track Length/mm");
  hist_length_vs_anode->Draw("COLZ");
  
  new TCanvas;
  hist_x_length_vs_y_length->SetTitle("Horizontal Length vs Vertical Length");
  hist_x_length_vs_y_length->SetXTitle("Horizontal Length / mm");
  hist_x_length_vs_y_length->SetYTitle("Vertical Length / mm");
  hist_x_length_vs_y_length->Draw("COLZ");

  new TCanvas;
  hist_meshpeak->SetTitle("Mesh Peak Value / Mesh Energy");
  hist_meshpeak->SetXTitle("Mesh Peak Value / Mesh Energy");
  hist_meshpeak->SetYTitle("Counts");
  hist_meshpeak->Draw();
  
  new TCanvas;
  hist_meshrise_0->SetTitle("Mesh Rise Times / sec");
  hist_meshrise_0->SetYTitle("Counts");
  hist_meshrise_0->SetXTitle("Mesh Rise Time / sec");
  hist_meshrise_0->Draw();
  hist_meshrise_10->SetLineColor(kRed);
  hist_meshrise_10->Draw("SAME");
  hist_meshrise_25->SetLineColor(kBlue);
  hist_meshrise_25->Draw("SAME");
  hist_meshrise_50->SetLineColor(kGreen);
  hist_meshrise_50->Draw("SAME");
  hist_meshrise_75->SetLineColor(kYellow);
  hist_meshrise_75->Draw("SAME");
  hist_meshrise_90->SetLineColor(kPink);
  hist_meshrise_90->Draw("SAME");

  new TCanvas;
  hist_meshfall_0->SetTitle("Mesh Fall Times / sec");
  hist_meshfall_0->SetYTitle("Counts");
  hist_meshfall_0->SetXTitle("Mesh Fall Time / sec");
  hist_meshfall_0->Draw();
  hist_meshfall_10->SetLineColor(kRed);
  hist_meshfall_10->Draw("SAME");
  hist_meshfall_25->SetLineColor(kBlue);
  hist_meshfall_25->Draw("SAME");
  hist_meshfall_50->SetLineColor(kGreen);
  hist_meshfall_50->Draw("SAME");
  hist_meshfall_75->SetLineColor(kYellow);
  hist_meshfall_75->Draw("SAME");
  hist_meshfall_90->SetLineColor(kPink);
  hist_meshfall_90->Draw("SAME");
  
  new TCanvas;
  hist_cutnum->SetTitle("Cut Number");
  hist_cutnum->SetXTitle("Cut Number");
  hist_cutnum->SetYTitle("Counts");
  hist_cutnum->Draw();

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  cout<<"Number of events cut= "<<nentries-neutronnumber<<endl;
  cout<<"Most used cut is cut number = "<<hist_cutnum->GetMaximumBin()<<endl;
  cout<<"Least used cut is cut number = "<<hist_cutnum->GetMinimumBin()<<endl;

  cout<<""<<endl;

  cout <<"Mean CCD Range= "<< hist_rangeccd->GetMean()<<endl;
  cout<<"Mean Track Length= "<<hist_tracklength->GetMean()<<endl;
  cout<<"Mean Maximum Pixel adu = "<<hist_maxpixel->GetMean()<<endl;

  cout<<""<<endl;

  cout<<"Mean Anode Energy= "<<hist_energy_trig->GetMean()<<endl;
  cout<<"Modal Anode Energy= "<<hist_energy_trig->GetMaximumBin()<<endl;
  cout<<"Mean CCD Energy = "<<hist_energy_track->GetMean()<<endl;
  cout<<"Modal CCD Energy= "<<hist_energy_trig->GetMaximumBin()<<endl;
  cout<<"Mean Mesh Energy= "<<hist_energy_mesh->GetMean()<<endl;
  cout<<"Modal Mesh Energy= "<<hist_energy_mesh->GetMaximumBin()<<endl;

  cout<<""<<endl;

  cout<<"Mean (Anode Energy - CCD Energy)= "<<gaus1->GetParameter(1)<<endl;
  cout<<"Anode Energy / CCD Energy= "<<func1->GetParameter(0)<<endl;
  cout<<"Mean (Mesh Energy - CCD Energy)= "<<gaus2->GetParameter(1)<<endl;
  cout<<"Mesh Energy / CCD Energy= "<<func2->GetParameter(0)<<endl;
  cout<<"Mean (Anode Energy - Mesh Energy)= "<<gaus3->GetParameter(1)<<endl;
  cout<<"Anode Energy / Mesh Energy= "<<func3->GetParameter(0)<<endl;

  cout<<""<<endl;

  cout<<"There were "<<neutronnumber<<" neutrons"<<endl;
  cout<<"Total Exposure time now "<<totalexposure/(60*60*24)<<" days"<<endl;
  cout<<"There were ~  "<<neutronnumber/(totalexposure/(60*60*24))<<" neutrons per day"<<endl;

  cuts.close();
  neutrons.close();
  
}

double fitPoly12Percent(double E)
{
  return sqrt(0.9741*0.9741*E*E+113.835*113.835)-62.675;
  //return -11.50+39.02*E+0.03459*E*E;                                                    
  //return 469 + 10*E + 0.548 *E*E -3.39e-3*E*E*E+7.5e-6*E*E*E*E;          
  
}               

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
