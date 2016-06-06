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


  //Plots to see where the anode and mesh saturate
  TH1D *hist_mesh_max = new TH1D("Mesh Max", "Mesh Max", 100, 0,  0.5);
  TH1D *hist_anode_max=new TH1D("Anode Max", "Anode Max", 50, 0, 0.25);
  TH1D *hist_mesh_peak = new TH1D("Mesh Peak", "Mesh Peak", 100, 0, 0.5);
  TH1D *hist_veto_peak = new TH1D("Veto Peak", "Veto Peak", 100, 0, 0.5);
  TH1D *hist_mesh_energy = new TH1D("Mesh Energy", "Mesh Energy", 1000, 0, 50000);
  TH1D *hist_anode_energy = new TH1D("Anode Energy", "Anode Energy", 200, 0, 2000);
  TH1D *hist_veto_energy=new TH1D("Veto Energy", "Veto Energy", 200, 0, 4000);
  
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
      
      if(rho_start>100 ||
	 rho_end>100)
	{
	  cutnum=8;
	  continue;
	}
      
      
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	       
     

      
      hist_mesh_max->Fill(aStep.Mesh_max_v);
      hist_anode_max->Fill(aStep.Anode_max_v);
      hist_mesh_peak->Fill(aStep.Mesh_peak_v);
      hist_veto_peak->Fill(aStep.Veto_peak_v);

      hist_mesh_energy->Fill(aStep.Mesh_max_v*VOLT_TO_ENERGY);
      hist_anode_energy->Fill(aStep.Anode_max_v*VOLT_TO_ENERGY);
      hist_veto_energy->Fill(aStep.Veto_peak_v*VOLT_TO_ENERGY);
         
      //above I am filling the histograms with data that has passed all the cuts - hopefully neutrons!
         
      cutnum=0;
      neutronnumber++;
  
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
  
  TF1 *func4 = new TF1("pol1", "pol1", 9000, 17000);//neutron rate

  TF1 *func5 = new TF1("pol1", "pol1", 1374000000 ,1386000000);//neutron no. vs time

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
  hist_mesh_max->SetTitle("Mesh Maximum / V");
  hist_mesh_max->SetYTitle("Counts");
  hist_mesh_max->SetXTitle("Mesh Voltage / V");
  hist_mesh_max->Draw();

  new TCanvas;
  hist_anode_max->SetTitle("Anode Maximum / V");
  hist_anode_max->SetYTitle("Counts");
  hist_anode_max->SetXTitle("Anode Voltage / V");
  hist_anode_max->Draw();

  new TCanvas;
  hist_mesh_peak->SetTitle("Mesh Peak / V");
  hist_mesh_peak->SetXTitle("Mesh Voltage / V");
  hist_mesh_peak->SetYTitle("Counts");
  hist_mesh_peak->Draw();

  new TCanvas;
  hist_veto_peak->SetTitle("Veto Peak / V");
  hist_veto_peak->SetXTitle("Veto Voltage / V");
  hist_veto_peak->SetYTitle("Counts");
  hist_veto_peak->Draw();
  
  new TCanvas;
  hist_mesh_energy->SetTitle("Raw Mesh Energy / keV");
  hist_mesh_energy->SetYTitle("Counts");
  hist_mesh_energy->SetXTitle("Mesh Energy / keV");
  hist_mesh_energy->Draw();

  new TCanvas;
  hist_anode_energy->SetTitle("Raw Anode Energy / keV");
  hist_anode_energy->SetYTitle("Counts");
  hist_anode_energy->SetXTitle("Anode Energy / keV");
  hist_anode_energy->Draw();

  new TCanvas;
  hist_veto_energy->SetTitle("Veto Energy / keV");
  hist_veto_energy->SetYTitle("Counts");
  hist_veto_energy->SetXTitle("Veto Energy / keV");
  hist_veto_energy->Draw();

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  cout<<"Number of events cut= "<<nentries-neutronnumber<<endl;
 

  cout<<""<<endl;

  cout<<"Mesh saturates at "<<hist_mesh_max->GetMaximumBin()<<"V / "<<hist_mesh_energy->GetMaximumBin()<<"keV"<<endl;
  cout<<"Anode saturates at "<<hist_anode_max->GetMaximumBin()<<"V / "<<hist_anode_energy->GetMaximumBin()<<"keV"<<endl;
  
  cout<<""<<endl;

  cout<<"There were "<<neutronnumber<<" neutrons"<<endl;
  cout<<"Total Exposure time now "<<totalexposure/(60*60*24)<<" days"<<endl;
  cout<<"There were ~  "<<neutronnumber/(totalexposure/(60*60*24))<<" neutrons per day"<<endl;

  double interestingevents=neutronnumber/83251.;
  cout<<"% of Interesting Events that were Neutrons: "<<interestingevents*100.<<endl;

}

double fitPoly12Percent(double E)
{
  return sqrt(0.9741*0.9741*E*E+113.835*113.835)-62.675;
  //return -11.50+39.02*E+0.03459*E*E;                                                    
  //return 469 + 10*E + 0.548 *E*E -3.39e-3*E*E*E+7.5e-6*E*E*E*E;          
  
}               

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
