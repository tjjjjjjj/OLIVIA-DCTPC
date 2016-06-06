#define DCTPCTree_cxx
#include "../input_files/DCTPCTree.h"
#include "../input_files/Calibration_constants_LittleDCTPC.h";
#include <TH2.h>
#include <TH1.h>
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
  
  TFile *outtree = new TFile("$LittleDCTPC_physics_input");
  TTree *dctreepc = (TTree*)outtree->Get("dctpc_eventinfo");
  DCTPCTree aStep(dctreepc);

   TH1D *hist_energy_length=new TH1D("Length Energy", "Length Energy", 125, 0, 12500);
   TH1D *hist_energy_track=new TH1D("CCD Energy", "CCD Energy", 125, 0, 12500);
  
//create event loop
  Long64_t nentries = dctreepc->GetEntries();
  cout << "Number of Entries: " << nentries << endl;

  int passccd=0;
  int passwfandccd=0;

  for (int event = 0; event<nentries; event++)
    {
      aStep.GetEntry(event);
      if(event%1000000==0)
      cout<<((double)event/(double)nentries)*100.<<"%"<<endl;

      r2 = (aStep.Track_x_pix*aStep.Track_x_pix+aStep.Track_y_pix*aStep.Track_y_pix)*pow(MMPERPIXEL/10.,2);
      weight=1./pow(RADIAL_PARAM0_CALIB/(RADIAL_PARAM1_CALIB+r2),RADIAL_PARAM2_CALIB);//radial calibration 
      rho_start=(sqrt(((aStep.Track_x_start_pix*MMPERPIXEL)*(aStep.Track_x_start_pix*MMPERPIXEL))+((aStep.Track_y_start_pix*MMPERPIXEL)*(aStep.Track_y_start_pix*MMPERPIXEL))));
      rho_end=(sqrt(((aStep.Track_x_end_pix*MMPERPIXEL)*(aStep.Track_x_end_pix*MMPERPIXEL))+((aStep.Track_y_end_pix*MMPERPIXEL)*(aStep.Track_y_end_pix*MMPERPIXEL))));
      verticallength=((aStep.Mesh_R0time_samp+aStep.Mesh_F10time_samp)*driftspeed);
      horizontallength=aStep.Track_range_pix*MMPERPIXEL;
      tracklength=(sqrt((verticallength*verticallength)+(horizontallength*horizontallength)));
      lengthenergy=El_0 + (El_1*tracklength) + (El_2*tracklength*tracklength);

      meshcalib=MESH_CALIB;
      ccdcalib=CCD_CALIB;
      anodecalib=ANODE_CALIB;
      delta_anodeccd=((aStep.Etrig_kev*anodecalib)-(aStep.Etrack_kev*weight*ccdcalib))/((aStep.Etrig_kev*anodecalib)+(aStep.Etrack_kev*weight*ccdcalib));
      delta_meshccd=((aStep.Emesh_kev*meshcalib)-(aStep.Etrack_kev*weight*ccdcalib))/((aStep.Emesh_kev*meshcalib)+(aStep.Etrack_kev*weight*ccdcalib));
      delta_anodemesh=((aStep.Etrig_kev*anodecalib)-(aStep.Emesh_kev*meshcalib))/((aStep.Etrig_kev*anodecalib)+(aStep.Emesh_kev*meshcalib));
      delta_meshlength=((lengthenergy)-(aStep.Emesh_kev*meshcalib))/((lengthenergy)+(aStep.Emesh_kev*meshcalib));
      delta_anodelength=((aStep.Etrig_kev*anodecalib)-(lengthenergy))/((aStep.Etrig_kev*anodecalib)+(lengthenergy));
      delta_ccdlength=((aStep.Etrack_kev*weight*ccdcalib)-(lengthenergy))/((aStep.Etrack_kev*weight*ccdcalib)+(lengthenergy));
      mesh_over_ccd=((aStep.Emesh_kev*meshcalib)/(aStep.Etrack_kev*weight*ccdcalib));
      ccdenergy=(Eccd_0)+(Eccd_1*aStep.Etrack_kev*weight)+(Eccd_2*aStep.Etrack_kev*weight*aStep.Etrack_kev*weight);

      if(aStep.Edge==1)
	{
	  continue;
	}
      
      if(aStep.EventNum-aStep.LastSpark<7)
	{
	  continue;
	}
      
      if(aStep.Track_range_pix*MMPERPIXEL<=0.)
	{
	  continue;
	}
      
      if(aStep.Etrig_kev<400)
	{
	  continue;
	}
      
      if(aStep.Emesh_kev<500)
	{
	  continue;
	}
      
      if(aStep.Etrack_kev*weight<300)
	{
	  continue;
	}
      
      if(aStep.Track_maxpixel_ccdadu<72 ||
				     aStep.Track_maxpixel_ccdadu>200)
	{
	  continue;
	}
      
      
      if(rho_start>100 ||
	 rho_end>100)
	{
 	  continue;
	}

        if(tracklength<72)//change to greater than 72 for neutrons
	{
	  continue;
	}
      

      //anode saturates at Etrig_keV=1280 

       hist_energy_length->Fill(lengthenergy);
	
       hist_energy_track->Fill(ccdenergy);

        
    }

  Double_t fitf(Double_t *x, Double_t *par)
  {
    Double_t fitval = par[0]*x[0];
    return fitval;
  }   

  // TF1 *func = new TF1("fitf",fitf, 400, 1900,1);
  //func->SetParameter(0,1);
  //hist_emesh->Fit(func, "r n");
  TF1 *func1 = new TF1("pol1", "pol1", 400, 1000);

  TF1 *func2 = new TF1("pol2", "pol2" , 0, 15000);

  TF1 *func3 = new TF1("pol1", "pol1", 600, 1100);

  TF1 *gaus1=new TF1("gaus", "gaus", -1, 1);
 
  TF1 *gaus2=new TF1("gaus", "gaus", -1, 1);

  TF1 *gaus3=new TF1("gaus", "gaus", 2000, 8000);

  TF1 *gaus4=new TF1("gaus", "gaus", 1000, 9000);

 
  new TCanvas;
  hist_energy_length->SetTitle("Calibration to Length Energy Peak");
  hist_energy_length->SetYTitle("Counts");
  hist_energy_length->SetXTitle("Energy / keV");
  hist_energy_length->Draw();
  hist_energy_length->Fit(gaus3, "r n");
  gaus3->SetLineColor(kRed);
  gaus3->Draw("SAME");
  hist_energy_track->Draw("SAME");
  hist_energy_track->Fit(gaus4, "r n");
  gaus4->SetLineColor(kBlue);
  gaus4->Draw("SAME");
  leg_energy = new TLegend(0.6, 0.7, 0.89, 0.89);
  leg_energy->AddEntry(gaus3, "Length Energy", "l");
  leg_energy->AddEntry(gaus4, "CCD Energy", "l");
  leg_energy->Draw("SAME");

  cout<<"Mean Length Energy : "<<gaus3->GetParameter(1)<<endl;
  cout<<"Sigma Length Energy : "<<gaus3->GetParameter(2)<<endl;
  cout<<"Length Peak Resolution : "<<((gaus3->GetParameter(2))/(gaus3->GetParameter(1)))*100.<<endl;
  cout<<"Length Energy Calibration Constant : "<<5300/gaus3->GetParameter(1)<<endl;

  cout<<""<<endl;

  cout<<"Mean CCD Energy : "<<gaus4->GetParameter(1)<<endl;
  cout<<"Sigma CCD Energy : "<<gaus4->GetParameter(2)<<endl;
  cout<<"CCD Peak Resolution : "<<((gaus4->GetParameter(2))/(gaus4->GetParameter(1)))*100.<<endl;
  cout<<"CCD Calibration Constant: "<<5300/(gaus4->GetParameter(1))<<endl;
}

double fitPoly12Percent(double E)
{
  return sqrt(0.9741*0.9741*E*E+113.835*113.835)-62.675;
  //return -11.50+39.02*E+0.03459*E*E;                                                    
  //return 469 + 10*E + 0.548 *E*E -3.39e-3*E*E*E+7.5e-6*E*E*E*E;          

}    
