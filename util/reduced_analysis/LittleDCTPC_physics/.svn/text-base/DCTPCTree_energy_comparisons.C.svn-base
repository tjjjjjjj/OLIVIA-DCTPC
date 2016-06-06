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

  //Comparing energy plots
  TH1D *hist_delta_anodeccd=new TH1D("Anode Energy Minus CCD Energy","Anode Energy Minus CCD Energy",50,-1.,1.);
  TH1D *hist_delta_meshccd=new TH1D("Mesh Energy Minus CCD Energy", "Mesh Energy Minus CCD Energy", 10, -1., 1.);
  TH1D *hist_delta_anodemesh=new TH1D("Anode Energy Minus Mesh Energy", "Anode Energy Minus Mesh Energy", 10, -1., 1.);
  TH1D *hist_delta_meshlength=new TH1D("Length Energy Minus Mesh Energy", "Length Energy Minus Mesh Energy", 10, -1., 1.);
  TH1D *hist_delta_anodelength=new TH1D("Anode Energy Minus Length Energy", "Anode Energy Minus Mesh Energy", 10, -1., 1.);
  TH1D *hist_delta_ccdlength=new TH1D("CCD Energy Minus Length Energy", "CCD Energy Minus Length Energy", 10, -1., 1.);
  
  TH2D *hist_energy_anodeccd=new TH2D("Anode Energy vs CCD Energy", "Anode Energy vs CCD Energy", 200, 0, 5000., 200, 0, 5000.);
  TH2D *hist_energy_meshccd = new TH2D ("Mesh Energy vs CCD Energy", "Mesh Energy vs CCD Energy", 200, 0, 10000., 200, 0, 10000.);
  TH2D *hist_energy_anodemesh = new TH2D ("Mesh Energy vs Anode Energy", "Mesh Energy vs Anode Energy", 200, 0, 5000, 100, 0, 2000);;
  TH2D *hist_energy_meshlength=new TH2D("Length Energy vs Mesh Energy", "Length Energy vs Mesh Energy", 200, 0, 5000, 200, 0, 5000);
  TH2D *hist_energy_anodelength=new TH2D("Length Energy vs Anode Energy", "Length Energy vs Anode Energy", 200, 0, 5000, 200, 0, 5000);
  TH2D *hist_energy_ccdlength=new TH2D("Length Energy vs CCD Energy", "Length Energy vs CCD Energy", 200, 0, 5000, 200, 0, 5000);
  
  
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
	{
	  cutnum=14;
	  continue;
	}
    */
      
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	  
      hist_delta_anodeccd->Fill(delta_anodeccd);
      hist_delta_anodemesh->Fill(delta_anodemesh);
      hist_delta_meshccd->Fill(delta_meshccd);
      hist_delta_meshlength->Fill(delta_meshlength);
      hist_delta_ccdlength->Fill(delta_ccdlength);
      hist_delta_anodelength->Fill(delta_anodelength);
      
      hist_energy_anodeccd->Fill(ccdenergy ,aStep.Etrig_kev*anodecalib);
      hist_energy_anodemesh->Fill(aStep.Emesh_kev*meshcalib, aStep.Etrig_kev*anodecalib);
      hist_energy_meshccd->Fill(ccdenergy, aStep.Emesh_kev*meshcalib); 
      hist_energy_meshlength->Fill(lengthenergy, aStep.Emesh_kev*meshcalib);
      hist_energy_anodelength->Fill(lengthenergy, aStep.Etrig_kev*anodecalib);
      hist_energy_ccdlength->Fill(lengthenergy, ccdenergy);
      
      //above I am filling the histograms with data that has passed all the cuts - hopefully neutrons!
      
      cutnum=0;
      neutronnumber++;

    }
  
  //scale after event loop, scale for each sequence and then scale for overall using 1/overall exposure
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  Double_t fitf(Double_t *x, Double_t *par)
  {
    Double_t fitval = par[0]*x[0];
    return fitval;
  }   
  
  TF1 *func1 = new TF1("pol1", "pol1", 200, 1100);//anode vs ccd
  
  TF1 *func2 = new TF1("pol2", "pol2", 100, 7000);//mesh vs ccd

  TF1 *func3 = new TF1("pol1", "pol1", 250, 1200);//anode vs mesh

  TF1 *func4 = new TF1("pol2", "pol2", 300, 4000);//mesh vs length

  TF1 *func5 = new TF1("pol1", "pol1", 300, 900);//anode vs length

  TF1 *func6 = new TF1("pol2", "pol2", 300, 5000);//ccd vs length

  TF1 *func7 = new TF1("pol1", "pol1", 300, 4000);//cut length energy

  TF1 *func8 = new TF1("expo", "expo", 300, 4410);
  
  //above func 1/2/3 is where I define line of best fit for 2D energy plots
  
  TF1 *gaus1 = new TF1("gaus", "gaus", -1, 1);
  
  TF1 *gaus2 = new TF1("gaus", "gaus", -1, 1);
  
  TF1 *gaus3 = new TF1("gaus", "gaus", -1, 1);

  TF1 *gaus4 = new TF1("gaus", "gaus", -1, 1);

  TF1 *gaus5 = new TF1("gaus", "gaus", -1, 1);
  
  TF1 *gaus6 = new TF1("gaus", "gaus", -1, 1);

  TF1 *gaus7 = new TF1("gaus", "gaus", 4410, 6000);
  
  //above gaus 1/2/3 is where I define the fits for the delta distributions
 
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
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
  hist_delta_meshccd->SetXTitle("Calibrated Mesh Energy/keV Minus Calibrated CCD Energy/keV");
  hist_delta_meshccd->SetYTitle("Counts");
  hist_delta_meshccd->Draw();
  hist_delta_meshccd->Fit(gaus2, "r n");
  gaus2->SetLineColor(kBlue);
  gaus2->Draw("SAME");

  new TCanvas;
  hist_delta_anodemesh->SetTitle("Anode Energy Minus Mesh Energy");
  hist_delta_anodemesh->SetXTitle("Calibrated Anode Energy/keV Minus Calibrated Mesh Energy/keV");
  hist_delta_anodemesh->SetYTitle("Counts");
  hist_delta_anodemesh->Draw();
  hist_delta_anodemesh->Fit(gaus3, "r n");
  gaus3->SetLineColor(kBlue);
  gaus3->Draw("SAME");

  new TCanvas;
  hist_delta_meshlength->SetTitle("Length Energy Minus Mesh Energy");
  hist_delta_meshlength->SetXTitle("Reconstructed Length Energy/keV Minus Calibrated Mesh Energy/keV");
  hist_delta_meshlength->SetYTitle("Counts");
  hist_delta_meshlength->Draw();
  hist_delta_meshlength->Fit(gaus4, "r n");
  gaus4->SetLineColor(kBlue);
  gaus4->Draw("SAME");

  new TCanvas;
  hist_delta_anodelength->SetTitle("Anode Energy Minus Length Energy");
  hist_delta_anodelength->SetXTitle("Calibrated Anode Energy/keV Minus Reconstructed Length Energy/keV");
  hist_delta_anodelength->SetYTitle("Counts");
  hist_delta_anodelength->Draw();
  hist_delta_anodelength->Fit(gaus5, "r n");
  gaus5->SetLineColor(kBlue);
  gaus5->Draw("SAME");

  new TCanvas;
  hist_delta_ccdlength->SetTitle("CCD Energy Minus Length Energy");
  hist_delta_ccdlength->SetXTitle("Calibrated CCD Energy/keV Minus Reconstructed Length Energy/keV");
  hist_delta_ccdlength->SetYTitle("Counts");
  hist_delta_ccdlength->Draw();
  hist_delta_ccdlength->Fit(gaus6, "r n");
  gaus6->SetLineColor(kBlue);
  gaus6->Draw("SAME");
  
  new TCanvas;
  hist_energy_anodeccd->SetTitle("Anode Energy vs CCD Energy");
  hist_energy_anodeccd->SetXTitle("Calibrated CCD Energy/keV");
  hist_energy_anodeccd->SetYTitle("Calibrated Anode Energy/keV");
  hist_energy_anodeccd->SetMarkerStyle(6);
  hist_energy_anodeccd->Draw("COLZ");
  hist_energy_anodeccd->Fit(func1, "r n");
  func1->SetLineColor(kRed);
  func1->Draw("SAME");

  new TCanvas;
  hist_energy_meshccd->SetTitle("Mesh Energy vs CCD Energy");
  hist_energy_meshccd->SetXTitle("Calibrated CCD Energy/keV");
  hist_energy_meshccd->SetYTitle("Calibrated Mesh Energy/keV");
  hist_energy_meshccd->SetMarkerStyle(6);
  hist_energy_meshccd->Draw("COLZ");
  hist_energy_meshccd->Fit(func2, "r n");
  func2->SetLineColor(kRed);
  func2->Draw("SAME");

  new TCanvas;
  hist_energy_anodemesh->SetTitle("Anode Energy vs Mesh Energy");
  hist_energy_anodemesh->SetXTitle("Calibrated Mesh Energy/keV");
  hist_energy_anodemesh->SetYTitle("Calibrated Anode Energy/keV");
  hist_energy_anodemesh->SetMarkerStyle(6);
  hist_energy_anodemesh->Draw("COLZ");
  hist_energy_anodemesh->Fit(func3, "r n");
  func3->SetLineColor(kRed);
  func3->Draw("SAME");

  new TCanvas;
  hist_energy_meshlength->SetTitle("Mesh Energy vs Length Energy");
  hist_energy_meshlength->SetYTitle("Mesh Energy / keV");
  hist_energy_meshlength->SetXTitle("Length Energy / keV");
  hist_energy_meshlength->Draw("COLZ");
  hist_energy_meshlength->Fit(func4, "r n");
  func4->SetLineColor(kRed);
  func4->Draw("SAME");

  new TCanvas;
  hist_energy_anodelength->SetTitle("Anode Energy vs Length Energy");
  hist_energy_anodelength->SetYTitle("Anode Energy / keV");
  hist_energy_anodelength->SetXTitle("Length Energy / keV");
  hist_energy_anodelength->Draw("COLZ");
  hist_energy_anodelength->Fit(func5, " r n");
  func5->SetLineColor(kRed);
  func5->Draw("SAME");

  new TCanvas;
  hist_energy_ccdlength->SetTitle("CCD Energy vs Length Energy");
  hist_energy_ccdlength->SetYTitle("CCD Energy / keV");
  hist_energy_ccdlength->SetXTitle("Length Energy / keV");
  hist_energy_ccdlength->Draw("COLZ");
  hist_energy_ccdlength->Fit(func6, "r n");
  func6->SetLineColor(kRed);
  func6->Draw("SAME");
 
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  cout<<"Mean (Anode Energy - CCD Energy)= "<<gaus1->GetParameter(1)<<endl;
  cout<<"Anode Energy / CCD Energy= "<<func1->GetParameter(1)<<endl;
  cout<<"Mean (Mesh Energy - CCD Energy)= "<<gaus2->GetParameter(1)<<endl;
  cout<<"Mesh Energy / CCD Energy= "<<func2->GetParameter(1)<<endl;
  cout<<"Mean (Anode Energy - Mesh Energy)= "<<gaus3->GetParameter(1)<<endl;
  cout<<"Anode Energy / Mesh Energy= "<<func3->GetParameter(1)<<endl;
  cout<<"Mean (Length Energy - Mesh Energy)= "<<gaus4->GetParameter(1)<<endl;
  cout<<"Mesh Energy / Length Energy= "<<func4->GetParameter(1)<<endl;
  cout<<"Mean (Anode Energy - Length Energy)= "<<gaus5->GetParameter(1)<<endl;
  cout<<"Anode Energy / Length Energy= "<<func5->GetParameter(1)<<endl;
  cout<<"Mean (CCD Energy - Length Energy)= "<<gaus6->GetParameter(1)<<endl;
  cout<<"CCD Energy / Length Energy= "<<func6->GetParameter(1)<<endl;

  cout<<""<<endl;

  cout<<"There were "<<neutronnumber<<" neutrons"<<endl;
  cout<<"Total Exposure time now "<<totalexposure/(60*60*24)<<" days"<<endl;
  cout<<"There were ~  "<<neutronnumber/(totalexposure/(60*60*24))<<" neutrons per day"<<endl;

}

double fitPoly12Percent(double E)
{
  return sqrt(0.9741*0.9741*E*E+113.835*113.835)-62.675;
  //return -11.50+39.02*E+0.03459*E*E;                                                    
  //return 469 + 10*E + 0.548 *E*E -3.39e-3*E*E*E+7.5e-6*E*E*E*E;          
  
}               

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
