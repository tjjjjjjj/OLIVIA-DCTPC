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
  
  ofstream file_passccd;
  file_passccd.open("Pass_CCD.txt");
  ofstream file_passwfandccd;
  file_passwfandccd.open("Pass_WF_and_CCD.txt");
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  double seqcalib [19] = {0.955114 , 1.06241 , 1.07105 , 1.07669 , 1.02827 , 1.15574 , 1.18828 , 1.08779 , 0.998239 , 1.05118 , 1.03128 , 1.02639 , 0.992746 , 0.925486 , 1.1664 , 1.05791, 0.989169};
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //create event loop
  Long64_t nentries = dctreepc->GetEntries();
  cout << "Number of Entries: " << nentries << endl;
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  int cutnum=-1;
  int neutronnumber=0;
  int passall=0;
  int passccd=0;
  int passwf=0;
  int passwfandccd=0;

  TH1D *hist_time_ccd=new TH1D("CCD Pass Rate", "CCD Pass Rate", 100, 1374000000, 1396000000);
  TH1D *hist_time_wf=new TH1D("WF Pass Rate", "WF Pass Rate", 100, 1374000000, 1396000000);
  TH1D *hist_time_wfandccd=new TH1D("CCD and WF Pass Rate", "CCD and WF Pass Rate", 100, 1374000000, 1396000000);

  TH1D *hist_time_ccd_percent = new TH1D("% Events Pass CCD", "% Events Pass CCD", 100, 1374000000, 1396000000);
  TH1D *hist_time_wfandccd_percent = new TH1D("% Events Pass WF and CCD", "% Events Pass WF and CCD", 100, 1374000000, 1396000000);

  TH2D *hist_ccdpass_rate=new TH2D("CCD Pass", "CCD Pass", 100, 1374000000, 1386000000, 2000, 0, 2000);
  TH2D *hist_wfpass_rate=new TH2D("WF Pass", "WF Pass", 100, 1374000000, 1386000000, 2000, 0, 2000);
  TH2D *hist_wfandccdpass_rate=new TH2D("CCD and WF Pass", "CCD and WF Pass", 100, 1374000000, 1386000000, 2000, 0, 2000);

  TH2D *hist_ccdpass_rate_percent=new TH2D("% CCD Pass", "% CCD Pass", 100, 1374000000, 1396000000, 100, 0, 100);
  TH2D *hist_wfandccdpass_rate_percent=new TH2D("% WF and CCD Pass", "% WF and CCD Pass", 100, 1374000000, 1396000000, 100, 0, 100);

  TH2D *hist_pass_rate=new TH2D("Pass", "Pass", 100, 1374000000, 1386000000, 100, 0, 1);

  TH1D *hist_eff_vs_eccd_ccd=new TH1D("Efficiency vs CCD Energy - CCD Only", "Efficiency vs CCD Energy - CCD Only", 200, 0, 6000);
  TH1D *hist_eff_vs_eccd_wfandccd=new TH1D("Efficiency vs CCD Energy - CCD and WF", "Efficiency vs CCD Energy - CCD and WF", 200, 0, 6000);

  TH1D *hist_eff_vs_elength=new TH1D("Efficiency vs Length Energy", "Efficiency vs Length Energy", 200, 0, 6000);
  TH1D *hist_eff_vs_elength_ccd=new TH1D("Efficiency vs Length Energy - CCD Only", "Efficiency vs Length Energy - CCD Only", 200, 0, 6000);
  TH1D *hist_eff_vs_elength_wf=new TH1D("Efficiency vs Length Energy - WF Only", "Efficiency vs Length Energy - WF Only", 200, 0, 6000);
  TH1D *hist_eff_vs_elength_wfandccd=new TH1D("Efficiency vs Length Energy - CCD and WF", "Efficiency vs Length Energy - CCD and WF", 200, 0, 6000);

  /*TH1D *hist_eff_vs_eccd_ccd_percent=new TH1D("%CCD vs Energy", "%CCD vs Energy", 1000, 0, 10000);
  TH1D *hist_eff_vs_eccd_wfandccd_percent = new TH1D("%WF and CCD vs Energy", "%WF and CCD vs Energy", 1000, 0, 10000);
  TH1D *hist_eff_vs_eccd=new TH1D("Efficiency vs CCD Energy", "Efficiency vs CCD Energy", 1000, 0, 10000);*/
  
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
	{continue;}
      
      if(aStep.EventNum-aStep.LastSpark<7)
	{continue;}

      // if(aStep.Etrig_kev<400
      //  || aStep.Emesh_kev<500)
      // {
      //   passwf++;
      //   hist_time_wf->Fill(aStep.Timenow_sec);
      //   hist_wfpass_rate->Fill(aStep.Timenow_sec, passwf);
      //   hist_eff_vs_elength->Fill(lengthenergy);
      // }
     
      
      if(aStep.Track_range_pix*MMPERPIXEL<=0.)
	{continue;}
      
      if(aStep.Etrack_kev*weight<300)
	{continue;}
      
      if(aStep.Track_maxpixel_ccdadu<72 ||
				     aStep.Track_maxpixel_ccdadu>200)
	{continue;}
      
      if(rho_start>100 ||
	 rho_end>100)
	{continue;}
      
      if(tracklength>72)
	{continue;}

      passccd++;
      file_passccd<<aStep.RunNum<<" "<<aStep.EventNum<<endl;
      hist_time_ccd->Fill(aStep.Timenow_sec);
      hist_time_ccd_percent->Fill(aStep.Timenow_sec);
      hist_ccdpass_rate->Fill(aStep.Timenow_sec, passccd);
      hist_ccdpass_rate_percent->Fill(aStep.Timenow_sec, (passccd/5577720.)*100.);
      hist_eff_vs_eccd_ccd->Fill(ccdenergy);
      hist_eff_vs_elength_ccd->Fill(lengthenergy);
    
      /////////all events above have passed ccd cuts only
      
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	 

      if(aStep.Etrig_kev<400)
	{continue;}
      
      if(aStep.Emesh_kev<500)
	{continue;}
          
      passwfandccd++;
      file_passwfandccd<<aStep.RunNum<<" "<<aStep.EventNum<<endl;
      hist_time_wfandccd->Fill(aStep.Timenow_sec);
      hist_time_wfandccd_percent->Fill(aStep.Timenow_sec);
      hist_wfandccdpass_rate->Fill(aStep.Timenow_sec, passwfandccd);
      hist_wfandccdpass_rate_percent->Fill(aStep.Timenow_sec, (passwfandccd/5577720.)*100.);
      hist_eff_vs_eccd_wfandccd->Fill(ccdenergy);
      hist_eff_vs_elength_wfandccd->Fill(lengthenergy);

      ////////////////////////////////////////////////////////////////////////////////////////////////////

      double detectorefficiency=((passwfandccd*1.0)/(passccd*1.0));
      
      hist_time=(TH1D*)hist_time_wfandccd_percent->Clone();
      hist_time->Divide(hist_time_ccd);
      hist_time_percent=(TH1D*)hist_time->Clone();
      hist_pass_rate->Fill(aStep.Timenow_sec, detectorefficiency);
      hist_eff_vs_eccd=(TH1D*)hist_eff_vs_eccd_wfandccd->Clone();
      hist_eff_vs_eccd->Divide(hist_eff_vs_eccd_ccd);
      hist_eff_vs_eccd_percent=(TH1D*)hist_eff_vs_eccd->Clone();
      hist_eff_vs_elength=(TH1D*)hist_eff_vs_elength_wfandccd->Clone();
      hist_eff_vs_elength->Divide(hist_eff_vs_elength_ccd);
      hist_eff_vs_elength2=(TH1D*)hist_eff_vs_elength_wfandccd->Clone();
      hist_eff_vs_elength2->Divide(hist_eff_vs_elength_wf);

      neutronnumber++;
     
    }
 
  hist_time_percent->Scale(100);
  hist_eff_vs_eccd_percent->Scale(100);
   
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  Double_t fitf(Double_t *x, Double_t *par)
  {
    Double_t fitval = par[0]*x[0];
    return fitval;
  }   
    
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  new TCanvas;
  hist_time_ccd_percent->SetTitle("No. Events vs  Time");
  hist_time_ccd_percent->SetXTitle("Unix Time");
  hist_time_ccd_percent->SetYTitle("No. Events");
  hist_time_ccd_percent->SetLineColor(kRed);
  hist_time_ccd_percent->Draw();
  // hist_time_wf->SetLineColor(kBlue);
  //hist_time_wf->Draw("SAME");
  hist_time_wfandccd_percent->SetLineColor(kGreen);
  hist_time_wfandccd_percent->Draw("SAME");
  leg_time_percent = new TLegend(0.6, 0.7, 0.89, 0.89);
  leg_time_percent->AddEntry(hist_time_ccd_percent, "No. Events Passed CCD Cuts Only", "l");
  leg_time_percent->AddEntry(hist_time_wfandccd_percent, "No. Events Passed Both CCD&WF Cuts", "l");
  leg_time_percent->AddEntry(hist_time_wf, "No. Events Passed WF Cuts Only", "l");
  leg_time_percent->Draw("SAME");

  new TCanvas;
  hist_time->SetTitle("Passed WF& CCD/Passed CCD vs Time");
  hist_time->SetXTitle("Unix Time");
  hist_time->SetYTitle("No. Passed WF&CCD/No. Passed CCD");
  hist_time->Draw();
 
  new TCanvas;
  hist_ccdpass_rate->SetTitle("Rate of Events");
  hist_ccdpass_rate->SetXTitle("Unix Time");
  hist_ccdpass_rate->SetYTitle("No. of Events");
  hist_ccdpass_rate->SetMarkerColor(2);
  hist_ccdpass_rate->SetMarkerStyle(2);
  hist_ccdpass_rate->Draw();
  //hist_wfpass_rate->SetMarkerStyle(2);
  //hist_wfpass_rate->SetMarkerColor(4);
  hist_wfpass_rate->Draw("SAME");
  hist_wfandccdpass_rate->SetMarkerStyle(2);
  hist_wfandccdpass_rate->Draw("SAME");
  leg_pass_rate = new TLegend(0.6, 0.7, 0.89, 0.89);
  leg_pass_rate->AddEntry(hist_ccdpass_rate_percent, "Passed CCD Cuts Only - Red", "l");
  leg_pass_rate->AddEntry(hist_wfpass_rate, "Passed WF Cuts Only - Blue", "l");
  leg_pass_rate->AddEntry(hist_wfandccdpass_rate_percent, "Passed Both CCD&WF Cuts", "l");
  leg_pass_rate->Draw("SAME");

  new TCanvas;
  hist_pass_rate->SetTitle("No. Passed CCD&WF/No. Passed CCD Rate");
  hist_pass_rate->SetXTitle("Unix Time");
  hist_pass_rate->SetYTitle("No. Passed CCD&WF/No. Passed CCD");
  hist_pass_rate->SetMarkerStyle(2);
  hist_pass_rate->Draw();

  new TCanvas;
  hist_eff_vs_eccd_ccd->SetTitle("No. Events vs CCD Energy");
  hist_eff_vs_eccd_ccd->SetYTitle("No.Events");
  hist_eff_vs_eccd_ccd->SetXTitle("CCD Energy/keV");
  hist_eff_vs_eccd_ccd->SetLineColor(kRed);
  hist_eff_vs_eccd_ccd->Draw();
  hist_eff_vs_eccd_wfandccd->SetLineColor(kBlue);
  hist_eff_vs_eccd_wfandccd->Draw("SAME");
  leg_eff_vs_eccd = new TLegend(0.6, 0.7, 0.89, 0.89);
  leg_eff_vs_eccd->AddEntry(hist_eff_vs_eccd_ccd, "Passed CCD Cuts Only", "l");
  leg_eff_vs_eccd->AddEntry(hist_eff_vs_eccd_wfandccd, "Pass Both CCD&WF Cuts", "l");
  leg_eff_vs_eccd->Draw("SAME");

  new TCanvas;
  hist_eff_vs_eccd->SetTitle("No. Events Passed CCD&WF/No. Events Passed CCD vs CCD Energy");
  hist_eff_vs_eccd->SetXTitle("CCD Energy/keV");
  hist_eff_vs_eccd->SetYTitle("No. Events Passed CCD&WF/No. Events Passed CCD");
  hist_eff_vs_eccd->Draw();

  new TCanvas;
  hist_eff_vs_elength_ccd->SetTitle("No. Events vs Length Energy");
  hist_eff_vs_elength_ccd->SetYTitle("No. Events");
  hist_eff_vs_elength_ccd->SetXTitle("Length Energy/keV");
  hist_eff_vs_elength_ccd->SetLineColor(kRed);
  hist_eff_vs_elength_ccd->Draw();
  //hist_eff_vs_elength_wf->SetLineColor(kBlue);
  //hist_eff_vs_elength_wf->Draw("SAME");
  hist_eff_vs_elength_wfandccd->SetLineColor(kGreen);
  hist_eff_vs_elength_wfandccd->Draw("SAME");
  leg_eff_vs_elength = new TLegend(0.6, 0.7, 0.89, 0.89);
  leg_eff_vs_elength->AddEntry(hist_eff_vs_elength_ccd, "Passed CCD Cuts Only", "l");
  leg_eff_vs_elength->AddEntry(hist_eff_vs_elength_wf, "Passed WF Cuts Only", "l");
  leg_eff_vs_elength->AddEntry(hist_eff_vs_elength_wfandccd, "Passed Both CCD&WF Cuts", "l");
  leg_eff_vs_elength->Draw("SAME");

  new TCanvas;
  hist_eff_vs_elength->SetTitle("No.Events Passed CCD&WF/No. Events Passed CCD vs Length Energy");
  hist_eff_vs_elength->SetYTitle("No.Events Passed CCD&WF/No. Events Passed CCD");
  hist_eff_vs_elength->SetXTitle("Length Energy/keV");
  //hist_eff_vs_elength2->SetLineColor(kBlue);
  hist_eff_vs_elength->SetLineColor(kRed);
  hist_eff_vs_elength->Draw();
  //hist_eff_vs_elength2->Draw("SAME");
  leg_eff=new TLegend(0.6, 0.7, 0.89, 0.89);
  leg_eff->AddEntry(hist_eff_vs_elength, "Pass WF&CCD/Pass CCD", "l");
  // leg_eff->AddEntry(hist_eff_vs_elength2, "Pass WF&CCD/Pass WF", "l");

  cout<<""<<endl;

  cout<<"There were "<<neutronnumber<<" neutrons"<<endl;
  cout<<"Total Exposure time now "<<totalexposure/(60*60*24)<<" days"<<endl;
  cout<<"There were ~  "<<neutronnumber/(totalexposure/(60*60*24))<<" neutrons per day"<<endl;
  
  cout<<""<<endl;

  cout<<passccd<<" events passed the CCD cuts only"<<endl;
  cout<<passwf<<" events passed the WF cuts only"<<endl;
  cout<<passwfandccd<<" events passed the waveform and CCD cuts"<<endl;
  cout<<"Detector Efficiency (Pass WFandCCD/Pass CCD): "<<(detectorefficiency)*100.<<"%"<<endl;
  cout<<"Detector Efficiency (Pass WFandCCD/Pass WF): "<<(passwfandccd/passwf)*100.<<"%"<<endl;
  
  file_passccd.close();
  file_passwfandccd.close();
  
}

double fitPoly12Percent(double E)
{
  return sqrt(0.9741*0.9741*E*E+113.835*113.835)-62.675;
  //return -11.50+39.02*E+0.03459*E*E;                                                    
  //return 469 + 10*E + 0.548 *E*E -3.39e-3*E*E*E+7.5e-6*E*E*E*E;          
  
}               

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
