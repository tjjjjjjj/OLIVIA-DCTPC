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
  
  ofstream neutrons;
  neutrons.open("Neutrons.txt");
  ofstream cuts;
  cuts.open("Cuts.txt");
  ofstream neutronruns;
  neutronruns.open("Neutron_Runs.txt");

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //Energy plots
  TH1D *hist_energy_trig=new TH1D("Anode Energy","Anode Energy",50,0,2000); 
  TH1D *hist_energy_track=new TH1D("CCD Energy", "CCD Energy", 100, 0, 10000);//this is going so high because i still need to cut the double alpha events which only the CCD picks up
  TH1D *hist_energy_mesh=new TH1D("Mesh Energy", "Mesh Energy", 100, 0, 10000);
  TH1D *hist_energy_trig_log=new TH1D ("Anode Energy Log Scale", "Anode Energy Log Scale", 20, 0, 2000);
  TH1D *hist_energy_track_log=new TH1D("CCD Energy Log Scale", "CCD Energy Log Scale", 100, 0, 10000);
  TH1D *hist_energy_mesh_log=new TH1D("Mesh Energy Log Scale", "Mesh Energy Log Scale", 100, 0, 10000);
  TH1D *hist_energy_all=new TH1D("All Energies", "All Energies", 250, 0, 10000);
  TH1D *hist_energy_trig_notsat=new TH1D("Non Saturated Anode", "Non Saturated Anode", 50, 0, 2000);

  //misc energy plots
  TH2D *hist_mesh_over_ccd=new TH2D("Mesh Energy/CCD Energy vs CCD Energy", "Mesh Energy/CCD Energy vs CCD Energy", 300, 0, 3, 2000, 0, 10000);
  TH1D *hist_mesh_over_time = new TH1D("Mesh Energy / Time", "Mesh Energy / Time", 125, 0, 2.5e9);

  //energy from length plot
  TH1D *hist_energy_length=new TH1D("Length Energy", "Length Energy", 100, 0, 10000);
  TH1D *hist_energy_length_all=new TH1D("Length Energy All", "Length Energy All", 100, 0, 10000);
  TH1D *hist_divide_lengthenergy=new TH1D("CCD Divider", "CCD Divider", 100, 0, 10000);
  TH1D *hist_divide_wf=new TH1D("WF Divider", "WF Divider", 100, 0, 10000);
  TH1D *hist_divide_ccdenergy=new TH1D("CCD Energy Divider", "CCD Energy Divider", 100, 0, 10000);

  //neutron number
  TH1D *hist_neutrons=new TH1D("Neutrons", "Neutrons",2000, 0, 2000);
  // TH2D *hist_neutronno=new TH2D("Neutron No.", "Neutron No.", 8000, 9000, 17000, 2000, 0, 2000);
  TH1D *hist_prediction=new TH1D("Prediction", "Prediction", 100, 0, 10000);
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
      
      cuts<<aStep.RunNum<<" "<<aStep.EventNum<<endl;
      
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

      hist_energy_length_all->Fill(lengthenergy);

      
      /*  if(tracklength>72)//change to greater than 72 for neutrons
	{
	  cutnum=14;
	  continue;
	}
      
      */
      
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	  
    
      /* hist_mesh_over_ccd->Fill(mesh_over_ccd, aStep.Etrack_kev*ccdcalib*weight);
      
	 hist_energy_trig->Fill(aStep.Etrig_kev*anodecalib);*/
      hist_energy_track->Fill(ccdenergy);
      /* hist_energy_mesh->Fill(aStep.Emesh_kev*meshcalib);
      hist_energy_trig_log->Fill(aStep.Etrig_kev*anodecalib);
      hist_energy_track_log->Fill(ccdenergy);
      hist_energy_mesh_log->Fill(aStep.Emesh_kev*(meshcalib));
      if(aStep.Mesh_max_v<0.39)
	hist_energy_all->Fill(aStep.Emesh_kev*meshcalib);
      if(aStep.Etrig_kev<1250)
	hist_energy_trig_notsat->Fill(aStep.Etrig_kev*anodecalib);

	hist_mesh_over_time->Fill((aStep.Emesh_kev*meshcalib)/aStep.Mesh_totaltime_samp);*/
      hist_energy_length->Fill(lengthenergy);

      //above I am filling the histograms with data that has passed all the cuts - hopefully neutrons!
      
      neutrons<<aStep.RunNum<<" "<<aStep.EventNum<<endl;
      neutronruns<<aStep.RunNum<<endl;

      cutnum=0;
      if(tracklength<72)
      neutronnumber++;

      hist_neutrons->Fill(neutronnumber);
     
    }

  TF1 *divide_length = new TF1("pol1", "pol1", 400, 6000);
  TF1 *divide_ccd = new TF1("pol1", "pol1", 400, 6000);
  divide_length->SetParameters(0.273859, -3.18967e-05);
  divide_ccd->SetParameters(0.273859, 9.17716e-06);

  hist_energy_all->SetMaximum(500);
  hist_energy_track->SetMaximum(80);

  //MIGHT NEED DELETING !!!!!!!
  hist_energy_length_all->Scale(1./detectorefficiency);
  hist_energy_track->Scale(1./detectorefficiency);
  
  //scale after event loop, scale for each sequence and then scale for overall using 1/overall exposure
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  hist_energy_length_ccd=(TH1D*)hist_energy_length_all->Clone();
  hist_energy_length_ccd->Divide(divide_length);
  hist_energy_ccd_ccd=(TH1D*)hist_energy_track->Clone();
  hist_energy_ccd_ccd->Divide(divide_ccd);

    // hist_energy_length_wf=(TH1D*)hist_energy_length_all->Clone();
  //hist_energy_length_wf->Divide(hist_divide_wf);

  Double_t fitf(Double_t *x, Double_t *par)
  {
    Double_t fitval = par[0]*x[0];
    return fitval;
  }   
  
  /* TF1 *func1 = new TF1("pol1", "pol1", 500, 1200);//anode vs ccd
  
  TF1 *func2 = new TF1("pol1", "pol1", 0, 3000);//mesh vs ccd

  TF1 *func3 = new TF1("pol1", "pol1", 250, 1200);//anode vs mesh

  TF1 *func4 = new TF1("pol1", "pol1", 300, 3000);//mesh vs length

  TF1 *func5 = new TF1("pol1", "pol1", 300, 1200);//anode vs length

  TF1 *func6 = new TF1("pol1", "pol1", 300, 3000);//ccd vs length

  TF1 *func7 = new TF1("pol1", "pol1", 300, 4000);//cut length energy

  TF1 *func8 = new TF1("expo", "expo", 300, 4410);

  TF1 *func9 = new TF1("expo", "expo", 300, 3510);

  TF1 *func10 = new TF1("func10", "gaus(0)+expo(3)", 5000, 6000);

  TF1 *func11 = new TF1("expo", "expo", 300, 4410);*/
 
  //above func 1/2/3 is where I define line of best fit for 2D energy plots
  
  /* TF1 *gaus1 = new TF1("gaus", "gaus", -1, 1);
  
  TF1 *gaus2 = new TF1("gaus", "gaus", -1, 1);
  
  TF1 *gaus3 = new TF1("gaus", "gaus", -1, 1);

  TF1 *gaus4 = new TF1("gaus", "gaus", -1, 1);

  TF1 *gaus5 = new TF1("gaus", "gaus", -1, 1);
  
  TF1 *gaus6 = new TF1("gaus", "gaus", -1, 1);

  TF1 *gaus7 = new TF1("gaus", "gaus", 4400, 7000);

  TF1 *gaus8 = new TF1("gaus", "gaus", 3500, 7000);

  TF1 *gaus9 = new TF1("gaus", "gaus", 4400, 7000);*/

  TF1 *expo1 = new TF1("expo", "expo", 300, 4410);
  TF1 *gaus1 = new TF1("gaus", "gaus", 4410, 7000);
  TF1 *total1 = new TF1("total", "gaus(0)+expo(3)", 300, 7000);

  TF1 *expo2 = new TF1("expo", "expo", 300, 3510);
  TF1 *gaus2 = new TF1("gaus", "gaus", 3510, 7000);
  TF1 *total2 = new TF1("total", "gaus(0)+expo(3)", 300, 7000);

  TF1 *expo3 = new TF1("expo", "expo", 300, 5000);
  TF1 *gaus3 = new TF1("gaus", "gaus", 5000, 7000);
  TF1 *total3 = new TF1("total", "gaus(0)+expo(3)", 300, 6000);

  TF1 *expo4 = new TF1("expo", "expo", 300, 3510);
  TF1 *gaus4 = new TF1("gaus", "gaus", 3510, 7000);
  TF1 *total4 = new TF1("total", "gaus(0)+expo(3)", 300, 7000);

  TF1 *pol1=new TF1("pol1", "pol1", 400, 3500);

  // double_t par[5];
  
  //above gaus 1/2/3 is where I define the fits for the delta distributions

  TArc *innerveto = new TArc(512*MMPERPIXEL,512*MMPERPIXEL,124);
  TArc *outerveto = new TArc(512*MMPERPIXEL,512*MMPERPIXEL,134);
  
  //above is where I define the veto ring shape for track position plot

  TLine *threshold = new TLine();
  threshold->SetX1(6500);
  threshold->SetX2(6500);
  threshold->SetY1(0);
  threshold->SetY2(100);
  TLine *alphasstart = new TLine();
  alphasstart->SetX1(4000);
  alphasstart->SetX2(4000);
  alphasstart->SetY1(0);
  alphasstart->SetY2(100);
  TLine *alphaspeak = new TLine();
  alphaspeak->SetX1(5300);
  alphaspeak->SetX2(5300);
  alphaspeak->SetY1(0);
  alphaspeak->SetY2(100);
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /* new TCanvas;
  c1->SetLogy();
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
  c1_n3->SetLogy();
  hist_energy_track_log->SetXTitle("CCD Energy/keV");
  hist_energy_track_log->SetYTitle("Log Counts");
  hist_energy_track_log->SetTitle("Calibrated CCD Energy Log Scale");
  hist_energy_track_log->Draw();
  */
  new TCanvas;
  hist_energy_track->SetTitle("Calibrated CCD Energy");
  hist_energy_track->SetXTitle("Calibrated CCD Energy/keV");
  hist_energy_track->SetYTitle("Counts");
  hist_energy_track->Fit(expo2, "r n");
  hist_energy_track->Fit(gaus2, "r n");
  //cout<<gaus2->GetParameter(0)<<endl;
  //cout<<gaus2->GetParameter(1)<<endl;
  //cout<<gaus2->GetParameter(2)<<endl;
  //cout<<expo2->GetParameter(0)<<endl;
  //cout<<expo2->GetParameter(1)<<endl;
  //expo2->SetLineColor(kBlue);
  //gaus2->SetLineColor(kBlue);
  hist_energy_track->Draw("E1");
  //expo2->Draw("SAME");
  //gaus2->Draw("SAME");
  total2->SetParameters(36.3646, 4848.43, 896.082, 2.74022, -1.31322e-05);
  hist_energy_track->Fit(total2, "r n");
  total2->SetLineColor(kBlue);
  total2->Draw("SAME");
  /*
  new TCanvas;
  c1_n5->SetLogy();
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
  hist_energy_all->SetTitle("No Saturated Mesh Energy");
  hist_energy_all->SetYTitle("Counts");
  hist_energy_all->SetXTitle("Energy / keV");
  hist_energy_all->Draw();

  new TCanvas;
  // hist_energy_trig_notsat->SetTitle("All Energies");
  // hist_energy_trig_notsat->SetXTitle("Energy / keV");
  //hist_energy_trig_notsat->SetYTitle("Counts");   
  //hist_energy_trig_notsat->SetLineColor(kBlue);
  //hist_energy_trig_notsat->Draw();
  hist_energy_track->SetTitle("CCD and Length Energy");
  hist_energy_track->SetXTitle("Energy / keV");
  hist_energy_track->SetYTitle("Counts");
  hist_energy_track->SetLineColor(kRed);
  hist_energy_track->Draw();
  hist_energy_length->Draw("SAME");
  //hist_energy_all->SetLineColor(kGreen);
  //hist_energy_all->Draw("SAME");
  //hist_energy_length->Draw("SAME");
  leg_energy=new TLegend(0.6, 0.7, 0.89, 0.89);
  //leg_energy->AddEntry(hist_energy_all, "Mesh Energy", "l");
  leg_energy->AddEntry(hist_energy_track, "CCD Energy", "l");
  //leg_energy->AddEntry(hist_energy_trig, "Anode Energy", "l");
  leg_energy->AddEntry(hist_energy_length, "Energy Reconstructed from Length", "l");
  leg_energy->Draw();

  new TCanvas;
  hist_mesh_over_ccd->SetTitle("Mesh Energy keV/CCD Energy keV vs CCD Energy keV");
  hist_mesh_over_ccd->SetXTitle("Mesh Energy/keV / CCD Energy/keV");
  hist_mesh_over_ccd->SetYTitle("CCD Energy/keV");
  hist_mesh_over_ccd->Draw("COLZ");
  
  new TCanvas;
  hist_energy_length->SetTitle("Energy Reconstrcuted from Track Length");
  hist_energy_length->SetXTitle("Energy from Length / keV");
  hist_energy_length->SetYTitle("Counts");
  hist_energy_length->Draw();
  hist_energy_length->Fit(func8, "r n");
  func7->SetLineColor(kRed);
  func7->Draw("SAME");*/
  
  new TCanvas;
  hist_energy_length_all->SetTitle("Energy for All Events Reconstructed from Track Length");
  hist_energy_length_all->SetYTitle("Counts");
  hist_energy_length_all->SetXTitle("Energy from Length / keV");
  hist_energy_length_all->Draw("E1");
  hist_energy_length_all->Fit(expo1, "r n");
  hist_energy_length_all->Fit(gaus1, "r n");
  // expo1->SetLineColor(kRed);
  //gaus1->SetLineColor(kRed);
  //expo1->Draw("SAME");
  //gaus1->Draw("SAME");
  // cout<<gaus1->GetParameter(0)<<endl;
  //cout<<gaus1->GetParameter(1)<<endl;
  //cout<<gaus1->GetParameter(2)<<endl;
  //cout<<expo1->GetParameter(0)<<endl;
  //cout<<expo1->GetParameter(1)<<endl;
  total1->SetParameters(2.88544, 5204.83, 599.398, 0.794064, -0.000107004);
  hist_energy_length_all->Fit(total1, "r n");
  total1->SetLineColor(kRed);
  total1->Draw("SAME");
  hist_energy_track->Fit(expo2, "r n");
  hist_energy_track->Fit(gaus2, "r n");
  hist_energy_track->SetLineColor(kGray);
  hist_energy_track->Draw("SAME E1");
  total2->SetLineColor(kBlue);
  total2->Draw("SAME");
  //cout<<gaus2->GetParameter(0)<<endl;
  //cout<<gaus2->GetParameter(1)<<endl;
  //cout<<gaus2->GetParameter(2)<<endl;
  //cout<<expo2->GetParameter(0)<<endl;
  //cout<<expo2->GetParameter(1)<<endl;
  //expo2->SetLineColor(kBlue);
  //gaus2->SetLineColor(kBlue);
  threshold->SetLineColor(kGreen-9);
  threshold->Draw("SAME");
  alphasstart->SetLineColor(kGreen-9);
  alphasstart->Draw("SAME");
  alphaspeak->SetLineColor(kGreen-9);
  alphaspeak->Draw("SAME");

  new TCanvas;
  hist_energy_length_ccd->SetTitle("Length Energy Scaled by Efficiency");
  hist_energy_length_ccd->SetYTitle("Counts");
  hist_energy_length_ccd->SetXTitle("Energy from Length / keV");
  hist_energy_length_ccd->Draw();
  hist_divide_lengthenergy->SetLineColor(kMagenta);
  hist_divide_lengthenergy->Draw();
  //hist_energy_length_ccd->Fit(expo3, "r n");
  //hist_energy_length_ccd->Fit(gaus3, "r n");
  //expo3->SetLineColor(kYellow);
  //gaus3->SetLineColor(kYellow);
  //expo3->Draw("SAME");
  //gaus3->Draw("SAME");
  //cout<<gaus3->GetParameter(0)<<endl;
  //cout<<gaus3->GetParameter(1)<<endl;
  //cout<<gaus3->GetParameter(2)<<endl;
  //cout<<expo3->GetParameter(0)<<endl;
  //cout<<expo3->GetParameter(1)<<endl;
  //total3->SetParameters(22.2347, 5073.87, 586.897, -0.910078, 0.000587304);
  //total3->SetLineColor(kYellow);
  //total3->Draw("SAME");
  //  hist_energy_length_ccd->SetLineColor(kBlue);
  //hist_energy_length_ccd->Fit(func11, "r n");
  //hist_energy_length_ccd->Fit(gaus9, "r n");
  //func11->SetLineColor(kGreen);
  //gaus9->SetLineColor(kGreen);
  //func11->Draw("SAME");
  //gaus9->Draw("SAME");
  //hist_energy_length_ccd->Draw("SAME");
  // hist_energy_length_wf->SetLineColor(kGreen);
  //hist_energy_length_wf->Draw("SAME");

  new TCanvas;
  hist_energy_ccd_ccd->SetTitle("CCD Energy Scaled by Efficiency");
  hist_energy_ccd_ccd->SetYTitle("Counts");
  hist_energy_ccd_ccd->SetXTitle("CCD Energy / keV");
  hist_energy_ccd_ccd->Draw();
  hist_divide_ccdenergy->SetLineColor(kYellow);
  hist_divide_ccdenergy->Draw();
  //hist_energy_ccd_ccd->Fit(expo4, "r n");
  //hist_energy_ccd_ccd->Fit(gaus4, "r n");
  //expo4->SetLineColor(kGreen);
  //gaus4->SetLineColor(kGreen);
  //expo4->Draw("SAME");
  //gaus4->Draw("SAME");
  //cout<<gaus4->GetParameter(0)<<endl;
  //cout<<gaus4->GetParameter(1)<<endl;
  //cout<<gaus4->GetParameter(2)<<endl;
  //cout<<expo4->GetParameter(0)<<endl;
  //cout<<expo4->GetParameter(1)<<endl;
  
  new TCanvas;
  hist_neutrons->SetTitle("Neutrons");
  hist_neutrons->SetYTitle("Neutron Number");
  hist_neutrons->SetXTitle("Length Energy");
  hist_neutrons->SetMarkerStyle(2);
  hist_neutrons->Draw();
  hist_neutrons->Fit(pol1, "r n");
  pol1->SetLineColor(kRed);
  pol1->Draw("SAME");
  cout<<pol1->GetParameter(0)<<" "<<pol1->GetParameter(1)<<endl;

  /* new TCanvas;
  hist_mesh_over_time->SetTitle("Mesh Energy / Mesh Total Time");
  hist_mesh_over_time->SetYTitle("Counts");
  hist_mesh_over_time->SetXTitle("Mesh Energy / Mesh Total Time, keV/sec");
  hist_mesh_over_time->Draw();*/
 

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /* cout<<"Mean Anode Energy= "<<hist_energy_trig->GetMean()<<endl;
  cout<<"Modal Anode Energy= "<<hist_energy_trig->GetMaximumBin()<<endl;
  cout<<"Mean CCD Energy = "<<hist_energy_track->GetMean()<<endl;
  cout<<"Modal CCD Energy= "<<hist_energy_trig->GetMaximumBin()<<endl;
  cout<<"Mean Mesh Energy= "<<hist_energy_mesh->GetMean()<<endl;
  cout<<"Modal Mesh Energy= "<<hist_energy_mesh->GetMaximumBin()<<endl;
  cout<<"Mean Length Energy= "<<hist_energy_length->GetMean()<<endl;
  cout<<"Modal Length Energy= "<<hist_energy_length->GetMaximumBin()<<endl;

  cout<<""<<endl;*/
  /*
  cout<<"Length Energy Linear Fit Parameter 0 = "<<func7->GetParameter(0)<<endl;
  cout<<"Length Energy Linear Fit Parameter 1 = "<<func7->GetParameter(1)<<endl;
  cout<<"Length Energy Gaussian Fit Parameter 0 = "<<gaus7->GetParameter(0)<<endl;
  cout<<"Length Energy Gaussian Fit Parameter 1 = "<<gaus7->GetParameter(1)<<endl;
  cout<<"Length Energy Gaussian Fit Parameter 2 = "<<gaus7->GetParameter(2)<<endl;
  cout<<"CCD Energy Gaussian Fit Parameter 0 = "<<gaus8->GetParameter(0)<<endl;
  cout<<"CCD Energy Gaussian Fit Parameter 1 = "<<gaus8->GetParameter(1)<<endl;
  cout<<"CCD Energy Gaussian Fit Parameter 2 = "<<gaus8->GetParameter(2)<<endl;
  */
  cout<<""<<endl;

  //cout<<"Length Energy Resolution: "<<((gaus7->GetParameter(2))/(gaus7->GetParameter(1)))*100.<<"%"<<endl;
  //cout<<"CCD Energy Resolution: "<<((gaus8->GetParameter(2))/(gaus8->GetParameter(1)))*100.<<"%"<<endl;

  cout<<""<<endl;

  //cout<<"Area under exponential curve : "<<func8->Integral(300, 7000)<<endl;
  //cout<<"Neutrons per day = : "<<(func8->Integral(300, 7000))/(totalexposure/(60*60*24))<<endl;
  cout<<"Length Energy"<<endl;
  cout<<"Parameter 0 : "<<total1->GetParameter(0)<<endl;
  cout<<"Parameter 1 : "<<total1->GetParameter(1)<<endl;
  cout<<"Parameter 2 : "<<total1->GetParameter(2)<<endl;
  cout<<"Parameter 3 : "<<total1->GetParameter(3)<<endl;
  cout<<"Parameter 4 : "<<total1->GetParameter(4)<<endl;
  cout<<"Parameter 5 : "<<total1->GetParameter(5)<<endl;

  cout<<"Length Energy Resolution : "<<((total1->GetParameter(2))/(total1->GetParameter(1)))*100.<<"%"<<endl;
  cout<<"Alpha mean : "<<total1->GetParameter(1)<<endl;
  cout<<"ChiSquare : "<<total1->GetChisquare()<<endl;
  cout<<"Degrees of Freedom : "<<total1->GetNDF()<<endl;
  cout<<"Probability that fit is correct: "<<(TMath::Prob(total1->GetChisquare(), total1->GetNDF()))*100.<<"%"<<endl;
  
  cout<<""<<endl;

  cout<<"Area Under graph : "<<total1->Integral(400, 6000)<<endl;
  cout<<"Sum of Y Axis : "<<hist_energy_length_all->GetSum()<<endl;
  cout<<"Neutrons per day : "<<hist_energy_length_all->GetSum()/(totalexposure/(60*60*24))<<endl;
  cout<<"Area Under exponential part: "<<expo1->Integral(400, 6000)<<endl;
  cout<<"Neutron number prediction : "<<(hist_energy_length_all->GetSum()/total1->Integral(400, 6000))*expo1->Integral(400, 6000)<<endl;
  cout<<"Neutrons per day prediction : "<<((hist_energy_length_all->GetSum()/total1->Integral(400, 6000))*expo1->Integral(400, 6000))/(totalexposure/(60*60*24))<<endl;

  cout<<""<<endl;

  cout<<"CCD Energy"<<endl;
  cout<<"Parameter 0 : "<<total2->GetParameter(0)<<endl;
  cout<<"Parameter 1 : "<<total2->GetParameter(1)<<endl;
  cout<<"Parameter 2 : "<<total2->GetParameter(2)<<endl;
  cout<<"Parameter 3 : "<<total2->GetParameter(3)<<endl;
  cout<<"Parameter 4 : "<<total2->GetParameter(4)<<endl;
  cout<<"Parameter 5 : "<<total2->GetParameter(5)<<endl;

  cout<<"CCD Energy Resolution : "<<((total2->GetParameter(2))/(total2->GetParameter(1)))*100.<<"%"<<endl;
  cout<<"Alpha mean : "<<total2->GetParameter(1)<<endl;
  cout<<"ChiSquare : "<<total2->GetChisquare()<<endl;
  cout<<"Degrees of Freedom : "<<total2->GetNDF()<<endl;
  cout<<"Probability that fit is correct: "<<(TMath::Prob(total2->GetChisquare(), total2->GetNDF()))*100.<<"%"<<endl;

  cout<<""<<endl;

  cout<<"Length Energy Scaled to Efficiency"<<endl;
  cout<<"Parameter 0 : "<<total3->GetParameter(0)<<endl;
  cout<<"Parameter 1 : "<<total3->GetParameter(1)<<endl;
  cout<<"Parameter 2 : "<<total3->GetParameter(2)<<endl;
  cout<<"Parameter 3 : "<<total3->GetParameter(3)<<endl;
  cout<<"Parameter 4 : "<<total3->GetParameter(4)<<endl;
  cout<<"Parameter 5 : "<<total3->GetParameter(5)<<endl;

  // cout<<"Length Energy Scaled to Efficiency Resolution : "<<((total3->GetParameter(2))/(total3->GetParameter(1)))*100.<<"%"<<endl;
  cout<<"Alpha mean : "<<total3->GetParameter(1)<<endl;
  cout<<"ChiSquare : "<<total3->GetChisquare()<<endl;
  cout<<"Degrees of Freedom : "<<total3->GetNDF()<<endl;
  cout<<"Probability that fit is correct: "<<(TMath::Prob(total3->GetChisquare(), total3->GetNDF()))*100.<<"%"<<endl;

  cout<<""<<endl;

  cout<<"There were "<<neutronnumber/detectorefficiency<<" neutrons"<<endl;
  cout<<"Total Exposure time now "<<totalexposure/(60*60*24)<<" days"<<endl;
  cout<<"There were ~  "<<(neutronnumber/(totalexposure/(60*60*24)))/detectorefficiency<<" neutrons per day"<<endl;

  cuts.close();
  neutrons.close();
  neutronruns.close();
}

double fitPoly12Percent(double E)
{
  return sqrt(0.9741*0.9741*E*E+113.835*113.835)-62.675;
  //return -11.50+39.02*E+0.03459*E*E;                                                    
  //return 469 + 10*E + 0.548 *E*E -3.39e-3*E*E*E+7.5e-6*E*E*E*E;          
  
}               

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
