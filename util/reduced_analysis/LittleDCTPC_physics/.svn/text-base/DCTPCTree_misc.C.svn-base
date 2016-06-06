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
  TH1D *hist_meshpeak=new TH1D("Mesh Peak Time", "Mesh Peak Time", 100, 0, 0.00001);
  TH1D *hist_meshtotaltime=new TH1D("Mesh Total Time", "Mesh Total Time", 100, 0, 0.00001);
  
  //Miscellaneous Plots
  TH1D *hist_maxpixel=new TH1D("Max Pixel", "Max Pixel", 100, 0, 1000);
  TH1D *hist_trackangle = new TH1D ("Track Angle", "Track Angle", 360, -180., 180.);
  TH1D *hist_cutnum=new TH1D("Cut Number","Cut Number",20,0,20);
  TH2D *hist_pos=new TH2D("Track Position","Track Position",100,-53.4,220,100,-53.4,220);
  TH2D *hist_startpos = new TH2D("Track Start Position", "Track Start Position", 100, -53.4., 220., 100, -53.4., 220.);
  TH2D *hist_endpos = new TH2D("Track End Position", "Track End Position", 100, -53.4., 220., 100, -53.4., 220.);
  TH1D *hist_rho_start=new TH1D("Rho Start", "Rho Start", 50, 0, 150);
  TH1D *hist_rho_end=new TH1D("Rho End", "Rho End", 50, 0, 150);
  TH2D *hist_rhostart_vs_rhoend=new TH2D("Rho Start vs Rho End", "Rho Start vs Rho End", 150, 0, 150, 150, 0, 150);
  TH1D *hist_rho_start_sq= new TH1D("Rho Start Squared", "Rho Start Squared", 250, 0, 22500);
  TH1D *hist_rho_end_sq= new TH1D("Rho End Squared", "Rho End Squared", 250, 0, 22500);
  TH1D *hist_skewness= new TH1D("Track Skewness", "Track Skewness", 100, 0, 1.5);

  //Plots to see where the anode and mesh saturate
  TH1D *hist_mesh_max = new TH1D("Mesh Max", "Mesh Max", 100, 0,  0.5);
  TH1D *hist_anode_max=new TH1D("Anode Max", "Anode Max", 50, 0, 0.25);

  //Neutron number plots
  TH2D *hist_neutronsperrun=new TH2D ("Number of Neutrons per Run", "Neutrons per Run", 8000, 9000, 17000, 1400, 0, 1400);
  TH1D *hist_time=new TH1D("Neutron Rate", "Neutron Rate", 100, 1374000000 ,1386000000);//this gives us a plot for the number of neutrons recorded at a certain time
  TH1D *hist_time_allevents=new TH1D("Event Rate", "Event Rate", 100, 1374000000 ,1386000000);
  TH2D *hist_neutronrate=new TH2D("Neutrons vs Time", "Neutrons vs Time", 100, 1374000000 ,1386000000, 600, 0, 600);
  TH1D *hist_ccd_neutrons=new TH1D("CCD Energy for Events Kept", "CCD Energy for Events Kept", 100, 0, 10000);
  TH1D *hist_ccd_allevents=new TH1D("CCD Energy for All Events", "CCD Energy for All Events", 100, 0, 10000);
  
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

  for (int event = 0; event<200000/*nentries*/; event++)
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
      
      hist_cutnum->Fill(cutnum);
      hist_time_allevents->Fill(aStep.Timenow_sec);
      hist_ccd_allevents->Fill(ccdenergy);
            
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
	  cutnum=13;
	  continue;
	}
      
      /* if(tracklength>72)
	{
	  cutnum=14;
	  continue;
	}
      */

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	       
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

      hist_meshpeak->Fill(aStep.Mesh_peaktime_samp*secondspersample);
      hist_meshtotaltime->Fill(aStep.Mesh_totaltime_samp*secondspersample);
      
      hist_maxpixel->Fill(aStep.Track_maxpixel_ccdadu);
      hist_trackangle->Fill(aStep.Track_phi_deg);
      hist_rho_start->Fill(rho_start);
      hist_rho_end->Fill(rho_end);
      hist_rhostart_vs_rhoend->Fill(rho_start, rho_end);
      hist_rho_start_sq->Fill(rho_start*rho_start);
      hist_rho_end_sq->Fill(rho_end*rho_end);
      hist_skewness->Fill(TMath::Abs(aStep.Track_skewness));

      hist_mesh_max->Fill(aStep.Mesh_max_v);
      hist_anode_max->Fill(aStep.Anode_max_v);

      hist_ccd_neutrons->Fill(ccdenergy);
         
      //above I am filling the histograms with data that has passed all the cuts - hopefully neutrons!
         
      cutnum=0;
      neutronnumber++;
     
      hist_neutronsperrun->Fill(aStep.RunNum, neutronnumber);
      hist_time->Fill(aStep.Timenow_sec);
      hist_neutronrate->Fill(aStep.Timenow_sec, neutronnumber);
      
    }

  hist_ccd_cuts=(TH1D*)hist_ccd_allevents->Clone();
  hist_ccd_cuts->Add(hist_ccd_neutrons, -1);
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  Double_t fitf(Double_t *x, Double_t *par)
  {
    Double_t fitval = par[0]*x[0];
    return fitval;
  }   
  
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
  hist_neutronsperrun->SetTitle("Number of Neutrons per Run");
  hist_neutronsperrun->SetXTitle("Run Number");
  hist_neutronsperrun->SetYTitle("Total Number of Neutrons");
  hist_neutronsperrun->SetMarkerStyle(2);
  hist_neutronsperrun->Draw();
  hist_neutronsperrun->Fit(func4, "r n");
  func4->SetLineColor(kRed);
  func4->Draw("SAME");

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
  hist_meshrise_90->SetLineColor(kMagenta);
  hist_meshrise_90->Draw("SAME");
  leg_meshrise = new TLegend(0.6, 0.7, 0.89, 0.89);
  leg_meshrise->AddEntry(hist_meshrise_0, "R0", "l");
  leg_meshrise->AddEntry(hist_meshrise_10, "R10", "l");
  leg_meshrise->AddEntry(hist_meshrise_25, "R25", "l");
  leg_meshrise->AddEntry(hist_meshrise_50, "R50", "l");
  leg_meshrise->AddEntry(hist_meshrise_75, "R75", "l");
  leg_meshrise->AddEntry(hist_meshrise_90, "R90", "l");
  leg_meshrise->Draw("SAME");
  
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
  hist_meshfall_90->SetLineColor(kMagenta);
  hist_meshfall_90->Draw("SAME");
  leg_meshfall = new TLegend(0.6, 0.7, 0.89, 0.89);
  leg_meshfall->AddEntry(hist_meshfall_0, "F0", "l");
  leg_meshfall->AddEntry(hist_meshfall_10, "F10", "l");
  leg_meshfall->AddEntry(hist_meshfall_25, "F25", "l");
  leg_meshfall->AddEntry(hist_meshfall_50, "F50", "l");
  leg_meshfall->AddEntry(hist_meshfall_75, "F75", "l");
  leg_meshfall->AddEntry(hist_meshfall_90, "F90", "l");
  leg_meshfall->Draw("SAME");

  new TCanvas;
  hist_meshpeak->SetTitle("Mesh Peak Time / sec");
  hist_meshpeak->SetYTitle("Counts");
  hist_meshpeak->SetXTitle("Mesh Peak Time / sec");
  hist_meshpeak->Draw();

  new TCanvas;
  hist_meshtotaltime->SetTitle("Mesh Total Time / sec");
  hist_meshtotaltime->SetXTitle("Mesh Total Time / sec");
  hist_meshtotaltime->SetYTitle("Counts");
  hist_meshtotaltime->Draw();
  
  new TCanvas;
  hist_time->SetTitle("Neutron Events Collected in Time");
  hist_time->SetXTitle("Unix Time");
  hist_time->SetYTitle("Neutron Events");
  hist_time->Draw();

  new TCanvas;
  hist_time_allevents->SetTitle("Events Collected in Time");
  hist_time_allevents->SetXTitle("Unix Time");
  hist_time_allevents->SetYTitle("Events");
  hist_time_allevents->Draw();
  cout<<"Check: "<<hist_time_allevents->GetSum()<<endl;

  new TCanvas;
  hist_neutronrate->SetTitle("Neutron Rate");
  hist_neutronrate->SetXTitle("Unix Time");
  hist_neutronrate->SetYTitle("Neutron Number");
  hist_neutronrate->SetMarkerStyle(2);
  hist_neutronrate->Draw();
  hist_neutronrate->Fit(func5, "r n");
  func5->SetLineColor(kBlue);
  func5->Draw("SAME");

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
  hist_cutnum->SetTitle("Cut Number");
  hist_cutnum->SetXTitle("Cut Number");
  hist_cutnum->SetYTitle("Counts");
  hist_cutnum->Draw();

  new TCanvas;
  hist_ccd_allevents->SetTitle("CCD Energy");
  hist_ccd_allevents->SetYTitle("Counts");
  hist_ccd_allevents->SetXTitle("CCD Energy/keV");
  hist_ccd_allevents->Draw();
  hist_ccd_cuts->SetLineColor(kRed);
  hist_ccd_cuts->Draw("SAME");
  hist_ccd_neutrons->SetLineColor(kBlue);
  hist_ccd_neutrons->Draw("SAME");
  leg_ccd=new TLegend(0.6, 0.7, 0.89, 0.89);
  leg_ccd->AddEntry(hist_ccd_allevents, "All Events", "l");
  leg_ccd->AddEntry(hist_ccd_cuts, "Cut Events", "l");
  leg_ccd->AddEntry(hist_ccd_neutrons, "Kept Events", "l");
  leg_ccd->Draw("SAME");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  cout<<"Number of events cut= "<<nentries-neutronnumber<<endl;
  cout<<"Most used cut is cut number = "<<hist_cutnum->GetMaximumBin()<<endl;
  cout<<"Least used cut is cut number = "<<hist_cutnum->GetMinimumBin()<<endl;
  cout<<"Mean Maximum Pixel adu = "<<hist_maxpixel->GetMean()<<endl;

  cout<<""<<endl;

  cout<<"Mesh saturates at "<<hist_mesh_max->GetMaximumBin()<<"V"<<endl;
  cout<<"Anode saturates at "<<hist_anode_max->GetMaximumBin()<<"V"<<endl;
  
  cout<<""<<endl;

  cout<<"There were "<<neutronnumber<<" neutrons"<<endl;
  cout<<"Total Exposure time now "<<totalexposure/(60*60*24)<<" days"<<endl;
  cout<<"There were ~  "<<neutronnumber/(totalexposure/(60*60*24))<<" neutrons per day"<<endl;
  cout<<"Parameter 0  in Number of Neutrons per Run Plot "<<func4->GetParameter(0)<<endl;
  cout<<"Parameter 1 in Number of Neutrons per Run Plot "<<func4->GetParameter(1)<<endl;
  cout<<"Parameter 0 in Neutron Number vs Time Plot "<<func5->GetParameter(0)<<endl;
  cout<<"Parameter 1 in Neutron Number vs Time Plot "<<func5->GetParameter(1)<<endl;

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
  
