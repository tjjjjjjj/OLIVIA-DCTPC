#define DCTPCTree_cxx
#include "../input_files/DCTPCTree.h"
#include "../input_files/Calibration_constants_BigDCTPC.h";
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TRandom3.h"
#include "TGraph.h"
#include "TMath.h"
#include "TArc.h"
#include "TVirtualFitter.h"
TGraph *gr;

//This macro demonstrates how the calibration constants are used in making some alpha and neutron plots

void DCTPCTree::Loop()
{

TStopwatch timer;
gROOT->SetStyle("Plain");
gStyle->SetEndErrorSize(3);
gStyle->SetPalette(1,0);
gStyle->SetLineWidth(2);
gStyle->SetHistLineWidth(2);
gStyle->SetOptStat(kTRUE);
// gStyle->SetOptFit(kFALSE);
TH1::AddDirectory(false);
  
TFile *outtree = new TFile("$BigDCTPC_physics_input");
TTree *dctreepc = (TTree*)outtree->Get("dctpc_eventinfo");
DCTPCTree aStep(dctreepc);

TH2D *hist_energy=new TH2D("", "", 50, 0, 15000, 50, 0, 15000);
TH2D *hist_energy2=new TH2D("", "", 50, 0, 15000, 50, 0, 15000);
TH2D *hist_energy3=new TH2D("", "", 50, 0, 15000, 50, 0, 15000);

TH1D *hist_energy4=new TH1D("", "",40,0.,15000);
TH1D *hist_energy44=new TH1D("", "", 40, 0.,15000);
TH1D *hist_energy444=new TH1D("", "", 40, 0.,15000);
TH1D *hist_energy4444=new TH1D("", "", 40, 0.,15000);
TH1D *hist_energy44444=new TH1D("", "", 40, 0.,15000);
TH1D *hist_energy444444=new TH1D("", "", 40, 0.,15000);
TH1D *hist_energy4444444=new TH1D("", "", 40, 0.,15000);
TH1D *hist_energy44444444=new TH1D("", "", 40, 0.,15000);
TH1D *hist_energy444444444=new TH1D("", "", 40, 0.,15000);
TH1D *hist_energy4444444444=new TH1D("", "", 40, 0.,15000);
TH1D *hist_energy5=new TH1D("", "", 50,0.,15000);
TH2D *hist_test=new TH2D("", "",100,0.,15,100,0.,15);

TH2D *hist_study1=new TH2D("", "",50,0.,2.,50,0.,15);
TH2D *hist_study2=new TH2D("", "",50,0.,3000.,50,0.,15);
TH2D *hist_study3=new TH2D("", "",50,0.,200.,50,0.,15);
TH2D *hist_study4=new TH2D("", "",50,2.,12.,50,0.,15);
TH2D *hist_study5=new TH2D("2d", "2d",50,0.,400.,50,0.,15);
TH2D *hist_study6=new TH2D("", "",50,0.,400.,50,0.,15);
TH2D *hist_study7=new TH2D("", "",50,0.9,1.1,50,0.,15);
TH2D *hist_study8=new TH2D("", "",50,0.5,1.5,50,0.,15);
TH2D *hist_study9=new TH2D("", "",50,0.,3000.,50,0.,15);
TH2D *hist_study10=new TH2D("", "",50,0.,200.,50,0.,15);
TH2D *hist_study11=new TH2D("", "",50,0.,.0015,50,0.,15);
TH2D *hist_study12=new TH2D("", "",50,0.,.0005,50,0.,15);
TH2D *hist_study13=new TH2D("", "",50,0.07,.12,50,0.,15);
TH2D *hist_study14=new TH2D("", "",50,1.015,1.03,50,0.,15);
TH2D *hist_study15=new TH2D("", "",50,1000.,4000.,50,0.,15);
TH2D *hist_study16=new TH2D("", "",50,1391000000,1407500000,50,0.,20);
TH2D *hist_study17=new TH2D("", "",50,0.,4000.,50,0.,15);
TH2D *hist_study18=new TH2D("", "",50,0.,200.,50,0.,15);
TH2D *hist_study19=new TH2D("", "",50,0.,1.,50,0.,15);
TH1D *hist_study20=new TH1D("", "",50,-10.,10.);

TH1D *hist_energy6=new TH1D("", "",40,0,2000);
TH1D *hist_energy66=new TH1D("", "", 40, 0,2000);
TH1D *hist_energy666=new TH1D("", "", 40, 0,2000);
TH1D *hist_energy6666=new TH1D("", "", 40, 0,2000);
TH1D *hist_energy66666=new TH1D("", "", 40, 0,2000);
TH1D *hist_energy666666=new TH1D("", "", 40, 0,2000);
TH1D *hist_energy6666666=new TH1D("", "", 40, 0,2000);
TH1D *hist_energy66666666=new TH1D("", "", 40, 0,2000);
TH1D *hist_energy666666666=new TH1D("", "", 40, 0,2000);
TH1D *hist_energy6666666666=new TH1D("", "", 40, 0,2000);
TH1D *hist_energy7=new TH1D("", "", 50,0,2000);

TH2D *hist_pos=new TH2D("","",200,0,1024,200,0,1024);
TH2D *hist_energyrange=new TH2D("","",200,0,10000,1000,0,1024);

TH1D *hist_time=new TH1D("","",100,1391000000,1407500000);
TH2D *hist_counter=new TH2D("","",1000,1391000000,1407500000,5000,0,5000);

TH1D *hist_length=new TH1D("", "",800,0.,15);

double rstart=0.;
double rend=0.;
int points=0;
int points2=0;
double x[1000000],y[1000000],x2[1000000],y2[1000000];
double Etrack_kev_calib,Etrig_kev_calib,Emesh_kev_calib,TrackX_pix_calib,TrackY_pix_calib,radius,TrackX_pix_calib_start,TrackY_pix_calib_start,TrackX_pix_calib_end,TrackY_pix_calib_end;

double Lvertical;
double Ltrack;
double deltaPoly;
double Elength;
double SRIMpar[3] = {0.14,0.9110,1.3};
double Vdrift=1.15e06;
double convertSample2sec = 4.e-09;//second per sample (250 MS/s)

int counter=0;

  Long64_t nentries = dctreepc->GetEntries();
  cout << "Number of Entries: " << nentries << endl;
 
for (int event = 0; event<nentries; event++)
{
      aStep.GetEntry(event); 
      
      if(event%100000==0)
	  cout<<((double)event/(double)nentries)*100.<<"%"<<endl;	



//calibration of variables 

//set-by-set position calibration
TrackX_pix_calib = ((aStep.Track_x_start_pix+aStep.Track_x_end_pix)/2.) - POSITION_X_CALIB[aStep.SetNum];
TrackY_pix_calib = ((aStep.Track_y_start_pix+aStep.Track_y_end_pix)/2.) - POSITION_Y_CALIB[aStep.SetNum];

 
 
TrackX_pix_calib_start = aStep.Track_x_start_pix - POSITION_X_CALIB[aStep.SetNum];
TrackY_pix_calib_start = aStep.Track_y_start_pix - POSITION_Y_CALIB[aStep.SetNum];
TrackX_pix_calib_end = aStep.Track_x_end_pix - POSITION_X_CALIB[aStep.SetNum];
TrackY_pix_calib_end = aStep.Track_y_end_pix - POSITION_Y_CALIB[aStep.SetNum];

//sequence-by-sequence energy calibration
Etrack_kev_calib = aStep.Etrack_kev * ENERGY_CALIB_ccd[aStep.SequenceNum];
Etrig_kev_calib = aStep.Etrig_kev * ENERGY_CALIB_wf[aStep.SequenceNum]; 
Emesh_kev_calib = aStep.Emesh_kev * ENERGY_CALIB_wf[aStep.SequenceNum];

//radial-based energy calibration
radius=(pow(TrackX_pix_calib,2)+pow(TrackY_pix_calib,2))*pow(MMPERPIXEL/10.,2);  
Etrack_kev_calib*=(pow(radius+RADIAL_PARAM1_CALIB,RADIAL_PARAM2_CALIB)/RADIAL_PARAM0_CALIB);

rstart=sqrt(pow(TrackX_pix_calib_start,2)+pow(TrackY_pix_calib_start,2));
rend=sqrt(pow(TrackX_pix_calib_end,2)+pow(TrackY_pix_calib_end,2));

if(aStep.Anode_base_v>1.03)
continue;

if(aStep.Mesh_base_v>0.2)
continue;

// if(aStep.Anode_max_v>1.84)
// continue;
// 
// if(aStep.Mesh_max_v>1.75)
// continue;

if(aStep.Mesh_peak_v/Emesh_kev_calib>0.00015)
continue;


if (
     aStep.Edge==0 && 
     rstart<450. && 
     rend<450.  &&
       aStep.SequenceNum>11 
     )
	{	
	 if(Etrig_kev_calib>300.&&Emesh_kev_calib>300.&&Etrack_kev_calib>300.&&
	 Etrig_kev_calib<30000.&&Emesh_kev_calib<30000.&&Etrack_kev_calib<30000. 
	 
// && TMath::Abs((Etrig_kev_calib-Emesh_kev_calib)/(Etrig_kev_calib+Emesh_kev_calib))<0.3
// && TMath::Abs((Etrack_kev_calib-Emesh_kev_calib)/(Etrack_kev_calib+Emesh_kev_calib))<0.3
	 )
	 {
	 



// cout<<aStep.RunNum<<" "<<aStep.EventNum<<endl;
// cout<<aStep.Etrack_kev<<" "<<aStep.Etrig_kev<<" "<<aStep.Emesh_kev<<endl;

	 hist_energy->Fill(Emesh_kev_calib,Etrig_kev_calib);
	 hist_energy2->Fill(Emesh_kev_calib,Etrack_kev_calib);
	 hist_energy3->Fill(Etrig_kev_calib,Etrack_kev_calib);


	hist_pos->Fill(512+aStep.Track_x_start_pix,512+aStep.Track_y_start_pix);



Lvertical = (aStep.Mesh_R0time_samp+aStep.Mesh_F10time_samp)*convertSample2sec*Vdrift;
Ltrack = sqrt((aStep.Track_range_pix*MMPERPIXEL/10.)*(aStep.Track_range_pix*MMPERPIXEL/10.)+(Lvertical)*(Lvertical));
deltaPoly = SRIMpar[1]*SRIMpar[1]-4*(SRIMpar[0]*(SRIMpar[2]-Ltrack));
Elength = (-1*SRIMpar[1]+sqrt(deltaPoly))/(2*SRIMpar[0]);

hist_test->Fill(Emesh_kev_calib/1000.,Elength);

hist_length->Fill(Elength);

cout<<aStep.RunNum<<" "<<aStep.EventNum<<endl;

if(Elength<8.0&&Elength>7.0)
{

counter++;
}
// if(Elength<6.0&&Elength>5.0){
if(Elength<6.0&&Elength>5.0){
	hist_study1->Fill(aStep.Track_rms_ccdadu/aStep.Track_mean_ccdadu,Emesh_kev_calib/1000.);
	hist_study2->Fill(aStep.Anode_R0time_samp,Emesh_kev_calib/1000.);
	hist_study3->Fill(aStep.NextSpark,Emesh_kev_calib/1000.);
	hist_study4->Fill(aStep.Track_width_pix,Emesh_kev_calib/1000.);
	hist_study5->Fill(aStep.Track_range_pix,Emesh_kev_calib/1000.);
	hist_study6->Fill(radius,Emesh_kev_calib/1000.);
	hist_study7->Fill(Emesh_kev_calib/Etrig_kev_calib,Emesh_kev_calib/1000.);
	hist_study8->Fill(Emesh_kev_calib/Etrack_kev_calib,Emesh_kev_calib/1000.);
	hist_study9->Fill(aStep.Mesh_F0time_samp,Emesh_kev_calib/1000.);
	hist_study10->Fill(aStep.Track_rms_ccdadu,Emesh_kev_calib/1000.);
	hist_study11->Fill(aStep.Mesh_rms_v,Emesh_kev_calib/1000.);
	hist_study12->Fill(aStep.Anode_rms_v,Emesh_kev_calib/1000.);
	hist_study13->Fill(aStep.Mesh_base_v,Emesh_kev_calib/1000.);
	hist_study14->Fill(aStep.Anode_base_v,Emesh_kev_calib/1000.);
	hist_study15->Fill(aStep.Mesh_totaltime_samp,Emesh_kev_calib/1000.);
	hist_study16->Fill(aStep.Timenow_sec,Emesh_kev_calib/1000.);
	hist_study17->Fill(aStep.Mesh_F10time_samp+aStep.Mesh_R0time_samp,Emesh_kev_calib/1000.);
	hist_study18->Fill(aStep.Anode_starttime_samp-aStep.Mesh_starttime_samp,Emesh_kev_calib/1000.);
	
	
	hist_study19->Fill(Lvertical/Ltrack,Emesh_kev_calib/1000.);

	hist_study20->Fill(Lvertical-(sqrt(pow(Ltrack,2)-pow(Lvertical,2))/sqrt(2)));
	//cout<<aStep.Anode_R0time_samp<<" "<<Emesh_kev_calib/1000.<<endl;
	}

     hist_counter->Fill(aStep.Timenow_sec,counter);

     if(aStep.SequenceNum==7)
     {
	 hist_energy4->Fill(Emesh_kev_calib);
	 hist_energy6->Fill(Emesh_kev_calib);
	 }
	 if(aStep.SequenceNum==8)
	 {
	 hist_energy44->Fill(Emesh_kev_calib);
	 hist_energy66->Fill(Emesh_kev_calib);
	 }
	 if(aStep.SequenceNum==9)
	 {
	 hist_energy444->Fill(Emesh_kev_calib);
	 hist_energy666->Fill(Emesh_kev_calib);
	 }
	 if(aStep.SequenceNum==10)
	 {
	 hist_energy4444->Fill(Emesh_kev_calib);
	 hist_energy6666->Fill(Emesh_kev_calib);
	 }
	 if(aStep.SequenceNum==11)
	 {
	 hist_energy44444->Fill(Emesh_kev_calib);
	 hist_energy66666->Fill(Emesh_kev_calib);
	 }
	 if(aStep.SequenceNum==12)
	 {
	 hist_energy444444->Fill(Emesh_kev_calib);
	 hist_energy666666->Fill(Emesh_kev_calib);
	 }
	 if(aStep.SequenceNum==13)
	 {
	 hist_energy4444444->Fill(Emesh_kev_calib);
	 hist_energy6666666->Fill(Emesh_kev_calib);
	 }
	 
	 if(aStep.SequenceNum==14)
	 {
	 hist_energy44444444->Fill(Emesh_kev_calib);
	 hist_energy66666666->Fill(Emesh_kev_calib);
	 }
	 
	 if(aStep.SequenceNum==15)
	 {
	 hist_energy444444444->Fill(Emesh_kev_calib);
	 hist_energy666666666->Fill(Emesh_kev_calib);
	 }
	 
	 if(aStep.SequenceNum==16)
	 {
	 hist_energy4444444444->Fill(Emesh_kev_calib);
	 hist_energy6666666666->Fill(Emesh_kev_calib);
	 }
	
     if(Emesh_kev_calib>400. && Emesh_kev_calib<2000.)
     {
     hist_time->Fill(aStep.Timenow_sec);
	 }	
	  
	 hist_energy5->Fill(Emesh_kev_calib);
	 hist_energy7->Fill(Emesh_kev_calib);
	 hist_energyrange->Fill(Etrig_kev_calib,aStep.Track_range_pix);
	}

	}
	   
}

hist_energy4->Scale(1./EXPOSURE_SEQ[7]);
hist_energy44->Scale(1./EXPOSURE_SEQ[8]);
hist_energy444->Scale(1./EXPOSURE_SEQ[9]);
hist_energy4444->Scale(1./EXPOSURE_SEQ[10]);
hist_energy44444->Scale(1./EXPOSURE_SEQ[11]);
hist_energy444444->Scale(1./EXPOSURE_SEQ[12]);
hist_energy4444444->Scale(1./EXPOSURE_SEQ[13]);
hist_energy44444444->Scale(1./EXPOSURE_SEQ[14]);
hist_energy444444444->Scale(1./EXPOSURE_SEQ[15]);
hist_energy4444444444->Scale(1./EXPOSURE_SEQ[16]);

hist_energy5->Scale(1./(EXPOSURE_SEQ[7]+EXPOSURE_SEQ[8]+EXPOSURE_SEQ[9]+EXPOSURE_SEQ[10]+EXPOSURE_SEQ[11]+EXPOSURE_SEQ[12]+EXPOSURE_SEQ[13]+EXPOSURE_SEQ[14]+EXPOSURE_SEQ[15]+EXPOSURE_SEQ[16]));

hist_energy6->Scale(1./EXPOSURE_SEQ[7]);
hist_energy66->Scale(1./EXPOSURE_SEQ[8]);
hist_energy666->Scale(1./EXPOSURE_SEQ[9]);
hist_energy6666->Scale(1./EXPOSURE_SEQ[10]);
hist_energy66666->Scale(1./EXPOSURE_SEQ[11]);
hist_energy666666->Scale(1./EXPOSURE_SEQ[12]);
hist_energy6666666->Scale(1./EXPOSURE_SEQ[13]);
hist_energy66666666->Scale(1./EXPOSURE_SEQ[14]);
hist_energy666666666->Scale(1./EXPOSURE_SEQ[15]);
hist_energy6666666666->Scale(1./EXPOSURE_SEQ[16]);

hist_energy7->Scale(1./(EXPOSURE_SEQ[7]+EXPOSURE_SEQ[8]+EXPOSURE_SEQ[9]+EXPOSURE_SEQ[10]+EXPOSURE_SEQ[11]+EXPOSURE_SEQ[12]+EXPOSURE_SEQ[13]));

hist_energy->SetXTitle("Mesh energy (keVee)");
hist_energy->SetYTitle("Anode energy (keVee)");
hist_energy->Draw("COLZ");

new TCanvas;
hist_energy2->SetXTitle("Mesh energy (keVee)");
hist_energy2->SetYTitle("CCD energy (keVee)");
hist_energy2->Draw("COLZ");

new TCanvas;
hist_energy3->SetXTitle("Anode energy (keVee)");
hist_energy3->SetYTitle("CCD energy (keVee)");
hist_energy3->Draw("COLZ");

new TCanvas;
hist_energy44444->SetMinimum(0.);
hist_energy44444->SetMaximum(0.0002);
hist_energy44444->SetTitle("Alpha energies, start and end r<400 pixels");
hist_energy44444->SetXTitle("Mesh energy (keVee)");
hist_energy44444->GetYaxis()->SetTitleOffset(1.2);
hist_energy44444->SetYTitle("events/second/bin");
hist_energy44444->SetLineColor(5);
hist_energy44444->Draw();
hist_energy44->SetLineColor(2);
hist_energy44->Draw("SAME");
hist_energy444->SetLineColor(3);
hist_energy444->Draw("SAME");
hist_energy4444->SetLineColor(4);
hist_energy4444->Draw("SAME");
hist_energy4->SetLineColor(1);
hist_energy4->Draw("SAME");
hist_energy444444->SetLineColor(6);
hist_energy444444->Draw("SAME");
hist_energy4444444->SetLineColor(7);
hist_energy4444444->Draw("SAME");
hist_energy44444444->SetLineColor(8);
hist_energy44444444->Draw("SAME");
hist_energy444444444->SetLineColor(9);
hist_energy444444444->Draw("SAME");
hist_energy4444444444->SetLineColor(41);
hist_energy4444444444->Draw("SAME");

new TCanvas;
hist_energy5->SetTitle("Alpha energies, start and end r<400 pixels");
hist_energy5->SetXTitle("Mesh energy (keVee)");
hist_energy5->GetYaxis()->SetTitleOffset(1.2);
hist_energy5->SetYTitle("events/second/bin");
hist_energy5->Draw();


new TCanvas;
hist_energy66666->SetMinimum(0.);
hist_energy66666->SetMaximum(0.00004);
hist_energy66666->SetTitle("Neutron energies, start and end r<400 pixels");
hist_energy66666->SetXTitle("Mesh energy (keVee)");
hist_energy66666->GetYaxis()->SetTitleOffset(1.2);
hist_energy66666->SetYTitle("events/second/bin");
hist_energy66666->SetLineColor(5);
hist_energy66666->Draw();
hist_energy66->SetLineColor(2);
hist_energy66->Draw("SAME");
hist_energy666->SetLineColor(3);
hist_energy666->Draw("SAME");
hist_energy6666->SetLineColor(4);
hist_energy6666->Draw("SAME");
hist_energy6->SetLineColor(1);
hist_energy6->Draw("SAME");
hist_energy666666->SetLineColor(6);
hist_energy666666->Draw("SAME");
hist_energy6666666->SetLineColor(7);
hist_energy6666666->Draw("SAME");
hist_energy66666666->SetLineColor(8);
hist_energy66666666->Draw("SAME");
hist_energy666666666->SetLineColor(9);
hist_energy666666666->Draw("SAME");
hist_energy6666666666->SetLineColor(41);
hist_energy6666666666->Draw("SAME");

new TCanvas;
hist_energy7->SetTitle("Neutron energies, start and end r<400 pixels");
hist_energy7->SetXTitle("Mesh energy (keVee)");
hist_energy7->GetYaxis()->SetTitleOffset(1.2);
hist_energy7->SetYTitle("events/second/bin");
hist_energy7->Draw();

new TCanvas;
hist_pos->SetXTitle("X coordinate (pixels)");
hist_pos->SetYTitle("Y coordinate (pixels)");
hist_pos->Draw("COLZ");

new TCanvas;
hist_time->SetTitle("Neutron events collected in time, w/o efficiency correction");
hist_time->SetXTitle("Unix time");
hist_time->SetYTitle("Neutron events");
hist_time->Draw();

new TCanvas;
hist_test->Draw("COLZ");

new TCanvas;
hist_counter->Draw();

new TCanvas;
hist_study1->Draw("COLZ");

new TCanvas;
hist_study2->Draw("COLZ");

new TCanvas;
hist_study3->Draw("COLZ");

new TCanvas;
hist_study4->Draw("COLZ");

new TCanvas;
hist_study5->Draw("COLZ");

new TCanvas;
hist_study6->Draw("COLZ");

new TCanvas;
hist_study7->Draw("COLZ");

new TCanvas;
hist_study8->Draw("COLZ");

new TCanvas;
hist_study9->Draw("COLZ");

new TCanvas;
hist_study10->Draw("COLZ");

new TCanvas;
hist_study11->Draw("COLZ");

new TCanvas;
hist_study12->Draw("COLZ");

new TCanvas;
hist_study13->Draw("COLZ");

new TCanvas;
hist_study14->Draw("COLZ");

new TCanvas;
hist_study15->Draw("COLZ");

new TCanvas;
hist_study16->Draw("COLZ");

new TCanvas;
hist_study17->Draw("COLZ");

new TCanvas;
hist_study17->Draw("COLZ");

new TCanvas;
hist_study18->Draw("COLZ");

new TCanvas;
hist_study19->Draw("COLZ");

new TCanvas;
hist_study20->Draw();

new TCanvas;
hist_length->Draw();

}
