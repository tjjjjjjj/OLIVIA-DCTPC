#define DCTPCTree_cxx
#include "../input_files/DCTPCTree.h"
#include "../input_files/DCTPC_runtree.h"
#include "../input_files/Calibration_constants_LittleDCTPC.h";
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <string>
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
  
  //Scaled energy plots
  TH1D *hist_energy_mesh2=new TH1D ("Scaled Mesh Energy", "Scaled Mesh Energy", 20, 0, 5000);
  TH1D *hist_energy_mesh3=new TH1D ("Scaled Mesh Energy", "Scaled Mesh Energy", 20, 0, 5000);
  TH1D *hist_energy_mesh4=new TH1D ("Scaled Mesh Energy", "Scaled Mesh Energy", 20, 0, 5000);
  TH1D *hist_energy_mesh5=new TH1D ("Scaled Mesh Energy", "Scaled Mesh Energy", 20, 0, 5000);
  TH1D *hist_energy_mesh6=new TH1D ("Scaled Mesh Energy", "Scaled Mesh Energy", 20, 0, 5000);
  TH1D *hist_energy_mesh7=new TH1D ("Scaled Mesh Energy", "Scaled Mesh Energy", 20, 0, 5000);
  TH1D *hist_energy_mesh8=new TH1D ("Scaled Mesh Energy", "Scaled Mesh Energy", 20, 0, 5000);
  TH1D *hist_energy_mesh9=new TH1D ("Scaled Mesh Energy", "Scaled Mesh Energy", 20, 0, 5000);
  TH1D *hist_energy_mesh10=new TH1D ("Scaled Mesh Energy", "Scaled Mesh Energy", 20, 0, 5000);
  TH1D *hist_energy_mesh11=new TH1D ("Scaled Mesh Energy", "Scaled Mesh Energy", 20, 0, 5000);
  TH1D *hist_energy_mesh12=new TH1D ("Scaled Mesh Energy", "Scaled Mesh Energy", 20, 0, 5000);
  TH1D *hist_energy_mesh13=new TH1D ("Scaled Mesh Energy", "Scaled Mesh Energy", 20, 0, 5000);
  TH1D *hist_energy_mesh14=new TH1D ("Scaled Mesh Energy", "Scaled Mesh Energy", 20, 0, 5000);
  TH1D *hist_energy_mesh15=new TH1D ("Scaled Mesh Energy", "Scaled Mesh Energy", 20, 0, 5000);
  TH1D *hist_energy_mesh16=new TH1D ("Scaled Mesh Energy", "Scaled Mesh Energy", 20, 0, 5000);
  TH1D *hist_energy_mesh17=new TH1D ("Scaled Mesh Energy", "Scaled Mesh Energy", 20, 0, 5000);
  TH1D *hist_energy_mesh18=new TH1D ("Scaled Mesh Energy", "Scaled Mesh Energy", 20, 0, 5000);
  TH1D *hist_energy_meshall=new TH1D ("Scaled Mesh Energy", "Scaled Mesh Energy", 100, 0, 5000);
  TH1D *hist_energy_mesh=new TH1D("Mesh Energy", "Mesh Energy", 100, 0, 10000);
  TH1D *hist_energy_ccd=new TH1D("CCD Energy", "CCD Energy", 100, 0, 10000);
  TH1D *hist_energy_anode=new TH1D("Anode Energy", "Anode Energy", 20, 0, 2000);
  //these plots give us a plot of neutron rate vs energy


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
      
      /*if(tracklength>72)
	{
	  cutnum=14;
	  continue;
	  }*/
      
                  
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	 

	  //above I am filling the histograms with data that has passed all the cuts - hopefully neutrons!
	       
      cutnum=0;
      neutronnumber++;
          
      if(aStep.SequenceNum==2)
	hist_energy_mesh2->Fill(aStep.Emesh_kev*meshcalib);
      
      if(aStep.SequenceNum==3)
	    hist_energy_mesh3->Fill(aStep.Emesh_kev*meshcalib);
      
      if(aStep.SequenceNum==4)
	hist_energy_mesh4->Fill(aStep.Emesh_kev*meshcalib);
      
      if(aStep.SequenceNum==5)
	hist_energy_mesh5->Fill(aStep.Emesh_kev*meshcalib);
      
      if(aStep.SequenceNum==6)
	hist_energy_mesh6->Fill(aStep.Emesh_kev*meshcalib);
      
      if(aStep.SequenceNum==7)
	hist_energy_mesh7->Fill(aStep.Emesh_kev*meshcalib);
      
      if(aStep.SequenceNum==8)
	hist_energy_mesh8->Fill(aStep.Emesh_kev*meshcalib);
	  
      if(aStep.SequenceNum==9)
	hist_energy_mesh9->Fill(aStep.Emesh_kev*meshcalib);
      
      if(aStep.SequenceNum==10)
	hist_energy_mesh10->Fill(aStep.Emesh_kev*meshcalib);
	  
      if(aStep.SequenceNum==11)
	hist_energy_mesh11->Fill(aStep.Emesh_kev*meshcalib);
      
      if(aStep.SequenceNum==12)
	hist_energy_mesh12->Fill(aStep.Emesh_kev*meshcalib);
      
      if(aStep.SequenceNum==13)
	hist_energy_mesh13->Fill(aStep.Emesh_kev*meshcalib);
      
      if(aStep.SequenceNum==14)
	hist_energy_mesh14->Fill(aStep.Emesh_kev*meshcalib);
      
      if(aStep.SequenceNum==15)
	hist_energy_mesh15->Fill(aStep.Emesh_kev*meshcalib);
      
      if(aStep.SequenceNum==16)
	hist_energy_mesh16->Fill(aStep.Emesh_kev*meshcalib);
      
      if(aStep.SequenceNum==17)
	hist_energy_mesh17->Fill(aStep.Emesh_kev*meshcalib);
      
      if(aStep.SequenceNum==18)
	hist_energy_mesh18->Fill(aStep.Emesh_kev*meshcalib);
      
      hist_energy_meshall->Fill(aStep.Emesh_kev*meshcalib);
      hist_energy_mesh->Fill(aStep.Emesh_kev*meshcalib);
      hist_energy_ccd->Fill(ccdenergy);
      hist_energy_anode->Fill(aStep.Etrig_kev*anodecalib);

    }
  
  hist_energy_mesh2->Scale(1./EXPOSURE_SEQ[1]);
  hist_energy_mesh3->Scale(1./EXPOSURE_SEQ[2]);
  hist_energy_mesh4->Scale(1./EXPOSURE_SEQ[3]);
  hist_energy_mesh5->Scale(1./EXPOSURE_SEQ[4]);
  hist_energy_mesh6->Scale(1./EXPOSURE_SEQ[5]);
  hist_energy_mesh7->Scale(1./EXPOSURE_SEQ[6]);
  hist_energy_mesh8->Scale(1./EXPOSURE_SEQ[7]);
  hist_energy_mesh9->Scale(1./EXPOSURE_SEQ[8]);
  hist_energy_mesh10->Scale(1./EXPOSURE_SEQ[9]);
  hist_energy_mesh11->Scale(1./EXPOSURE_SEQ[10]);
  hist_energy_mesh12->Scale(1./EXPOSURE_SEQ[11]);
  hist_energy_mesh13->Scale(1./EXPOSURE_SEQ[12]);
  hist_energy_mesh14->Scale(1./EXPOSURE_SEQ[13]);
  hist_energy_mesh15->Scale(1./EXPOSURE_SEQ[14]);
  hist_energy_mesh16->Scale(1./EXPOSURE_SEQ[15]);
  hist_energy_mesh17->Scale(1./EXPOSURE_SEQ[16]);
  hist_energy_mesh18->Scale(1./EXPOSURE_SEQ[17]);

  hist_energy_meshall->Scale(1./EXPOSURE_SEQ[1]+EXPOSURE_SEQ[2]+EXPOSURE_SEQ[3]+EXPOSURE_SEQ[4]+EXPOSURE_SEQ[5]+EXPOSURE_SEQ[6]+EXPOSURE_SEQ[7]+EXPOSURE_SEQ[8]+EXPOSURE_SEQ[9]+EXPOSURE_SEQ[10]+EXPOSURE_SEQ[11]+EXPOSURE_SEQ[12]+EXPOSURE_SEQ[13]+EXPOSURE_SEQ[14]+EXPOSURE_SEQ[15]+EXPOSURE_SEQ[16]+EXPOSURE_SEQ[17]);

  hist_energy_mesh->Scale(1./EXPOSURE_SEQ[1]+EXPOSURE_SEQ[2]+EXPOSURE_SEQ[3]+EXPOSURE_SEQ[4]+EXPOSURE_SEQ[5]+EXPOSURE_SEQ[6]+EXPOSURE_SEQ[7]+EXPOSURE_SEQ[8]+EXPOSURE_SEQ[9]+EXPOSURE_SEQ[10]+EXPOSURE_SEQ[11]+EXPOSURE_SEQ[12]+EXPOSURE_SEQ[13]+EXPOSURE_SEQ[14]+EXPOSURE_SEQ[15]+EXPOSURE_SEQ[16]+EXPOSURE_SEQ[17]);

  hist_energy_ccd->Scale(1./EXPOSURE_SEQ[1]+EXPOSURE_SEQ[2]+EXPOSURE_SEQ[3]+EXPOSURE_SEQ[4]+EXPOSURE_SEQ[5]+EXPOSURE_SEQ[6]+EXPOSURE_SEQ[7]+EXPOSURE_SEQ[8]+EXPOSURE_SEQ[9]+EXPOSURE_SEQ[10]+EXPOSURE_SEQ[11]+EXPOSURE_SEQ[12]+EXPOSURE_SEQ[13]+EXPOSURE_SEQ[14]+EXPOSURE_SEQ[15]+EXPOSURE_SEQ[16]+EXPOSURE_SEQ[17]);

  hist_energy_anode->Scale(1./EXPOSURE_SEQ[1]+EXPOSURE_SEQ[2]+EXPOSURE_SEQ[3]+EXPOSURE_SEQ[4]+EXPOSURE_SEQ[5]+EXPOSURE_SEQ[6]+EXPOSURE_SEQ[7]+EXPOSURE_SEQ[8]+EXPOSURE_SEQ[9]+EXPOSURE_SEQ[10]+EXPOSURE_SEQ[11]+EXPOSURE_SEQ[12]+EXPOSURE_SEQ[13]+EXPOSURE_SEQ[14]+EXPOSURE_SEQ[15]+EXPOSURE_SEQ[16]+EXPOSURE_SEQ[17]);

  //scale after event loop, scale for each sequence and then scale for overall using 1/overall exposure

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
  hist_energy_mesh2->SetTitle("Scaled Mesh Energy Sequence 2");
  hist_energy_mesh2->SetXTitle("Mesh Energy / keV");
  hist_energy_mesh2->SetYTitle("Events/Second/Bin");
  hist_energy_mesh2->Draw();

  new TCanvas;
  hist_energy_mesh3->SetTitle("Scaled Mesh Energy Sequence 3");
  hist_energy_mesh3->SetXTitle("Mesh Energy / keV");
  hist_energy_mesh3->SetYTitle("Events/Second/Bin");
  hist_energy_mesh3->Draw();
  
  new TCanvas;
  hist_energy_mesh4->SetTitle("Scaled Mesh Energy Sequence 4");
  hist_energy_mesh4->SetXTitle("Mesh Energy / keV");
  hist_energy_mesh4->SetYTitle("Events/Second/Bin");
  hist_energy_mesh4->Draw();
  
  new TCanvas;
  hist_energy_mesh5->SetTitle("Scaled Mesh Energy Sequence 5");
  hist_energy_mesh5->SetXTitle("Mesh Energy / keV");
  hist_energy_mesh5->SetYTitle("Events/Second/Bin");
  hist_energy_mesh5->Draw();
  
  new TCanvas;
  hist_energy_mesh6->SetTitle("Scaled Mesh Energy Sequence 6");
  hist_energy_mesh6->SetXTitle("Mesh Energy / keV");
  hist_energy_mesh6->SetYTitle("Events/Second/Bin");
  hist_energy_mesh6->Draw();
  
  new TCanvas;
  hist_energy_mesh7->SetTitle("Scaled Mesh Energy Sequence 7");
  hist_energy_mesh7->SetXTitle("Mesh Energy / keV");
  hist_energy_mesh7->SetYTitle("Events/Second/Bin");
  hist_energy_mesh7->Draw();
  
  new TCanvas;
  hist_energy_mesh8->SetTitle("Scaled Mesh Energy Sequence 8");
  hist_energy_mesh8->SetXTitle("Mesh Energy / keV");
  hist_energy_mesh8->SetYTitle("Events/Second/Bin");
  hist_energy_mesh8->Draw();
  
  new TCanvas;
  hist_energy_mesh9->SetTitle("Scaled Mesh Energy Sequence 9");
  hist_energy_mesh9->SetXTitle("Mesh Energy / keV");
  hist_energy_mesh9->SetYTitle("Events/Second/Bin");
  hist_energy_mesh9->Draw();
  
  new TCanvas;
  hist_energy_mesh10->SetTitle("Scaled Mesh Energy Sequence 10");
  hist_energy_mesh10->SetXTitle("Mesh Energy / keV");
  hist_energy_mesh10->SetYTitle("Events/Second/Bin");
  hist_energy_mesh10->Draw();
  
  new TCanvas;
  hist_energy_mesh11->SetTitle("Scaled Mesh Energy Sequence 11");
  hist_energy_mesh11->SetXTitle("Mesh Energy / keV");
  hist_energy_mesh11->SetYTitle("Events/Second/Bin");
  hist_energy_mesh11->Draw();
  
  new TCanvas;
  hist_energy_mesh12->SetTitle("Scaled Mesh Energy Sequence 12");
  hist_energy_mesh12->SetXTitle("Mesh Energy / keV");
  hist_energy_mesh12->SetYTitle("Events/Second/Bin");
  hist_energy_mesh12->Draw();
  
  new TCanvas;
  hist_energy_mesh13->SetTitle("Scaled Mesh Energy Sequence 13");
  hist_energy_mesh13->SetXTitle("Mesh Energy / keV");
  hist_energy_mesh13->SetYTitle("Events/Second/Bin");
  hist_energy_mesh13->Draw();
  
  new TCanvas;
  hist_energy_mesh14->SetTitle("Scaled Mesh Energy Sequence 14");
  hist_energy_mesh14->SetXTitle("Mesh Energy / keV");
  hist_energy_mesh14->SetYTitle("Events/Second/Bin");
  hist_energy_mesh14->Draw();
  
  new TCanvas;
  hist_energy_mesh15->SetTitle("Scaled Mesh Energy Sequence 15");
  hist_energy_mesh15->SetXTitle("Mesh Energy / keV");
  hist_energy_mesh15->SetYTitle("Events/Second/Bin");
  hist_energy_mesh15->Draw();
  
  new TCanvas;
  hist_energy_mesh16->SetTitle("Scaled Mesh Energy Sequence 16");
  hist_energy_mesh16->SetXTitle("Mesh Energy / keV");
  hist_energy_mesh16->SetYTitle("Events/Second/Bin");
  hist_energy_mesh16->Draw();
  
  new TCanvas;
  hist_energy_mesh17->SetTitle("Scaled Mesh Energy Sequence 17");
  hist_energy_mesh17->SetXTitle("Mesh Energy / keV");
  hist_energy_mesh17->SetYTitle("Events/Second/Bin");
  hist_energy_mesh17->Draw();
  
  new TCanvas;
  hist_energy_mesh18->SetTitle("Scaled Mesh Energy Sequence 18");
  hist_energy_mesh18->SetXTitle("Mesh Energy / keV");
  hist_energy_mesh18->SetYTitle("Events/Second/Bin");
  hist_energy_mesh18->Draw();
  
  new TCanvas;
  hist_energy_mesh2->SetTitle("Mesh Energy");
  hist_energy_mesh2->SetXTitle("Mesh Energy / keV");
  hist_energy_mesh2->SetYTitle("Events/Second/Bin");
  hist_energy_mesh2->Draw();
  hist_energy_mesh3->SetLineColor(kRed);
  hist_energy_mesh3->Draw("SAME");
  hist_energy_mesh4->SetLineColor(kBlue);
  hist_energy_mesh4->Draw("SAME");
  hist_energy_mesh5->SetLineColor(kGreen);
  hist_energy_mesh5->Draw("SAME");
  hist_energy_mesh6->SetLineColor(kYellow);
  hist_energy_mesh6->Draw("SAME");
  hist_energy_mesh7->SetLineColor(kMagenta);
  hist_energy_mesh7->Draw("SAME");
  hist_energy_mesh8->SetLineColor(kTeal-6);
  hist_energy_mesh8->Draw("SAME");
  hist_energy_mesh9->SetLineColor(kOrange-3);
  hist_energy_mesh9->Draw("SAME");
  hist_energy_mesh10->SetLineColor(kCyan);
  hist_energy_mesh10->Draw("SAME");
  hist_energy_mesh11->SetLineColor(kSpring-6);
  hist_energy_mesh11->Draw("SAME");
  hist_energy_mesh12->SetLineColor(kAzure);
  hist_energy_mesh12->Draw("SAME");
  hist_energy_mesh13->SetLineColor(kViolet-6);
  hist_energy_mesh13->Draw("SAME");
  hist_energy_mesh14->SetLineColor(kPink-9);
  hist_energy_mesh14->Draw("SAME");
  hist_energy_mesh15->SetLineColor(kGray);
  hist_energy_mesh15->Draw("SAME");
  hist_energy_mesh16->SetLineColor(kYellow+2);
  hist_energy_mesh16->Draw("SAME");
  hist_energy_mesh17->SetLineColor(kViolet+1);
  hist_energy_mesh17->Draw("SAME");
  hist_energy_mesh18->SetLineColor(kRed-10);
  hist_energy_mesh18->Draw("SAME");

  new TCanvas;
  hist_energy_meshall->SetTitle("Mesh Energy");
  hist_energy_meshall->SetXTitle("Mesh Energy / keV");
  hist_energy_meshall->SetYTitle("Events/Second/Bin");
  hist_energy_meshall->Draw();

  new TCanvas;
  hist_energy_mesh->SetTitle("All Energies");
  hist_energy_mesh->SetXTitle("Energy / keV");
  hist_energy_mesh->SetYTitle("Events/Second/Bin");
  hist_energy_mesh->Draw();
  hist_energy_ccd->SetLineColor(kRed);
  hist_energy_ccd->Draw("SAME");
  hist_energy_anode->SetLineColor(kBlue);
  hist_energy_anode->Draw("SAME");
  leg_energy = new TLegend(0.6, 0.7, 0.89, 0.89);
  leg_energy->AddEntry(hist_energy_mesh, "Mesh Energy", "l");
  leg_energy->AddEntry(hist_energy_ccd, "CCD Energy", "l");
  leg_energy->AddEntry(hist_energy_anode, "Anode Energy", "l");
  leg_energy->Draw();
 
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  cout<<"Sequence 2 Mean Mesh Energy "<<hist_energy_mesh2->GetMean()<<endl;
  cout<<"Sequence 3 Mean Mesh Energy "<<hist_energy_mesh3->GetMean()<<endl;
  cout<<"Sequence 4 Mean Mesh Energy "<<hist_energy_mesh4->GetMean()<<endl;
  cout<<"Sequence 5 Mean Mesh Energy "<<hist_energy_mesh5->GetMean()<<endl;
  cout<<"Sequence 6 Mean Mesh Energy "<<hist_energy_mesh6->GetMean()<<endl;
  cout<<"Sequence 7 Mean Mesh Energy "<<hist_energy_mesh7->GetMean()<<endl;
  cout<<"Sequence 8 Mean Mesh Energy "<<hist_energy_mesh8->GetMean()<<endl;
  cout<<"Sequence 9 Mean Mesh Energy "<<hist_energy_mesh9->GetMean()<<endl;
  cout<<"Sequence 10 Mean Mesh Energy "<<hist_energy_mesh10->GetMean()<<endl;
  cout<<"Sequence 11 Mean Mesh Energy "<<hist_energy_mesh11->GetMean()<<endl;
  cout<<"Sequence 12 Mean Mesh Energy "<<hist_energy_mesh12->GetMean()<<endl;
  cout<<"Sequence 13 Mean Mesh Energy "<<hist_energy_mesh13->GetMean()<<endl;
  cout<<"Sequence 14 Mean Mesh Energy "<<hist_energy_mesh14->GetMean()<<endl;
  cout<<"Sequence 15 Mean Mesh Energy "<<hist_energy_mesh15->GetMean()<<endl;
  cout<<"Sequence 16 Mean Mesh Energy "<<hist_energy_mesh16->GetMean()<<endl;
  cout<<"Sequence 17 Mean Mesh Energy "<<hist_energy_mesh17->GetMean()<<endl;
  cout<<"Sequence 18 Mean Mesh Energy "<<hist_energy_mesh18->GetMean()<<endl;
  cout<<"Overall Mean Mesh Energy "<<hist_energy_meshall->GetMean()<<endl;
  
}

double fitPoly12Percent(double E)
{
  return sqrt(0.9741*0.9741*E*E+113.835*113.835)-62.675;
  //return -11.50+39.02*E+0.03459*E*E;                                                    
  //return 469 + 10*E + 0.548 *E*E -3.39e-3*E*E*E+7.5e-6*E*E*E*E;          
  
}               

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
