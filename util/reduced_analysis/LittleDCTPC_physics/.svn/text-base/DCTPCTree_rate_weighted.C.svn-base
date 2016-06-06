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
  
  //Scaled to Exposure Neutron Number vs Time plots
  TH2D *hist_time_2=new TH2D("Rate", "Rate", 100, 1374000000, 1386000000, 600, 0, 600);
  TH2D *hist_time_3=new TH2D("Rate", "Rate", 100, 1374000000, 1386000000, 600, 0, 600);
  TH2D *hist_time_4=new TH2D("Rate", "Rate", 100, 1374000000, 1386000000, 600, 0, 600);
  TH2D *hist_time_5=new TH2D("Rate", "Rate", 100, 1374000000, 1386000000, 600, 0, 600);
  TH2D *hist_time_6=new TH2D("Rate", "Rate", 100, 1374000000, 1386000000, 600, 0, 600);
  TH2D *hist_time_7=new TH2D("Rate", "Rate", 100, 1374000000, 1386000000, 600, 0, 600);
  TH2D *hist_time_8=new TH2D("Rate", "Rate", 100, 1374000000, 1386000000, 600, 0, 600);
  TH2D *hist_time_9=new TH2D("Rate", "Rate", 100, 1374000000, 1386000000, 600, 0, 600);
  TH2D *hist_time_10=new TH2D("Rate", "Rate", 100, 1374000000, 1386000000, 600, 0, 600);
  TH2D *hist_time_11=new TH2D("Rate", "Rate", 100, 1374000000, 1386000000, 600, 0, 600);
  TH2D *hist_time_12=new TH2D("Rate", "Rate", 100, 1374000000, 1386000000, 600, 0, 600);
  TH2D *hist_time_13=new TH2D("Rate", "Rate", 100, 1374000000, 1386000000, 600, 0, 600);
  TH2D *hist_time_14=new TH2D("Rate", "Rate", 100, 1374000000, 1386000000, 600, 0, 600);
  TH2D *hist_time_15=new TH2D("Rate", "Rate", 100, 1374000000, 1386000000, 600, 0, 600);
  TH2D *hist_time_16=new TH2D("Rate", "Rate", 100, 1374000000, 1386000000, 600, 0, 600);
  TH2D *hist_time_17=new TH2D("Rate", "Rate", 100, 1374000000, 1386000000, 600, 0, 600);
  TH2D *hist_time_18=new TH2D("Rate", "Rate", 100, 1374000000, 1386000000, 600, 0, 600);

   //Scaled to Exposure Neutron Number vs Run Number plots
  TH2D *hist_neutronsperrun_2=new TH2D("Neutrons per Run", "Neutrons per Run", 800, 9000, 17000, 600, 0, 600);
  TH2D *hist_neutronsperrun_3=new TH2D("Neutrons per Run", "Neutrons per Run", 800, 9000, 17000, 600, 0, 600);
  TH2D *hist_neutronsperrun_4=new TH2D("Neutrons per Run", "Neutrons per Run", 800, 9000, 17000, 600, 0, 600);
  TH2D *hist_neutronsperrun_5=new TH2D("Neutrons per Run", "Neutrons per Run", 800, 9000, 17000, 600, 0, 600);
  TH2D *hist_neutronsperrun_6=new TH2D("Neutrons per Run", "Neutrons per Run", 800, 9000, 17000, 600, 0, 600);
  TH2D *hist_neutronsperrun_7=new TH2D("Neutrons per Run", "Neutrons per Run", 800, 9000, 17000, 600, 0, 600);
  TH2D *hist_neutronsperrun_8=new TH2D("Neutrons per Run", "Neutrons per Run", 800, 9000, 17000, 600, 0, 600);
  TH2D *hist_neutronsperrun_9=new TH2D("Neutrons per Run", "Neutrons per Run", 800, 9000, 17000, 600, 0, 600);
  TH2D *hist_neutronsperrun_10=new TH2D("Neutrons per Run", "Neutrons per Run", 800, 9000, 17000, 600, 0, 600);
  TH2D *hist_neutronsperrun_11=new TH2D("Neutrons per Run", "Neutrons per Run", 800, 9000, 17000, 600, 0, 600);
  TH2D *hist_neutronsperrun_12=new TH2D("Neutrons per Run", "Neutrons per Run", 800, 9000, 17000, 600, 0, 600);
  TH2D *hist_neutronsperrun_13=new TH2D("Neutrons per Run", "Neutrons per Run", 800, 9000, 17000, 600, 0, 600);
  TH2D *hist_neutronsperrun_14=new TH2D("Neutrons per Run", "Neutrons per Run", 800, 9000, 17000, 600, 0, 600);
  TH2D *hist_neutronsperrun_15=new TH2D("Neutrons per Run", "Neutrons per Run", 800, 9000, 17000, 600, 0, 600);
  TH2D *hist_neutronsperrun_16=new TH2D("Neutrons per Run", "Neutrons per Run", 800, 9000, 17000, 600, 0, 600);
  TH2D *hist_neutronsperrun_17=new TH2D("Neutrons per Run", "Neutrons per Run", 800, 9000, 17000, 600, 0, 600);
  TH2D *hist_neutronsperrun_18=new TH2D("Neutrons per Run", "Neutrons per Run", 800, 9000, 17000, 600, 0, 600);

  //Scaled Neutron Number Plots
  TH2D *hist_time_counter=new TH2D("Neutrons vs Time", "Neutrons vs Time", 1000, 1374000000, 1386000000, 600, 0, 600);
  TH2D *hist_run_counter=new TH2D("Neutrons vs Run", "Neutrons vs Run", 8000, 9000, 17000, 600, 0, 600);
  
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
  double counter_time=0;
  double counter_run=0;
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for (int event = 0; event</*200000*/nentries; event++)
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
      
      if(tracklength>72)
	{
	  cutnum=9;
	  continue;
	}
    

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	      
         
      //above I am filling the histograms with data that has passed all the cuts - hopefully neutrons!
         
      cutnum=0;
      neutronnumber++;
      counter_time=counter_time+1+(1./EXPOSURE_SEQ[aStep.SequenceNum-1]);
      counter_run=counter_run+1+(1./EXPOSURE_SEQ[aStep.SequenceNum-1]);

      if(aStep.SequenceNum==2)
	{
	  hist_time_2->Fill(aStep.Timenow_sec, neutronnumber, 1./345130.);//1./EXPOSURE_SEQ[1]);
	  hist_neutronsperrun_2->Fill(aStep.RunNum, neutronnumber, 1./345130.);
	}

      if(aStep.SequenceNum==3)
	{
	  hist_time_3->Fill(aStep.Timenow_sec, neutronnumber, 1./87757.);// 1./EXPOSURE_SEQ[2]);
	  hist_neutronsperrun_3->Fill(aStep.RunNum, neutronnumber, 1./87757.);
	}

      if(aStep.SequenceNum==4)
	{
	  hist_time_4->Fill(aStep.Timenow_sec, neutronnumber, 1./248539.);//1./EXPOSURE_SEQ[3]);
	  hist_neutronsperrun_4->Fill(aStep.RunNum, neutronnumber, 1./248539.);
	}

      if(aStep.SequenceNum==5)
	{
	  hist_time_5->Fill(aStep.Timenow_sec, neutronnumber, 1./367407.);// 1./EXPOSURE_SEQ[4]);
	  hist_neutronsperrun_5->Fill(aStep.RunNum, neutronnumber, 1./367407.);
	}

      if(aStep.SequenceNum==6)
	{
	  hist_time_6->Fill(aStep.Timenow_sec, neutronnumber, 1./364018.);//1./EXPOSURE_SEQ[5]);
	  hist_neutronsperrun_6->Fill(aStep.RunNum, neutronnumber, 1./364018.);
	}

      if(aStep.SequenceNum==7)
	{
	  hist_time_7->Fill(aStep.Timenow_sec, neutronnumber, 1./89584.);//1./EXPOSURE_SEQ[6]);
	  hist_neutronsperrun_7->Fill(aStep.RunNum, neutronnumber, 1./89584.);
	}

      if(aStep.SequenceNum==8)
	{
	  hist_time_8->Fill(aStep.Timenow_sec, neutronnumber, 1./251442.);//1./EXPOSURE_SEQ[7]);
	  hist_neutronsperrun_8->Fill(aStep.RunNum, neutronnumber, 1./251442.);
	}

      if(aStep.SequenceNum==9)
	{
	  hist_time_9->Fill(aStep.Timenow_sec, neutronnumber, 1./324605.);// 1./EXPOSURE_SEQ[8]);
	  hist_neutronsperrun_9->Fill(aStep.RunNum, neutronnumber, 1./324605.);
	}

      if(aStep.SequenceNum==10)
	{
	  hist_time_10->Fill(aStep.Timenow_sec, neutronnumber, 1./358127);//1./EXPOSURE_SEQ[9]);
	  hist_neutronsperrun_10->Fill(aStep.RunNum, neutronnumber, 1./358127.);
	}

      if(aStep.SequenceNum==11)
	{
	  hist_time_11->Fill(aStep.Timenow_sec, neutronnumber, 1./359194.);//1./EXPOSURE_SEQ[10]);
	  hist_neutronsperrun_11->Fill(aStep.RunNum, neutronnumber, 1./359194.);
	}

      if(aStep.SequenceNum==12)
	{
	  hist_time_12->Fill(aStep.Timenow_sec, neutronnumber, 1./367272.);//1./EXPOSURE_SEQ[11]);
	  hist_neutronsperrun_12->Fill(aStep.RunNum, neutronnumber, 1./367272.);
	}

      if(aStep.SequenceNum==13)
	{
	  hist_time_13->Fill(aStep.Timenow_sec, neutronnumber, 1./365276);//1./EXPOSURE_SEQ[12]);
	  hist_neutronsperrun_13->Fill(aStep.RunNum, neutronnumber, 1./365276.);
	}

      if(aStep.SequenceNum==14)
	{
	  hist_time_14->Fill(aStep.Timenow_sec, neutronnumber, 1./354349.);//1./EXPOSURE_SEQ[13]);
	  hist_neutronsperrun_14->Fill(aStep.RunNum, neutronnumber, 1./354349.);
	}

      if(aStep.SequenceNum==15)
	{
	  hist_time_15->Fill(aStep.Timenow_sec, neutronnumber, 1./370056.);//1./EXPOSURE_SEQ[14]);
	  hist_neutronsperrun_15->Fill(aStep.RunNum, neutronnumber, 1./370056.);
	}

      if(aStep.SequenceNum==16)
	{
	  hist_time_16->Fill(aStep.Timenow_sec, neutronnumber, 1./125736.);//1./EXPOSURE_SEQ[15]);
	  hist_neutronsperrun_16->Fill(aStep.RunNum, neutronnumber, 1./125736.);
	}

      if(aStep.SequenceNum==17)
	{
	  hist_time_17->Fill(aStep.Timenow_sec, neutronnumber, 1./103122.);//1./EXPOSURE_SEQ[16]);
	  hist_neutronsperrun_17->Fill(aStep.RunNum, neutronnumber, 1./103122.);
	}

      if(aStep.SequenceNum==18)
	{
	  hist_time_18->Fill(aStep.Timenow_sec, neutronnumber, 1./1111687.);
	  hist_neutronsperrun_18->Fill(aStep.RunNum, neutronnumber, 1./1111687);
	}

      hist_time_counter->Fill(aStep.Timenow_sec, counter_time);
      hist_run_counter->Fill(aStep.RunNum, counter_run);
      
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
 
 
  new TCanvas;
  hist_time_2->SetTitle("Neutron Number vs Time Weighted by Exposure");
  hist_time_2->SetYTitle("Neutron Number");
  hist_time_2->SetXTitle("Unix Time / sec");
  hist_time_2->SetMarkerStyle(7);
  hist_time_2->Draw("COLZ");
  hist_time_3->SetMarkerStyle(7);
  hist_time_3->Draw("COLZ SAME");
  hist_time_4->SetMarkerStyle(7);
  hist_time_4->Draw("SAME COLZ");
  hist_time_5->SetMarkerStyle(7);
  hist_time_5->Draw("SAME COLZ");
  hist_time_6->SetMarkerStyle(7);
  hist_time_6->Draw("SAME COLZ");
  hist_time_7->SetMarkerStyle(7);
  hist_time_7->Draw("SAME COLZ");
  hist_time_8->SetMarkerStyle(7);
  hist_time_8->Draw("SAME COLZ");
  hist_time_9->SetMarkerStyle(7);
  hist_time_9->Draw("SAME COLZ");
  hist_time_10->SetMarkerStyle(7);
  hist_time_10->Draw("SAME COLZ");
  hist_time_11->SetMarkerStyle(7);
  hist_time_11->Draw("SAME COLZ");
  hist_time_12->SetMarkerStyle(7);
  hist_time_12->Draw("SAME COLZ");
  hist_time_13->SetMarkerStyle(7);
  hist_time_13->Draw("SAME COLZ");
  hist_time_14->SetMarkerStyle(7);
  hist_time_14->Draw("SAME COLZ");
  hist_time_15->SetMarkerStyle(7);
  hist_time_15->Draw("SAME COLZ");
  hist_time_16->SetMarkerStyle(7);
  hist_time_16->Draw("SAME COLZ");
  hist_time_17->SetMarkerStyle(7);
  hist_time_17->Draw("SAME COLZ");
  hist_time_18->SetMarkerStyle(7);
  hist_time_18->Draw("SAME COLZ");
  
  new TCanvas;
  hist_neutronsperrun_2->SetTitle("Number of Neutrons per Run Weighted By Exposure");
  hist_neutronsperrun_2->SetYTitle("Neutron Number");
  hist_neutronsperrun_2->SetXTitle("Run Number");
  hist_neutronsperrun_2->SetMarkerStyle(7);
  hist_neutronsperrun_2->Draw("COLZ");
  hist_neutronsperrun_3->SetMarkerStyle(7);
  hist_neutronsperrun_3->Draw("SAME COLZ");
  hist_neutronsperrun_4->SetMarkerStyle(7);
  hist_neutronsperrun_4->Draw("SAME COLZ");
  hist_neutronsperrun_5->SetMarkerStyle(7);
  hist_neutronsperrun_5->Draw("SAME COLZ");
  hist_neutronsperrun_6->SetMarkerStyle(7);
  hist_neutronsperrun_6->Draw("SAME COLZ");
  hist_neutronsperrun_7->SetMarkerStyle(7);
  hist_neutronsperrun_7->Draw("SAME COLZ");
  hist_neutronsperrun_8->SetMarkerStyle(7);
  hist_neutronsperrun_8->Draw("SAME COLZ");
  hist_neutronsperrun_9->SetMarkerStyle(7);
  hist_neutronsperrun_9->Draw("SAME COLZ");
  hist_neutronsperrun_10->SetMarkerStyle(7);
  hist_neutronsperrun_10->Draw("SAME COLZ");
  hist_neutronsperrun_11->SetMarkerStyle(7);
  hist_neutronsperrun_11->Draw("SAME COLZ");
  hist_neutronsperrun_12->SetMarkerStyle(7);
  hist_neutronsperrun_12->Draw("SAME COLZ");
  hist_neutronsperrun_13->SetMarkerStyle(7);
  hist_neutronsperrun_13->Draw("SAME COLZ");
  hist_neutronsperrun_14->SetMarkerStyle(7);
  hist_neutronsperrun_14->Draw("SAME COLZ");
  hist_neutronsperrun_15->SetMarkerStyle(7);
  hist_neutronsperrun_15->Draw("SAME COLZ");
  hist_neutronsperrun_16->SetMarkerStyle(7);
  hist_neutronsperrun_16->Draw("SAME COLZ");
  hist_neutronsperrun_17->SetMarkerStyle(7);
  hist_neutronsperrun_17->Draw("SAME COLZ");
  hist_neutronsperrun_18->SetMarkerStyle(7);
  hist_neutronsperrun_18->Draw("SAME COLZ");
 
  new TCanvas;
  hist_time_counter->SetTitle("Neutron Number vs Unix Time Scaled");
  hist_time_counter->SetYTitle("Neutron Number");
  hist_time_counter->SetXTitle("Unix Time / sec");
  hist_time_counter->SetMarkerStyle(2);
  hist_time_counter->Draw();

  new TCanvas;
  hist_run_counter->SetTitle("Neutron Number vs Run Number Scaled");
  hist_run_counter->SetYTitle("Neutron Number");
  hist_run_counter->SetXTitle("Run Number");
  hist_run_counter->SetMarkerStyle(2);
  hist_run_counter->Draw();

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  cout<<"Number of events cut= "<<nentries-neutronnumber<<endl;
 
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
  
