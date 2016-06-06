#define DCTPC_runtree_cxx
#include "../input_files/DCTPC_runtree.h"
#include "../input_files/Calibration_constants_LittleDCTPC.h";
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TRandom3.h"
#include "TGraph.h"
#include "TMath.h"
#include "TArc.h"
#include "TVirtualFitter.h"
TGraph *gr;

//This macro is used to find the exposure in seconds for every sequence

void DCTPC_runtree::Loop()
{
  
  TStopwatch timer;
  gROOT->SetStyle("Plain");
  gStyle->SetEndErrorSize(3);
  gStyle->SetPalette(1,0);
  gStyle->SetLineWidth(2);
  gStyle->SetHistLineWidth(2);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptFit(kFALSE);
  TH1::AddDirectory(false);
  
  TFile *outtree = new TFile("$LittleDCTPC_physics_input");
  TTree *dctreepc = (TTree*)outtree->Get("dctpc_runinfo");
  
  DCTPC_runtree aStep(dctreepc);
  
  Long64_t nentries = dctreepc->GetEntries();
  cout << "Number of Entries: " << nentries << endl;
  
  TH2D *hist_triggers_vs_time2=new TH2D("No. Trigs vs Time", "No. Trigs vs Time", 1000, 1374000000, 1386000000, 2000, 0, 10000);
  TH2D *hist_triggers_vs_time3=new TH2D("No. Trigs vs Time", "No. Trigs vs Time", 1000, 1374000000, 1386000000, 2000, 0, 10000);
  TH2D *hist_triggers_vs_time4=new TH2D("No. Trigs vs Time", "No. Trigs vs Time", 1000, 1374000000, 1386000000, 2000, 0, 10000);
  TH2D *hist_triggers_vs_time5=new TH2D("No. Trigs vs Time", "No. Trigs vs Time", 1000, 1374000000, 1386000000, 2000, 0, 10000);
  TH2D *hist_triggers_vs_time6=new TH2D("No. Trigs vs Time", "No. Trigs vs Time", 1000, 1374000000, 1386000000, 2000, 0, 10000);
  TH2D *hist_triggers_vs_time7=new TH2D("No. Trigs vs Time", "No. Trigs vs Time", 1000, 1374000000, 1386000000, 2000, 0, 10000);
  TH2D *hist_triggers_vs_time8=new TH2D("No. Trigs vs Time", "No. Trigs vs Time", 1000, 1374000000, 1386000000, 2000, 0, 10000);
  TH2D *hist_triggers_vs_time9=new TH2D("No. Trigs vs Time", "No. Trigs vs Time", 1000, 1374000000, 1386000000, 2000, 0, 10000);
  TH2D *hist_triggers_vs_time10=new TH2D("No. Trigs vs Time", "No. Trigs vs Time", 1000, 1374000000, 1386000000, 2000, 0, 10000);
  TH2D *hist_triggers_vs_time11=new TH2D("No. Trigs vs Time", "No. Trigs vs Time", 1000, 1374000000, 1386000000, 2000, 0, 10000);
  TH2D *hist_triggers_vs_time12=new TH2D("No. Trigs vs Time", "No. Trigs vs Time", 1000, 1374000000, 1386000000, 2000, 0, 10000);
  TH2D *hist_triggers_vs_time13=new TH2D("No. Trigs vs Time", "No. Trigs vs Time", 1000, 1374000000, 1386000000, 2000, 0, 10000);
  TH2D *hist_triggers_vs_time14=new TH2D("No. Trigs vs Time", "No. Trigs vs Time", 1000, 1374000000, 1386000000, 2000, 0, 10000);
  TH2D *hist_triggers_vs_time15=new TH2D("No. Trigs vs Time", "No. Trigs vs Time", 1000, 1374000000, 1386000000, 2000, 0, 10000);
  TH2D *hist_triggers_vs_time16=new TH2D("No. Trigs vs Time", "No. Trigs vs Time", 1000, 1374000000, 1386000000, 2000, 0, 10000);
  TH2D *hist_triggers_vs_time17=new TH2D("No. Trigs vs Time", "No. Trigs vs Time", 1000, 1374000000, 1386000000, 2000, 0, 10000);
  TH2D *hist_triggers_vs_time18=new TH2D("No. Trigs vs Time", "No. Trigs vs Time", 1000, 1374000000, 1386000000, 2000, 0, 10000);
  
  TH2D *hist_tracks_vs_time2=new TH2D("No. Tracks vs Time", "No. Tracks vs Time", 100, 1374000000, 1386000000, 700, 0, 700);
  TH2D *hist_tracks_vs_time3=new TH2D("No. Tracks vs Time", "No. Tracks vs Time", 100, 1374000000, 1386000000, 700, 0, 700);
  TH2D *hist_tracks_vs_time4=new TH2D("No. Tracks vs Time", "No. Tracks vs Time", 100, 1374000000, 1386000000, 700, 0, 700);
  TH2D *hist_tracks_vs_time5=new TH2D("No. Tracks vs Time", "No. Tracks vs Time", 100, 1374000000, 1386000000, 700, 0, 700);
  TH2D *hist_tracks_vs_time6=new TH2D("No. Tracks vs Time", "No. Tracks vs Time", 100, 1374000000, 1386000000, 700, 0, 700);
  TH2D *hist_tracks_vs_time7=new TH2D("No. Tracks vs Time", "No. Tracks vs Time", 100, 1374000000, 1386000000, 700, 0, 700);
  TH2D *hist_tracks_vs_time8=new TH2D("No. Tracks vs Time", "No. Tracks vs Time", 100, 1374000000, 1386000000, 700, 0, 700);
  TH2D *hist_tracks_vs_time9=new TH2D("No. Tracks vs Time", "No. Tracks vs Time", 100, 1374000000, 1386000000, 700, 0, 700);
  TH2D *hist_tracks_vs_time10=new TH2D("No. Tracks vs Time", "No. Tracks vs Time", 100, 1374000000, 1386000000, 700, 0, 700);
  TH2D *hist_tracks_vs_time11=new TH2D("No. Tracks vs Time", "No. Tracks vs Time", 100, 1374000000, 1386000000, 700, 0, 700);
  TH2D *hist_tracks_vs_time12=new TH2D("No. Tracks vs Time", "No. Tracks vs Time", 100, 1374000000, 1386000000, 700, 0, 700);
  TH2D *hist_tracks_vs_time13=new TH2D("No. Tracks vs Time", "No. Tracks vs Time", 100, 1374000000, 1386000000, 700, 0, 700);
  TH2D *hist_tracks_vs_time14=new TH2D("No. Tracks vs Time", "No. Tracks vs Time", 100, 1374000000, 1386000000, 700, 0, 700);
  TH2D *hist_tracks_vs_time15=new TH2D("No. Tracks vs Time", "No. Tracks vs Time", 100, 1374000000, 1386000000, 700, 0, 700);
  TH2D *hist_tracks_vs_time16=new TH2D("No. Tracks vs Time", "No. Tracks vs Time", 100, 1374000000, 1386000000, 700, 0, 700);
  TH2D *hist_tracks_vs_time17=new TH2D("No. Tracks vs Time", "No. Tracks vs Time", 100, 1374000000, 1386000000, 700, 0, 700);
  TH2D *hist_tracks_vs_time18=new TH2D("No. Tracks vs Time", "No. Tracks vs Time", 100, 1374000000, 1386000000, 700, 0, 700);
  
  for (int event = 0; event<nentries; event++)
    {
      aStep.GetEntry(event);
      
      if(aStep.RunNum>=9535 && aStep.RunNum<=9906)
	{
	  hist_triggers_vs_time2->Fill(aStep.Time_startofrun_sec, aStep.Totaltrigs);
	  hist_tracks_vs_time2->Fill(aStep.Time_startofrun_sec, aStep.Totaltracks);
	}

      if(aStep.RunNum>=9907 && aStep.RunNum<=10001)
	{
	  hist_triggers_vs_time3->Fill(aStep.Time_startofrun_sec, aStep.Totaltrigs);
	  hist_tracks_vs_time3->Fill(aStep.Time_startofrun_sec, aStep.Totaltracks);
	}

      if(aStep.RunNum>=10004 && aStep.RunNum<=10273)
	{
	  hist_triggers_vs_time4->Fill(aStep.Time_startofrun_sec, aStep.Totaltrigs);
	  hist_tracks_vs_time4->Fill(aStep.Time_startofrun_sec, aStep.Totaltracks);
	}
      
      if(aStep.RunNum>=10274 && aStep.RunNum<=10671)
	{
	  hist_triggers_vs_time5->Fill(aStep.Time_startofrun_sec, aStep.Totaltrigs);
	  hist_tracks_vs_time5->Fill(aStep.Time_startofrun_sec, aStep.Totaltracks);
	}	

      if(aStep.RunNum>=10672 && aStep.RunNum<=11064)
	{
	  hist_triggers_vs_time6->Fill(aStep.Time_startofrun_sec, aStep.Totaltrigs);
	  hist_tracks_vs_time6->Fill(aStep.Time_startofrun_sec, aStep.Totaltracks);
	}

      if(aStep.RunNum>=11065 && aStep.RunNum<=11161)
	{
	  hist_triggers_vs_time7->Fill(aStep.Time_startofrun_sec, aStep.Totaltrigs);
	  hist_tracks_vs_time7->Fill(aStep.Time_startofrun_sec, aStep.Totaltracks);
	}

      if(aStep.RunNum>=11162 && aStep.RunNum<=11435)
	{
	  hist_triggers_vs_time8->Fill(aStep.Time_startofrun_sec, aStep.Totaltrigs);
	  hist_tracks_vs_time8->Fill(aStep.Time_startofrun_sec, aStep.Totaltracks);
	}

      if(aStep.RunNum>=11436 && aStep.RunNum<=11788)
	{
	  hist_triggers_vs_time9->Fill(aStep.Time_startofrun_sec, aStep.Totaltrigs);
	  hist_tracks_vs_time9->Fill(aStep.Time_startofrun_sec, aStep.Totaltracks);
	}

      if(aStep.RunNum>=11789 && aStep.RunNum<=12175)
	{
	  hist_triggers_vs_time10->Fill(aStep.Time_startofrun_sec, aStep.Totaltrigs);
	  hist_tracks_vs_time10->Fill(aStep.Time_startofrun_sec, aStep.Totaltracks);
	}

      if(aStep.RunNum>=12176 && aStep.RunNum<=12561)
	{
	  hist_triggers_vs_time11->Fill(aStep.Time_startofrun_sec, aStep.Totaltrigs);
	  hist_tracks_vs_time11->Fill(aStep.Time_startofrun_sec, aStep.Totaltracks);
	}

      if(aStep.RunNum>=12562 && aStep.RunNum<=12947)
	{
	  hist_triggers_vs_time12->Fill(aStep.Time_startofrun_sec, aStep.Totaltrigs);
	  hist_tracks_vs_time12->Fill(aStep.Time_startofrun_sec, aStep.Totaltracks);
	}

      if(aStep.RunNum>=12948 && aStep.RunNum<=13340)
	{
	  hist_triggers_vs_time13->Fill(aStep.Time_startofrun_sec, aStep.Totaltrigs);
	  hist_tracks_vs_time13->Fill(aStep.Time_startofrun_sec, aStep.Totaltracks);
	}

      if(aStep.RunNum>=13341 && aStep.RunNum<=13723)
	{
	  hist_triggers_vs_time14->Fill(aStep.Time_startofrun_sec, aStep.Totaltrigs);
	  hist_tracks_vs_time14->Fill(aStep.Time_startofrun_sec, aStep.Totaltracks);
	}

      if(aStep.RunNum>=13724 && aStep.RunNum<=14120)
	{
	  hist_triggers_vs_time15->Fill(aStep.Time_startofrun_sec, aStep.Totaltrigs);
	  hist_tracks_vs_time15->Fill(aStep.Time_startofrun_sec, aStep.Totaltracks);
	}

      if(aStep.RunNum>=14121 && aStep.RunNum<=14256)
	{
	  hist_triggers_vs_time16->Fill(aStep.Time_startofrun_sec, aStep.Totaltrigs);
	  hist_tracks_vs_time16->Fill(aStep.Time_startofrun_sec, aStep.Totaltracks);
	}

      if(aStep.RunNum>=14262 && aStep.RunNum<=14374)
	{
	  hist_triggers_vs_time17->Fill(aStep.Time_startofrun_sec, aStep.Totaltrigs);
	  hist_tracks_vs_time17->Fill(aStep.Time_startofrun_sec, aStep.Totaltracks);
	}

      if(aStep.RunNum>=14375 && aStep.RunNum<=15579)
	{
	  hist_triggers_vs_time18->Fill(aStep.Time_startofrun_sec, aStep.Totaltrigs);
	  hist_tracks_vs_time18->Fill(aStep.Time_startofrun_sec, aStep.Totaltracks);
	}

      
    }
  
  new TCanvas;
  hist_triggers_vs_time2->SetTitle("No. of Triggers vs Time");
  hist_triggers_vs_time2->SetXTitle("Unix Time / sec");
  hist_triggers_vs_time2->SetYTitle("No. of Triggers");
  hist_triggers_vs_time2->SetMarkerStyle(2);
  hist_triggers_vs_time2->Draw();
  hist_triggers_vs_time3->SetMarkerColor(2);
  hist_triggers_vs_time3->SetMarkerStyle(2);
  hist_triggers_vs_time3->Draw("SAME");
  hist_triggers_vs_time4->SetMarkerStyle(2);
  hist_triggers_vs_time4->SetMarkerColor(3);
  hist_triggers_vs_time4->Draw("SAME");
  hist_triggers_vs_time5->SetMarkerStyle(2);
  hist_triggers_vs_time5->SetMarkerColor(4);
  hist_triggers_vs_time5->Draw("SAME");
  hist_triggers_vs_time6->SetMarkerStyle(2);
  hist_triggers_vs_time6->SetMarkerColor(5);
  hist_triggers_vs_time6->Draw("SAME");
  hist_triggers_vs_time7->SetMarkerStyle(2);
  hist_triggers_vs_time7->SetMarkerColor(6);
  hist_triggers_vs_time7->Draw("SAME");
  hist_triggers_vs_time8->SetMarkerStyle(2);
  hist_triggers_vs_time8->SetMarkerColor(7);
  hist_triggers_vs_time8->Draw("SAME");
  hist_triggers_vs_time9->SetMarkerStyle(2);
  hist_triggers_vs_time9->SetMarkerColor(8);
  hist_triggers_vs_time9->Draw("SAME");
  hist_triggers_vs_time10->SetMarkerStyle(2);
  hist_triggers_vs_time10->SetMarkerColor(9);
  hist_triggers_vs_time10->Draw("SAME");
  hist_triggers_vs_time11->SetMarkerStyle(5);
  hist_triggers_vs_time11->SetMarkerColor(1);
  hist_triggers_vs_time11->Draw("SAME");
  hist_triggers_vs_time12->SetMarkerStyle(5);
  hist_triggers_vs_time12->SetMarkerColor(2);
  hist_triggers_vs_time12->Draw("SAME");
  hist_triggers_vs_time13->SetMarkerStyle(5);
  hist_triggers_vs_time13->SetMarkerColor(3);
  hist_triggers_vs_time13->Draw("SAME");
  hist_triggers_vs_time14->SetMarkerStyle(5);
  hist_triggers_vs_time14->SetMarkerColor(4);
  hist_triggers_vs_time14->Draw("SAME");
  hist_triggers_vs_time15->SetMarkerStyle(5);
  hist_triggers_vs_time15->SetMarkerColor(5);
  hist_triggers_vs_time15->Draw("SAME");
  hist_triggers_vs_time16->SetMarkerStyle(5);
  hist_triggers_vs_time16->SetMarkerColor(6);
  hist_triggers_vs_time16->Draw("SAME");
  hist_triggers_vs_time17->SetMarkerStyle(5);
  hist_triggers_vs_time17->SetMarkerColor(7);
  hist_triggers_vs_time17->Draw("SAME");
  hist_triggers_vs_time18->SetMarkerStyle(5);
  hist_triggers_vs_time18->SetMarkerColor(8);
  hist_triggers_vs_time18->Draw("SAME");

   new TCanvas;
  hist_tracks_vs_time2->SetTitle("No. of Tracks vs Time");
  hist_tracks_vs_time2->SetXTitle("Unix Time / sec");
  hist_tracks_vs_time2->SetYTitle("No. of Tracks");
  hist_tracks_vs_time2->SetMarkerStyle(2);
  hist_tracks_vs_time2->Draw();
  hist_tracks_vs_time3->SetMarkerColor(2);
  hist_tracks_vs_time3->SetMarkerStyle(2);
  hist_tracks_vs_time3->Draw("SAME");
  hist_tracks_vs_time4->SetMarkerColor(3);
  hist_tracks_vs_time4->SetMarkerStyle(2);
  hist_tracks_vs_time4->Draw("SAME");
  hist_tracks_vs_time5->SetMarkerColor(4);
  hist_tracks_vs_time5->SetMarkerStyle(2);
  hist_tracks_vs_time5->Draw("SAME");
  hist_tracks_vs_time6->SetMarkerColor(5);
  hist_tracks_vs_time6->SetMarkerStyle(2);
  hist_tracks_vs_time6->Draw("SAME");
  hist_tracks_vs_time7->SetMarkerColor(6);
  hist_tracks_vs_time7->SetMarkerStyle(2);
  hist_tracks_vs_time7->Draw("SAME");
  hist_tracks_vs_time8->SetMarkerColor(7);
  hist_tracks_vs_time8->SetMarkerStyle(2);
  hist_tracks_vs_time8->Draw("SAME");
  hist_tracks_vs_time9->SetMarkerColor(8);
  hist_tracks_vs_time9->SetMarkerStyle(2);
  hist_tracks_vs_time9->Draw("SAME");
  hist_tracks_vs_time10->SetMarkerColor(9);
  hist_tracks_vs_time10->SetMarkerStyle(2);
  hist_tracks_vs_time10->Draw("SAME");
  hist_tracks_vs_time11->SetMarkerColor(1);
  hist_tracks_vs_time11->SetMarkerStyle(5);
  hist_tracks_vs_time11->Draw("SAME");
  hist_tracks_vs_time12->SetMarkerColor(2);
  hist_tracks_vs_time12->SetMarkerStyle(5);
  hist_tracks_vs_time12->Draw("SAME");
  hist_tracks_vs_time13->SetMarkerColor(3);
  hist_tracks_vs_time13->SetMarkerStyle(5);
  hist_tracks_vs_time13->Draw("SAME");
  hist_tracks_vs_time14->SetMarkerColor(4);
  hist_tracks_vs_time14->SetMarkerStyle(5);
  hist_tracks_vs_time14->Draw("SAME");
  hist_tracks_vs_time15->SetMarkerColor(5);
  hist_tracks_vs_time15->SetMarkerStyle(5);
  hist_tracks_vs_time15->Draw("SAME");
  hist_tracks_vs_time16->SetMarkerColor(6);
  hist_tracks_vs_time16->SetMarkerStyle(5);
  hist_tracks_vs_time16->Draw("SAME");
  hist_tracks_vs_time17->SetMarkerColor(7);
  hist_tracks_vs_time17->SetMarkerStyle(5);
  hist_tracks_vs_time17->Draw("SAME");
  hist_tracks_vs_time18->SetMarkerColor(8);
  hist_tracks_vs_time18->SetMarkerStyle(5);
  hist_tracks_vs_time18->Draw("SAME");
 
}
