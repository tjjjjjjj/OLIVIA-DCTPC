#include "TROOT.h"
#include "TApplication.h"
#include "TPluginManager.h"
#include "TColor.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TSpectrum.h"
#include "TGraph.h"
#include "TMath.h"
#include "TF1.h"
#include "calibTools.hh"
#include "../../../MaxCam/DmtpcSkimDataset.hh"
#include "../../../MaxCam/MaxCamClusterImage.hh"
#include "../../../MaxCam/MaxCamImageTools.hh"
#include "../../../MaxCam/MaxCamConfig.hh"
#include "../../../MaxCam/DmtpcGainMap.hh"

#include <iostream>

using namespace std;

/* Main function */ 

int main(int argc, char *argv[])
{

   const Int_t NRGBs = 7;
   const Int_t NCont = 104;
   
   Double_t stops[NRGBs] = { 0.00, 0.10, 0.30, 0.50, 0.65, 0.90, 1.00 };
   Double_t red[NRGBs]   = { 0.31, 0.00, 0.00, 0.10, 0.87, 1.00, 0.51 };
   Double_t green[NRGBs] = { 0.00, 0.00, 0.81, 0.90, 1.00, 0.20, 0.00 };
   Double_t blue[NRGBs]  = { 0.48, 0.51, 1.00, 0.10, 0.12, 0.00, 0.00 };
      
   TStyle *st1 = new TStyle("st1","my style");
   st1->SetPalette(1);
   //st1->CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
   //st1->SetNumberContours(NCont);
   st1->SetPadColor(0);
   st1->SetCanvasColor(0);
   st1->SetOptStat(0);
   st1->cd(); 

   TApplication theApp("App",&argc,argv);

   gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
					 "*", "TStreamerInfo",
					 "RIO","TStreamerInfo()");


   TCanvas* c0 = new TCanvas("c0","c0",0,0,1200,600);
   c0->Divide(2,1);
   c0->cd(1);

   DmtpcSkimDataset d;
   d.openRootFile(theApp.Argv(1));
   d.getEvent(0);

   TH2F* hSum = (TH2F*)d.event()->cluster(0)->getImage()->Clone("hSum");
   hSum->Reset();

   cout << "tag a" << endl;

   TString serial = "A80333";
   
   cout << "tag b" << endl;
   
   DmtpcGainMap* gm = new DmtpcGainMap(serial);
   
   cout << "tag c" << endl;
   
   for(int i=1; i<theApp.Argc(); i++)
   {
      DmtpcSkimDataset ds;
      ds.openRootFile(theApp.Argv(i));
      ds.getEvent(0);
      const int ncam = ds.event()->ncamera();
      
      TObjArray* sumOfImages = calibTools::sumImages(ds,"gamma");
      c0->cd(1);
      sumOfImages->At(0)->Draw("colz");

      if(i==0)
      {
	 hSum = (TH2F*)sumOfImages->At(0)->Clone("hSum");
      }
      else
      {
	 hSum->Add((TH2F*)sumOfImages->At(0));
      }

      cout << "Finished processing file number" << i << " of " << theApp.Argc()-1 << endl;

   }
   
   for(int i=1; i<=hSum->GetNbinsX(); i++)
   {
      for(int j=1; j<=hSum->GetNbinsY(); j++)
      {
	 hSum->SetBinContent(i,j,
			     hSum->GetBinContent(i,j)/double(theApp.Argc()));
      }
   }

      c0->cd(1);
//   hSum->SetMaximum(5*(theApp.Argc()-1));
//   hSum->SetMinimum(-5*(theApp.Argc()-1));
   hSum->Draw("colz");
   c0->Update();

   getchar();
      
   TH1D* hProjection = hSum->ProjectionX("hProjection");
   double maxval = hProjection->GetMaximum();
   for(int i=1; i<=hProjection->GetNbinsX(); i++)
   {
      hProjection->SetBinContent(i,maxval-hProjection->GetBinContent(i));
   }
   
   
   c0->cd(2);
   hProjection->Draw();
   c0->Update();
   
   TSpectrum* spacers = new TSpectrum(10,1.0);
   spacers->Search(hProjection,1.00,"",0.4);

   const int nspacers = spacers->GetNPeaks();
   double spacerpts[nspacers][8];
   double ypos[8];
   double width[nspacers];
   for(int i=0; i<nspacers; i++){width[i]=0;}
  
   for(int i=0; i<8; i++)
   {
      TH1D* hProjRange = hSum->ProjectionX("hProjRange",i*32+1,(i+1)*32);
      ypos[i]=(i*32+16)*4;
      double maxval = hProjRange->GetMaximum();
      for(int j=1; j<=hProjRange->GetNbinsX(); j++)
      {
	 hProjRange->SetBinContent(j,maxval-hProjRange->GetBinContent(j));
      }
      
      
      c0->cd(2);
      hProjRange->Draw();
      c0->Update();
      
      spacers->Search(hProjRange,1.00,"",0.4);
      
      float* xpos = spacers->GetPositionX();
      int xposind[nspacers];
      TMath::Sort(nspacers,xpos,xposind,kFALSE);
      for(int j=0; j<nspacers; j++)
      {
	 //spacerpts[j][i]=xpos[xposind[j]];
	 double lowfitrange = (xpos[xposind[j]]-15>0) ? xpos[xposind[j]]-15 : -8;
	 double highfitrange = (xpos[xposind[j]]+15<1024) ? xpos[xposind[j]]+15 : 1032;

	 // cout << lowfitrange << "," << highfitrange << endl;
	 hProjRange->Fit("gaus","Q","SAME",lowfitrange,highfitrange);
	 if (hProjRange->GetFunction("gaus")->GetParameter(2)<15)
	    width[j]+=hProjRange->GetFunction("gaus")->GetParameter(2);
	 else
	    width[j]+=15;
	 spacerpts[j][i]=hProjRange->GetFunction("gaus")->GetParameter(1);
	 if(spacerpts[j][i]>1024) spacerpts[j][i]=1024;
	 c0->Update();
		 
      }

   }
   for(int i=0; i<nspacers; i++) {width[i]=width[i]/8.0;}

   hProjection->Draw();
   c0->Update();

   TCanvas* c1 = new TCanvas("c1","c1",0,0,800,800);
   c1->cd();

   double center[nspacers],diffs[nspacers-1];
   double slope[nspacers],intercept[nspacers];

   TGraph* g = new TGraph(8,spacerpts[0],ypos);
   g->SetMarkerStyle(21);
   g->GetXaxis()->Set(1024,0,1024);
   g->Draw("AP");
   g->Fit("pol1","Q");
   TF1* pol1 = g->GetFunction("pol1");
   center[0]=pol1->GetX(512);
   intercept[0]=pol1->GetParameter(0); slope[0] = pol1->GetParameter(1);
   cout << center[0] << endl;
   cout << pol1->GetParameter(0) << "," << pol1->GetParameter(1) 
	<< "," << width[0] << endl;
   gm->addSpacer(pol1->GetParameter(1),pol1->GetParameter(0),width[0]);
   c1->Update();

   for(int i=1; i<nspacers; i++)
   {
      TGraph* gg = new TGraph(8,spacerpts[i],ypos);
      gg->SetMarkerStyle(21);
      gg->Draw("P SAME");
      gg->Fit("pol1","Q","SAME");
      pol1 = gg->GetFunction("pol1");
      center[i]=pol1->GetX(512);
      intercept[i]=pol1->GetParameter(0); slope[i] = pol1->GetParameter(1);
      cout << center[i] << endl;
      cout << pol1->GetParameter(0) << "," << pol1->GetParameter(1)
	   << "," << width[i] << endl;
      gm->addSpacer(pol1->GetParameter(1),pol1->GetParameter(0),width[i]);
      c1->Update();
   }

   double avgdiff=0;
   for(int i=0; i<nspacers-1; i++)
   {
      diffs[i]=center[i+1]-center[i];
      cout << diffs[i] << endl;
      avgdiff+=diffs[i];
   }
   avgdiff=avgdiff/double(nspacers-1);
   
   cout << "The average spacer difference is: " << avgdiff << endl;

   TH2F* spacerMask = calibTools::makeSpacerMask(hSum,nspacers,slope,intercept,width);


   //hSum = MaxCamImageTools::blur(hSum,1,0.3);

   TH2F* hBlur = calibTools::blurWithoutSpacers(hSum,spacerMask,0.3);
   c0->cd(1);
   hBlur->Draw("colz");
   gm->setGainMap(hBlur);
   c0->Update();


   TString s;
   cout << "Save to file? (y/n)" << endl;;
   cin >> s;
   if(s=="y")
   {
      TString outfilename;
      cout << "Filename?" << endl;
      cin >> outfilename;
      TFile* outfile = new TFile(outfilename,"RECREATE");
      gm->Write();
      //Blur->SetName("gainMap");
      //hBlur->Write();
      //hProjection->Write();
   }



   return 0;
}
