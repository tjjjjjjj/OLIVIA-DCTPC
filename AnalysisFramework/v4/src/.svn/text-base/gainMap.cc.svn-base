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
   //st1->SetPalette(1);
   TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
   st1->SetNumberContours(NCont);
   st1->SetPadColor(0);
   st1->SetCanvasColor(0);
   st1->SetOptStat(0);
   st1->cd(); 

   TApplication theApp("App",&argc,argv);

   gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
					 "*", "TStreamerInfo",
					 "RIO","TStreamerInfo()");


   TCanvas* c0 = new TCanvas("c0","c0",0,0,1200,600);
   c0->Divide(3,1);
   c0->cd(1);

   TString skimpath = "/net/zwicky/dmtpc/akaboth/projects/DarkMatter/AnalysisFramework/v4/skim/dmtpc_10L_";
   //   TString skimpath = "/net/zwicky/dmtpc/production/skimout_4sh/skim/dmtpc_4sh_";
   TString filename = skimpath;
   filename+=theApp.Argv(1);
   filename+="skim.root";

   cout << filename << endl;
   
   DmtpcSkimDataset d;
   d.openRootFile(filename);
   d.getEvent(0);
   const int ncam = d.event()->ncamera();
   
   TH2F* hSum[ncam];
   TString serial[ncam];
   DmtpcGainMap* gm[ncam];

   for(int i=0; i<ncam; i++)
   {
      TString name="hSum_";
      name+=i;
      hSum[i] = (TH2F*)d.event()->cluster(i)->getImage()->Clone(name);
      hSum[i]->Reset();
      cout << hSum[i] << endl;
      serial[i] = d.event()->cameraSerialNumber(i);
      cout << "Serial: " << serial[i] << endl;
      gm[i] = new DmtpcGainMap(serial[i]);
      
   }


   TObjArray* datasets = new TObjArray(theApp.Argc());
   for(int i=1; i<theApp.Argc(); i++)
   {
      DmtpcSkimDataset* ds = new DmtpcSkimDataset();
      filename = skimpath;
      filename+=theApp.Argv(i);
      filename+="skim.root";
      cout << filename << endl;
      ds->openRootFile(filename);
      datasets->Add(ds);
      
   }

   TObjArray* sumOfImages = calibTools::sumImagesWithN(datasets);

   for(int i=0; i<ncam; i++)
   {
      hSum[i]=(TH2F*)sumOfImages->At(i);
      for(int j=1; j<=hSum[i]->GetNbinsX(); j++)
      {
	 for(int k=1; k<=hSum[i]->GetNbinsY(); k++)
	 {
	    float bin = hSum[i]->GetBinContent(j,k);
	    if(bin>4 || bin <-4)
	    {
	       int nbins=0;
	       double avg=0;
	       for(int m=j-1; m<=j+1; m++)
	       {
		  for(int n=k-1; n<=k+1; n++)
		  {
		     if(m>0 && m<= hSum[i]->GetNbinsX() && 
			n>0 && n<= hSum[i]->GetNbinsY() &&
			!(m==j && n==k))
		     {
			nbins++;
			avg+=hSum[i]->GetBinContent(m,n);
		     }
		  }
	       }
	       avg=avg/double(nbins);
	       hSum[i]->SetBinContent(j,k,avg);
	    }
	 }
      }
   }
   
  
   TCanvas* c1 = new TCanvas("c1","c1",0,0,800,800);
   
   for(int n=0; n<ncam; n++)
   {
      c0->cd(1);
      hSum[n]->DrawCopy("colz");
      c0->Update();

      /////////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////
      ///
      /// Start of spacer finding and fitting; needs the full field of view to be active, and 
      /// needs the spacers to be more-or-less vertical
      ///
      /////////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////      
      int ns;
      cout << "How many spacers are there?" << endl;
      cin >> ns;

      TH1D* hProjection = hSum[n]->ProjectionX("hProjection");
      double maxval = hProjection->GetMaximum();
      for(int i=1; i<=hProjection->GetNbinsX(); i++)
      {
	 hProjection->SetBinContent(i,maxval-hProjection->GetBinContent(i));
      }
      
   
      c0->cd(2);
      hProjection->Draw();
      c0->Update();
      
      TH1F* hProjectionSmooth = (TH1F*)hProjection->Clone("hProjectionSmooth");
      hProjectionSmooth->SetLineColor(kRed);
      hProjectionSmooth->Smooth(500);
      hProjectionSmooth->Draw("SAME");
      c0->Update();

      TH1F* hProjectionSub = (TH1F*)hProjection->Clone("hProjectionSub");
      hProjectionSub->Add(hProjectionSmooth,-1);
      
      c0->cd(3);
      hProjectionSub->Draw();
      c0->Update();
      
      TSpectrum* spacers = new TSpectrum(ns,1.0);
//   spacers->Search(hProjection,1.00,"",0.4);
      spacers->Search(hProjectionSub);

      const int nspacers = spacers->GetNPeaks();
      double spacerpts[nspacers][4];
      double ypos[4];
      double width[nspacers];
      for(int i=0; i<nspacers; i++){width[i]=0;}
      
      for(int i=0; i<4; i++)
      {
	 TH1D* hProjRange = hSum[n]->ProjectionX("hProjRange",i*64+1,(i+1)*64);
	 ypos[i]=(i*64+32)*4;
	 double maxval = hProjRange->GetMaximum();
	 for(int j=1; j<=hProjRange->GetNbinsX(); j++)
	 {
	 hProjRange->SetBinContent(j,maxval-hProjRange->GetBinContent(j));
	 }
	 
	 
	 TH1F* hProjRangeSmooth = (TH1F*)hProjRange->Clone("hProjRangeSmooth");
	 hProjRangeSmooth->SetLineColor(kRed);
	 hProjRangeSmooth->Smooth(500);
	 c0->Update();
      
	 TH1F* hProjRangeSub = (TH1F*)hProjRange->Clone("hProjRangeSub");
//	 hProjRangeSub->Add(hProjRangeSmooth,-1);
	 
	 c0->cd(2);
	 hProjRange->Draw();
	 hProjRangeSmooth->Draw("SAME");
	 c0->Update();
	 
	 c0->cd(3);
	 hProjRangeSub->Draw();
	 c0->Update();
      
	 spacers->Search(hProjRange,1.00);
//      spacers->Search(hProjRangeSub,1.00,"",0.35);
      
	 float* xpos = spacers->GetPositionX();
	 int xposind[nspacers];
	 TMath::Sort(nspacers,xpos,xposind,kFALSE);

	 
	 for(int j=0; j<nspacers; j++)
	 {
	    //spacerpts[j][i]=xpos[xposind[j]];
	    double lowfitrange = (xpos[xposind[j]]-15>0) ? xpos[xposind[j]]-15 : -8;
	    double highfitrange = (xpos[xposind[j]]+15<1024) ? xpos[xposind[j]]+15 : 1032;
	    int lowbin = hProjRangeSub->FindBin(lowfitrange);
	    int highbin = hProjRangeSub->FindBin(highfitrange);
	    double maxval=0; double maxxval=0;
	    for(int k=lowbin; k<=highbin; k++)
	    {
	       if(hProjRangeSub->GetBinContent(k)>maxval)
	       {
		  maxval=hProjRangeSub->GetBinContent(k);
		  maxxval=hProjRangeSub->GetBinCenter(k);
	       }
	    }
	    
	    cout << lowfitrange << "," << highfitrange << ":" << maxval << "," << maxxval << endl;
	    hProjRangeSub->Fit("gaus","Q","SAME",lowfitrange,highfitrange);
	    if (hProjRangeSub->GetFunction("gaus")->GetParameter(2)<15)
	       width[j]+=hProjRangeSub->GetFunction("gaus")->GetParameter(2);
	    else
	       width[j]+=15;
	    spacerpts[j][i]=hProjRangeSub->GetFunction("gaus")->GetParameter(1);
	    if(spacerpts[j][i]>1024) spacerpts[j][i]=1024;
	    c0->Update();
	    
	    cout << xpos[xposind[j]] << ":" << spacerpts[j][i] << ", " << ypos[i] << endl;
	 }

      }
      for(int i=0; i<nspacers; i++) {width[i]=width[i]/4.0;}
      
      hProjection->Draw();
      c0->Update();

      c1->cd();
      double center[nspacers],diffs[nspacers-1];
      double slope[nspacers],intercept[nspacers];
      
      TGraph* g = new TGraph(4,ypos,spacerpts[0]);
      g->SetMarkerStyle(21);
      g->GetXaxis()->Set(1024,0,1024);
      g->SetMinimum(0);
      g->SetMaximum(1024);
      g->Draw("AP");
      g->Fit("pol1","");
      TF1* pol1 = g->GetFunction("pol1");
      center[0]=pol1->Eval(512);
      intercept[0]=-1*pol1->GetParameter(0)/pol1->GetParameter(1); 
      slope[0] = 1/pol1->GetParameter(1);
      cout << center[0] << endl;
      cout << intercept[0] << "," << slope[0]
	   << "," << width[0] << endl;

      gm[n]->addSpacer(slope[0],intercept[0],width[0]);
      c1->Update();

      for(int i=1; i<nspacers; i++)
      {
	 TGraph* gg = new TGraph(4,ypos,spacerpts[i]);
	 gg->SetMarkerStyle(21);
	 gg->Draw("P SAME");
	 gg->Fit("pol1","","SAME");
	 pol1 = gg->GetFunction("pol1");
	 center[i]=pol1->Eval(512);
	 intercept[i]=-1*pol1->GetParameter(0)/pol1->GetParameter(1); 
	 slope[i] = 1/pol1->GetParameter(1);
	 cout << center[i] << endl;
	 cout << intercept[i] << "," << slope[i]
	      << "," << width[i] << endl;
	 gm[n]->addSpacer(slope[i],intercept[i],width[i]);
	 c1->Update();
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////
      ///
      /// End of spacer finding and fitting
      ///
      /////////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////      

      double avgdiff=0;
      for(int i=0; i<nspacers-1; i++)
      {
	 diffs[i]=center[i+1]-center[i];
	 cout << diffs[i] << endl;
	 avgdiff+=diffs[i];
      }
      avgdiff=avgdiff/double(nspacers-1);
      
      cout << "The average spacer difference is: " << avgdiff << endl;
      
      TH2F* spacerMask = calibTools::makeSpacerMask(hSum[n],nspacers,slope,intercept,width);
      TH2F* hBlur = calibTools::blurWithoutSpacers(hSum[n],spacerMask,0.3);
      TH2F* gmregion = (TH2F*)hBlur->Clone("gmregion");
      gmregion->Multiply(spacerMask);
      double minval = gmregion->GetMinimum();
      for(int m=1; m<=hBlur->GetNbinsX(); m++)
      {
 	 for(int n=1; n<=hBlur->GetNbinsY(); n++)
 	 {
 	    hBlur->SetBinContent(m,n,hBlur->GetBinContent(m,n)-minval);
	    if(hBlur->GetBinContent(m,n)<0)
	       hBlur->SetBinContent(m,n,0);
	 }
      }
      c0->cd(1);
      hBlur->Draw("colz");
      gm[n]->setGainMap(hBlur);
      c0->Update();
   }
   
   TString s;
   cout << "Save to file? (y/n)" << endl;;
   cin >> s;
   if(s=="y")
   {
      TString outfilename;
      cout << "Filename?" << endl;
      cin >> outfilename;
      TFile* outfile = new TFile(outfilename,"RECREATE");
      for(int n=0; n<ncam; n++)
      {
	 gm[n]->Write();
      }
   }



   return 0;
}
