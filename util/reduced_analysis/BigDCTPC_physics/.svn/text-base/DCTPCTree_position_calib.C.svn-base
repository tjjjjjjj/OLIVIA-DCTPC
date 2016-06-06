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

//This macro is used to calibrate the position of the field cage with respect to the camera

void DCTPCTree::Loop()
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
  
TFile *outtree = new TFile("$BigDCTPC_physics_input");
TTree *dctreepc = (TTree*)outtree->Get("dctpc_eventinfo");
DCTPCTree aStep(dctreepc);
  
TH2D *hist_pos=new TH2D("","",200,0,1024,200,0,1024);

double rstart=0.;
double rend=0.;
int points=0;
int points2=0;
double x[1000000],y[1000000],x2[1000000],y2[1000000];

//create event loop
  Long64_t nentries = dctreepc->GetEntries();
  cout << "Number of Entries: " << nentries << endl;
 
for (int event = 0; event<nentries; event++)
{
      aStep.GetEntry(event); 
      
      if(event%10000==0)
	  cout<<((double)event/(double)nentries)*100.<<"%"<<endl;	
 


//constants are added/subtracted here (before the fit even takes place) because I want the cuts to isolate the rings. That is, these are used just to make sure that the cuts are most effective. 
rstart=sqrt(pow(aStep.Track_x_start_pix-3.,2)+pow(aStep.Track_y_start_pix-19.,2));
rend=sqrt(pow(aStep.Track_x_end_pix-3.,2)+pow(aStep.Track_y_start_pix-19.,2));
 
 //these cuts are meant to isolate clean alphas coming from the rings
     if (
     aStep.Edge==0  && 
     aStep.Ntrig <=2 && 
     rstart>505 && 
     rstart<660 &&
     rend<500 
     )
	{	

	 hist_pos->Fill(512+aStep.Track_x_start_pix,512+aStep.Track_y_start_pix);
	 x[points]=512+(aStep.Track_x_start_pix);
	 y[points]=512+(aStep.Track_y_start_pix);
	 points++;
	}
	
	 if ( 
     aStep.Ntrig <=2 
     )
	{	
	 x2[points2]=512+(aStep.Track_x_start_pix);
	 y2[points2]=512+(aStep.Track_y_start_pix);
	 points2++;
	}
	   
}



TH2D *dummy =new TH2D("","",100,-200,1250,100,-200,1250);
dummy->Draw();
dummy->GetYaxis()->SetTitleOffset(1.2);
dummy->SetXTitle("X coordinate (pixels)");
dummy->SetYTitle("Y coordinate (pixels)");
gr = new TGraph(points,x,y);
gr2 = new TGraph(points2,x2,y2);


gr->Draw("p SAME");

TVirtualFitter::SetDefaultFitter("Minuit");  //default is Minuit
TVirtualFitter *fitter = TVirtualFitter::Fitter(0, 3);
fitter->SetFCN(myfcn);

fitter->SetParameter(0, "x0",   0, 0.1, 0,0);
fitter->SetParameter(1, "y0",   0, 0.1, 0,0);
fitter->SetParameter(2, "R",    1, 0.1, 0,0);

Double_t arglist[1] = {0};
fitter->ExecuteCommand("MIGRAD", arglist, 0);

//Draw the circle on top of the points
TArc *arc = new TArc(fitter->GetParameter(0),
fitter->GetParameter(1),fitter->GetParameter(2));
arc->SetLineColor(kRed);
arc->SetLineWidth(1);
arc->Draw("SAME");
gr->Draw("p SAME");

new TCanvas;
dummy->Draw();
gr2->Draw("p SAME");
arc->Draw("SAME"); 
gr2->Draw("p SAME");
}


void myfcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
   //minimisation function computing the sum of squares of residuals
   Int_t np = gr->GetN();
   f = 0;
   Double_t *x = gr->GetX();
   Double_t *y = gr->GetY();
   for (Int_t i=0;i<np;i++) {
      Double_t u = x[i] - par[0];
      Double_t v = y[i] - par[1];
      Double_t dr = par[2] - TMath::Sqrt(u*u+v*v);
      f += dr*dr;
   }
}
