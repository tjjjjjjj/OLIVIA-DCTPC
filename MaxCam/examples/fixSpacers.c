#include "../../MaxCam/MaxCamMC.hh"
#include "../../MaxCam/MaxCamImageTools.hh"
#include "../../MaxCam/MaxCamClusterImage.hh"

#include "TF1.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TVector3.h"
#include "TNtuple.h"

#include "TFile.h"
#include "TSystem.h"
#include "TTimeStamp.h"

#include "TRandom.h"
#include "TRandom3.h"
#include "TMath.h"

#include "TTree.h"
#include "TCanvas.h"
#include "TPolyMarker.h"
#include "TClonesArray.h"

#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>
//typedef vector<int> vector_int;

#include "../../MaxCam/DmtpcEvent.hh"

using namespace std;

using std::vector;

TCanvas c = TCanvas("c");




vector<double> fixImages(TString date, int runNum, double mindist)
{

  double gain=30; // counts/keV

  MaxCamClusterImage* clust;
  TTree* Corrected = new TTree("Corrected","Corrected");
  TTree* Simulation;
  TTree* dmtpc;

  TH2F* savedImg = 0;
  Corrected->Branch("fixedImage",&savedImg);

  TString fname = "SimData/dmtpc_run" + date + "_" + runNum + "_sim.root"; // open skim file
  TFile *f =  new TFile(fname,"READ");
  dmtpc = (TTree*) f->Get("dmtpc");
  Simulation = (TTree*) f->Get("Simulation");

  DmtpcEvent* d;
  dmtpc->SetBranchAddress("event",&d);
  float simD;
  //Simulation->SetBranchAddress("simdata",simD); // for old version of runSimulations
  Simulation->SetBranchAddress("E",&simD);
  ////Simulation->GetEntry(entryNum);

  TH2F* image;
  TH2F* tempimg;
  TH2F* baseimage;
  TH2F* copy;
  TH2F* clusterImage;
  
  // Get spacer info
  fname = "SimData/dmtpc_run" + date + "_" + runNum + "_spacers.root"; // open spacer file
  TFile *file1 = new TFile(fname,"UPDATE");
  TTree* spacer = (TTree*) file1->Get("spacer");
  double spacerRef;
  spacer->SetBranchAddress("spacerY",&spacerRef);
  spacer->GetEntry(0);
  
  int nbinsx;
  int nbinsy;
  
  // PARAMETERS TO SET
  //double dfitRange=20.0; // distance from spacer to boundary of fit region, in pixels (total range is 2x+1)
  const int fitRangeB=12; // in bins
  int sumRange = int(fitRangeB/2.3);
  int bin2pix=4; // number of pixels per bin
  
  TH1D* profile; 
  TH1D* profile1;
  TH1D* midpt;
  TH1D* pNew;

  int sN = int(spacer->GetEntries());
  const int spacerNum = sN;
  double spacer_array[spacerNum][fitRangeB*2+1];
  
  // Find bin corrections from profile
  vector<double> spacerY;
  double spacerMid;
  double binVal,expVal; // actual bin value and expected value
  double val_left, val_right, sum_left, sum_right;
  int spacerBin;

  // Initializations for variables in big loop
  double offset;
  TTimeStamp* ttime = new TTimeStamp();
  int ntracks;

  int spacerOffset;
  double xsum,counter,xmid;

  int spacerOffset;
  double tWidth = 2.3; //2.3 <- best value for parameter
  double xIntg = 1.5; //1.5
  double xfactor;

  //for (int eN=0; eN<dmtpc->GetEntries(); eN++)
  for (int eN=0; eN<10; eN++)
    {

      dmtpc->GetEntry(eN);
      Simulation->GetEntry(eN);

      image = (TH2F*)d->ccdData(0);

      tempimg = (TH2F*)image->Clone("tempimg");
      baseimage = (TH2F*)tempimg->Clone("baseimage");
      copy = (TH2F*)tempimg->Clone("copy");
      clusterImage = (TH2F*)tempimg->Clone("clusterImage");
      clusterImage->Scale(0.0);
      
      nbinsx = copy->GetNbinsX();
      nbinsy = copy->GetNbinsY();

      profile = image->ProjectionY("profile"); // for now project whole image; later can project just clusters
      profile1 = image->ProjectionY("profile1");
      midpt = image->ProjectionY("midpt");
      midpt->Scale(0.0);

      profile->Smooth(7);
  
      double holInt = profile->Integral();

      //Find bounds
      double prof_min = profile1->GetMinimum()-1.0;
      double min = profile1->Integral()/nbinsx/nbinsy;
      double max = image->GetMaximum();
      max = .99*max+.01*min;

      for (int i=0; i<spacerNum; i++)
	{
	  spacer->GetEntry(i);
	  spacerMid = spacerRef;
	  spacerBin = int(spacerMid/bin2pix);

	  sum_left=0.0;
	  sum_right=0.0;
	  for(int ii=0; ii<sumRange; ii++)
	    {
	      sum_left = sum_left + profile->GetBinContent(spacerBin-fitRangeB+ii);
	      sum_right = sum_right + profile->GetBinContent(spacerBin+fitRangeB-ii);
	    }
	  val_left = sum_left/double(sumRange);
	  val_right = sum_right/double(sumRange);
 
	  // For spacers too close to edge:
	  if ((nbinsy-spacerBin)<(fitRangeB*.75))    val_right = val_left;
	  elseif ((spacerBin)<(fitRangeB*.75))    val_left = val_right;
	  double localMax = (val_left + val_right)/2.0;
	  double slope = (val_right - val_left)/2.0/fitRangeB;
      
	  // Correct for bin values from spacerBin-fitRangeB to spacerBin+fitRangeB
	  for (int j=0; j<=2*fitRangeB; j++)
	    { 
	      binVal = profile1->GetBinContent(spacerBin+j-fitRangeB);
	      expVal = localMax+slope*(j-fitRangeB);
	      spacer_array[i][j]=expVal-binVal;
	    }

	  spacerY.push_back(spacerMid);      
      	}

      baseimage = MaxCamImageTools::blur(baseimage,1,0.5);
      clust = new MaxCamClusterImage(baseimage,ttime);
      ntracks = MaxCamImageTools::findClustersNG(baseimage,clust,min+10,max,10,mindist);
  

      for(int v =0; v<ntracks; v++)
	{
	  vector<int> bins = clust->getCluster(v);
	  for(int i=0; i<int(bins.size()); i++)
	    {
	      tempimg->SetBinContent(bins[i],1.0);
	    }
	}

      for(int j=1; j<=tempimg->GetNbinsY(); j++) 
	{
	  xsum=0;
	  counter=0;
	  for(int i=1; i<=tempimg->GetNbinsX(); i++)
	    {
	      if(tempimg->GetBinContent(i,j)==1.0)
		{
		  xsum = xsum + i;
		  counter = counter +1;
		}
	    }
	  if (counter==0) xmid=0.0;
	  else xmid=xsum/counter;
	  midpt->SetBinContent(j,xmid);
	}

      for(int v =0; v<ntracks; v++)
	{
	  vector<int> bins = clust->getCluster(v);
	  for(int i=0; i<int(bins.size()); i++)
	    {
	      int ii = double((bins[i]-nbinsx)%(nbinsx+2))-1;
	      int j= int(ceil((bins[i]-nbinsx)/(nbinsx+2)))+1;
	      binVal = image->GetBinContent(bins[i]);
	      xmid = midpt->GetBinContent(j);
	      xfactor = exp(-pow((fabs(ii-xmid)*2.4/tWidth),2.3))/xIntg*2.4/tWidth;
	      if(1)
		{
		  for(int k=0; k<spacerNum; k++)
		    {	 
		      spacerOffset = j-int(spacerY[k]/double(bin2pix));
		      if(fabs(spacerOffset)<=fitRangeB)
			{
			  expVal = binVal+xfactor*spacer_array[k][spacerOffset+fitRangeB];
			  if(expVal<binVal) expVal=binVal;
			  image->SetBinContent(bins[i],expVal);
			  clusterImage->SetBinContent(bins[i],1);
			}
		    }
		}
	    }
	}

      pNew = image->ProjectionY("pNew");
      double fixInt = pNew->Integral();

      //cout << "Holey Intg.: " << holInt/gain << endl;
      //cout << "Fixed Int: " << fixInt/gain << endl;
      //cout << "Orig Int: " << simD << endl;

      if (eN%100==0) cout <<eN<< " ";

      savedImg = image;
      Corrected->Fill();

    }

  f->Close();

  Corrected->Write();
  file1->Close();

}




vector<double> fixImage(TString date, int runNum, int entryNum, double mindist)
     // Returns visible energy, corrected visible energy, and energy input for MC
{

  double gain=30; // counts/keV

  MaxCamClusterImage* clust;
  TTree* Simulation;

  TString fname = "SimData/dmtpc_run" + date + "_" + runNum + "_sim.root"; // open skim file
  TFile *f =  new TFile(fname,"READ");TTree* 
  TTree* dmtpc = (TTree*) f->Get("dmptc");
  Simulation = (TTree*) f->Get("Simulation");
  DmtpcEvent* d;
  dmtpc->SetBranchAddress("event",&d);
  dmtpc->GetEntry(entryNum);
  float simD;
  //Simulation->SetBranchAddress("simdata",simD); // for old version of runSimulations
  Simulation->SetBranchAddress("E",&simD); //Integral //Escint // Choose which form of energy to read out
  Simulation->GetEntry(entryNum);

  TH2F* image = (TH2F*)d->ccdData(0);

  TH2F* tempimg = (TH2F*)image->Clone("tempimg");
  TH2F* baseimage = (TH2F*)tempimg->Clone("baseimage");
  TH2F* copy = (TH2F*)tempimg->Clone("copy");
  TH2F* clusterImage = (TH2F*)tempimg->Clone("clusterImage");
  clusterImage->Scale(0.0);
  
  // Get spacer info
  fname = "SimData/dmtpc_run" + date + "_" + runNum + "_spacers.root"; // open spacer file
  TFile *file1 = new TFile(fname,"READ");
  TTree* spacer = (TTree*) file1->Get("spacer");
  double spacerRef;
  spacer->SetBranchAddress("spacerY",&spacerRef);
  spacer->GetEntry(0);
  
  int nbinsx = copy->GetNbinsX();
  int nbinsy = copy->GetNbinsY();
  
  // PARAMETERS TO SET
  //double fitRange=20.0; // distance from spacer to boundary of fit region, in pixels (total range is 2x+1)
  const int fitRangeB=12; // in bins
  int sumRange = int(fitRangeB/2.3);
  int bin2pix=4; // number of pixels per bin
  
  TH1D* profile = new TH1D("profile","",image->GetNbinsY(),image->GetYaxis()->GetXmin(),image->GetYaxis()->GetXmax());
  TH1D* profile1 = new TH1D("profile1","",image->GetNbinsY(),image->GetYaxis()->GetXmin(),image->GetYaxis()->GetXmax());
  profile = image->ProjectionY("profile"); // for now project whole image; later can project just clusters
  profile1 = image->ProjectionY("profile1");
  TH1D* midpt = image->ProjectionY("midpt");
  midpt->Scale(0.0);

  double holInt = profile->Integral();

  // GET SPACERS
  profile->Smooth(7);
  
  int sN = int(spacer->GetEntries());
  const int spacerNum = sN;
  double spacer_array[spacerNum][fitRangeB*2+1];
  
  // Find bounds
  double prof_min = profile1->GetMinimum()-1.0;
  double min = profile1->Integral()/nbinsx/nbinsy;
  double max = image->GetMaximum();
  max = .99*max+.01*min;
  //cout << "Pixel bounds: " << min << " " << max << endl;

  // Find bin corrections from profile
  vector<double> spacerY;
  double spacerMid;
  double binVal,expVal; // actual bin value and expected value
  double val_left,val_right, sum_left, sum_right;
  int spacerBin;

  double offset;

  cout << "spacers: ";
  for (int i=0; i<spacerNum; i++)
    {
      spacer->GetEntry(i);
      spacerMid = spacerRef;
      spacerBin = int(spacerMid/bin2pix);

      sum_left=0.0;
      sum_right=0.0;
      for(int ii=0; ii<sumRange; ii++)
	{
	  sum_left = sum_left + profile->GetBinContent(spacerBin-fitRangeB+ii);
	  sum_right = sum_right + profile->GetBinContent(spacerBin+fitRangeB-ii);
	}
      val_left = sum_left/double(sumRange);
      val_right = sum_right/double(sumRange);
 
      // For spacers too close to edge:
      if ((nbinsy-spacerBin)<(fitRangeB*.75)) val_right = val_left;
      elseif ((spacerBin)<(fitRangeB*.75)) val_left = val_right;
      double localMax = (val_left + val_right)/2.0;
      double slope = (val_right - val_left)/2.0/fitRangeB;
      
      // Correct for bin values from spacerBin-fitRangeB to spacerBin+fitRangeB
      for (int j=0; j<=2*fitRangeB; j++)
	{ 
	  binVal = profile1->GetBinContent(spacerBin+j-fitRangeB);
	  expVal = localMax+slope*(j-fitRangeB);
	  spacer_array[i][j]=expVal-binVal;
	}

      spacerY.push_back(spacerMid);      
      
      cout << spacerMid << " ";
    }

  cout << endl;


  baseimage = MaxCamImageTools::blur(baseimage,1,0.5);
  TTimeStamp* ttime = new TTimeStamp();
  clust = new MaxCamClusterImage(baseimage,ttime);
  int ntracks = MaxCamImageTools::findClustersNG(baseimage,clust,min+10,max,10,mindist);

  int spacerOffset;
  for(int v =0; v<ntracks; v++)
    {
      vector<int> bins = clust->getCluster(v);
      for(int i=0; i<int(bins.size()); i++)
	{
		      tempimg->SetBinContent(bins[i],1.0);
	}
    }

  double xsum,counter,xmid;
  for(int j=1; j<=tempimg->GetNbinsY(); j++) 
    {
      xsum=0;
      counter=0;
      for(int i=1; i<=tempimg->GetNbinsX(); i++)
	{
	  if(tempimg->GetBinContent(i,j)==1.0)
	    {
	      xsum = xsum + i;
	      counter = counter +1;
	    }
	}
      if (counter==0) xmid=0.0;
      else xmid=xsum/counter;
      midpt->SetBinContent(j,xmid);
    }

  int spacerOffset;
  double tWidth = 2.3;
  double xIntg = 1.5;
  double xfactor;
  for(int v =0; v<ntracks; v++)
    {
      vector<int> bins = clust->getCluster(v);
      for(int i=0; i<int(bins.size()); i++)
	{
	  int ii = double((bins[i]-nbinsx)%(nbinsx+2))-1;
	  int j= int(ceil((bins[i]-nbinsx)/(nbinsx+2)))+1;
	  binVal = image->GetBinContent(bins[i]);
	  xmid = midpt->GetBinContent(j);
	  xfactor = exp(-pow((fabs(ii-xmid)*2.4/tWidth),2.3))/xIntg*2.4/tWidth;
	  if(1) // binVal>min && binVal<max <- no longer need min requirement
	    {
	      for(int k=0; k<spacerNum; k++)
		{	 
		  spacerOffset = j-int(spacerY[k]/double(bin2pix));
		  if(fabs(spacerOffset)<=fitRangeB)
		    {
		      expVal = binVal+xfactor*spacer_array[k][spacerOffset+fitRangeB];
		      //if(expVal>max+25) expVal=max;
		      if(expVal<binVal) expVal=binVal;
		      image->SetBinContent(bins[i],expVal);
		      clusterImage->SetBinContent(bins[i],1);
		    }
		}
	    }
	}
    }


  TH1D* pNew = image->ProjectionY("pNew");
  double fixInt = pNew->Integral();
  cout << "Holey Intg.: " << holInt/gain << endl;
  cout << "Fixed Intg.: " << fixInt/gain << endl;
  cout << "Orig. Intg.: " << simD << endl;

  
  c.Clear();
  c.Divide(3,2);

  c.cd(1);
  //  Original Image
  copy->Draw("colz");

  c.cd(2);
  //  Corrected Image
  image->Draw("colz");

  c.cd(3);
  //  Corrected pixels
  clusterImage->Draw("colz");

  c.cd(4);
  //  Original profile
  profile1->Draw();

  c.cd(5);
  //  Corrected profile
  TH1D* pNew = image->ProjectionY("pNew");
  pNew->Draw();

  c.cd(6);
  //  Smoothed image (used by findClustersNG to find clusters)
  baseimage->Draw("colz");

  vector<double> E;
  E.push_back(fixInt/gain);
  E.push_back(simD);
  E.push_back(holInt/gain);

//     f->Close();
//     file1->Close();

  return E;

}




void Ecompare(TString date, int runNum)
     // Use to compare energy before and after spacer correction
{

  TFile* file1;
  DmtpcEvent* dm1;

  TString fname = "SimData/dmtpc_run" + date + "_" + runNum + "_sim.root"; // open skim file
  file1 =  new TFile(fname,"READ");

  TH2F* Ehist = new TH2F("Ehist","",1000,0,25000,700,-5000,30000);
  TH2F* holHist = new TH2F("holhist","",1000,0,25000,700,-5000,30000);
  TH2F* Hist = new TH2F("Hist","",1000,0,25000,700,-5000,30000);
  TProfile* Epro = new TProfile("Corrected","",20,0,2000,-5000,30000);
  TProfile* Hpro = new TProfile("Spacers","",20,0,2000,-5000,30000);
  TH1F* s1 = new TH1F("Corrected 0.5 MeV","",13,-500,6000);
  TH1F* s2 = new TH1F("Corrected 1.0 MeV","",13,-500,6000);
  TH1F* s3 = new TH1F("Corrected 5.0 MeV","",13,-500,6000);
  TH1F* s1h = new TH1F("Spacers 0.5 MeV","",13,-500,6000);
  TH1F* s2h = new TH1F("Spacers 1.0 MeV","",13,-500,6000);
  TH1F* s3h = new TH1F("Spacers 5.0 MeV","",13,-500,6000);
  
  vector<double> E;
  for(int i=0; i<30; i++) // Sum skimmed images
    {
      E = fixImage(date,runNum,i,10);
      if(E[1]>20.0)
	{
//       Ehist->Fill(E[1],E[0]);
//       holHist->Fill(E[1],E[2]);
      Hist->Fill(18000,E[0]);
      Hist->Fill(18000,E[2]);
//       Epro->Fill(E[1],E[0]);
//       Hpro->Fill(E[1],E[2]); 

      Ehist->Fill(E[1],(E[0]-E[1])/E[1]);
      holHist->Fill(E[1],(E[2]-E[1])/E[1]);
      Epro->Fill(E[1],(E[0]-E[1])/E[1]);
      Hpro->Fill(E[1],(E[2]-E[1])/E[1]); 

      if(E[0]>400 && E[0]<500) s1->Fill(E[0],1);
      elseif(E[0]>900 && E[0]<1000) s2->Fill(E[0],1);
      elseif(E[0]>5000 && E[0]<5200) s3->Fill(E[0],1);

      if(E[2]>400 && E[0]<500) s1h->Fill(E[2],1);
      elseif(E[2]>900 && E[2]<1000) s2h->Fill(E[2],1);
      elseif(E[2]>5000 && E[2]<5200) s3h->Fill(E[2],1);
    }
    }

  c.Clear();

  if (1)
    {
      c.Divide(2,1);
      c.cd(1);
      Epro->GetXaxis()->SetTitle("E_{MC} (keV)");
      Epro->GetYaxis()->SetTitle("(E_{C}-E_{MC})/E_{MC}");
      Epro->Draw();
      c.cd(2);
      Hpro->GetXaxis()->SetTitle("E_{MC} (keV)");
      Hpro->GetYaxis()->SetTitle("(E_{S}-E_{MC})/E_{MC}");
      Hpro->Draw();
      //c.cd(3);
      //Hist->Draw();
    }

  else
    {
      c.Divide(3,2);
      c.cd(1);
      s1->Draw();
      c.cd(2);
      s3->Draw();
      c.cd(3);
      s3->Draw();
      c.cd(4);
      s1h->Draw();
      c.cd(5);
      s2h->Draw();
      c.cd(6);
      s3h->Draw();
    }
}




void sumImages(TH2F* ImgSum, TString date, int runNum)
{
  TFile* file1;
  DmtpcEvent* dm1;
  TH2F* tempImg;
  int ntracks;

  TString fname = "SimData/dmtpc_run" + date + "_" + runNum + "_sim.root"; // open skim file
  file1 =  new TFile(fname,"READ");
  TTree* dmtpc = (TTree*) file1->Get("dmtpc");
  dmtpc->SetBranchAddress("event",&dm1);
  dmtpc->GetEntry(0);
 
     cout<< "Entries: " << dmtpc->GetEntries()<<endl;
     for(int i=1; i<dmtpc->GetEntries(); i++) // Sum skimmed images
       {
         dmtpc->GetEntry(i);
	 tempImg = dm1->ccdData(0);
	 ImgSum->Add(ImgSum,tempImg,1,1);
       }

}




void fixSpacers(TString date, int runNum, int ntimes=7, double width=5.6, double amp=.25)
     // Use to find and save spacer locations; last three inputs are optimized peak fitting parameters
{  
  
  TH2F* ImgSum;
  DmtpcEvent* dm1;
  TH2F* h,tempImg;

  TString fname = "SimData/dmtpc_run" + date + "_" + runNum + "_sim.root"; // open skim file
  TFile* file1 = new TFile(fname,"READ");
  TTree* dmtpc = (TTree*) file1->Get("dmtpc");
  dmtpc->SetBranchAddress("event",&dm1);
  dmtpc->GetEntry(0);
 
  ImgSum = dm1->ccdData(0);

  sumImages(ImgSum, date, runNum);
  
  c.Clear();
  c.Divide(2,1);
  c.cd(1);
  ImgSum->Draw("colz");
  double ymax = double(ImgSum->GetYaxis()->GetXmax());
  int nbinsy = ImgSum->GetNbinsY();
  double y2bin = nbinsy/ymax;

  c.cd(2);
  TH1D* profSmooth = ImgSum->ProjectionY("profSmooth");
  profSmooth->Scale(-1);
  profSmooth->Smooth(ntimes);
  profSmooth->ShowPeaks(width,"",amp); // may have to modify peak fitting parameters;

  TList *fxns = profSmooth->GetListOfFunctions(); // get spacers from polymarker list
  TPolyMarker *pm = (TPolyMarker*) fxns->FindObject("TPolyMarker");
  int spacerNum = pm->GetN();

  vector<double> Spacers;
  double tempSpacer;
  for (int i=0; i<spacerNum; i++)
    {
      tempSpacer = *(pm->GetX()+i);
      Spacers.push_back(tempSpacer);
    }
  sort(Spacers.begin(),Spacers.end());

  vector<double> profK0;
  vector<double> profK0a;
  int tRange= 2; // int(spacerNum/2)-1; // 2 is best value for tRange
  int i,j,k;
  for (i=0; i<spacerNum; i++)
    {
      for (j=1; j<=tRange; j++)
	{
	  // if (j==0) continue;
	  if ((i+j)>=spacerNum) continue;
	  for (k=1; k<=tRange; k++)
	    {
		  profK0.push_back(double((Spacers[i+j]-Spacers[i])/k));
		  profK0a.push_back(1);
	    }
	}
    }

  sort(profK0.begin(),profK0.end());
  double kernel=10;
  double val;
  for (i=0; i<profK0.size()-1; i++)
    {
      for (j=i+1; j<=profK0.size(); j++)
	{
	  if (profK0[j]-profK0[i]>kernel) break;
	  val = exp(-pow(((profK0[i]-profK0[j])/kernel),2)*2);
	  profK0a[i]=profK0a[i]+val;
	  profK0a[j]=profK0a[j]+val;
	}
    }

  double k0 = *(max_element(profK0a.begin(),profK0a.end())-profK0a.begin()+profK0.begin());

  vector<double> profK1(spacerNum);
  vector<double> profK1a(spacerNum);
  for (int i=0; i<spacerNum; i++)
    {
      profK1[i]=int(Spacers[i])%int(k0);
      profK1a[i]=1.0;
    }

  sort(profK1.begin(),profK1.end());
  kernel=10;
  for (i=0; i<spacerNum-1; i++)
    {
      for (j=i+1; j<=spacerNum; j++)
	{
	  if (profK1[j]-profK1[i]>kernel) break;
	  val = exp(-pow(((profK1[i]-profK1[j])/kernel),2)*2);
	  profK1a[i]=profK1a[i]+val;
	  profK1a[j]=profK1a[j]+val;
	}
    }

  double k1 = *(max_element(profK1a.begin(),profK1a.end())-profK1a.begin()+profK1.begin());

  cout << "Spacing: " << k0 << "  Offset: " << k1 << endl;

  fname = "SimData/dmtpc_run" + date + "_" + runNum + "_spacers.root"; // open bias file to save spacers
  TFile *file2 = new TFile(fname,"RECREATE");
  TTree *spacerTree = new TTree("spacer","spacers");
  double spacerY;
  spacerTree->Branch("spacerY",&spacerY,"spacerY/D");

  int check;
  vector<double> spacersFin;
  for (double spcr=k1; spcr<=ymax; spcr=spcr+k0)
    {
      check=0;
      for (int i=0; i<Spacers.size(); i++)
	{
	  if (fabs(Spacers[i]-spcr)<(0.1*k0))
	    {
	      spacerY = Spacers[i];
	      spacerTree->Fill();
	      Spacers.erase(Spacers.begin()+i-1);
	      spacersFin.push_back(spacerY);
	      check=1;
	      break;
	    }
	}
      if (check==0)
	{
	  spacerY = spcr;
	  spacersFin.push_back(spacerY);
	  spacerTree->Fill();
	}
    }

  sort(spacersFin.begin(),spacersFin.end());
  //for (int i=0; i<spacersFin.size(); i++)
  //{
  //  cout << spacersFin[i] <<endl;
  //}

  file2->Write();
  file2->Close();

  spacersFin.~vector<double>();
  profK0.~vector<double>();
  profK0a.~vector<double>();
  profK1.~vector<double>();
  profK1a.~vector<double>();
  Spacers.~vector<double>();
  cout << "done" <<endl;

}
