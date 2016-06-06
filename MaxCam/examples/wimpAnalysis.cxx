#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1.h"
#include "TH2.h"
#include "TNtuple.h"
#include "TROOT.h"
#include "TSystem.h"
#include "/net/cockroach/data02/jocelyn/data_processing_14May2009/projects/DarkMatter/MaxCam/MaxCamCluster.hh"
#include "/net/cockroach/data02/jocelyn/data_processing_14May2009/projects/DarkMatter/MaxCam/MaxCamImageTools.hh"
#include "/net/cockroach/data02/jocelyn/data_processing_14May2009/projects/DarkMatter/MaxCam/DmtpcDataset.hh"
#include "/net/cockroach/data02/jocelyn/data_processing_14May2009/projects/DarkMatter/MaxCam/DmtpcEvent.hh"
#include "/net/cockroach/data02/jocelyn/data_processing_14May2009/projects/DarkMatter/MaxCam/MaxCamChannel.hh"
#include <fstream>
#include <iostream>
#include <vector>

#include "TF1.h"
#include "TH3.h"
#include "TSystem.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TProfile.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TROOT.h"


void wimpAnalysis2(int run1,int run2){

  TString outFileName="/net/cockroach/data02/jocelyn/data_processing_14May2009/analysis_ntuples/analysis_ntuples";
  outFileName+=run1;
  outFileName+="_";
  outFileName+=run2;
  outFileName+=".root";
  TFile *outFile=new TFile(outFileName,"RECREATE");
  TTree *outTree=new TTree("skim","Skimmed Events");
  double E[2][15],range[2][15],x[2][15],y[2][15],phi[2][15],theta[2][15];
  double integral[2],skewness[2][15], worm[2][15], maxpixel[2][15], energy_density[2][15], npixel[2][15], cluster_mean[2][15], cluster_rms[2][15], neighbors[2][15];
  bool edge[2][15],spark[2];
  int ntracks[2],date,time,eventnum,run;
  outTree->Branch("run",&run,"run/I");
  outTree->Branch("eventnum",&eventnum,"eventnum/I");
  outTree->Branch("date",&date,"date/I"); // empty for now
  outTree->Branch("time",&time,"time/I"); // empty for now
  outTree->Branch("spark",&spark,"spark[2]/O");
  outTree->Branch("worm",&worm,"worm[2][15]/D");
  outTree->Branch("ntracks",&ntracks,"ntracks[2]/I");
  outTree->Branch("edge",&edge,"edge[2][15]/O");
  outTree->Branch("E",&E,"E[2][15]/D");
  outTree->Branch("range",&range,"range[2][15]/D");
  outTree->Branch("x",&x,"x[2][15]/D");
  outTree->Branch("y",&y,"y[2][15]/D");
  outTree->Branch("theta",&theta,"theta[2][15]/D");
  outTree->Branch("phi",&phi,"phi[2][15]/D");
  outTree->Branch("maxpixel",&maxpixel,"maxpixel[2][15]/D");
  outTree->Branch("energy_density",&energy_density,"energy_density[2][15]/D");
  outTree->Branch("npixel",&npixel,"npixel[2][15]/D");
  outTree->Branch("integral",&integral,"integral[2]/D");
  outTree->Branch("cluster_mean",&cluster_mean,"cluster_mean[2][15]/D");
  outTree->Branch("cluster_rms",&cluster_rms,"cluster_rms[2][15]/D");
  outTree->Branch("neighbors",&neighbors,"neighbors[2][15]/D");

  int binno=0;  int countNeighbors=0;
  double clusmean=0; double clusrms=0; 
  for (int l=0; l<2; l++) {
    for (int p=0; p<15; p++) {
      maxpixel[l][p] = 0.;
      energy_density[l][p] = 0.;
      neighbors[l][p] = 0.;
    }
  }

  for(int i=run1;i<run2;i++){
    cout<<"run:"<<"\t"<<i<<endl;
    TString fname="/net/cockroach/data02/jocelyn/data_processing_14May2009/skim/dmtpc_run";
    TString fname3="/net/cockroach/data02/jocelyn/data_processing_14May2009/sort/dmtpc_run";
    if (i<100){fname+="000";fname3+="000";}
    else if (i<1000){fname+="00";fname3+="00";}
    else if (i<10000){fname+="0";fname3+="0";}
    fname+=i;
    fname3+=i;
    fname+="skim.root";
    fname3+="sksort.root";
    // add a check that the file exists
    FILE *ftest=fopen(fname,"r");
    FILE *f2test=fopen(fname3,"r");

    // need to exclude runs that don't have both trees
    if ( (ftest!=NULL) && (f2test!=NULL) ) {
      TFile *f=TFile::Open(fname);
      TFile *f2=TFile::Open(fname3);
      TTree *skimIn=(TTree*)f->Get("skim");

      TObjArray* clustarr=0;
      //     TClonesArray* cleanimg=0;
      skimIn->SetBranchAddress("clusters",&clustarr);
      //      skimIn->SetBranchAddress("cleanimage",&cleanimg);

      TTree *sort=(TTree*)f2->Get("sksort");
      Int_t nev=skimIn->GetEntries();
      skimIn->SetBranchAddress("eventnum",&eventnum);
      skimIn->SetBranchAddress("date",&date);
      skimIn->SetBranchAddress("time",&time);
      skimIn->SetBranchAddress("E",E);
      skimIn->SetBranchAddress("range",range);
      skimIn->SetBranchAddress("x",x);
      skimIn->SetBranchAddress("y",y);
      skimIn->SetBranchAddress("phi",phi);
      skimIn->SetBranchAddress("theta",theta);
      skimIn->SetBranchAddress("integral",integral);
      skimIn->SetBranchAddress("skewness",skewness);
      skimIn->SetBranchAddress("edge",edge);
      skimIn->SetBranchAddress("spark",spark);
      skimIn->SetBranchAddress("ntracks",ntracks);
      sort->SetBranchAddress("worm",worm);
      for(int j=0;j<nev;j++){
	skimIn->GetEvent(j);
	sort->GetEvent(j);

	run = i;

	// try to cut double worms with extra variables
	//      loop over cameras
	for (int m=0; m<clustarr->GetEntries(); m++) {
	  if ( (spark[m]==0) ) {
	    MaxCamClusterImage* c = (MaxCamClusterImage*)clustarr->At(m);
	    TH2F* img = c->getImage();
	    // loop over tracks
 	    for(int k=0; k<c->getNCluster(); k++) {
	      if ( (edge[m][k]==0) ) {
		vector <int> px = c->getCluster(k);
		npixel[m][k]=c->getCluster(k).size();
		energy_density[m][k] = c->getEnergyDensity(k);

		// loop over pixels in the cluster
		vector <float> binvals;
		double maxx=0; 
		for(int l=0; l<int(px.size()); l++) {
		  binvals.push_back(img->GetBinContent(px[l]));
		  if (img->GetBinContent(px[l])>maxx) {
		    maxx=img->GetBinContent(px[l]);
		    binno=px[l];
		  }
		  clusmean+=img->GetBinContent(px[l]);
		}
		clusmean /= double(px.size());
		for(int l=0; l<int(px.size()); l++) {
		  clusrms+=TMath::Power((img->GetBinContent(px[l])-clusmean),2);
		}
		clusrms=TMath::Sqrt((clusrms/double(px.size())));
		cluster_mean[m][k] = clusmean;
		cluster_rms[m][k] = clusrms;
	 
		double maxval =  *max_element(binvals.begin(),binvals.end());
		maxpixel[m][k] = maxval;

		// get neighbors around the max pixel
		if (npixel[m][k]<64) {
		  double mean, rms;
		  MaxCamImageTools::meanRMSNoOutliers(img,mean,rms);
		  if (k==0) integral[m] = rms;
		  double threshold=3.7*rms;
		  
		  // get maxpixel position
		  int nbx=img->GetNbinsX()+2;
		  int nby=img->GetNbinsY()+2;
		  int pi = binno%nbx;
		  int pj = ((binno-pi)/nbx)%nby;
		
		  countNeighbors = 0;
		  for (int ii=pi-1; ii<=pi+1; ii++) {
		    for (int jj=pj-1; jj<=pj+1; jj++) {
		      if (ii==pi && jj==pj) continue;
		      if (img->GetBinContent(ii,jj)<threshold) continue;
		      countNeighbors++;
		    }
		  }
		  //		  cout << pi << "\t" << pj << "\t" << threshold << "\t" << countNeighbors << "\t" << clusmean << "\t" << clusrms << endl;
		} // end if npixel<64
		neighbors[m][k] = countNeighbors;

	      } // end if edge==0
	    } // end loop over cluster entries
	  } // end if spark == 0
	} // end loop over cameras
	
	outFile->cd();
	outTree->Fill();
	gROOT->cd();
      }
      f->Close();
      f2->Close();
    } // end if (files) are open
  } // end loop over runs
  outFile->cd();
  outTree->Write();
  delete outTree;
  outFile->Close();
}


