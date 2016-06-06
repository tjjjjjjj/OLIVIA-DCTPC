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
#include "../MaxCamCluster.hh"
#include "../MaxCamImageTools.hh"
#include "../DmtpcDataset.hh"
#include <fstream>
#include <iostream>
#include <vector>

void tracks(int run1, int run2, int speed);
TString makeTree(int run1,int run2);
void plots(TString fileName);
void go();

TString basedir = "/Users/jbattat/research/dmtpc/cameras/pixis1024b/";

void tracks(int run1,int run2,int speed){
  TCanvas *c1=new TCanvas("c1","c1",0,0,800,700);
  c1->cd(1);
  for(int i=run1;i<run2;i++){
    TString fname="/export/dmtpcdata01/analysis/v1/skim/dmtpc_run";
    TString fname2="data/dmtpc_run";
    TString fname3="/export/dmtpcdata01/analysis/v1/sort/dmtpc_run";
    if (i<100){fname+="000";fname2+="000";fname3+="000";}
    else if (i<1000){fname+="00";fname2+="00";fname3+="00";}
    else if (i<10000){fname+="0";fname2+="0";fname3+="0";}
    fname+=i;
    fname2+=i;
    fname3+=i;
    fname+="skim.root";
    fname2+=".root";
    fname3+="sksort.root";
    TFile *f=TFile::Open(fname);
    TFile *f2=TFile::Open(fname3);
    TTree *skim=(TTree*)f->Get("skim");
    TTree *sort=(TTree*)f2->Get("sksort");
    Int_t nev=skim->GetEntries();
    cout << "number of events: " << nev << endl;
    cout<<endl;cout<<"run:"<<"\t"<<i<<endl;cout<<endl;
    double E[2][15],x[2][15],y[2][15],range[2][15],phi[2][15],worm[2][15];
    int ntracks[2],eventnum;
    bool edge[2][15],spark[2];
    skim->SetBranchAddress("E",E);
    skim->SetBranchAddress("ntracks",ntracks);
    skim->SetBranchAddress("edge",edge);
    skim->SetBranchAddress("x",x);
    skim->SetBranchAddress("y",y);
    skim->SetBranchAddress("phi",phi);
    skim->SetBranchAddress("range",range);
    skim->SetBranchAddress("eventnum",&eventnum);
    skim->SetBranchAddress("spark",spark);
    sort->SetBranchAddress("worm",worm);
    DmtpcDataset *d=new DmtpcDataset;
    d->openRootFile(fname2);
    DmtpcDataset *bd=new DmtpcDataset;
    bd->openRootFile("data/dmtpc_run00320.root");
    TH2F *bias=(TH2F*)bd->getBiasFrame(1)->Clone("bias");
    TH2F *image;
    for(Int_t j=0;j<nev;j++){
      skim->GetEvent(j);
      sort->GetEvent(j);
      if(spark[0]==0&&worm[0][0]==0&&edge[0][0]==0&E[0][0]>0&&E[0][0]<100000&&range[0][0]>0&&ntracks[0]==1){
      	cout<<eventnum<<":"<<"\t"<<"(x,y)="<<"\t"<<"("<<x[0][0]<<","<<y[0][0]<<")"<<endl;
	cout<<"\t"<<"E="<<"\t"<<E[0][0]<<endl;
	cout<<"\t"<<"range="<<"\t"<<range[0][0]<<endl;
	cout<<"\t"<<"phi="<<"\t"<<phi[0][0]<<endl;
	cout<<"\t"<<"worm="<<"\t"<<worm[0][0]<<endl;
	cout<<"\t"<<"edge="<<"\t"<<edge[0][0]<<endl;
	cout<<"\t"<<"spark="<<"\t"<<spark[0]<<endl;
	cout<<endl;
	for(int k=0;k<speed;k++){
	  d->getEvent(eventnum);
	  image=d->event()->ccdData(0);
	  image->Add(bias,-1);
	  image->SetMinimum(-100);
	  image->SetMaximum(100);
	  image->Draw("colz");
	  c1->Update();
	}
      }
    }
    f->Close();
    f2->Close();
  }
}



TString makeTree(int run1,int run2){
  //TString outFileName="/export/dmtpcdata01/analysis/v1/skim/skimFile_runs";
  TString outFileName = basedir+"skim/skimFile_runs";
  outFileName+=run1;
  outFileName+="_";
  outFileName+=run2;
  outFileName+=".root";
  cout << "outFileName = " << outFileName << endl;
  TFile *outFile=new TFile(outFileName,"RECREATE");
  TTree *outTree=new TTree("skim","Skimmed Events");
  double E[2][15],range[2][15],x[2][15],y[2][15],phi[2][15],theta[2][15];
  double integral[2],skewness[2][15], worm[2][15];
  bool edge[2][15],spark[2];
  int ntracks[2],date,time,eventnum,run;
  outTree->Branch("run",&run,"run/I");
  outTree->Branch("eventnum",&eventnum,"eventnum/I");
  outTree->Branch("date",&date,"date/I");
  outTree->Branch("time",&time,"time/I");
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
  outTree->Branch("integral",&integral,"integral[2]/D");
  for(int i=run1;i<=run2;i++){
    cout<<"run:"<<"\t"<<i<<endl;
    TString fname  = basedir+"skim/dmtpc_run";
    TString fname3 = basedir+"sort/dmtpc_run";
    //if (i<100){fname+="000";fname3+="000";}
    //else if (i<1000){fname+="00";fname3+="00";}
    //else if (i<10000){fname+="0";fname3+="0";}
    //fname+=i;
    //fname3+=i;
    TString iStr = TString::Format("%05d",i);
    fname  += iStr;
    fname3 += iStr;
    fname+="skim.root";
    fname3+="sksort.root";
    TFile *f=TFile::Open(fname);
    TFile *f2=TFile::Open(fname3);
    TTree *skimIn=(TTree*)f->Get("skim");
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
    cout << "event number  ";
    for(int j=0;j<nev;j++){
      cout << j << " ... ";
      skimIn->GetEvent(j);
      sort->GetEvent(j);
      outFile->cd();
      outTree->Fill();
      gROOT->cd();
    }
    cout << endl;
    f->Close();
    f2->Close();
  }
  outFile->cd();
  outTree->Write();
  delete outTree;
  outFile->Close();
  return outFileName;
}



void plots(TString fileName){
  TFile *file=TFile::Open(fileName);
  TTree *tree=(TTree*)file->Get("skim");
  int nev=tree->GetEntries();
  cout << "number of events: " << nev << endl;
  double E[2][15],range[2][15],x[2][15],y[2][15],phi[2][15],theta[2][15];
  double integral[2],skewness[2][15], worm[2][15];
  bool edge[2][15],spark[2];
  int ntracks[2],date,time,eventnum,run;
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("eventnum",&eventnum);
  tree->SetBranchAddress("date",&date);
  tree->SetBranchAddress("time",&time);
  tree->SetBranchAddress("spark",spark);
  tree->SetBranchAddress("worm",worm);
  tree->SetBranchAddress("ntracks",ntracks);
  tree->SetBranchAddress("edge",edge);
  tree->SetBranchAddress("E",E);
  tree->SetBranchAddress("range",range);
  tree->SetBranchAddress("x",x);
  tree->SetBranchAddress("y",y);
  tree->SetBranchAddress("theta",theta);
  tree->SetBranchAddress("phi",phi);
  tree->SetBranchAddress("integral",integral);
  TH1F *dNdE=new TH1F("dNdE","",100,0,10000);
  dNdE->SetXTitle("Energy (ADU)");
  dNdE->SetYTitle("Number of Events");
  dNdE->GetYaxis()->SetTitleOffset(1.2);
  TH1F *dNdE=new TH1F("h2","",100,0,10000);
  h2->SetXTitle("Number of Tracks");
  h2->SetYTitle("Number of Events");
  h2->GetYaxis()->SetTitleOffset(1.2);
  for(int i=0;i<nev;i++){
    cout << i << "  ";
    tree->GetEvent(i);
    for(int j=0;j<2;j++){
      for(int k=0;k<15;k++){
	//if(spark[j]==0&&E[j][k]>0&&E[j][k]<100000&&edge[j][k]==0&&range[j][k]>2&&worm[j][k]==0){
	if (E[j][k]>0) {
	  dNdE->Fill(E[j][k]);
	  h2->Fill(
	//h3->Fill(E[j][k]/2.5);
	  //h4->Fill(E[j][k]/2.5,phi[j][k]);
	  //h5->Fill(phi[j][k]);
	  //h6->Fill(skewness[j][k]*0.143,E[j][k]/2.5);
	  //h7->Fill(x[j][k]*0.143,y[j][k]*0.143);
	}
      }
    }
  }
  //cout<<h1->GetEntries()<<endl;
  TCanvas *c1=new TCanvas("c1","c1",0,0,1200,1000);
  //c1->Divide(3,3);
  //c1->cd(1);
  dNdE->Draw("p");
  c1->Update();
}


void go() {
  //TString ff = makeTree(0,0);
  TString ff = basedir+"skim/skimFile_runs0_0.root";
  plots(ff);
  return;
}
