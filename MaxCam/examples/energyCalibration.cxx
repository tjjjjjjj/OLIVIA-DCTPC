#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "../DmtpcDataset.hh"
#include <fstream>
#include <iostream>

void E(){
  int n=0;
  DmtpcDataset *d1=new DmtpcDataset;
  d1->openRootFile("dmtpc_run00005.root");
  TH2F *bias1=(TH2F*)d1->getBiasFrame(0)->Clone("bias");
  TH2F *image;
  TH1D *projX1,*projX3,*projX4,*projY1;
  TH1D *projX2=new TH1D("projX2","projX2",12,0,12);
  TH1D *projY2=new TH1D("projY2","projY2",256,0,256);
  TH1F *hE=new TH1F("hE","",100,-1000,5000);
  TH1F *hE2=new TH1F("hE2","",50,5000,25000);
  TH1F *hE3=new TH1F("hE3","",50,-10000,10000);
  for(int i=0;i<3000;i++){
    Double_t integral=0,mean=0,bg=0;
    Int_t ntracks=-1;
    d1->getEvent(i);
    image=d1->event()->ccdData(0);
    image->Add(bias1,-1);
    projX1=image->ProjectionX("projX1",27,47);
    projX3=image->ProjectionX("projX3",3,22);
    for(int j=0;j<101;j++){mean+=projX3->GetBinContent(150+j);}
    bg=mean/(100*20);
    for(int j=0;j<12;j++){projX2->SetBinContent(j,-21*bg+projX1->GetBinContent(223+j));}
    for(int j=0;j<12;j++){integral+=projX2->GetBinContent(j);}
    cout<<i<<"\t"<<integral<<endl;
    hE->Fill(integral);
    if(integral<200){
      Double_t E0=0;
      ntracks=0;
      projX4=image->ProjectionX("projX4",0,255);
      for(int j=0;j<50;j++){E0+=-256*bg+projX4->GetBinContent(188+j);}
      hE3->Fill(E0);
      cout<<E0<<endl;
    }
    else if(integral>500&&integral<900){
      Double_t totalE=0;
      ntracks=1;
      projX4=image->ProjectionX("projX4",0,255);
      for(int j=0;j<78;j++){totalE+=-256*bg+projX4->GetBinContent(178+j);}
      hE2->Fill(totalE);
      cout<<totalE<<endl;
      projY1=image->ProjectionY("projY1",188,238);
      if(n==0){
	for(int j=0;j<256;j++){projY2->SetBinContent(j,-50*bg+projY1->GetBinContent(j));}
      }
      else if(n>0){
	for(int j=0;j<256;j++){
	  Double_t bin=projY2->GetBinContent(j);
	  projY2->SetBinContent(j,bin-50*bg+projY1->GetBinContent(j));
	}
      }
      n+=1;
    }
  }
  TCanvas *c1=new TCanvas("c1","",0,0,1200,800);
  c1->Divide(3,2);
  c1->cd(1);
  hE->Draw();
  c1->cd(2);
  hE2->Draw();
  c1->cd(3);
  hE3->Draw();
  c1->cd(4);
  projY2->Draw();
  c1->Update();
  Double_t Esum=0;
  for(int j=20;j<248;j++){Esum+=projY2->GetBinContent(j);}
  cout<<Esum<<endl;
  projY2->SetBinContent(45,45000);
  projY2->SetBinContent(46,47000);
  projY2->SetBinContent(47,49000);
  projY2->SetBinContent(48,51000);
  projY2->SetBinContent(49,53000);
  projY2->SetBinContent(50,54000);
  projY2->SetBinContent(51,55000);
  projY2->SetBinContent(78,75000);
  projY2->SetBinContent(79,76000);
  projY2->SetBinContent(80,77000);
  projY2->SetBinContent(81,78000);
  projY2->SetBinContent(82,79000);
  projY2->SetBinContent(83,80000);
  projY2->SetBinContent(84,81000);
  projY2->SetBinContent(85,82500);
  projY2->SetBinContent(86,82500);
  projY2->SetBinContent(116,89000);
  projY2->SetBinContent(117,90000);
  projY2->SetBinContent(118,91000);
  projY2->SetBinContent(119,92000);
  projY2->SetBinContent(120,93000);
  projY2->SetBinContent(121,94000);
  projY2->SetBinContent(122,95000);
  projY2->SetBinContent(151,90000);
  projY2->SetBinContent(152,93000);
  projY2->SetBinContent(153,91000);
  projY2->SetBinContent(154,94000);
  projY2->SetBinContent(155,92000);
  projY2->SetBinContent(156,93000);
  projY2->SetBinContent(157,95000);
  projY2->SetBinContent(185,88000);
  projY2->SetBinContent(186,87000);
  projY2->SetBinContent(187,85000);
  projY2->SetBinContent(188,84000);
  projY2->SetBinContent(189,82000);
  projY2->SetBinContent(190,81000);
  projY2->SetBinContent(191,80000);
  c1->cd(5);
  projY2->Draw();
  Double_t Esum2=0;
  for(int j=20;j<248;j++){Esum2+=projY2->GetBinContent(j);}
  cout<<Esum2<<endl;
  c1->Update();
  cout<<n<<endl;
}


void E2(Int_t run,Int_t nevent,Int_t x1,Int_t x2,Int_t y1,Int_t y2){
  TString fname="data/dmtpc_run";
  if (run<100) fname+="000";
  else if (run<1000) fname+="00";
  else if (run<10000) fname+="0";
  fname+=run;
  fname+=".root";
  TH2F *bias,*bias2,*bias3,*bias4,*bias5;
  TH2F *image;
  TH1D *projX,*projX2;
  DmtpcDataset *d=new DmtpcDataset;
  DmtpcDataset *d2=new DmtpcDataset;
  DmtpcDataset *d3=new DmtpcDataset;
  DmtpcDataset *d4=new DmtpcDataset;
  DmtpcDataset *d5=new DmtpcDataset;
  Double_t integral=0;
  TH2F *image2=new TH2F("image2","",256,0,256,256,0,256);
  d->openRootFile(fname);
  d2->openRootFile(fname);
  d3->openRootFile(fname);
  d4->openRootFile(fname);
  d5->openRootFile(fname);
  bias=d->getBiasFrame(1);
  d->getEvent(nevent);
  d2->getEvent(nevent+1);
  d3->getEvent(nevent-1);
  d4->getEvent(nevent+2);
  d5->getEvent(nevent-2);
  image=d->event()->ccdData(0);
  bias2=d2->event()->ccdData(0);
  bias3=d3->event()->ccdData(0);
  bias4=d4->event()->ccdData(0);
  bias5=d5->event()->ccdData(0);
  image->Add(bias,-1);
  bias2->Add(bias,-1);
  bias3->Add(bias,-1);
  bias4->Add(bias,-1);
  bias5->Add(bias,-1);
  for(int i=0;i<257;i++){
    for(int j=0;j<257;j++){
      Double_t bin=0;
      bin=image->GetBinContent(i,j);
      bin+=-0.25*bias2->GetBinContent(i,j);
      bin+=-0.25*bias3->GetBinContent(i,j);
      bin+=-0.25*bias4->GetBinContent(i,j);
      bin+=-0.25*bias5->GetBinContent(i,j);
      image2->SetBinContent(i,j,bin);
    }
  }
  TCanvas *c1=new TCanvas("c1","",0,0,1200,500);
  c1->Divide(2,1);
  c1->cd(1);
  image2->SetMaximum(150);
  image2->SetMinimum(-50);
  image2->Draw("colz");
  c1->Update();
  projX=image2->ProjectionX("projX",y1,y2);
  for(int i=x1;i<x2+1;i++){integral+=projX->GetBinContent(i);}
  cout<<"Run"<<"\t"<<fname<<endl;
  cout<<"\t"<<nevent<<endl;
  cout<<"\t"<<integral<<endl;
  cout<<endl;
  projX2=image2->ProjectionX("projX2",0,256);
  c1->cd(2);
  projX2->Draw();
  c1->Update();
}


void E3(){
  E2(115,1599,13,100,137,175);
  E2(115,2675,181,206,150,187);
  E2(118,1588,162,237,87,200);
  E2(120,1933,25,225,150,237);
  E2(123,2509,25,225,12,150);
  E2(125,469,187,225,25,212);
  E2(125,536,12,175,150,175);
  E2(128,1387,50,87,87,107);
  E2(133,442,62,131,100,200);
  E2(148,666,100,200,150,225);
  E2(145,1804,75,100,175,237);
  E2(142,1506,125,225,50,100);
  E2(138,942,175,250,175,225);
  E2(144,2069,175,212,125,150);
}


void E4(){
  TH2F *h=new TH2F("h","",50,0,70000,50,0,70000);
  h->Fill(5038,16000);
  h->Fill(1627.75,7300);
  h->Fill(2209,7800);
  h->Fill(3422,7500);
  h->Fill(21322.2,62400);
  h->Fill(14923.8,59100);
  h->Fill(14245.8,45700);
  h->Fill(905,5500);
  h->Fill(9592,39700);
  h->Fill(6587.75,33800);
  h->Fill(1317.5,9600);
  h->Fill(8102,26600);
  h->Fill(383.75,6000);
  h->Fill(2859.75,9700);
  h->SetXTitle("E' (counts)");
  h->SetYTitle("E (counts)");
  TCanvas *c1=new TCanvas("c1","c1",0,0,800,800);
  c1->cd(1);
  h->GetXaxis()->SetTitleColor(1);
  h->GetYaxis()->SetTitleOffset(1.2);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.8);
  h->Draw("p");
  c1->Update();
}
