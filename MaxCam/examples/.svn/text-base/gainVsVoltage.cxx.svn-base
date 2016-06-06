/************************************************************
//
//
*************************************************************/
#include "TGraphErrors.h"
#include "../MaxCamRead.hh"
#include "../MaxCamChannel.hh"
#include "../MaxCamConfig.hh"
#include "TMath.h"
#include "TF1.h"
#include "TH2F.h"
#include "TTree.h"
#include "TProfile.h"
#include "TPolyMarker.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TLegend.h"
#include <iostream>
#include <vector>


MaxCamRead *ana;
int i, n, nparam;
TF1 *fun;
float x[9];
TGraphErrors *gr=0;


TGraph *gcal;


void calibrateWidth() {
    int i0=1;
    float xx[6], yy[6];
    for (int ii=i0+1; ii<wireList.size(); ii++) {
	int iii=ii-i0-1;
	xx[iii]=wireList[ii]-wireList[i0];
	yy[iii]=(ii-i0)*5000; // spacing in microns
	cout << iii << "  " << xx[iii] << "  " << yy[iii] << endl;
    }
    cout << wireList.size()-i0 << endl;
    gcal=new TGraph(wireList.size()-i0-1,xx,yy);
    gcal->Draw("APL");
    gcal->Fit("pol1");
}



TProfile *hprof[100];
TLegend *leg;
void plotVoltage(TString opt="", int iw) {

    hprof[iw]=new TProfile(TString("hprof")+iw,"",9,2.525,2.975);
    hprof[iw]->SetMarkerColor(iw);
    hprof[iw]->SetMarkerStyle(20+iw);
    hprof[iw]->SetLineColor(iw);

    TString selectWire=TString("wire.y>0&&wire.y<7e6&&wire.width>8&&wire.width<16&&wire.id=="+iw);

    nt->Draw( TString("wire.y:conf.wirehv>>hprof") + iw, selectWire, "prof"+opt);

    leg->AddEntry( hprof[iw], TString("wire ")+(iw+1), "p");
}

void plotVoltageScan() {
  leg = new TLegend(0.3,0.6,0.5,0.8);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  plotVoltage(    "",1);
  //plotVoltage("same",2);
  //plotVoltage("same",3);
  //plotVoltage("same",4);
  //plotVoltage("same",5);

  leg->Draw();
}



void plotTime(TString opt, int iw) {
    hprof[iw]=new TProfile(TString("hprof")+iw,"",20,0,10);
    hprof[iw]->SetMarkerColor(iw+1);
    hprof[iw]->SetMarkerStyle(20+iw);
    hprof[iw]->SetLineColor(iw+1);

    TString selectWire=TString("wire.y>2500&&wire.y<7000&&id==")+iw;
    selectWire+="&&wire.width>8&&wire.width<16");
    selectWire+="&&conf.wirehv>2.82&&conf.wirehv<2.88");

    nt->Draw( TString("wire.y:conf.time/3600>>hprof"), selectWire, "prof"+opt);

    leg->AddEntry( hprof[iw], TString("wire ")+(iw+1), "p");
}

void plotTimeScan() {
  leg = new TLegend(0.3,0.6,0.5,0.8);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  plotTime(    "",0);
  plotTime("same",1);
  plotTime("same",2);
  plotTime("same",3);
  plotTime("same",4);
  
  leg->Draw();
}



void findPressureAndWire(float press, int wid) {
    int nn=nt->GetEntries();
    int ii=0;
    for (; ii<nn; ii++) {
	nt->GetEvent(ii);
	if (fabs(expconf.setpress-press)>1e-3) continue;
	if (wire.id!=wid) continue;
	return;
    }
    if (ii=nn) cout <<"could not find wire=" << wid << " pressure="<< press<<endl; 
}

void findPressure(float press) {
    i=0;
    for (; i<n; i++) {
	ana->getEvent(i);
	if (fabs(ana->pressure()->setValue-press)>1e-3) continue;
	return;
    }
    if (i==n) cout <<"could not find pressure="<< press<<endl; 
}

TGraph *srim;
void stopping(int pressure) {
    TString fileName="/Users/ddujmic/Desktop/DCH/SRIM/data/ioniz_CF4_";
    fileName+=pressure;
    fileName+=TString("Torr_150mrad.dat");
    cout <<fileName<<endl;
    vector<double> xs;
    vector<double> ys;
    double xx, yy, zz;
    ifstream fin(fileName);
    while (!fin.eof()) {
    	fin >> xx >> yy >> zz;
    	xs.push_back(xx*1e-8);
    	ys.push_back(yy*10);
	//cout << xx << "  " << yy << endl;
    }
    srim= new TGraph( xs.size(), &xs[0], &ys[0]);
    TString title="P=";
    title+=pressure;
    title+="Torr";
    srim->SetTitle(title);
}
TGraphErrors *bragg;
vector<float> calib;
float xoffset=2.2;
void plotBraggX(float press, int scaleWire=-1, float sp=-1, int useCalib=0) {
    int nw=5;
    float *xw=new float[nw];
    float *xwerr=new float[nw];
    float *yw=new float[nw];
    float *ywerr=new float[nw];
    if (useCalib<0) calib.clear();
    for (int iw=0; iw<nw; iw++) {
	xw[iw] = iw*0.5+xoffset;
	xwerr[iw]=0;
	findPressureAndWire(press,iw);
	yw[iw] = wire.y;
	ywerr[iw] = wire.yerr;
	if (scaleWire>=0) { 
	    findPressureAndWire(press,scaleWire);
	    if (sp>0) wire.y/=sp;
	    if (useCalib>0) wire.y/=calib[iw]; // apply calibration
	    yw[iw]/= wire.y;
	    ywerr[iw]/=wire.y;
	    if (useCalib<0) {
		calib.push_back( srim->Eval(xw[iw])/yw[iw] ); // compute calibration
		cout <<"wire="<<iw<<" calib="<<srim->Eval(xw[iw])/yw[iw]<<endl;
	    }
	}
    }
    bragg = new TGraphErrors(nw, xw, yw, xwerr, ywerr);
}

TGraph *braggy;
void plotBraggY(float press, float sp=-1) {
    findPressure(press);
    TH1D *yproj = ana->ccdImage()->ProjectionY("yproj",60,65);
    int np=yproj->GetNbinsX();
    float *xp=new float[np];
    float *xperr=new float[np];
    float *yp=new float[np];
    float *yperr=new float[np];
    for (int ii=0; ii<np; ii++) {
	xp[ii] = 36e-4*ii*8+xoffset;
	xperr[ii]=0;
	yp[ii] = yproj->GetBinContent(ii+1);
	yperr[ii]=0;
        if (sp>0) { yp[ii]/=sp; yperr[ii]/=sp; }
    }
    braggy = new TGraph(np, xp, yp);
}

void plotAllBragg() {
    float p[]={135, 180, 195, 210, 280, 345, 410, 455, 500};
    int np=sizeof(p)/sizeof(float);
    c= new TCanvas("c","",650,800);   
    c->Divide(3,3);
    gStyle->SetOptTitle();
    // calibrate
    stopping(135);
    plotBraggX(135, 0, srim->Eval(xoffset), -1);
    
    for (int ii=0; ii<np; ii++) {
    	c->cd(ii+1); 
	stopping(int(p[ii])); 
	srim->Draw("AL");
	srim->GetHistogram()->SetXTitle("depth (cm)");
	srim->GetHistogram()->SetYTitle("dE/dx (MeV/mm)");
	plotBraggX(p[ii],0,srim->Eval(xoffset), 1);
	bragg->Draw("P");
    }
}


void tmpplot98_100() {
    TProfile *hp98=new TProfile("hp98","",5, 2.275,2.525); hp98->SetMarkerStyle(21);
    nt->Draw("wire.y:conf.wirehv>>hp98", "wire.id==0&&wire.y>0&&wire.y<10000&&wire.mean>300&&wire.mean<340&&conf.setpress==49","prof");
    TProfile *hp99=new TProfile("hp99","",5, 2.575,2.825); hp99->SetMarkerStyle(24);
    nt->Draw("wire.y:conf.wirehv>>hp99","wire.id==0&&wire.y>0&&wire.y<10000&&wire.mean>350&&wire.mean<390&&wire.width>6&&wire.width<13&&conf.setpress==110","prof");   
    TProfile *hp100=new TProfile("hp100","",4, 2.625,2.825); hp100->SetMarkerStyle(20);
    nt->Draw("wire.y:conf.wirehv>>hp100","wire.id==0&&wire.y>0&&wire.y<10000&&wire.mean>350&&wire.mean<390&&wire.width>6&&wire.width<13&&conf.setpress==180","prof");
    hh=new TH2F("hh","",100,2.2,2.9, 100, 500,3000);hh->Draw();
    hp100->Draw("same");
    hp99->Draw("same");
    hp98->Draw("same");
}

TF1* dcalib=new TF1("dcalib","2.6-x*36e-4*0.93"); 
TGraphErrors* profile(int nbin, float xmin, float xmax, int wireid=-1) {

    TH1F hpw("hpw","",100,0,pow(15*36e-3,2));
    TH1F hpx("hpx","",100,0.35,2.35);
    float xx[100], yy[100], xerr[100], yerr[100];
    int i=0;
    float dx=(xmax-xmin)/nbin;
    for (float x=xmin; x<xmax; x+=dx) {
    	TString xcut="y>1200&&y<4000&&width<15&&(2.6-mean*36e-4*0.93>";
    	xcut += x;
    	xcut += ")&&(2.6-mean*36e-4*0.93<";
    	xcut += x+dx;
	xcut += ")";
        if (wireid>=0) { xcut += "&&id=="; xcut += wireid; }
	cout << xcut;
 	ana->nt->Draw("2.6-mean*36e-4*0.93>>hpx",xcut);  //hpx->Fit("gaus","L"); gPad->Update(); gSystem->Sleep(1000);
	ana->nt->Draw("(width*36e-3)**2>>hpw",xcut); hpw->Fit("gaus","L"); //hpw->Fit("gaus","L"); gPad->Update(); gSystem->Sleep(1000);
	xx[i]=hpx->GetMean();
	xerr[i]=hpx->GetRMS()/sqrt(hpx->Integral());
	yy[i]=hpw->GetMean();
	yerr[i]=hpw->GetRMS()/sqrt(hpw->Integral());
	i++;
    }
    return new TGraphErrors(i,xx,yy,xerr,yerr);
}


void tmpplot102_105(float depth=1) {
  init();
  ana->findWires(4,0.1,"y");
  ana->makeNtuple(0,3,3,0,"y");

    gStyle->SetTitleOffset(1.6,"y");
    c= new TCanvas("c","",800,400); c->Divide(2,1); c->cd(1);
    TGraphErrors *gr =profile(4, 0.35,2.35);
    TGraphErrors *gr0=profile(4, 0.35,2.35,0); gr0->SetMarkerColor(1);
    TGraphErrors *gr1=profile(4, 0.35,2.35,1); gr1->SetMarkerColor(2);
    TGraphErrors *gr2=profile(4, 0.35,2.35,2); gr2->SetMarkerColor(3);
    TGraphErrors *gr3=profile(4, 0.35,2.35,3); gr3->SetMarkerColor(4);
     gr0->Draw("AP");
    gr1->Draw("P");
    gr2->Draw("P");
    gr3->Draw("P");
    gr0->GetHistogram()->SetXTitle("drift distance (cm)");
    gr0->GetHistogram()->SetYTitle("width_{T}^{2} (mm^{2})");
    c->cd(2);
    gr->Draw("AP");
    gr->GetHistogram()->SetXTitle("drift distance (cm)");
    gr->GetHistogram()->SetYTitle("width_{T}^{2} (mm^{2})");
    gStyle->SetOptFit();
    gr->Fit("pol1");
    float p0=gr->GetFunction("pol1")->GetParameter(0);
    float p0err=gr->GetFunction("pol1")->GetParError(0);
    float p1=gr->GetFunction("pol1")->GetParameter(1);
    float p1err=gr->GetFunction("pol1")->GetParError(1);
    float width = sqrt( gr->GetFunction("pol1")->Eval(depth) );
    float werr = 0.5/width*sqrt( pow(p0err,2) + pow(p1err*depth,2) + pow(p1*0.1,2) );
    cout <<" width at " << depth << "cm  = " << width << " +/- " << werr << endl;
}

void tmpplot116(int wireid=-1, int pressure=-1, TString opt="") {
    //TProfile *hp=new TProfile("hp","",6, 2.675,2.975);
    TProfile *hp=new TProfile("hp","",14, 2.425,3.125);
    TString cut("y>0&&y<10000&&mean>250&&mean<350&&width<20&&width>6");
    if (pressure>0) {
    	cut += TString("&&abs(setpress-");
    	cut += pressure;
    	cut += ")<1e-3";
	hp->SetMarkerColor(pressure/50);
    }
    if (wireid>=0) {
	cut += "&&id==";
	cut += wireid;
	hp->SetMarkerStyle(wireid+20);
	hp->SetMarkerColor(wireid+1);
    }
    cout << cut << endl;
    ana->nt->Draw("y:wirehv>>hp",cut,"goff");
    hp->Draw(opt);
}


void init() {
  ana = new MaxCamRead("data/ccdrun_00139.root");
  ana->readBiasFrame();
  //ana->addRun("data/ccdrun_00103.root");
  //ana->addRun("data/ccdrun_00104.root");
  //ana->addRun("data/ccdrun_00105.root");
  //ana->addRun("data/ccdrun_00120.root");
  //ana->addRun("data/ccdrun_00121.root");

  i=0; 

  fun = new TF1("fun","[0]+[1]*exp(-0.5*(x-[2])**2/[3]**2)"); 

  n=ana->tree()->GetEntries(); 
  cout << "TOTAL = " << n << endl;

  ana->getEvent(0);

  //ana->findHotPixels(500,3,0.05);
  ana->findHotPixels("hotpixels.dat");

}



void gainVsVoltage() {

  init();

  ana->findWires(4,0.1,"y");

  ana->makeNtuple(0,2,3,0,"y");
}
