//
//  Measurement of transverse diffusion
//
//  Analysis:
//    .L diffMeasurement.cxx++
//    init("run_file.root")
//    run()
//  or
//    .L diffMeasurement.cxx++
//    all()
// 
//
#include "../DmtpcEvent.hh"
#include "../DmtpcDataset.hh"
#include "../MaxCamSRIM.hh"

#include "TH2F.h"
#include "TH1D.h"
#include "TF1.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TChain.h"
#include "TLegend.h"
#include "TDatime.h"
#include "TProfile.h"
#include "TRandom2.h"

#include <iostream>
using std::cout;
using std::endl;

#include <vector>

////////////////////////////////////////////////////////////////////////////
//
// Analysis of diffusion data
//
////////////////////////////////////////////////////////////////////////////



DmtpcDataset *d=0; // dataset
TH2F* biasFrame=0;  // bias frame for subtraction
TF1* fun; // fit function for segment profile
TF1* fDiff; // fit function for D/mu
TNtuple *nt; // ntuple with results
TCanvas *call;

// length calibration
float lengthCalibration=1.4e-2; // pixel -> cm

// energy calibration
float energyCalibration=20; // ADC->keV

//
// Range-of-interest (ROI)
//
// ROI: vertical strip
int minXBin=210;
int maxXBin=211;
// down stream segement
int minXBin_ff=130;
int maxXBin_ff=131; 
// ROI: y-bin ranges for individual sources 
int sourceRanges[]={
    20,   60,
    60,  100,
    90,  130,
    140, 180,
    180, 220
};
float *sourceWidth, *sourceWidthErr;
float *sourceWidthSq, *sourceWidthSqErr, *sourceWidthReso;
float *sourceYield, *sourceYieldErr, *sourceYieldReso, *sourceYieldResoErr;

// Coordinate of fishing line intersecting track
int fishingLineLocation[]={
    212,
    212,
    213,
    213,
    213
};

// specify heights of sources (in cm)
float sourceZ[]={ 5.1,  13.1, 17.1,    9.1,  2.1};
float sourceZErr[]={0.1, 0.1, 0.1, 0.1, 0.1};
    
// number of sources
int nSources=5;

// length of drift volume
float driftL=19; // cm

// drift voltage (to be read in from data file)
float driftHV=0;

// pressure (to be read from data file)
float pressureCF4=0;

// temperature (tp be read from data file)
float temperatureCF4=0; 

TDatime *t0=0;

// select straight track only?
bool doStraightSelect=false;

// maximum yield for source
float sourceYieldMax=5000;
int sourceYieldBins=200;

const float sqrtTwoPi=sqrt(TMath::TwoPi());

void
init(TString filename) {
    if (d) { delete d; d=0; }
    d = new DmtpcDataset;
    d->openRootFile(filename);
    biasFrame=d->getBiasFrame(1);
    if (nt) { delete nt; nt=0; }
    nt= new TNtuple("nt", "source activities", "source:y:yErr:mean:meanErr:width:widthErr:anodeHV:driftHV:press:temp0:mean_ff");
    if (fun) { delete fun; fun=0; }
    fun=new TF1("fun","[0]*exp(-0.5*((x-[1])/[2])**2)+[3]"); // segment profile
    if (fDiff) { delete fDiff; fDiff=0; }
    fDiff = new TF1("fDiff","[0]**2+2*[1]*[2]*x");
    fDiff->SetParName(0, "#sigma_{0}");
    fDiff->SetParName(1, "D/#mu");
    fDiff->SetParName(2, "L/V");
    if (call) { delete call; call=0; }
    call= new TCanvas("call","",800,600);
    call->Divide(5,3);
    sourceWidth=new float[nSources];
    sourceWidthErr=new float[nSources];
    sourceWidthSq=new float[nSources];
    sourceWidthSqErr=new float[nSources];
    sourceWidthReso=new float[nSources];
    sourceYield=new float[nSources];
    sourceYieldErr=new float[nSources];
    sourceYieldReso=new float[nSources];
    sourceYieldResoErr=new float[nSources];
}

void getEvent(int i) {
    d->getEvent(i);
}


void fitEvent() {
    
    // read amplification, drift voltages (in kV) and pressure (in Torr)
    float anode=d->event()->experimentConfig("anodeHV")->currentValue;
    float drift=d->event()->experimentConfig("driftHV")->setValue;
    float press=d->event()->experimentConfig("pressure")->currentValue;
    float temp0=d->event()->experimentConfig("temp0")->currentValue;

    if (fabs(driftHV-drift)>1e-3) {
        cout << "Change in setup !!!!" << endl;
        cout << "NEW DRIFT VOLTAGE = " << drift << endl;
        driftHV=drift;
    }
    if (fabs(pressureCF4-press)>1e-3) {
        cout << "Change in setup !!!!" << endl;
        cout << "NEW PRESSURE = " << press << endl;
        pressureCF4=press;        
    }

    
    // get and analyze event image
    TH2F* image = d->event()->ccdData(0);
    if (biasFrame) image->Add(biasFrame, -1);
    TH1D* projY_ff = doStraightSelect ? image->ProjectionY("projY_ff", minXBin_ff, maxXBin_ff) : 0;
    TH1D* projY = image->ProjectionY("projY", minXBin, maxXBin);

    int nEdges=nSources*2;
    for (int i=0; i<nEdges; i+=2) {
        projY->GetXaxis()->SetRange(sourceRanges[i], sourceRanges[i+1]);
        fun->SetParameter(0, projY->GetMaximum());
        fun->SetParameter(1, projY->GetBinCenter(projY->GetMaximumBin()));
        fun->SetParameter(2, 2);
        fun->SetParameter(3, projY->GetMinimum());

        if (doStraightSelect) { // select only tracks that are parallel with x axis
            projY_ff->Fit("fun", "0Q");
        }
        projY->Fit("fun", "0Q");
        //projY->Fit("fun"); call->Update(); getchar();
        float y    = sqrtTwoPi * fun->GetParameter(0) * fun->GetParameter(2) / projY->GetBinWidth(1);
        float yErr = y*sqrt(pow(fun->GetParError(0)/fun->GetParameter(0),2)+
                          pow(fun->GetParError(2)/fun->GetParameter(2),2));
        nt->Fill( i/2,
                  y, yErr,
                  fun->GetParameter(1),fun->GetParError(1),
                  fun->GetParameter(2),fun->GetParError(2),
                  anode, drift, press, temp0,
                  doStraightSelect ? projY_ff->GetFunction("fun")->GetParameter(1) : fun->GetParameter(1)
                  );
    }
    
}




void sourceEval(int iSource) {
    TH2F *hWY=new TH2F("hWY","", 50,0,15, sourceYieldBins,0,sourceYieldMax );
    
    call->cd(iSource+1);
    
    // find maximum bin - assume it's 1-track bin
    TH1F hmean("hmean","",100,0,1024);
    TString cut="source=="; cut+=iSource;
    cut+="&&y*width>"; cut+=300*(maxXBin-minXBin+1)*pressureCF4/75.;
    nt->Draw("mean>>hmean", cut);
    float mean=hmean.GetBinCenter( hmean.GetMaximumBin() );
    cut+="&&mean>"; cut+=mean-20; cut+="&&mean<"; cut+=mean+20;
    cut+="&&abs(mean-mean_ff)<4";
    
    nt->Draw("y:width>>hWY",cut,"goff");
    hWY->DrawCopy("colz");
    
    TH1D *hy=hWY->ProjectionY();
    int  ybin=hy->GetMaximumBin();

    call->cd(nSources+iSource+1);
    TH1D *hw=hWY->ProjectionX("hw", ybin-3, ybin+2);
    hw->Fit("gaus", "QL"); hw->DrawCopy(); // Q
    
    call->cd(nSources*2+iSource+1);
    int  wbin=hw->GetMaximumBin();
    hy=hWY->ProjectionY("hy", wbin-7, wbin+7);
    hy->GetXaxis()->SetRange(ybin-3, ybin+2);
    hy->Fit("gaus", "Q"); // Q
    hy->GetXaxis()->SetRange(1, hy->GetNbinsX());
    hy->DrawCopy();

    TF1 *gfit=hw->GetFunction("gaus");
    sourceWidthSq[iSource]=pow(gfit->GetParameter(1),2);
    sourceWidthSqErr[iSource]=2*gfit->GetParameter(1)*gfit->GetParError(1);
    sourceWidthReso[iSource]=gfit->GetParameter(2)/gfit->GetParameter(1);
    delete hw;
    
    gfit=hy->GetFunction("gaus");
    sourceYield[iSource]    =gfit->GetParameter(1);
    sourceYieldErr[iSource] =gfit->GetParError(1);
    sourceYieldReso[iSource]=gfit->GetParameter(2)/gfit->GetParameter(1);
    sourceYieldResoErr[iSource] = sourceYieldReso[iSource]*sqrt( pow(gfit->GetParError(1)/gfit->GetParameter(1),2) +
                                                                 pow(gfit->GetParError(2)/gfit->GetParameter(2),2) );
    delete hy;
    
    delete hWY;
}



void run() {

    int nev=d->tree()->GetEntries();
    for (int i=0; i<nev; i++) {
        getEvent(i);
        fitEvent();
    }
    for (int i=0; i<nSources; i++) sourceEval(i);

    TH1F *hpress=new TH1F("hpress","",100,0,250);
    nt->Draw("press>>hpress","","goff");
    pressureCF4=hpress->GetMean();
    delete hpress;
    TH1F *htemp0=new TH1F("htemp0","",100,10,30);
    nt->Draw("temp0>>htemp0","","goff");
    temperatureCF4=htemp0->GetMean()+273.15;
    delete htemp0;
    TDatime *t=d->event()->timeStamp();
    if (!t0) { t0=new TDatime; *t0=*t; }
    float dtime = (t->Convert()-t0->Convert())/(24.*3600.);
    TH1F *hanode=new TH1F("hanode","",1000,0.5,1);
    nt->Draw("anodeHV>>hanode","","goff");
    float anodeHV=hanode->GetMean();
    delete hanode;
    
    
    int runno=d->event()->runNumber();

    // Fit D/mu
    driftHV*=1e3; // kV -> V
    float EN = driftHV/driftL  /  MaxCamSRIM::numberDensity(pressureCF4, temperatureCF4); // in V cm2

    // calibrate length
    for (int i=0; i<nSources; i++) {
        sourceWidthSq[i]   *=lengthCalibration*lengthCalibration;
        sourceWidthSqErr[i]*=lengthCalibration*lengthCalibration;
    }

    // set fit parameters
    fDiff->SetParameter(2, driftL/driftHV);
    fDiff->SetParLimits(2, driftL/driftHV, driftL/driftHV);
    
    TCanvas *cdiff=new TCanvas; cdiff->cd();
    TGraphErrors *gdiff = new TGraphErrors(nSources, sourceZ, sourceWidthSq, sourceZErr, sourceWidthSqErr);
    gdiff->Draw("AP");
    gdiff->GetHistogram()->SetXTitle("z (cm)");
    gdiff->GetHistogram()->SetYTitle("#sigma^{2} (cm^{2})");
    gdiff->Fit("fDiff", "QME");
    float Dmu   =gdiff->GetFunction("fDiff")->GetParameter(1);
    float DmuErr=gdiff->GetFunction("fDiff")->GetParError(1);
    float sigma0   =gdiff->GetFunction("fDiff")->GetParameter(0);
    float sigma0Err=gdiff->GetFunction("fDiff")->GetParError(0);
    float corr=0;//gdiff->GetFunction("fDiff")->GetCorrelation(0,1);
   
    cout << "Run number ...... " << runno << endl;
    cout << "Anode voltage ... " << anodeHV << "kV" << endl;
    cout << "Drift voltage ... " << driftHV << " kV" << endl;
    cout << "Pressure ........ " << pressureCF4 << " Torr" << endl;
    cout << "Temperature ..... " << temperatureCF4 << " K" << endl;
    cout << "E/N ............. " << EN << "  V cm2" << endl; 
    cout << "D/mu ............ " << Dmu << "+/-" << DmuErr << " V" << endl;
    cout << "Sigma0 .......... " << sigma0 << "+/-" << sigma0Err << " cm" << endl;
    cout << "Correlation ..... " << corr << endl;
    

    // output results to file
    TFile* fout = new TFile("diffMeasurement.root","UPDATE");
    TNtuple *nout=(TNtuple*)fout->Get("nout");
    if (!nout) nout = new TNtuple("nout","","EN:Dmu:DmuErr:sig0:sig0Err:press:temp0:minx:maxx:run:dtime:anode:corr");
    nout->Fill(EN, Dmu, DmuErr, sigma0, sigma0Err, pressureCF4, temperatureCF4, minXBin, maxXBin,runno, dtime, anodeHV, corr);
    fout->Write();
    fout->Close();
    delete fout;
}



void all() {
    for (int i=692; i<=694; i++) {
        TString fname="data/dmtpc_run00";
        fname += i;
        fname += ".root";
        init(fname);
        run();
    }
}


void nominal() {
    vector<int> runList;
    //for (int i=671; i<=681; i++) runList.push_back(i); // 75Torr
    //for (int i=683; i<=686; i++) runList.push_back(i); // 150Torr
    for (int i=695; i<=704; i++) runList.push_back(i); // 50Torr

    for (unsigned int i=0; i<runList.size(); i++) {
        TString fname="data/dmtpc_run00";
        fname += runList[i];
        fname += ".root";
        init(fname);
        run();        
    }
    
}

void roiSize(TString fname) {
    for (int i=160; i<=180; i+=5) { //160-240
        minXBin=i;
        for (int j=0; j<10; j++) {
            maxXBin=i+j;
            init(fname);
            run();
        }
    }
}


void runSingleSource(TString fname) {
    nSources=1;
    init(fname);
    TFile* fout = new TFile("diffMeasurement.root","UPDATE");
    if (nt) { delete nt; nt=0; }
    nt= new TNtuple("nt", "source activities", "source:y:yErr:mean:meanErr:width:widthErr:anodeHV:driftHV:press:temp0:mean_ff");
    int nev=d->tree()->GetEntries();
    for (int i=0; i<nev; i++) {
        getEvent(i);
        fitEvent();
    }
    sourceEval(0);
    nt->Write();
    fout->Write();
}


void runCollimationStudy() {
    doStraightSelect=false; init("data/dmtpc_run00685.root"); run();
    doStraightSelect=true;  init("data/dmtpc_run00685.root"); run();
}





////////////////////////////////////////////////////////////////////////////
//
// Plots
//
////////////////////////////////////////////////////////////////////////////


void plotDmu(TString opt="") {

// Recommended values from Christophorou et al.
    float ENRec[]={0.14, 0.17, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.80,
                   1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 35, 40, 45, 50,
                   60, 70, 80, 90, 100, 150, 200, 250, 300};
    float DmuRec[]={0.025, 0.026, 0.027, 0.028, 0.029, 0.030, 0.031, 0.032, 0.033,
                    0.034, 0.035, 0.035, 0.036, 0.037, 0.039, 0.041, 0.046, 0.052,
                    0.062, 0.084, 0.155, 0.293, 0.492, 0.736, 1.01, 1.29, 1.58, 2.16,
                    2.68, 3.12, 3.49, 3.81, 4.78, 5.16, 5.29, 5.39 };
    TGraph *gDmuRec=new TGraph(36, ENRec, DmuRec);
    
    
// Measured values by Curtis et al. 1988
    float ENCur[]={0.14, 0.17, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 1.0,
                   1.2, 1.4, 1.7, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10, 11,
                   12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};
    float DmuCur[]={0.0258, 0.0256, 0.0263, 0.0274, 0.0287, 0.0296, 0.0312, 0.0324, 0.0333,
                    0.0338, 0.0343, 0.0348, 0.0352, 0.0346, 0.0343, 0.0348, 0.0354, 0.0362,
                    0.0365, 0.0382, 0.0406, 0.0431, 0.0468, 0.0479, 0.0494, 0.0538, 0.0567,
                    0.0617, 0.0650, 0.0723, 0.0790, 0.0846, 0.0952, 0.105,  0.114,  0.130,
                    0.144,  0.160};
    TGraph *gDmuCur=new TGraph(38, ENCur, DmuCur);
    gDmuCur->SetMarkerStyle(24);
    gDmuCur->SetMarkerColor(kBlack);
    gDmuCur->SetMarkerSize(1);

    // Measured values by Schmidt&Polenz
    float SchPol[]={0.175119,	0.0292471,	
        0.284223,	0.0356001,	
        0.351872,	0.0350513,	
        0.717083,	0.0389358,	
        0.887799,	0.0395175,	
        1.08364,	0.0407224,	
        1.26737,	0.0413347,	
        1.38034,	0.0407062,	
        1.44032,	0.0360481,	
        1.4823,	0.0432497,	
        1.50339,	0.0400872,	
        1.56909,	0.0419522,	
        1.61419,	0.0382979,	
        1.75804,	0.0371472,	
        1.83491,	0.0394702,	
        1.99844,	0.0382844,	
        2.20793,	0.0394582,	
        2.5456,	0.0382691,	
        3.24287,	0.0406489,	
        4.37289,	0.0412504,	
        5.11461,	0.0438217,	
        6.51555,	0.0465468,	
        8.06848,	0.0549876,	
        10.7277,	0.0611244};

    int nSchPol=sizeof(SchPol)/sizeof(float);
    float *ENSch=new float[nSchPol];
    float *DmuSch=new float[nSchPol];
    for (int i=0; i<nSchPol; i+=2) { ENSch[i]=SchPol[i]; DmuSch[i]=SchPol[i+1]; }
    
    
    /*float ENSch[]={7.72581, 3.90586, 2.27017, 1.24459, 0.96299, 0.634162,
                   0.494341, 0.402653, 0.250125, 0.129396, 0.100148, 0.0613063 };
    float DmuSch[]={0.65729, 0.266194, 0.141528, 0.0836252, 0.0683541, 0.049855,
    0.0444936, 0.0397054, 0.03363, 0.0297741, 0.0303169, 0.025454}; */
    TGraph *gDmuSch=new TGraph(nSchPol, ENSch, DmuSch);
    gDmuSch->SetMarkerStyle(25);
    gDmuSch->SetMarkerSize(1);
    gDmuSch->SetMarkerColor(kBlack);
     

    
// Measured vales by Lakshminarasimha et al. 1973
    float ENLak[]={15.5,  31};
    float DmuLak[]={0.1, 0.5 };
    TGraph *gDmuLak=new TGraph(2, ENLak, DmuLak);
    gDmuLak->SetMarkerStyle(26);
    gDmuLak->SetMarkerSize(1);
    gDmuLak->SetMarkerColor(kBlack);
    
    
// Our measured values from diffMeasurement.root file
    TChain nout("nout");
    nout.Add("DIFFUSION/diffMeasurement_run695-704_1.15.root"); // 50Torr
    nout.Add("DIFFUSION/diffMeasurement_run671-681_1.15.root"); // 75Torr
    nout.Add("DIFFUSION/diffMeasurement_run683-686_1.15.root"); // 150Torr

    TH2F* frame=new TH2F("frame","",100,0.5,20, 100,0,0.15);
    frame->SetXTitle("E/N (10^{-17} V cm^{2})");
    frame->SetYTitle("D_{T}/#mu (V)");
    frame->Draw();
    gPad->SetLogx();
    gDmuCur->Draw("P");
    //gDmuRec->Draw("C");
    gDmuLak->Draw("P");
    gDmuSch->Draw("P");
    if (opt.Contains("prof")) {
        TProfile *hDmu=new TProfile("hDmu","", 14, 0.5, 14.5);
        hDmu->SetMarkerStyle(21);
        hDmu->SetLineColor(kRed);
        hDmu->SetMarkerColor(kRed);
        //nout.Draw("Dmu:EN*1e17>>hDmu","","prof same");
        nout.Draw("Dmu:EN*1e17>>hDmu","","prof");
        vector <float> xDmu,yDmu;
        for (int i=1; i<=hDmu->GetNbinsX(); i++) {
            xDmu.push_back( hDmu->GetBinCenter(i) );
            yDmu.push_back( hDmu->GetBinContent(i) );
        }
        TGraph *gDmu = new TGraph(xDmu.size(), &xDmu[0], &yDmu[0]);
        gDmu->SetMarkerStyle(21);
        gDmu->SetLineColor(kRed);
        gDmu->SetMarkerColor(kRed);
        frame->Draw();
        gDmuCur->Draw("P");
        //gDmuRec->Draw("C");
        gDmuLak->Draw("P");
        gDmuSch->Draw("P");
        gDmu->Draw("P");
        
        TLegend *leg=new TLegend(0.3,0.6, 0.75,0.85);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->AddEntry(gDmuLak, "Lakshminarasimha 1973", "P");
        leg->AddEntry(gDmuCur, "Curtis 1988", "P");
        leg->AddEntry(gDmuSch, "Schmidt 1988", "P");
        //leg->AddEntry(gDmuRec, "Recommended (Cristophorou et al.)", "L");
        leg->AddEntry(hDmu, "DMTPC", "P");
        leg->Draw();

        for (int ip=1; ip<=hDmu->GetNbinsX(); ip++) {
            cout << hDmu->GetBinCenter(ip) << "  &  " << hDmu->GetBinContent(ip) << endl;
        }
        
        if (opt.Contains("sigma")) {
            gPad->Clear();
            float z=20; // cm
            vector<float> sEN;
            vector<float> sig;
            vector<float> sig_dmtpc;
            float N=MaxCamSRIM::numberDensity(75)*1e-17; //Torr
            for (float EN=0.5; EN<40.5; EN+=1.0) {
                float E=EN*N;
                float sig_d_sq = gDmuRec->Eval(EN) * 2 * z / E; // cm2
                float sig_0_sq = 0.05*0.05; // cm2
                float sig_d_sq_dmtpc = hDmu->GetBinContent( hDmu->FindBin(EN) ) * 2 * z / E; // cm2
                sEN.push_back( EN );
                sig.push_back      ( sqrt(sig_0_sq + sig_d_sq       ) );
                sig_dmtpc.push_back( sqrt(sig_0_sq + sig_d_sq_dmtpc ) );
                //cout << EN << "  " << E << "  " << sig_d_sq << "  " << sig_0_sq << "  " << sig_d_sq_dmtpc << endl;
            }
            TGraph *gSig=new TGraph(sEN.size(), &sEN[0], &sig[0]);
            TGraph *gSig_dmtpc=new TGraph(15, &sEN[0], &sig_dmtpc[0]);
            gSig_dmtpc->SetMarkerColor(kRed);
            gSig->SetLineWidth(4);
            gSig->SetLineStyle(2);
            gSig->Draw("AC");
            gSig->GetHistogram()->SetXTitle("E/N (10^{-17} V cm^{2})");
            gSig->GetHistogram()->SetYTitle("#sigma  (cm)");
            gSig_dmtpc->Draw("P");
        }

    } else if (opt.Contains("press")) {
        nout.SetMarkerStyle(25);
        nout.SetMarkerColor(kBlack);
        nout.Draw("Dmu:EN*1e17","abs(press-75)<0.5");
        TGraph *gDmu_75 = (TGraph*)gPad->GetPrimitive("Graph")->Clone("gDmu_75");
        gDmu_75->SetMarkerStyle(24);
        nout.Draw("Dmu:EN*1e17","abs(press-150)<0.5");
        TGraph *gDmu_150 = (TGraph*)gPad->GetPrimitive("Graph")->Clone("gDmu_150");
        gDmu_150->SetMarkerStyle(25);
        nout.Draw("Dmu:EN*1e17","abs(press-50)<0.5");
        TGraph *gDmu_50 = (TGraph*)gPad->GetPrimitive("Graph")->Clone("gDmu_50");
        gDmu_50->SetMarkerStyle(26);
        frame->Draw();
        gDmu_50->Draw("P");
        gDmu_75->Draw("P");
        gDmu_150->Draw("P");
        TLegend *leg=new TLegend(0.3,0.6, 0.85,0.85);
        leg->SetFillColor(0);
        leg->SetBorderSize(0);
        leg->AddEntry(gDmu_50, "50 Torr", "P");
        leg->AddEntry(gDmu_75, "75 Torr", "P");
        leg->AddEntry(gDmu_150, "150Torr", "P");
        leg->Draw();
    }

    
}



void plotSigma0() {
    TChain nout("nout");
    nout.Add("DIFFUSION/diffMeasurement_run695-704_1.15.root"); // 50Torr
    nout.Add("DIFFUSION/diffMeasurement_run671-681_1.15.root"); // 75Torr
    nout.Add("DIFFUSION/diffMeasurement_run683-686_1.15.root"); // 150Torr

    TH2F* frame=new TH2F("frame","",100,0,18, 100,0,0.2);
    frame->SetXTitle("E/N (10^{-17} V cm^{2})");
    frame->SetYTitle("#sigma_{0} (cm)");

    

    nout.SetMarkerStyle(25);
    nout.SetMarkerColor(kBlack);
    nout.Draw("sig0:EN*1e17","abs(press-75)<0.5");
    TGraph *gSig0_75 = (TGraph*)gPad->GetPrimitive("Graph")->Clone("gSig0_75");
    gSig0_75->SetMarkerStyle(24);
    nout.Draw("sig0:EN*1e17","abs(press-150)<0.5");
    TGraph *gSig0_150 = (TGraph*)gPad->GetPrimitive("Graph")->Clone("gSig0_150");
    gSig0_150->SetMarkerStyle(25);
    nout.Draw("sig0:EN*1e17","abs(press-50)<0.5");
    TGraph *gSig0_50 = (TGraph*)gPad->GetPrimitive("Graph")->Clone("gSig0_50");
    gSig0_50->SetMarkerStyle(26);
    frame->Draw();
    gSig0_50->Draw("P");
    gSig0_75->Draw("P");
    gSig0_150->Draw("P");
    TLegend *leg=new TLegend(0.3,0.6, 0.85,0.85);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(gSig0_50, "50 Torr", "P");
    leg->AddEntry(gSig0_75, "75 Torr", "P");
    leg->AddEntry(gSig0_150, "150Torr", "P");
    leg->Draw();
   
}




void plotRoiSize(TString what="Dmu") {
    TChain nout("nout");
    nout.Add("DIFFUSION/diffMeasurement_150Torr_run685_roi.root");
    TCanvas *c1 = new TCanvas( "c1" ); 
    TH2F * frame= new TH2F("frame","",10,-0.5,9.5,  17,160-2.5,240+2.5);
    frame->SetMinimum(0);
    frame->SetXTitle("Segment length (bin)");
    frame->SetYTitle("Segment position (bin)");
    if (what=="Dmu") { nout.Draw("minx:maxx-minx>>frame","Dmu*(maxx-minx>-1)", "", nout.GetEntries()-1, 1);  c1->SetPhi(60);}
    else { nout.Draw("minx:maxx-minx>>frame","sig0*(maxx-minx>-1)", "",
                     nout.GetEntries()-1, 1);  c1->SetPhi(210); c1->SetTheta(10); }
   
    frame->Draw("lego2");

    float sumDmu=0, sumDmu2=0;
    int nDmu=0;
    for (int i=230; i>=190; i-- ) {
        for (int j=0; j<4; j++) {
            float val = frame->GetBinContent( frame->FindBin(j,i) );
            sumDmu += val;
            sumDmu2 += val*val;
            nDmu++;
        }
    }
    sumDmu /= nDmu;
    sumDmu2 /= nDmu;
    float rmsDmu=sqrt( sumDmu2 - sumDmu*sumDmu );
    cout << "mean ...." << sumDmu << endl;
    cout << "RMS ....." << rmsDmu << endl;
}


void plotFocusing(TString what="points") {
    TChain nout("nout");
    nout.Add("DIFFUSION/diffMeasurement_focus.root");
    nout.Add("DIFFUSION/diffMeasurement_run00665.root");

    if (what=="points") {
        TH2F * frame= new TH2F("frame","",  5,-0.5,4.5,  100,0.02, 0.06);
        frame->GetXaxis()->SetBinLabel(1,"3");
        frame->GetXaxis()->SetBinLabel(2,"3|1");
        frame->GetXaxis()->SetBinLabel(3,"3|2");
        frame->GetXaxis()->SetBinLabel(4,"3|3");
        frame->GetXaxis()->SetBinLabel(5,"3|4");
        frame->SetXTitle("Lens focus dial");
        frame->SetYTitle("D/#mu (V), or #sigma_{0} (cm)");
        nout.Draw("Dmu:run");
        TGraph *gDmu = (TGraph*)gPad->GetPrimitive("Graph")->Clone("gDmu");
        gDmu->SetMarkerStyle(25);
        // change run number -> number on focus dial
        gDmu->GetX()[0]=4;
        gDmu->GetX()[1]=0;
        gDmu->GetX()[2]=1;
        gDmu->GetX()[3]=3;
        gDmu->GetX()[4]=2;
        gDmu->GetX()[5]=2;
        
        nout.Draw("sig0:run");
        TGraph *gSigma0 = (TGraph*)gPad->GetPrimitive("Graph")->Clone("gSigma0");
        gSigma0->GetX()[0]=4;
        gSigma0->GetX()[1]=0;
        gSigma0->GetX()[2]=1;
        gSigma0->GetX()[3]=3;
        gSigma0->GetX()[4]=2;
        gSigma0->GetX()[5]=2;
        
        frame->Draw();
        gDmu->Draw("P");
        gSigma0->Draw("P");
        
        TLegend *leg=new TLegend(0.5,0.3, 0.75,0.4);
        leg->AddEntry(gDmu, "D/#mu", "P");
        leg->AddEntry(gSigma0, "#sigma_{0}", "P");
        leg->SetFillColor(0);
        leg->SetBorderSize(0);
        leg->Draw();
    }
    else {
        TH1F *hDmu=new TH1F("hDmu","",100,0.03,0.08);
        nout.Draw("Dmu>>hDmu");
        hDmu->Draw();
        hDmu->Fit("gaus", "L");
    }
}


void plotTimeDependence() {
    TChain nout("nout");
    nout.Add("DIFFUSION/diffMeasurement_time_dependence.root");
    TH2F * frame= new TH2F("frame","",  100,-0.5,3,  100,0.02, 0.06);
    frame->SetXTitle("Time (days)");
    frame->SetYTitle("D/#mu (V), or #sigma_{0} (cm)");
    nout.Draw("Dmu:dtime");
    TGraph *gDmu = (TGraph*)gPad->GetPrimitive("Graph")->Clone("gDmu");
    gDmu->SetMarkerStyle(25);
    nout.Draw("sig0:dtime");
    TGraph *gSigma0 = (TGraph*)gPad->GetPrimitive("Graph")->Clone("gSigma0");
    
    frame->Draw();
    gDmu->Draw("P");
    gSigma0->Draw("P");

    TLegend *leg=new TLegend(0.5,0.3, 0.75,0.4);
    leg->AddEntry(gSigma0, "#sigma_{0}", "P");
    leg->AddEntry(gDmu, "D/#mu", "P");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();
    
}


void plotAnodeDependence() {
    TChain nout("nout");
    nout.Add("DIFFUSION/diffMeasurement_anode.root");
    TH2F * frame= new TH2F("frame","",  6,0.695,0.755,  100,0.02, 0.06);
    frame->SetXTitle("Anode (kV)");
    frame->SetYTitle("D/#mu (V), or #sigma_{0} (cm)");
    nout.Draw("Dmu:anode");
    TGraph *gDmu = (TGraph*)gPad->GetPrimitive("Graph")->Clone("gDmu");
    gDmu->SetMarkerStyle(25);
    gDmu->Print();
    nout.Draw("sig0:anode");
    TGraph *gSigma0 = (TGraph*)gPad->GetPrimitive("Graph")->Clone("gSigma0");
    frame->Draw();
    gDmu->Draw("P");
    gSigma0->Draw("P");

    TLegend *leg=new TLegend(0.5,0.3, 0.75,0.4);
    leg->AddEntry(gSigma0, "#sigma_{0}", "P");
    leg->AddEntry(gDmu, "D/#mu", "P");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();
}



void plotStruggling() {
    // data points from Timur
    float x[]={0, 0.5, 1, 1.5};
    float y[]={0, 171, 226, 258};
    TGraph *g= new TGraph(5, x, y);
    g->Draw("AC");
}




void energyResolutionStudy(TString fname) {

    vector<float> x, xerr, y, yerr;
    MaxCamSRIM he("SRIM_He_in_CF4_100Torr");
    he.setPressure(75);

    float sourceX=243;

    sourceYieldMax=20000;
    sourceYieldBins=200;
    
    sourceRanges[0]=110; sourceRanges[1]=180; // 4x4 binning
    doStraightSelect=true;
    minXBin_ff=55;
    maxXBin_ff=70; 

    
    maxXBin=220;
    for (int i=217; i>216/*180*/; i-=5) {   
        minXBin=i;
        runSingleSource(fname);
        x.push_back( he.calcEnergyLoss( 5350,
                                        (sourceX-maxXBin)*40*lengthCalibration, 
                                        (sourceX-minXBin+1)*40*lengthCalibration) );
        xerr.push_back( 0 );
        y.push_back( sourceYieldReso[0] );
        yerr.push_back( sourceYieldResoErr[0] );
        
    }

    TCanvas *cEres=new TCanvas("cEres");
    cEres->cd();
    TFile* fout = new TFile("diffMeasurement.root","RECREATE");
    TGraphErrors *g=new TGraphErrors(x.size(), &x[0], &y[0], &xerr[0], &yerr[0] );
    g->Draw("AP");
    g->GetHistogram()->SetXTitle("E (keV)");
    g->GetHistogram()->SetYTitle("#DeltaE/E");
    g->Write();
    fout->Write();
    TF1 *fres = new TF1("fres","sqrt([0]**2+[1]**2/x)", 10, 500);
    fres->SetParNames("E_{0}", "E_1");
    g->Fit("fres");
}

void plotEnergyResolution() {
    TH2F *frame=new TH2F("frame","",100, 0, 500, 100, 0, 0.17);
    frame->SetXTitle("E (keV)");
    frame->SetYTitle("#sigma_{E}/E");
    frame->Draw();
    TFile *f=new TFile("DIFFUSION/diffMeasurement_run778_energyReso.root");
    TGraphErrors *g=(TGraphErrors*)f->Get("Graph");
    g->Draw("P");

    float xray_E[]={5.9, 5350};
    float xray_reso[]={0.1, 0.02};
    TGraph *gxray=new TGraph(2,xray_E, xray_reso);
    gxray->SetMarkerColor(kRed);
    gxray->SetMarkerStyle(29);
    gxray->SetMarkerSize(2);
    gxray->Draw("P");

}


void spatialResolutionStudy(TString fname) {

    vector<float> x, xerr, y, yerr;
    MaxCamSRIM he("SRIM_He_in_CF4_100Torr");
    he.setPressure(200);

    float sourceX=980;

    sourceRanges[0]=530; sourceRanges[1]=570; // 1x1 binning
    
    
    maxXBin=801;
    minXBin=801;
    maxXBin_ff=900;
    minXBin_ff=905;
    doStraightSelect=true;
    
    float deltaE=he.calcEnergyLoss( 5350,
                                    (sourceX-maxXBin)*40*lengthCalibration, 
                                    (sourceX-minXBin+1)*40*lengthCalibration);
    runSingleSource(fname);

    x.push_back( deltaE );
    xerr.push_back( 0 );
    y.push_back( sourceWidth[0] );
    yerr.push_back( sourceWidthErr[0] );
    

    TCanvas *cEres=new TCanvas("cEres");
    cEres->cd();
    TGraphErrors *g=new TGraphErrors(x.size(), &x[0], &y[0], &xerr[0], &yerr[0] );
    g->Draw("AP");
    g->GetHistogram()->SetXTitle("E (keV)");
    g->GetHistogram()->SetYTitle("#DeltaE/E");
    TF1 *fres = new TF1("fres","sqrt([0]**2+[1]**2/x)", 10, 500);
    fres->SetParNames("E_{0}", "E_1");
    g->Fit("fres");
}


void plotSpatialResolution() {

    TChain chnt("nt");
    chnt.Add("DIFFUSION/diffMeasurement_run785.root");
    chnt.Add("DIFFUSION/diffMeasurement_run787.root");
    chnt.Add("DIFFUSION/diffMeasurement_run788.root");

    vector<float> x,y, xerr,yerr;
    TH1F *hwidth=new TH1F("hwidth","",100,0,5);
    chnt.Draw("width>>hwidth","width*y>150&&y<500&&width<3.5&&y>0&&abs(press-200)<1");
    hwidth->Fit("gaus");
    x.push_back(200);
    xerr.push_back(1);
    y.push_back(hwidth->GetFunction("gaus")->GetParameter(1)*lengthCalibration*10);
    yerr.push_back(hwidth->GetFunction("gaus")->GetParError(1)*lengthCalibration*10);

    chnt.Draw("width>>hwidth","width*y>600&&y<1200&&width<3.5&&y>0&&abs(press-75)<1");
    hwidth->Fit("gaus");
    x.push_back(75);
    xerr.push_back(1);
    y.push_back(hwidth->GetFunction("gaus")->GetParameter(1)*lengthCalibration*10);
    yerr.push_back(hwidth->GetFunction("gaus")->GetParError(1)*lengthCalibration*10);
    
    chnt.Draw("width>>hwidth","width*y>80&&y<250&&width<3.5&&y>0&&abs(press-250)<1");
    hwidth->Fit("gaus");
    x.push_back(250);
    xerr.push_back(1);
    y.push_back(hwidth->GetFunction("gaus")->GetParameter(1)*lengthCalibration*10);
    yerr.push_back(hwidth->GetFunction("gaus")->GetParError(1)*lengthCalibration*10);

    gPad->Clear();
    TGraphErrors *g=new TGraphErrors(x.size(), &x[0],&y[0], &xerr[0],&yerr[0]);
    g->Draw("AP");
}


void fitAttenuation(int sig_run=761, int bkg_run=776) {
    
    TChain ntWF("ntWF");
    TString sigrun=TString("DIFFUSION/plotMCA_run"); sigrun+=sig_run; sigrun+=TString(".root");
    ntWF.Add(sigrun);
    TString bkgrun=TString("DIFFUSION/plotMCA_run"); bkgrun+=bkg_run; bkgrun+=TString(".root");
    ntWF.Add(bkgrun);

    TH1F* hwf_all = new TH1F("hwf_all","", 128,    0-0.039,   10-0.039);
    hwf_all->SetMarkerStyle(1);
    TH1F* hwf_bkg = new TH1F("hwf_bkg","", 128,    0-0.039,   10-0.039);
    hwf_bkg->SetLineStyle(2);
    TH1F* hwf_sig = new TH1F("hwf_sig","", 128,    0-0.039,   10-0.039);

    TCanvas *cWF = new TCanvas("cWF");
    cWF->Divide(2,1);

    cWF->cd(1);

    TString sigcut=TString("run=="); sigcut+=sig_run;
    ntWF.Draw("peak>>hwf_all",sigcut);
    TString bkgcut=TString("run=="); bkgcut+=bkg_run;
    ntWF.Draw("peak>>hwf_bkg",bkgcut);
    
    hwf_bkg->Scale(1./4);
    hwf_all->GetXaxis()->SetRange(10,80);
    hwf_all->Draw("E");
    hwf_bkg->Draw("same");

    cWF->cd(2);
    hwf_sig->Add(hwf_all, hwf_bkg, 1, -1);
    hwf_sig->GetXaxis()->SetRange(10,80);
    hwf_sig->SetMarkerStyle(1);
    hwf_sig->Draw("E");

    TF1 *funFe55=new TF1("funFe55","[0]*exp(-0.5*(x-[1])**2/[2]**2)+[3]*exp(-0.5*(x-[4])**2/[5]**2)");
    funFe55->SetParameters(1320, 1.75, 0.17, 300, 1.2, 0.56);
    hwf_sig->Fit("funFe55");
}



void plotAttenuation() {
    TCanvas *cWF = new TCanvas("cWF","",800, 400);
    cWF->Divide(2,1);

    cWF->cd(1);
    float N=MaxCamSRIM::numberDensity(75);

    // -2600V
    float x[]={2.1, 9.1, 13.1, 17.1};
    float xerr[]={0.1,0.1,0.1};
    float y[]={1.752, 1.750, 1.731, 1.740};
    float yerr[]={0.003, 0.003, 0.003, 0.003};
    TGraphErrors *g= new TGraphErrors(4,x,y, xerr,yerr);
    g->Draw("AP");
    g->GetHistogram()->SetXTitle("E/N (Td)");
    TF1 *fatt=new TF1("fatt","[0]*exp([1]*x)");
    fatt->SetParNames("M_{0}", "#eta");
    g->Fit("fatt");
    vector<float> EN, ENErr, alphaN, alphaNErr;
    EN.push_back(2600./19.5/N*1e17);
    ENErr.push_back(0.1);
    alphaN.push_back(g->GetFunction("fatt")->GetParameter(1)/N*1e18);
    alphaNErr.push_back(g->GetFunction("fatt")->GetParError(1)/N*1e18);
        
    // -5000V
    //     765  1.79107e+00   3.03049e-03
    //     768  1.78179e+00   2.92083e-03
    //     763    1.80040e+00   2.86867e-03
    float y5000[]={1.79107e+00, 1.78179e+00,  1.80040e+00};
    float yerr5000[]={0.003, 0.003, 0.003};
    TGraphErrors *g5000= new TGraphErrors(3,&x[1],y5000, &xerr[1],yerr5000);
    g5000->Draw("P");
    g5000->Fit("fatt");
    EN.push_back(5000./19.5/N*1e17);
    ENErr.push_back(0.1);
    alphaN.push_back(g5000->GetFunction("fatt")->GetParameter(1)/N*1e18);
    alphaNErr.push_back(g5000->GetFunction("fatt")->GetParError(1)/N*1e18);    
    
    
    cWF->cd(2);
    TGraphErrors *gan = new TGraphErrors(2, &EN[0], &alphaN[0], &ENErr[0], &alphaNErr[0]);
    gan->SetMarkerColor(kRed);
    
    float ENrec[]={8, 9, 10, 15, 20, 25, 30, 35, 40};
    float alphaNrec[]={-0.052, -0.056, -0.058, -0.051, -0.043, -0.143, -0.281, -0.438, -0.730};
    TGraph *ganrec= new TGraph( sizeof(ENrec)/sizeof(float) ,ENrec, alphaNrec );
    ganrec->SetLineWidth(3);
    ganrec->SetLineStyle(2);
    ganrec->SetMarkerStyle(28);

    for (unsigned int i=0; i<EN.size(); i++) {
        cout << EN[i] << "  " << alphaN[i] << "+-" << alphaNErr[i] << endl;
    }

    // Masured values by Hunter et al
    float ENhun[]={20, 22, 25, 30, 35, 40};
    float alphaNhun[]={-0.001, -0.0027, -0.0134, -0.093, -0.306, -0.636};
    TGraph *ganhun= new TGraph( sizeof(ENhun)/sizeof(float), ENhun, alphaNhun );
    ganhun->SetMarkerStyle(26);


    // Measured values by Datskos, Carted, Christophorou
    float ENDat[]={2.5, 5, 7.5, 10, 15, 20};
    int nDat=sizeof(ENDat)/sizeof(float);
    for (int i=0; i<nDat; i++) ENDat[i] *=  ( 100 * 760./101325 * 1.38e-23*300 * 1e21); // ->Vcm2
    float alphaNDat[]={-0.0496278, -0.0496278, -0.0496278, -0.248139, -1.21588, -2.48139};
    TGraph *gandat= new TGraph( nDat, ENDat, alphaNDat );
    gandat->SetMarkerStyle(25);
    gandat->SetMarkerColor(kBlue);

    
    TCanvas *cmea=new TCanvas;
    cmea->cd();
    TH2F *hatt=new TH2F("hatt","",100,0, 45,  100,-0.35,0.05);
    hatt->SetXTitle("E/N (Td)");
    hatt->SetYTitle("#eta/N (#times10^{18} cm^{2})");
    hatt->Draw();
    //ganrec->Draw("C");
    gandat->Draw("P");
    gan->Draw("P");
    ganhun->Draw("P");

    TLegend *leg=new TLegend(0.20,0.3, 0.85,0.5);
    leg->SetFillStyle(0);
    leg->AddEntry(ganhun, "Hunter 1987", "P");
    leg->AddEntry(gandat, "Datskos 1992", "P");
    //leg->AddEntry(ganrec, "Recommended (Christophorou et al)", "L");
    leg->AddEntry(gan, "DMTPC", "P");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();


}



void calcRMS(int n=1000, float rangl=0) {
    gStyle->SetOptStat(1);

    TRandom2 *rand= new TRandom2;
    float xsize=256;
    float ysize=256;
    rangl*=TMath::DegToRad();
    float cosr=cos(rangl);
    float sinr=sin(rangl);
    TH2F *hxy=new TH2F("hxy","", int(2*xsize),-xsize,xsize,  int(2*ysize),-ysize,ysize);
    while (n>0) {
        float x = (rand->Uniform()-0.5)*xsize;
        float y = (rand->Uniform()-0.5)*ysize;
        hxy->Fill(x*cosr+y*sinr,y*cosr-x*sinr);
        n--;
    }
    TCanvas *c=new TCanvas;
    c->Divide(2,1);
    c->cd(1); hxy->Draw("box");
    c->cd(2); TH1D* hx=hxy->ProjectionX(); hx->Draw();

    cout <<"x RMS ............. " << hx->GetRMS() << endl;
    cout <<"xsize/sqrt(12) .... " << xsize/sqrt(12) << endl;
}

