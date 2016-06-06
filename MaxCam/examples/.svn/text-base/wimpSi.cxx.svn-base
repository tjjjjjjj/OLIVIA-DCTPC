//
//
//
//
//
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooArgSet.h"
#include "RooProdPdf.h"
#include "RooMCStudy.h"
#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooFitResult.h"

#include "TCanvas.h"
#include "TH2F.h"
#include "TFeldmanCousins.h"
#include "TF1.h"
#include "TStyle.h"
#include "TNtuple.h"

#include <vector>
using std::vector;
using std::sort;


//////////////////////////////////////////////////////////////
//
//
//     LIKELIHOOD
//
//
//////////////////////////////////////////////////////////////


// model
RooAddPdf *model=0;
// recoil energy
RooRealVar *eRecoil=0;
// yield in CCD
RooRealVar *eCCD=0;
// simulated data
RooDataSet *data=0;
// fit result
RooFitResult *fitres=0;
// toy manager
RooMCStudy *toy=0;

// WIMP mass
RooRealVar *MD=0;
// WIMP rate
RooRealVar *nSig=0;
// WIMP pdf
RooGenericPdf *sigPdf=0;
RooGenericPdf *sigCCDPdf=0;

// Background rate
RooRealVar *nBkg=0;



void buildModel() {

    // exposure
    //double years=1.;
    //double kilos=100.;

    // energy range
    double eMin=20; // keV
    double eMax=100; //keV    

    // observable: recoil energy
    eRecoil = new RooRealVar("eRecoil", "Recoil energy", 100, eMin, eMax, "keV");
    // observable: recoil energy in CCD
    eCCD = new RooRealVar("eCCD", "Recoil energy in CCD", 100, 0, 65536, "keV");


    
    // signal PDF
    //
    //  dR / dE = R0/(E0*r) * exp( -E / (E0*r) )
    //
    //
    // where
    //        
    //  R0 =  (503/MD*MT) * (sigma/1pb) * (rho/0.4 GeV/cm3) * ( v0/230km/s ) 
    //
    //  r = 4*MD*MT/(MD+MT)^2  
    //
    //  E0 = 1/2*MD*v0^2  
    //
    //
    //
    
    MD = new RooRealVar("MD","WIMP mass", 60e6,  "keV/c^{2}"); 
    RooRealVar *MT = new RooRealVar("MT","Target mass", 28e6, "keV/c^{2}"); 
    RooRealVar *v0 = new RooRealVar("v0","WIMP velocity spread", 220, "km/s");
    RooFormulaVar *r = new RooFormulaVar("r", "Kinematic scattering variable", "4*@0*@1/((@0+@1)*(@0+@1))", RooArgList( *MD, *MT ) );
    RooFormulaVar *E0 = new RooFormulaVar("E0", "WIMP energy spread", "0.5*@0/9e10*@1*@1", RooArgList( *MD, *v0 ) );
    sigPdf    = new RooGenericPdf("sigPdf",    "signal PDF",     "exp(-@0/(@1*@2))", RooArgSet( *eRecoil, *E0, *r));

    
    // quenching factor
    
    RooFormulaVar *Q = new RooFormulaVar("Q", "Quenching in Si", "0.25", RooArgList() ); // need to put realistic values
    RooFormulaVar *W = new RooFormulaVar("W", "Work function in Si in keV", "0.0036", RooArgList() ); // need to put realistic values
    sigCCDPdf = new RooGenericPdf("sigCCDPdf", "CCD signal PDF", "exp(-@0/(@1*@2)*@3/@4)", RooArgSet( *eCCD, *E0, *r, *W, *Q));

    E0->Print();
    r->Print();
    
    // signal rate
    //RooRealVar *sigma=new RooRealVar("sigma","WIMP-nucleon scattering cross-section", 1e-8, "pb");
    //RooRealVar *rho=new RooRealVar("rho","WIMP density", 0.4, "GeV/cm3");
    //RooFormulaVar *R0= new RooFormulaVar("R0", "WIMP recoil rate", "503/(@0*@1)*(@2/1)*(@3/0.4)*(@4/230)", RooArgList(*MD, *MT, *sigma, *rho, *v0) );
    nSig=new RooRealVar("nSig","Number of WIMP recoils", 100 );


    
    
    // background PDF
    //
    // Assume exponential in recoil energy 
    // and flat in cosine of recoil angle
    //
    //  dR / dE = B*beta * exp( -beta*E )
    //
    //  where B = backgroundRatePerKeV * deltaE * years * targetMass
    //
    
    RooRealVar *beta = new RooRealVar("beta","", -0.2); 
    RooExponential *bkgPdf= new RooExponential("bkgPdf","", *eRecoil, *beta);
    RooExponential *bkgCCDPdf= new RooExponential("bkgCCDPdf","", *eCCD, *beta);

    // background rate
    nBkg = new RooRealVar("nBkg","Number of background events", 0);


    

    //
    // Likelihood model with signal and background PDF's 
    //
    model=new RooAddPdf ("model",  "Recoil PDF", RooArgSet(*sigCCDPdf, *bkgCCDPdf), RooArgList(*nSig, *nBkg) );    
}


void plotRecoilModel() {
    //
    // Plot the model
    //
    TCanvas *cp=new TCanvas; cp->SetTitle("PDF Canvas");
    RooPlot *p2=eRecoil->frame(10);   model->plotOn( p2 );  p2->Draw();
}

void plotCCDModel() {
    //
    // Plot the model
    //
    TCanvas *cp=new TCanvas; cp->SetTitle("PDF Canvas");
    RooPlot *p2=eCCD->frame(100);   model->plotOn( p2 );  p2->Draw();
}




void generateData(TString fname="wimpSi.dat", int nSamples=1) {
    //
    // Create simulated sample based on the model
    //
    if (!model) return;
    if (toy) delete toy;
    toy=new RooMCStudy( *model, RooArgSet(*eRecoil), RooFit::Extended(), RooFit::FitOptions("rhem") );
    cout <<  "generating = " <<nBkg->getVal()<< "  " <<nSig->getVal() << endl;
    toy->generate(nSamples, int(nBkg->getVal()+nSig->getVal()), kTRUE);
    data = (RooDataSet*)toy->genData(0);
    if (fname!="none") data->write(fname.Data());
}


void readData(TString fname="wimpSi.dat") {
    //
    // Read ascii file
    //
    if (!model) return;
    data=RooDataSet::read( fname.Data(), RooArgList(*eRecoil) );
}




void fitData() {
    //
    // Fit data using a model for recoil energy distribution
    // for signal and background events
    //    
    if (!data || ! model) return;
    if (fitres) delete fitres;
    //fitres = model->fitTo(*data, RooFit::Extended(kTRUE), RooFit::Hesse(kTRUE), RooFit::Save(kTRUE), RooFit::PrintLevel(-1) );
    fitres = model->fitTo(*data, RooFit::Extended(kTRUE), RooFit::Hesse(kTRUE), RooFit::Save(kTRUE) );
}



void plotFit() {
    //
    // Plot the model
    //
    TCanvas *cp=new TCanvas; cp->SetTitle("PDF Canvas");
    RooPlot *p1=data->plotOn( eRecoil->frame(20) );
    model->plotOn( p1, RooFit::ProjWData(*data) );
    model->plotOn( p1, RooFit::Components("bkgPdf"), RooFit::LineStyle(kDashed));
    p1->Draw("h");

    TCanvas *cl=new TCanvas; cl->SetTitle("NLL Canvas");
    RooPlot *p2=nSig->frame(-3, 30, 50); model->plotNLLOn( p2, data, "L", 0.01, kTRUE );
    p2->Draw();

    TCanvas *cc=new TCanvas; cc->SetTitle("NLL countours");
    RooPlot *p3=new RooPlot(*nBkg, *nSig, 0, 25,  -3, 30);
    fitres->plotOn(p3, "nBkg", "nSig", "ME12");
    p3->Draw();
    
}


double likelihoodScan(bool doPlot=false) {

    assert (model && data);
    
    
    vector<double> x,y,L;
    double lkSum=0;
    for (double ns=0; ns<35; ns+=0.25) {
        x.push_back(ns);
        nSig->setVal(ns);
        nSig->setConstant();
        fitData();
        y.push_back(fitres->minNll());
        double lk=exp(-fitres->minNll());
        L.push_back(lk);
        if (ns>=0) { lkSum+=lk; }
    }
    cout << "lkSum = " << lkSum << endl;
    if (doPlot) {
        TCanvas *cc=new TCanvas; cc->SetTitle("NLL scan");
        TGraph *g=new TGraph(x.size(), &x[0], &y[0]);
        g->SetLineColor(kBlue);
        g->SetLineWidth(2);
        g->Draw("APL");
        g->GetHistogram()->SetXTitle("#mu");
        g->GetHistogram()->SetYTitle("-log(L)");
        TCanvas *cL=new TCanvas; cL->SetTitle("Likelihood scan");
        TGraph *gL=new TGraph(x.size(), &x[0], &L[0]);
        gL->SetLineColor(kBlue);
        gL->SetLineWidth(2);
        gL->Draw("APL");
        gL->GetHistogram()->SetXTitle("#mu");
        gL->GetHistogram()->SetYTitle("L");
    }
    
    unsigned int i=0;
    double sum=0;
    for ( ;i<L.size(); i++) {
        if (x[i]<0) continue;
        //cout << x[i]<< "  " << L[i] << "  " << lkSum<<endl;
        sum += L[i];
        if (sum/lkSum>0.9) {  cout << "mu90="<< x[i] << endl; break;}
    }

    return x[i];
}




