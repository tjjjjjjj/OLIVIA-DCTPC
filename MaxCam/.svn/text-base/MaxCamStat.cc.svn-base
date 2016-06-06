#include "MaxCamStat.hh"

#include <iostream>
#include <fstream>
using std::cout;
using std::endl;
using std::cerr;
#include <vector>
using std::vector;


#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TRandom.h"
#include "TMath.h"

//NamespaceImp(MaxCamStat)



void
MaxCamStat::calculateError(int Nneg, int Ntot, double &eup, double &elo, ValueType type) {
    // Given Nneg out of Ntot observed events, compute Poisson errors
    // on fraction Nneg/Ntot
    
    double Npos=Ntot-Nneg;
    int nn=10000;

    double nnegMin = Nneg-10*sqrt(Nneg);
    if (nnegMin<0) nnegMin=0;
    double nnegMax = Nneg+10*sqrt(Nneg);
    if (nnegMax<10) nnegMax=10;
    
    double nposMin = Npos-10*sqrt(Npos);
    if (nposMin<0) nposMin=0;
    double nposMax = Npos+10*sqrt(Npos);
    if (nposMax<10) nposMax=10;
    
    TH1F *htmp;
    if (type==MaxCamStat::Fraction) htmp = new TH1F("htmp","",1000,0,1);
    else if (type==MaxCamStat::Asymmetry) htmp = new TH1F("htmp","",1000,-1,1);
    else assert(0);
    for (int i=0;i<nn;i++) {
        double npos = nposMin + (nposMax-nposMin) * gRandom->Rndm();
        while (TMath::Poisson(Npos, npos)/TMath::Poisson(Npos,Npos)<gRandom->Rndm()) 
            npos = nposMin + (nposMax-nposMin) * gRandom->Rndm();
        
        double nneg = nnegMin + (nnegMax-nnegMin) * gRandom->Rndm();
        while (TMath::Poisson(Nneg, nneg)/TMath::Poisson(Nneg,Nneg)<gRandom->Rndm()) 
            nneg = nnegMin + (nnegMax-nnegMin) * gRandom->Rndm();
        
        if (type==MaxCamStat::Fraction) htmp->Fill( nneg/(npos+nneg) );
        else if (type==MaxCamStat::Asymmetry)  htmp->Fill( (npos-nneg)/(npos+nneg) );
        else assert(0);
        
    }
    
    double f;
    if (type==MaxCamStat::Fraction) f = Nneg/(Npos+Nneg);
    else if (type==MaxCamStat::Asymmetry)  f = (Npos-Nneg)/(Npos+Nneg);
    else assert(0);
    int cen = htmp->FindBin( f );
    int ilow=1, ihigh=htmp->GetNbinsX();
    while ( ilow<cen && htmp->Integral(1,ilow)/htmp->Integral(1,cen) < 0.317 ) ilow++;
    while ( ihigh>cen && htmp->Integral(ihigh,htmp->GetNbinsX())/htmp->Integral(cen,htmp->GetNbinsX())<0.317 ) ihigh--;

    eup = htmp->GetBinLowEdge(++ihigh)-f;
    elo= f-htmp->GetBinLowEdge(ilow);

    cout << "Value type  :   " << type << endl;
    cout << "Value       :   " << f << endl;
    cout << "N1 vs N2        " << Nneg << ", " << Npos << endl;
    cout << "This error  :  +" << eup << " / -" << elo << endl;
    cout << "Binomial:   :   " << sqrt(f*(1-f)/(Npos+Nneg)) << endl;
    //eup = elo = sqrt(f*(1-f)/(Npos+Nneg));

    htmp->Draw();
    delete htmp;
}



TGraphAsymmErrors*
MaxCamStat::fractionError(TH1F* h1, TH1F* h2) {
    int n=h1->GetNbinsX();
    double eu, el;
    vector<double> x,y,xl,xu,yl,yu;
    for (int i=1; i<=n; i++) {
        calculateError((int)h1->GetBinContent(i), (int)h2->GetBinContent(i), eu, el, MaxCamStat::Fraction);
        x.push_back( h1->GetBinCenter(i) );
        double f=h1->GetBinContent(i)/h2->GetBinContent(i);
        y.push_back( f );
        xl.push_back( h1->GetBinWidth(i)*0.5 );
        xu.push_back( h1->GetBinWidth(i)*0.5 );
        yl.push_back( el );
        yu.push_back( eu );
    }
    return new TGraphAsymmErrors(n, &x[0],&y[0], &xl[0],&xu[0], &yl[0],&yu[0]);
}


TGraphAsymmErrors*
MaxCamStat::asymmetryError(TH1F* h1, TH1F* h2) {
    int n=h1->GetNbinsX();
    double eu, el;
    vector<double> x,y,xl,xu,yl,yu;
    for (int i=1; i<=n; i++) {
        calculateError((int)h1->GetBinContent(i), (int)h2->GetBinContent(i), eu, el, MaxCamStat::Asymmetry);
        x.push_back( h1->GetBinCenter(i) );
        double f=1-2*h1->GetBinContent(i)/h2->GetBinContent(i);
        y.push_back( f );
        xl.push_back( h1->GetBinWidth(i)*0.5 );
        xu.push_back( h1->GetBinWidth(i)*0.5 );
        yl.push_back( el );
        yu.push_back( eu );
    }
    return new TGraphAsymmErrors(n, &x[0],&y[0], &xl[0],&xu[0], &yl[0],&yu[0]);
}



