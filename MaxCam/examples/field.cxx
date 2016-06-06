
#include "TF2.h"
#include "TGraph.h"

#include <iostream>
using std::cout;
using std::endl;

#include "../MaxCamUnits.hh"
using MaxCamUnits::cm;
using MaxCamUnits::mm;
using MaxCamUnits::um;
using MaxCamUnits::volt;
#include "../MaxCamElectricField.hh"


MaxCamElectricField *elf=0;

double V(double *var, double *par) {
    if (par);
    return elf->V(var);
}


TH1F* hM=0;
void plotField(double xmin, double xmax, double xstep, double y) {
    hM=new TH1F("hM","",100,0,40);
    for (double x=xmin; x<xmax; x+=xstep) {
        TGraph *M=(TGraph*)elf->fieldLine(x,y,50)->Clone("M");
        M->Draw("L");
        cout << "Int(alpha*dx) = " <<elf->getAlphaSum() << endl;
        hM->Fill( elf->getAlphaSum() );
    }
}



TF2 *fV;
void doublemesh() {
    elf=new MaxCamElectricField;
    elf->setDriftPlane( 2.6*cm );  // plane
    elf->setWirePlane ( 0.3*mm, 1500*volt,  25.0*um,  0.25*mm,  240, -30*mm); // mesh
    elf->setWirePlane ( 0.0*mm, (1500+820)*volt,  7.50*um,  0.25*mm,  240, -30*mm+0.25*mm/2); // mesh
    elf->setFieldPlane(-0.3*mm, 1500*volt ); // plane

    fV = new TF2("fV",V,  -1.1*mm,1.1*mm,  -0.2*mm,0.5*mm,  1);
    fV->SetNpx(200);
    fV->SetNpy(200);
    fV->Draw("cont1");
    fV->GetHistogram()->SetXTitle("X (m)");
    fV->GetHistogram()->SetYTitle("Y (m)");
    
    plotField(-1*mm, 1.1*mm, 0.02*mm, 1.8*mm);
}


void mesh() {   
    elf=new MaxCamElectricField;
    elf->setDriftPlane( 1.6*cm ); // plane
    elf->setWirePlane ( 0.15*mm, 1500*volt,  25.0*um,  0.25*mm,  240, -30*mm); // mesh
    elf->setFieldPlane(-0.15*mm, (1500+900)*volt ); //plane

    fV = new TF2("fV",V,  -1.1*mm,1.1*mm,  -0.2*mm,0.5*mm,  1);
    fV->SetNpx(200);
    fV->SetNpy(200);
    fV->Draw("cont1");
    fV->GetHistogram()->SetXTitle("X (m)");
    fV->GetHistogram()->SetYTitle("Y (m)");
    
    plotField(-1*mm, 1.1*mm, 0.02*mm, 1.55*mm);
}



void wires() {   
    elf=new MaxCamElectricField;
    elf->setDriftPlane( 2.6*cm ); //plane
    elf->setWirePlane ( 0.3*cm, 1500*volt,  25.0*um,  0.2*cm,  80, -30*mm); // wires
    elf->setWirePlane ( 0.0*cm, 3700*volt,  7.50*um,  0.5*cm,  20, -30*mm); // wires
    elf->setFieldPlane(-0.3*cm, 1500*volt ); // plane

    fV = new TF2("fV",V,  -1.1*cm,1.1*cm,  -0.4*cm,2.0*cm,  1);
    fV->SetNpx(200);
    fV->SetNpy(200);
    fV->Draw("cont1");
    fV->GetHistogram()->SetXTitle("X (m)");
    fV->GetHistogram()->SetYTitle("Y (m)");
    
    plotField(-1*cm, 1.1*cm, 0.05*cm, 1.8*cm);
}



