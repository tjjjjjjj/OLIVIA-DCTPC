#include "McDarkElectricFieldMWPC.hh"
#include "G4PhysicalConstants.hh"
#include "globals.hh"

#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMatrixDLazy.h"
#include "TString.h"
#include "TF2.h"
#include "TGraph.h"


// ------ voltages:

double Vc = 1500*volt;
double Va = 2500*volt + Vc;

// ------ geometry:

double dr = 26*mm;
double h1 = 3*mm;
double h2 = 3*mm;
double sc = 2*mm;
double sa = 5*mm;
double rc = 50*um/2;
double ra = 15*um/2;

int nc = 50;
int na = 20;
double xa = -50*mm;
double xc = -50*mm;

// ------- derived:

double d = h1-h2+dr;
double D = h1+h2+dr;
double a = pi/D;
double VR = Vc/Va;
TVectorD *Qa=0;
TVectorD *Qc=0;


McDarkElectricFieldMWPC::McDarkElectricFieldMWPC() {
    Qa = charge("a");
    Qc = charge("c");
}


McDarkElectricFieldMWPC::~McDarkElectricFieldMWPC() {}


McDarkElectricFieldMWPC::McDarkElectricFieldMWPC(const McDarkElectricFieldMWPC &p)
   : G4ElectricField(p) {
    G4cout << "NOT implemented" << G4endl;
}


McDarkElectricFieldMWPC&
McDarkElectricFieldMWPC::operator = (const McDarkElectricFieldMWPC &p) {
    G4cout << "NOT implemented" << G4endl;
    return *this;
}



double
McDarkElectricFieldMWPC::calcL(double x, double y, double xn, double yn, double r) const {

    if ( (x-xn)*(x-xn) + (y-yn)*(y-yn) < r*r )
         return  log( a * r * 0.5 / cos(a*(yn-d*0.5)) );
    
    double A = cosh(a*(x-xn)) - cos(a*(y-yn));
    double B = cosh(a*(x-xn)) + cos(a*(y+yn-d));
    return 0.5*log(A/B);
}


TVectorD*
McDarkElectricFieldMWPC::charge(TString opt) {
    
    bool isAnode=opt=="a" ? true : false;

    int n = isAnode ? na : nc;
    double s = isAnode ? sa : sc;
    double r = isAnode ? ra : rc;
    double y = isAnode ? 0.0: h1;
    
    // inverse capacitance matrix
    TMatrixD L = THilbertMatrixT<double>( n, n);
    for (int i=0; i<n; i++) {
        double xi = s*i;
        for (int j=0; j<n; j++) {
            double xj =  s*j;
            L(i,j) = calcL(xi, y, xj, y, r);
        }
    }
    TMatrixD L0 = L;
    double det=0;
    L.Invert(&det);
    //L0*=L;
    //L0.Print(); // check inversion
    
    // voltage vector
    double Pj = isAnode ? 1 : VR;
    double PC = isAnode ? VR*(h1+dr)/D : VR*dr/D;
    TVectorT<double> P(n);
    for (int i=0; i<n; i++) P(i) = PC-Pj;
    
    // compute charge
    P *= L;
    
    return new TVectorD(P);
}



double
McDarkElectricFieldMWPC::V(double x, double y) const {

    double PW=0;
    for (int i=0; i<na; i++) PW  -= (*Qa)[i] * calcL(x,y, xa+i*sa,  0.0, ra);

    double PC1=0;
    for (int i=0; i<nc; i++) PC1 -= (*Qc)[i] * calcL(x,y, xc+i*sc, h1, rc);


    double PC2=VR*(h1+dr-y)/D;

    return Va*(PW - PC1 + PC2);
}



void
McDarkElectricFieldMWPC::GetFieldValue( const G4double Point[4], G4double *field ) const {

    double x=Point[0];
    double y=Point[1];
    double dx=1*um;
    double dy=1*um;

    
    field[0]=field[1]=field[2]=0;

    double Vxy = V(x,y);
    field[3] = (V(x+dx,y)-Vxy) / dx;
    field[4] = (V(x,y+dy)-Vxy) / dy;
    field[5] = 0;


}
