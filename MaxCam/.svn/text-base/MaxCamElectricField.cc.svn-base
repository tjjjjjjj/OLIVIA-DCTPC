#include "MaxCamElectricField.hh"

#include "TMatrixD.h"
#include "TMatrixDLazy.h"
#include "TF2.h"
#include "TGraph.h"

#include <iostream>
using std::cout;
using std::endl;

#include <vector>
using std::vector;

using MaxCamUnits::pi;
using MaxCamUnits::um;
using MaxCamUnits::mm;
using MaxCamUnits::cm;






MaxCamElectricField::MaxCamElectricField()  {
    //
    //   MWPC with drift region (from Mathieson, 4.3.5):
    //
    //     _____________________________________________ ground
    //
    //
    // dr 
    //      
    //     o          o          o          o          o Vc
    //
    //
    // h1
    //
    //     o     o     o     o     o     o     o     o   Va
    //
    // h2
    //     _____________________________________________ 
    //                                                   Vc
    //
    
    
    _Vc = 1500;
    _Va = 3700;
    _Vp = _Vc;
    
// ------ geometry:
    
    _dr = 26*mm;
    _h1 = 5*mm;
    _h2 = 5*mm;

    
// ------- derived:

    _Qa=0;
    resetCharges();

    _fieldLine=0;

// ------- Townsend & attachemnt (from Garfield/Magboltz)

    double press2N=101325./760/(1.38e-23*300)*1e-6;

    // christophorou
    double ep[]={8, 9, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 200, 250, 300, 350, 400, 450, 500, 600, 620};
    double ap[]={-0.052, -0.056, -0.058, -0.051, -0.043, -0.143, -0.281, -0.438, -0.730, -1.09, -1.47, -2.24, -2.92, -3.47, -3.73,
                 -3.66, -3.43, -2.68, -1.54, -0.033, 1.61, 12.0, 23.5, 35.0, 47.2, 59.4, 71.4, 83.6, 107.4, 111.7};
    int nN=sizeof(ep)/sizeof(double);
    double att[30];
    for (int i=0; i<nN; i++) {
        ep[i]*=1e-17*press2N;
        ap[i]*=1e-18*press2N;
        att[i]=0; // ap is effective;
    }

    
    // 50Torr
    /*double ep[]={2, 2.8769, 4.13828, 5.9527, 8.56266, 12.317, 17.7173, 25.4855, 36.6596, 52.733,
                 75.8538, 109.112, 156.952, 225.768, 324.755, 467.144, 671.964, 966.586, 1390.39, 2000};

    double ap[]={0, 0, 0, 0, 0, 0, 0, 0.0010662, 0.038151, 0.24281, 0.7088, 1.4702, 2.5809,
                 4.09658, 6.14396, 8.70694, 11.873, 15.4921, 19.4514, 23.4873};

    double att[]={0, 0, 0, 0, 0.01674, 0.50717, 2.85437, 6.07039, 7.79018, 6.45795, 4.25768,
    2.74698, 1.73707, 1.10161, 0.67724, 0.41392, 0.2565, 0.15875, 0.095986, 0.053093};*/

    // 75 Torr
    /*double ep[]={1.33333, 1.91793, 2.75885, 3.96847, 5.70844, 8.21131, 11.8116, 16.9903, 24.4397, 35.1553,
                  50.5692, 72.7413, 104.635, 150.512, 216.504, 311.43, 447.976, 644.391, 926.924, 1333.33};

    double ap[]={0, 0, 0, 0, 0, 0, 0, 0, 0.00061402, 0.028361, 0.20817, 0.6404, 1.35982, 2.42219,
                 3.9046, 5.85874, 8.3866, 11.4566, 15.0599, 19.0303};

    double att[]={0, 0, 0, 0, 0, 0.013145, 0.60384, 3.68083, 8.56442, 11.5585,
    10.0296, 6.74724, 4.3172, 2.75925, 1.70817, 1.05968, 0.66746, 0.40269, 0.25343, 0.14685};*/
    
    // 100 Torr
    /*double ep[]={1, 1.43845, 2.06914, 2.97635, 4.28133, 6.15848, 8.85867, 12.7428, 18.3298, 26.3665,
                 37.9269, 54.5559, 78.476, 112.884, 162.378, 233.572, 335.982, 483.293, 695.193, 1000};

    double ap[]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0.001864, 0.046279, 0.27022, 0.76001, 1.55887,
                 2.70397, 4.27503, 6.3567, 9.00146, 12.1739, 15.829};

    double att[]={0, 0, 0, 0, 0, 0, 0.033837, 1.29909, 6.26403, 13.0678, 15.4957, 11.9844,
    8.29852, 5.20542, 3.29469, 2.08159, 1.28679, 0.79609, 0.50365, 0.2871};*/
    

    
    // 200Torr
    /*
    double ep[]={0.5, 0.71922, 1.03457, 1.48818, 2.14067, 3.07924, 4.42933, 6.37138, 9.1649, 13.1833,
                 18.9635, 27.278, 39.238, 56.4419, 81.1888, 116.786, 167.991, 241.647, 347.596, 500};

    double ap[]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.000015633, 0.0024984, 0.062939, 0.30352,
                 0.82138, 1.65068, 2.82334, 4.45546, 6.57182, 9.25806};
    
    double att[]={0, 0, 0, 0, 0, 0, 0, 0, 0.15502, 3.30997, 13.6512, 27.2128, 32.9839,
                  23.8228, 15.8653, 10.2734, 6.25024, 3.91459, 2.47421, 1.53363};
    
    */            
    _ap_vs_ep = new TGraph(30, ep, ap);
    _at_vs_ep = new TGraph(30, ep, att);
}





void
MaxCamElectricField::resetCharges() {
    _d = _h1-_h2+_dr;
    _D = _h1+_h2+_dr;
    _a = pi/_D;
    _VR = _Vc/_Va;
    _VP = _Vp/_Va;
    if (_Qa) delete _Qa;
    _Qa=0;
}



double
MaxCamElectricField::calcL(double x, double y, double xn, double yn, double r) {

    if ( (x-xn)*(x-xn) + (y-yn)*(y-yn) < r*r ) {
        _statFlag |= 1;
        return  log( _a * r * 0.5 / cos(_a*(yn-_d*0.5)) );
    }
    
    double A = cosh(_a*(x-xn)) - cos(_a*(y-yn));
    double B = cosh(_a*(x-xn)) + cos(_a*(y+yn-_d));
    return 0.5*log(A/B);
}



TVectorD*
MaxCamElectricField::calcCharge() {
    cout << GetName() << ": Compute charges for "
         << wireList.size()
         << " wires" << endl;

    int n=wireList.size();
    
    TMatrixD L = THilbertMatrixT<double>( n, n);
    for (int i=0; i<n; i++) {
        double xi = wireList[i].x;
        double yi = wireList[i].y;
        double ri = wireList[i].radius;
        for (int j=0; j<n; j++) {
            double xj = wireList[j].x;
            double yj = wireList[j].y;
            L(i,j) = calcL(xi, yi, xj, yj, ri);
        }
    }
    TMatrixD L0 = L;
    double det=0;
    L.Invert(&det);
    //L0*=L;
    //L0.Print(); // check inversion
    
    // voltage vector
    TVectorT<double> P(n);
    for (int i=0; i<n; i++) {
        P(i) = wireList[i].V - _Vp*(_h1+_dr - wireList[i].y)/_D;
    }
    
    // compute charge
    P *= L;
    
    return new TVectorD(P);    
}


double
MaxCamElectricField::V(double *var) {
    
    double x=var[0];
    double y=var[1];

    if (!_Qa) { _Qa=calcCharge(); _statFlag=0; } // reset flag about being too close to wire
    double PW=0;
    for (unsigned int i=0; i<wireList.size(); i++) PW -= (*_Qa)[i] * calcL(x,y, wireList[i].x,  wireList[i].y,  wireList[i].radius);

    double PC2 = _Vp*(_h1+_dr-y)/_D;

    double V = (-PW + PC2);
    
    if (y<-_h2) { V=_Vp; _statFlag|=1; }
    else if (y>_h1+_dr) { V=0; _statFlag|=1; }

    
    return V;
}


double
MaxCamElectricField::V(double x, double y) {
    double var[]={x,y};
    return V(var);
}


double
MaxCamElectricField::Ex(double x, double y) {
    double dx=0.1*um;
    double var[]={x,y};
    double dvar[]={x+dx,y};
    return (V(var)-V(dvar))/dx;
}

double
MaxCamElectricField::Ey(double x, double y) {
    double dy=0.1*um;
    double var[]={x,y};
    double dvar[]={x,y+dy};
    return (V(var)-V(dvar))/dy;
}

void
MaxCamElectricField::Print() {
    cout << _dr+_h1 << "  " << 0 << endl;
    for (unsigned int i=0; i<wireList.size(); i++) {
        cout << wireList[i].x <<", "<< wireList[i].y <<";   "
             << wireList[i].radius <<", "<< wireList[i].V <<", "<< endl;
    }
    cout << -_h2    << "  " << _Vp << endl;
}

TGraph*
MaxCamElectricField::fieldLine(double x, double y, double pressure, double step, int charge) {

    if (_fieldLine) delete _fieldLine;
    vector<double> xx,yy;
    _electricFieldSum = _alphaSum = 0;
    _statFlag = 0;
    double ex,ey,el,al;
    while (!_statFlag) {
        ex = Ex(x,y)*charge;
        ey = Ey(x,y)*charge;
        el = sqrt(ex*ex+ey*ey);
        al = alpha( el, pressure);
        _electricFieldSum+=el;
        _alphaSum+=al;
        x += ex/el*step;
        y += ey/el*step;
        if (!_statFlag) {
            xx.push_back(x);
            yy.push_back(y);
            //cout << "E("<<x<<","<<y<<")="<<ex<<","<<ey<<"   " << step << "  " << al << endl;
        }
    }
    _electricFieldSum *= step;
    _alphaSum *= step;
    _fieldLine = new TGraph(xx.size(),&xx[0],&yy[0]);
    return _fieldLine;
}


void
MaxCamElectricField::driftStep(double &x, double &y, double step) {
    double ex= -Ex(x,y);
    double ey= -Ey(x,y);
    double el = sqrt(ex*ex+ey*ey);
    
    x += ex/el*step;
    y += ey/el*step;
    //cout << "E("<<x<<","<<y<<")="<<ex/el<<","<<ey/el<<"   " << step << endl;
}



double
MaxCamElectricField::alpha(double E, double P) {
    assert(P>0);
    E*=cm; // V/m -> V/cm
    double alpha_over_p = _ap_vs_ep->Eval( E/P );
    double attachment = _at_vs_ep->Eval( E/P );

    //cout << E << "  " << alpha_over_p*P << "  " << attachment << endl;
    
    return (alpha_over_p*P - attachment )/cm  ; // 1/cm -> 1/m
}


double
MaxCamElectricField::attachment(double E, double P) {
    assert(P>0);
    E*=cm; // V/m -> V/cm
    double attachment = _at_vs_ep->Eval( E/P );
    
    return attachment/cm  ; // 1/cm -> 1/m
}
