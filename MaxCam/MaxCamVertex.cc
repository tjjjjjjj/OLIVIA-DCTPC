#include "MaxCamVertex.hh"

#include "MaxCamImageTools.hh"

#include "TF1.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

#include <math.h>
#include <iostream>
using std::cout;
using std::endl;

    MaxCamVertex::MaxCamVertex(){;}



MaxCamVertex::MaxCamVertex(const MaxCamVertex &other) {
    // Not yet done
}



int
MaxCamVertex::fit2D() {


    //TMinuit
    return 0;
}



double
MaxCamVertex::calcChiSquare2D() {

    double chi2=0;
    unsigned int n=_trackList2D.size();
    for (unsigned int i=0; i<n; i++) {
        chi2 += distanceSquare2D( _xfit, _yfit, _trackList2D[i]);
    }
    return chi2;
}


double
MaxCamVertex::distanceSquare2D(double x, double y, TF1* f) {

    double a0 = f->GetParameter(0);
    double a1 = f->GetParameter(1);
//     cout << "a0, a1     :" << a0 << ", " << a1 << endl;
    
    double b1 = -1/a1;
    double b0 = y-b1*x;
//     cout << "b0, b1     :" << b0 << ", " << b1 << endl;

    double yint = (a0*b1-a1*b0)/(b1-a1);
    double xint = (yint-b0)/b1;
//     cout << "xint, yint :" << xint << ", " << yint << endl;

    return (x-xint)*(x-xint) + (y-yint)*(y-yint);
}
