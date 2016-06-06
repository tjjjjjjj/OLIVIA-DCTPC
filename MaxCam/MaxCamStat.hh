#ifndef MAXCAM_STAT_HH
#define MAXCAM_STAT_HH

#include <assert.h>

class TH1F;
class TGraphAsymmErrors;

namespace MaxCamStat {

    enum ValueType {	Fraction, Asymmetry    };
    enum ErrorType {	Poisson, Binomial    };
    
    void calculateError(int Nneg, int Ntot, double &eplus, double &eminus, ValueType type);
    
    TGraphAsymmErrors* fractionError(TH1F* h1, TH1F* h2);
    TGraphAsymmErrors* asymmetryError(TH1F* h1, TH1F* h2);

};

#endif

