#ifndef MAXCAM_SRIM_HH
#define MAXCAM_SRIM_HH
#include "DmtpcMath.hh" 

#include <assert.h>

#include "TString.h"
#include "TGraph.h"
#include "TNtuple.h"
#include "TH2.h"
#include "math.h"

class MaxCamSRIM  {

public:

    // Ctors
    MaxCamSRIM( const char* fileName="SRIM_F_in_CF4_100Torr", TString energyUnitString="keV" );
    
    MaxCamSRIM(const MaxCamSRIM &other);

    virtual ~MaxCamSRIM();
    
    void setStopping(float pressure);

    void setPressure(float pressure) { setStopping(pressure); }
    
    TGraph *readSRIM(const char* fileName, TString what, float minE, float maxE, float dEdxFactor);
    
    void   fillSrimTable( const char* fileName );
    
    TGraph *getStoppingVsEnergy(bool isTotal=true) { return isTotal ? _stopping : _stoppingElec; }
  
    TGraph *getStoppingVsRange(bool isTotal=true) { return isTotal ? _range : _rangeElec; }
    
    TGraph *getRangeVsEnergy(bool isTotal=true) { return isTotal ? _rangenergy : _rangenergyElec; }
    
    TGraph *getEnergyVsRange(bool isTotal=true) { return isTotal ? _enerange : _enerangeElec; }

    TGraph *getLStraggleVsEnergy() { return _lstraggle; }

    TGraph *getTStraggleVsEnergy() { return _tstraggle; }

    TGraph *getRangeStep() { return _stepSize; }
    
    static double density(double pTorr, double amu=88, double T=300);

    static double numberDensity(double pTorr, double T=300);
    
        double getdedx(double energy);
    
    double calcEnergyLoss(double ionEnergy, double distance0=0, double distance1=1e10, bool isTotal=true, double stepsize_in_um = 50);
    
    virtual TString GetName() { return TString("MaxCamSRIM"); }
    
    void Print();

    TNtuple *getTable() { return _srimTable; }

    float getCurrentPressure() { return _press; }

    static double X0(int A, int Z);

    float getDefaultDensity() { return _targetDensityAt600Torr; }


    // added 4/10/2010 by Shawn; creates TGraphs of the raw numbers in the SRIM file
    TGraph *getRawTotalStoppingVsEnergy() { return _rawtotalstopping; }
    TGraph *getRawElectricStoppingVsEnergy() { return _rawelectricstopping; }
    TGraph *getRawNuclearStoppingVsEnergy() { return _rawnuclearstopping; }

    TGraph *getRawRangeVsEnergy() { return _rawrangenergy; }

    void SetEnergyUnits(TString input_string);
    
    TString GetEnergyUnits(void){ return _energyUnitString; }

    float GetEnergyConversionFactor(void);

    Double_t ROOTv24TGraphEval(TGraph* gr,Double_t x);


private:
    
    TFile *f;
    TH2D *hist;
    TH1D *proj;
    
    TNtuple *_srimTable;
    TGraph  *_stopping; // dE/dx = dE/dx(E) (total=elec+nucl)
    TGraph  *_range;    // dE/dx = dE/dx(x) (total=elec+nucl)
    TGraph  *_rangenergy; // x = x(E)  (total=elec+nucl)
    TGraph  *_enerange; // E = E(x)  (total=elec+nucl)
    TGraph  *_lstraggle; // sigmaR = sigmaR(E)
    TGraph  *_tstraggle; // sigmaT = sigmaT(E)

    TGraph  *_stoppingElec; // dE/dx = dE/dx(E)  (elec)
    TGraph  *_rangeElec;    // dE/dx = dE/dx(x)  (elec)
    TGraph  *_rangenergyElec; // x = x(E)   (elec)
    TGraph  *_enerangeElec; // E = E(x)  (elec)
    TGraph  *_stepSize; 

    // added 4/10/2010 by Shawn; TGraphs of the raw numbers in the SRIM file
    TGraph *_rawtotalstopping; // dE/dx = dE/dx(E) (raw SRIM file total=elec+nucl)
    TGraph *_rawelectricstopping; // dE/dx = dE/dx(E) (raw SRIM file elec)
    TGraph *_rawnuclearstopping; // dE/dx = dE/dx(E) (raw SRIM file nucl)

    TGraph *_rawrangenergy; // x = x(E)  (total expected range, raw in SRIM file)

    // added 4/10/2010 by Shawn to allow all energies to be in different units that keV
    // although keV is still default
    TString _energyUnitString;
  
    float _targetDensityAt600Torr;
    float _projectileMass;
    float _press;
    double _dedx;
    
    ClassDef(MaxCamSRIM,0)
};

#endif

