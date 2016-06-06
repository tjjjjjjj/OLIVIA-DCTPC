#ifndef MAXCAM_ENDF_HH
#define MAXCAM_ENDF_HH

#include "TString.h"

class TGraph;
class TGraph2D;
class TF1;

#include "math.h"

class MaxCamENDF  {

public:

    // Ctors
    MaxCamENDF( const char* fileName="ENDF_DCS_n_on_19F", TString what="scattering" );
    
    MaxCamENDF(const MaxCamENDF &other);
    
    virtual ~MaxCamENDF();
    
    void   fillScatteringTable( const char* fileName );

    void fillCSTable(const char* fileName);
    
    void   fillFissionTable( const char* fileName );
    
    virtual TString GetName() { return TString("MaxCamENDF"); }
    
    TGraph2D* getDifferentialCrossSection() { return _DCS; }

    double generateCosAngleCMS(double energy);
    
    TGraph* getCrossSection() { return _CS; }
    bool acceptEnergy(double E);

    
    double generateEnergy(double emin=0.1, double emax=10000);

    TF1 *legendre() { return _fDCS; }

    double maxCS() { return _maxCS; }

    static double generateRockNeutronEnergy(TString location="generic");
    static double generateCosmicNeutronEnergy();
    
private:

    TF1       *_fDCS;
    TGraph2D  *_DCS; 
    TGraph    *_CS;
    
    double _maxCS;
    
  ClassDef(MaxCamENDF,0)
};

#endif

