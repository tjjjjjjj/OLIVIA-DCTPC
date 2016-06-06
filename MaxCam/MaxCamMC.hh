#ifndef MAXCAM_MC_HH
#define MAXCAM_MC_HH

#include <assert.h>

#include "TString.h"
#include "MaxCamSRIM.hh"
#include "MaxCamDataset.hh"
#include "MaxCamElectricField.hh"

#include "TH2F.h"
#include "TRandom2.h"
#include "TGraph.h"
#include "TLorentzVector.h"
#include "TVector3.h"

class TRotation;

#include "math.h"
#include <vector>
using std::vector;


class MaxCamMC : public MaxCamSRIM, public MaxCamDataset {

public:
    
    // Ctors
    MaxCamMC(TString fname="");
    
    MaxCamMC(const MaxCamMC &other);
    
    virtual ~MaxCamMC();
    
    void setProjectile(double Ex, double Ey, double Ez, double m);
    TLorentzVector* getProjectile() { return _projectile; }
    void setIsotropicProjectile(double Ekin, double m);
    
    void setRecoil(double Ex, double Ey, double Ez, double M) {
        _recoil = new TLorentzVector(Ex,Ey,Ez, sqrt(Ex*Ex+Ey*Ey+Ez*Ez+M*M));
    }
    int makeRecoil(double E, double phi);
    int setRecoilEnergy(double E);
    int setRecoilEnergyAngle(double E, double phi);
    void setRecoilMass(double M);

    
    double calcCosRecoil(double E);
    
    TLorentzVector* getRecoil() { return _recoil; }
    
    bool propagateRecoil();
    
    void setPressure(double PressureInTorrs);
    
    void setRecoilCoord(double x, double y, double z) { 
        if (_recoilCoord) _recoilCoord->SetXYZ(x,y,z);
        else _recoilCoord=new TVector3(x,y,z); 
    }
    TVector3 *getRecoilCoord() { return _recoilCoord; }
    void setRandomRecoilCoord( double x1=0, double x2=0, double y1=0, double y2=0, double z1=0, double z2=0);
    
    void event(bool resetImage=true);
    
    void setStepSize(double step) { _stepSize=step; }
    double getStepSize() { return _stepSize; }

    void setStraggling(bool value) {_stragglingOn=value;}
    bool getStraggling() {return _stragglingOn;}  
    
    
    TH2F *getTrackImage() { return _trackImage; }
    void setTrackImage(int nx=400, float minx=-20, float maxx=20, int ny=400, float miny=-20, float maxy=20) {
        _trackImage=new TH2F("trackImage","",nx,minx,maxx, ny,miny,maxy);
    }
    
    TH2F *getWireImage() { return getAnodeImage(); }
    void setWireImage(int nx=330, float minx=-20, float maxx=20, int ny=330, float miny=-20, float maxy=20) {
        setAnodeImage(nx,minx,maxx, ny,miny,maxy);
    }
    TH2F *getAnodeImage() { return _wireImage; }
    void setAnodeImage(int nx=330, float minx=-20, float maxx=20, int ny=330, float miny=-20, float maxy=20) {
        _wireImage=new TH2F("wireImage","",nx,minx,maxx, ny,miny,maxy);
    }
    
    TH2F *getCCDImage() { return _ccdImage; }
    void setCCDImage(int nx=96, float minx=0, float maxx=768, int ny=64, float miny=0, float maxy=512) {
        _ccdImage=new TH2F("ccdImage","",nx,minx,maxx, ny,miny,maxy);
    }
    
    double applyDiffusion(double x);
    
    void applyBiasADC(TH2F* h);

      void applyRadialEffect(TH2F* h);
    
    double getNumPhotons() { return _totalPhotons; }

    double applyAvalancheWidth(double y);
    
    vector<double> getAnodeWireListX() const { return _anodeWireListX; }
    vector<double> getAnodeWireListY() const { return _anodeWireListY; }

    vector<double> getSpacerListX() const { return _spacerListX; }
    vector<double> getSpacerListY() const { return _spacerListY; }
    
    void addWire(double r, TString which="y") {
        if (which=="x") _anodeWireListX.push_back(r);
        else            _anodeWireListY.push_back(r);
    }
    void clearWireList(TString which="y") {
        if (which=="x") if (getAnodeWireListX().size()) getAnodeWireListX().clear();
        else            if (getAnodeWireListY().size()) getAnodeWireListY().clear();
    }

    void addSpacer(double r, TString which="y") {
        if (which=="x") _spacerListX.push_back(r);
        else            _spacerListY.push_back(r);
    }
    void clearSpacerList(TString which="y") {
        if (which=="x") getSpacerListX().clear();
        else            getSpacerListY().clear();
    }

    double getSpacerDiameter() const {return _spacerDiameter;}
    void setSpacerDiameter( double diam) { _spacerDiameter = diam;}

    void applySpacers(TH2* image);

    virtual TString GetName() { return TString("MaxCamMC"); }
    
    TGraph *getRangeVsEnergyProject(TString opt="long");
    TGraph *getEnergyVsRangeProject(TString opt="long");
    
    void setNoiseADC(double noise) { _calibration._noiseADC=noise; }
    double getNoiseADC() { return _calibration._noiseADC; }
    
    void setElecScintillation(double f) { _calibration._f_elec=f; }
    double getElecScintillation() { return _calibration._f_elec; }
    
    void setNuclScintillation(double f) { _calibration._f_nucl=f; }
    double getNuclScintillation() { return _calibration._f_nucl; }
    
    void setPixelsPermm(double f) { _calibration._pixelsPermm=f; }
    double getPixelsPermm() { return _calibration._pixelsPermm; }
    
    void setPhotonsPerkeV(double val) { _calibration._photonsPerkeV=val; }
    double getPhotonsPerkeV() { return  _calibration._photonsPerkeV; }

    void setAvalancheWidth(double val) { _calibration._avalancheWidth=val; }
    double getAvalancheWidth() { return _calibration._avalancheWidth; }

    void setDiffusionConstTerm(double val) { _calibration._diffusionConstTerm=val; }
    double getDiffusionConstTerm() { return _calibration._diffusionConstTerm; }

    void setDiffusionDzTerm(double val) { _calibration._diffusionDzTerm=val; }
    double getDiffusionDzTerm() { return _calibration._diffusionDzTerm; }


    
private:
    TRandom2 *_rnd;
    
    TLorentzVector *_projectile;
    TLorentzVector *_recoil;
    TVector3 *_recoilCoord;
    TRotation *_rotation;
    
    TH2F *_trackImage;
    TH2F *_wireImage;
    TH2F *_ccdImage;
    
    double _minEnergy; // energy cutoff in simulation
    double _stepSize;  // step size in SRIM simulation
    double _totalPhotons; // count of photons per track

    bool _stragglingOn;
    
    struct Calibration {
        double _photonsPerkeV; // # of photons per energy (from alpha calibration)
        double _pixelsPermm;
        
        double _f_elec; // fraction of electronic energy losss that scintillates
        double _f_nucl; // fraction of nuclear energy loss that scintillates
        
        double _noiseADC; // ADC noise

        double _avalancheWidth; // width of scintillation perpendicular to wires
        double _diffusionConstTerm; // constant diffusion term 
        double _diffusionDzTerm; // diffusion term that depends on drift length

        
    } _calibration;
        
    int findClosestWire(double r, vector<double> &wireList) {
        double min=1e10;
        int    iw=-1;
        if (!wireList.size()) return -1;
        for (unsigned int i=0; i<wireList.size(); i++) {
            if (fabs(wireList[i]-r)<min) { 
                min=fabs(wireList[i]-r); iw=i; 
            }
        }
        return iw;
    }

    vector<double> _spacerListX;
    vector<double> _spacerListY;

    double _spacerDiameter;
    
    vector<double> _anodeWireListX;
    vector<double> _anodeWireListY;
    bool findAnodeWire(double x, double y, double &xw, double &yw) {
        int ix=findClosestWire(x, _anodeWireListX ); 
        int iy=findClosestWire(y, _anodeWireListY );
        assert (ix>-1 || iy>-1);
        bool isX=true;
        if (ix>-1 && iy<0)  isX=true; // only x-wires
        else if (ix<0  && iy>-1) isX=false; // only y-wires
        else isX = fabs(getAnodeWireListX()[ix]-x)<fabs(getAnodeWireListY()[iy]<y) ? true : false; // both x&y
        if (isX) {
            xw=getAnodeWireListX()[ix]; yw=y;
        }
        else {
            xw=x; yw=getAnodeWireListY()[iy];
        }
        return isX;
    }

    bool findAnodeWire(double &x, double &y);

    
    MaxCamElectricField *_elf;
    
    void driftElectrons(double npho);

    ClassDef(MaxCamMC,0)

        };

#endif

