#ifndef McDarkCamera_H
#define McDarkCamera_H 1

#include "globals.hh"

#include "TH2.h"
#include "MaxCamConfig.hh"
#include "G4ThreeVector.hh"

class G4PVPlacement;

class McDarkCamera {

public:

    McDarkCamera(G4String camName, G4String chipName,
                 G4double bias, G4double noise,
                 G4int binX, G4int binY );

    ~McDarkCamera();

    void processHit(G4ThreeVector &hitPos, G4double hitPhot, G4int &chanID,  G4double &weight);
    // Turn hits in sensitive volume to CCD pixels.
 
    void processWorm(G4ThreeVector &hitPos, G4double quenchedEnergy, G4int &chanID,  G4double &weight);
    // Turn worms in sensitive volume to CCD pixels.
    

    void setViewfield(G4ThreeVector posInFocalPlane,
                      G4ThreeVector xaxis,
                      G4ThreeVector yaxis,
                      G4double avgAcceptance);
    void setAngularPar(G4double par) { _angularPar=par; }
    // Sets viewfield parameters: origin and axes of the viewfield,
    // average acceptance. Gaussian decrease of acceptance with radius.
    // (default sigma = maximum radius)
    
    
    TH2F* image() { return _image; }

    TH2S* rawImage();

    MaxCamConfig* config() { return _config; }

    void clearImage() { _image->Reset(); }

    void addNoise();
    // Adds readout noise to current image.
    // The mean value, RMS are taken from MaxCamConfig::meanForExposedPixels
    // and MaxCamConfig::rmsForExposedPixels, respectively.

    void setActiveVolume(G4PVPlacement* pl) { _activeVolume=pl; }
    // Set active volume for worms.
    
    G4PVPlacement* getActiveVolume() { return _activeVolume; }
    // Get active volume


protected:
   

    G4double calcAcceptance(G4double x, G4double y);
    //
    // Total optical efficiency:
    // _maxAcceptance = total efficiency for center of viewfield 
    // Gaussian falloff with radius with sigma2=_angularPar
    // i.e.
    //    acceptance = _maxAcceptance * exp( -0.5*r^2/_angularPar )
    //
    //  
  
private:


    TH2F *_image;
    // Pointer to image histogram
    
    MaxCamConfig *_config;
    // Pointer to CCD configuration

    G4ThreeVector _origin, _xaxis, _yaxis;
    // viewfield in TPC coordinates
    // origin is the corner of CCD viewfield, and
    // xaxis and yaxis are 2-dimensional vectors (z=0).

    G4double _maxAcceptance;
    // Maximum photon acceptance (center of lens);
    // includes all possible losses (mesh, window, lens, ccd QE)

    G4double _angularPar;
    // Relative decrease of lens efficiency for larger viewing angles.
    // It's 1.0 at the center of the lens with Gaussian decrease for
    // larger angles. The sigma of the Gaussian is _angularPar.
    // i.e. for maximum angle r=sqrt(x*x+y*y)/2, _angularPar=r
    // so  efficincy = exp(-0.5/_angularPar)

    G4PVPlacement *_activeVolume;
    // Active volume for worm creation 

};
    
#endif
