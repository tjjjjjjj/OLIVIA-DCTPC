//-------------------------------------------------------------------
//
// Created by: Timur Sahin (tcsahin@MIT.EDU), Denis Dujmic (ddujmic@mit.edu)
// Date:       Oct 15, 2007
// Copyright:  MIT 2007
//
//
//
// $Id: McDarkElectricFieldMWPC.hh,v 1.2 2007/11/03 14:35:14 ddujmic Exp $

#ifndef McDarkElectricFieldMWPC_H
#define McDarkElectricFieldMWPC_H 1


#include "G4Types.hh"
#include "G4ElectricField.hh"

#include "TVectorD.h"
#include "TString.h"

class McDarkElectricFieldMWPC : public G4ElectricField {
    
public:  // with description
    
    McDarkElectricFieldMWPC();
    virtual ~McDarkElectricFieldMWPC();
    
    McDarkElectricFieldMWPC(const McDarkElectricFieldMWPC &r);
    McDarkElectricFieldMWPC& operator = (const McDarkElectricFieldMWPC &p);
    
    virtual void  GetFieldValue( const G4double Point[4],
                                 G4double *field ) const;

    void setAnodePlane(double &x0, double &y0, double &diameter, double &xpitch) {
        _anodeX0=&x0;
        _anodeY0=&y0;
        _anodeDiameter=&diameter;
        _anodePitch=&xpitch;
    }

    void setAnodeVoltage(double &v) { _anodeVoltage=&v; }
    
    void setGroundPlane(double &uppery, double &lowery) {
        _upperGroundY = &uppery;
        _lowerGroundY = &lowery;
    }

    void setCathodePlane(double &y) {
        _cathodeY = &y;
    }

    void setCathodeVoltage(double &v) { _cathodeVoltage=&v; }
    
private:

    double V(double x, double y) const;

    double calcL(double x, double y, double xn, double yn, double r) const;

    TVectorD* charge(TString opt="a");
    
    double *_anodeVoltage;
    double *_anodeX0;
    double *_anodeY0;
    double *_anodeDiameter;
    double *_anodePitch;

    double *_upperGroundY;
    double *_lowerGroundY;

    double *_cathodeVoltage;
    double *_cathodeY;
};


#endif
