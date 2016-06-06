#ifndef DMTPC_DETECTORPART_HH
#define DMTPC_DETECTORPART_HH

#include "TClonesArray.h"
#include "TObject.h"

#include "DmtpcLocation.hh"

#define DDP_NPAR 10

class DmtpcDetectorPart : public TObject {

public:
    
    enum DmtpcPartShape {
        kCylinder, 
        kCuboid, 
        kSphere,
        kOther
    };
    
    
    //
    // Ctors
    DmtpcDetectorPart();
    
    DmtpcDetectorPart(DmtpcDetectorPart::DmtpcPartShape geo);
    
    DmtpcDetectorPart(const DmtpcDetectorPart &other);
    
    
    virtual ~DmtpcDetectorPart();
    
    virtual const char* GetName() const { return "DmtpcDetectorPart"; }  
    
    void  setParameter(int i, float val) { if ( i>-1 && i<DDP_NPAR)  _param[i]=val; }
    float getParameter(int i) { return i>-1 && i<DDP_NPAR ? _param[i] : 0; }
    
    void setLocalCoordinate( DmtpcLocation loc ) { _localCoord=loc; }
    DmtpcLocation getLocalCoordinate() { return _localCoord; }
    
    void setComment(TString com) { _comment=com; }
    TString getComment() { return _comment; }
    
    void setMaterial(TString mat) { _material=mat; }
    TString getMaterial() { return _material; }

    void setShape(DmtpcDetectorPart::DmtpcPartShape sh) { _shape=sh; }
    DmtpcDetectorPart::DmtpcPartShape getShape() { return _shape; }

    int getModuleID() { return _moduleID; }
    void setModuleID(int id) { _moduleID=id; }

    int getPartID() { return _partID; }
    void setPartID(int id) { _partID=id; }

    TString getSN() { return _serial; }
    void setSN(TString ser) { _serial=ser; }
    
private:  
  
    DmtpcLocation _localCoord;
    
    float   _param[DDP_NPAR]; // parameters to describe the part
    
    TString _comment; // Comment about detector part, shape, purpose, etc.
    
    TString _material; // Description of material; e.g. SS104, Cu101, etc
    
    DmtpcPartShape _shape; // geometrical shape of object
    
    int _moduleID; // module identification for multi-module run
    
    int _partID; // part number inside a module
    
    TString _serial; // alphanumeric serial tag by manufacturer
    
    ClassDef(DmtpcDetectorPart,1)
        
};

#endif
