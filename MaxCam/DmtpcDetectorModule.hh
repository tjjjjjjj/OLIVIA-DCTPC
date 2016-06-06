#ifndef DMTPC_DETECTORMODULE_HH
#define DMTPC_DETECTORMODULE_HH

#include "TClonesArray.h"
#include "TObject.h"

//class DmtpcDetectorLocation;



class DmtpcDetectorModule : public TObject {

public:
  //
  // Ctors
  DmtpcDetectorModule();
  
  DmtpcDetectorModule(const DmtpcDetectorModule &other);
  
  virtual ~DmtpcDetectorModule();
  
  virtual const char* GetName() const { return "DmtpcDetectorModule"; }  
  
  TClonesArray*  listOfParts() { return _listOfParts; }
  
  int DetectorModuleSerialNumber() { return _DetectorModuleSerialNumber; }
  void setDetectorModuleSerialNumber(int iev) { _DetectorModuleSerialNumber=iev; }

  //DmtpcLocation* DetectorModuleLocation() { return _DetectorModuleLocation; }
  //void setDetectorModuleLocation(DmtpcLocation *loc) { _DetectorModuleLocation=loc; }

private:    
    
  TClonesArray *_listOfParts; // list of DetectorModule parts

  int _DetectorModuleSerialNumber;
    
  //DmtpcLocation *_DetectorModuleLocation;

  ClassDef(DmtpcDetectorModule,1)
        
};

#endif
