#include "DmtpcDetectorPart.hh"

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;




ClassImp(DmtpcDetectorPart)

//____________________
//
// Class describing a detector part. Relevant detector parts with
// brief description is used to store experimental setup. 
// 
//
//
DmtpcDetectorPart::DmtpcDetectorPart()  {
  
}


DmtpcDetectorPart::DmtpcDetectorPart(DmtpcDetectorPart::DmtpcPartShape shape) :
  _shape(shape)  {
}


DmtpcDetectorPart::DmtpcDetectorPart(const DmtpcDetectorPart &other) : TObject() { 
  
  _localCoord=other._localCoord;
  _comment=other._comment;
  _material=other._material;
  _shape=other._shape;
  _moduleID=other._moduleID;
  _partID=other._partID;
  _serial=other._serial;
  
  for (int i=0; i<DDP_NPAR; i++) _param[i]=other._param[i];
}

DmtpcDetectorPart::~DmtpcDetectorPart() {
}

