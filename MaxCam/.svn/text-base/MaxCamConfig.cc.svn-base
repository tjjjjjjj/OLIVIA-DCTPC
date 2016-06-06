#include "MaxCamConfig.hh"
//#include "MaxCam.hh"

#include <iostream>
#include <fstream>
#include <vector>
using std::cout;
using std::endl;
using std::cerr;
using std::vector;


ClassImp(MaxCamConfig)

//____________________
// Class contains camera setup for output. 
//
//
MaxCamConfig::MaxCamConfig(const char *name, const char *title) : 
  TNamed(name, title), 
  cameraID(-1),
  row_width(-1), img_rows(-1),
  hbin(1), vbin(1),
  ul_x(-1), ul_y(-1), lr_x(-1), lr_y(-1),
  CCDTemp(666), CCDTempSet(-20.0), exposureTime(-1),
  frameType(normal),
  cleanType(none),
  nFlushes(-1),
  bitDepth(65535),
  daqTime(-1),
  serialNumber(""),
  overscanColumns(0),
  overscanRows(0),
  digitizeOverscan(false)
{}


MaxCamConfig::MaxCamConfig(const MaxCamConfig &other) :
  TNamed(other), 
  cameraID(other.cameraID),
  row_width(other.row_width), img_rows(other.img_rows),
  hbin(other.hbin), vbin(other.vbin),
  ul_x(other.ul_x), ul_y(other.ul_y), lr_x(other.lr_x), lr_y(other.lr_y),
  CCDTemp(other.CCDTemp), CCDTempSet(other.CCDTempSet), 
  exposureTime(other.exposureTime),
  frameType(other.frameType),
  cleanType(other.cleanType),
  nFlushes(other.nFlushes),
  bitDepth(other.bitDepth),
  daqTime(other.daqTime),
  serialNumber(other.serialNumber),
  overscanColumns(other.overscanColumns),
  overscanRows(other.overscanRows),
  digitizeOverscan(other.digitizeOverscan),
  meanForExposedPixels(other.meanForExposedPixels),
  rmsForExposedPixels(other.rmsForExposedPixels),
  meanForOverscanPixels(other.meanForOverscanPixels),
  rmsForOverscanPixels(other.rmsForOverscanPixels)
{}


