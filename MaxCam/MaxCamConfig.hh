#ifndef MAXCAM_CONFIG_HH
#define MAXCAM_CONFIG_HH

#include "TNamed.h"


class MaxCamConfig : public TNamed {

public:

  // Ctors
  MaxCamConfig() {};

  MaxCamConfig(const char *name, const char *title);

  MaxCamConfig(const MaxCamConfig &other);

  virtual ~MaxCamConfig() {};


  long cameraID; // Number for camera currently in use.

  long row_width; // Number of bins in ccd chip row. (does NOT include overscan)
  long img_rows;  // Number of bins in ccd chip col. (does NOT include overscan)

  long hbin; // Horizontal bin factor for pixels.
  long vbin; // Vertical bin factor for pixels.

  long ul_x; // upper-left x  in pixels (does NOT include overscan)
  long ul_y; // upper-left y  in pixels (does NOT include overscan)
  long lr_x; // lower-right x in pixels (does NOT include overscan)
  long lr_y; // lower-right y in pixels (does NOT include overscan)

  double CCDTemp; // Actual CCD temperature
  double CCDTempSet; // Goal CCD temperature
  
  long exposureTime; // Duration of exposure in ms
  
  enum FrType { normal, dark };
  int  frameType; // Normal (shutter open) or dark frame

  enum ClType { none, flush, bg};
  int  cleanType; //  Camera cleaning type 

  int nFlushes; // number of flushes before exposure.

  int bitDepth; // 8-bit or 16-bit pixels are supported.

  double daqTime; // total time for DAQ time ( = flushing + exposure + readout).

  TString serialNumber;

  unsigned short overscanColumns; // numer of column bins that are not exposed but digitized
  unsigned short overscanRows;    // number of row bins that are not exposed but digitized
  bool  digitizeOverscan;         // apply digitization to overscan pixels

  float meanForExposedPixels;   // mean for exposed vixels
  float rmsForExposedPixels;    // rms for exposed vixels
  float meanForOverscanPixels;  // mean for masked vixels
  float rmsForOverscanPixels;   // rms for masked vixels


  ClassDef(MaxCamConfig,2)
};

#endif

  
//For example, if you set the binning to be 4x4, then you'd have:
//
//  cam->getConfiguration()->ul_x, ul_y = 0, 0
//  cam->getConfiguration()->lr_x, lr_y = 1024, 1024
//  cam->getConfiguration()->hbin, vbin = 4, 4
//  cam->getConfiguration()->row_width, img_rows = 256, 256
//  cam->getConfiguration()->overscanColumns = 2
//  cam->getConfiguration()->overscanRows = 0
