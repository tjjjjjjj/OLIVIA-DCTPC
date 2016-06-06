#include "MaxCamAlta.hh"
#include "MaxCamConfig.hh"
#include <iostream>
#include <fstream>
#include <vector>
using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::flush;
#include <cmath>
#include "TString.h"
#include "TFile.h"
#include "TH2.h"
#include "TStopwatch.h"

//// begin camera block
#include "fitsio.h"
#include "png.h"
#include "ApnCamera.h"
#include "ApogeeUsb/ApogeeUsb.h"

// libccd:
#include "tcl.h"
#include "ccd.h"
extern int bias_start, bias_end, bcols;
typedef struct {
     unsigned short *pixels;
     int            size;
     short          xdim;
     short          ydim;
     short          zdim;
     short          xbin;
     short          ybin;
     short          type;
     char           name[64];
     int            shmid;
     size_t         shmsize;
     char           *shmem;
} CCD_FRAME;
typedef void *PDATA;
#define MAX_CCD_BUFFERS  1000
PDATA CCD_locate_buffer(char *name, int idepth, short imgcols, short imgrows, short hbin, short vbin);
int   CCD_free_buffer();
int   CCD_locate_buffernum(char *name);
extern CCD_FRAME CCD_Frame[MAX_CCD_BUFFERS];
extern int CCD_free_buffer(char *name);
/////////////////////////////////////////////////// end camera block



ClassImp(MaxCamAlta)

//____________________
// Class to run Apogee Alta camera. 
//
// *** Example:
//
//
MaxCamAlta::MaxCamAlta(int debug)  {
  // Constructor

  _alta = (CApnCamera*)new CApnCamera();
  //_alta->debugFlag=debug;
}



int 
MaxCamAlta::openCamera(int i) {
  // Prepare camera for taking data: estalish a connection and set the image area;
  cout << "openCamera():  i = " << i << endl;
  //cout << "_alta->debugFlag = " << _alta->debugFlag << endl;
  assert( _alta->InitDriver(i,0,0) ); 
  assert(_alta->ResetSystem() );
  //	Single image per download 
  _alta->write_ImageCount(1); ;
  setFanSpeed(2);

  copyConfiguration();

  return 0;
}


int
MaxCamAlta::copyConfiguration() {

  config->cameraID = _alta->m_CamIdA;

  config->hbin=_alta->read_RoiBinningH(); 
  config->vbin=_alta->read_RoiBinningV();

  config->ul_x=_alta->read_RoiStartX();
  config->ul_y=_alta->read_RoiStartY(); 
  config->lr_x=_alta-> read_RoiPixelsH();
  config->lr_y=_alta->read_RoiPixelsV(); 

  config->row_width=(config->lr_x-config->ul_x)/config->hbin;
  config->img_rows=(config->lr_y-config->ul_y)/config->vbin;  

  config->CCDTemp=getTemperature();
  config->CCDTempSet=getGoalTemperature(); 
  
  config->exposureTime=1000;
  
  config->frameType=MaxCamConfig::dark;

  //enum ClType { none, flush, bg};
  //int  config->cleanType; //  Camera cleaning type 

  //int config->nFlushes; // number of flushes before exposure.

  config->bitDepth=_alta->read_DataBits() ? 0xFFF : 0xFFFF; // 16bit=0, 12bit=1

  //double config->daqTime; // total time for DAQ time ( = flushing + exposure + readout).

  //config->overscanColumns=_alta->m_ApnSensorInfo->m_OverscanColumns;
  //config->overscanRows=_alta->m_ApnSensorInfo->m_OverscanRows;
  //config->digitizeOverscan=_alta->m_pvtDigitizeOverscan;
  // by default, no overscan region is recorded
  config->overscanColumns=0; //N_OVERSCAN_COLUMNS;
  config->overscanRows=0;    //should eventually be N_OVERSCAN_ROWS;
  config->digitizeOverscan=false;

  char sernumber[256];
  long buflength=256;
  _alta->GetCameraSerialNumber( sernumber, &buflength);
  cout << "sernumber = " << sernumber << endl;
  cout << "_alta->m_SysDeviceName = " << _alta->m_SysDeviceName << endl;
  //config->serialNumber=TString(_alta->m_SysDeviceName);
  config->serialNumber=TString(sernumber);

  cout <<" +++++++++++++++++++++ SETUP ++++++++++++++++++" << endl;
  cout << "[" << _alta->m_pvtProductId<<"]  [" << _alta->m_pvtDeviceId<<"]  [" <<_alta->m_SysDeviceName << "]"
       <<endl;

  cout << "[";
  for (int i=0; i<buflength; i++) cout << int(sernumber[i]);
  cout << "]" << endl;

  return 0;
}



int 
MaxCamAlta::closeCamera() {
  // Close camera handle.

  _alta->CloseDriver();

  return 0;
}


int
MaxCamAlta::clickCamera() {

  int status = expose();	    

  status = grabImage();

  return status;
}


int
MaxCamAlta::expose() {
  bool shutter = 1;
  if (config->frameType==MaxCamConfig::dark) shutter=0;
 return _alta->Expose( config->exposureTime*0.001, shutter);	    
}


int
MaxCamAlta::grabImage() {

  /*	    Readout the image and save in a named buffer (tempobs) */
  int status = _alta->BufferImage("tempobs"); 
    
  /*	    Use the libccd routine to find the corresponding buffer index */
  _bnum = CCD_locate_buffernum("tempobs");
  
  /*	    Print details about the buffer for debug purposes */
  //printf("Buffer %4d %s = %d bytes cols=%d rows=%d depth=%d\n",_bnum,CCD_Frame[_bnum].name,CCD_Frame[_bnum].size,CCD_Frame[_bnum].xdim,CCD_Frame[_bnum].ydim,CCD_Frame[_bnum].zdim);

  return status;
}


TH2S*
MaxCamAlta::createHisto(TString histoName) {
  // Create a histrogram from the image array in the memory.

  // After exposing the Alta U6 camera, the image data is "downloaded" from
  // camera to computer using MaxCamAlta::grabImage().  This puts the image
  // data into a structure called CCD_Frame.  There is a 1-d vector of
  // pixel values that represents the entire image.  This is called
  // CCD_Frame[_bnum].pixels.  
  // 
  // The data must then be put into a TH2S.  This happens in
  // MaxCamAlta::createHisto().  The overscan region is put into a TH2S in
  // createOverscanHisto().  If you want to read out the entire vector
  // (exposed pixels plus overscan) into a TH2S, then use createFullHisto()
  // 
  //   unsigned short *image = CCD_Frame[_bnum].pixels;
  // 
  // The 1D array is indexed such that pixels in a row are the fast-moving
  // index.  For example, if you have a 4x4 array (row x col) of exposed
  // pixels and a 4x2 array of overscan pixels, then the indexes are:
  // 
  // 
  //     <--- exposed -------><-- o.s.-->
  // 
  //      ------------------------------
  //     |  0 |  1 |  2 |  3 ||  4 |  5 | 
  //     |------------------------------
  //     |  6 |  7 |  8 |  9 || 10 | 11 |
  //     |------------------------------
  //     | 12 | 13 | 14 | 15 || 16 | 17 |
  //     |------------------------------
  //     | 18 | 19 | 20 | 21 || 22 | 23 |
  //      ------------------------------
  // 
  // In the early days, it was decided (by Denis likely) that the displayed
  // image in the TH2S should match the view that you would see if you
  // looked through the camera from the back side (fans) forward.  This
  // means that the x-axis (cols) on the chip become the y-axis (rows) in
  // the TH2S and vice versa.  It also requires that the TH2S x-axis be
  // mirrored.
  // 
  // So, when you populate the TH2S with data from the 1D image array, you need to compute three values:
  //   1. the index into the 1D image array (called "pindex" in our code)
  //   2. the xbin number for the TH2S (1-indexed)
  //   3. the ybin number for the TH2S (1-indexed)
  //
  //
  // If your version of ROOT does not have TH2S::SetBinContent(xbin, ybin, value), 
  // then try: SetBinContent(GetBin(xbin,ybin), value)

  unsigned short *image = CCD_Frame[_bnum].pixels;
  int nx_bin = CCD_Frame[_bnum].xdim;   // number of x bins (includes overscan columns, if active)
  int ny_bin = CCD_Frame[_bnum].ydim;   // number of y bins 

  config->meanForExposedPixels=config->rmsForExposedPixels=0;

  img_histo = new TH2S(histoName, "", 
		       ny_bin, 0, config->lr_y-config->ul_y,  
		       nx_bin-config->overscanColumns, 0, config->lr_x-config->ul_x );

  //cout << "nx_bin, ny_bin = " << nx_bin << ", " << ny_bin << endl;
  //cout << "config->overscanColumns = " << config->overscanColumns << endl;
  //cout << "config->lr_y, ul_y, lr_y - ul_y = " << config->lr_y << ", " << config->ul_y << ", " << config->lr_y-config->ul_y << endl;
  //cout << "config->lr_x, ul_x, lr_x - ul_x = " << config->lr_x << ", " << config->ul_x << ", " << config->lr_x-config->ul_x << endl;
  //cout << "img_histo->GetNbinsX() = " << img_histo->GetNbinsX() << endl;
  //cout << "img_histo->GetNbinsY() = " << img_histo->GetNbinsY() << endl;

  int nimgcols = nx_bin-config->overscanColumns;
  int pindex=0;
  for (int icol=0; icol<nimgcols; icol++) {
    for (int jrow=0; jrow<ny_bin; jrow++) {
      pindex=nx_bin * (ny_bin-jrow-1) + icol;
      if (image[pindex]>=config->bitDepth) cerr << "MaxCamAlta: saturated point "<<icol<<","<<jrow<<"  maximum="<<config->bitDepth<<endl;
      //img_histo->SetBinContent((ny_bin+2)*(icol+1) + ny_bin-jrow,(float)image[pindex]);
      img_histo->SetBinContent(ny_bin-jrow, icol+1, short(image[pindex]));

      config->meanForExposedPixels+=(float)image[pindex];
      config->rmsForExposedPixels+=(float)(image[pindex]*image[pindex]);
    }
  }

  //=_meanForMaskedPixels=_rmsForMaskedPixels

  config->meanForExposedPixels /= (ny_bin*nx_bin);
  config->rmsForExposedPixels /= (ny_bin*nx_bin);
  config->rmsForExposedPixels =sqrt( config->rmsForExposedPixels - config->meanForExposedPixels*config->meanForExposedPixels );
  return img_histo;
}

TH2S*
MaxCamAlta::createOverscanHisto(TString histoName) {

  /* 
     Code to read out part of the overscan region of the ccd chip
     The bin numbering of resulting image is chosen to match that of the 
     imaging pixels.
     For example, if the imaging pixels cover 0,1024 in x and y, then 
     the overscan region will cover 1024-1032 in y and 0-1024 in x.
     See WARNING #2 about the meaning of x and y w.r.t. ccd chip rows and columns.
     
     ************* WARNING #0 ***************
     This fxn has not been tested with on-chip binning....
     ****************************************

     ************* WARNING #1 ***************
     only 8 overscan columns are used here (cols 1025-1032)
     The first 4 "under"scan columns are not used.
     Nor are any of the underscan or overscan rows
     ****************************************
     
     ************* WARNING #2 ***************
     A row on the ccd chip corresponds to a 
     column in the TH2S.  While I am not thrilled
     with this convention, I followed it since
     all historical data was taken in this configuration.
     see, e.g. createHisto()
     ****************************************

     ************* WARNING #3 ***************
     Getting the pixels-to-bin mapping correct 
     is an absolute nightmare.  Proceed at your
     own risk.
     ****************************************

     ************* WARNING #4 ***************
     Should you want to extend this fxn to include 
     more overscan pixels, be aware that you will
     have to not only edit this routine, but also
     setDigitizeOverscan() and, likely, 
     createHisto()
     ****************************************
     
     The Apogee Alta U6 uses a Kodak KAF-1001 chip
     (actually, it's the KAF-10011105, see:
     apogee/src/apogee/ApnCamData_KAF10011105.cpp)
     
     The chip has several rows and columns of covered (unexposed) pixels
     that can be used to estimate the pedestal, read noise and dark noise
     on an image-by-image basis.  
     
     The geometry of ccd chip in the Apogee Alta U6
     including overscan region, is shown below.
     
     The inner box is the imaging pixels
     The outer box is the overscan region
     The dimensions are:
     imaging pixels:  1024 x 1024
     unexposed region:
       4 rows above imaging pixels
       4 rows below imaging pixels
       4 columns to the left  of imaging pixels
       8 columns to the right of imaging pixels
     So there are 12 overscan columns (well, 4 underscan and 8 overscan).
     And there are 8 overscan rows (well, 4 underscan and 4 overscan).
    ________________________
   |                        |
   |  ------------------    |
   | |                  |   |
   | |                  |   |
   | |                  |   |
   | |     imaging      | <-|---- non-imaging pixels 
   | |     pixels       |   |
   | |                  |   |
   | |                  |   |
   | |                  |   |
   | |__________________|   |
   |                        |
    ------------------------

  */

  unsigned short *image = CCD_Frame[_bnum].pixels;
  int ncol_bin = CCD_Frame[_bnum].xdim; //# of col bins on chip --> y of TH2S  (includes imaging region and overscan cols if active)
  int nrow_bin = CCD_Frame[_bnum].ydim; //# of row bins on chip --> x of TH2S  (only imaging region)

  config->meanForOverscanPixels=config->rmsForOverscanPixels=0;

  //cout << "\n\n" << endl;
  //cout << "createOverscanHisto()" << endl;
  //cout << "nrow_bin, ncol_bin = " << nrow_bin << ", " << ncol_bin << endl;
  //cout << "config->overscanColumns = " << config->overscanColumns << endl;
  //cout << "config->overscanRows = " << config->overscanRows << endl;
  //cout << "config->row_width = " << config->row_width << endl;
  //cout << "config->img_rows  = " << config->img_rows << endl;
  //cout << "bias_start = " << bias_start << endl;
  //cout << "bias_end   = " << bias_end << endl;
  //cout << "_alta->m_pvtRoiPixelsH = " << _alta->m_pvtRoiPixelsH << endl;
  //cout << "_alta->m_ApnSensorInfo->m_PostRoiSkipColumns = " << _alta->m_ApnSensorInfo->m_PostRoiSkipColumns << endl;
  //cout << "config->lr_y = " << config->lr_y << endl;
  //cout << "config->ul_y = " << config->ul_y << endl;
  //cout << "config->lr_y-config->ul_y = " << config->lr_y-config->ul_y << endl;

  int n_overscan_col_bins     = config->overscanColumns;
  int overscan_col_bin_start  = config->row_width;  // 0-indexed bin number of first overscan column (FAU6BB4=256)
  //cout << "n_overscan_col_bins    = " << n_overscan_col_bins << endl;
  //cout << "overscan_col_bin_start = " << overscan_col_bin_start << endl;

  img_histo_overscan = new TH2S(histoName, "", 
				nrow_bin, 0, config->lr_y-config->ul_y, 
				n_overscan_col_bins, bias_start-1, bias_end);

  // FAU6BB4 = "For Alta U6 Binned By 4"
  int pindex=0;
  int xbin, ybin;  
  for (int icol=0; icol<n_overscan_col_bins; icol++) {  // chip cols, histogram y
    ybin = icol+1;  // y-index into TH2S
    for (int irow=0; irow<nrow_bin; irow++) {  // chip rows, histogram x
      xbin = nrow_bin-irow;  // x-index into TH2S

      // 1-D index into Alta image array
      pindex = ncol_bin * (nrow_bin-irow-1) + icol+overscan_col_bin_start;
      img_histo_overscan->SetBinContent(xbin, ybin, short(image[pindex]));

      config->meanForOverscanPixels+=(float)image[pindex];
      config->rmsForOverscanPixels+=(float)(image[pindex]*image[pindex]);
    }
  }

  config->meanForOverscanPixels /= (config->overscanColumns*nrow_bin);
  config->rmsForOverscanPixels /= (config->overscanColumns*nrow_bin);
  config->rmsForOverscanPixels =sqrt( config->rmsForOverscanPixels - config->meanForOverscanPixels*config->meanForOverscanPixels );

  return img_histo_overscan;
}


TH2S*
MaxCamAlta::createFullHisto(TString histoName) {

  unsigned short *image = CCD_Frame[_bnum].pixels;
  int ncol = CCD_Frame[_bnum].xdim;  //cols of chip --> y of TH2F
  int nrow = CCD_Frame[_bnum].ydim;  //rows of chip --> x of TH2F

  int hist_y_max = config->lr_x-config->ul_x;
  if (bias_end) hist_y_max = bias_end;

  img_histo_full = new TH2S(histoName, "", 
			    nrow, 0, config->lr_y-config->ul_y, 
			    ncol, 0, hist_y_max);
			    //ncol, 0, (config->lr_x-config->ul_x)+config->overscanColumns);
  int pindex=0;
  int xbin, ybin;  
  for (int icol=0; icol<ncol; icol++) {    // histogram y
    ybin = icol+1;
    for (int irow=0; irow<nrow; irow++) {  // histogram x
      xbin = nrow-irow;
      pindex = ncol * (nrow-irow-1) + icol;
      img_histo_full->SetBinContent(xbin, ybin, short(image[pindex]));
    }
  }

  return img_histo_full;
}

int MaxCamAlta::setHBin(long hb) { 
  // Set horizontal binning.
  config->hbin=hb; 
  config->row_width=(config->lr_x-config->ul_x)/config->hbin;

  _alta->m_pvtRoiPixelsH /= hb;
  _alta->write_RoiBinningH(hb);

  config->overscanColumns /= hb;

  return 0; 
}

int MaxCamAlta::setVBin(long vb) { 
  // Set vertical binning.
  config->vbin=vb; 
  config->img_rows=(config->lr_y-config->ul_y)/config->vbin;  

  _alta->m_pvtRoiPixelsV /= vb;
  _alta->write_RoiBinningV(vb);


  return 0; 
}

int MaxCamAlta::resetWithFlush() {
  _alta->ResetSystem();

  //cout << "_alta->m_pvtRoiPixelsH = " << _alta->m_pvtRoiPixelsH << endl;
  //cout << "_alta->m_pvtRoiPixelsV = " << _alta->m_pvtRoiPixelsV << endl;
  //cout << "_alta->read_RoiBinningH() = " << _alta->read_RoiBinningH() << endl;
  //cout << "_alta->read_RoiBinningV() = " << _alta->read_RoiBinningV() << endl;
  //cout << "_alta->read_DigitizeOverscan() = " << _alta->read_DigitizeOverscan() << endl;

  return 0;
}

int MaxCamAlta::reset() {
  _alta->ResetSystemNoFlush();
  return 0;
}

int 
MaxCamAlta::setExposureTime(long exptime) { 
  // Set exposure time (in msec).
  config->exposureTime=exptime;
  return 0; 
}

int
MaxCamAlta::cancelExposure() {
  // Stop the current exposure.
  return _alta->StopExposure(false);
}


double
MaxCamAlta::getTemperature() {
  // Read the actual temperature of the CCD (in deg C).
  config->CCDTemp= _alta->read_TempCCD();
  return config->CCDTemp;
}


int
MaxCamAlta::setTemperature(double temp) {
  // Set CCD temperature (in deg C).
  config->CCDTempSet=temp;
  _alta->write_CoolerEnable(1);
  _alta->write_CoolerSetPoint( temp );
  return 0;
}

double 
MaxCamAlta::getGoalTemperature() { 
  // Return the desired temperature (in deg C).
  return _alta->read_CoolerSetPoint(); 
}


int
MaxCamAlta::setDataBits(int bits) {
  // Set depth of ADC converter on camera.
  // only 12 and 16 bits are allowed
  switch (bits) {
  case 12: _alta->write_DataBits( Apn_Resolution_TwelveBit ); break;
  case 16: _alta->write_DataBits( Apn_Resolution_SixteenBit ); break;
  default: cout << "MaxCamAlta::setDataBits: wrong bit depth" << endl; return -1;
  }
  return 0;
}

int
MaxCamAlta::getDataBits() {
  // Get bith depth
  Apn_Resolution bitdepth=_alta->read_DataBits();
  return bitdepth ? 12 : 16;
}

int
MaxCamAlta::getGain() {
  // ADC gain
  int bits=getDataBits();
  switch(bits) {
  case 12: return  _alta->read_TwelveBitGain();
  case 16: return 0;//return  _alta->read_Alta2ADGainSixteen();
  default: assert(0);
  }
  return 0;
}

int
MaxCamAlta::setGain(unsigned short gain) {
  // ADC gain
  int bits=getDataBits();
  switch(bits) {
  case 12:   _alta->write_TwelveBitGain(gain); break;
  case 16:   break;//_alta->write_Alta2ADGainSixteen(gain); break;
  default: assert(0);
  }
  return 0;
}


int
MaxCamAlta::setDarkFrame() {
  // Do not open shutter during exposure.
  config->frameType=MaxCamConfig::dark;
  _alta->write_CameraMode(Apn_CameraMode_Normal);
  return 0;
}

int
MaxCamAlta::setNormalFrame() {
  // Open shutter during exposure.
  config->frameType=MaxCamConfig::normal;
  _alta->write_CameraMode(Apn_CameraMode_Normal);
  return 0; 
}

int
MaxCamAlta::setTestFrame() {
  // Flag camera as test.
  _alta->write_CameraMode(Apn_CameraMode_Test);
  return 0;
}



int
MaxCamAlta::setFanSpeed(int speed) {
  if (speed > Apn_FanMode_High)     speed=Apn_FanMode_High;
  else if (speed < Apn_FanMode_Off) speed=Apn_FanMode_Off;
  _alta->write_FanMode(speed);
  return 0;
}




int 
MaxCamAlta::flushCamera(long repeat) {
  // Flush ccd image...

  // Start flushing
  //_alta->Flush();
     
  return 0;
}

void
MaxCamAlta::print() {
  //_alta->Print();
}


unsigned short  MaxCamAlta::nCamera=0;

#include "usb.h" 
#define USB_ALTA_VENDOR_ID      0x125c
#define USB_ALTA_PRODUCT_ID     0x0010


void
MaxCamAlta::discoverCameras() {
  
  usb_init();
  
  usb_find_busses();
  usb_find_devices();
  
  nCamera = 0;
  
  /* find ALTA device */
  for(struct usb_bus *bus = usb_busses; bus; bus = bus->next) {
    cout << "found bus" << endl;
    for(struct usb_device *dev = bus->devices; dev; dev = dev->next) {
      cout << "found device" << endl;
      if (dev->descriptor.idVendor == USB_ALTA_VENDOR_ID && 
	  dev->descriptor.idProduct == USB_ALTA_PRODUCT_ID) {
	cout << "found apogee"<<endl;
	cout << "file name .......  " << dev->filename << endl; 
	nCamera++;
      }
    }
  }
  cout << "Total ccd found = " << nCamera << endl;

}


bool
MaxCamAlta::imageReady() {
  // Is camera image ready?

  return _alta->ImageReady();
}


int
MaxCamAlta::openShutter() {
  _alta->write_ForceShutterOpen(true);
  return 0;
}


int
MaxCamAlta::closeShutter() {
  _alta->write_ForceShutterOpen(false);
  return 0;
}


int
MaxCamAlta::setDigitizeOverscan(bool digitize) {

  config->digitizeOverscan=digitize;
  _alta->write_DigitizeOverscan( digitize );

  int n_bias_cols = digitize ? N_OVERSCAN_COLS : 0;
  int n_bias_rows = digitize ? 0 : 0;  // shoudl eventually use N_OVERSCAN_ROWS...
  cout << "n_bias_cols = " << n_bias_cols << endl;

  config->overscanColumns = n_bias_cols;
  config->overscanRows = n_bias_rows;   
  
  // the *pixel* number of the start and end of the overscan region
  // this is "1" indexed.
  // For the Apogee Alta U6, if the overscan is active then:
  //   bias_start = 1025
  //   bias_end   = 1032
  bias_start = digitize ? (_alta->m_ApnSensorInfo->m_ImagingColumns+1) : 0;
  bias_end   = digitize ? (_alta->m_ApnSensorInfo->m_ImagingColumns+n_bias_cols) : 0;

  _alta->m_pvtRoiPixelsH = _alta->m_ApnSensorInfo->m_ImagingColumns + n_bias_cols;
  _alta->m_ApnSensorInfo->m_PostRoiSkipColumns = n_bias_cols;


  cout << "bias_start, bias_end = " << bias_start << ", " << bias_end << endl;
  cout << "_alta->m_pvtRoiPixelsH = " << _alta->m_pvtRoiPixelsH << endl;
  cout << "_alta->m_ApnSensorInfo->m_PostRoiSkipColumns = " << _alta->m_ApnSensorInfo->m_PostRoiSkipColumns << endl;


  cout << " *******************************************************" << endl;
  cout << "  WARNING:  Need to propagate overscan info into config" << endl;
  cout << "     setDigitizeOverscan() " << endl;
  cout << " *******************************************************" << endl;

  return 0;
}
