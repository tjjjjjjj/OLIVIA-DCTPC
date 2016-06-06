#include "MaxCam.hh"
#include "MaxCamConfig.hh"
#include <iostream>
#include <fstream>
#include <vector>
using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::flush;

#include "TString.h"
#include "TFile.h"
#include "TH2.h"
#include "TStopwatch.h"

//// begin libfli block - move to header file 

#include "libfli.h"
#include "png.h"
#include "fitsio.h"

typedef struct {
  flidomain_t domain;
  char *dname;
  char *name;
} cam_t;

vector<cam_t*> camList;
flidev_t device;

extern const char *__progname;

#define LIBVERSIZ 1024

////// end libfli block


ClassImp(MaxCam)

//____________________
// Class to run MaxCam camera. 
// (Finger Lakes Instrumentation, http://www.fli-cam.com)
//
// *** Example:
//  
//  // Camera setup:
//  MaxCam cam;
//  cam.openCamera();
//  cam.setNumberOfFlushes(1); // cleaning of CCD before exposure
//  cam.setTemperature(-20);   // CCD temperature
//  cam.setExposureTime(1000); // exposure in msec
//  cam.setHBin(16);           // x binning
//  cam.setVBin(16);           // y binning
//  cam.setDarkFrame();        // or setNormalFrame (default)
//
//  // take data:
//  cam.clickCamera();
//
//  // access CCD histogram:
//  TH2F* hCCD=cam.createHisto("hCCD");
//  hCCD->Draw("colz");
//
//  // close connection
//  cam.closeCamera();
//
//
  MaxCam::MaxCam(int debugLevel) : numcams(0) {
  // Constructor

  config = new MaxCamConfig("ccdconfig","CCD Camera configuration");;
  img_histo = new TH2S;
  img=0;

  char  libver[LIBVERSIZ];  

  FLISetDebugLevel(NULL, debugLevel ? FLIDEBUG_ALL : FLIDEBUG_NONE);

  FLIGetLibVersion(libver, LIBVERSIZ);
  //cout << "MaxCam: Library version " << libver << endl;

  findCams();

  totalExposures=0;

  _sw = new TStopwatch;
}


void 
MaxCam::findCams() {
  // Search USB ports for attached cameras.

  char **tmplist;
  FLIList(FLIDOMAIN_USB | FLIDEVICE_CAMERA, &tmplist);
  cout << "+++++++ FLIList found # of cameras = " << tmplist << endl;

  int ncams=0;
  cout << "+++++++ FLIList found # of cameras = " << tmplist[0][0] << endl;
  if (tmplist && tmplist[0] ) {
  
    for (int i = 0; tmplist[i]; i++) ncams++;
    cout << "+++++++ FLIList found # of cameras = " << ncams << "    " << tmplist[0]<< endl;
        
    for (int i = 0; tmplist[i]; i++) {


      for (int j = 0; tmplist[i][j] != '\0'; j++) {
	if (tmplist[i][j] == ';') {
	  tmplist[i][j] = '\0';
	  break;
	}
      }
      
      cam_t *tmpcam = new cam_t;
      tmpcam->domain = FLIDOMAIN_USB;
      tmpcam->dname = "USB";
      tmpcam->name = strdup(tmplist[i]);

      camList.push_back(tmpcam);
    }

    numcams += ncams;
  }

  cout << "Found # of cameras = " << numcams << endl;

  FLIFreeList(tmplist);
  
}


int 
MaxCam::openCamera(int i) {
  // Prepare camera for taking data: estalish a connection and set the image area;

  if (!numcams) {
    cout << "MaxCam: no camera attached to the system" << endl;
    return 1;
  }
  if (i>=numcams || i<0) {
    cout << "MaxCam: illegal camera number " << i
	 << " allowed range " << "0," << numcams << endl;
    return 2;
  }
  FLIOpen( &device, camList[i]->name, FLIDEVICE_CAMERA | camList[i]->domain);

  setImageArea();

  return 0;
}

int 
MaxCam::closeCamera() {
  // Close camera handle.
  FLIClose( device);
  return 0;
}


char* 
MaxCam::getCameraModel() {
  // Returns camera model.
  if (!device) return 0;
#define BUFF_SIZ 4096
  char *buff= new char[BUFF_SIZ];
  FLIGetModel( device, buff, BUFF_SIZ);
  return buff;
}

long 
MaxCam::getHWRevision() {
  // Returns camera hardware revision
  if (!device) return 0;
  long tmp;
  FLIGetHWRevision( device, &tmp);
  return tmp;
}

long 
MaxCam::getFWRevision() {
  // Returns camera firmware revision
  if (!device) return 0;
  long tmp;
  FLIGetFWRevision( device, &tmp);
  return tmp;
}


int 
MaxCam::setExposureTime(long exptime) { 
  // Set exposure time (in msec).
  config->exposureTime=exptime;
  return FLISetExposureTime( device, exptime); 
}

int
MaxCam::cancelExposure() {
  // Stop the current exposure.
  return FLICancelExposure( device );
}

int
MaxCam::cleanCCD() {
  // Make one "quick" readout to clean camera before real exposure.
  long vbin=config->vbin;
  long expotime = config->exposureTime;
  setExposureTime(0);
  setVBin(512);
  clickCamera();
  setVBin(vbin);
  setExposureTime(expotime);
  return 0;
}

int 
MaxCam::setNumberOfFlushes(long nflush) { 
  // Set number of times the CCD is readout before exposure.
  // The cleaning of the CCD before exposure is not necessary if
  // background is continuously flushed (see startBkgFlush and stopBkgFlush).
  config->nFlushes=nflush;
  config->cleanType=MaxCamConfig::flush;
  return FLISetNFlushes(device, nflush); 
}



int
MaxCam::startBkgFlush() {
  // Start continuous cleaning of the CCD.
  config->cleanType=MaxCamConfig::bg;
  return FLIControlBackgroundFlush( device, FLI_BGFLUSH_START);
}


int
MaxCam::stopBkgFlush() {
   // Stop continuous cleaning of the CCD.
  config->cleanType=MaxCamConfig::none;
  return FLIControlBackgroundFlush( device, FLI_BGFLUSH_STOP);
}


int
MaxCam::flushCamera(long repeat) {
  long rows=config->img_rows;
  return FLIFlushRow( device, rows, repeat);
}

int
MaxCam::openShutter() {
  // Open shutter.
  return FLIControlShutter( device, FLI_SHUTTER_OPEN);
}

int 
MaxCam::closeShutter() {
  // Close shutter.
  return FLIControlShutter( device, FLI_SHUTTER_CLOSE); 
}


int
MaxCam::triggerShutter() {
  return FLIControlShutter( device, FLI_SHUTTER_EXTERNAL_TRIGGER);  
}


int
MaxCam::setDarkFrame() {
  // Do not open shutter during exposure.
  config->frameType=MaxCamConfig::dark;
  return FLISetFrameType( device, FLI_FRAME_TYPE_DARK);
}

int
MaxCam::setNormalFrame() {
  // Open shutter during exposure.
  config->frameType=MaxCamConfig::normal;
  return FLISetFrameType( device, FLI_FRAME_TYPE_NORMAL); 
}


double
MaxCam::getTemperature() {
  // Read the actual temperature of the CCD (in deg C).
  double temp;
  FLIGetTemperature(device, &temp);
  config->CCDTemp=temp;
  return temp;
}

// double
// MaxCam::readTemperature() {
//   // Get CCD temperature.
//   double temp;
//   FLIReadTemperature(device, &temp);
//   return 0;
// }

int
MaxCam::setTemperature(double temp) {
  // Set CCD temperature (in deg C).
  config->CCDTempSet=temp;
  return FLISetTemperature(device, temp);
}


double 
MaxCam::getGoalTemperature() { 
  // Return the desired temperature (in deg C).
  return config->CCDTempSet; 
}


double
MaxCam::getPixelSize(TString what) {
  // Returns x, y pixel sizes. 
  // Arguments are "x", "y".
  what.ToUpper();
  if (!device) return -1;
  double xsize, ysize;
  FLIGetPixelSize( device, &xsize, &ysize);
  return what=="X" ? xsize : ysize;
}

long
MaxCam::getArrayArea(TString what) {
  // Returns coordinate of active CCD area.
  // Arguments are "ULX", "ULY", "LRX", "LRY".
  what.ToUpper();
  long tmp1, tmp2, tmp3, tmp4;
  FLIGetArrayArea( device, &tmp1, &tmp2, &tmp3, &tmp4);
  long ret=-1;
  if      (what=="ULX") ret=tmp1;
  else if (what=="ULY") ret=tmp2;
  else if (what=="LRX") ret=tmp3;
  else if  (what="LRY") ret=tmp4;
  else cerr << "MaxCam: Wrong coordinate type" << endl;
  return ret;
}

long
MaxCam::getVisibleArea(TString what) {
  // Returns coordinate of active CCD area.
  // Arguments are "ULX", "ULY", "LRX", "LRY".
  what.ToUpper();
  long tmp1, tmp2, tmp3, tmp4;
  FLIGetVisibleArea( device, &tmp1, &tmp2, &tmp3, &tmp4);
  long ret=-1;
  if      (what=="ULX") ret=tmp1;
  else if (what=="ULY") ret=tmp2;
  else if (what=="LRX") ret=tmp3;
  else if  (what="LRY") ret=tmp4;
  else cerr << "MaxCam: Wrong coordinate type" << endl;
  return ret;
}

int MaxCam::setHBin(long hb) { 
  // Set horizontal binning.
  config->hbin=hb; 
  FLISetHBin( device, hb);
  setImageArea();
  return 0; 
}

int MaxCam::setVBin(long vb) { 
  // Set vertical binning.
  config->vbin=vb; 
  FLISetVBin( device, vb); 
  setImageArea();
  return 0; 
}


int 
MaxCam::setImageArea(long ul_x, long ul_y, long lr_x, long lr_y) {
  // Set readout area. This function is called during initialization
  // without arguments (full frame). Setting smaller frame speeds up
  // the CCD readout.

  if (ul_x<0 || ul_y<0 || lr_x<0 || lr_y<0) 
    FLIGetVisibleArea( device, &ul_x, &ul_y, &lr_x, &lr_y);  

  config->ul_x=ul_x; config->ul_y=ul_y; config->lr_x=lr_x; config->lr_y=lr_y;

  config->row_width = (lr_x-ul_x) / config->hbin;
  config->img_rows  = (lr_y-ul_y) / config->vbin;

  FLISetImageArea( device, 
		   ul_x,                  ul_y, 
		   ul_x+(lr_x-ul_x)/config->hbin, ul_y+(lr_y-ul_y)/config->vbin);
  return 0;
}





int
MaxCam::clickCamera() {
  // Take a picture using current camera setup.

  _sw->Start();

  int ret = expose();
  ret = grabImage();

  _sw->Stop();
  config->daqTime=_sw->RealTime();
  
  return ret;
}


int
MaxCam::expose() {
  totalExposures++;

  if (config->frameType==MaxCamConfig::dark) {
    if (config->exposureTime<1) cout <<"b"<<flush;
    else cout << "d"<<flush;
  }
  else cout << "*" << flush;


  int ret=0;
  FLIExposeFrame( device );
  long timeLeft;
  do {
    ret = FLIGetExposureStatus( device, &timeLeft);
    if (ret) break;
    usleep(timeLeft * 1000);
  } while (timeLeft);

  return ret;
}


int
MaxCam::grabImage() {

  if (img) delete img;
  img = new u_int16_t[config->img_rows * config->row_width];

  int ret=0;
  for (int row = 0; row < config->img_rows; row++) {
    ret = FLIGrabRow( device, &img[row * config->row_width], config->row_width);
    if (ret) break;
  }
  return ret;
}


TH2S*
MaxCam::createHisto(TString histoName) {
  // Create a histrogram from the image array in the memory.

  img_histo = new TH2S(histoName, "", 
		       config->row_width, 0, config->lr_x-config->ul_x,  
		       config->img_rows,  0, config->lr_y-config->ul_y );

  for (int i=0; i<config->row_width; i++) {
    for (int j=0; j<config->img_rows; j++) {
      if (img[config->row_width * j + i]>=config->bitDepth) 
	cerr << "MaxCam: saturated point "<<i<<","<<j<<endl;
      img_histo->SetBinContent((config->row_width+2)*(j+1) + i + 1,short(img[config->row_width * (config->img_rows-j-1) + i]));
    }
    //cout << endl;
  }

  return img_histo;
}



int
MaxCam::writePNG(TString fileName) {
  // Write a png file. [note yet done]

  FILE *fp = NULL;
  if ((fp = fopen((const char*)fileName, "wb")) == NULL) {
    cerr << "MaxCam: Cannot open file " << fileName << endl;
    return -1;
  }

  png_structp pngptr = NULL;
  if ((pngptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
					NULL, NULL, NULL)) == NULL)
  {
    cerr << "MaxCam: Cannot create png write structure" << endl;
    return -1;
  }

  png_infop infoptr = NULL;
  if ((infoptr = png_create_info_struct(pngptr)) == NULL)
  {
    cerr << "MaxCam: Cannot create png info" << endl;
    return -1;
  }

  png_init_io(pngptr, fp);

  png_set_compression_level(pngptr, Z_BEST_COMPRESSION);

  png_set_IHDR(pngptr, infoptr, config->row_width, config->img_rows, 16, PNG_COLOR_TYPE_GRAY,
	       PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
	       PNG_FILTER_TYPE_DEFAULT);

  png_write_info(pngptr, infoptr);

  png_set_swap(pngptr);

  png_byte *row;
  long height=config->img_rows;
  for (row = (png_byte*)img; height > 0; row += config->row_width * sizeof(u_int16_t), height--)
    png_write_row(pngptr, row);

  png_write_end(pngptr, infoptr);

  if (fp) fclose(fp);

  if (pngptr) png_destroy_write_struct(&pngptr, &infoptr);

  return 0;
}


int
MaxCam::writeFITS(TString fileName) {
  // Write a fits file. 

  int status = 0;
  long naxes[2] = {config->row_width, config->img_rows};
  fitsfile *fp;

  fits_create_file(&fp, (const char*)fileName, &status);
  if (status)
  {
    fits_report_error(stderr, status);
    return -1;
  }

  fits_create_img(fp, SHORT_IMG, 2, naxes, &status);
  if (status)
  {
    fits_report_error(stderr, status);
    return -1;
  }

  fits_write_img(fp, TUSHORT, 1, config->row_width * config->img_rows, img, &status);
  if (status)
  {
    fits_report_error(stderr, status);
    return -1;
  }

  fits_close_file(fp, &status);
  if (status)
  {
    fits_report_error(stderr, status);
    return -1;
  }

  return 0;
}


int
MaxCam::writeROOT(TString fileName) {
  // Write a ROOT file.
  TFile f(fileName,"RECREATE");
  createHisto("tmpCCD");
  img_histo->Write();
  f.Write();
  f.Close();
  return 0;
}


void
MaxCam::printCameraInfo(int iCam) {
  // Print basic information about the camera and the
  // current setup.
  //

  if (openCamera(iCam)) return;

  cout << "MaxCam: Trying camera " << camList[iCam]->name
       << " from " << camList[iCam]->dname 
       << " domain" << endl; 


  cout << "MaxCam: Model " << getCameraModel() << endl;
  cout << "MaxCam: Hardware Rev. " << getHWRevision() << endl;
  cout << "MaxCam: Firmware Rev. " << getFWRevision() << endl;
  cout << "MaxCam: Pixel Size " << getPixelSize("x") 
       << "x" << getPixelSize("y") << endl; 
  cout << "MaxCam: Image area ("<<getVisibleArea("ULX")<<", "
       <<getVisibleArea("ULY")<<")("
       <<getVisibleArea("LRX")<<", "
       <<getVisibleArea("LRY")<<")"<<endl;



}


TH1F* 
MaxCam::createYieldHisto(TString histoName) {
  // Make 1D histogram of pixel yields from current image.
  //

  if (!img_histo) return 0;

  float min = img_histo->GetMinimum();
  float max = img_histo->GetMaximum();
  int nbin = int(max-min);

  TH1F* hY=new TH1F(histoName,"", nbin, min,max);

  int x0=1, x1=img_histo->GetXaxis()->GetNbins();
  int y0=1, y1=img_histo->GetYaxis()->GetNbins();

  for (int i=x0; i<=x1; i++) {
    for (int j=y0; j<=y1; j++) {
      hY->Fill( img_histo->GetBinContent(i,j) );
    }
  }

  return hY;
}
