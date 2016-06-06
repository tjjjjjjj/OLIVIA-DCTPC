#include "MaxCamAndor.hh"
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
#include "TH2F.h"
#include "TStopwatch.h"

#include "atmcdLXd.h"



ClassImp(MaxCamAndor)

//____________________
// Class to run Andor camera. 
//
// *** Example:
//
//
MaxCamAndor::MaxCamAndor(int debug)  {
  // Constructor
}



int 
MaxCamAndor::openCamera(int i) {
  // Prepare camera for taking data: estalish a connection and set the image area;
  cout << endl << endl 
       << "--------------------------------------------------------------"
       <<endl;
  cout << GetName() << ": ANDOR CAMERA" << endl;
  cout << GetName() << ": openCamera():  i = " << i << endl;

  config->cameraID=i;

  long int lCameraHandle;
  GetCameraHandle(i, &lCameraHandle);
  SetCurrentCamera(lCameraHandle);

  //Initialize CCD
  unsigned long error = Initialize("/usr/local/etc/andor");
  if(error!=DRV_SUCCESS){
    cout << "Initialisation error...exiting" << endl;
    assert(0);
  }

  //sleep to allow initialization to complete
  sleep(20); 
  
  //Set Read Mode to --Image--
  SetReadMode(4);

  //Set Acquisition mode to --Single scan--
  SetAcquisitionMode(1);
  
  //Set initial exposure time
  SetExposureTime(1);

  //Get Detector dimensions
  int width, height;
  GetDetector(&width, &height);

  // set gain to 1x
  setPreAmpGain(0);

  // set ADC to 16-bit readout
  setDataBits(16);

  // maintain temp after camera shutdown
  SetCoolerMode(1);

  //Initialize Shutter
  SetShutter(1, 0, 0, 0);
        
  //Setup Image dimensions  
  SetImage(1,1,1,width,1,height);

  // Rotate image
  SetImageRotate(1);

  // maximum buffer size (no binning)
  imageData = new at_32[width*height];

  // set minimum gain
  setGain(0);

  // set fan to maximum
  setFanSpeed(1);

  copyConfiguration();

  return 0;
}


int
MaxCamAndor::copyConfiguration() {

  int width, height;
  GetDetector(&width, &height);


  config->hbin=1; 
  config->vbin=1;

  config->ul_x=1;
  config->ul_y=1; 
  config->lr_x=width;
  config->lr_y=height; 

  config->row_width=(config->lr_x-config->ul_x+1)/config->hbin;
  config->img_rows =(config->lr_y-config->ul_y+1)/config->vbin;  


  config->CCDTemp=getTemperature();
  config->CCDTempSet=getGoalTemperature(); 
  
  calcDaqTime();
  
  config->frameType=MaxCamConfig::dark;

  config->bitDepth=0xFFFF;

  config->overscanColumns=0; //N_OVERSCAN_COLUMNS;
  config->overscanRows=0;    //should eventually be N_OVERSCAN_ROWS;
  config->digitizeOverscan=false;

  int serialNumber;
  GetCameraSerialNumber( &serialNumber );
  config->serialNumber=TString( serialNumber );
  cout << GetName() <<": Camera serial number = " << serialNumber << endl;

  return 0;
}



int 
MaxCamAndor::closeCamera() {
  // Close camera handle.
  ShutDown();
  return 0;
}


int
MaxCamAndor::clickCamera() {

  expose();

  grabImage();

  return 0;
}


int
MaxCamAndor::expose() {
  StartAcquisition();
  return 0;
}


int
MaxCamAndor::grabImage() {
  int status;
  GetStatus(&status);
  while(status==DRV_ACQUIRING) GetStatus(&status);

  unsigned long dataSize=config->row_width*config->img_rows;

  status=GetAcquiredData( imageData, dataSize);
  if (status!=DRV_SUCCESS) {
    cout << GetName() <<": Acq. Error Code = " << status << endl;
  }
  //SaveAsFITS("image.fits", 0);


  return 0;
}


TH2S*
MaxCamAndor::createHisto(TString histoName) {
  // Create a histrogram from the image array in the memory.

  int nx_bin = config->row_width;   // number of x bins (includes overscan columns, if active)
  int ny_bin = config->img_rows;   // number of y bins 

  config->meanForExposedPixels = config->rmsForExposedPixels = 0;

  img_histo = new TH2S(histoName, "", 
		       ny_bin, 0, config->lr_y-config->ul_y+1,  
		       nx_bin, 0, config->lr_x-config->ul_x+1 );

  int pindex=0;
  for (int icol=0; icol<nx_bin; icol++) {
    for (int jrow=0; jrow<ny_bin; jrow++) {

      pindex=nx_bin * (ny_bin-jrow-1) + icol;

      if (imageData[pindex]>config->bitDepth) 
	cerr << "MaxCamAndor: saturated point "
	     <<icol<<","<<jrow
	     <<"  maximum="<<config->bitDepth<<endl;

      //if (icol>nx_bin-10 && jrow>ny_bin-2) cout << "  " << imageData[pindex] << flush;
      img_histo->SetBinContent(jrow+1, icol+1, (short)imageData[pindex]);

      config->meanForExposedPixels+=(float)imageData[pindex];
      config->rmsForExposedPixels+=(float)(imageData[pindex]*imageData[pindex]);
    }
  }


  config->meanForExposedPixels /= (ny_bin*nx_bin);
  config->rmsForExposedPixels /= (ny_bin*nx_bin);
  config->rmsForExposedPixels =sqrt( config->rmsForExposedPixels - 
				     config->meanForExposedPixels*config->meanForExposedPixels );

  return img_histo;
}




TH2S*
MaxCamAndor::createOverscanHisto(TString histoName) {
  return img_histo_overscan;
}


void
MaxCamAndor::calcDaqTime() {
  float exposure, accumulate, kinetics;
  GetAcquisitionTimings( &exposure, &accumulate, &kinetics);
  config->exposureTime=(long int)(exposure)*1000;
  config->daqTime=exposure+accumulate+kinetics;
}


void
MaxCamAndor::setImage() {
    SetImage(config->hbin, config->vbin, 
	     config->ul_x, config->lr_x, config->ul_y, config->lr_y);
    calcDaqTime();
}


TH2S*
MaxCamAndor::createFullHisto(TString histoName) {
  return img_histo_full;
}

int MaxCamAndor::setHBin(long hb) { 
  // Set horizontal binning.
  config->hbin=hb; 
  config->row_width=(config->lr_x-config->ul_x+1)/config->hbin;
  setImage();
  return 0; 
}

int MaxCamAndor::setVBin(long vb) { 
  // Set vertical binning.
  config->vbin=vb; 
  config->img_rows=(config->lr_y-config->ul_y+1)/config->vbin;  
  setImage();
  return 0; 
}

int MaxCamAndor::resetWithFlush() {
  cout << GetName() << "::resetWithFlush not done" << endl;
  return 0;
}

int MaxCamAndor::reset() {
  cout << GetName() << "::reset not done" << endl;
  return 0;
}

int 
MaxCamAndor::setExposureTime(long exptime) { 
  // Set exposure time (in msec).
  cout << GetName() << ": Setting exposure time to " 
       << (float(exptime)*1e-3) << "sec" << endl;
  SetExposureTime(float(exptime)*1e-3); // sec
  calcDaqTime();
  return 0; 
}

int
MaxCamAndor::cancelExposure() {
  // Stop the current exposure.
  AbortAcquisition();
  return 0;
}


double
MaxCamAndor::getTemperature() {
  // Read the actual temperature of the CCD (in deg C).
  int temp;
  GetTemperature(&temp);
  config->CCDTemp=double(temp);
  return config->CCDTemp;
}


int
MaxCamAndor::setTemperature(double temp) {
  // Set CCD temperature (in deg C).

  int minTemp, maxTemp;
  GetTemperatureRange(&minTemp, &maxTemp);
  if (temp<minTemp) temp=minTemp;
  else if (temp>maxTemp) temp=maxTemp;
  cout << GetName() <<": Setting temperature to " << temp << "C" << endl;

  // make sure the cooling is on
  CoolerON();
  config->CCDTemp=temp;
  SetTemperature( int(temp) );
  return 0;
}

int
MaxCamAndor::coolerOff() {
  return CoolerOFF();
}

int
MaxCamAndor::coolerOn() {
  return CoolerON();
}


int
MaxCamAndor::baselineClampOn() {
  return SetBaselineClamp(1);
}

int
MaxCamAndor::baselineClampOff() {
  return SetBaselineClamp(0);
}



double 
MaxCamAndor::getGoalTemperature() { 
  // Return the desired temperature (in deg C).
  float AmbientTemp, CoolerVolts, currTemp, goalTemp;
  GetTemperatureStatus(&currTemp, &goalTemp, &AmbientTemp, &CoolerVolts);
  return goalTemp;
}


int
MaxCamAndor::setDataBits(int bits) {
  // Set depth of ADC converter on camera.

  int num, ivalue;
  GetNumberADChannels(&num);
  int i=0;
  for (; i<num; i++) {
    GetBitDepth(i, &ivalue);
    if (ivalue==bits) break;
  }
  assert(i<num);
  SetADChannel(i);
  cout << GetName() << ": Set bit depth = " << i << "  value="<<ivalue<<endl;

  config->bitDepth=0xFFFF >> (16 - bits);

  return 0;
}

int
MaxCamAndor::getDataBits() {
  // Get bith depth
  cout << GetName() << "::getDataBits not done" << endl;
  return 16;
}

int
MaxCamAndor::getGain() {
  // ADC gain
  int gain;
  if (GetEMCCDGain(&gain) != DRV_SUCCESS) return -1;
  return gain;
}

int
MaxCamAndor::setGain(unsigned short gain) {
  // EM gain

  // for gain==0 use standard CCD amplifier
  if (gain==0) {
    SetOutputAmplifier(1);
    return 0;
  } 

  // for gain>0 use EMCCD amplifier
  SetOutputAmplifier(0);  
  int lowGain, hiGain;
  GetEMGainRange(&lowGain, &hiGain);      
  if (gain<lowGain) gain=lowGain;
  else if (gain>hiGain) gain=hiGain;

  cout << GetName() << ": Setting EM gain to " << gain 
       << "   possible values = 0 and (" << lowGain << ", " << hiGain << ")" <<  endl;

  if (SetEMCCDGain( (int)gain ) != DRV_SUCCESS) return -1;
  return 0;
}


int
MaxCamAndor::setDarkFrame() {
  // Do not open shutter during exposure.
  config->frameType=MaxCamConfig::dark;
  SetShutter(1,2,0,0);
  return 0;
}

int
MaxCamAndor::setNormalFrame() {
  // Open shutter during exposure.
  config->frameType=MaxCamConfig::normal;
  SetShutter(1,0,50,50);
  return 0; 
}

int
MaxCamAndor::setTestFrame() {
  // Flag camera as test.
  cout << GetName() << "::setTestFrame not done" << endl;
  return 0;
}



int
MaxCamAndor::setFanSpeed(int speed) {
  SetFanMode( speed );
  return 0;
}




int 
MaxCamAndor::flushCamera(long repeat) {
  // Flush ccd image...
  cout << GetName() << "::flushCamera not done" << endl;     
  return 0;
}

void
MaxCamAndor::print() { }



bool
MaxCamAndor::imageReady() {
  // Is camera image ready?
  int status;
  GetStatus(&status);
  return status!=DRV_ACQUIRING;
}


int
MaxCamAndor::openShutter() {
  SetShutter(1,1,0,0);
  return 0;
}


int
MaxCamAndor::closeShutter() {
  SetShutter(1,2,0,0);
  return 0;
}


int
MaxCamAndor::setDigitizeOverscan(bool digitize) {
  cout << GetName() << "::setDigitizeOverscan not done" << endl;     
  return 0;
}


float
MaxCamAndor::getReadoutTime() {
  float rot=0;
  GetReadOutTime(&rot);
  return rot;
}


int
MaxCamAndor::setPreAmpGain(int igain) {
  int num;
  GetNumberPreAmpGains(&num);
  assert (igain>-1 && igain <num);
  cout << GetName() << ": Set gain index to " << igain 
       << " range=(0, " << num << ")" << endl;
  float value;
  GetPreAmpGain(igain, &value);
  cout << GetName() << ": Set gain value to " << value << endl; 
  return SetPreAmpGain(igain);
}
