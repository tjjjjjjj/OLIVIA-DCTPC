#include "TH2.h"
#include "SimCamera.hh"
#include "TObject.h"
#include "TF1.h" 
#include "TFile.h"
#include "TRandom3.h"
#include <iostream>
#include <math.h>
//________________________________
/*Begin_Html
<center><h2>The SimCamera Class </h2></center>
The SimCamera class is designed to hold the various parameters of the 
CCD camera used in the detector.  A separate class is used so that in
future versions of the simulation, several instances of the class may
be added to a TObjectArray or TClonesArray to allow for generation of
multiple images, each having its own properties.
End_Html*/
//________________________________

ClassImp(SimCamera)

SimCamera::SimCamera()
{
  std::cout << "default constructor" << std::endl;
  //<<<<<Default Constructor>>>>>
  fSerialNumber = "";
  fCameraNumber = 0;
  fRnd = new TRandom3(0);
  fX = 0;
  fY = 0;
  fXBins = 256;
  fYBins = fXBins;
  fPixelsPerBin = 4;
  fXWidth = 0.179*fXBins*fPixelsPerBin;
  fYWidth = fXWidth;
  fNoiseFactor = 1;
  fEmGain = 1;
  fGain = 9.3;
  fGainMap = NULL; 
  fDarkCurrent = 0;
  fActiveRegion = 0; 
  fActiveRegionUnits = "mm";
  fReadNoise = 7.3;
  fBias = 0;

}
SimCamera::SimCamera(double x, double y,  int xbin, int ybin,int pixperbin,
		     double xwidth,double ywidth,int camnum)
{
  std::cout << "alternate constructor" << std::endl;
  //<<<<<Alternate Constructor>>>>
  fSerialNumber = "";
  fRnd = new TRandom3(0);
  fCameraNumber = camnum;
  fX = x;
  fY = y;
  fXBins = xbin;
  fYBins = ybin;
  fPixelsPerBin = pixperbin;
  fXWidth = xwidth;
  fYWidth = ywidth;
  fNoiseFactor = 1;
  fEmGain = 1;
  fDarkCurrent = 0;
  fReadNoise = 7;
  fGainMap = NULL; 
  fGain = 9.3;
  fActiveRegion = 0; 
  fActiveRegionUnits = "mm";
  setImage();

}

SimCamera::~SimCamera()
{
  //<<<<<Destructor>>>>>
  delete fRnd;
  if (fccdImage) delete fccdImage;
  if (fGainMap) delete fGainMap; 
  if (fActiveRegion) delete fActiveRegion; 
}

void
SimCamera::setPosition(double x, double y)
{
  //Sets the x-y position of the image center
  fX = x;
  fY = y;
}

void
SimCamera::setWidths(double xwidth, double ywidth)
{
  //Sets the x and y widths in mm of the image 
  fXWidth = xwidth;
  fYWidth = ywidth;
}
void
SimCamera::setBins(int xbin,int ybin)
{
  //Sets the number of x and y bins of the image
  fXBins = xbin;
  fYBins = ybin;
}

void
SimCamera::setGain_Landau(double landau_1,double landau_2)
{
  fGainLandau_1 = landau_1;
  fGainLandau_2 = landau_2;
}  


void
SimCamera::setRadialParams(double radial_1,double radial_2,double radial_3)
{
  fRadialParams_1 = radial_1;
  fRadialParams_2 = radial_2;
  fRadialParams_3 = radial_3;
}

void
SimCamera::setPhotonsADU(double photonsadu)
{
  fPhotonsADU = photonsadu;
}

void
SimCamera::setImage()
{
  //Sets the image according to the current parameters.
  //Uses fCameraNumber in the name to avoid having multiple
  //images with the same name.
  TString name = "fccd_";
  name += fCameraNumber;
  fccdImage = new TH2F(name,name,fXBins,0,fXBins*fPixelsPerBin,
		       fYBins,0,fYBins*fPixelsPerBin);
}

void
SimCamera::resetImage()
{
  //Resets the CCD image
  fccdImage->Reset();
}

void
SimCamera::applyNoise()
{
  //Applies bias and Gaussian read noise to image
  for(int i = 1; i<=fXBins;i++){
    for(int j = 1; j<= fYBins;j++){
      double signal = fccdImage->GetBinContent(i,j);
      if (fDarkCurrent > 0){
        double dc = fRnd->Poisson(fDarkCurrent);
        signal += fRnd->Poisson(dc*fEmGain);
      }
      //add 1000 here because the eventual stored TH2S has trouble with negative values
      fccdImage->SetBinContent(i,j,fRnd->Gaus(fBias+signal,fReadNoise));
      // fccdImage->SetBinContent(i,j,fRnd->Gaus(fBias+signal,fReadNoise));
      //std::cout<<i<<" "<<j<<" "<<fccdImage->GetBinContent(i,j)<<std::endl;
    }
  }

  //for(int i=0;i<256;i++)
    //for(int j=0;j<256;j++)
      // {
	//cout<<i<<" "<<j<<" "<<fccdImage->GetBinContent(i,j)<<endl;

	//}

}

void
SimCamera::emptyImage()
{
  //Creates an empty image with only dark current and read noise
  resetImage();
  double sigma = sqrt( fDarkCurrent*fNoiseFactor*fNoiseFactor + 
		       fReadNoise*fReadNoise);
  for(int i = 1; i<=fXBins;i++){
    for(int j = 1; j<= fYBins;j++){
      fccdImage->SetBinContent(i,j,fRnd->Gaus(fBias,sigma));
      //std::cout<<i<<" "<<j<<" "<<fRnd->Gaus(fBias,sigma)<<std::endl;
    }
  }
}

double 
SimCamera::getBinPositionX(int i, TString opt)
{
  //Gets x position in mm of a given x bin i.
  double bin = i -0.5;//Using convention of taking midpoint as bin position
  if (opt == "lo") bin = i - 1;
  else if (opt == "hi") bin = i; 
  double position = bin / fXBins * fXWidth + fX - fXWidth / 2;
  return position;
}

double 
SimCamera::getBinPositionY(int i, TString opt)
{
  //Gets y position in mm of a given y bin i
  double bin = i-0.5;
  if (opt == "lo") bin = i - 1;
  else if (opt == "hi") bin = i; 
  double position = bin / fYBins * fYWidth + fY - fYWidth/2;
  return position;
}

double
SimCamera::getPixelPositionX(int i, TString opt)
{
  //=====Returns the position of the given x pixel: Option for high or low edge
  double pixel = i - 0.5;
  if (opt == "lo") pixel = i - 1;
  else if (opt == "hi") pixel = i;
  double position = pixel / (fXBins*fPixelsPerBin) * fXWidth + fX - fXWidth/2;
  return position;
}
double
SimCamera::getPixelPositionY(int i, TString opt)
{
  //=====Returns the position of the given y pixel: Option for high or low edge
  double pixel = i - 0.5;
  if (opt == "hi") pixel = i - 1;
  else if (opt == "lo") pixel = i;
  double position = pixel / (fYBins*fPixelsPerBin) * fYWidth + fY - fYWidth/2;
  return position;
}

double 
SimCamera::getPixelX(double x)
{
  //======Returns the value of the x-pixel with the given x-coordinate==
  double pos = x + fXWidth/2 - fX;//distance from far left edge
  double pixel = pos / fXWidth * fXBins*fPixelsPerBin;
  return pixel;
}

double 
SimCamera::getPixelY(double y)
{
  //======Returns value of y-pixel with the given y-coordinate===========
  double pos = y + fYWidth / 2 - fY; //distance from lower edge
  double pixel = pos / fYWidth * fYBins * fPixelsPerBin;
  return pixel;
}

bool
SimCamera::isInImage(double x, double y)
{
  //======Test if (x,y) is in the image described by this object==========
  bool inImage = false;
  if (x >= fX - fXWidth/2. && y > fY - fYWidth/2. && x < fX+fXWidth/2. && y < fY+ fYWidth/2.) inImage = true;

  return inImage;
}

void SimCamera::normalizeGainMap()
{
	if (fGain==0 || fGainMap==NULL) return; 
//  std::cout <<fGainMap->Integral() << std::endl; 

  if (fGainMap->Integral() != fGainMap->GetNbinsX() * fGainMap->GetNbinsY())
  {
    double min = fGainMap->GetMinimum(); 

    if (min < 0)
    {
      for (int i = 1; i <= fGainMap->GetNbinsX(); i++)
      {
        for (int j = 1; j <= fGainMap->GetNbinsY(); j++)
        {
          fGainMap->SetBinContent(i,j,fGainMap->GetBinContent(i,j) - min);
	  //std::cout<<fGainMap->GetBinContent(i,j) - min<<std::endl;        
        }
      }
    }

    fGainMap->Scale(fGainMap->GetNbinsX() * fGainMap->GetNbinsY() /fGainMap->Integral());

    TString fname = "normalized_"; 
    fname += fGainMap->GetName(); 
    fname +=".root"; 
    TFile f(fname,"RECREATE"); 
    ((TH2F*) (fGainMap->Clone("normalized")))->Write();
    f.Close(); 
  }

	fGainMap->Scale(fGain); 
}

// input is vixels
double SimCamera::getGain(int u, int v)
{
  if (fGainMap==NULL && fActiveRegion == NULL) return getGain(); 

  if (fActiveRegion)
  {
    
    bool is_inside_active_region=true;

    // active region is specified in units of camera pixels
    // have to convert vixels into hardware pixels
    if(fActiveRegionUnits.Contains("pix")){
     
      // 
      is_inside_active_region=fActiveRegion->IsInside(u*fPixelsPerBin,v*fPixelsPerBin);

    } else { // default is mm

      double xcal = fXWidth / (fXBins); 
      double ycal = fYWidth / (fYBins); 
      double x = fX + (u - fXBins/2) * xcal; 
      double y = fY + (v - fYBins/2) * ycal; 

      is_inside_active_region=fActiveRegion->IsInside(x,y);

    }
    
    if(!is_inside_active_region)
    {
        return 0; 
    }
  }

  if (fGainMap)
  {
    double g =  fGainMap->GetBinContent(u,v); 
    return g; 
  }

  return getGain(); 
}

TH2S* SimCamera::getRawCCDImage()
{
  TString name = "ccd_";
  name += fCameraNumber;
  TH2S* rawIm = new TH2S(name,name,fXBins,0,fXBins*fPixelsPerBin,
		         fYBins,0,fYBins*fPixelsPerBin);

  // cout<<fXBins<<" 0 "<<fXBins*fPixelsPerBin<<" "<<fYBins<<" 0 "<<fYBins*fPixelsPerBin<<endl;

  for (int i = 1; i<=fXBins;i++){
    for (int j = 1; j<=fYBins; j++){
      rawIm->SetBinContent(i,j,(short)fccdImage->GetBinContent(i,j));
      //std::cout<<i<<" "<<j<<" "<<(short)fccdImage->GetBinContent(i,j)<<std::endl;
    }
  }

  return rawIm;
}
