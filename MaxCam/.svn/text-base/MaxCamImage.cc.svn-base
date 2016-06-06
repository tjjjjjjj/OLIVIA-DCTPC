#include "MaxCamImage.hh"

ClassImp(MaxCamImage)



    
MaxCamImage::MaxCamImage(TH2F *image, unsigned int camID) : TH2F(*image) {

        _cameraID=camID;
    
}



MaxCamImage::MaxCamImage(MaxCamImage &image) : TH2F(), 
    _cameraID(image.getCameraID())
    
{        
    ((TH2F&)image).Copy(*this);
  
}


MaxCamImage::~MaxCamImage() {
}
