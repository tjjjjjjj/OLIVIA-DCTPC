#ifndef MAX_CAM_IMAGE
#define MAX_CAM_IMAGE

#include "TH2F.h"


class MaxCamImage : public TH2F {

public:
    MaxCamImage(TH2F* image, unsigned int cameraID);

  MaxCamImage() {}

    MaxCamImage(MaxCamImage &image);

    ~MaxCamImage();
    
    unsigned int getCameraID() { return _cameraID; }

    
private:
    
    unsigned int _cameraID;

    ClassDef(MaxCamImage, 1)
};



#endif
