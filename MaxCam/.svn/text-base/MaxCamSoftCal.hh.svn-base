#ifndef MAX_CAM_SOFT_CAL_HH
#define MAX_CAM_SOFT_CAL_HH
#ifdef DM_DAQ 
#include <comedilib.h>
#endif


class MaxCamSoftCal
{

#ifdef DM_DAQ 
  public:
    MaxCamSoftCal(comedi_t * dev, char * path = NULL);
    ~MaxCamSoftCal(); 

    double toPhysical(unsigned int value, int subdev, int chan, int range);
    unsigned int fromPhysical(double value, int subdev, int chan, int range); 

  private:  
    char * calibpath; 
    comedi_calibration_t * calib; 
    comedi_polynomial_t poly; 
  #endif
};

#endif
