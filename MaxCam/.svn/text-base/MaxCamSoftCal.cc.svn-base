#include "MaxCamSoftCal.hh"
#include <stdlib.h>
#include <stdio.h>
//#include <malloc.h>

static const unsigned int TO_PHYS = 0; 
static const unsigned int FROM_PHYS = 1; 

#ifdef DM_DAQ
MaxCamSoftCal::MaxCamSoftCal(comedi_t * dev, char * path)
{
  if (path==NULL) calibpath = comedi_get_default_calibration_path(dev); 
  else calibpath = path; 
  calib = comedi_parse_calibration_file(calibpath); 
  comedi_close(dev); 
}

MaxCamSoftCal::~MaxCamSoftCal()
{
  free(calibpath); 
//  comedi_cleanup_calibration_file(calib); 
}


double MaxCamSoftCal::toPhysical(unsigned int value, int subdev, int chan, int range)
{
    comedi_get_softcal_converter(subdev,chan,range, COMEDI_TO_PHYSICAL, calib, &poly); 
    return comedi_to_physical(value,&poly);
}


unsigned int MaxCamSoftCal::fromPhysical(double value, int subdev, int chan, int range)
{
    comedi_get_softcal_converter(subdev,chan,range, COMEDI_FROM_PHYSICAL, calib, &poly); 
    return comedi_from_physical(value,&poly);
}

#endif 
