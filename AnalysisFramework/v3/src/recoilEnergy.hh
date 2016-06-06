#ifndef RECOIL_ENERGY_HH
#define RECOIL_ENERGY_HH

#include <string>
#include "../../../MaxCam/DmtpcSkimEvent.hh"

/* recoilEnergy.hh 

   Utility function for getting the recoil energy from the 
   observerd energy. 

*/


namespace RecoilEnergy
{
    double visibleToRecoil(const double E);
    double getRecoilEnergy(double E, int cam);
    double getRecoilEnergy(double E, int cam, double x, double y);
    double getRecoilEnergy(DmtpcSkimEvent * ev, int cam, int track,
                           bool useGainMap = true); 
    double getVisibleEnergy(DmtpcSkimEvent * ev, int cam = 0, int track = 0, bool useGainMap = true); 
    void setGainMapDirectory(std::string dir); 


};

#endif

