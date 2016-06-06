#ifndef SIM_RINGS_HH
#define SIM_RINGS_HH

#include <vector> 
#include "TRandom3.h"
#include "TObject.h"

class SimRings : public TObject
{

  public:

    /**
     * all coords in cm 
     *
     * @param n number of rings
     * @param x0  x coordinate of center of rings
     * @param y0  y coordinate of center of rings
     * @param rmin  minimum radius of rings
     * @param rmax  maximum radius of rings
     * @param zmin the z coordinate of the bottom of the bottom ring 
     * @param zmax the z coordinate of the top of the topmost ring 
     * @param dz  the thickness of the ring 
     */

    SimRings(int n, double x0, double y0, double rmin, double rmax, double zmin, double zmax, double dz, double ztop, double surftop, double zbottom, double surfbottom); 
    SimRings(){}; 
    void generateDecay(double *pos, double * dir); 

  private: 
    
    std::vector<double> zs; 
    double zgap; 
    double dz; 
    double x0;
    double y0; 
    double rin; 
    double rout; 
    double dr; 
    double ztop; 
    double zbottom; 
    double surftop;
    double surfbottom;
    int nrings; 
    double ring_fraction; 
    double top_fraction; 
    double bottom_fraction; 
    TRandom3 rand; 

    ClassDef(SimRings,1); 

};

#endif
