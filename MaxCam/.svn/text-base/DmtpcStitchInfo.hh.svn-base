#ifndef DMTPC_STITCH_INFO_HH
#define DMTPC_STITCH_INFO_HH

#include "TNamed.h" 
#include <vector>
#include "TH2F.h"
#include "DmtpcLensCorrection.hh" 
#include "TH2D.h"

class DmtpcStitchInfo : public TNamed
{

  public:

    DmtpcStitchInfo (const char * name, unsigned n, const DmtpcLensCorrection * lens,  const double * xmins, const double * ymins, const double * xmaxs, const double * ymaxs,
                     const int * nbinsx, const int * nbinsy,  const double * xorigin = 0, const double * yorigin = 0, const double * rotation = 0, 
                     const double * scale = 0, const double * weight = 0, unsigned binning = 0); 

    DmtpcStitchInfo() {weight_sum = 0;}

    virtual ~DmtpcStitchInfo(); 

    unsigned nimages;
    std::vector<std::vector<std::vector<double> > > corners;
    std::vector<double> xmins;
    std::vector<double> xmaxs;
    std::vector<double> ymins;
    std::vector<double> ymaxs;
    std::vector<int> nbinsx;
    std::vector<int> nbinsy;
    std::vector<double> xorigin;
    std::vector<double> yorigin;
    std::vector<double> rotation;
    std::vector<double> sintheta;
    std::vector<double> costheta;
    std::vector<double> scale; 
    std::vector<double> weight; 
    std::vector<double>xwidths;
    std::vector<double>ywidths;
    std::vector<TH2F*> frac;
    double stitched_xmax; 
    double stitched_ymax; 
    double stitched_xmin; 
    double stitched_ymin; 
    unsigned stitched_nbinsx; 
    unsigned stitched_nbinsy; 
    TH2F * weight_sum; 

    ClassDef(DmtpcStitchInfo,4); 
};


#endif
