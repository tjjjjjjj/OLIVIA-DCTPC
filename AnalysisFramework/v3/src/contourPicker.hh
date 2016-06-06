#ifndef CONTOUR_PICKER_HH
#define CONTOUR_PICKER_HH
#include "TH2.h"
#define CNTRDBG 0

static double cntr_defaultLevels[] = {0, 0.68, 0.95, 0.99, 1} ;

class ContourPicker
{

  public:
    

    ContourPicker(int n = 5, 
        double * levels = cntr_defaultLevels , 
        double tol = 1); 
    ~ContourPicker(); 

    void setLevels(int n, double * levels) { _n = n; _levels = levels; }
    void setTol(double tol){ _tol = tol;} 
    double * getContours(const TH2 * hist); 
    void drawHistogramContours(TH2 * hist, const bool same); 
    unsigned int nLevels(){ return _n;} 

    TH2F * makeContourHist(const TH2 * hist, char * name); 
    
  private:
    double _tol; 
    int _n; 
    double * _levels; 

};

#endif
