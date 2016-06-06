#ifndef MAXCAM_VERTEX_HH
#define MAXCAM_VERTEX_HH

#include "TObject.h"
#include "TF1.h"
#include <vector>
using std::vector;



class MaxCamVertex : public TObject {

public:

    // Ctors
    
    MaxCamVertex();
    
    MaxCamVertex(const MaxCamVertex &other);

    virtual ~MaxCamVertex() {}

    void addTrack2D(TF1* f) {
        _trackList2D.push_back( f );
    }

    double calcChiSquare2D();
    
    static double distanceSquare2D(double x, double y, TF1* f);

    int fit2D();

    double getXfit() { return _xfit; }
    double getYfit() { return _yfit; }
    double getXerr() { return _xerr; }
    double getYerr() { return _yerr; }
    
private:

    double _xfit, _yfit; // fitted coordinates
    double _xerr, _yerr; // fitted errors
    double _chi2; // chi2
    
    vector<TF1*> _trackList2D; // list of tracks

    ClassDef(MaxCamVertex,1)
};

#endif

