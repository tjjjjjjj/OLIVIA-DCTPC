#ifndef MAXCAM_SEGMENT_HH
#define MAXCAM_SEGMENT_HH

#include "TObject.h"

class MaxCamSegment : public TObject {

public:

    // Ctors
    
    MaxCamSegment();
    
    MaxCamSegment(const MaxCamSegment &other);
    
    virtual ~MaxCamSegment() {}

    void setEnergy(Float_t E, Float_t EErr) { _energy=E; _energyErr=EErr; }
    Float_t E() { return _energy; }
    Float_t EErr() { return _energyErr; }

    void setNPixels(Int_t n) { _npixels=n; }
    Int_t  n() { return _npixels; }

    void setNSegments(Int_t n) { _nsegments=n; }
    Int_t  nseg() { return _nsegments; }
    
    void setCorrelation(Float_t c) { _correlation=c; }
    Float_t correlation() { return _correlation; }
    
    void setLength(Float_t E, Float_t EErr) { _length=E; _lengthErr=EErr; }
    Float_t L() { return _length; }
    Float_t LErr() { return _lengthErr; }

    void setWidth(Float_t w, Float_t we) { _width=w; _widthErr=we; }
    Float_t W() { return _width; }
    Float_t WErr() { return _widthErr; }

    void setAmplitude(Float_t a, Float_t ae) { _amplitude=a; _amplitudeErr=ae;    }
    Float_t A() { return _amplitude; }
    Float_t AErr() { return _amplitudeErr; }

    void setBackground(Float_t b, Float_t be) { _bgLevel=b; _bgLevelErr=be; }
    Float_t B() { return _bgLevel; }
    Float_t BErr() { return _bgLevelErr; }
    
    void setMean(Float_t m, Float_t me) { _mean=m; _meanErr=me; }
    Float_t Mean() { return _mean; }
    Float_t MeanErr() { return _meanErr; }

    void print();

    void setDischarge(Bool_t b=true) { _isDischarge=b; }
    Bool_t isDischarge() { return _isDischarge; }
    
    void setBorderline(Bool_t b=true) { _isBorderline=b; }
    Bool_t isBorderline() { return _isBorderline; }

    void setSkewness(Float_t s) {  _skewness=s; }
    Float_t skewness() { return _skewness; }

    void reset() {
        _energy=_energyErr=_correlation=_length=_lengthErr=_width=_widthErr=0;
        _amplitude=_amplitudeErr=_bgLevel=_bgLevelErr=_mean=_meanErr;
        _npixels=_nsegments=0;
        _skewness=-999;
        _isDischarge=_isBorderline=false;
        _Ix=_Iy=_cosRecoil=0;
    }

    void setIxy(float ix, float iy) { _Ix=ix; _Iy=iy; }
    Float_t Ix() { return _Ix; }
    Float_t Iy() { return _Iy; }

    void setCosRecoil( float c) { _cosRecoil=c; }
    Float_t cosRecoil() { return _cosRecoil; }
    
private:

    Float_t _energy;
    Float_t _energyErr;
    Int_t   _npixels;
    Int_t   _nsegments;
    Float_t _correlation;
    Float_t _length;
    Float_t _lengthErr;
    Float_t _width;
    Float_t _widthErr;
    Float_t _amplitude;
    Float_t _amplitudeErr;
    Float_t _bgLevel;
    Float_t _bgLevelErr;
    Float_t _mean;
    Float_t _meanErr;
    Bool_t  _isDischarge;
    Bool_t  _isBorderline;
    Float_t  _skewness;
    Float_t _Ix;
    Float_t _Iy;
    Float_t _cosRecoil;
    
    ClassDef(MaxCamSegment,1)
};

#endif

