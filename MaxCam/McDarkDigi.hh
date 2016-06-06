#ifndef MCDARK_DIGI_HH
#define MCDARK_DIGI_HH

#include "TObject.h"
#include "TVector3.h"


class McDarkDigi : public TObject {

public:

    // Ctors
    McDarkDigi();    
    
    McDarkDigi(const McDarkDigi &other);
    
    virtual ~McDarkDigi() {};
    
    void setModuleID(int m) { _moduleID=m; }
    int getModuleID() { return _moduleID; }

    void setChannelID(int c) { _channelID=c; }
    int getChannelID() { return _channelID; }

    void setWeight(double w) { _weight=w; }
    double getWeight() { return _weight; }    

    void setTrackID(int index) { _trackIndex=index; }
    int getTrackID() { return _trackIndex; }
    
  private:
    
    int _moduleID;  // module number, e.g. camera number or PMT number for this TPC
    int _channelID; // pixel-bin in CCD, time-bin in charge/PMT readout
    double _weight; // ADC value
    
    int _trackIndex; // index of track causing this digi

    ClassDef(McDarkDigi,1)
};

#endif

