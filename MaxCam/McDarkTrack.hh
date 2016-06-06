#ifndef MCDARK_TRACK_HH
#define MCDARK_TRACK_HH

#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector3.h"


class McDarkTrack : public TObject {

public:

    // Ctors
    McDarkTrack();    
    
    McDarkTrack(const McDarkTrack &other);
    
    virtual ~McDarkTrack() {};
    
    UInt_t index() { return _index; }
    void setIndex(UInt_t idx) { _index=idx; }

    UInt_t motherIndex() { return _motherIndex; }
    void setMotherIndex(UInt_t idx) { _motherIndex=idx; }

    ULong_t PDGIndex() { return _PDGIndex; }
    void setPDGIndex(ULong_t idx) { _PDGIndex=idx; }
    
    TLorentzVector P4() { return _P4; }
    void setP4(TLorentzVector v);
    void setP4(Double_t px, Double_t py, Double_t pz, Double_t E);

    void setQuenchedEnergy( double qEn ) { _quenchedEnergy=qEn; }
    double quenchedEnergy() { return _quenchedEnergy; }
    
    Double_t mass() { return _P4.M(); }
    Double_t Ekine() { return _Ekine; }

    TVector3 begPoint() { return _begPoint; }
    void setBegPoint(TVector3 point) { _begPoint=point; }

    TVector3 endPoint() { return _endPoint; }
    void setEndPoint(TVector3 point) { _endPoint=point; }
    
    void print();
  
private:

    UInt_t   _index; // index in Geant's track table
    UInt_t   _motherIndex; // points to mother's index

    TLorentzVector _P4; // 4-momentum

    double _quenchedEnergy; // ionization energy in senstive volume
    
    ULong_t _PDGIndex; // PDG particle index
    
    Double_t _Ekine; // precomputed kinetic energy (from 4-vector)

    TVector3 _begPoint; // 1st point on the track
    TVector3 _endPoint; // last point on the track

    ClassDef(McDarkTrack,1)
};

#endif

