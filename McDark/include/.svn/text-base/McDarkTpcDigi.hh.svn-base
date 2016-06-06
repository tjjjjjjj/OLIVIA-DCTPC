#ifndef McDarkTpcDigi_h
#define McDarkTpcDigi_h 1

#include "G4VDigi.hh"
#include "G4TDigiCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"



class McDarkTpcDigi : public G4VDigi {
    
  public:

      McDarkTpcDigi();
      ~McDarkTpcDigi();
      McDarkTpcDigi(const McDarkTpcDigi&);
      const McDarkTpcDigi& operator=(const McDarkTpcDigi&);
      int operator==(const McDarkTpcDigi&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      void Draw();
      void Print();

  public:

    void setModuleID(G4int m) { _moduleID=m; }
    G4int getModuleID() { return _moduleID; }

    void setChannelID(G4int c) { _channelID=c; }
    G4int getChannelID() { return _channelID; }

    void setWeight(G4double w) { _weight=w; }
    G4double getWeight() { return _weight; }
    
    void setHitID(G4int index) { _hitIndex=index; }
    G4int getHitID() { return _hitIndex; }
    
    void setDigiID(G4int index) { _digiIndex=index; }
    G4int getDigiID() { return _digiIndex; }

    void setTrackID(G4int index) { _trackIndex=index; }
    G4int getTrackID() { return _trackIndex; }
    
  private:
    
    G4int _moduleID;  // module number, e.g. camera number or PMT number for this TPC
    G4int _channelID; // pixel-bin in CCD, time-bin in charge/PMT readout
    G4double _weight; // ADC value
    
    G4int _hitIndex;  // index of hit causing this digi 
    G4int _digiIndex; // index of this digi
    G4int _trackIndex; // index of track causing this digi
};


typedef G4TDigiCollection<McDarkTpcDigi> McDarkTpcDigiCollection;

extern G4Allocator<McDarkTpcDigi> McDarkTpcDigiAllocator;


inline void*
McDarkTpcDigi::operator new(size_t) {
  void* aDigi;
  aDigi = (void*) McDarkTpcDigiAllocator.MallocSingle();
  return aDigi;
}


inline void
McDarkTpcDigi::operator delete(void* aDigi)  {
    McDarkTpcDigiAllocator.FreeSingle((McDarkTpcDigi*) aDigi);
}

#endif

