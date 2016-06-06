#ifndef McDarkTpcHit_h
#define McDarkTpcHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"



class McDarkTpcHit : public G4VHit
{
  public:

      McDarkTpcHit();
      ~McDarkTpcHit();
      McDarkTpcHit(const McDarkTpcHit&);
      const McDarkTpcHit& operator=(const McDarkTpcHit&);
      int operator==(const McDarkTpcHit&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      void Draw();
      void Print();

  public:
  
    void setDepositEnergy(G4double de) { _edep = de; }
    void setDepositPosition(G4ThreeVector xyz) { _pos = xyz; }
    void setDepositTime(G4double t) { _time = t; }
    G4double getDepositEnergy() { return _edep; } 
    G4ThreeVector getDepositPosition() { return _pos; }
    G4double getDepositTime() { return _time; }

    void setTrackID(G4int index) { _trackIndex=index; }
    G4int getTrackID() { return _trackIndex; }
    
    void setHitID(G4int index) { _hitIndex=index; }
    G4int getHitID() { return _hitIndex; }

    void setPDGIndex(G4int pdg) { _pdgIndex=pdg; }
    G4int PDGIndex() { return _pdgIndex; }

    void setQuenchingFactor(G4double qf) { _qFactor=qf; }
    G4double getQuenchingFactor() { return _qFactor; }
    
  private:
    
    G4double      _edep; // deposited energy
    G4ThreeVector _pos; // mean position of segment
    G4double      _time; // mean time of segment

    G4int        _trackIndex; // index of track causing ionization in G4Tracks table
    G4int        _hitIndex; // index of this hit
    G4int        _pdgIndex; // PDG index of track causing ionization
    G4double _qFactor; // quenching factor = ionization loss/total energy loss
};


typedef G4THitsCollection<McDarkTpcHit> McDarkTpcHitCollection;

extern G4Allocator<McDarkTpcHit> McDarkTpcHitAllocator;


inline void* McDarkTpcHit::operator new(size_t) {
  void* aHit;
  aHit = (void*) McDarkTpcHitAllocator.MallocSingle();
  return aHit;
}


inline void
McDarkTpcHit::operator delete(void* aHit) {
  McDarkTpcHitAllocator.FreeSingle((McDarkTpcHit*) aHit);
}

#endif

