#include "McDarkTpcHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include <iomanip>


G4Allocator<McDarkTpcHit> McDarkTpcHitAllocator;


McDarkTpcHit::McDarkTpcHit() :
    _edep(0),
    _pos(0., 0., 0.),
    _time(0),
    _trackIndex(-1),
    _hitIndex(-1),
    _pdgIndex(-1),
    _qFactor(1)
{}


McDarkTpcHit::~McDarkTpcHit() {;}


McDarkTpcHit::McDarkTpcHit(const McDarkTpcHit& right)
  : G4VHit(right) {
  _edep = right._edep;
  _pos = right._pos;
  _time = right._time;
  _trackIndex = right._trackIndex;
  _hitIndex = right._hitIndex;
  _pdgIndex = right._pdgIndex;
  _qFactor = right._qFactor;
}


const McDarkTpcHit&
McDarkTpcHit::operator=(const McDarkTpcHit& right) {
  _edep = right._edep;
  _pos = right._pos;
  _time = right._time;
  _trackIndex = right._trackIndex;
  _hitIndex = right._hitIndex;
  _pdgIndex = right._pdgIndex;
  _qFactor = right._qFactor;
  return *this;
}


int
McDarkTpcHit::operator==(const McDarkTpcHit& right) const {
  return (this==&right) ? 1 : 0;
}


void
McDarkTpcHit::Draw() {;}


void
McDarkTpcHit::Print() {
  G4cout << "      TPC hit: " << std::setw(5) << G4BestUnit(_edep,"Energy") 
	 << ", at " << G4BestUnit(_pos,"Length") << G4endl;
}



