#include "McDarkTpcDigi.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include <iomanip>

G4Allocator<McDarkTpcDigi> McDarkTpcDigiAllocator;


McDarkTpcDigi::McDarkTpcDigi() : 
    _moduleID(-1),
    _channelID(-1),
    _weight(0.0),
    _hitIndex(-1),
    _digiIndex(-1),
    _trackIndex(-1)
{;}


McDarkTpcDigi::~McDarkTpcDigi() {;}


McDarkTpcDigi::McDarkTpcDigi(const McDarkTpcDigi& right) : G4VDigi(right) {
    _moduleID = right._moduleID;
    _channelID = right._channelID;
    _weight = right._weight;
    _hitIndex = right._hitIndex;
    _digiIndex = right._digiIndex;
    _trackIndex = right._trackIndex;
}


const McDarkTpcDigi&
McDarkTpcDigi::operator=(const McDarkTpcDigi& right) {
    _moduleID = right._moduleID;
    _channelID = right._channelID;
    _weight = right._weight;
    _hitIndex = right._hitIndex;
    _digiIndex = right._digiIndex;
    _trackIndex = right._trackIndex;
    return *this;
}


int
McDarkTpcDigi::operator==(const McDarkTpcDigi& right) const {
  return (this==&right) ? 1 : 0;
}


void
McDarkTpcDigi::Draw() {;}


void
McDarkTpcDigi::Print() {;}



