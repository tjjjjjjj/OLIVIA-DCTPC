#include "McDarkTrack.hh"
#include <iostream>
using std::endl;
using std::cout;
#include <iomanip>
using std::setw;

ClassImp(McDarkTrack)

//
    McDarkTrack::McDarkTrack() : TObject() {;}

McDarkTrack::McDarkTrack(const McDarkTrack &other) : TObject(other) {

    _index=other._index;; 
    _motherIndex=other._motherIndex;
    
    setP4( other._P4 );

    _quenchedEnergy=other._quenchedEnergy;
    
    _PDGIndex=other._PDGIndex;

    _begPoint=other._begPoint;
    _endPoint=other._endPoint;

}

    

void
McDarkTrack::setP4(TLorentzVector v) {
    _P4=v;
    _Ekine=v.E()-v.M();
}

void
McDarkTrack::setP4(Double_t px, Double_t py, Double_t pz, Double_t E) {
    TLorentzVector v(px,py,pz,E);
    setP4(v);
}


void
McDarkTrack::print() {
// Print MC track info.

    cout << setw(6)  << index() << " | "
         << setw(10) << PDGIndex() << " | "                   
         << setw(6)  << motherIndex() << " | "
         << setw(12) << Ekine() << " | "
         << setw(12) << (begPoint()-endPoint()).Mag() << " | "
         << setw(12) << quenchedEnergy()
         << endl; 
}
