#include "WaveformVector.hh"
#include "SkimWaveform.hh"
#include "DmtpcPulse.hh"
#include <vector>
#include <string>
#include "TObject.h"
using std::vector;

ClassImp(WaveformVector)

WaveformVector::WaveformVector(int n,int bid, int cid, const char* nm) :
name(nm), boardID(bid), chanID(cid), N(n), wf(n){;}

WaveformVector::WaveformVector(const WaveformVector& w) : TObject(w)
{
  name = w.name;
  boardID = w.boardID;
  chanID = w.chanID;
  N = w.N;
  wf = w.wf;
}


void
WaveformVector::add(const SkimWaveform& w)
{
  N++;
  wf.push_back(w);

}

void
WaveformVector::insert(int i, const SkimWaveform& w)
{
  N++;
  wf.insert(wf.begin()+i,w);
}

void
WaveformVector::swap(int i, const SkimWaveform& w)
{
  wf[i] = w;
}

void
WaveformVector::rm(int i)
{
  N--;
  wf.erase(wf.begin()+i);
}

void
WaveformVector::clearWf()
{
  N=0;
  wf.clear();
}

void
WaveformVector::resize(int i)
{
  N=i;
  wf.resize(i);
}

WaveformVector&
WaveformVector::operator=(const WaveformVector& v)
{
  if (&v==this) return *this;
  N = v.N;
  name=v.name;
  wf=v.wf;
  return *this;
}


const SkimWaveform& 
WaveformVector::at(int i) const{return wf[i];}

const SkimWaveform& 
WaveformVector::operator()(int i) const{return wf[i];}

const SkimWaveform& 
WaveformVector::operator[](int i) const{return wf[i];}

SkimWaveform& 
WaveformVector::at(int i) {return wf[i];}

SkimWaveform& 
WaveformVector::operator()(int i) {return wf[i];}

SkimWaveform& 
WaveformVector::operator[](int i) {return wf[i];}

const DmtpcPulse& 
WaveformVector::at(int i, int j) const{return wf[i][j];}

const DmtpcPulse& 
WaveformVector::operator()(int i, int j) const{return wf[i][j];}

DmtpcPulse& 
WaveformVector::at(int i, int j){return wf[i][j];}

DmtpcPulse& 
WaveformVector::operator()(int i, int j){return wf[i][j];}

