#include "CspWfVector.hh"
#include "CspPulse.hh"
#include "CspWaveform.hh"
#include <vector>
#include <string>
#include "TObject.h"
using std::vector;


ClassImp(CspWfVector)

CspWfVector::CspWfVector(int n, int bid, int cid,const char* nm) :
name(nm), boardID(bid), chanID(cid), N(n), wf(n){;}

CspWfVector::CspWfVector(const CspWfVector& w) : TObject(w)
{
  name = w.name;
  boardID = w.boardID;
  chanID = w.chanID;
  N = w.N;
  wf = w.wf;
}


void
CspWfVector::add(const CspWaveform& w)
{
  N++;
  wf.push_back(w);

}

void
CspWfVector::insert(int i, const CspWaveform& w)
{
  N++;
  wf.insert(wf.begin()+i,w);
}

void
CspWfVector::swap(int i, const CspWaveform& w)
{
  wf[i] = w;
}

void
CspWfVector::rm(int i)
{
  N--;
  wf.erase(wf.begin()+i);
}

void
CspWfVector::clearWf()
{
  N=0;
  wf.clear();
}

void
CspWfVector::resize(int i)
{
  N=i;
  wf.resize(i);
}

CspWfVector&
CspWfVector::operator=(const CspWfVector& v)
{
  if (&v==this) return *this;
  N = v.N;
  name=v.name;
  wf=v.wf;
  return *this;
}


const CspWaveform& 
CspWfVector::at(int i) const{return wf[i];}

const CspWaveform& 
CspWfVector::operator()(int i) const{return wf[i];}

const CspWaveform& 
CspWfVector::operator[](int i) const{return wf[i];}

CspWaveform& 
CspWfVector::at(int i) {return wf[i];}

CspWaveform& 
CspWfVector::operator()(int i) {return wf[i];}

CspWaveform& 
CspWfVector::operator[](int i) {return wf[i];}

const CspPulse& 
CspWfVector::at(int i, int j) const{return wf[i][j];}

const CspPulse& 
CspWfVector::operator()(int i, int j) const{return wf[i][j];}

CspPulse& 
CspWfVector::at(int i, int j){return wf[i][j];}

CspPulse& 
CspWfVector::operator()(int i, int j){return wf[i][j];}

