#include "PMTWfVector.hh"
#include "PMTPulse.hh"
#include "PMTWaveform.hh"
#include <vector>
#include <string>
#include "TObject.h"
using std::vector;


ClassImp(PMTWfVector)

PMTWfVector::PMTWfVector(int n, int bid, int cid,const char* nm) :
name(nm), boardID(bid), chanID(cid), N(n), wf(n){;}

PMTWfVector::PMTWfVector(const PMTWfVector& w) : TObject(w)
{
  name = w.name;
  boardID = w.boardID;
  chanID = w.chanID;
  N = w.N;
  wf = w.wf;
}


void
PMTWfVector::add(const PMTWaveform& w)
{
  N++;
  wf.push_back(w);

}

void
PMTWfVector::insert(int i, const PMTWaveform& w)
{
  N++;
  wf.insert(wf.begin()+i,w);
}

void
PMTWfVector::swap(int i, const PMTWaveform& w)
{
  wf[i] = w;
}

void
PMTWfVector::rm(int i)
{
  N--;
  wf.erase(wf.begin()+i);
}

void
PMTWfVector::clearWf()
{
  N=0;
  wf.clear();
}

void
PMTWfVector::resize(int i)
{
  N=i;
  wf.resize(i);
}

PMTWfVector&
PMTWfVector::operator=(const PMTWfVector& v)
{
  if (&v==this) return *this;
  N = v.N;
  name=v.name;
  wf=v.wf;
  return *this;
}


const PMTWaveform& 
PMTWfVector::at(int i) const{return wf[i];}

const PMTWaveform& 
PMTWfVector::operator()(int i) const{return wf[i];}

const PMTWaveform& 
PMTWfVector::operator[](int i) const{return wf[i];}

PMTWaveform& 
PMTWfVector::at(int i) {return wf[i];}

PMTWaveform& 
PMTWfVector::operator()(int i) {return wf[i];}

PMTWaveform& 
PMTWfVector::operator[](int i) {return wf[i];}

const PMTPulse& 
PMTWfVector::at(int i, int j) const{return wf[i][j];}

const PMTPulse& 
PMTWfVector::operator()(int i, int j) const{return wf[i][j];}

PMTPulse& 
PMTWfVector::at(int i, int j){return wf[i][j];}

PMTPulse& 
PMTWfVector::operator()(int i, int j){return wf[i][j];}

