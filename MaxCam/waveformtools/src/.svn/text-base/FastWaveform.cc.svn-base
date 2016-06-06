#include "FastWaveform.hh"
#include "FastPulse.hh"

FastWaveform::FastWaveform(int n) : TObject(),
N(n), pulse(n),
base(0), rms(0), 
wfMax(0), wfMaxTime(0),wfMaxBin(0),
wfMin(0), wfMinTime(0),wfMinBin(0)
{;}

FastWaveform::FastWaveform(const FastWaveform& w) : TObject(w)
{
  N=w.size();
  pulse=w.pulse;
  base=w.base;
  rms=w.rms;
  wfMax=w.wfMax;
  wfMaxTime=w.wfMaxTime;
  wfMaxBin=w.wfMaxBin;
  wfMin=w.wfMin;
  wfMinTime=w.wfMinTime;
  wfMinBin=w.wfMinBin;
}

FastWaveform&
FastWaveform::operator=(const FastWaveform& w) 
{
  if (&w==this) return *this;
  N=w.size();
  pulse=w.pulse;
  base=w.base;
  rms=w.rms;
  wfMax=w.wfMax;
  wfMaxTime=w.wfMaxTime;
  wfMaxBin=w.wfMaxBin;
  wfMin=w.wfMin;
  wfMinTime=w.wfMinTime;
  wfMinBin=w.wfMinBin;
  return *this;
}

const FastPulse& 
FastWaveform::at(int i) const{return pulse[i];}

const FastPulse& 
FastWaveform::operator()(int i)const{return pulse[i];}

const FastPulse& 
FastWaveform::operator[](int i) const
{return pulse[i];}

FastPulse& 
FastWaveform::at(int i){return pulse[i];}

FastPulse& 
FastWaveform::operator()(int i){return pulse[i];}

FastPulse& 
FastWaveform::operator[](int i){return pulse[i];}

void
FastWaveform::add(const FastPulse& p)
{
  N++;
  pulse.push_back(p);
}

void
FastWaveform::insert(int i, const FastPulse& p)
{
  N++;
  pulse.insert(pulse.begin()+i,p);
}

void
FastWaveform::swap(int i, const FastPulse& p)
{
  pulse[i] = p;
}

void
FastWaveform::rm(int i)
{
  N--;
  pulse.erase(pulse.begin()+i);
}

void
FastWaveform::clear()
{
  clearPulse();
  base=0;
  rms=0;
  wfMin=0;
  wfMinTime=0;
  wfMinBin=0;
  wfMax = 0;
  wfMaxTime=0;
  wfMaxBin=0;
}

void
FastWaveform::clearPulse()
{
  N=0;
  pulse.clear();
}

void
FastWaveform::resize(int i)
{
  N=i;
  pulse.resize(i);
}
