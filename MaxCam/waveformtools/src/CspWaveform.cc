#include "CspWaveform.hh"
#include "CspPulse.hh"

CspWaveform::CspWaveform(int n) : TObject(),
N(n), pulse(n),
base(0), rms(0), 
wfMax(0), wfMaxTime(0),wfMaxBin(0),
wfMin(0), wfMinTime(0),wfMinBin(0)
{;}


CspWaveform::CspWaveform(const CspWaveform& w) : TObject(w)
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

CspWaveform&
CspWaveform::operator=(const CspWaveform& w) 
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


const CspPulse& 
CspWaveform::at(int i) const{return pulse[i];}

const CspPulse& 
CspWaveform::operator()(int i)const{return pulse[i];}

const CspPulse& 
CspWaveform::operator[](int i) const
{return pulse[i];}

CspPulse& 
CspWaveform::at(int i){return pulse[i];}

CspPulse& 
CspWaveform::operator()(int i){return pulse[i];}

CspPulse& 
CspWaveform::operator[](int i){return pulse[i];}

void
CspWaveform::add(const CspPulse& p)
{
  N++;
  pulse.push_back(p);
}

void
CspWaveform::insert(int i, const CspPulse& p)
{
  N++;
  pulse.insert(pulse.begin()+i,p);
}

void
CspWaveform::swap(int i, const CspPulse& p)
{
  pulse[i] = p;
}

void
CspWaveform::rm(int i)
{
  N--;
  pulse.erase(pulse.begin()+i);
}

void
CspWaveform::clear()
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
CspWaveform::clearPulse()
{
  N=0;
  pulse.clear();
}

void
CspWaveform::resize(int i)
{
  N=i;
  pulse.resize(i);

}
