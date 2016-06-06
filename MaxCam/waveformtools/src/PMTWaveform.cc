#include "PMTWaveform.hh"
#include "PMTPulse.hh"

PMTWaveform::PMTWaveform(int n) : TObject(),
N(n), pulse(n),
base(0), rms(0), 
wfMax(0), wfMaxTime(0),wfMaxBin(0),
wfMin(0), wfMinTime(0),wfMinBin(0)
{;}


PMTWaveform::PMTWaveform(const PMTWaveform& w) : TObject(w)
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

PMTWaveform&
PMTWaveform::operator=(const PMTWaveform& w) 
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


const PMTPulse& 
PMTWaveform::at(int i) const{return pulse[i];}

const PMTPulse& 
PMTWaveform::operator()(int i)const{return pulse[i];}

const PMTPulse& 
PMTWaveform::operator[](int i) const
{return pulse[i];}

PMTPulse& 
PMTWaveform::at(int i){return pulse[i];}

PMTPulse& 
PMTWaveform::operator()(int i){return pulse[i];}

PMTPulse& 
PMTWaveform::operator[](int i){return pulse[i];}

void
PMTWaveform::add(const PMTPulse& p)
{
  N++;
  pulse.push_back(p);
}

void
PMTWaveform::insert(int i, const PMTPulse& p)
{
  N++;
  pulse.insert(pulse.begin()+i,p);
}

void
PMTWaveform::swap(int i, const PMTPulse& p)
{
  pulse[i] = p;
}

void
PMTWaveform::rm(int i)
{
  N--;
  pulse.erase(pulse.begin()+i);
}

void
PMTWaveform::clear()
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
PMTWaveform::clearPulse()
{
  N=0;
  pulse.clear();
}

void
PMTWaveform::resize(int i)
{
  N=i;
  pulse.resize(i);

}
