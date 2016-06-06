#include "SkimWaveform.hh"
#include "DmtpcPulse.hh"


SkimWaveform::SkimWaveform(int n) : TObject(),
N(n), pulse(n),
base(0), rms(0), 
wfMin(0), wfMinTime(0),wfMinBin(0),
wfMax(0), wfMaxTime(0),wfMaxBin(0)

{;}

SkimWaveform::SkimWaveform(const SkimWaveform& w) : TObject(w)
{
  N=w.size();
  pulse=w.pulse;
  base=w.base;
  rms=w.rms;

  wfMin=w.wfMin;
  wfMinTime=w.wfMinTime;
  wfMinBin=w.wfMinBin;

  wfMax=w.wfMax;
  wfMaxTime=w.wfMaxTime;
  wfMaxBin=w.wfMaxBin;
}

SkimWaveform&
SkimWaveform::operator=(const SkimWaveform& w) 
{
  if (&w==this) return *this;
  N=w.size();
  pulse=w.pulse;
  base=w.base;
  rms=w.rms;

  wfMin=w.wfMin;
  wfMinTime=w.wfMinTime;
  wfMinBin=w.wfMinBin;

  wfMax=w.wfMax;
  wfMaxTime=w.wfMaxTime;
  wfMaxBin=w.wfMaxBin;

  return *this;
}

const DmtpcPulse& 
SkimWaveform::at(int i) const{return pulse[i];}

const DmtpcPulse& 
SkimWaveform::operator()(int i)const{return pulse[i];}

const DmtpcPulse& 
SkimWaveform::operator[](int i) const
{return pulse[i];}

DmtpcPulse& 
SkimWaveform::at(int i){return pulse[i];}

DmtpcPulse& 
SkimWaveform::operator()(int i){return pulse[i];}

DmtpcPulse& 
SkimWaveform::operator[](int i){return pulse[i];}

void
SkimWaveform::add(const DmtpcPulse& p)
{
  N++;
  pulse.push_back(p);
}

void
SkimWaveform::insert(int i, const DmtpcPulse& p)
{
  N++;
  pulse.insert(pulse.begin()+i,p);
}

void
SkimWaveform::swap(int i, const DmtpcPulse& p)
{
  pulse[i] = p;
}

void
SkimWaveform::rm(int i)
{
  N--;
  pulse.erase(pulse.begin()+i);
}

void
SkimWaveform::clear()
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
SkimWaveform::clearPulse()
{
  N=0;
  pulse.clear();
}

void
SkimWaveform::resize(int i)
{
  N=i;
  pulse.resize(i);
}


