#include "FastWaveform.hh"
#include "FastPulse.hh"

FastWaveform::FastWaveform(int n) : TObject(),
N(n), pulse(n),
base(0), rms(0), 
wfMax(0), wfMaxTime(0),wfMaxBin(0),
				    wfMin(0), wfMinTime(0),wfMinBin(0),
				    leftint(0), leftint2(0), rightint(0), rightint2(0),
				    jmin(0), jbragg(0), jterm(0), jterm2(0), jterm3(0), origin(0),
				    rec_SD(0),wfd_delta(0),
				    peak1(0),peak2(0),peak1val(0),peak2val(0),
				    half1(0),half2(0),
				    termdist(0),maxdevloc(0),
				    rms_left(0),rms_right(0),rms_outer(0),rms_full(0),
				    dt(0)
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
  leftint=w.leftint;
  leftint2=w.leftint2;
  rightint=w.rightint;
  rightint2=w.rightint2;
  jmin=w.jmin;
  jbragg=w.jbragg;
  jterm=w.jterm;
  jterm2=w.jterm2;
  jterm3=w.jterm3;
  origin=w.origin;
  rec_SD=w.rec_SD;
  wfd_delta=w.wfd_delta;
  peak1=w.peak1;
  peak2=w.peak2;
  peak1val=w.peak1val;
  peak2val=w.peak2val;
  half1=w.half1;
  half2=w.half2;
  termdist=w.termdist;
  maxdevloc=w.maxdevloc;
  rms_left=w.rms_left;
  rms_right=w.rms_right;
  rms_outer=w.rms_outer;
  rms_full=w.rms_full;
  dt=w.dt;
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
  leftint=w.leftint;
  leftint2=w.leftint2;
  rightint=w.rightint;
  rightint2=w.rightint2;
  jmin=w.jmin;
  jbragg=w.jbragg;
  jterm=w.jterm;
  jterm2=w.jterm2;
  jterm3=w.jterm3;
  origin=w.origin;
  rec_SD=w.rec_SD;
  wfd_delta=w.wfd_delta;
  peak1=w.peak1;
  peak2=w.peak2;
  peak1val=w.peak1val;
  peak2val=w.peak2val;
  half1=w.half1;
  half2=w.half2;
  termdist=w.termdist;
  maxdevloc=w.maxdevloc;
  rms_left=w.rms_left;
  rms_right=w.rms_right;
  rms_outer=w.rms_outer;
  rms_full=w.rms_full;
  dt=w.dt;
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
  leftint=0;
  leftint2=0;
  rightint=0;
  rightint2=0;
  jmin=0;
  jbragg=0;
  jterm=0;
  jterm2=0;
  jterm3=0;
  origin=0;
  rec_SD=0;
  wfd_delta=0;
  peak1=0;
  peak2=0;
  peak1val=0;
  peak2val=0;
  half1=0;
  half2=0;
  termdist=0;
  maxdevloc=0;
  rms_left=0;
  rms_right=0;
  rms_outer=0;
  rms_full=0;
  dt=0;
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
