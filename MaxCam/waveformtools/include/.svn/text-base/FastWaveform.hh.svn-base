#ifndef FAST_WAVEFORM
#define FAST_WAVEFORM

class FastPulse;
#include "TObject.h"
#include <vector>
class FastWaveform : public TObject
{

  public:
    FastWaveform(int n=0);
    FastWaveform(const FastWaveform& w);
    virtual ~FastWaveform(){;}

    void add(const FastPulse& p);
    void insert(int i, const FastPulse& p);
    void swap(int i, const FastPulse& p);
    void rm(int i);
    void clearPulse();
    void clear();
    void resize(int i);    
   
    int size() const {return N;}

    double getBase() const {return base;}
    double getRMS() const {return rms;}
    double getWfMax() const {return wfMax;}
    double getWfMaxTime() const {return wfMaxTime;}
    int getWfMaxBin() const {return wfMaxBin;}
    double getWfMin() const {return wfMin;}
    double getWfMinTime() const {return wfMinTime;}
    int getWfMinBin() const {return wfMinBin;}

    void setBase(double x){base = x;}
    void setRMS(double x){rms=x;}
    void setWfMax(double x){wfMax=x;}
    void setWfMaxTime(double x){wfMaxTime=x;}
    void setWfMaxBin(int x){wfMaxBin=x;}
    void setWfMin(double x){wfMin=x;}
    void setWfMinTime(double x){wfMinTime=x;}
    void setWfMinBin(int x){wfMinBin=x;}

    FastWaveform& operator=(const FastWaveform& w);

    const FastPulse& at(int i) const;
    const FastPulse& operator()(int i) const;
    const FastPulse& operator[](int i) const;
    FastPulse& at(int i);
    FastPulse& operator()(int i);
    FastPulse& operator[](int i);
    

  protected:
    int N;
    std::vector<FastPulse> pulse;
    double base;
    double rms;
    double wfMax;
    double wfMaxTime;
    int wfMaxBin;
    double wfMin; 
    double wfMinTime; 
    int wfMinBin;

    ClassDef(FastWaveform,1)
};


#endif
