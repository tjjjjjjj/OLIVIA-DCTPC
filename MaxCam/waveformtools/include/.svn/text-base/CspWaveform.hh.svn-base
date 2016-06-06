#ifndef CSP_WAVEFORM
#define CSP_WAVEFORM

#include "TObject.h"
#include <vector>
class CspPulse;

class CspWaveform : public TObject
{

  public:
    CspWaveform(int n=0);
    CspWaveform(const CspWaveform& w);
    virtual ~CspWaveform(){;}

    void add(const CspPulse& p);
    void insert(int i, const CspPulse& p);
    void swap(int i, const CspPulse& p);
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

    CspWaveform& operator=(const CspWaveform& w);

    const CspPulse& at(int i) const;
    const CspPulse& operator()(int i) const;
    const CspPulse& operator[](int i) const;
    CspPulse& at(int i);
    CspPulse& operator()(int i);
    CspPulse& operator[](int i);
    

  protected:
    int N;
    std::vector<CspPulse> pulse;
    double base;
    double rms;
    double wfMax;
    double wfMaxTime;
    int wfMaxBin;
    double wfMin; 
    double wfMinTime; 
    int wfMinBin;

    ClassDef(CspWaveform,1)
};


#endif
