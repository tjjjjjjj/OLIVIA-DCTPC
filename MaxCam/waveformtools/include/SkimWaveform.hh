#ifndef DMTPC_SKIM_WAVEFORM
#define DMTPC_SKIM_WAVEFORM

#include "TObject.h"
#include <vector>
class DmtpcPulse;



class SkimWaveform : public TObject
{

  public:
    SkimWaveform(int n=0);
    SkimWaveform(const SkimWaveform& w);
    virtual ~SkimWaveform(){;}

    void add(const DmtpcPulse& p);
    void insert(int i, const DmtpcPulse& p);
    void swap(int i, const DmtpcPulse& p);
    void rm(int i);
    void clearPulse();
    void clear();
    void resize(int i); 
   
    int size() const {return N;}

    
    double getTime() const {return time;}
    double getBase() const {return base;}
    double getRMS() const {return rms;}
    double getWfMax() const {return wfMax;}
    double getWfMaxTime() const {return wfMaxTime;}
    int getWfMaxBin() const {return wfMaxBin;}
    double getWfMin() const {return wfMin;}
    double getWfMinTime() const {return wfMinTime;}
    int getWfMinBin() const {return wfMinBin;}

    void setTime(double x){time = x;}
    void setBase(double x){base = x;}
    void setRMS(double x){rms=x;}
    void setWfMax(double x){wfMax=x;}
    void setWfMaxTime(double x){wfMaxTime=x;}
    void setWfMaxBin(int x){wfMaxBin=x;}
    void setWfMin(double x){wfMin=x;}
    void setWfMinTime(double x){wfMinTime=x;}
    void setWfMinBin(int x){wfMinBin=x;}

    SkimWaveform& operator=(const SkimWaveform& w);

    const DmtpcPulse& at(int i) const;
    const DmtpcPulse& operator()(int i) const;
    const DmtpcPulse& operator[](int i) const;
    DmtpcPulse& at(int i);
    DmtpcPulse& operator()(int i);
    DmtpcPulse& operator[](int i);
    

  protected:
    int N;
    std::vector<DmtpcPulse> pulse;
    double base;
    double rms;
 
    double time;

    double wfMin; 
    double wfMinTime; 
    int wfMinBin;

    double wfMax;
    double wfMaxTime;
    int wfMaxBin;

    ClassDef(SkimWaveform,2)
};


#endif
