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

  double getLeftInt() const {return leftint;}
  double getLeftInt2() const {return leftint2;}
  double getRightInt() const {return rightint;}
  double getRightInt2() const {return rightint2;}
  int getjMin() const {return jmin;}
  int getjBragg() const {return jbragg;}
  int getjTerm() const {return jterm;}
  int getjTerm2() const {return jterm2;}
  int getjTerm3() const {return jterm3;}
  int getOrigin() const {return origin;}
  double getRecSD() const {return rec_SD;}
  int getWfdDelta() const {return wfd_delta;}
  int getPeak1() const {return peak1;}
  int getPeak2() const {return peak2;}
  double getPeak1Val() const {return peak1val;}
  double getPeak2Val() const {return peak2val;}
  int getHalf1() const {return half1;}
  int getHalf2() const {return half2;}
  double getTermDist() const {return termdist;}
  int getMaxDevLoc() const {return maxdevloc;}
  double getRMSLeft() const {return rms_left;}
  double getRMSRight() const {return rms_right;}
  double getRMSOuter() const {return rms_outer;}
  double getRMSFull() const {return rms_full;}
  double getdt() const {return dt;}

  void setLeftInt(double x){leftint=x;}
  void setLeftInt2(double x){leftint2=x;}
  void setRightInt(double x){rightint=x;}
  void setRightInt2(double x){rightint2=x;}
  void setjMin(int n){jmin=n;}
  void setjBragg(int n){jmin=n;}
  void setjTerm(int n){jterm=n;}
  void setjTerm2(int n){jterm2=n;}
  void setjTerm3(int n){jterm3=n;}
  void setOrigin(int n){origin=n;}
  void setRecSD(double x){rec_SD=x;}
  void setWfdDelta(int n){wfd_delta=n;}
  void setPeak1(int n){peak1=n;}
  void setPeak2(int n){peak2=n;}
  void setPeak1Val(double x){peak1val=x;}
  void setPeak2Val(double x){peak2val=x;}
  void setHalf1(int n){half1=n;}
  void setHalf2(int n){half2=n;}
  void setTermDist(double x){termdist=x;}
  void setMaxDevLoc(int n){maxdevloc=n;}
  void setRMSLeft(double x){rms_left=x;}
  void setRMSRight(double x){rms_right=x;}
  void setRMSOuter(double x){rms_outer=x;}
  void setRMSFull(double x){rms_full=x;}
  void setdt(double x){dt=x;}
  

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

  double leftint;
  double leftint2;
  double rightint;
  double rightint2;
  int jmin;
  int jbragg;
  int jterm;
  int jterm2;
  int jterm3;
  int origin;
  double rec_SD;
  int wfd_delta;
  int peak1;
  int peak2;
  double peak1val;
  double peak2val;
  int half1;
  int half2;
  double termdist;
  double maxdevloc;
  double rms_left;
  double rms_right;
  double rms_outer;
  double rms_full;
  double dt;

    ClassDef(FastWaveform,1)
};


#endif
