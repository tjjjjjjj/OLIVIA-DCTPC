#ifndef WAVEFORM_VECTOR
#define WAVEFORM_VECTOR

#include "SkimWaveform.hh"
#include "DmtpcPulse.hh"

class SkimWaveform;
class DmtpcPulse;

#include "TObject.h"
#include <vector>
#include <string>

class WaveformVector : public TObject
{

  public:
    WaveformVector(int n=0,int bid=0, int cid=0, const char* nm="Pulse");
    WaveformVector(const WaveformVector& w);
    virtual ~WaveformVector(){;}
 
    const char* GetName() const {return name.c_str();}
    const char* getName() const {return GetName();}

    void setBoard(int id){boardID = id;}
    void setChan(int id){chanID = id;}

    int getBoard(){return boardID;}
    int getChan(){return chanID;}

    int size() const {return N;}
    void add(const SkimWaveform& w);
    void insert(int i, const SkimWaveform& w);
    void swap(int i, const SkimWaveform& w);
    void rm(int i);
    void clearWf();
    void clear(){ clearWf();}
    void resize(int i);

    WaveformVector& operator=(const WaveformVector& w);

    const SkimWaveform& at(int i) const;
    const SkimWaveform& operator()(int i) const;
    const SkimWaveform& operator[](int i) const;
    SkimWaveform& at(int i) ;
    SkimWaveform& operator()(int i) ;
    SkimWaveform& operator[](int i) ;

    const DmtpcPulse& at(int i, int j) const;
    const DmtpcPulse& operator()(int i, int j) const;
    DmtpcPulse& at(int i, int j);
    DmtpcPulse& operator()(int i, int j);
  

  protected:
    std::string name;
    int boardID;
    int chanID;
    int N;
    std::vector<SkimWaveform> wf;
    ClassDef(WaveformVector,1)

};


#endif
