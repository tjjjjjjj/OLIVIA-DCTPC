#ifndef __SIMMESH__
#define __SIMMESH__

//ROOT includes
#include "TMath.h"
#include "TObject.h"
#include "TString.h"

//code includes
#include "EventGenerator.hh"
#include "SimScope.hh"
#include "SimChamber.hh"

//class EventGenerator; I might have no idea what I'm doing
//class SimScope;

class SimMesh : public TObject{

public:

  SimMesh();

  SimMesh(double Initial_Energy_EG, double StepSize_EG, SimScope* scope_EG, int nscope_EG, double fMinEnergy_EG, vector<double>& fRecoilZ_EG, vector<double>& fRecoilEn_EG, SimChamber* fChamber);
  virtual ~SimMesh();
  void simulate();

private:

  double noise(double error) {return rand->Gaus(0,error);}
  double timenoise(double sz) {return (sz>0) ? a1 + a2 * pow(sz,0.5) : a1;}
  double j(double sz){return (t + sz/vd)/dt;}
  double vel(){return (En<0) ? 0 : pow(2*En/mass,0.5);}
  void step(int i);
  double vecsum(vector<double>& B, int risesteps);
  void decay();
  double round(double x, double res);
  void plot();
  void noisify();
  void allocate();

  double t;
  double Tr;
  double Td;
  double vd;
  double Scale;
  double scoperes;
  double rate;
  double dt;
  double a1;
  double a2;
  double SD;
  double mass;
  int zero_offset;
  double zero_offset_var;
  int Nsamples;
  double vertical_offset;
  int need_space;
  int jmax_init;
  int imin;
  int imax;
  int j_init;
  int jmax;
  int jmin;
  double En;
  vector<double> Q; //electrons hit mesh               
  vector<double> R; //signal recorded              
  double Initial_Energy;
  double dx;
  int plotmax;
  int presize;
  int plotmin;
  int nscope;
  SimScope* scope;
  TRandom3* rand;
  vector<double> fRecoilEn;
  vector<double> fRecoilZ;
  double fMinEnergy;
  SimChamber* fChamber;

  ClassDef(SimMesh,1)
    };

#endif