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



class SimMesh : public TObject{ /*stores variables and defines functions necessary for producing a waveform.*/

public:


  SimMesh();
  //empty constructor

  SimMesh(int plot_offset_EG, int imin_EG, double Initial_Energy_EG, double StepSize_EG, SimScope* scope_EG, double fMinEnergy_EG, vector<double>& fRecoilZ_EG, vector<double>& fRecoilEn_EG, SimChamber* fChamber_EG);
  //full constructor

  virtual ~SimMesh();
  //destructor

  void simulate();
  //runs a waveform simulation from an initialized instance of SimMesh



private:



  double noise(double error) {return rand->Gaus(0,error);} 
  //Mathematical function. Returns gaussian random number.
  
  double timenoise(double sz) {return (sz>0) ? a1 + a2 * pow(sz,0.5) : a1;} 
  //smears waveform on the time axis
  
  double j(double sz){return (t + sz/vd)/dt;} 
  //arrival time of electron
  
  double vel(){return (En<0) ? 0 : pow(2*En/mass,0.5);} 
  //converts energy to velocity
  
  void step(int i);
  //runs one step of simulation
  
  double vecsum(vector<double>& B, int risesteps);
  //Mathematical function. Computes sum of vector's elements.
  
  void decay();
  //simulates rise,decay times
  
  double round(double x, double res);
  //Mathematical function. Rounds x to nearest multiple of res.
  
  void plot();
  //plots waveform
  
  void plotfinish();
  //adds noise, vertically shifts, and discretizes waveform

  void noisify();
  //adds noise to waveform (Has been replaced by plotfinish).
  
  void allocate();
  //allocates space to Q,R

  int SetFinalPosition(ifin){final_position = ifin;}
  int GetFinalPosition(){return final_position;}
  //returns last used index of fRecoilZ. Used for compatibility with EventGenerator's double alphas.

  int SetPlotOffset(poff){plot_offset = poff;}
  int GetPlotOffset(){return plot_offset;}

  //see v1/runParameters/LittleDCTPC_far/scopeProperties.temp for more info on those of the variables here that are tune-able
  double t;  //current time
  double Tr; //rise time
  double Td; //decay time
  double vd; //drift velocity
  double Scale; //final vertical scaling
  double scoperes; //voltage resolution
  double rate; //sampling rate
  double dt; //scope time resolution
  double a1; //smearing (constant)
  double a2; //smearing (position-dependent)
  double SD; // voltage noise
  double mass; //alpha mass
  int zero_offset; //horizontal centering
  double zero_offset_var; //hor. centering (position-dependent)
  int Nsamples; //accuracy of smearing
  double vertical_offset; //vertical centering
  int need_space; //Q, R almost out of space
  int jmax_init; //amount of space allocated to Q,R at a time
  int imin; //beginning of position vector
  int imax; //end of position vector
  int j_init; //first expected electron arrival time
  int jmax; //current size of Q,R
  int jmin; //earliest detected electron
  double En; //current energy
  vector<double> Q; //electrons hit mesh               
  vector<double> R; //signal recorded              
  double Initial_Energy; //it does what it says. Comes from EventGenerator (a.k.a. EG)
  double dx; //step size
  int plotmax; //limits of plot
  int presize; //limits of plot
  int plotmin; //limits of plot
  int nscope; //currently unused
  SimScope* scope; //scope from (EG)
  TRandom3* rand; //random number generator
  vector<double> fRecoilEn; //energy loss from (EG)
  vector<double> fRecoilZ; //position from (EG)
  double fMinEnergy; //signals end of simulation (EG)
  SimChamber* fChamber; //Contains TPC parameters (EG)
  int final_position; //the last index of fRecoilZ that was called before exiting step() loop.
  int plot_offset; //keeps track of the plotmin value used in plotting. Needed for correct double-alpha waveforms.

  ClassDef(SimMesh,1)
    };

#endif
