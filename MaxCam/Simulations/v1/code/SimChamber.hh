#ifndef SIMCHAMBER_HH
#define SIMCHAMBER_HH

#include <vector>
#include "TString.h"
#include "TObject.h"
using std::vector;

class SimChamber : public TObject{

public:

  SimChamber();
  virtual ~SimChamber(){}

  void setHeight(double height){fHeight = height;}
  void setDriftLength(double len){fDriftLength = len;}
  void setDriftVoltage(double volt){fDriftVoltage = volt;}
  void setAnodeVoltage(double volt){fAnodeVoltage = volt;}
  void setTemperature(double temp){fTemperature = temp;}
  void setPressure(double pres){fPressure = pres;}
  void setDOverMu(double dom){fDOverMu = dom;}
  void setDOverMu();
  void setDiffusionConstantTerm(double diffcon){fDiffusionConstTerm = diffcon;}
  void setDiffusionDzTerm(double diffdz){fDOverMu = diffdz*fDriftVoltage/(2*fDriftLength);}
  void setElectronLifetime(double elifetime){fElectronLifetime = elifetime;}
  void setElectricScintillation(double scint){fElecScintillation = scint;}
  void setNuclearScintillation(double scint){fNuclScintillation = scint;}
  void setAttenuation(double atten){fAttenuation = atten;}
  void setElectronPerkeV(double elperkeV){fElectronPerkeV = elperkeV;} 
  void setDriftVelocity(double v){fDriftVelocity = v;}
  void setOrientationAngle(double ang){fOrientationAngle = ang;}

  double getHeight(){return fHeight;}
  double getDriftLength(){return fDriftLength;}
  double getDriftVoltage(){return fDriftVoltage;}
  double getAnodeVoltage(){return fAnodeVoltage;}
  double getTemperature(){return fTemperature;}
  double getPressure(){return fPressure;}
  double getDOverMu(){return fDOverMu;}
  double getDiffusionConstantTerm(){return fDiffusionConstTerm;}
  double getDiffusionDzTerm(){return 2*fDOverMu*fDriftLength / fDriftVoltage;}
  double getElectronLifetime(){return fElectronLifetime;}
  double getElectricScintillation(){return fElecScintillation;}
  double getNuclearScintillation(){return fNuclScintillation;}
  double getAttenuation(){return fAttenuation;}
  double getElectronPerkeV(){return fElectronPerkeV;}
  double getDriftVelocity(){return fDriftVelocity;}
  double getOrientationAngle(){return fOrientationAngle;}

  void setSpacerAxis(TString ax = "x"){fSpacerAxis = ax;}
  TString getSpacerAxis(){return fSpacerAxis;}

  void setLongDiffusionConstantTerm(double width){fLongDiffConstTerm = width;}
  void setLongDiffusionDzTerm(double dz){fLongDiffDzTerm = dz;}
  double getLongDiffusionConstantTerm(){return fLongDiffConstTerm;}
  double getLongDiffusionDzTerm(){return fLongDiffDzTerm;}

  void setSpacerWidth(double width){fSpacerWidth = width;}
  double getSpacerWidth(){return fSpacerWidth;}

  void addSpacer(double pos){fSpacerPositions.push_back(pos);}
  void clearSpacers(){fSpacerPositions.clear();}
  vector<double> getSpacers(){return fSpacerPositions;}

private:

  double fHeight;                    //Height in mm = |z(mesh)-z(lens)|
  double fDriftLength;               //Drift length in mm
  double fDriftVoltage;              //Drift voltage in V
  double fAnodeVoltage;              //Anode voltage in V
  double fTemperature;               //Gas temperature in K
  double fPressure;                  //Gas pressure in torr
  double fDOverMu;                   //D/mu in V
  double fDiffusionConstTerm;        //in mm^2
  double fElectronLifetime;          //electron lifetime in units TBD
  double fElecScintillation;         //% of e&m en loss scintillated
  double fNuclScintillation;         //% of nuclear en loss scintillated
  double fAttenuation;               //Attenuation of electrons per mm
  double fElectronPerkeV;            //Electrons at mesh/keV energy loss (work function * electron gain)
  double fDriftVelocity;             //Drift velocity in mm/ns
  vector<double> fSpacerPositions;   //spacer positions in mm
  TString fSpacerAxis;               //Axis along which spacers are aligned ("x" or "y")
  double fSpacerWidth;               //Spacer width in mm
  double fLongDiffConstTerm;         //Constant term in longitudinal diffusion
  double fLongDiffDzTerm;            //Dz term in longitudinal diffusion
  double fOrientationAngle;          //Degrees west (counterclockwise) from north that the detector +y axis is oriented.

  ClassDef(SimChamber,1)
};

#endif
