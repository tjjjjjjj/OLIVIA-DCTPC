#ifndef EVENTGENERATOR_HH
#define EVENTGENERATOR_HH

#include "../../../DmtpcDataset.hh"
#include "../../../MaxCamSRIM.hh"
#include "../../../MaxCamElectricField.hh"
#include "../../../MaxCamClusterImage.hh"
#include "ParticleGenerator.hh"
#include "SimCamera.hh"
#include "SimPMT.hh"
#include "SimScope.hh"
#include "SimChamber.hh"

#include "SimMesh.hh"

#include "TH1D.h"
#include "TH2F.h"
#include "TString.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TVector3.h"
#include "TObject.h"
#include "TVectorT.h"
#include "TF1.h"
#include "TH1.h"
#include "TObjArray.h"
#include <math.h>
#include <vector>
using std::vector;

class SimMesh;
class EventGenerator : public TObject {

public:

  EventGenerator();
  ~EventGenerator();

  ParticleGenerator* 	particleGen(){return fPartGen;}
  MaxCamSRIM* 		srim(){return fSrim;}
  SimChamber* 		chamber(){return fChamber;}
  TObjArray* 		camera(){return fCamera;}
  SimCamera* 		camera(int i);  

  TObjArray*            pmt() {return fPMT;}
  SimPMT*               pmt(int i); 

  TObjArray*            scope() {return fScope;}
  SimScope*             scope(int i); 

  SimMesh*              simMesh(){return fMesh;}

  void 		reset();
  void 		setStepSize(double step){fStepSize = step;}
  double 	getStepSize(){return fStepSize;}
  void 		setMinEnergy(double minen){fMinEnergy = minen;}
  double 	getMinEnergy(){return fMinEnergy;}

  void 		setSrimName(TString name){fSrimName = name;}
  void 		setSrim();
  void 		setPressure(double pres);
  void 		setSkew(double skew){fSkew=skew;}
  TString 	getSrimName(){return fSrimName;}

//  TH1F*		getElecSignal(){return fSignal;}
//  TF1*		getStraggleProfile(){return fLongStraggle;}

  double 	getLength(){return fLength;}
  double        getZLength(){return fZLength;}
  double	getTrackTheta(){return fTrackTheta;}
  double	getTrackPhi(){return fTrackPhi;}

  void 		setTimeResolution(double tr){fTimeResolution = tr;}
  double 	getTimeResolution(){return fTimeResolution;}
  double    getSkew(){return fSkew;}
  //void 		setTimeResolutionPMT(double tr){fTimeResolutionPMT = tr;}
  //double 	getTimeResolutionPMT(){return fTimeResolutionPMT;}

  TObjArray*    getTrueClusterArray(){ return fTrueClusterArray; }

  void 		useSpacers(bool use = true){fUseSpacers = use;}
  void 		useRadialEffect(bool use = true){fUseRadialEffect = use;}
  void          useLandauEffect(bool use = true){fUseLandauEffect = use;}
  bool 		propagateRecoil();
  void 		driftElectrons();
  void 		runEvent(bool resetImage = true);
  double 	applyDiffusion(double x,double z);
  double    applyElectronLifetime(double x,double z);
  Int_t     getSequence(); 
  void 		applyRadialEffect(int ncam);
  void      applyTPCImageEffect(int ncam);
  void      applyGainLandauEffect(int ncam);
  void 		applySpacers(int ncam);
  void 		makeImage(int ncam);
//  void 		makeElectronicSignal();
//  void 		setElectronicSignal();
//  TH1F* 	getElectronicSignal(){return fSignal;}

  void          makePMTSignal(Int_t ipmt, Int_t itrig);


  Double_t acceptance(Double_t D, Double_t h);
  Double_t solidAngleSphericalCap(Double_t D, Double_t h);

private:
  
  SimMesh *fMesh;//TJ's

  TFile *f;
  TH1D *hist;
  ParticleGenerator *fPartGen;    //ParticleGenerator Object
  DmtpcDataset *fData;
  MaxCamSRIM *fSrim;              //MaxCamSRIM Object
  SimChamber *fChamber;           //SimChamber Object
  TObjArray *fCamera;             //SimCamera Object
  TObjArray *fPMT;                //SimPMT Object
  TObjArray *fScope;              //SimScope Object
//  TH1F *fSignal;                  //Electronic signal
  //TF1 *fLongStraggle;             //Longitudinal straggling profile
  double fStepSize;               //Step size in mm
  double fMinEnergy;              //Min particle energy before stopping loop [keV]
  double fTempEnergy;             //Holds current energy in loop [keV]
  double fTempEnergy2;             //Holds current energy in loop [keV] 
  double fTimeResolution;         //Electronics time resolution in ns
  //TH1F  *fPMTSignal;              //PMT signal in time (x-axis=ns)
  //double fTimeResolutionPMT;      //PMT electronics time resolution in ns
  double fLength;                 //length of track in mm (dist from start to finish)
  double fZLength;                //length of z part of track in in mm 
  double fSkew;
  double fLengthScaleFactor;      //Multiply all lengths in the recoil by this to account for straggling
  double fTrackPhi;               //Azimuth angle on (-pi,pi] where +x=0, +y=pi/2, -x=pi [rad]
  double fTrackTheta;             //Polar angle on [0,pi] where +z is 0 and -z is pi. [rad]
  double fStraggle;
  TString fSrimName;              //Title of SRIM file
  vector<double> fRecoilX;        //x of each step in track [mm]
  vector<double> fRecoilY;        //y of each step in track [mm]
  vector<double> fRecoilZ;        //z of each step in track [mm]
  vector<double> fRecoilEn;       //en released at each step [keV]
  double fRecoilZ_min;            //minimum z-coordinate of a track [mm]
  double fRecoilZ_max;            //maximum z-coordinate of a track [mm]
  bool fUseSpacers;               //True: add spacer effects
  bool fUseRadialEffect;          //True: add cos^3 radial effect
  bool fUseTPCImageEffect;         //True: add TPC image effect
  bool fUseLandauEffect;          //True: add gain variation according to a Landau distribution
  bool fUseDoubleAlpha;           //True: propagate two alphas, instead of one nuclear recoil
  TVector3 fTempRecoil;           //Current particle direction in loop
  TVector3 fTempRecoil2;           //Current particle direction in loop 
  TVector3 fTempPosition;         //Current particle position in loop
  TRandom3 *fRnd;                 //Random number generator
  TObjArray *fTrueClusterArray;
  ClassDef(EventGenerator,1)
};

#endif
