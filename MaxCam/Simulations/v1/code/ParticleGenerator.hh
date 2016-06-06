#ifndef __PARTICLEGENERATOR__
#define __PARTICLEGENERATOR__

#include "../../../MaxCamENDF.hh"
#include "../../../DmtpcDecayChain.hh"
#include "SimTools.hh"
#include "SimRings.hh"

#include "TString.h"
#include "TMath.h"
#include "TVector3.h"
#include "TObject.h"
#include "TTimeStamp.h"
#include "TDatime.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"


class ParticleGenerator : public TObject {

public:

  ParticleGenerator();
  virtual ~ParticleGenerator();

  //Source position vector
  void setSourcePosition(TVector3* pos) {fSourcePos = pos;}
  void setSourcePosition(TVector3 pos) 
  {
    fSourcePos->SetXYZ(pos.X(),pos.Y(),pos.Z());
  }
  void setSourcePosition(Double_t x, Double_t y, Double_t z)
  {
    fSourcePos->SetXYZ(x,y,z);
  }
  TVector3* getSourcePosition(){return fSourcePos;}
  TVector3 sourcePosition();

  //Recoil position vector
  void setRecoilPosition(TVector3* pos) {fRecoilPos = pos;}
  void setRecoilPosition(TVector3 pos) 
  {
    fRecoilPos->SetXYZ(pos.X(),pos.Y(),pos.Z());
  }
  void setRecoilPosition(Double_t x, Double_t y, Double_t z)
  {
    fRecoilPos->SetXYZ(x,y,z);
  }
  TVector3* getRecoilPosition(){return fRecoilPos;}
  TVector3 recoilPosition();
  
  //Recoil position limits
  void setPositionMaxValues(Double_t x, Double_t y, Double_t z);
  void setPositionMinValues(Double_t x, Double_t y, Double_t z);
  void setRingValues(int n, Double_t x0, Double_t y0, Double_t rin, Double_t rout, 
                     Double_t zMin, Double_t zMax, Double_t dz, Double_t topz,
                     Double_t fractop, Double_t bottomz, Double_t frac_bot ); 

  //Projectile and recoil energies
  Double_t getProjectileEnergy(){return fProjEn;}
  void setProjectileEnergy(Double_t en){fProjEn = en;}
  Double_t getProjectileGaussianEnergy(){return fProjGaussianEn;}
  void setProjectileGaussianEnergy(Double_t gausen){fProjGaussianEn = gausen;}
  Double_t getProjectileGaussianEnergySpread(){return fProjGaussianEnSpread;}
  void setProjectileGaussianEnergySpread(Double_t en){fProjGaussianEnSpread = en;}
  Double_t getRecoilEnergy(){return fRecoilEn; }
  Double_t getRecoilEnergy2(){return fRecoilEn2;}
  void setRecoilEnergy(Double_t en){fRecoilEn = en;}
  Double_t getCosCygnus() {return fCosCygnus;}
  void setCosCygnus(Double_t coscyg) {fCosCygnus=coscyg;}

  //Projectile direction vector
  void setProjectileVector(TVector3* vect);
  void setProjectileVector(TVector3 vect){
    fProjVector->SetMagThetaPhi(1,vect.Theta(),vect.Phi());
  }
  void setProjectileVector(Double_t theta, Double_t phi)
  {
    fProjVector->SetMagThetaPhi(1,theta,phi);
  }
  void setProjectileVector(Double_t x, Double_t y, Double_t z);
  TVector3* getProjectileVector(){return fProjVector;}
  TVector3 projectileVector();

  //Recoil direction vector
  void setRecoilVector(TVector3* vect);
  void setRecoilVector(TVector3 vect){
    fRecVector->SetMagThetaPhi(1,vect.Theta(),vect.Phi());
  }
  void setRecoilVector(Double_t theta, Double_t phi){
    fRecVector->SetMagThetaPhi(1,theta,phi);
  }
  void setRecoilVector(Double_t x, Double_t y, Double_t z);
  TVector3* getRecoilVector(){return fRecVector;}
  
  TVector3* getRecoilVector2(){return fRecVector2;}
  TVector3 recoilVector();
  TVector3 recoilVector2();
  
  //Source Direction vector
  void setSourceDir(TVector3* vect);
  void setSourceDir(TVector3 vect)
  {
    fSourceDir->SetMagThetaPhi(1,vect.Theta(),vect.Phi());
  }
  void setSourceDir(Double_t theta, Double_t phi)
  {
    fSourceDir->SetMagThetaPhi(1,theta,phi);
  }
  void setSourceDir(Double_t x, Double_t y, Double_t z);
  TVector3* getSourceDir(){return fSourceDir;}
  TVector3 sourceDir();


  //Energy Limits
  void setMaxEnergy(Double_t en){fMaxEn = en;}
  void setMinEnergy(Double_t en){fMinEn = en;}
  Double_t getMaxEnergy(){return fMaxEn;}
  Double_t getMinEnergy(){return fMinEn;}
  void setEnergySpectrum(TH1* spec){fEnSpec = spec;}
  TH1* getEnergySpectrum(){return fEnSpec;}

  void setSeqNum(Int_t seqnum){fSeq=seqnum;}
  Int_t getSeqNum(){return fSeq;}	
  //Times
  void setTime(Double_t time){fTime = time; fTimeSet = time;}
  void setTime(TDatime *time){fTime = SimTools::getTimeDouble(time); fTimeSet = fTime;}
  void setTime(TTimeStamp *time){fTime = SimTools::getTimeDouble(time); fTimeSet = fTime;}
  Double_t getTime(){return fTime;}
  TTimeStamp *getTimeStamp(){return SimTools::makeTimeStamp(fTime);}
  TDatime *getTDatime(){return SimTools::makeTDatime(fTime);}
  void setTimeLimits(Double_t min, Double_t max){fMinTime = min;fMaxTime = max;}
  void setBeginTime(Int_t year,Int_t month,Int_t day,Int_t hour, Int_t minute, Int_t sec);
  void setBeginTime(TTimeStamp *time){fMinTime = SimTools::getTimeDouble(time);}
  void setBeginTime(TDatime *time){fMinTime = SimTools::getTimeDouble(time);}
  void setEndTime(Int_t year, Int_t month, Int_t day, Int_t hour,Int_t minute, Int_t sec);
  void setEndTime(TTimeStamp *time){fMaxTime = SimTools::getTimeDouble(time);}
  void setEndTime(TDatime *time){fMaxTime = SimTools::getTimeDouble(time);}
  Double_t getBeginTime(){return fMinTime;}
  Double_t getEndTime(){return fMaxTime;}  

  Double_t getTimeStep(){return fTimeStep;}
  void setTimeStep(Double_t ts){fTimeStep = ts;}

  //Particle Masses
  void setRecoilMass(Double_t mass = 19e6){fRecoilMass = mass;}
  void findRecoilMass(){fRecoilMass = SimTools::findParticleMass(fRecoilParticle);}
  Double_t getRecoilMass(){return fRecoilMass;}  
  void setRecoilA(Double_t A = 19){fRecoilA = A;}
  void findRecoilA(){fRecoilA = SimTools::findParticleA(fRecoilParticle);}
   Double_t getRecoilA(){return fRecoilA;}
  void setProjectileMass(Double_t mass = 1.e6){fProjMass = mass;}
  void findProjectileMass(){fProjMass = SimTools::findParticleMass(fProjParticle);}
  Double_t getProjectileMass(){return fProjMass;}


  //Collimation maximum angle
  void setThetaMax(Double_t theta){fThetaMax = theta;}
  Double_t getThetaMax(){return fThetaMax;}

  //Particle Names
  void setRecoilParticle(TString name){fRecoilParticle = name;}
  
  void setProjectileParticle(TString name){fProjParticle = name;}
  TString getRecoilParticle(){return fRecoilParticle;}
  TString getProjectileParticle(){return fProjParticle;}

  //Location
  void setRunType(TString name){fRunType = name;}
  void setConfigType(TString name){fConfigType = name;}
  void setLocation(TString loc);
  void setLatLong(Double_t lat, Double_t lon){fLat = lat; fLong = lon;}
  Double_t getLatitude(){return fLat;}
  Double_t getLongitude(){return fLong;}
  void setRotationAngle(Double_t ang){fRotationAngle = ang;}
  Double_t getRotationAngle(){return fRotationAngle;}
  void setTopOrBottomTPC(Bool_t tb){fTopBottom = tb;}//1 = top, 0 = bottom
  Bool_t getTopOrBottomTPC(){return fTopBottom;}

  //Options
  void setRecoilType(TString type){type.ToLower(); fRecoilType = type;}
  void setEnergyOption(TString opt){opt.ToLower(); fEnOption = opt;}
  void setPositionOption(TString opt){opt.ToLower(); fPosOption = opt;}
  void setDirectionOption(TString opt){opt.ToLower(); fDirOption = opt;}
  void setTimeOption(TString opt){opt.ToLower(); fTimeOption = opt;}
  void setSpecialOption(TString opt){opt.ToLower(); fSpecialOption = opt;}
  void setTheoryEnergyOption(TString opt){opt.ToLower(); fTheoryEnergyOption = opt;}
  void setDecayChain(TString s,TString key); 
  TString getConfigType(){return fConfigType;}
  TString getRecoilType(){return fRecoilType;}
  TString getRunType(){return fRunType;}
  TString getEnergyOption(){return fEnOption;}
  TString getPositionOption(){return fPosOption;}
  TString getDirectionOption(){return fDirOption;}
  TString getTimeOption(){return fTimeOption;}
  TString getSpecialOption(){return fSpecialOption;}
   TString getTheoryEnergyOption(){return fTheoryEnergyOption;}

  //ENDF files
  void setENDFfiles();
  void setSpectrum(TString filename){fSpecName = filename;}
  void setScattering(TString filename){fScatDCSname = filename;}
  void setCrossSection(TString filename){fScatCSname = filename;}
  MaxCamENDF *spectrum(){return fSpectrum;}
  MaxCamENDF *scattering(){return fScatteringDCS;}
  MaxCamENDF *crossSection(){return fScatteringCS;}
  TString getSpectrum(){return fSpecName;}
  TString getScattering(){return fScatDCSname;}
  TString getCrossSection(){return fScatCSname;}

  TVector3 getWimpWindDirection();
 
  //Decay Chain
  //


  //Generate projectiles & recoils
  void generateTime();
  void generatePosition();
  void generateDirection();
  void generateEnergy();
  void generateDirectionAndEnergyDoubleAlpha();
  void generateWimpProjectile();
   void generateWimpRecoil();
  void generateRandomRecoil();
  void generateENDF();

  //Generate the 3 types of particles included so far
  void generateWIMP();
  void generateAlphaType();
  void generateDoubleAlphaType();
  void generateNeutronType();

  //Pick the type of recoil and generate events
  void generateRecoil();



private:

  TRandom3 *fRandom;             //Random number generator

  TVector3 *fSourcePos;          //Source Position in mm
  TVector3 *fSourceDir;          //Source direction (unit vector)
  TVector3 *fRecoilPos;          //Recoil Position in mm
  TVector3 *fProjVector;         //Projectile direction (unit vector)
  TVector3 *fRecVector;          //Recoil direction (unit vector)
  TVector3 *fRecVector1;          //Recoil direction (unit vector)
  TVector3 *fRecVector2;          //Recoil direction (unit vector)
  Double_t fProjEn;              //Projectile energy in keV
  Double_t fProjEn2;              //Projectile energy in keV  
  Double_t fProjGaussianEn;      //Projectile energy in keV
  Double_t fProjGaussianEnSpread;//Projectile energy spread in kev
  Double_t fRecoilEn;            //Recoil energy in keV
  Double_t fRecoilEn2;            //Recoil energy in keV
  Double_t fProjMass;            //Projectile mass in keV
  Double_t fRecoilMass;          //Recoil Mass in keV
  Double_t fRecoilA;             //Recoil A
  Int_t fSeq;
  Double_t fTime;                //Time in days from 1200 31/12/1999
  Double_t fTimeSet;                //Time in days from 1200 31/12/1999
  Double_t fThetaMax;            //Maximum extent of collimation in degrees
  Double_t fXMax;                //Maximum recoil x (mm)
  Double_t fXMin;                //Minumum recoil x (mm)
  Double_t fYMin;                //Minimum recoil y (mm)
  Double_t fYMax;                //Maximum recoil y (mm)
  Double_t fZMin;                //Minimum recoil z (mm)
  Double_t fZMax;                //Maximum recoil z (mm)
  Double_t fMinEn;               //Minimum energy (keV)
  Double_t fMaxEn;               //Maximum energy (keV)
  Double_t fMinTime;             //Earliest allowed time (days)
  Double_t fMaxTime;             //Latest allowed time (days)
  Double_t fTimeStep;	         //Time between events for "series" time option runs (s)
  Double_t fCosCygnus;            //Angle with Cygnus
  Double_t fLat;                 //Latitude (degrees)
  Double_t fLong;                //Longitude (degrees)
  Double_t fRotationAngle;       //Degrees west (counterclockwise) of north of the +y axis
  Bool_t fTopBottom;             //If 1, top TPC, looking down 2: bottom TPC (looking up)

  TString fRecoilParticle;       //Name of recoil particle
  TString fProjParticle;         //Name of projectile particle
  TString fRunType;
  TString fConfigType;
  TString fRecoilType;           //Recoil type: Alpha, neutron or wimp
  TString fEnOption;             //Energy Options: Fix, Random
  TString fPosOption;            //Position Options: Fix, Random
  TString fDirOption;            //Direction Options: Fix, Source, Isotropic, Collimate
  TString fTimeOption;           //Time Options: Fix, Random, Current
  TString fSpecialOption;        //Special options: wimp or endf run
  TString fTheoryEnergyOption;   //Theoretical Energy Options: flat, vesc, infinity
  
  

  DmtpcDecayChain * fDecay; 
  SimRings * fRing; 
  TH1* fEnSpec;			 //Projectile energy spectrum
  MaxCamENDF *fSpectrum;         //ENDF Neutron emission spectrum
  MaxCamENDF *fScatteringDCS;    //Scattering total cross section
  MaxCamENDF *fScatteringCS;     //Scattering cross section
  TString fSpecName;             //Neutron spectrum filename
  TString fScatDCSname;          //Scattering total cross section filename
  TString fScatCSname;           //Scattering cross section filename

   TF1* fdRdEInf;                //Recoil Energy Spectrum with vesc=infinity
   TF1* fdRdEVesc;               //Recoil Energy Spectrum with vesc=600
   TF1* fdRdTh;                  //Recoil Angle Spectrum with set energy
//   TTree *t1;
//   TFile *f;

  Int_t nentries;


  ClassDef(ParticleGenerator,2);
};

#endif
