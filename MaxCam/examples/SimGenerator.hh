#ifndef SIMGENERATOR_HH
#define SIMGENERATOR_HH
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
////Includes various particle generators to be called                  //
////by makeSimulations.cc                                              //
/////////////////////////////////////////////////////////////////////////
#include "../MaxCamImageTools.hh"
#include "../MaxCamWIMP.hh"
#include "../McWimp.hh"
#include "../MaxCamENDF.hh"
#include "../MaxCamTwoBodyKinematics.hh"

#include "TH2.h"
#include "TSystem.h"
#include "TString.h"
#include "TMath.h"
#include "TVector3.h"

#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

class SimGenerator {

public:

  SimGenerator();
  virtual ~SimGenerator();

  //Source position vector
  void setSourcePosition(TVector3* pos) {_sourcePos = pos;}
  void setSourcePosition(TVector3 pos) {
    _sourcePos->SetXYZ(pos.X(),pos.Y(),pos.Z());
  }
  void setSourcePosition(double x, double y, double z){
    _sourcePos->SetXYZ(x,y,z);
  }
  TVector3* getSourcePosition(){return _sourcePos;}
  TVector3 sourcePosition();

  //Recoil position vector
  void setRecoilPosition(TVector3* pos) {_recoilPos = pos;}
  void setRecoilPosition(TVector3 pos) {
    _recoilPos->SetXYZ(pos.X(),pos.Y(),pos.Z());
  }
  void setRecoilPosition(double x, double y, double z){
    _recoilPos->SetXYZ(x,y,z);
  }
  TVector3* getRecoilPosition(){return _recoilPos;}
  TVector3 recoilPosition();
  
  //Recoil position limits
  void setPositionMaxValues(double x, double y, double z);
  void setPositionMinValues(double x, double y, double z);

  //Projectile and recoil energies
  double getProjectileEnergy(){return _projEn;}
  void setProjectileEnergy(double en){_projEn = en;}
  double getRecoilEnergy(){return _recoilEn; }
  void setRecoilEnergy(double en){_recoilEn = en;}

  //Particle type
  void setType(TString type){_type = type;}
  TString getType(){return _type;}

  //Projectile direction vector
  void setProjectileVector(TVector3* vect);
  void setProjectileVector(TVector3 vect){
    _projVector->SetMagThetaPhi(1,vect.Theta(),vect.Phi());
  }
  void setProjectileVector(double theta, double phi){
    _projVector->SetMagThetaPhi(1,theta,phi);
  }
  void setProjectileVector(double x, double y, double z);
  TVector3* getProjectileVector(){return _projVector;}
  TVector3 projectileVector();

  //Recoil direction vector
  void setRecoilVector(TVector3* vect);
  void setRecoilVector(TVector3 vect){
    _recVector->SetMagThetaPhi(1,vect.Theta(),vect.Phi());
  }
  void setRecoilVector(double theta, double phi){
    _recVector->SetMagThetaPhi(1,theta,phi);
  }
  void setRecoilVector(double x, double y, double z);
  TVector3* getRecoilVector(){return _recVector;}
  TVector3 recoilVector();

  void setMaxEnergy(double en){_maxEn = en;}
  void setMinEnergy(double en){_minEn = en;}
  double getMaxEnergy(){return _maxEn;}
  double getMinEnergy(){return _minEn;}

  void setTime(double time){_time = time;}
  double getTime(){return _time;}

  void setRecoilMass(double mass = 19e6){_recoilMass = mass;}
  double getRecoilMass(){return _recoilMass;}

  void setProjectileMass(double mass = 1.e6){_projMass = mass;}
  double getProjectileMass(){return _projMass;}

  void setThetaMax(double theta){_thetaMax = theta;}
  double getThetaMax(){return _thetaMax;}

  TVector3 wimpGun();

  void generateFixedAlpha();
  void generateCf252Neutron(TString recoiltype ="");
  void generateWimp();
  void generateRandomAlpha();
  void generateRandomNeutron();
  bool generateRecoil(TString recoilpart = "");
  bool generateRecoil(TString type, TString type2);

private:

  TVector3 *_sourcePos;
  //Source Position applies to fixed alpha and Cf-252 events
  TVector3 *_recoilPos;
  TVector3 *_projVector;
  TVector3 *_recVector;
  double _projEn;
  //projectile energy not included in wimp events
  double _recoilEn;
  //set recoil energy for fixed-source alpha events
  double _projMass;
  //projectile mass must be set for wimp events
  double _recoilMass;
  //recoil mass must be set for neutron and wimp events
  double _time;
  double _thetaMax;
  //ThetaMax set for fixed alpha runs (collimation in degrees)
  double _xMax;
  double _xMin;
  double _yMin;
  double _yMax;
  double _zMin;
  double _zMax;
  double _minEn;
  //Min Energy applies to wimps, neutrons and random alphas.
  double _maxEn;
  //Max Energy applies to neutrons and random alphas
  TString _type;
  //Types supported: "ar": random alpha, "af": fixed alpha, "nr": random neutron, "cf": Californium neutron, "w": wimp

};



SimGenerator::SimGenerator(){

  _recoilMass = 19e6;
  _projMass = 1e6;
  _time = 365.*gRandom->Rndm();
  _xMin = 0;
  _yMin = 0;
  _zMin = 0;
  _xMax = 0.143*1024;
  _yMax = 0.143*1024;
  _zMax = 250;
  _thetaMax = 5;
  _recoilEn = 0;
  _projEn = 1;
  _minEn = 0;
  _maxEn = 1000;
  _sourcePos = new TVector3(1,0,0);
  _recoilPos = new TVector3(1,0,0);
  _projVector = new TVector3(1,0,0);
  _recVector = new TVector3(1,0,0);

}
SimGenerator::~SimGenerator(){

  delete _sourcePos;
  delete _recoilPos;
  delete _projVector;
  delete _recVector;

}

TVector3 
SimGenerator::sourcePosition(){

  TVector3 vect(_sourcePos->X(),_sourcePos->Y(),_sourcePos->Z());
  return vect;

}

TVector3 
SimGenerator::recoilPosition(){

  TVector3 vect(_recoilPos->X(),_recoilPos->Y(),_recoilPos->Z());
  return vect;

}

void
SimGenerator::setProjectileVector(TVector3 *vect){
  vect->SetMag(1);
  _projVector = vect;
}

TVector3 
SimGenerator::projectileVector(){

  TVector3 vect(_projVector->X(),_projVector->Y(),_projVector->Z());
  return vect;

}

void
SimGenerator::setRecoilVector(TVector3 *vect){
  vect->SetMag(1);
  _recVector = vect;
}

TVector3 
SimGenerator::recoilVector(){

  TVector3 vect(_recVector->X(),_recVector->Y(),_recVector->Z());
  return vect;

}

void 
SimGenerator::setProjectileVector(double x, double y, double z){
  _projVector->SetXYZ(x,y,z);
  _projVector->SetMag(1);
}

void
SimGenerator::setRecoilVector(double x, double y, double z){
  _recVector->SetXYZ(x,y,z);
  _recVector->SetMag(1);
}


void 
SimGenerator::setPositionMaxValues(double x, double y, double z){
  _xMax = x;
  _yMax = y;
  _zMax = z;
}

void 
SimGenerator::setPositionMinValues(double x, double y, double z){
  _xMin = x;
  _yMin = y;
  _zMin = z;
}



bool
SimGenerator::generateRecoil(TString type,TString type2){

  setType(type);
  bool isDone = generateRecoil(type2);
  return isDone;
}

bool
SimGenerator::generateRecoil(TString recoilpart){

  recoilpart.ToLower();
  if (recoilpart == "xe") setRecoilMass(132e6);
  else if (recoilpart == "c") setRecoilMass(12e6);
  else if (recoilpart == "he") setRecoilMass(4e6);
  else setRecoilMass(19e6);

  bool isGenerated = true;
  //Fixed gun Alpha
  if (getType() == "af")generateFixedAlpha();
  //Random Alpha
  else if (getType() == "ar")generateRandomAlpha();
  //Neutron
  else if (getType() == "nr")generateRandomNeutron();
  //Wimp
  else if (getType() == "w")generateWimp();
  //Cf-252
  else if (getType() == "cf" || getType() == "nf")generateCf252Neutron(recoilpart);
  else {
    cout << "Type not found. No recoil generated."<<endl;
    isGenerated = false;
  }

  return isGenerated;
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
////Generate Alpha tracks within range of energies (emin and emax)          //
////Random direction                                                        //
//////////////////////////////////////////////////////////////////////////////

void 
SimGenerator::generateRandomAlpha(){

  setProjectileMass(4e6);
  setRecoilMass(4e6);
  double x = _xMin + (_xMax-_xMin)*gRandom->Rndm();
  double y = _yMin + (_yMax-_yMin)*gRandom->Rndm();
  double z = _zMin + (_zMax-_zMin)*gRandom->Rndm();
  setRecoilPosition(x,y,z);

  double phi = 2*TMath::Pi()*gRandom->Rndm() - TMath::Pi();
  double theta = 2.*asin(sqrt(gRandom->Rndm()));//Inverse cumulative dist. function
  setProjectileVector(theta,phi);
  setRecoilVector(theta,phi);

  double recoilEn = _minEn + (_maxEn-_minEn)*gRandom->Rndm();
  setRecoilEnergy(recoilEn);
  setProjectileEnergy(recoilEn);

}//End generateAlpha()

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
////Generate Alpha tracks with known energy and start position              //
////                                                                        //
//////////////////////////////////////////////////////////////////////////////

void 
SimGenerator::generateFixedAlpha(){

  setProjectileMass(4.e6);
  setRecoilMass(4.e6);

  double thetaMax = _thetaMax*TMath::Pi() / 180.; //Collimated to deltaTheta degrees
  //Phi and theta measured from y-axis: shooting alphas roughly along y;
  //Theta distribution goes like sin(theta)
  //Want to cut off at theta max
  //Probability density p(theta) = sin(theta)/Integral(sin(theta), 0 < theta < thetaMax)
  //sin(theta) = p(theta)*(1-cos(thetaMax)) = p * 2*sin^2(thetaMax/2)
  //Cumulative probability function
  //1 - sin(theta) = P * 2 * sin^2(thetaMax/2) = 2 * sin^2(theta/2), P is random from 0 to 1
  //==> theta = 2*arcsin{sqrt(P) * sin(thetaMax/2)}

  double phi = 2*TMath::Pi()*gRandom->Rndm() - TMath::Pi();
  double theta = 2.*asin(sqrt(gRandom->Rndm())*sin(thetaMax/2.));//Inverse cumulative dist. function
  TVector3 projVector = projectileVector().Unit();

  TVector3 orthoVector = projVector.Orthogonal().Unit();
  orthoVector.Rotate(phi,projVector);

  setRecoilVector(projVector*cos(theta)+orthoVector*sin(theta));
  setProjectileVector(getRecoilVector());
  setRecoilPosition(getSourcePosition());

}//End generateAlpha()

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
////Generate neutrons with random energies within emin and emax            //
////Random neutron directions, so yields random recoil direction           //
/////////////////////////////////////////////////////////////////////////////

void 
SimGenerator::generateRandomNeutron(){

  setProjectileMass(1.e6);

  //Generate interaction position
  double x = _xMin + (_xMax-_xMin)*gRandom->Rndm();
  double y = _yMin + (_yMax-_yMin)*gRandom->Rndm();
  double z = _zMin + (_zMax-_zMin)*gRandom->Rndm();
  
  double recoilEn, projEn;
  TVector3 projVector(1,0,0);
  
  //Generate neutron direction
  double phi = 2*TMath::Pi()*gRandom->Rndm();
  double theta = 2.*asin(sqrt(gRandom->Rndm()));//Inverse cumulative dist. function

  //Generate neutron energy
  bool done = false;
  while(!done){
    projEn = _minEn + (_maxEn-_minEn)*gRandom->Rndm();
    //Check if any recoil is possible in [minEn,maxEn]
    double cosmin = MaxCamTwoBodyKinematics::
                    calcCosRecoilFromRecoilEnergy(_minEn,projEn,_recoilMass,_projMass);
    if (cosmin >= -1 && cosmin <= 1) done = true;  
  }

  projVector.SetMagThetaPhi(1.,theta,phi);
  TVector3 orthoVector = projVector.Orthogonal().Unit();

  //Generate recoil energy and angle
  done = false;
  double cosRecoil = -2;
  while(!done){
    recoilEn = _minEn + (projEn-_minEn)*gRandom->Rndm();
    //See MaxCamTwoBodyKinematics::calcCosRecoilFromRecoilEnergy
    cosRecoil = MaxCamTwoBodyKinematics::
                calcCosRecoilFromRecoilEnergy(recoilEn,projEn,_recoilMass,_projMass);
    //Test for valid cosine
    if (cosRecoil >= -1 && cosRecoil <= 1) done = true;
  }
  double sinRecoil = sqrt(1.-cosRecoil*cosRecoil);
  //cosRecoil is cosine of angle (theta) between neutron vector and recoil vector
  //Get random phi (still in neutron direction coords)
  double phi1 = 2.*TMath::Pi()*gRandom->Rndm();
  orthoVector.Rotate(phi1,projVector);
  //Set recoil vector
  TVector3 recVector = cosRecoil*projVector.Unit() + sinRecoil*orthoVector.Unit();

  setProjectileEnergy(projEn);
  setRecoilEnergy(recoilEn);
  setRecoilPosition(x,y,z);
  setRecoilVector(recVector);
  setProjectileVector(projVector);

}//End generateNeutron()

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
////Generate Cf252 neutrons                                            //
/////////////////////////////////////////////////////////////////////////

void 
SimGenerator::generateCf252Neutron(TString recoiltype){

  setProjectileMass(1.e6);
  //Source at known position sourcePos
  //neutrons generated at random position within window
  double x = gRandom->Rndm()*(_xMax - _xMin) + _xMin;
  double y = gRandom->Rndm()*(_yMax - _yMin) + _yMin;
  double z = gRandom->Rndm()*(_zMax - _zMin) + _zMin;
  //Neutron direction gotten from difference in source and interaction positions.
  TVector3 projVector(x-_sourcePos->X(),y-_sourcePos->Y(),z-_sourcePos->Z());
  projVector = projVector.Unit();
  TVector3 orth = projVector.Orthogonal().Unit();

  //Now extract energy from known radiation spectrum
  //Also check to be sure a recoil is possible within desired range.
  MaxCamENDF *neutronEnergy = new MaxCamENDF("ENDF_Cf-252_n_spectrum","fission");
  recoiltype.ToLower();
  MaxCamENDF *scatteringDCS, *scatteringCS;
  //  if (recoiltype == "xe"){
  //    scatteringDCS = new MaxCamENDF("ENDF_DCS_n_on_132Xe","elastic scattering");
  //    scatteringCS = new MaxCamENDF("ENDF_CS_n_on_132Xe","cs");
  //  }
  //  else if (recoiltype == "c"){
  //scatteringDCS = new MaxCamENDF("ENDF_DCS_n_on_12C","elastic scattering");
  // scatteringCS = new MaxCamENDF("ENDF_CS_n_on_12C","cs");
  // }
  // else if (recoiltype == "he"){
  //  scatteringDCS = new MaxCamENDF("ENDF_DCS_n_on_4He","elastic scattering");
  //  scatteringCS = new MaxCamENDF("ENDF_CS_n_on_4He","cs");
  // }
  //else{
    scatteringDCS = new MaxCamENDF("ENDF_DCS_n_on_19F","elastic scattering");
    scatteringCS = new MaxCamENDF("ENDF_CS_n_on_19F","cs");
    //}
  double projEn;
  bool done = false;
  if (getType() == "cf"){
    while(!done){
      projEn = neutronEnergy->generateEnergy(_minEn,_maxEn);
      if (!scatteringCS->acceptEnergy(projEn)) continue; //check cross section
      //Check if any recoil is possible in [minEn,maxEn]
      double cosmin = MaxCamTwoBodyKinematics::
	calcCosRecoilFromRecoilEnergy(_minEn,projEn,_recoilMass,_projMass);
      if (cosmin >= -1 && cosmin <= 1) done = true;  
    }
  }
  else {
    projEn = _projEn;
    double cosmin = MaxCamTwoBodyKinematics::
      calcCosRecoilFromRecoilEnergy(_minEn,projEn,_recoilMass,_projMass);
    if (!(cosmin >= -1 && cosmin <= 1)){
      cerr<<"Minimum energy higher than maximum recoil energy for given projectile energy. Quitting.";
      assert(0);
    }

  }

  //Generate recoil energy and angle
  done = false;
  double cosRecoil = -2;
  double recoilEn;
  //  while(!done){
  //    recoilEn = _minEn + (projEn - _minEn)*gRandom->Rndm();
  //See MaxCamTwoBodyKinematics::calcCosRecoilFromRecoilEnergy
  //    cosRecoil = MaxCamTwoBodyKinematics::
  //                calcCosRecoilFromRecoilEnergy(recoilEn,projEn,_recoilMass,_projMass);
  //Test for valid cosine
  //    if (cosRecoil >= -1 && cosRecoil <= 1) done = true;
  //  }
  while(!done){
    double cosScatterCMS = scatteringDCS->generateCosAngleCMS( projEn);
    recoilEn = MaxCamTwoBodyKinematics::calcRecoilEnergyFromCosScatterCMS(cosScatterCMS,projEn,_recoilMass,_projMass);
    if (recoilEn >= _minEn && recoilEn <= _maxEn) done = true;
  }
  cosRecoil = MaxCamTwoBodyKinematics::calcCosRecoilFromRecoilEnergy(recoilEn,projEn,_recoilMass,_projMass);
  double sinRecoil = sqrt(1.-cosRecoil*cosRecoil);
  //cosRecoil is cosine of angle (theta) between neutron vector and recoil vector
  //Get random phi (still in neutron direction coords)
  double phi1 = 2.*TMath::Pi()*gRandom->Rndm();
  orth.Rotate(phi1,projVector);
  //Set recoil vector
  TVector3 recVector = cosRecoil*projVector.Unit() + sinRecoil*orth.Unit();

  setRecoilVector(recVector);
  setRecoilPosition(x,y,z);
  setRecoilEnergy(recoilEn);
  setProjectileEnergy(projEn);
  setProjectileVector(projVector);

  delete neutronEnergy;
  delete scatteringDCS;
  //  delete scatteringCS;
}//end generateCf252

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
////Generate WIMP using McWimp class                                   //
/////////////////////////////////////////////////////////////////////////
void 
SimGenerator::generateWimp(){

  McWimp *simwimp = new McWimp();
  simwimp->setTime(getTime());
  simwimp->setWimpMass(getProjectileMass());
  simwimp->setOpt("lv");
  simwimp->setEmin(getMinEnergy());
  simwimp->genWimp();

  double x = _xMin + (_xMax-_xMin)*gRandom->Rndm();
  double y = _yMin + (_yMax-_yMin)*gRandom->Rndm();
  double z = _zMin + (_zMax-_zMin)*gRandom->Rndm();
  setRecoilPosition(x,y,z);
  setRecoilVector(simwimp->getRecoilV());
  setProjectileVector(simwimp->getWimpV());
  setProjectileEnergy(simwimp->getWimpV().Mag2()/2*_projMass);
  setRecoilEnergy(simwimp->getRecoilE());

}//end generateWimp()

// void
// SimGenerator::generateWimp2(){

//   double lat = 44.355;
//   double long = -103.772;
//   TVector3 wimpV(-240,0,0);//in km/s
//   double c = 3e5;//in km/s
//   TVector3 vEarth = wimpGun();
//   double vx, vy, bz;
//   TVector3 wimpVLab;
//   double vEsc = 600*10;//max wimp velocity
//   double v0 = 230;//average WIMP velocity

//   bool done = false;
//   while(!done){

//     vx = (gRandom->Rndm()-0.5)*vEsc;
//     vy = (gRandom->Rndm()-0.5)*vEsc;
//     vz = (gRandom->Rndm()-0.5)*vEsc;
//     wimpVLab = TVector3(vx,vy,vz);

//     double rnd = gRandom->Rndm();
//     fmax = exp(-(vEarth+wimpVLab).Mag2()/(v0*v0));
//     if (fmax < rnd) continue;
//     double projEn = 0.5 * _projMass*wimpVLab.Mag2()/(c*c);
//     double mRatio = _recoilMass / _projMass;
//     double cosMax = MaxCamTwoBodyKinematics::
//                     calcCosRecoilFromRecoilEnergy(_minEn,projEn,_recoilMass,_projMass);
//     if (abs(cosrecoil) <= 1 ) done = true;
  
//   }

//   _projEn = projEn;
//   _projVect->SetXYZ(vx,vy,vz);
//   _projVect->SetMag(1);

//   //Generate

// }

// TVector3 
// SimGenerator::wimpGun() {
//   //Copied from McWimp
//     // Made all parameters into data members 7/17
//   //DUSEL (Lead, SD)
//   //    double lat = 44.355;
//   //    double lon = -103.772;
//   //Cambridge, MA
//   double lat = 42+22/60. + 25/3600.;
//   double lon = -71 - 6/60. - 38/3600.;
//   //WIPP (Carlsbad, NM)
//   //double lat = 32 + 24/60.+43/3600.;
//   //double lon = -104 - 14/60.-11/3600.;


//     TVector3 vS(9,242,7); // Sun's velocity relative to WIMP cloud
//     //// Set constants to calculate Earth's velocity

//     double
//       sE_avg = 29.79, // Earth's average orbital speed
//       e = 0.016722,   // ellipticity of Earth's orbit
//       L0 = 13 *TMath::Pi()/180,  // longitude of orbit's minor axis, +- 1 degree uncertainty
    
//       Bx = -5.5303 *TMath::Pi()/180,  //lat & lon of ecliptic's xyz axes - from Lewin-Smith
//       Lx = 266.141 *TMath::Pi()/180,
//       By = 59.575 *TMath::Pi()/180,
//       Ly = -13.3485 *TMath::Pi()/180,
//       Bz = 29.812 *TMath::Pi()/180,
//       Lz = 179.3212 *TMath::Pi()/180;


//     double l = 280.460 + 0.9856474*_time;
//     double g = (357.528 + 0.9856003*_time)*TMath::Pi()/180;

//     double L = ( l + 1.915*sin(g) + 0.020*sin(2*g) ) *TMath::Pi()/180;
//     double sE = sE_avg*( 1-e*sin(L-L0) ); // Earth's speed at position L


//     double vE_x = sE * cos(Bx) * sin(L-Lx);
//     double vE_y = sE * cos(By) * sin(L-Ly);
//     double vE_z = sE * cos(Bz) * sin(L-Lz);

//     TVector3 vE(vE_x, vE_y, vE_z); // Earth's velocity relative to Sun
//     TVector3 vT( vS + vE); // Earth's velocity relative to WIMP cloud
//     //    if (opt=="g") {return vT;} //galactic
//     // If do not want galactic coord, convert to equatorial
//     // More accurate numbers may be available - check to update soon
//     vT.Rotate(62.8716639*TMath::Pi()/180, TVector3(.83867,.54464, 0));
//     vT.RotateZ(249.75*TMath::Pi()/180);
//     //    if (opt=="e") {return vT;} //equatorial
//     double side = 1.00273790935;    // Number of sidereal days per day
//     double t0 = (lon + 10.0)*TMath::Pi()/180; // Relative global long. to eq. long. set to 0 at t=0
//     vT.RotateZ(-(t*side)*2.0*pi - t0);
//     vT.RotateX((lat-90)*TMath::Pi()/180);
//     //    if (opt=="lv") {return vT;} //lab, oriented vertically
//     return vT;
//     //    vT.RotateX(90*TMath::Pi()/180);
//     //    if (opt=="lh") {return vT;} //lab, oriented horizontally    
// }
#endif
