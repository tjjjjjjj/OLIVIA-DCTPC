#ifndef MC_WIMP_HH
#define MC_WIMP_HH

#include "TROOT.h"
#include "TRandom.h"

#include "TF1.h"
#include "TVector3.h"

//#include "Wimp.hh"
#include "DetObj.hh"

#include <iostream>

using namespace std;


/** DEPRECATED!  \deprecated Use DmtpcTheory and DmtpcAstro */

class McWimp {


public:

  // McWimp();

  McWimp(double wM = 100e6,
	 TVector3 wV = TVector3(-240,0,0),
	 double t = 0, 
	 double rE = 0, 
	 TVector3 sP = TVector3(12.5,0,0), 
	 TVector3 rV = TVector3(-100,0,0))

  {
    _wimpMass = wM;
    _wimpV = wV;
    _time = t; 
    
    _recoilE = rE;  

    _startPos = sP;
    _recoilV = rV;

    _type = "w";

    _fxn = new TF1("","",0,0);

    // set detector config. currently set for 5 CCD tiling model 7/16
    DetObj* tempDet = new DetObj(5,10,(400.0/256.0),(400.0/256.0));
    _det = tempDet;
    _det->addCCD(-200,200,256,256,0,1);
    _det->addCCD(200,200,256,256,0,2);
    _det->addCCD(-200,-200,256,256,0,3);
    _det->addCCD(200,-200,256,256,0,4);
    _det->addCCD(0,0,256,256,1,5);

    //default parameters for genWimp and cosmicGun
    _opt = "e";
    _Emin = 0;
    _lat = 44.355;
    _lon = -103.772;
    _ccd = 0;
  }



  void setAlpha() {
    setWimpMass(4e6);
    _type="a";

    _wimpV = TVector3(0,0,0);
    _time = 0;    
    _opt = "g"; // set if you want to know direction of cygnus even while running backgrounds
    _Emin = 0;
    _lat = 0;
    _lon = 0;
    _ccd = 0;
  }


  void setWimp() {
    _wimpMass = 100e6;
    _type = "w";

    _wimpV = TVector3(-240,0,0);
    _opt = "e";
    _lat = 44.355;
    _lon = -103.772;
    _ccd = 0;
  }





  //////////////////////////////
  //// Data Member Handlers
  //////////////////////////////

  double getWimpMass() {return _wimpMass;}
  void setWimpMass(double wM) { _wimpMass = wM;}
  
  TVector3 getWimpV() {return _wimpV;}
  void setWimpV(TVector3 wV) { _wimpV = wV;}
  
  double getTime() {return _time;}
  void setTime(double t) { _time = t;}
  
  double getRecoilE() {return _recoilE;}
  void setRecoilE(double rE) { _recoilE = rE;}
  
  TVector3 getStartPos() {return _startPos;}
  void setStartPos(TVector3 sP) { _startPos = sP;}
  
  TVector3 getRecoilV() {return _recoilV;}
  void setRecoilV(TVector3 rV) { _recoilV = rV;}
  
  TH2F getSimRecoil() {return _simRecoil;}
  void setSimRecoil(TH2F sR) { _simRecoil = sR;}

  string getType() {return _type;}
  // void setType(string t) {_type = t;}
  
  DetObj* getDet() {return _det;}            // Det handlers 7/16
  void setDet(DetObj* det) {_det = det;}

  TF1* getFxn() {return _fxn;}
  void setFxn(TF1* f) {_fxn = f;}

  void setOpt(string opt) {_opt = opt;}
  void setEmin(double Em) {_Emin = Em;}
  void setLat(double lat) {_lat = lat;}
  void setLon(double lon) {_lon = lon;}
  void setCCD(int ccd) {
    _ccd = ccd;
    _det->imageCCDs(ccd);
  }

  void startInside() {_det->startInside();}
  void startAnywhere() {_det->startAnywhere();}
  
  




  ////////////////////////////////
  //// Data Analysis
  ////////////////////////////////

  
  double Lxy() {
    TVector3 rV = getRecoilV();
    double lxy = sqrt( rV.X()*rV.X() + rV.Y()*rV.Y() ) * getRecoilE();
    return lxy; 
  }


  double Lz() {
    TVector3 rV = getRecoilV();
    double lz = rV.Z() * getRecoilE();
    return lz;
  }


  double thetaRecoil() {
    double cR = cosRecoil();
    return acos(cR);
  }


  double cosRecoil() {
    TVector3 rV = getRecoilV().Unit();
    TVector3 wV = getWimpV().Unit();
    double cR = rV * wV;
    return cR;
  }


  double sinRecoil() {
    TVector3 rV = getRecoilV().Unit();
    TVector3 wV = getWimpV().Unit();
    TVector3 crossV = rV.Cross(wV);
    return crossV.Mag();
  }


  double thetaXY() {
    TVector3 rV = (getRecoilV()-getRecoilV().Z()).Unit();
    TVector3 wV = (getWimpV()-getWimpV().Z()).Unit();
    //return abs((int) (rV.Phi()-wV.Phi())*1000)/1000.0;
    double cR = rV * wV;
    return acos(cR);
  }

  
  double sinXY() {
    TVector3 rV = (getRecoilV()-getRecoilV().Z()).Unit();
    TVector3 wV = (getWimpV()-getWimpV().Z()).Unit();
    TVector3 crossV = rV.Cross(wV);
    return crossV.Z();
  }


  double sinZ() { // This is not very useful...
    TVector3 rV = getRecoilV().Unit();
    TVector3 wV = getWimpV().Unit();
    return rV.Theta()-wV.Theta();
  }
    

  double thetaCy() {
    TVector3 rV = (getRecoilV()-getRecoilV().Z()).Unit();
    TVector3 cT = (cosmicGun()-cosmicGun().Z()).Unit();    
    // double rT = getRecoilV().Phi();
    // double cT = cosmicGun().Phi();
    // return abs((int) (rT - cT)*1000)/1000.00;
    double cR = rV * cT;
    return acos(cR);
  }






  void simulate() {
    if (_type == "w") {
      genWimp();
    }
    if (_type == "a") {
      genAlpha();
    }
  }



  /////////////////////////////////////////////////////////////////////
  ////  genWimp()
  /////////////////////////////////////////////////////////////////////
  
  void genWimp() { 
    //  Made all parameters data members 7/16

    if(_type == "a") {
      cout << "Can't run genWimp on alphas" << endl;
      return;
    }

    double t = getTime();
    string opt = _opt;   // May want "get" fxns for these
    double Emin = _Emin;


    double wimpM = getWimpMass();
    double recoilM = 19e6;  // Set to fluorine nucleus
    // double R = 4*wimpM*recoilM/(wimpM+recoilM)/(wimpM+recoilM);
    
    double recoilE;
    double Emax = 300; // max WIMP energy
    TVector3 vEarth = cosmicGun();  
    
    double c = 3e5;  // speed of light
    
    
    //// Generate WIMP velocity in lab frame
    
    double vx,vy,vz;
    TVector3 wimpVLab; // velocity of WIMP in lab frame
    
    double vEsc = 600*10; // limit for velocity
    double v0 = 230;      // average WIMP speed
    
    double rnd, fmax;
    double cosRecoil;
    
    while (1) {
      
      vx = (gRandom->Rndm()-0.5)*vEsc;
      vy = (gRandom->Rndm()-0.5)*vEsc;
      vz = (gRandom->Rndm()-0.5)*vEsc;
      wimpVLab = TVector3(vx,vy,vz);
         
      
      //// Check if WIMP speed fits Max_Boltz Distr.
      rnd = gRandom->Rndm();
      fmax = exp( - (vEarth+wimpVLab).Mag2()/(v0*v0) );
      if (fmax<rnd) continue;

      //// Make sure recoil is possible w/ current energy cut 7/16
      double wimpE = 0.5*wimpM*wimpVLab*wimpVLab/(c*c); // KE of WIMP 
      double mRatio = recoilM/wimpM; // mass ratio: target over projectile
      cosRecoil = sqrt(Emin/(mRatio*wimpE))*(1+mRatio)*0.5;
      if (cosRecoil>-1 && cosRecoil<1) break;

    }
      
      
    while (1) {

      recoilE = gRandom->Rndm() * (Emax-Emin) + Emin;

      //// Calculate cosine of recoil relative to WIMP direction
      double wimpE = 0.5*wimpM*wimpVLab*wimpVLab/(c*c); // KE of WIMP 
      double mRatio = recoilM/wimpM; // mass ratio: target over projectile
      cosRecoil = sqrt(recoilE/(mRatio*wimpE))*(1+mRatio)*0.5;     
      
      
      //// check if kinematics is valid   
      if (cosRecoil>-1 && cosRecoil<1) break;
      
    }
    
    
    //// Construct recoil with random aziumuthal angle & random starting point (mm)
    
    double pi = TMath::Pi();
    
    double phi = 2*pi*gRandom->Rndm(); 
    double sinRecoil = sqrt(1-cosRecoil*cosRecoil);
    
    TVector3 orth = wimpVLab.Orthogonal().Unit() * sinRecoil;
    TVector3 parallel = wimpVLab.Unit() * cosRecoil;
    
    TVector3 recoilV = orth + parallel;
    recoilV.Rotate(phi+180, wimpVLab);
    
    TVector3 startPos(_det->randPos(recoilV,recoilE));

    
    //// Set new parameters for Wimp

    setWimpV(wimpVLab);
    //_wimpV = wimpVLab;
    setTime(t);
    setRecoilE(recoilE);
    setStartPos(startPos);
    setRecoilV(recoilV);

    
    // return phi;
    
  }





  //////////////////////////////////////////////////////////
  ////  genAlpha()
  //////////////////////////////////////////////////////////

  void genAlpha() {

    if(_type == "w") {
      cout << "Can't run genAlpha on wimps" << endl;
      return;
    }
    
    double alphaM = getWimpMass();
    // double Emax = 300; // max Alpha energy
    double c = 3e5;  // speed of light

    double pi = 3.14159265358979;
    

    //// Generate Alpha velocity in lab frame
    
    double vx,vy,vz;
    TVector3 alphaV; // velocity of alpha in lab frame
    double alphaE; // alpha kinetic energy
    
    
    if(_fxn->GetXmin()==0 && _fxn->GetXmax()==0) {
      double vEsc = 600*10; // limit for velocity
      double v0 = 15000;      // average Alpha speed // MAKE SURE SPEED IS RIGHT
    
      double rnd, fmax;
    
      while (1) {
      
	vx = (gRandom->Rndm()-0.5)*vEsc;
	vy = (gRandom->Rndm()-0.5)*vEsc;
	vz = (gRandom->Rndm()-0.5)*vEsc;
	alphaV = TVector3(vx,vy,vz);
         
      
	//// Check if Alpha speed fits Max_Boltz Distr. // MAKE SURE THIS IS RIGHT DISTR.
	rnd = gRandom->Rndm();
	fmax = exp( - (alphaV).Mag2()/(v0*v0) );
	if (fmax>rnd) break;

      }

      alphaE = 0.5*alphaM*alphaV*alphaV/(c*c); 
      alphaV.SetMag(1.0);
    }



    else {
      alphaE = _fxn->GetRandom();
      alphaV.SetX(1);
      alphaV.SetPhi( gRandom->Rndm()*2*pi );
      alphaV.SetTheta( gRandom->Rndm()*pi );
    }


    //// Construct track with random starting point (mm)
    
    TVector3 startPos, normal;

    // The logical check below will be changed when det. model is changed.
    if( (gRandom->Rndm() > 0.5) && (_det->Inside()!=1) && (_det->marginsFit(alphaE))) {
      startPos = _det->randWall(alphaV,alphaE);
      normal = _det->getNormal();
      
      if (normal.Mag()==0) startPos = _det->randPos(alphaV,alphaE);

      // Fold alpha velocity distribution over so it doesn't point outside detector
      else if((normal * alphaV)<0) {
	TVector3 parallel = (normal.Unit() * alphaV) * normal.Unit();
	alphaV = -alphaV; //  - 2*parallel;
      }
    }

    else startPos = _det->randPos(alphaV,alphaE);



    //// Set simulated McAlpha parameters
    setRecoilV(alphaV.Unit());
    setStartPos(startPos);

    setRecoilE(alphaE);
  
  }





  ///////////////////////////////////////////////////////
  ////  cosmicGun()   
  ///////////////////////////////////////////////////////

  TVector3 cosmicGun() {

//     if(_type == "a") {
//       cout << "Can't run cosmicGun on alphas" << endl;
//       return TVector3(0,0,0);
//     }

    // Made all parameters into data members 7/17

    double t = getTime();
    string opt = _opt;    // May want "get methods for these
    double lat = _lat;
    double lon = _lon;

    double pi = 3.14159265358979; // TMath::Pi();
    double conv = pi/180;

  
    TVector3 vS(9,242,7); // Sun's velocity relative to WIMP cloud


    //// Set constants to calculate Earth's velocity

    double
      sE_avg = 29.79, // Earth's average orbital speed
      e = 0.016722,   // ellipticity of Earth's orbit
      L0 = 13 *conv,  // longitude of orbit's minor axis, +- 1 degree uncertainty
    
      Bx = -5.5303 *conv,  //lat & lon of ecliptic's xyz axes - from Lewin-Smith
      Lx = 266.141 *conv,
      By = 59.575 *conv,
      Ly = -13.3485 *conv,
      Bz = 29.812 *conv,
      Lz = 179.3212 *conv;


    double l = 280.460 + 0.9856474*t;
    double g = (357.528 + 0.9856003*t)*conv;

    double L = ( l + 1.915*sin(g) + 0.020*sin(2*g) ) *conv;
    double sE = sE_avg*( 1-e*sin(L-L0) ); // Earth's speed at position L


    double vE_x = sE * cos(Bx) * sin(L-Lx);
    double vE_y = sE * cos(By) * sin(L-Ly);
    double vE_z = sE * cos(Bz) * sin(L-Lz);

    TVector3 vE(vE_x, vE_y, vE_z); // Earth's velocity relative to Sun

    TVector3 vT( vS + vE); // Earth's velocity relative to WIMP cloud



    if (opt=="g") {return vT;}


    // If do not want galactic coord, convert to equatorial
    // More accurate numbers may be available - check to update soon

    vT.Rotate(62.8716639*conv, TVector3(.83867,.54464, 0));
    vT.RotateZ(249.75*conv);
  
    if (opt=="e") {return vT;}


    double side = 1.00273790935;    // Number of sidereal days per day
    double t0 = (lon + 10.0)*conv; // Relative global long. to eq. long. set to 0 at t=0
    vT.RotateZ(-(t*side)*2.0*pi - t0);
    vT.RotateX((lat-90)*conv);
  
    if (opt=="lv") {return vT;}

  
    vT.RotateX(90*conv);

    if (opt=="lh") {return vT;}
    

    else {
      cout << "Invalid option: _opt must be 'g' 'e' or 'l'" << endl;
      return TVector3(0,0,0);
    }

  }
  



  // virtual ~McWimp() {};


  private:

  double _wimpMass;

  TVector3 _wimpV;

  double _time;

  double _recoilE;

  TVector3 _startPos;

  TVector3 _recoilV;

  TH2F _simRecoil;

  string _type;

  DetObj* _det; // detector class added to handle starting position 7/16


  //// Parameters for genWimp and cosmicGun - May want "get" methods for these

  TF1* _fxn;

  string _opt;

  double _Emin;
  
  double _lat;
  double _lon;

  int _ccd; // May be able to delete this


};


#endif
