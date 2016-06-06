//// Need to find better method for data storage. Ran into bugs with arrays


#ifndef DET_OBJ
#define DET_OBJ

#include "TROOT.h"
#include "TMath.h"
//#include "TList.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom.h"

#include <iostream>
#include <vector>

using namespace std;



class DetObj {


public:

  DetObj(int numCCD=5,
	 int numWires=10,
	 double pX=595.7e-3,
	 double pY=595.7e-3)

  {
    _ccdList = TH2D("","", numCCD,0,numCCD, 5,0,5);
    _wireList = TH1D("","",numWires,0,numWires);

    for(int i=1; i<=numCCD; i++) {_ccdList.SetBinContent(i,5,-1);}

    for(int i=1; i<=numWires; i++) {_wireList.SetBinContent(i,9001e3);}
    // wire position in mm. We won't have 9km detector, so 9001e3 means wire not set

    _pixelX = pX;
    _pixelY = pY;

    _ccd = 0;
    _startInside = 0;

    _height = 250;
    _topMargin = 55;
    _bottomMargin = 920;

    double pi = TMath::Pi();
    _radius = 1000/sqrt(pi);
    _sideMargin = 60;

  }



  //// Methods for CCD configuration
  //// x,y : position of CCD center, in mm in detector frame (center of detector at 0,0)
  //// dX,dY : dimensions of CCD in pixels
  //// imaged : 1 means CCD on, 0 means off, -1 means CCD not set
  //// You must remove the CCD in slot number "n" before adding a CCD in that slot

  void addCCD(double x, double y, double dX, double dY, int imaged, int ccdNum) {

    if(_ccdList.GetBinContent(ccdNum,5) >= 0) {
      cout << "CCD #" << ccdNum << " already set." <<endl;
      return;
    }

    _ccdList.SetBinContent(ccdNum, 1, x);
    _ccdList.SetBinContent(ccdNum, 2, y);
    _ccdList.SetBinContent(ccdNum, 3, dX);
    _ccdList.SetBinContent(ccdNum, 4, dY);
    _ccdList.SetBinContent(ccdNum, 5, imaged);
  }


  double getPosX(int ccdNum) {return _ccdList.GetBinContent(ccdNum,1);}
  double getPosY(int ccdNum) {return _ccdList.GetBinContent(ccdNum,2);}
  double getDimX(int ccdNum) {return _ccdList.GetBinContent(ccdNum,3);}
  double getDimY(int ccdNum) {return _ccdList.GetBinContent(ccdNum,4);}
  int imaged(int ccdNum) {return (int) _ccdList.GetBinContent(ccdNum,5);}


  void removeCCD(int ccdNum) {_ccdList.SetBinContent(ccdNum,5,-1);}

  void ccdOff(int ccdNum) {
    if(imaged(ccdNum)==-1) {
      cout << "CCD not set: can't turn off." << endl;
      return;}
    _ccdList.SetBinContent(ccdNum,5,0);
  }

  void ccdOn(int ccdNum) {
    if(imaged(ccdNum)==-1) {
      cout << "CCD not set: can't turn on." << endl;
      return;}
    _ccdList.SetBinContent(ccdNum,5,1);
  }



  //// Methods for wire configuration
  //// x : position of wire in mm in detector frame (assume all wires parallel to y axis)
  //// if x > 9000e3, wire not set

  void addWire(double x, int wireNum) {
    if (_wireList.GetBinContent(wireNum) < 9000e3) {
      cout << "Wire #" << wireNum << "already set." << endl;
      return;
    }
    _wireList.SetBinContent(wireNum,x);
  }

  void removeWire(int wireNum) {_wireList.SetBinContent(wireNum,9001e3);}



  //// Handle data members

  TH2D cList() {return _ccdList;}
  //void setList(TH2D list) {_ccdList = list;}

  TH1D wList() {return _wireList;}

  void setPixelX(double pX) {_pixelX = pX;}
  double getPixelX() {return _pixelX;}

  void setPixelY(double pY) {_pixelY = pY;}
  double getPixelY() {return _pixelY;}

  void setNormal(TVector3 n) {_normal = n;}
  TVector3 getNormal() {return _normal;}

  void startInside() {_startInside = 1;}
  void startAnywhere() {_startInside = 0;}
  int Inside() {return _startInside;}



  //// Returns position of pixel's center in detector frame
  TVector2 truePos(int ccdNum, double x, double y) {
    double dimX = getDimX(ccdNum);
    double dimY = getDimY(ccdNum);
    double trueX = getPosX(ccdNum) - (dimX*_pixelX/2) + (x*_pixelX) - _pixelX/2;
    double trueY = getPosY(ccdNum) - (dimY*_pixelY/2) + (y*_pixelY) - _pixelY/2;
    return TVector2(trueX,trueY);
  }


  //// Methods for finding an imaged point in xy plane

  int inCCD(int ccdNum, double truX, double truY) {
    double mmX = _pixelX*getDimX(ccdNum);
    double mmY = _pixelY*getDimY(ccdNum);
    double cntrX = getPosX(ccdNum);
    double cntrY = getPosY(ccdNum);
    double minX = cntrX - mmX/2;
    double maxX = cntrX + mmX/2;
    double minY = cntrY - mmY/2;
    double maxY = cntrY + mmY/2;
    if (truX>=minX && truX<=maxX && truY>=minY && truY<=maxY) {return 1;}
    return 0;
  }


  vector<int> inCCDs(double truX, double truY) {
    int off = 0;
    vector<int> list;
    for(int i=1; i<=(_ccdList.GetNbinsX()); i++) {
	if(imaged(i)<0) {continue;}
	if(imaged(i)==0) off++;
	if(inCCD(i,truX,truY) && (imaged(i)==1)) {list.push_back(i);}
    }
    if(off!=0) cout << "Warning: some CCDs are off" <<  endl;
    return list;
  }


  void imageCCDs(int ccd) {_ccd = ccd;}



  TVector3 randPos(TVector3 vec, double E) {
    if (_startInside==1) return randPos1();
    else return randPos2(vec, E);
  }


  TVector3 randPos1() {

    double diam = 2*_radius;
    
    double x,y,z;
    int counter;
    int on;

    while(1) {
      x = (gRandom->Rndm()-0.5)*diam;
      y = (gRandom->Rndm()-0.5)*diam;
      if ( sqrt(x*x + y*y) > _radius ) continue;
      
      if(_ccd!=1) {
	on=-1;
	break;
      }

      counter=0;
      on=0;
      for(int i=1; i<=(_ccdList.GetNbinsX()); i++) {
	if(imaged(i)<=0) {continue;}
	on++;
	if(inCCD(i,x,y)) {counter++;}
      }

      if(on==0) {break;}
      if(counter>0) {break;}
    }

    if(on==0) {
      cout << "CCDs are off: no detection" << endl;
      return TVector3(9e6,9e6,9e6);
    }

    z = (gRandom->Rndm())*_height;

    setNormal(TVector3(0,0,0));
    return TVector3(x,y,z);
  }




  int trackFits(TVector3 start, TVector3 sect, double E) {
    double length = (sect - start).Mag();
    if (length >= (E/10)) return 0;
    else return 1;
  }



  int intersectDet(double x, double y, double z, TVector3 vec, double E) {
    TVector3 tempV;
    TVector3 sect;
    TVector3 start(x,y,z);
    while (1) {
      // check if path intersects with
      // top circular face
      sect = start + vec*((_height-z)/vec.Z());
      if ( (sect.Perp() < _radius) && trackFits(start,sect,E) ) break;
      // bottom circular face
      sect = start + vec*(-z/vec.Z());
      if ( (sect.Perp() < _radius) && trackFits(start,sect,E) ) break;
      // curved face
      tempV = ( vec.Cross(TVector3(0,0,1)) ).Unit();
      tempV = (start * tempV) * tempV;
      double cos = tempV.Perp() / _radius;
      double phi;
      if ( (tempV.Cross(start)).Z() > 0 ) {phi = tempV.Phi()+acos(cos);}
      else {phi = tempV.Phi()-acos(cos);}
      sect.SetPhi(phi);
      sect.SetPerp(_radius);
      double dXY = (start - sect).Perp();
      sect = start + (dXY/vec.Perp())*vec;
      if ( (sect.Z() < _height) && (sect.Z() > 0) && trackFits(start,sect,E) ) break;
      // if all checks fail...
      return 0;
    }
    return 1;
  }



  int intersectCCD(double x, double y, double z, TVector3 vec, double E, int CCD) {
    TVector3 sect;
    TVector3 start(x,y,z);
    double mmX = _pixelX*getDimX(CCD);
    double mmY = _pixelY*getDimY(CCD);
    double cntrX = getPosX(CCD);
    double cntrY = getPosY(CCD);
    double minX = cntrX - mmX/2;
    double maxX = cntrX + mmX/2;
    double minY = cntrY - mmY/2;
    double maxY = cntrY + mmY/2;
    while (1) {
      // check if path intersects with
      // top face
      sect = start + vec*((_height-z)/vec.Z());
      if ( (sect.X()>=minX && sect.X()<=maxX && sect.Y()>=minY && sect.Y()<=maxY) &&
	   trackFits(start,sect,E) ) break;
      // bottom face
      sect = start + vec*((-z)/vec.Z());
      if ( (sect.X()>=minX && sect.X()<=maxX && sect.Y()>=minY && sect.Y()<=maxY) &&
	   trackFits(start,sect,E) ) break;
      // right face
      sect = start + vec*((maxX-x)/vec.X());
      if ( (sect.Z()>=0 && sect.Z()<=_height && sect.Y()>=minY && sect.Y()<=maxY) &&
	   trackFits(start,sect,E) ) break;
      // left face
      sect = start + vec*((minX-x)/vec.X());
      if ( (sect.Z()>=0 && sect.Z()<=_height && sect.Y()>=minY && sect.Y()<=maxY) && 
	   trackFits(start,sect,E) ) break;
      // back face
      sect = start + vec*((maxY-y)/vec.Y());
      if ( (sect.Z()>=0 && sect.Z()<=_height && sect.X()>=minX && sect.X()<=maxX) &&
	   trackFits(start,sect,E) ) break;
      // front face
      sect = start + vec*((minY-y)/vec.Y());
      if ( (sect.Z()>=0 && sect.Z()<=_height && sect.X()>=minX && sect.X()<=maxX) &&
	   trackFits(start,sect,E) ) break;
      // if all checks fail...
      return 0;
    }
    return 1;
  }



  TVector3 randPos2(TVector3 vec, double E) {

    double diam = 2*_radius;
    double margin = 2*_sideMargin;
    
    double x,y,z;
    int counter;
    int on;
    TVector3 end;


    while(1) {
      x = (gRandom->Rndm()-0.5)*(margin+diam);
      y = (gRandom->Rndm()-0.5)*(margin+diam);
      z = (gRandom->Rndm())*(_topMargin+_bottomMargin+_height)-_bottomMargin;
      // z will be changed to reflect double triple mesh config. later

      if ( (sqrt(x*x + y*y) > _radius) || ((z>_height) || (z<0)) ) {
	if (!intersectDet(x,y,z,vec,E)) continue;}
      
      if(_ccd!=1) {
	on=-1;
	break;
      }

      counter=0;
      on=0;
      for(int i=1; i<=(_ccdList.GetNbinsX()); i++) {
	if(imaged(i)<=0) {continue;}
	on++;
	if(inCCD(i,x,y)) {counter++;}
      }

      if(on==0) {break;}

      for(int i=1; i<=(_ccdList.GetNbinsX()); i++) {
	if(imaged(i)<=0) {continue;}
	if(intersectCCD(x,y,z,vec,E,i)) {counter++;}
      }

      if(counter>0) {break;}
    }

    if(on==0) {
      cout << "CCDs are off: no detection" << endl;
      return TVector3(9e6,9e6,9e6);
    }
    
    setNormal(TVector3(0,0,0));
    return TVector3(x,y,z);
  }
  



  TVector3 randWall(TVector3 vec, double E) {
    // If you're only taking events that start within the working region
    // alphas/recoils that start from the detector walls are excluded, so
    // randWall1() is currently obsolete.  If Hayk/Denis's rings/wires or the
    // cathodes/meshes end up being a source of background, randWall1()
    // will be activated again.
    // if (_startInside==1) return randWall1();
    // else return randWall2(vec,E);
    return randWall2(vec,E);
  }



  int marginsFit(double E) {
    if (_topMargin < (E/11)) return 1;
    if (_bottomMargin < (E/11)) return 1;
    if (_sideMargin < (E/11)) return 1;
    return 0;
  }




//   TVector3 randWall1() {
    
//     double pi = TMath::Pi();
//     double diam = 2*_radius;
//     double areaCir = 2*pi*_radius*_radius;
//     double areaRec = 2*pi*_radius*_height;

//     double x,y,z;
//     int counter, on;

//     while(1) {

//       //////// If on the curved wall ...
//       if(gRandom->Rndm() > areaCir/(areaCir+areaRec)) {
// 	z = (gRandom->Rndm())*_height;
// 	double theta = gRandom->Rndm()*2*pi;
// 	x = _radius*cos(theta);
// 	y = _radius*sin(theta);
    
// 	setNormal(TVector3(-x,-y,0));
//       }

//       //////// If on the circular top/bottom ... 
//       else {
// 	int zz;
// 	if(gRandom->Rndm() > 0.5) {
// 	  z = _height;
// 	  zz = -1;
// 	}
// 	else {
// 	  z = 0;
// 	  zz = 1;
// 	}

// 	while(1) {
// 	  x = (gRandom->Rndm()-0.5)*diam;
// 	  y = (gRandom->Rndm()-0.5)*diam;
// 	  if( sqrt(x*x + y*y) < _radius ) break;
// 	}
      
// 	setNormal(TVector3(0,0,zz));
//       }

//       //////// Make sure within the range of ccds
//       if(_ccd!=1) {
// 	on=-1;
// 	break;
//       }

//       counter=0;
//       on=0;
//       for(int i=1; i<=(_ccdList.GetNbinsX()); i++) {
// 	if(imaged(i)<=0) {continue;}
// 	on++;
// 	if(inCCD(i,x,y)) {counter++;}
//       }

//       if(on==0) {break;}
//       if(counter>0) {break;}
//     }


//     if(on==0) {
//       cout << "CCDs are off: no detection" << endl;
//       return TVector3(9e6,9e6,9e6);
//     }

//     return TVector3(x,y,z);

//   }



  TVector3 randWall2(TVector3 vec, double E) {
    
    double pi = TMath::Pi();
    double detRad = _radius+_sideMargin;
    double diam = 2*detRad;
    double detHgt = _height + _topMargin + _bottomMargin;
    double areaCir = 2*pi*detRad*detRad;
    double areaRec = 2*pi*detRad*detHgt;

    double x,y,z;
    int counter, on, misses;

    misses=0;
    while(1) {

      //////// If on the curved wall ...
      if(gRandom->Rndm() > areaCir/(areaCir+areaRec)) {
 	z = (gRandom->Rndm())*detHgt - _bottomMargin;
	double theta = gRandom->Rndm()*2*pi;
	x = detRad*cos(theta);
	y = detRad*sin(theta);
   
	setNormal(TVector3(-x,-y,0));
      }

      //////// If on the circular top/bottom ... 
      else {
	int zz;
	if(gRandom->Rndm() > 0.5) {
	  z = detHgt - _bottomMargin;
	  zz = -1;
	}
	else {
	  z = -_bottomMargin;
	  zz = 1;
	}

	while(1) {
	  x = (gRandom->Rndm()-0.5)*diam;
	  y = (gRandom->Rndm()-0.5)*diam;
	  if( sqrt(x*x + y*y) < detRad ) break;
	}
      
	setNormal(TVector3(0,0,zz));
      }
 
      misses++;
      if (misses>100) {
	setNormal(TVector3(0,0,0));
	break;
      }

      //// Make sure within working region of detector
      if (1) { // ( (sqrt(x*x + y*y) > _radius) || ((z>_height) || (z<0)) ) {
	if (!intersectDet(x,y,z,vec,E)) continue;}

      //////// Make sure within the range of ccds
      if(_ccd!=1) {
	on=-1;
	break;
      }

     counter=0;
      on=0;
      for(int i=1; i<=(_ccdList.GetNbinsX()); i++) {
	if(imaged(i)<=0) {continue;}
	on++;
	if(inCCD(i,x,y)) {counter++;}
      }

      if(on==0) {break;}

      for(int i=1; i<=(_ccdList.GetNbinsX()); i++) {
	if(imaged(i)<=0) {continue;}
	if(intersectCCD(x,y,z,vec,E,i)) {counter++;}
      }

      if(counter>0) {break;}
    }


    if(on==0) {
      cout << "CCDs are off: no detection" << endl;
      return TVector3(9e6,9e6,9e6);
    }

    return TVector3(x,y,z);

  }



private:
 
  double _height; // height of detector working region (distance between cathode & ground)

  double _radius; // radius of detector working regoin

  double _topMargin; // distance from cathode to top wall of detector

  double _bottomMargin; // distance from ground to bottom wall
                        // This will be corrected to be the distance from the bottom
                        // cathode to the bottom wall after initial simulations.

  double _sideMargin; // difference between detector wall radius and working regoin radius

  double _pixelX; // mm per pixel

  double _pixelY;

  TH2D _ccdList;

  TH1D _wireList;

  int _ccd; // Consider CCDs in generating starting positions?

  int _wire; // Consider wires in generating starting positions?
  
  int _startInside; // Only generate events that start inside FoV?

  TVector3 _normal;

};


#endif

