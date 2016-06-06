#ifndef DMTPC_LOCATION_HH
#define DMTPC_LOCATION_HH

#include "TROOT.h"
#include "TVector3.h"

class DmtpcLocation : public TObject {

/**This describes the position and orientation of an object by its x position, y position, z position and three euler angles for the
orientation, where in the defined coordinate system the x axis points East
the y axis points North, and the z axis is normal to the surface of the earth
pointing away from the center of the earth \author Asher Kaboth*/

public:

      /** 
	  Constructor
	  \param x x-position
	  \param y y-position
	  \param z z-position
	  \param alpha first euler angle
	  \param beta second euler angle
	  \param gamma third euler angle
      */
      DmtpcLocation(double x, double y, double z, double alpha, double beta, double gamma);
      /** Constructor 
	  \param position TVector3 of the position as defined above
	  \param angles TVector3 of the euler angles as defined above
      */
   DmtpcLocation(TVector3* position, TVector3* angles);
  DmtpcLocation() {}
    
   virtual ~DmtpcLocation();
    
   virtual const char* GetName() const { return "DmtpcLocation"; }

   double getX() {return _position->X();}
   double getY() {return _position->Y();}
   double getZ() {return _position->Z();}
   double getEuler1() {return _angles->X();}
   double getEuler2() {return _angles->Y();}
   double getEuler3() {return _angles->Z();}

   void setX(double lat) { _position->SetX(lat);}
   void setY(double lon) { _position->SetY(lon);}
   void setZ(double z) { _position->SetZ(z);}
   void setEuler1(double alpha) { _angles->SetX(alpha);}
   void setEuler2(double beta) { _angles->SetY(beta);}
   void setEuler3(double gamma) { _angles->SetZ(gamma);}

   TVector3* getPosition() {return _position;}
   TVector3* getEulerAngles() {return _angles;}
   
   void setPosition(TVector3* position) {_position=position;}
   void setAngles(TVector3* angles) {_angles=angles;}

  
private:
    
   TVector3* _position;
   TVector3* _angles;
    
  ClassDef(DmtpcLocation,1)
};

#endif

