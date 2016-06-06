#ifndef DMTPC_ABS_LOCATION_HH
#define DMTPC_ABS_LOCATION_HH

#include "TROOT.h"
#include "TVector3.h"

/**This describes the position and orientation of an object by its latitude,
longitude, z position wrt sea level and three euler angles for the
orientation, where in the defined coordinate system the x axis points East
the y axis points North, and the z axis is normal to the surface of the earth
pointing away from the center of the earth \author Asher Kaboth*/

class DmtpcAbsLocation : public TObject {



public:

      /** Constructor
	  \param lat latitude in degrees
	  \param lon longitude in degrees
	  \param z height in meters from sea level
	  \param alpha first euler angle
	  \param beta second euler angle
	  \param gamma third euler angle
      */
   DmtpcAbsLocation(double lat, double lon, double z, double alpha, double beta, double gamma);
      /**
	 Constructor
	 \param positition TVector3 of lat,long,z as defined above
	 \param angles TVector3 of alpha, beta, gamma as defined above
      */
   DmtpcAbsLocation(TVector3* position, TVector3* angles);
    
   virtual ~DmtpcAbsLocation();
    
   virtual const char* GetName() const { return "DmtpcAbsLocation"; }
    
      
      /** 
	  getLatitude
	  \return the set latitude
      */
   double getLatitude() {return _position->X();}
      /** 
	  getLongitude
	  \return the set longitude
      */
   double getLongitude() {return _position->Y();}
      /** getZ
	  \return the set z postion
      */
   double getZ() {return _position->Z();}
      /** getEuler1
	  \return the first euler angle
      */
   double getEuler1() {return _angles->X();}
      /** getEuler2
	  \return the second euler angle
      */
   double getEuler2() {return _angles->Y();}
      /** getEuler3
	  \return the third euler angle
      */
   double getEuler3() {return _angles->Z();}


   void setLatitude(double lat) { _position->SetX(lat);}
   void setLongitude(double lon) { _position->SetY(lon);}
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
    
  ClassDef(DmtpcAbsLocation,1)
};

#endif

