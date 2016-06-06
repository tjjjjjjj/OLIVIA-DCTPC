
#include "TObject.h"

class TVector3;




class RobinHoodTriangle : public TObject {

public:

  RobinHoodTriangle(TVector3 *a, TVector3 *b, TVector3 *c, double charge);

  RobinHoodTriangle(RobinHoodTriangle& t);

  ~RobinHoodTriangle();


  double getCharge() { return _charge; }
  void   setCharge(double charge) { _charge=charge; }

  double getArea() { return _area; }
  
  TVector3* getVertex(int i) { return  (i<0 || i>2) ? 0 : _v[i]; }


private:

  // basic parameters:

  TVector3* _v[3]; // corners of triangle
  double _charge; // charge deposited on triangle
  


  // derived values:
  
  TVector3* _center; //center of triangle
  double _area; // area of triangle


   ClassDef(RobinHoodTriangle, 0)

};
