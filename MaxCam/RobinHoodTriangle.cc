#include "RobinHoodTriangle.hh"
#include "TVector3.h"

ClassImp(RobinHoodTriangle)


RobinHoodTriangle::RobinHoodTriangle(TVector3 *a, TVector3 *b, TVector3 *c, double charge) {

  _v[0]=a;
  _v[1]=b;
  _v[2]=c;

  _charge=charge;

  _center = new TVector3(  (*a)+(*b)+(*c) );
  (*_center) *= 1./3;


  _area = 0.5 * ( (*b)-(*a) ).Cross( (*c)-(*a) ).Mag();
}



RobinHoodTriangle::RobinHoodTriangle( RobinHoodTriangle &t ) : TObject(t) {
  
  _v[0] = t._v[0];
  _v[1] = t._v[1];
  _v[2] = t._v[2];
  
  _charge = t._charge;
  _center = t._center;
  _area   = t._area; 
}




RobinHoodTriangle::~RobinHoodTriangle() {

  delete _v[0];
  delete _v[1];
  delete _v[2];

  delete _center;
}
