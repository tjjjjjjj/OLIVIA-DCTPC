#include "TObject.h"
#include "RobinHoodTriangle.hh"

#include "TList.h"

class TVector3;




class RobinHoodTriangleMaker : public TObject {

public:

  RobinHoodTriangleMaker() {}

  ~RobinHoodTriangleMaker() {}

  //TList &getList() const { return &_list; }
  
  
  // Detector geometry - volume:
  void cylinder( TVector3 *center, TVector3 *norm, 
		 double ID=-1.0, double OD=-1.0, double thickness=-1.0, int nTriangles=0, 
		 TString opt="");
  
  // area:
  void sphere( TVector3 *center, TVector3 *norm, 
	       double radius=-1.0, double height=-1.0, int nTriangles=0);
  
  void disk( TVector3 *center, TVector3 *norm, 
	     double Diameter=0, int nTriangles=0,
	     TString opt="");
  
  void rectangle( TVector3 *center, 
		  TVector3 *side1, TVector3 *side2, int nTriangles=0, 
		  TString opt="");



  // RH input files:
  void dumpTriangles(const char* geoFile);

  void addConductorToRhc(const char* rhcFile, double volts, const char * geoFile);


private:

  TList _list;
				  
  ClassDef(RobinHoodTriangleMaker,0)

};
