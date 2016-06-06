#ifndef MAXCAM_CLUSTER_HH
#define MAXCAM_CLUSTER_HH

#include "TMath.h"
#include "TH2F.h"
#include "TDatime.h"
#include <vector>

using std::vector;

/** DEPRECATED! \deprecated Predecessor to MaxCamClusterImage. \author Asher Kaboth */

class MaxCamCluster : public TObject{
      
public:
   
   //constructors
   MaxCamCluster();
   MaxCamCluster(vector<int> pixels, TH2* image);
   //destructor
   ~MaxCamCluster();
   
   void setTime(TDatime* time);
   
   //returns the amount of light
   double getIntegral();
   //finds the length of a cluster; returns endpoints (x1,y1) and (x2,y2)
   double getLength(double &x1, double &y1, double &x2, double &y2);
   //returns the angle in the x-y plane
   double getPhi();
   //returns the x and y position of the vertex
   void getXY(double &x, double &y);
   //finds out if bin is in the cluster
   bool isInCluster(int bin);
   //returns the list of bins
   vector<int> getCluster() const {return _pixels;}
   //returns theta (=0 for the moment)
   double getTheta() const {return 0;}
   
   TH2* getImage() {return _image;}
   
   void changeImage(TH2* newimage);
   
   double getEnergyDensity() {return getIntegral()/double(_pixels.size());}
   
   double getPhi2();
   double getPhi3();
   double getPhi4();
   double getLength2(double theta, double pixpermm = 1);
   double getSkewness(double theta);
   bool hitsEdge();
   
   double getCygnusAngle(double yangN, double theta);
   
   
private:
   
   TH2* _image;
   vector< int > _pixels;
   int _nbinsx;
   int _nbinsy;
   int _nbins;
   double _time; //in julian days past the J2000 epoch
   
   ClassDef(MaxCamCluster,2)
   
};

#endif
