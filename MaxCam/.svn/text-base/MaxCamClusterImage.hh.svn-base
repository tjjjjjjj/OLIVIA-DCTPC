#ifndef MAXCAM_CLUSTERIMAGE_HH
#define MAXCAM_CLUSTERIMAGE_HH

#include "TMath.h"
#include "TH2.h"
#include "TDatime.h"
#include "DmtpcStitchInfo.hh"
#include "TTimeStamp.h"
#include "DmtpcGainMap.hh"
#include <vector>
#include <list>
#include "TGraph.h"

using std::vector;

/** Holds information from track finding algorithm and basic reconstructed quanities for the found tracks. Tracks are found through MaxCamImageTools::findClustersCI. \author Asher Kaboth */


class TCutG;

class MaxCamClusterImage : public TObject{
      
public:
   
   MaxCamClusterImage();
   /**
      Constructor
      \param image ccd image
      \param time TDatime for the time image was taken
   */
   MaxCamClusterImage(TH2* image, TDatime* time);
   /**
      Constructor
      \param image ccd image
      \param time TTimeStamp for the time image was taken
   */
   MaxCamClusterImage(TH2* image, TTimeStamp* time);
   //destructor
   ~MaxCamClusterImage();
   
      /**Add new cluster 
	 \param newpx vector of bin numbers of new cluster
      */
   void addCluster(vector < int > newpx, vector<int> redpx = vector<int>(),  vector<vector<char> > * frac = 0 );
      /**get number of clusters
	 \return number of clusters
      */
   int getNCluster() const {return _pixels.size();}

      /**returns the amount of light
	 \param i ith cluster
      */
   double getIntegral(int i) const ;
      /**returns the amount of light, weighted by a gain map
	 \param i ith cluster
      */
  double getIntegralWithGainMap(int i, const DmtpcGainMap* gainMap, double mingain = -1) const ;
      /**Finds the length of a cluster from the reduced cluster, by finding the two pixels furthest apart.
	 \param i ith cluster
	 \param x1 xcoordinate of one end
	 \param y1 ycoordinate of one end
	 \param x2 xcoordinate of other end
	 \param y2 ycoordinate of other end
	 \return length between (x1,y1) and (x2,y2)
	 
      */
   double getLength(int i, double &x1, double &y1, double &x2, double &y2) const;
      /**Finds the length of a cluster from the full cluster, by finding the two pixels furthest apart.
	 \param i ith cluster
	 \param x1 xcoordinate of one end
	 \param y1 ycoordinate of one end
	 \param x2 xcoordinate of other end
	 \param y2 ycoordinate of other end
	 \return length between (x1,y1) and (x2,y2)
	 
      */
   double getDiffusedLength(int i, double &x1, double &y1, double &x2, double &y2) const;
      /** returns the angle in the x-y plane by calculating the angle of the line returned by getLength. 
	  \param i ith cluster
	  \return phi
      */
   double getPhi(int i) const;

      /**returns the x and y position of the track as a an average of the x and y positions, weighted by the amount of light in each pixel in the cluster
	 \param i ith cluster
	 \param x xposition output
	 \param y yposition output

      */
   void getXY(int i, double &x, double &y) const;
      /**finds out if bin is in the cluster
	 \param i ith cluster
	 \param bin test bin
	 \return true if in, false if out
      */
   bool isInCluster(int i, int bin, bool red = false) const;
      /**returns the list of bins
	 \param i ith cluster
	 \return vector of ints containing the bins
      */

  /** is the cluster contained in a region?
      \param i ith cluster
      \param region a TCutG region
      \return true if all cluster pixels fall inside the TCutG
      .       false otherwise
   */
  bool isContainedInRegion(int itrk, TCutG* region) const;


  bool hitsVeto(int itrk, double rin, double rout) const; 
  bool crossesCameras(int itrk, const DmtpcStitchInfo * stitch) const; 
  bool hitsInactive(int itrk, double rout) const; 


  /** Apply a threshold to the reduced cluster pixels.
   * 
   *  @param threshold the threshold, in ADU. 
   *
   * **/
  void applyRedThreshold(double threshold);

  /** is any part of the cluster inside the region?
      \param i -- ith cluster
      \param region -- a TCutG region
      \return true if any cluster pixel falls inside the TCutG
      .       false otherwise
   */
  bool intersectsRegion(int itrk, TCutG* region) const;

  vector< vector< int> > getClusters() const {return _pixels;}

   vector<int> getCluster(int i) const {return _pixels[i];}
      /** returns the reduced cluster 
	  \param i ith cluster
	  \return vector of ints containing the bins
      */
   vector <int> getClusterRed(int i) const {return _pixelsred[i];}

      /**returns theta (=0 for the moment)
	 \param i ith cluster
      */
    double getTheta(int /*i*/) { return 0; }
      /** returns the time
       */
   double getTime() const {return _time;}
   
      /** returns the image the clusters were found in 
       */
   const TH2* getImage() const {return _image;}


   /** Sets the image to NULL */ 
   void forgetImage() { if (_image) delete _image; _image = 0; } 

   /** Sets the image the the given value without changing binning */ 
   void setImage(TH2 * img)  { _image =  img; } 

      /** changes the image. New image must have same or finer binning.
	  \param newimage the new image
      */
   void changeImage(TH2* newimage);
      /** changes the image and applys a threshold to make the reduced cluster
	  \param newimage the new image
	  \param threshold the threshold pixels must be above to make it into the reduced cluster
      */
   void changeImageWithThreshold(TH2* newimage, double threshold);
      /** sets the time \param time a TDatime with the new time */
   void setTime(TDatime* time);
      /** sets the time \param time a TTimeStamp with the new time */
   void setTime(TTimeStamp* time);

      /** calulates the major and minor axes of the best fit ellipse
	  \param i track number
	  \param a destination for major axis
	  \param b destination for minor axis
      */
      void getEllipseAxes(int i, double& a, double& b);
      /** calculates the energy density
	  \param i ith cluster
      */

   double getEnergyDensity(int i) {return getIntegral(i)/double(_pixels[i].size());}

   void morphologicalOperation(int i, int nerode , int ndilate, bool red  = false); 
   
      /**
	 finds a vector angle based on elliptical analysis and getAsymmetry
	 \param i ith cluster
	 \return phi
       */
   double getPhi2(int i, int nerode = 0, int ndilate = 0) const;
      /** finds an axial  angle based on the angle with maximum sigma
       \param i ith cluster
       \return phi
      */
   double getPhi3(int i ) const;
      /** finds an axial angle by fitting a line to the cluster
	  \param i ith cluster
	  \return phi
      */
   double getPhi4(int i) const;


    /** 
	  calculates the length along an angle
	  \param i ith cluster
	  \param theta angle to calculate along
	  \param pxpermm conversion factor to real distance 
     */
   double getLength2(int i, double theta, double pixpermm = 1);
      /** calculates third moment along an angle 
	  \param i ith cluster
	  \param theta angle to calculate along
	  \return skewness
      */
   double getSkewness(int i, double theta);
      /** calculates two bin "skewness" to determine direction of track
	  \param i ith cluster
	  \param phi angle to calculate along
	  \return 
       */
   double getAsymmetry(int i, double phi) const;
      /** calcualtes whether a track crosses the edge of the image
	  \param i ith cluster
      */
   bool hitsEdge(int i);
 
      /** enum to categorize camera orientation (10L only)
       */
   typedef enum
   {
      BOTTOM,
      TOP
   }CAMERA_ORIENTATION; 

   /* Returns the cygnus angle of the ith cluster.
        yangN: Angle of the y axis IN DEGREES with respect to North
        cam_orientation: (see enum above)
        lat,lon: latitude and longitude of the detector. 
          Default values for Cambridge. 
        phi: phi coordinate of the cluster or DBL_MAX to use the 
           value from getPhi2()
        theta: theta coordinate of the cluster or DBL_MAX to use 
           the value from getTheta() 
   */
   double getCygnusAngle(int i, double yangN, 
                         CAMERA_ORIENTATION cam_orientation, 
                         double lat=42.373611,
                         double lng = -71.110556,
                         double phi=DBL_MAX,
                         double theta=DBL_MAX );
   
   void getRADec(double phi, double theta, TDatime * time, double lat, double lon, double nang, 
                   CAMERA_ORIENTATION cam_orientation, double & ra, double & dec, 
                   double & l, double & b); 
                  

      /**Returns the mean value of the light of the cluster 
	 \param i ith cluster*/
   double getMean(int i) { return getEnergyDensity(i); }

      /**Returns the rms of the ith cluster 
       \param i ith cluster
       \param mean mean of the cluster */
   double getRMS(int i,double mean); 

      /**Returns the maximum value of the ith cluster. If maxBin is not NULL, the int at that value will be populated with the bin number of the maximum value.
	 \param i ith cluster
	 \param maxBin loction for bin number of maximum value
      */
   double getMax(int i, int * maxBin = NULL); 

      /**Returns the neighbors around the maximum pixel above threshold of the ith cluster. The threshold is to be given in units of the rms (calculated with meanRMSNoOutliers()  If the third argument is positive, then that will be used as the maxBin. Otherwise this method will calculate it using getMax (you may want to populate it for performance reasons)
	 \param i ith cluster
	 \param threshold for neighbor finding
    */
   int getNumNeighbors(int i, double thresholdInRMSUnits, int maxBin =-1);
   
   /** Draws the ith cluster in the MaxCamImage with a borderpx of 20. Makes copy of image so there should be no worries...
    \param i ith cluster
    \param borderpx number of pixels to draw around the cluster
   */
   void drawCluster(int i, int borderpx=20);

   /** Returns a list of TGraphs of a cluster's boundary. Each TGraph is a simple line segment. You
    *  should delete the TGraphs when you're done with them. 
    *  \param i ith cluster
    *  \param color the color of the boundary
    *  \param linewidth the linewidth of the boundary
    *  \param draw if true, the boundary will get drawn 
    * **/
      std::list<TGraph*> getClusterBoundary(int i, int color = 15, int linewidth = 2, bool draw = false, bool red = false) const; 
    

   /** Returns the perimeter of the cluster
    *  \param i ith cluster
    * **/
      int getClusterPerimeter(int i, bool red = false) const; 
    

   /** Draws a region of the image from xmin to xmax and ymin to ymax
    */
   void drawRegion(int xmin, int xmax, int ymin, int ymax); 
   /**
      finds the extremes of the cluster, plus a border
      \param i ith cluster
      \param xmin return variable for the lower x edge of the range
      \param xmax return variable for the upper x edge of the range
      \param ymin return variable for the lower y edge of the range
      \param ymax return variable for the upper y edge of the range
      \param borderpx size of the border for the region

    */
   void clusterBounds(int i, int * xmin, int * xmax, int * ymin, int * ymax, int borderpx=20) const; 
   
   /** finds the x and y bins or the bin centers of a pixel from its master bin number
       \param bin master bin
       \param x return variable for xcoordinate
       \param y return variable for ycoordinate
       \param undobinning false to return the x and y bin numbers; true to return the x and y bin centers
   */
   void getXYFromBinNo(int bin, int * x, int * y, bool undobinning = false) const; 

   /** finds the min and the max of the expression x*cos(phi)+y*sin(phi)
       \param i ith cluster
       \param phi angle to calculate from
       \param minVal return value for the minimum 
       \param maxVal return value for the maximum
 
    */
   void getMinMaxPosition(int i, double phi, double &minVal, double &maxVal) const;

   /**
   Finds any integer moment of a track along a given axis.
   \param i ith cluster
   \param n nth central moment
   \param phi Angle to calculate moments along
   \param nbins binning set from the binType variable. 0 for unbinned moment.
   \param binType
   \return nthe moment of ith vector.
   */
   double getMoment(int i,int n, double phi, int nbins = 4, TString binType = "totalBins");

   /**
   Finds any integer moment of a track along a given axis.  Uses arrays instead of histograms for better memory management.
   \param i ith cluster
   \param n nth central moment
   \param phi Angle to calculate moments along
   \param nbins binning set from binType variable: either pixels per bin or total number of bins.
   \param binType pixelPerBin or totalBins
   \return nth moment of ith vector
   */
   double getMoment2(int i, int n, double phi, int nbins = 4, char* binType = "pixelPerBin");

   /**
   Finds the values for a given set of moments we wish to calculate without doing any binning.
   \param i ith cluster
   \param phi Angle to calculate the moments along
   \param n the number of moments to calculate
   \param moments an array of moments to calculate
   \param useClustRed true: use reduced cluster, false: use full cluster
   \return the values of these moments in values
   */
   void getMomentsUnbinned(int i,double phi, int n, int* moments, double* values,bool useClustRed = false);

   /**
   Finds the values for a given set of moments we wish to calculate after projecting the cluster along an angle \f$\phi\f$ and binning the values
   \param i ith cluster
   \param phi Angle to calculate the moments along
   \param n the number of moments to calculate
   \param moments an array of moments to calculate
   \param binning determines the bin sizes for the total bins or
   \return the values of these moments in values
   */
   void getMoments(int i,double phi,int n, int* moments, double* values,int binning = 8, char* opt = "");

   /**
   Calculates the value of a Rayleigh vector of a given cluster.  The Rayleigh angle is defined as

\f$\vec{R} = \left(\sum\limits_i w_i\right)^{-1}\sum\limits_iw_i\frac{\vec{x}_i-\vec{x}_0}{|\vec{x}_i-\vec{x}_0|} \f$

  where \f$i\f$ represents each bin.  The current weighting scheme is to remove all negative valued bins and take the square of the value in the remaining bins to be the weight.
Thus, the Rayleigh angle is a weighted mean of the unit vector of each bin in the cluster and some center position \f$\vec{x}_0 \f$.

For a nuclear recoil with head-tail, the Rayleigh angle should point in the direction of the "head" part of the track, so \f$-\vec{R}\f$ gives us a rough estimate of the angle with head-tail.  The head-tail discrimination with this algorithm works significantly better than the moments algorithm, but the axial angle resolution is worse.

  There are several options that can be set.
 
   Center position option:

   'a' or '': Find the center of mass, weighting each bin the same

   'b': Find the center of the smallest box surrounding the cluster

   'c': Find the center of mass, weighting each bin by its value in the histogram

   Cluster option:

   '': Use full cluster

   'r': Use reduced cluster
 
   Other options:

   'u': Use unweighted bin values.  This will find asymmetry in the cluster shape

   'v': Use verbose mode for more output.  Good for debugging.

   \param i ith cluster
   \param opt the options to choose
   \return x the vector component along the x axis
   \return y the vector component along the y axis
   */
   void getRayleigh(int i, double& x, double& y, char* opt="c");


   /**
   Project a cluster along a given axis.  For binType, set to totalBins if binning gives the total number of bins in the projection.  Set to pixelPerBin if binning represents the projection bin width in pixels.
   \param i ith cluster
   \param phi projection angle/axis.
   \param binning the value for the binning option
   \param the method to determine bin sizes
   */
   TH1* projectCluster(int i, double phi, int binning, TString binType="totalBins") const;

   /** Get the portion of an image around the cluster as a new cluster **/ 
   TH2* getClusterHist(int i, int padding=1, bool red = false, const DmtpcGainMap * map = 0, double min_gain = 0, const char * name = 0, bool setzero = true) const; 

   /** Project a cluster along a given axis. phi=0 is x. Use the given interpolation method to rotate the pixels prior to projecting them. 
    */ 
   TH1* projectClusterInterpolate(int i, double phi, const char * interpolation = "bilinear", const DmtpcGainMap * map = 0, double min_gain = 0, bool reset_xaxis = true, const char * name = 0, double * startpos = 0) const; 


   /** Get the projection on the y axis of the pixels to the left of cluster, potentially excluding clusters **/  
   TH1 * getLeftProj(int i, int border = 8, bool ignoreClusters = true, double outlier_pct = 0.01) const; 

   /** merge two tracks
       \param i ith track
       \param j jth track
    */
   void mergeTracks(unsigned int i, unsigned int j);
   /** alternate energy finding method; finds energy in a box around the cluster
       \param i ith cluster
       \param nbin number of bins around the cluster to also integrate over
    */
   double getIntegral2(int i, int nbin);

   void meanRmsNoClustersNoHotSingle(double * mean, double * rms, double outlierFactor = 3) const; 
  
   void changeHistType(char type='S'); 
   void roundValues(bool onlyOutsideClusters = true, double roundTo=1.); 

//   vector<vector<char> > getCamFrac(int cluster) const { return _camfrac[cluster]; }

private:
   
   TH2* _image;
   vector< vector <int> > _pixels;
//   vector< vector< vector<char> > > _camfrac; //store as char out of fraction of 255 
   vector< vector <int> > _pixelsred;
   int _nbinsx;
   int _nbinsy;
   int _nbins;
   double _time; //in julian days past the J2000 epoch
   double image_rms; //Will be calculated once per image if needed.  
   
  
   
   ClassDef(MaxCamClusterImage,4)
   
};

#endif
