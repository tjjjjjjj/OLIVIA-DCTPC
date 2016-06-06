#ifndef MAXCAM_IMAGE_TOOLS_HH
#define MAXCAM_IMAGE_TOOLS_HH

#include "TH2F.h"
#include "TH1F.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TNtuple.h"
#include "TString.h"
#include "TArray.h"
#include "MaxCamCluster.hh"
#include "MaxCamClusterImage.hh"
#include "DmtpcGainMap.hh"
#include "DmtpcStitchInfo.hh"
#include "DmtpcDataset.hh" //for TH2Fkappasigmaclip()
#include "DmtpcRootTools.hh"
#include "MaxCamImageTools_Template.hh"
#include <cmath>

class TVector2;

#include <vector>
#include <set>
using std::vector;
class TGraph;
class TCanvas;
class TPoints;
class TCutG;

/** Macro to create a new TH2 with the same binning as another */
#define NEW_HIST2D_WITH_SAME_SIZE(IMG, TYPE, NAME) (new (TYPE) ((NAME),(NAME), (IMG)->GetNbinsX(), (IMG)->GetXaxis()->GetXmin(), (IMG)->GetXaxis()->GetXmax(), (IMG)->GetNbinsY(), (IMG)->GetYaxis()->GetXmin(), (IMG)->GetYaxis()->GetXmax())) 


/** Namespace of functions for performing functions on images. 
 * Typically take a TH2, plus some parameters to perform functions
 *
 *  The implementations of these functions are divided among the following source files: 
 *  
 *
 *  **MaxCamImageTools_ClusterFinding.cc: Cluster finding algorithms and exclusive helper functions
 *
 *  **MaxCamImageTools_ImageCalculations.cc: Calculations on an image or regions of an image but NOT calculations on groups of pixels
 *
 *  **MaxCamImageTools_ImageProcessing.cc: Filters, convolutions, and other things that modify an image (including stuff like killing pixel lists, pedestal subtraction, gradient, blurring, median filter, edge detection, thresholds)
 *
 *  **MaxCamImageTools_ImageTransforms.cc: Things that transform an image. Distinction with Processing is not too great. This includes stuff like rotations, hough transforms, polar coordinates, montage and padding.
 *
 *  **MaxCamImageTools_Interpolation.cc: Bicubic interpolation implementation
 *
 *  **MaxCamImageTools_IO.cc: FITS input and output and display stuff
 *
 *  **MaxCamImageTools_PixelProcessing.cc: Calculations and operations on groups of pixels 
 *
 *
 * */

namespace MaxCamImageTools {

  typedef enum
  {
    BILATERAL_GAUSSIAN,
    BILATERAL_CAUCHY, 
    BILATERAL_BOX, 
    BILATERAL_TRIANGLE, 
    BILATERAL_TUKEY,  
  } BILATERAL_VALUE_FN;

 
  /**
     Set the minimum and maximum of a TH2F for plotting purposes
     
     Implementation in MaxCamImageTools_IO.cc

   */
  void setFriendlyRange(TH2* img, float nsigma=3.0);




   /** Returns a histogram of the distribution of bin yields; e.g. 1 bin with value 0, 3 bins with value 1, etc.
       \param h input histogram
       \param min minumum bin value to histogram
       \param max maximum bin value to histogram
       \return a TH1F containting the yield histogram
      
       Implementation in MaxCamImageTools_ImageCalculations.cc
   */
  TH1F* createYieldHisto(TH1* h, float min=-1e10, float max=1e10);
  

   /** Applies an overall threshold to an image; sets all bins not above the threshold to zero
       \param image input image
       \param threshold threshold      

       Implementation in MaxCamImageTools_ImageProcessing.cc
    */
  void applyThreshold(TH2* image, float threshold);
   /** Applies an overall ceiling to an image; sets all bins not below ceiling to zero
       \param image input image
       \param ceiling ceiling      

       Implementation in MaxCamImageTools_ImageProcessing.cc
    */
  void applyCeiling(TH2* image, float ceiling);


   /** Determines the number of neighbors of a particular pixel which are above a given threshold
       \param image input image (1D)
       \param xbin test bin
       \param threshold threshold which neighbors must be above to be counted
       \param minNeighbors number of neighbors which the bin must have to return true
       \returns true if the number of bins above threshold is greater than or equal to minNeighbors; false otherwise

       Implementation in MaxCamImageTools_ImageCalculations.cc
   */
  bool hasNeighbor(TH1F* image, int xbin, float threshold, int minNeighbors=2);


  /** Determines the number of neighbors of a particular pixel which are above a given threshold
       \param image input image 
       \param xbin test bin
       \param threshold threshold which neighbors must be above to be counted
       \param minNeighbors number of neighbors which the bin must have to return true
       \returns true if the number of bins above threshold is greater than or equal to minNeighbors; false otherwise

       Implementation in MaxCamImageTools_ImageCalculations.cc
   */
  bool hasNeighbor(TH2* image, int xbin, int ybin, float threshold, int minNeighbors=2);


   /** Determines the neighbors of a particular pixel which are above a given threshold
       \param image input image 
       \param ibin test bin
       \param threshold threshold which neighbors must be above to be counted
       \returns a vector of the bin numbers that are above threshold

       Implementation in MaxCamImageTools_PixelProcessing.cc
   */
  vector<int> findNeighbors(TH2 *image, float threshold, int ibin);
  
   /** Set to average all pixels that are above threshold
       \param image input image
       \param threshold threshold above which all pixels will be set to the image mean
       \return number of pixels killed

       Implementation in MaxCamImageTools_ImageProcessing.cc
    */
   int killPixels(TH2* image, float threshold,vector<int>* killedPixels=NULL);


  /**Set to zero all pixels that are below the threshold, or have no neighbors above the threshold.
     \param image input image
     \param threshold threshold
     \param minNeighbors minimum number of neighbors a pixel must have to survive
     \return number of pixels killed

       Implementation in MaxCamImageTools_ImageProcessing.cc
   */
  int killLonePixels(TH2* image, float threshold, int minNeighbors=2);


  /**Set to the image mean all pixels that are below the threshold, or have no neighbors above the threshold.
     \param image input image
     \param threshold threshold
     \return number of pixels killed

       Implementation in MaxCamImageTools_ImageProcessing.cc
   */
   int killLonePixels2(TH2* image, float threshold,vector<int>* killedPixels=NULL);


   /** takes an image, usually an averaged image, and finds arrays of pixels that are hot or cold
       \param image input image
       \param low threshold below which pixels are "cold"
       \param high threshold above which pixels are "hot"
       \param cold return variable for array of cold pixels
       \param hot return variable for array of hot pixels
       \param ncold return varible for number of cold pixels
       \param nhot return variable for number of hot pixels

       Implementation in MaxCamImageTools_PixelProcessing.cc
   */
  void findHotCold(TH2* image, float low, float high, int* cold, int* hot, int &ncold, int &nhot);
  

   /** set a vector of pixels to zero
       \param image input image
       \param pixellist list of pixels to kill
       \return number of pixels killed

       Implementation in MaxCamImageTools_ImageProcessing.cc
   */
  int killPixelList(TH2* image, vector<int> pixellist);
  
   /* kills pixels outside of a designated region of interest. The region is rectangular.
      \param image input image 
      \param roi vector of bins containing pixels to be saved
      \param d number of bins to preserve outside of the roi
      \param opt choose x or y to kill the bins in

       Implementation in MaxCamImageTools_ImageProcessing.cc
   */
  int killUnusedPixels(TH2* image, vector<int> &roi, int d, TString opt);


   /** Sets a row in a 2d image to zero
       \param image input image
       \param irow the row to kill

       Implementation in MaxCamImageTools_ImageProcessing.cc
   */
  int killRow(TH2* image, int irow);


   /** Sets a column in a 2d image to zero
       \param image input image
       \param icolumn the column to kill

       Implementation in MaxCamImageTools_ImageProcessing.cc
   */
  int killColumn(TH2* image, int icolumn);


   /** subtracts a constant from the value of every bin in the image
       \param image input image
       \param pedestal constant pedestal to be subtracted

       Implementation in MaxCamImageTools_ImageProcessing.cc
   */
  void subtractPedestal(TH2* image,float pedestal);
  

   /**
       Implementation in MaxCamImageTools_ClusterFinding.cc
    */
  vector<int> killSecondaryClusters(TH2 *image, float threshold);
  
  /** Implementation in MaxCamImageTools_ImageCalculations.cc **/ 
  double calcSkewness(TH1* strip, int minx, int maxx, double bkg=0);

  /** Implementation in MaxCamImageTools_ImageCalculations.cc **/ 
  double calcSkewness2D(TH2* image, int minbin, int maxbin, double bkg=0);
  
  /** Implementation in MaxCamImageTools_ImageCalculations.cc **/ 
  double calcPixelDistance(TH2* image, int i, int j);

  /** Implementation in MaxCamImageTools_ImageCalculations.cc **/ 
  double calcPixelCorrelation(TH2* image, int i, int j);
  
  /** Implementation in MaxCamImageTools_ImageCalculations.cc **/ 
  double calcIntensity2D(TH2 *image, int imax, int jmax, int sidebandBins, float &bg);
  
  /** Implementation in MaxCamImageTools_ImageCalculations.cc **/ 
  void countSegmentLength(TH1 *project, float threshold, int pixelsBelowThreshold, int &iR, int &iL);

  /** Implementation in MaxCamImageTools_ImageCalculations.cc **/ 
  void countSegmentLength2D(TH2 *project, float threshold, int &iR, int &iL);
  
  /** Implementation in MaxCamImageTools_ClusterFinding.cc **/ 
  int findClusters(TH2* image, MaxCamCluster** array, 
		   double minsig, double maxsig, int minsize, double mindist);

  /** Implementation in MaxCamImageTools_ClusterFinding.cc **/ 
  int findClustersCI(TH2* image, MaxCamClusterImage* clust, 
		   double minsig, double maxsig, int minsize, double mindist);

  /** Implementation in MaxCamImageTools_ClusterFinding.cc **/ 
  int findClustersHysteresisGM(const TH2 *image, MaxCamClusterImage * clust, 
                             double high_thresh, double low_thresh, 
                             unsigned int min_neighbors, unsigned int minsize, 
                             double mindist, 
                             double minrxy_global, double minrxy_cluster, double maxjoin_residual,
                             double least_squares_weight, const DmtpcGainMap * map, double spacer_width = 2,
                             bool allow_edge_seed = false, bool allow_outliers = true); 

 int findClustersGMRing( const TH2 * image, MaxCamClusterImage * clust, const DmtpcStitchInfo * sinfo, 
                         const double * image_rms, const double * image_means, double space_sigma, 
                         double rms_sigma, double core_thresh_high, double core_thresh_low,  
                         double ring_thresh, double ring_nsigma, unsigned int min_size, BILATERAL_VALUE_FN fn, 
                         unsigned ncleanup, 
                         double minJoinDistance, double minrxy_global, double minrxy_cluster, 
                         double maxjoin_residual, double ls_weight, const DmtpcGainMap * map = 0,
                         double spacer_width = 2, const char * debug = 0, int fast_bilateral_scale_exp = 7); 


  /** Implementation in MaxCamImageTools_ClusterFinding.cc **/ 
  int findClustersGMSeed(const TH2 * image, MaxCamClusterImage * clust,
                           double seed_thresh, double thresh_pct,
                           double max_wrong_p, double min_thresh,
                           int blur_n, double blur_sigma, 
                           unsigned int neighbors_thresh, int min_neigh, 
                           unsigned int minsize, double minDistance, 
                           double minrxy_global, double minrxy_cluster, 
                           double maxjoin_residual, double ls_weight, const DmtpcGainMap * map, double spacer_width=2, 
                           bool reproduce_v4_bug = true); 
  
  /** Implementation in MaxCamImageTools_ClusterFinding.cc **/ 
  int findClustersGMSeedStitch(const TH2 * image, MaxCamClusterImage * clust,
                           const DmtpcStitchInfo * sinfo, 
                           const double * image_mean, const double * blurred_mean, 
                           const double * image_rms, const double * blurred_rms,  
                           double seed_thresh, double red_thresh, 
                           double thresh_pct,
                           double max_wrong_p, double min_thresh,
                           int blur_n, double blur_sigma, 
                           unsigned int neighbors_thresh, int min_neigh, 
                           unsigned int minsize, double minDistance, 
                           double minrxy_global, double minrxy_cluster, 
                           double maxjoin_residual, double ls_weight, 
                           const DmtpcGainMap * map, double spacer_width=2, const char * debug_outfile = 0, 
                           bool reproduce_v4_bug = false 
                           ); 


  /** Find clusters using GainMap. Only joins across spacers. . 
   *  Implementation in MaxCamImageTools_ClusterFinding.cc 
   */
  int findClustersGM(const TH2* image, MaxCamClusterImage* clust,
      double minsig, double maxsig, unsigned int minsize, double mindist, 
      double minrxy_global, double minrxy_cluster, double maxjoin_residual,
      double least_squares_weight, const DmtpcGainMap * map, double spacer_width = 2); 
  
  /** Implementation in MaxCamImageTools_ClusterFinding.cc **/ 
  int findClustersNG(TH2* image, MaxCamCluster** array, 
		     double min, double max, int minsize, double mindist);
  
  /** Implementation in MaxCamImageTools_ClusterFinding.cc **/ 
  int findClusterSE(TH2* image, MaxCamCluster** array, 
		    double minsig, double maxsig, int minsize, double mindist);

   TH2* makeClusterHist(TH2* image, vector<int> clusta, int padding=3);
  
  /** Implementation in MaxCamImageTools_PixelProcessing.cc **/ 
  int clustWidth(vector<int> clusta, TH2* image);
  
  /** Implementation in MaxCamImageTools_PixelProcessing.cc **/ 
  vector<int> halfway(vector<int> clusta, vector<int> clustb, TH2* image, double mindist);
  
  /** Implementation in MaxCamImageTools_PixelProcessing.cc **/ 
  double minDist(vector<int> clusta, vector<int> clustb, TH2* image);

  /** Implementation in MaxCamImageTools_ImageProcessing.cc **/ 
  TH2* blur(TH2* image, int blurn, double blurfrac);

  typedef enum
  {
    RENORMALIZE,
    MIRROR_EDGES, 
    EXTEND_EDGES,
    ZERO_EDGES
  }BLUR_EDGE_BEHAVIOR;

  /** 
   * Return the image blurred witha  gaussian kernel of size (2*blurn+1) X (2*blurn+1) and sigma sigma. 
   * The fourth parameter allows you to store the resulting kernel so that it needn't be calculated 
   * each time. 
   *
   * Note that the last specifed, non-saved kernel is cached, so repeated calls with the same parameters
   * should be quick. 
   *
   * \param image the image to blur
   * \param the number of bins on each side used in the gaussia kernel
   * \param sigma the standard deviation, in bins, of the gaussian kernel
   * \param kernel if this parameter is NULL, the kernel is calculated and not saved. 
   * If it is a pointer to a double * that is set to NULL, the kernel is calculated and 
   * saved to that location to that memory location (it is the responsibility of the user
   * to free it using the free(). If this parameter is a pointer to a double * that is not equal to NULL, then 
   * the existing kernel is used. It is assumed that it is of size 2*blurn+1 otherwise you may experience a segfault. 
   * \param edge_action specifies what to do when the kernel hits an edge
   *  
   * Implementation in MaxCamImageTools_ImageProcessing.cc
   *
   */ 
  TH2 * gaussianBlur(const TH2* image, int blurn, double sigma, double ** kernel = NULL, BLUR_EDGE_BEHAVIOR edge_action = RENORMALIZE) ;


  /** Implementation in MaxCamImageTools_ImageProcessing.cc **/ 
  TH2 * bilateralFilter(const TH2 * image, double space_sigma, double value_sigma, double nsigma=3, BILATERAL_VALUE_FN fn = BILATERAL_GAUSSIAN); 

  // same as above but faster and less precise if non-integer values (sounds to nearest integer and uses lookup table).  
  TH2 * fastBilateralFilter(const TH2 * image, double space_sigma, double value_sigma, double nsigma=3, BILATERAL_VALUE_FN fn = BILATERAL_GAUSSIAN, int scale_exp = 8, bool cache_lookup = false); 

  /** Implementation in MaxCamImageTools_ImageCalculations.cc **/ 
  int findMaximum(TH1 *project, float threshold, int minNeighbors);

  /** Implementation in MaxCamImageTools_ImageCalculations.cc **/ 
  int findMaximum2D(TH2 *project, float threshold, int minNeighbors);
  
  /** Implementation in MaxCamImageTools_PixelProcessing.cc **/ 
  int countPixelsAboveThreshold(TH2* image, float threshold);
  
  /**
   * Rotate all points in a TCutG object by an angle, with respect to a point
   *
   *  cut     The input TCutG object to be rotated (it is altered in-place)
   *  ang_rad The angle by which to rotate cut [radians]
   *  pt      The pivot point (x,y) for the rotation (in data coordinates)
   *          By default, the pivot point is the average of the extreme
   *          x and y values in the TCutG object
   *
   *  Implementation in MaxCamImageTools_ImageCalculations.cc
   *
   **/
  void rotateTCutG(TCutG *cut, float ang_rad, TPoints *pt=0);

  /**
   *  Embed the input TH2F within a larger TH2F
   * 
   *  Useful to pad before rotating to avoid cropping the image
   * 
   *  image      Input image to be padded
   *  padfactor  Factor by which to pad the image.  3 means the padded image
   *             will be 3 times larger in each dimension (9x as many pixels)
   *  CURRENTLY, THIS CODE ONLY WORKS FOR padfactor=3
   *
   *  Implementation in MaxCamImageTools_ImageTransforms.cc
   **/
  TH2* padImage(TH2* image, int padfactor=3);

  /** Implementation in MaxCamImageTools_ImageTransforms.cc */
  TH2* rotatePerfect(TH2 *hin, TString dir);

  /** Implementation in MaxCamImageTools_ImageTransforms.cc */
  TH2* rotateRight(TH2* image);

  /** Implementation in MaxCamImageTools_ImageTransforms.cc */
  TH2* rotateLeft(TH2* image);

  /** Implementation in MaxCamImageTools_ImageTransforms.cc */
  TH2* rotate180(TH2* image);


  /**
   *  Rotate an image or a subset of an image by an arbitrary amount
   *  
   *  WARNING: Currently only coded to rotate the entire image
   *  WARNING: algorithm seems to do worst when rotation is multiples 
   *           of pi/4
   *  
   *  The rotation is counter-clockwise and the pivot point is set by the
   *  input parameter "pt"
   *  
   *  Algorithm is to loop over all bins.  Sub-divide the bins into
   *  N=ndiv^2 sub-bins and then rotate each sub-bin. This spreads the
   *  "flux" in each source bin among several destination bins.  There
   *  is a tradeoff between accuracy and speed here.  Typically if ndiv
   *  = 8 or 16 the quality is fine (at least to my eye).  More
   *  detailed investigations needed to really quantify the quality of
   *  the rotation.
   *  
   *  @param[in] img  = Input image to be rotated.
   *  @param[in] ang  = Angle of counter-clockwise rotation [radians]
   *  @param[in] pt   = Point around which the rotation is applied.  
   *                    In "real" units, not bin number.
   *                    If omitted or zero, the center of the image is used.
   *  @param[in] ndiv = Each source pixel will be divided into ndiv^2 
   *                    sub-pixels.
   *  @param[in] reg  = A sub-region of the image to rotate.  If not specified
   *                    or zero, then the entire image is rotated.
   *                    [NOT YET IMPLEMENTED]
   *  @return    roti = The rotated image.
   *
   *  Implementation in MaxCamImageTools_ImageTransforms.cc
   */
  TH2* rotateImg(TH2* img, Float_t ang, TPoints* pt=0, Int_t ndiv=4, TCutG* reg=0);

  TH2 * rotateInterpolate(const TH2 *img, double ang, double xorig = -1, double yorig = -1, TH2 *out = 0, const char * interpolation = "bilinear"); 

  /**
     from the days of yore (see MaxCamRead::setRectangularCut())
   *  Implementation in MaxCamImageTools_ImageTransforms.cc
  */
  TH2* rotate(TH2* image, int nrot, TVector2 *xaxis, TVector2 *yaxis, float x0, float y0);

  /**
   *  Get the center of an image in data coordinates (not bins)
   *  Implementation in MaxCamImageTools_ImageCalculations.cc 
   **/
  TPoints* getCenter(TH2* img);

  /**
   *  Get the center of a TCutG object in data coordinates
   *  Implementation in MaxCamImageTools_ImageCalculations.cc 
   **/
  TPoints* getCenter(TCutG* cut);

  /**
   * Get the min, max values in x and y for a TCutG object
   *
   *  cut     The TCutG object to characterize
   *  xrng    (xmin, xmax) of the TCutG
   *  yrng    (ymin, ymax) of the TCutG
   *  Implementation in MaxCamImageTools_ImageCalculations.cc 
   **/
  void getExtent(TCutG* cut, TPoints* xrng, TPoints* yrng);
  

  /**
   * Count the number of pixels in an image that fall within a TCutG region
   *
   *  @param[in] img    = An image
   *  @param[in] cut    = A TCutG in that image
   *  @param[in] cutint = Include pixels interior or exterior to the cut 
   *                      (default is interior)
   *  Implementation in MaxCamImageTools_ImageCalculations.cc 
   **/
  Int_t countPix(TH2* img, TCutG* cut=0, bool cutint=true);

  /**
   *  code to handle the creation of rectangular TCutG objects
   *  Implementation in MaxCamImageTools_ImageCalculations.cc 
   **/
  TCutG* makeTCutGRect(TString name, float x0, float y0, float x1, float y1);

  /**
   * Make a annulus to approximate the ring or veto regions in the 4shooter
   * Does not assume that 
   *
   * @param[in] x0   -- center of circle, x-coordinate (data units)
   * @param[in] y0   -- center of circle, y-coordinate (data units)
   * @param[in] r1   -- Inner radius
   * @param[in] r2   -- Outer radius
   * @param[in] img  -- TH2 of image (used to get boundaries of sector)
   *
   *  Implementation in MaxCamImageTools_ImageCalculations.cc 
   **/
  TCutG* makeAnnulusSector(TString name, Float_t x0, Float_t y0, Float_t r1, Float_t r2, TH2* img);


  typedef enum
  {
    SOBEL,
    PREWITT, 
    SCHARR
  }GRADIENT_OPERATOR; 

  /** Canny edge detector. 
   *
   *  Uses nxn gaussian kernel with the gradient operator of your choice and hysteresis thresholding. 
   *  
   *  \param img An image to undergo edge detection 
   *  \param blurlevel Sigma of gaussian blur filter
   *  \param thresh_low Lower hysteresis threshold
   *  \param thresh_high Higher hysteresis threshold 
   *  \param g 3x3 gradient operator to use. 
   *  \param kernel_size The radius of the gaussian kernel to convolve with 
   *  \returns a TH2C with 0 for no edge and 1 for an edge. The histogram name is the input hist's name
   *   with a _thresh suffix. 
   *
   *  Implementation in MaxCamImageTools_ImageProcessing.cc 
   * */
  TH2C* edgeDetect(const TH2* img, double blurlevel, double thresh_low, double thresh_high, GRADIENT_OPERATOR g = SOBEL, unsigned int kernel_size = 7); 


  /** 2D Convolution
   *
   *  Creates a new image convolved with the specified kernel  
   *  
   *  \param img The image to convolve    
   *  \param kernel A collapsed array specifying the kernel. kernel[row][column] should be kernel + (width * row) + column. Note that if 
   *                the convolution is separable, the height is ignored and the kernel will be treated as one dimensional 
   *  \param width the width of the kernel 
   *  \param height the height of the kernel 
   *  \param linear treat the convolution as separable (i.e. equivalent to applying in x followed by y)
   *  \param edge_behavior what to do with the kernel by an edge
   *  \param new_name The name of the the resultant image. If NULL, the name will be the current name with _convolved appended
   *
   *  \return the convolved image
   *  Implementation in MaxCamImageTools_ImageProcessing.cc 
   */
  TH2 * convolve(const TH2* img, const double *kernel, int width, int height, bool separable = false , BLUR_EDGE_BEHAVIOR edge_behavior = RENORMALIZE, const char * new_name = 0); 
  

 /** Gradient, as used by canny edge detector
  *
  *  \param img the image to analyze
  *  \param magnitude this histogram will be filled with the magnitude of the gradient 
  *  \param orientation this histogram will be filled with the rounded axial orientation. 
  *           the orientation will be rounded to 0, 45, 90, and 135 degrees. 
  *  \param blurlevel sigma of gaussian blur filter 
  *  \param g 3x3 gradient operator to use 
  *  \param kernel_size radius of gaussian kernel to convolve with 
  *
  *  Implementation in MaxCamImageTools_ImageProcessing.cc 
  */ 
  void gradient(const TH2 * img, TH2 * magnitude, TH2S * orientation,  double blurlevel, GRADIENT_OPERATOR g = SOBEL, unsigned int kernel_size = 7); 
 
  /** Implementation in MaxCamImageTools_ImageCalculations.cc */ 
  double maxProjectedStep(const TH2 * img, char axis =0);


  typedef enum 
  {
    GAUSSIAN, 
    LORENTZIAN,
    TUKEY
  } ANISOTROPIC_DIFFUSION_FN; 

 
  /** Anisotropic Diffusion using the algorithm suggested by Perona and Malik for smoothing while preserving features. 
   *  See http://dx.doi.org/10.1109%2F34.56205 for paper 
   *
   *  \param img The image to perform the operation on 
   *  \param lambda Lambda parameter. Related to amount of diffusion. 0 < lambda <=1 
   *  \param K K parameter, related to strength of edge sensitivity
   *  \param fn The edge sensitivity function to use. Right now, choices are Gaussian, Lorentzian and Tukey Bilevel. 
   *  \param gradient_sigma the blur amount for computing the sigma 
   *  \param gradop The gradient operator for performing the gradient 
   *
   *  Implementation in MaxCamImageTools_ImageProcessing.cc 
   **/
  TH2* anisotropicDiffusion(const TH2* img, double lambda, double K, ANISOTROPIC_DIFFUSION_FN fn = TUKEY, double gradient_sigma = 1, GRADIENT_OPERATOR gradop = SOBEL, bool use_diagonals = false); 

  //TH2F* convertFitsIntoROOT(TString fitsFile, TString histName="imageHisto");
  /** Implementation in MaxCamImageTools_IO.cc */ 
  TH2* convertFitsIntoROOT(TString fitsFile, TString histName="imageHisto", int imageNumber=0, char type='F');

  /** Implementation in MaxCamImageTools_IO.cc */ 
  int getNumberOfImagesInFitsFile(TString fitsFile);

  /** Implementation in MaxCamImageTools_IO.cc */ 
  const char* getTimeStamp(TString fitsFile);

  /** Implementation in MaxCamImageTools_IO.cc */ 
  time_t getTimeStampUnix(TString fitsFile);

  /** Implementation in MaxCamImageTools_IO.cc */ 
  char* getFitsHeaderValue(TString fitsFile, char* key);

  /** Implementation in MaxCamImageTools_IO.cc */ 
  // Returns NULL if keyword not found in FITS file 
  // Does NOT quit FITS functions fail (status)
  char* getFitsHeaderValueIfExists(TString fitsFile, char* key);
  
  /** Implementation in MaxCamImageTools_IO.cc */ 
  int convertIntoFits(TH2* image, TString fileName);
  



  /** Implementation in MaxCamImageTools_ClusterFinding.cc */ 
   int findClustersADHysteresisGM(const TH2 *image, MaxCamClusterImage * clust, 
                             double K, double Lambda, ANISOTROPIC_DIFFUSION_FN f, int ntimes, double gradient_sigma, 
                             GRADIENT_OPERATOR g, double high_thresh, double low_thresh, 
                             unsigned int min_neighbors, unsigned int minsize, 
                             double mindist, 
                             double minrxy_global, double minrxy_cluster, double maxjoin_residual,
                             double least_squares_weight, const DmtpcGainMap * map, double spacer_width = 2); 

 
  /** Implementation in MaxCamImageTools_ImageTransforms.cc */ 
  TH2* resizeImage(TH2* image, float xmin, float xmax, float ymin, float ymax);

  /** Implementation in MaxCamImageTools_ImageTransforms.cc */ 
  TH1* resizeImage(TH1* image, float xmin, float xmax);

  /** Implementation in MaxCamImageTools_ImageTransforms.cc */ 
  TH1*  resizeImage(TH1* image, TH1* resize);
  
  /** Implementation in MaxCamImageTools_ImageCalculations.cc */ 
  void principalAxes(TH2 *image, float threshold, float &Ix, float &Iy);
  
  /** Implementation in MaxCamImageTools_ImageCalculations.cc */ 
  double cosRecoil2D( TH2 *image, int imax, int jmax);
  
  /** Implementation in MaxCamImageTools_ImageCalculations.cc */ 
  TH1F* makeProfile(TGraph *graph, int nx, float minx, float maxx, float miny, float maxy, TCanvas *c=0);
  
  /** Implementation in MaxCamImageTools_ClusterFinding.cc */ 
  double findPedestalWithTracks(TH2* image, MaxCamCluster** array, int ntracks);

  TH2F *rmsScale(const TH2 * image, const DmtpcStitchInfo * sinfo, const double * image_rms, const double * image_mean); 

  /**
   *  Return the mean value of the pixels in an image
   *
   *  @param[in] image  - Image to be analyzed
   *  @param[in] cut    - TCutG to be applied to image
   *  @param[in] cutint - Include pixels interior or exterior to the cut 
   *                      (default is interior)
   *
   *  Implementation in MaxCamImageTools_ImageCalculations.cc 
   **/
  double getMean(const TH2* image, TCutG* cut=0, bool cutint=true);
  /**
   *  Return the RMS value of the pixels in an image
   *
   *  @param[in] image  - Image to be analyzed
   *  @param[in] cut    - TCutG to be applied to image
   *  @param[in] cutint - Include pixels interior or exterior to the cut 
   *                      (default is interior)
   *
   *  Implementation in MaxCamImageTools_ImageCalculations.cc 
   *
   **/
  double getRMS(const TH2* image, TCutG* cut=0, bool cutint=true);

  /** Implementation in MaxCamImageTools_ImageCalculations.cc */ 
  void meanRMSNoOutliers(const TH2* image, double& mean, double& rms);

  /** Implementation in MaxCamImageTools_ImageCalculations.cc */ 
  TH1F* pixelHist(TH2* image, Int_t nbins=100, Float_t minVal=9999, Float_t maxVal=-9999, TCutG* cut=0, bool cutint=true);
  
  /** Implementation in MaxCamImageTools_ImageCalculations.cc */ 
  int distanceToImageEdge(TH2 *image, int ibin);

  /** Implementation in MaxCamImageTools_ImageProcessing.cc */ 
  TH2* TH2Fkappasigmaclip(DmtpcDataset* d, int cameraNum, int start, int stop, float sigmas, int iter, float change);

  /** Performs a hough transform on the input image looking for a line 
   * 
   *  The output format is a 2-dimensional histogram of line votes where
   *  the X axis is the angle of the line and the Y axis is the distance from 0,0 of the line (y = -cot(#theta) + r / sin(#theta))
   *
   *  You may pass an output histogram (assumed to be filled with zeros) with the range and binning you want or one may 
   *  be generated if you pass 0 with nangbins (ranging from 0 to 2 pi) and nrbins (ranging from 0 to the diagonal distance). 
   *  A TH2I or a TH2D will be generated depending on if the input is float or integer type. 
   *
   *   \param in the histogram to transform
   *   \param out either a histogram to use for the output (that defines the parameter ranges) or 0 to allocate a new one. 
   *   \param min The minimum threshold in the input image to add a curve in hough space
   *   \param nangbins if we allocate a new histogram, this is the number of angle bins, otherwise this is ignored
   *   \param nrbins if we allocate a new histogram, this is the number of r bins, otherwise this is ignored
   *   \return the hough transform of the image, with angle on the x axis and distance from smallest corner on the y axis
   *
   *   Implementation in MaxCamImageTools_ImageTransforms.cc 
   */ 
  TH2* houghTransformLine(TH2* in, TH2* out = 0, double min = 1,  unsigned nangbins = 512, unsigned nrbins=512); 

  // order: xcenter, ycenter, minor, major, rotation_angle (phi)
  typedef enum 
  {
    ELL_X0,ELL_Y0,ELL_A,ELL_B,ELL_PHI
  } ELLIPSE_PARAMETERS; 

  /** Naive implementation of ellipse hough transform. Use may be impractical due to storage
   *  and computation requirements. The accumulator is implemented as a 5-dimensional sparse 
   *  histogram with B (one of the axes) as the "independent" coordinate (the others will always be filled
   *  if there exists a valid value of b)
   *
   *  \warning This hasn't been tested yet...
   *
   *  The output array must be prepared ahead of time with the number of bins and limits for each parameters. 
   *  The parameter order is defined the ELLIPSE_PARAMETERS enum. For example, you may do:
   *
   *  double nbins[5], min[5], max[5]; 
   *  nbins[ELL_X0] = ...  
   *  min[ELL_X0] = ...  
   *  max[ELL_X0] = ...  
   *  ... 
   *  THnSparse * out = new THnSparse{F,D,I,S,C}("hist","hist",5,nbins,min,max); 
   *  houghTransformEllipse(in,out); 
   *
   *  \param in The histogram to find ellipses in
   *  \param out a prepared output 5d sparse histogram that includes the number of bins and limits you want
   *  \param min The minimum input histogram value to fill out a manifold in hough space
   *  \return pointer to out (the accumulator)
   *
   *   Implementation in MaxCamImageTools_ImageTransforms.cc 
   */
  THnSparse * houghTransformEllipse(TH2* in, THnSparse * out, double min = 1); 


  /** Circle hough transform. Accumulator is implemented as a 3-dimensional histogram with
   *  r as the "independent" coordinate. 
   *
   * The output histogram must be prepared ahead of time with the number of bins and limits for each paramter.
   * The parameter order (x,y,z) is x0, y0, r. 
   *
   *  ... 
   *  TH3 * out = new TH3IF,D,I,S,C}("hist","hist",nbinsx,xmin,xmax,nbinsy,ymin,ymax,nbinsr,rmin,rmax);  
   *  houghTransformCircle(in,out); 
   *
   *  \param in the histogram to find cirlces in 
   *  \param out the prepared output 3d 
   *  \param min The minimum input histogram value to fill out a manifold in hough space
   *  \return pointer to out (the accumulator) 
   *
   *   Implementation in MaxCamImageTools_ImageTransforms.cc 
   */
  TH3 * houghTransformCircle(TH2* in, TH3 * out, double min = 1, void (* progressFn)(int row) = 0); 

  /** Suppress all values that are not local maxima 
   * 
   * \param in input histogram
   * \param n  the number of pixels in neighborhood to check
   * \param setval the value to set nonmaximal values to 
   * \param suffix the suffix appended to the name of the original
   * \return a suppressed histogram
   *
   *  Implementation in MaxCamImageTools_ImageProcessing.cc 
   * */
  TH2* nonMaximumSuppress(const TH2 * in, int n =1, double setval =0, const char * suffix = "_suppressed"); 



  /** Implementation in MaxCamImageTools_ImageProcessing.cc*/
  TH2 * medianFilter(const TH2 * in, unsigned width =3, unsigned niter = 1); 

  /** Implementation in MaxCamImageTools_ImageProcessing.cc*/
  TH2 * laplacian (const TH2 * in); 
  /** Implementation in MaxCamImageTools_ImageProcessing.cc*/
  TH2 * laplacianOfGaussian (const TH2 * in, double sigma, int width); 

  /** Implementation in MaxCamImageTools_ImageTransforms.cc*/
  TH2 * toPolarCoordinates(TH2 * in, unsigned nrbins = 256, unsigned nthetabins = 256, double centerx = 0,
                           double centery = 0, char type = 'D', const char * interpolation = "bilinear", const char * name = 0);

  /** Build a composite image (TH2) from multiple images (TH2's). 
   *  
   *  @param in array of pointers to histograms to use
   *  @param info DmtpcStitchInfo with stitching parameters 
   *  @param interpolation the interpolation method. Currently only bilinear,bicubic and nearest neighbor are supported. 
   *  @param histType the output histogram type. 'C',S','I','F' and 'D' are supported. 
   *
   *   Implementation in MaxCamImageTools_ImageTransforms.cc 
   */
  TH2* montage(const TH2 * const* in, DmtpcStitchInfo * info, const char * name, 
               const char * interpolation = "bilinear", char histType ='F');

 /** Implementation in MaxCamImageTools_PixelProcessing.cc*/
 int selectPixelsAbove(const TH2 * img, std::set<int> * pixels, double thresh); 

 /** Implementation in MaxCamImageTools_PixelProcessing.cc*/
 int selectPixelsBelow(const TH2 * img, std::set<int> * pixels, double thresh); 

 /** Implementation in MaxCamImageTools_PixelProcessing.cc*/
 int selectPixelsBetween(const TH2 * img, std::set<int> * pixels, double thresh_low, double thresh_high); 

 /** 
  * Adds all pixels in the image but not inside in to out
  *
  * Returns number of pixels added to out; 
  *
  * Implementation in MaxCamImageTools_PixelProcessing.cc*/
 int getInversePixels(const TH2 * img, const std::set<int> * in, std::set<int> * out);  

 /** Adds all pixels that are at least partiall contained within the box to the set *
  *
  * Implementation in MaxCamImageTools_PixelProcessing.cc 
  */ 
 int selectPixelsInsideBox(const TH2 * img, std::set<int> * pixels, double xmin, double ymin, double xmax, double ymax); 

 /** Adds all pixels whose centers are contained within the circle to the set
  *
  * Implementation in MaxCamImageTools_PixelProcessing.cc 
  */ 
 int selectPixelsInsideCircle(const TH2 * img, std::set<int> * pixels, double x0, double y0, double r); 

 int selectPixelsAboveWithAtLeastNNeighbors(const TH2 * img, std::set<int> * pixels, double thresh, int n); 
 int getNeighborsAboveThreshold(const TH2 * img, int bin, std::set<int> * pixels, double thresh); 


 /** Implementation in MaxCamImageTools_PixelProcessing.cc*/
 int dilate(const TH2 *img, std::set<int> * pixels); 

 /** Implementation in MaxCamImageTools_PixelProcessing.cc*/
 int erode(const TH2 *img, std::set<int> * pixels); 

 /** Implementation in MaxCamImageTools_PixelProcessing.cc*/
 void fillPixels(TH2 *img, const std::set<int> * pixels, double val); 

 /** Implementation in MaxCamImageTools_PixelProcessing.cc*/
 int innerBorder(const TH2 *img, const std::set<int> * pixels, std::set<int> * out); 

 /** Implementation in MaxCamImageTools_PixelProcessing.cc*/
 int outerBorder(const TH2 *img, const std::set<int> * pixels, std::set<int> * out); 

 /** Implementation in MaxCamImageTools_ImageCalculations.cc*/
 double rectIntersectionArea(double * xcoords1, double * ycoords1, double * xcoords2, double * ycoords2); 

 /** Implementation in MaxCamImageTools_ImageCalculations.cc*/
 double polygonArea(int n, const double * x, const double * y, bool ordered = true); 

 /** Assumes same size bins!!! */ 
 /** Implementation in MaxCamImageTools_Interpolation.cc*/
 double interpolateBicubic(const TH2* in, double x, double y); 

 /** Implementation in MaxCamImageTools_ImageProcessing.cc*/
 TH2C * binaryThreshold (const TH2 * in, double thresh); 

 /** Implementation in MaxCamImageTools_PixelProcessing.cc*/ 
 TH2I * distanceTransform(const TH2* in, const set<int> * pixels); 

 /** Implementation in MaxCamImageTools_ImageProcessing.cc*/
 void fillEdges(TH2 * in, unsigned width = 1, double setval = 0); 

 /** Interpolates a value in the given histogram with the given method. 
  *
  * Current methods are :
  *   bilinear
  *   bicubic
  *   nearest 
  *
  * Use testInterpolate to check a method. 
  *
  * Implementation in MaxCamImageTools_Interpolation.cc
  */ 
 double interpolate(const TH2* in, double x, double y, const char * method = "bilinear"); 

 /** Implementation in MaxCamImageTools_Interpolation.cc*/
 bool testInterpolate(const char * method); 

 /** The output of each pixel is the pixel value divided  (or subtracted) by the mean (or median) of the neighbor values. Implementation in MaxCamImageTools_ImageProcessing **/
 TH2 * neighborRatio(const TH2 * in, bool abs = false, bool median = false, bool difference = false); 



 /** 
  * 
  *
  * Implementation in MaxCamImageTools_ImageProcessing
  *
  * **/
 int killLonePixelsMedianDifference(TH2 * tokill, double threshold); 

 void projectAlongLine (const TH2 * img, TH1 ** longi, TH1 ** tranv, double x0, double y0, double x1, double y1, double width, const char * interpolation_method = "bilinear"); 

 TH2D * radonTransform (const TH2 * img, int nbinstheta = 360, int nprojbins = -1);


 /** Computes the median of the image **
  *
  * Implemented in MaxCamImageTools_ImageCalculations.cc 
  */ 
 double median(const TH2 * img); 

 TH2 * crop(const TH2 * img, int xbinmin, int xbinmax, int ybinmin, int ybinmax); 

 TH2 * WienerGaussDeconvolve(const TH2 * in, double gauss_sigma, double noise_sigma); 

 /** Implementation in MaxCamImageTools_ImageTransforms.cc  **/ 
 // These all modify the input, and then return it for convenience
 
 TH2* hist_apply (TH2 * in, double (*fn) (double)); 
 TH2* hist_apply (TH2 * in, double (*fn) (double,double), double arg); 
 TH2* hist_sqrt (TH2 * in); 
 TH2* hist_abs (TH2* in);
 TH2* hist_pow (TH2* in, double b);

 /** Shifts the quadrants of a 2D FFT so that the 0 frequency bin is in the center
  * */
 TH2 * fftshift (TH2 * in, double scalex = 1, double scaley = 1, bool normalize = false); 


 /** Implementation in MaxCamImageTools_ImageProcessing **/ 
 /* Add noise from distribution in  f */ 
 void addGaussianNoise(TH2 * in, double sigma = 1); 


 
 //defined in MaxCamImageTools_ImageTransforms.cc
 //Zero pad the input histogram to the size newx, newy; 
 //The newxpos and newypos define the position of the original within the new image
 //where -1 = left align, 1 = right align, 0 = centered
 TH2* zeroPad(TH2 * in, int newx, int newy, float newxpos = 0, float newypos = 0, const char * name = "_zeropad"); 

 TH2* zeroPadSquare(TH2 * in, float newpos = 0, const char * name = "_zeropad"); 

 int countUniqueValues(const TH1 * in); 

};

#endif
