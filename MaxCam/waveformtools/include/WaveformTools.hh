/**\file WaveformTools.hh
\author Jeremy Lopez
\date 19 April 2011

Header file for signal processing methods

*/

#ifndef DMTPCTRACETOOLS
#define DMTPCTRACETOOLS
class TH1F;
class TH2F;
class TObjArray;
class WaveformVector;
class TH1; 
#include "IirFilter.hh"
#include "FirFilter.hh"
namespace waveform{




/**\brief Contains methods for signal processing

  This namespace contains various methods used in signal processing of 1-dimensional scope traces.
  This is meant to contain only methods that transform the signal in some way tomake it easier to analyze.

*/
  namespace tools{


  //Convolution
  /**Return the convolution of a histogram with the given array.  This also normalizes the array to 1.  The convolution is
     of the form
 
     \f$\tilde H(x_i) = (H\ast f)(x_i)=\left(\sum\limits_{k=-N}^{N}{f(x_k)}\right)^{-1} \sum\limits_{j=-N}^{N}{H(x_{i+j})f(x_j)} \f$

     where \f$H(x_i)\f$ is the histogram, \f$\tilde H(x_i)\f$ is the convolution,  and \f$f(x_i)\f$ is the convolution kernel of length \f$2N+1\f$

     \param hist the histogram
     \param N the half-width of the array (total length 2N+1)
     \param arr the array to convolve
     \param removeEdges if true removes the N bins on the edges of the histogram to remove all edge effects
     \param name the name and title of the new histogram

  */
  TH1F* convolution(const TH1F* hist, int N, const double* arr, 
		    bool removeEdges = true, const char* name = "convolution");

  /**Return the convolution of a histogram with the given array without doing any normalization

     \f$\tilde H(x_i) = H\ast f(x_i)=\sum\limits_{j=-N}^{N}{H(x_{i+j})f(x_j)} \f$

     where \f$H(x_i)\f$ is the histogram, \f$\tilde H(x_i)\f$ is the convolution,  and \f$f(x_i)\f$ is the convolution kernel of length \f$2N+1\f$
     \param hist the histogram
     \param N the half-width of the array (total length 2N+1)
     \param arr the array to convolve
     \param removeEdges if true removes the N bins on the edges of the histogram
     \param name the name and title of the new histogram
  */
  TH1F* convNoNorm(const TH1F* hist, int N, const double* arr, 
		   bool removeEdges=true,const char* name="convolution");
  //Window smoothing

  /**Convolve with a function of the form

   \f$f(x_i) =  a_0+a_1\cos\frac{i\pi}{N}+a_2\cos\frac{2i\pi}{N}+a_3\cos\frac{3i\pi}{N}+a_4\cos\frac{4i\pi}{N}\f$

   where \f$i\in[-N,N]\f$.

   \param hist the histogram
   \param N the width of the kernel
   \param a0 the constant term
   \param a1 the coefficient of the first cosine term
   \param a2 the coefficient of the second cosine term
   \param a3 the coefficient of the third cosine term
   \param a4 the coefficient of the fourth cosine term
   \param rmEdge if true remove edge bins
   \param name the new histogram name
  */
  TH1F* raisedCosWin(const TH1F* hist, int N, double a0, double a1, 
			double a2, double a3, double a4, bool rmEdge = true,
			const char* name="raisedCosWin");
  /**Convolve with the Hann window.  This is the raised cosine (raisedCosWin()) window with

     \f$a_i = {0.5,0.5,0,0,0}\f$
     \param hist the histogram
     \param N the width of the convolution function
     \param rmEdge if true remove edge bins
     \param name the new histogram name
  */
  TH1F* hannWin(const TH1F* hist, int N, bool rmEdge = true, 
		const char* name="hannWindow");
  /**Convolve with the Hamming window.  This is the raised cosine (raisedCosWin()) window with

     \f$a_i = {0.54,0.46,0,0,0}\f$
     \param hist the histogram
     \param N the width of the convolution function
     \param rmEdge if true remove edge bins
     \param name the new histogram name
  */
  TH1F* hammingWin(const TH1F* hist, int N, bool rmEdge=true, 
		      const char* name="hammingWindow");
  /**Convolve with the Blackman window.  This is the raised cosine (raisedCosWin()) window with

     \f$a_i = {0.5(1-\alpha),0.5,0.5\alpha,0,0}\f$
     \param hist the histogram
     \param N the width of the convolution function
     \param alpha the window parameter
     \param rmEdge if true remove edge bins
     \param name the new histogram name
  */
  TH1F* blackmanWin(const TH1F* hist, int N, double alpha=0.16, 
		    bool rmEdge=true, 
		    const char* name = "blackmanWindow");
  /**Convolve with the Blackman-Harris window.  This is the raised cosine (raisedCosWin()) window with

     \f$a_i = {0.35875,0.48829,0.14128,0.01168,0}\f$
     \param hist the histogram
     \param N the width of the convolution function
     \param rmEdge if true remove edge bins
     \param name the new histogram name
  */
  TH1F* blackmanHarrisWin(const TH1F* hist, int N, bool rmEdge=true, 
			  const char* name="blackmanHarris");
  /**Convolve with the Blackman-Nutall window.  This is the raised cosine (raisedCosWin()) window with

     \f$a_i = {0.3635819,0.4891775,0.1365995,0.0106411,0}\f$
     \param hist the histogram
     \param N the width of the convolution function
     \param rmEdge if true remove edge bins
     \param name the new histogram name
  */
  TH1F* blackmanNutallWin(const TH1F* hist, int N, bool rmEdge=true, 
			  const char* name="blackmanNutall");
  /**Convolve with the Flat top window.  This is the raised cosine (raisedCosWin()) window with

     \f$a_i = {1,1.93,1.29,0.388,0.032}\f$
     \param hist the histogram
     \param N the width of the convolution function
     \param rmEdge if true remove edge bins
     \param name the new histogram name
  */
  TH1F* flatTopWin(const TH1F* hist, int N, bool rmEdge=true,
		   const char* name="flatTopeWindow");
  /**Convolve with the Nutall window.  This is the raised cosine (raisedCosWin()) window with

     \f$a_i = {0.355768,0.487396,0.144232,0.012604,0}\f$
     \param hist the histogram
     \param N the width of the convolution function
     \param rmEdge if true remove edge bins
     \param name the new histogram name
  */
  TH1F* nutallWin(const TH1F* hist, int N, bool rmEdge=true,
		  const char* name="nutallWindow");
  /**Convolve with the Bohman window:

     \f$f(x_k) = \left(1-\left|\frac{k}{N}\right|\right)\cos\frac{k\pi}{N} + \frac{1}{\pi}\sin\frac{k\pi}{N}\f$
     \param hist the histogram
     \param N the width of the convolution function
     \param rmEdge if true remove edge bins
     \param name the new histogram name
  */
  TH1F* bohmanWin(const TH1F* hist, int N, bool rmEdge=true,
		  const char* name="bohmanWindow");
  /**Convolve with the Kaiser window:

     \f$f(x_k) = \frac{I_0\left(\pi a\sqrt{1-k^2/N^2}\right)}{I_0(\pi a)}\f$
     where \f$I_0(x) \f$ is the 0th order modified Bessel function of the first kind.
     \param hist the histogram
     \param N the width of the convolution function
     \param a the parameter
     \param rmEdge if true remove edge bins
     \param name the new histogram name
  */
  TH1F* kaiserWin(const TH1F* hist,int N,double a=3, bool rmEdge=true,
		  const char* name="kaiserWindow");
  /**Convolve with the Lanczos window:

     \f$f(x_k) = \frac{N}{k\pi}\sin\frac{k\pi}{N}\f$
     \param hist the histogram
     \param N the width of the convolution function
     \param rmEdge if true remove edge bins
     \param name the new histogram name
  */
  TH1F* lanczosWin(const TH1F* hist, int N, bool rmEdge=true,
		   const char* name="lanczosWindow");
  /**Convolve with the Lanczos filter:

     \f$f(x_k) = \frac{N}{k^2\pi^2}\sin\frac{k\pi}{N}\sin k\pi\f$
     \param hist the histogram
     \param N the width of the convolution function
     \param rmEdge if true remove edge bins
     \param name the new histogram name
  */
  TH1F* lanczosFilter(const TH1F* hist, int N, bool rmEdge=true,
		      const char* name="lanczosFilter");

  /**Convolve with the cosine window:

     \f$f(x_k) = \cos\frac{k\pi}{2N}\f$
     \param hist the histogram
     \param N the width of the convolution function
     \param rmEdge if true remove edge bins
     \param name the new histogram name
  */
  TH1F* cosineWin(const TH1F* hist, int N, bool rmEdge=true,
		  const char* name = "cosineWindow");
  /**Convolve with a rectangular window of half-width N.
     \param hist the histogram
     \param N the width of the convolution function
     \param rmEdge if true remove edge bins
     \param name the new histogram name
  */
  TH1F* rectWin(const TH1F* hist, int N, bool rmEdge=true,
		const char* name="rectWindow");
  /**Convolve with a Gaussian window of half-width N.

     \f$f(x_k) = \exp\frac{-k^2}{2\sigma^2}\f$
     \param hist the histogram
     \param N the width of the convolution function
     \param sigma the width in bins
     \param rmEdge if true remove edge bins
     \param name the new histogram name
  */
  TH1F* gausWin(const TH1F* hist, int N, double sigma,bool rmEdge=true,
		const char* name="gausWin");
  /**Convolve with a triangular window of half-width N.
     \param hist the histogram
     \param N the width of the convolution function
     \param rmEdge if true remove edge bins
     \param name the new histogram name
  */
  TH1F* triangleWin(const TH1F* hist, int N, bool rmEdge=true,
		    const char* name="triWin");
  /**Convolve with an exponential window of half-width N.

     \f$f(x_k) = \exp\frac{-|k|}{x_0}\f$
     \param hist the histogram
     \param N the width of the convolution function
     \param x0 the width in bins
     \param rmEdge if true remove edge bins
     \param name the new histogram name
  */
  TH1F* expoWin(const TH1F* hist, int N, double x0, bool rmEdge=true,
		const char* name="expoWin");

  //Convolution smoothing
  /**Convolve with a Gaussian function (Gaussian low-pass filter in frequency domain)

     \f$f(x_k) =\int\limits_{k-0.5}^{k+0.5}\!\!dx \exp\frac{-x^2}{2\sigma^2}\f$
     \param hist the histogram
     \param sigma the width of the convolution function
     \param cut the value of the integrand at which to cut off the convolution kernel
     \param rmEdge if true remove edge bins
     \param name the new histogram name
  */
  TH1F* gausConv(const TH1F* hist, double sigma, bool rmEdge=true, 
		 double cut=0.01, const char* name="gausConv");
  /**Convolve with a triangle function (\f$\sin^2(x)/x^2\f$ low-pass filter in frequency domain)

     \f$f(x_k) =\int\limits_{k-0.5}^{k+0.5}{\!\!dx\left(1-\frac{|x|}{N}\right) }\f$
     \param hist the histogram
     \param N the width of the convolution function
     \param rmEdge if true remove edge bins
     \param name the new histogram name
  */
  TH1F* triangleConv(const TH1F* hist, int N,bool rmEdge=true,
		     const char* name="triConv");
  /**Convolve with an exponential function (Breit-Wigner low-pass filter in frequency domain)

     \f$f(x_k) =\int\limits_{k-0.5}^{k+0.5}\!\!dx \exp\frac{-|x|}{x_0}\f$
     \param hist the histogram
     \param x0 the width of the convolution function
     \param cut the value of the integrand at which to cut off the convolution kernel
     \param rmEdge if true remove edge bins
     \param name the new histogram name
  */
  TH1F* expoConv(const TH1F* hist, double x0,bool rmEdge=true, 
		 double cut=0.01, const char* name="expoConv");
  /**Convolve with a Fermi-Dirac distribution. In the frequency domain, this reduces to a Breit-Wigner low pass filter for very negative x0/t and a sinc filter for very small t and positive x0. 

     \f$f(x_k) =\int\limits_{k-0.5}^{k+0.5}\!\!dx\left(\exp\frac{x-x_0}{t}+1\right)^{-1}\f$
     \param hist the histogram
     \param x0 the cutoff value of the filter
     \param t the width parameter of the filter
     \param cut the value of the integrand at which to cut off the convolution kernel
     \param rmEdge if true remove edge bins
     \param name the new histogram name
  */
  TH1F* fermiConv(const TH1F* hist, double x0,double t, bool rmEdge=true,
		  double cut=0.01,const char* name="fermiConv");

  //Other smoothing
  /**Use a Savitzky-Golay filter (polynomial interpolation).  This is done using a least squares fit in the following manner:

    Suppose we wish to have an nth order polynomial fit to 2N+1 points.  
   Then we define \f$a_k\f$ to be the coefficients we wish to solve for and \f$y_k\f$ to be the values, 
 equally spaced in \f$x_k\f$ and centered about \f$y_0\f$.  
We can then define the matrix \f$A_{jk} = j^k\f$, where \f$j\in[-N,N], k\in[0,n], 0^0=1\f$.  The coefficients may then be obtained through the following relation:

   \f$a_k = (A^T A)^{-1}_{ki}A^T_{ij}y_j \f$
     
   Here, \f$a_0\f$ is the interpolated value of the histogram at \f$x_0\f$ and \f$k!a_k\f$ is the kth order derivative at that point.
     \param hist the histogram
     \param order the order of the polynomial
     \param size the half-width of the filter
     \param rmEdge if true, remove edge bins
     \param name the new histogram name
  */
  TH1F* savitzkyGolay(const TH1F* hist, int order, int size, 
		      bool rmEdge=true,const char* name = "SGFilter");

  /**Use a median filter.  This takes 2size+1 points around each bin and replaces that bin with the median value of those points.
    \param hist the histogram
    \param size the half-width of the filter
    \param name the new histogram name
  */
  TH1F* medianFilter(const TH1F* hist, int size, 
		     const char* name="medianFilter");


  //Other functions
  /**Numerically integrate the histogram using linear interpolation
    \param hist the histogram
    \param useXwidth if true, include the width of each bin in determining the integral
    \param name the new histogram name
  */
  TH1F* integ(const TH1F* hist, bool useXwidth = false,
	      const char* name ="Integral");
  /**Use the Savitzky-Golay smoothing algorithm, savitzkyGolay(), to calculate derivatives. This is also known as Lanczos differentiation.
    \param hist the histogram
    \param order the order of the polynomial
    \param size the half-width of the window
    \param nder the number of the derivative to calculate
    \param rmEdge if true remove edge bins
    \param name the new histogram name
  */
  TH1F* sgDeriv(const TH1F* hist, int order, int size, 
		int nder = 1,bool rmEdge =true, const char* name = "Derivative");

  /**Compute the value of the unipolar z-transform at a given point on the complex plane.  If rectangular coordinates are used, x and a are the real parts and y and b are imaginary.  Otherwise, x and a are magnitutes and y and b are arguments.
    \param hist the histogram
    \param x one part of the complex number 
    \param y the other part of the complex number
    \param a part of the value of the z-transform
    \param b the other part of the value of the z-transform
    \param isReIm if true use rectangular coordinates, else use polar
  */
  void zTransform(const TH1F* hist, double x, double y, 
		  double& a, double& b, bool isReIm=true);
    
  /**Add a WaveformVector object to a TObjArray
    \param wfvlist a pointer to the TObjArray 
    \param index the index in the TObjArray at which to put the WaveformVector
    \param name the name to give the WaveformVector instance
    \param type the type of WaveformVector that corresponds to this channel
  */
  void addWaveformVectorToTObjArray(TObjArray * wfvlist, int index, 
				    const char * name, const char * type);
  
  /** 
   * Compute the IIR coefficients for a (2*n) order band stop butterworth filter. 
   * \param w0 the center of the stopband in normalized frequency (w/wnyquist)
   * \param delta width of the stopband in normalized frequency (w/wnyquist)
   * \param the order of the filter (actually a band stop is two filters, so it comes out to 2*n)
   * \param acoeffs pointer to double *, memory will be allocated for this array
   * \param bcoeffs pointer to double *, memory will be allocated for this array
   */
  void butterworthBandStop(double w0, double delta, unsigned int n, double ** acoeffs, double ** bcoeffs);

  /** Compute an IIR filter of the given histogram. acoeffs and bcoeffs should be arrays of size order+1
   *  \param hist the histogram to compute the filter of 
   *  \param order the order of the filter
   *  \param acoeffs array of a coefficients 
   *  \param bcoeffs array of b coefficients 
   *  \return the filtered histogram
   * */ 
  TH1 *IIRFilter(const TH1* hist, unsigned int order, const double * acoeffs, const double * bcoeffs); 

  /** Compute an FIR filter of the given histogram. bcoeffs should be array of size order+1
   *  \param hist the histogram to compute the filter of 
   *  \param order the order of the filter
   *  \param bcoeffs array of b coefficients 
   *  \return the filtered histogram
   * */ 
  TH1 *FIRFilter(const TH1* hist, unsigned int order, const double * bcoeffs); 



  //FIR Filter Tools
  /** Create an output histogram for an FIR filter given the input histogram and a kernel size.
   *  The type is either double for TH1D or float for all others
   *  \param in The input histogram
   *  \param nkernel The size of the convolution kernel
   *  \return An output histogram with the correct number of bins, axis titles, and bin widths/positions
   * */
  TH1* createFirOutputHist(TH1* in, int nkernel);
  /** Create an output histogram for an FIR filter given the number of input bins, a kernel size, and a
   *  histogram type. The type should be "C" for char, "S" for short, "I" for int, "F" for float and "D"
   *  for double. In practice only float or double should be used.
   * \param nin the size of the input histogram
   * \param nkernel The size of the convolution kernel
   * \param type the type of output to be used
   * \return An output histogram with the correct number of bins. 
   * */
  TH1* createFirOutputHist(int nin, int nkernel,char type);
  /** Use an FIR filter on a ROOT histogram. If out points to NULL, a new 
      histogram will be created.
      \param in The input histogram
      \param out The output histogram
      \param fir The filter object
      \return The output (either out or a new histogram if out is null)
   * */
  TH1* FirHistFilter(TH1* in, TH1* out, FirFilter<double>* fir);

  //IIR Filter Tools
  /** Create an output histogram for an IIR filter given the input histogram and a kernel size.
   *  The type is either double for TH1D or float for all others
   *  \param in The input histogram
   *  \param nkernel The number of FIR coefficients in the filter
   *  \return An output histogram with the correct number of bins, axis titles, and bin widths/positions
   * */
  TH1* createIirOutputHist(TH1* in, int nkernel);
  /** Create an output histogram for an IIR filter given the number of input bins, a kernel size, and a
   *  histogram type. The type should be "C" for char, "S" for short, "I" for int, "F" for float and "D"
   *  for double. In practice only float or double should be used.
   * \param nin the size of the input histogram
   * \param nkernel The number of FIR coefficients in the filter
   * \param type the type of output to be used
   * \return An output histogram with the correct number of bins. 
   * */
  TH1* createIirOutputHist(int nin,int nkernel, char type);
  /** Use an IIR filter on a ROOT histogram. If out points to NULL, a new 
      histogram will be created.
      \param in The input histogram
      \param out The output histogram
      \param iir The filter object
      \return The output (either out or a new histogram if the output is null
   * */
  TH1* IirHistFilter(TH1* in, TH1* out, IirFilter<double>* iir);

  //uses libsamplerate for now
  //resample_ratio is input/output
  //float precision is used by libsamplerate, so if you need double precision you're out of luck. 
  TH1 * resample(const TH1 *in, double resample_ratio = 2., TH1* out = 0); 


  unsigned zeroSuppress(TH1 *in, double zero_level, double noise_ceiling, unsigned min, unsigned window_size, double * mean = 0, double * rms =0); 

  TH1* cropZeros(const TH1* in, TH1* out = 0); 

  }

}

#endif
