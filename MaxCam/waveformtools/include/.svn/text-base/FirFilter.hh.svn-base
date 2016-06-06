#ifndef __FIR_FILTER_H__
#define __FIR_FILTER_H__
#include <vector>
class TH1;
namespace waveform{

/**
  \brief A class to perform finite-impulse response (FIR) filters
  on one-dimensional data. 

  The classes FirFilter<float>, FirFilter<double>
  and FirFilter<long double> are compiled automatically. The filter 
  works on single and double precision input. These are provided rather
  than relying on templates in order to get everything to work with CINT.
  
  \section FIR Finite-Impulse Response Filters
  Finite-impulse response (or FIR) filters are a simple type of digital filter
  used to extract local properties of a discrete signal. They can be used to
  smooth data using the values of a set of nearby points and can perform derivatives,
  as well as other functions. Am FIR filter is defined as a series of coefficients
  (we will call this the kernel and denote its elements by \f$ \{b[n]\}\f$ ) and the 
  filtering action is a discrete convolution of an input signal and the FIR kernel.
  Given an input signal x, the output signal is given by
  \f[ y[n] = \sum\limits_{m=0}^{N_k-1}b[m]x[n+m] = b[N_k-1-n] \ast x[n] \f] 
  
  Note that the order of the coefficients is the reverse of the definition
  used in some other sources such as Wikipedia. This definition uses the actual
  array indices used in the calculation. The last equality is ignoring some time
  shifting since y[n] should really be y[N-1+n] for the convolution definition to
  be correct.

  \section ZT Z-Transforms and the Transfer Function  
 
  A filter is typically characterized by its transfer function, related to its z-transform.
  The z-transform of a digital signal is given by
  \f[ X(z) = Z\{x[n]\} = \sum\limits_{-\infty}^\infty x[n] z^{-n} \f]
  The z-transform is the discrete analog to the Laplace transform. 

  One useful property of FIR filters is the convolution identity
  \f[ Z(x[n]\ast y[n]) = Z(x[n])Z(y[n]) \f]
  Applying this to our definition of the filter, 
  \f[Y(z) = X(z) \sum\limits_{k=0}^{N-1} z^{-k} b[N-1-k] \f]
  We then define the transfer function to be
  \f[ T(z) = \frac{Y(z)}{X(z)} = \sum\limits_{k=0}^{N-1} z^{-k} b[N-1-k] \f]
  
  \section FR Frequency Response
  The frequency response is a special case of the z-transform that is particularly useful in 
  characterizing the filter for actual signals. It is just the transfer function
  of the filter on the unit circle, \f$ z = \exp(2\pi i w) \f$, where w is the
  frequency. This is equivalent to taking a discrete Fourier transform. We let
 \f[ F(w) = T(e^{2\pi iw}) = F_r(w)+ i F_i(w) = \sum\limits_{k=0}^{N-1} \cos (2\pi kw) b[N-1-k] -i \sum\limits_{k=0}^{N-1} \sin(2\pi kw) b[N-1-k] \f]

  The power attenuation gives the ratio of the output signal power to the input signal power and is given by
  \f[ P(w) = F_r(w)^2 + F_i(w)^2 \f]
  The magnitude attenuation is the ratio of the output signal magnitude to the input signal magnitude and is given by
  \f[ M(w) = \sqrt{P(w)} \f]

  The phase shift gives the phase change between the input and output signals. 
  \f[\phi(w) = \arg(F(w)) = \tan^{-1}\frac{F_i}{F_r}  \f]
  
  The phase and group delays are additional ways to measure time delays or phase shifts between the input and output signals.
  the phase delay is 
  \f[\tau_{\phi} = -\frac{\phi(w)}{2\pi w}\f]
  and the group delay is defined as
  \f[\tau_g = -\frac{1}{2\pi}\frac{d\phi}{d\omega}\f]
  If we let \f$F_{r,n} = \sum\limits{k=0}^{N-1} k\cos(2\pi kw) b[N-1-k]\f$ and
  \f$F_{i,n} = -\sum\limits_{k=0}^{N-1} k\sin(2\pi kw)b[N-1-k]\f$ then the group delay
  can be shown to be
  \f[\tau_g = \frac{1}{1+\left(\frac{F_i}{F_r}\right)^2}\left(\frac{F_{r,n}}{F_r}+\frac{F_{i,n}F_i}{F_r^2}\right)\f]

  \author Jeremy Lopez
  \date August 2012
 */
template <class REALTYPE>
class FirFilter
{
  public:
    /** Default Constructor */
    FirFilter(){}
    /** 
     * Constructs a filter object with kernel size n. If ker is not null,
     * its contents are copied to the fKernel
     * \param n The kernel size
     * \param ker an array holding the kernel
     */
    FirFilter(int n,REALTYPE ker[]=0);
    /**
     * Constructs a filter using the contents of ker as the kernel
     * \param ker the vector to be copied into fKernel
     */
    FirFilter(std::vector<REALTYPE> ker);
    /** Destructor */
    virtual ~FirFilter(){}

    /** Returns the current size of the output vector. This is equal to 
      * \f$N_{out} = N_{in} - N_{k} + 1 \f$. If no input has
      * been specified, the size will be 0.
      */
    int output_size() const {return fOutput.size();}
    /** Returns the size of the kernel (number of coefficients used in the filter */
    int kernel_size() const {return fKernel.size();}
    /** Calculates the output size for a given input size. This is given by
      * \f$N_{out} = N_{in} - N_{k} + 1 \f$/
      * \param insize the input size
      * \return the output array size
      */
    int calculate_output_size(int insize) const {return insize + 1 - kernel_size();}

    /** Sets a new kernel
      * \param n the size of the new kernel
      * \param ker an array holding the new coefficients
      */
    void set_kernel(int n, const REALTYPE* ker);
    /** Sets a new kernel
      * \param ker a vector holding the new coefficients
      */
    void set_kernel(std::vector<REALTYPE> ker){set_kernel(ker.size(),&ker[0]);}

    /** Get a reference to the output vector */
    std::vector<REALTYPE>& output() {return fOutput;}
    /** Get a const reference to the output vector */
    const std::vector<REALTYPE>& output() const{return fOutput;}
    /** Get a pointer to the output array */
    REALTYPE* getOutput(){return &fOutput[0];}   
    /** Get a const pointer to the output array */
    const REALTYPE* getOutput() const{return &fOutput[0];}

    /** Normalize the filter coefficients. This divides all coefficients by
      * \f[ B =  \sum\limits_0^{N_k-1} b_n \f]
      * so that \f$ 1 = \sum b_n\f$. Be careful to avoid calling this 
      * on filters that are not meant to preserve signal integral.
      * Many basic derivative filters may have \f$\sum b_n = 0 \f$ and this
      * function will generate floating point errors.
      */
    void normalize();

//Float
    /** Copy the output into a float array
      * \param out the array
      */
    void copy_output(float* out) const;
    /** Run the filter on a float array and save to the output vector fOutput.
      * \param n the input size
      * \param in the input array
      */
    void run(int n,float in[]);
    /** Run the filter on a float vector and save to the output vector fOutput.
      * \param in the input vector
      */
    void run(std::vector<float> in) {return run(in.size(),&in[0]);}
    /** Run the filter on a float array and save to a given output array.
      * Note that the output vector is not filled in this function.
      * \param n the input size
      * \param in the input array
      * \param out the output array
      */
    void run(int n,float in[],float out[]) const;
    /** Run the filter on a float vector and save to a given output vector.
      * Note that the output vector fOutput, saved internally in the class is
      * not filled in this function.
      * \param in the input vector
      * \param out the output vector
      */
    void run(std::vector<float> in, std::vector<float> out) {return run(in.size(),&in[0],&out[0]);}

//Double
    /** Copy the output into a float array
      * \param out the array
      */
    void copy_output(double* out) const;
    /** Run the filter on a float array and save to the output vector fOutput.
      * \param n the input size
      * \param in the input array
      */
    void run(int n,double in[]);
    /** Run the filter on a double vector and save to the output vector fOutput.
      * \param in the input vector
      */
    void run(std::vector<double> in) {return run(in.size(),&in[0]);}
    /** Run the filter on a double array and save to a given output array.
      * Note that the output vector is not filled in this function.
      * \param n the input size
      * \param in the input array
      * \param out the output array
      */
    void run(int n,double in[],double out[]) const;
    /** Run the filter on a double vector and save to a given output vector.
      * Note that the output vector fOutput, saved internally in the class is
      * not filled in this function.
      * \param in the input vector
      * \param out the output vector
      */
    void run(std::vector<double> in,std::vector<double> out) const{return run(in.size(),&in[0],&out[0]);}

//Transforms
    /** Calculates the transfer function of the filter
      * \param zr the real part of the argument
      * \param zi the imaginary part of the argument
      * \param ztr the real part of the transfer function
      * \param zti the imaginary part of the transfer function
      */
    void transfer_function(double zr,double zi,double& ztr,double& zti) const;
    /** Calculates the frequency response of the filter. The frequency is in
      * inverse bins so that \f$w=0.5\f$ is the Nyquist frequency, and the response
      * is the same for all integral w.
      * \param w the frequency, in inverse bins.
      * \param re the real part of the frequency response
      * \param im the imaginary part of the frequency response
      */
    void frequency_response(double w, double& re,double& im) const;
    /** Calculate the attenuation in signal power at a given frequency. See
      * frequency_response() for an explanation of the frequency units. 
      * \param w the frequency
      * \return the attenuation in signal power
      */
    double power_attenuation(double w) const;
    /** Calculate the attenuation in signal magnitude at a given frequency. See
      * frequency_response() for an explanation of the frequency units. 
      * \param w the frequency
      * \return the attenuation in signal magnitude
      */
    double magnitude_attenuation(double w) const;
    /** Calculate the signal phase shift at a given frequency. See
      * frequency_response() for an explanation of the frequency units. 
      * \param w the frequency
      * \return the phase shift
      */
    double phase_shift(double w) const;
    /** Calculate the signal phase delay at a given frequency. See
      * frequency_response() for an explanation of the frequency units. 
      * \param w the frequency
      * \return the phase delay
      */
    double phase_delay(double w) const;
    /** Calculate the signal group delay at a given frequency. See
      * frequency_response() for an explanation of the frequency units. 
      * \param w the frequency
      * \return the group delay
      */
    double group_delay(double w) const; 

  protected:
    /** The FIR coefficients */
    std::vector<REALTYPE> fKernel; 
    /** The output vector */
    std::vector<REALTYPE> fOutput;
};


}


//Definitions of template functions
#ifndef __CINT__
#include <cmath>
#ifndef M_PI
  #define M_PI 3.14159265358979323846L
#endif
namespace
{
  const static long double PI = M_PI;
}

template <class REALTYPE>
waveform::FirFilter<REALTYPE>::FirFilter(int n, REALTYPE* kernel) : fKernel(n)
{
  if (kernel!=0)
    for (int i = 0; i < n; i++)
      fKernel[i] = kernel[i];
}

template <class REALTYPE>
waveform::FirFilter<REALTYPE>::FirFilter(std::vector<REALTYPE> kernel) : fKernel(kernel.size())
{
  for (unsigned int i = 0; i<kernel.size(); i++)
    fKernel[i] = kernel[i];
}

template <class REALTYPE>
void
waveform::FirFilter<REALTYPE>::set_kernel(int n, const REALTYPE* k)
{
  fKernel.resize(n);
  for (int i = 0; i<n; i++)
    fKernel[i] = k[i];


}

template <class REALTYPE>
void
waveform::FirFilter<REALTYPE>::normalize()
{
  double norm = 0;
  for (int i = 0 ; i<(int)fKernel.size(); i++)
    norm+=fKernel[i];
  for (int i = 0 ; i<(int)fKernel.size(); i++)
    fKernel[i]/=norm;
}

template <class REALTYPE>
void
waveform::FirFilter<REALTYPE>::transfer_function(double zr,double zi,double& ztr, double& zti) const
{
  ztr = 0;
  zti = 0;
  if (zr==0&&zi==0) return;
  double norm = sqrt(zr*zr+zi*zi);
  double phi = atan2(zi,zr);
  double power = 1;
  for (int i = 0 ; i<(int) fKernel.size(); i++)
  {
    ztr += fKernel[fKernel.size()-1-i]/power*cos(i*phi);
    zti -= fKernel[fKernel.size()-1-i]/power*sin(i*phi);
    power/=norm;
  }

}

template <class REALTYPE>
void
waveform::FirFilter<REALTYPE>::frequency_response(double w, double& re,double& im) const
{

  re = 0; 
  im = 0; 
  for (int i =0; i <(int) fKernel.size(); i++)
  {
    re += fKernel[fKernel.size()-1-i] * cos(2*::PI*i*w);
    im -= fKernel[fKernel.size()-1-i] * sin(2*::PI*i*w);

  }
}

template<class REALTYPE>
double
waveform::FirFilter<REALTYPE>::power_attenuation(double w) const
{
  double r,i;
  frequency_response(w,r,i);
  return r*r+i*i;
}

template<class REALTYPE>
double
waveform::FirFilter<REALTYPE>::magnitude_attenuation(double w) const
{return sqrt(power_attenuation(w));}

template<class REALTYPE>
double
waveform::FirFilter<REALTYPE>::phase_shift(double w) const
{
  double r,i;
  frequency_response(w,r,i);
  return atan2(i,r);
}

template<class REALTYPE>
double
waveform::FirFilter<REALTYPE>::phase_delay(double w) const
{
  if (w==0) return 0;
  return -phase_shift(w)/w;
}

template<class REALTYPE>
double
waveform::FirFilter<REALTYPE>::group_delay(double w) const
{
  if (w==0) return 0;
  double Cf=0,Cfn=0,Sf=0,Sfn=0;
  for (unsigned int i = 0; i < fKernel.size(); i++)
  {
    Cf += fKernel[fKernel.size()-1-i]*cos(2*::PI*w);
    Cfn += i*fKernel[fKernel.size()-1-i]*cos(2*::PI*w);
    Sf -= fKernel[fKernel.size()-1-i]*sin(2*::PI*w);
    Sfn -= i*fKernel[fKernel.size()-1-i]*sin(2*::PI*w);
  }
  double Rf = Sf/Cf;
  return -(Cfn/Cf+Rf*Sfn/Cf)/(1+Rf*Rf);
}

template <class REALTYPE> 
void
waveform::FirFilter<REALTYPE>::copy_output(float out[]) const
{
  for (int i = 0; i < output_size(); i++)
    out[i] = fOutput[i];
}


template <class REALTYPE>
void
waveform::FirFilter<REALTYPE>::run(int n, float in[])
{
  int outsize = n + 1 - fKernel.size();
  if (outsize<1)
  {
    fOutput.resize(0);
    return;
  }else if (outsize!=(int)fOutput.size())
    fOutput.resize(outsize);

  for (int i = 0; i < (int)fOutput.size(); i++)
  {
    fOutput[i] = 0;
    for (int j = 0; j < (int)fKernel.size(); j++)
    {
      fOutput[i] += fKernel[j] * in[i+j];
    }
  }
}


template <class REALTYPE>
void
waveform::FirFilter<REALTYPE>::run(int n, float in[], float out[]) const
{
  int outsize = n + 1 - fKernel.size();
    
  for (int i = 0; i < outsize; i++)
  {
    out[i] = 0;
    for (int j = 0; j < (int)fKernel.size(); j++)
    {
      out[i] += fKernel[j] * in[i+j];
    }
  }
}

template <class REALTYPE> 
void
waveform::FirFilter<REALTYPE>::copy_output(double out[]) const
{
  for (int i = 0; i < output_size(); i++)
    out[i] = fOutput[i];
}


template <class REALTYPE>
void
waveform::FirFilter<REALTYPE>::run(int n, double in[])
{
  int outsize = n + 1 - fKernel.size();
  if (outsize<1)
  {
    fOutput.resize(0);
    return;
  }else if (outsize!=(int)fOutput.size())
    fOutput.resize(outsize);

  for (int i = 0; i < (int)fOutput.size(); i++)
  {
    fOutput[i] = 0;
    for (int j = 0; j < (int)fKernel.size(); j++)
    {
      fOutput[i] += fKernel[j] * in[i+j];
    }
  }
}


template <class REALTYPE>
void
waveform::FirFilter<REALTYPE>::run(int n, double in[], double out[]) const
{
  int outsize = n + 1 - fKernel.size();
    
  for (int i = 0; i < outsize; i++)
  {
    out[i] = 0;
    for (int j = 0; j < (int)fKernel.size(); j++)
    {
      out[i] += fKernel[j] * in[i+j];
    }
  }
}
#endif


#endif
