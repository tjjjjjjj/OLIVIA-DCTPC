#ifndef __IIRFILTER_HH__
#define __IIRFILTER_HH__
#include <vector>
class TH1;
namespace waveform{

template <class A>
class IirFilter
{
  public:
    IirFilter(){}
    IirFilter(int nf,int ni,const A fir[]=0,const A iir[] =0,const A init[] = 0);
    IirFilter(std::vector<A> fir,std::vector<A> iir, std::vector<A> init);
    virtual ~IirFilter(){}

    int output_size() const {return fOutput.size();}
    int fir_size() const {return fFir.size();}
    int iir_size() const {return fIir.size();}
    int initial_vals_size() const {return fInit.size();}

    int calculate_output_size(int insize) const {return insize + 1 - fir_size();}

    void set_filter(int nf,const A* fir, int ni, const A* iir, const A* init=0);
    void set_filter(const std::vector<A>& fir, const std::vector<A>& iir,const std::vector<A>& init);
    void set_fir(int n,const A* fir);
    void set_iir(int n,const A* fir);
    void set_init(int n,const A* init);
    void set_fir(const std::vector<A>& fir){set_fir(fir.size(),&fir[0]);}
    void set_iir(const std::vector<A>& iir){set_iir(iir.size(),&iir[0]);}
    void set_init(const std::vector<A>& init){set_init(init.size(),&init[0]);}

    const std::vector<A>& output() const{return fOutput;}
    void output(A* out){out = &fOutput[0];}
    A* getOutput(){return &fOutput[0];}   
    const A* getOutput() const{return &fOutput[0];}

    void normalize();

//Float
    void copy_output(float* out) const;
    void run(int n,const float in[]);
    void run(const std::vector<float>& in) {return run(in.size(),&in[0]);}
    void run(int n,const float in[],float out[]) const;
    void run(const std::vector<float>& in, std::vector<float>& out) const 
      {return run(in.size(),&in[0],&out[0]);}

    void run_reverse(int n,const float in[]);
    void run_reverse(const std::vector<float>& in) {return run(in.size(),&in[0]);}
    void run_reverse(int n,const float in[],float out[]) const;
    void run_reverse(const std::vector<float>& in, std::vector<float>& out) const
      {return run_reverse(in.size(),&in[0],&out[0]);}


//Double
    void copy_output(double* out) const;
    void run(int n,const double in[]);
    void run(const std::vector<double>& in) {return run(in.size(),&in[0]);}
    void run(int n,const double in[],double out[]) const;
    void run(const std::vector<double>& in, std::vector<double>& out) {return run(in.size(),&in[0],&out[0]);}

    void run_reverse(int n,const double in[]);
    void run_reverse(const std::vector<double>& in) {return run(in.size(),&in[0]);}
    void run_reverse(int n,const double in[],double out[]) const;
    void run_reverse(const std::vector<double>& in, std::vector<double>& out) const
      {return run_reverse(in.size(),&in[0],&out[0]);}


//Transforms
    void transfer_function(double zr,double zi,double& ztr,double& zti) const;
    void frequency_response(double w, double& re,double& im) const;
    double magnitude_attenuation(double w) const;
    double power_attenuation(double w) const;
    double phase_shift(double w) const;
    double phase_delay(double w) const;
    double group_delay(double w) const;

  protected:
    std::vector<A> fFir;
    std::vector<A> fIir;
    std::vector<A> fInit;
    std::vector<A> fOutput;
};

//  TH1* createIirOutputHist(TH1* in, int nkernel);
//  TH1* createIirOutputHist(int nin, int nkernel, char type);
//  void IirHistFilter(TH1* in, TH1* out,IirFilter<double>* fir);

}


//Definitions of template functions
#ifndef __CINT__
#include <cmath>
#ifndef M_PI
  #define M_PI 3.14159265358979323846L


namespace
{
  const static long double PI = M_PI;
}
#else
#include "FirFilter.hh"
#endif

template <class A>
waveform::IirFilter<A>::IirFilter(int nf, int ni,const A fir[],const A iir[],const A init[]) : fFir(nf), fIir(ni), fInit(ni-1)
{
  if (fir!=0)
    for (int i = 0; i < nf; i++)
      fFir[i] = fir[i];

  if (iir!=0)
    for (int i = 0; i < ni; i++)
      fIir[i] = iir[i];

  for (int i = 0; i < ni-1; i++)
    fInit[i] = init ? init[i] : 0;    

}

template <class A>
waveform::IirFilter<A>::IirFilter(std::vector<A> fir,std::vector<A> iir, std::vector<A> init) : 
  fFir(fir.size()), fIir(iir.size()), fInit(iir.size()-1)
{
  for (unsigned int i = 0; i<fir.size(); i++)
    fFir[i] = fir[i];
  for (unsigned int i = 0; i<iir.size(); i++)
    fIir[i] = iir[i];
  for (unsigned int i = 0; i<init.size(); i++)
    fInit[i] = init[i];
}

template <class A>
void
waveform::IirFilter<A>::set_filter(int nf, const A* fir, int ni, const A* iir, const A* init)
{
  if (fir!=0)
    for (int i = 0; i < nf; i++)
      fFir[i] = fir[i];

  if (iir!=0)
    for (int i = 0; i < ni; i++)
      fIir[i] = iir[i];

  for (int i = 0; i < ni-1; i++)
    fInit[i] = init ? init[i] : 0;    

}

template <class A>
void
waveform::IirFilter<A>::set_filter(const std::vector<A>& fir,const std::vector<A>& iir,const std::vector<A>& init)
{
  fFir.resize(fir.size());
  fIir.resize(iir.size());
  fInit.resize(iir.size()-1,(A)0);
  for (unsigned int i = 0; i<fir.size(); i++)
    fFir[i] = fir[i];
  for (unsigned int i = 0; i<iir.size(); i++)
    fIir[i] = iir[i];
  for (unsigned int i = 0; i<init.size(); i++)
    fInit[i] = init[i];


}

template <class A>
void
waveform::IirFilter<A>::set_fir(int n, const A* k)
{
  fFir.resize(n);
  for (int i = 0; i<n; i++)
    fFir[i] = k[i];
}

template <class A>
void
waveform::IirFilter<A>::set_iir(int n, const A* k)
{
  fIir.resize(n);
  for (int i = 0; i<n; i++)
    fIir[i] = k[i];
}

template <class A>
void
waveform::IirFilter<A>::set_init(int n, const A* k)
{
  fInit.resize(n);
  for (int i = 0; i<n; i++)
    fInit[i] = k==0?0:k[i];
}

template <class A>
void
waveform::IirFilter<A>::normalize()
{
  //Normalize FIR part to unity
  double norm = 0;
  for (int i = 0; i<(int) fFir.size(); i++)
    norm+=fFir[i];
  for (int i = 0; i<(int) fFir.size(); i++)
    fFir[i]/=norm;

}

template <class A>
void
waveform::IirFilter<A>::transfer_function(double zr,double zi,double& ztr, double& zti) const
{
  ztr=0, zti=0;
  if (zr==0&&zi==0) return;
  double norm = sqrt(zr*zr+zi*zi);
  double phi = atan2(zi,zr);
  double Nr=0,Ni=0;
  double power = 1;
  for (int i = 0 ; i<(int) fFir.size(); i++)
  {
    Nr += fFir[fFir.size()-1-i]/power*cos(i*phi);
    Ni -= fFir[fFir.size()-1-i]/power*sin(i*phi);
    power *= norm;
  }
  double Dr=0,Di=0;
  power = 1;
  for (int i = 0; i<(int) fIir.size();i++)
  {
    Dr += fIir[i]/power*cos(i*phi);
    Di -= fIir[i]/power*sin(i*phi);
    power *= norm;
  }

  ztr = (Dr*Nr+Di*Ni)/sqrt(Dr*Dr+Di*Di);
  zti = (Ni*Dr-Nr*Di)/sqrt(Dr*Dr+Di*Di);


}

template <class A>
void
waveform::IirFilter<A>::frequency_response(double w, double& re,double& im) const
{
  double Nr=0,Ni=0, Dr=0,Di=0;
  for (int i =0; i< (int) fFir.size(); i++)
  {
    Nr += fFir[fFir.size()-1-i] * cos(2*::PI*i*w);
    Ni -= fFir[fFir.size()-1-i] * sin(2*::PI*i*w);
  }
  Nr /= fIir[0];
  Ni /= fIir[0];
  for (int i = 0; i < (int) fIir.size(); i++)
  {
    Dr += fIir[i]*cos(2*::PI*i*w);
    Di -= fIir[i]*sin(2*::PI*i*w);
  }
  Dr /= fIir[0];
  Di /= fIir[0];
  
  re = (Nr*Dr + Ni*Di)/sqrt(Dr*Dr+Di*Di);
  im = (Ni*Dr - Nr*Di)/sqrt(Dr*Dr+Di*Di); 
}

template<class A>
double
waveform::IirFilter<A>::power_attenuation(double w) const
{
  double re,im;
  frequency_response(w,re,im);
  return re*re+im*im;
}

template<class A>
double
waveform::IirFilter<A>::magnitude_attenuation(double w) const
{return sqrt(power_attenuation(w));}

template<class A>
double
waveform::IirFilter<A>::phase_shift(double w) const
{
  double re,im;
  frequency_response(w,re,im);
  return atan2(im,re);
}

template<class A>
double
waveform::IirFilter<A>::phase_delay(double w) const
{
  if (w==0) return 0;
  return -phase_shift(w)*2*::PI/w;
}

template<class A>
double
waveform::IirFilter<A>::group_delay(double w) const
{
  if (w==0) return 0;
  double Cf=0,Cfn=0,Sf=0,Sfn=0;
  for (unsigned int i = 0; i < fFir.size(); i++)
  {
    Cf += fFir[fFir.size()-1-i]*cos(2*::PI*w);
    Cfn += i*fFir[fFir.size()-1-i]*cos(2*::PI*w);
    Sf -= fFir[fFir.size()-1-i]*sin(2*::PI*w);
    Sfn -= i*fFir[fFir.size()-1-i]*sin(2*::PI*w);
  }
  double Ci=0,Cin=0,Si=0,Sin=0;
  for (unsigned int i = 0; i < fIir.size(); i++)
  {
    Ci += fIir[i]*cos(2*::PI*w);
    Cin += i*fIir[i]*cos(2*::PI*w);
    Si -= fIir[i]*sin(2*::PI*w);
    Sin -= i*fIir[i]*sin(2*::PI*w);
  }
  double Rf = Sf/Cf;
  double Ri = Si/Ci;

  return (Cin/Ci + Ri * Sin/Ci)/(1+Ri*Ri) - (Cfn/Cf+Rf*Sfn/Cf)/(1+Rf*Rf);
  

}

template <class A> 
void
waveform::IirFilter<A>::copy_output(float out[]) const
{
  for (int i = 0; i < output_size(); i++)
    out[i] = fOutput[i];
}


template <class A>
void
waveform::IirFilter<A>::run(int n, const float in[])
{
  int outsize = n + 1 - fFir.size();
  if (outsize<1)
  {
    fOutput.resize(0);
    return;
  }else if (outsize!=(int)fOutput.size())
    fOutput.resize(outsize);

  for (int i = 0; i < (int)fOutput.size(); i++)
  {
    fOutput[i] = 0;
    //FIR part
    for (int j = 0; j < (int)fFir.size(); j++)
    {
      fOutput[i] += fFir[j] * in[i+j];
    }
    //IIR part
    for (int j = 1; j <(int) fIir.size(); j++)
    {
      double x = fIir[j];
      x*= j>i? fInit[fIir.size() - (j-i)]: fOutput[i-j];
      fOutput[i] -= x;
    }
    fOutput[i] /= fIir[0];

  }
}

template <class A>
void
waveform::IirFilter<A>::run_reverse(int n, const float in[])
{
  int outsize = n + 1 - fFir.size();
  if (outsize<1)
  {
    fOutput.resize(0);
    return;
  }else if (outsize!=(int)fOutput.size())
    fOutput.resize(outsize);

  int outmax = outsize-1;
  int inmax = n-1;
  //Loop over all output bins
  for (int i = 0; i < outsize; i--)
  {
    fOutput[outmax-i] = 0;
    //FIR part
    for (int j = 0; j < (int)fFir.size(); j++)
    {
      fOutput[outmax-i] += fFir[j] * in[inmax-i-j];
    }
    //IIR part
    for (int j = 1; j <(int) fIir.size(); j++)
    {
      double x = fIir[j];
      x*= j>i? fInit[fIir.size() - (j-i)]: fOutput[outmax-i+j];
      fOutput[outmax-i] -= x;
    }
    fOutput[outmax-i] /= fIir[0];

  }
}



template <class A>
void
waveform::IirFilter<A>::run(int n, const float in[], float out[]) const
{
  int outsize = n + 1 - fFir.size();
    
  for (int i = 0; i < outsize; i++)
  {
    out[i] = 0;
    //FIR part
    for (int j = 0; j < (int)fFir.size(); j++)
    {
      out[i] += fFir[j] * in[i+j];
    }
    //IIR part
    for (int j = 1; j <(int) fIir.size(); j++)
    {
      double x = fIir[j];
      x*= j>i? fInit[fIir.size() - (j-i)]: out[i-j];
      out[i] -= x;
    }
    out[i] /= fIir[0];
  }
}

template <class A>
void
waveform::IirFilter<A>::run_reverse(int n, const float in[],float out[]) const
{
  int outsize = n + 1 - fFir.size();

  int outmax = outsize-1;
  int inmax = n-1;
  //Loop over all output bins
  for (int i = 0; i < outsize; i--)
  {
    out[outmax-i] = 0;
    //FIR part
    for (int j = 0; j < (int)fFir.size(); j++)
    {
      out[outmax-i] += fFir[j] * in[inmax-i-j];
    }
    //IIR part
    for (int j = 1; j <(int) fIir.size(); j++)
    {
      double x = fIir[j];
      x*= j>i? fInit[fIir.size() - (j-i)]: out[outmax-i+j];
      out[outmax-i] -= x;
    }
    out[outmax-i] /= fIir[0];

  }
}



template <class A> 
void
waveform::IirFilter<A>::copy_output(double out[]) const
{
  for (int i = 0; i < output_size(); i++)
    out[i] = fOutput[i];
}


template <class A>
void
waveform::IirFilter<A>::run(int n, const double in[])
{
//Set the output array if not already done
//Return if output size is zero
  int outsize = n + 1 - fFir.size();
  if (outsize<1)
  {
    fOutput.resize(0);
    return;
  }else if (outsize!=(int)fOutput.size())
    fOutput.resize(outsize);

//Loop over all output bins
  for (int i = 0; i < (int)fOutput.size(); i++)
  {
    fOutput[i] = 0;
    //FIR part
    for (int j = 0; j < (int)fFir.size(); j++)
    {
      fOutput[i] += fFir[j] * in[i+j];
    }
    //IIR part
    for (int j = 1; j <(int) fIir.size(); j++)
    {
      double x = fIir[j];
      x*= j>i? fInit[fIir.size() - (j-i)]: fOutput[i-j];
      fOutput[i] -= x;
    }
    fOutput[i] /= fIir[0];

  }
}


template <class A>
void
waveform::IirFilter<A>::run(int n,const double in[], double out[]) const
{
  int outsize = n + 1 - fFir.size();
    
  for (int i = 0; i < outsize; i++)
  {
    out[i] = 0;
    //FIR part
    for (int j = 0; j < (int)fFir.size(); j++)
    {
      out[i] += fFir[j] * in[i+j];
    }
    //IIR part
    for (int j = 1; j <(int) fIir.size(); j++)
    {
      double x = fIir[j];
      x*= j>i? fInit[fIir.size() - (j-i)]: out[i-j];
      out[i] -= x;
    }
    out[i] /= fIir[0];
  }
}

template <class A>
void
waveform::IirFilter<A>::run_reverse(int n, const double in[])
{
  int outsize = n + 1 - fFir.size();
  if (outsize<1)
  {
    fOutput.resize(0);
    return;
  }else if (outsize!=(int)fOutput.size())
    fOutput.resize(outsize);

  int outmax = outsize-1;
  int inmax = n-1;
  //Loop over all output bins
  for (int i = 0; i < outsize; i--)
  {
    fOutput[outmax-i] = 0;
    //FIR part
    for (int j = 0; j < (int)fFir.size(); j++)
    {
      fOutput[outmax-i] += fFir[j] * in[inmax-i-j];
    }
    //IIR part
    for (int j = 1; j <(int) fIir.size(); j++)
    {
      double x = fIir[j];
      x*= j>i? fInit[fIir.size() - (j-i)]: fOutput[outmax-i+j];
      fOutput[outmax-i] -= x;
    }
    fOutput[outmax-i] /= fIir[0];

  }
}
template <class A>
void
waveform::IirFilter<A>::run_reverse(int n, const double in[],double out[]) const
{
  int outsize = n + 1 - fFir.size();

  int outmax = outsize-1;
  int inmax = n-1;
  //Loop over all output bins
  for (int i = 0; i < outsize; i--)
  {
    out[outmax-i] = 0;
    //FIR part
    for (int j = 0; j < (int)fFir.size(); j++)
    {
      out[outmax-i] += fFir[j] * in[inmax-i-j];
    }
    //IIR part
    for (int j = 1; j <(int) fIir.size(); j++)
    {
      double x = fIir[j];
      x*= j>i? fInit[fIir.size() - (j-i)]: out[outmax-i+j];
      out[outmax-i] -= x;
    }
    out[outmax-i] /= fIir[0];

  }
}


#endif


#endif
