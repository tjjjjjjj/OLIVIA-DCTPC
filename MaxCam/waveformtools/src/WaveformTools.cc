

#include "WaveformVector.hh"
#include "SkimWaveform.hh"
#include "FastWfVector.hh"
#include "FastWaveform.hh"
#include "FastPulse.hh"
#include "CspWfVector.hh"
#include "CspWaveform.hh"
#include "CspPulse.hh"
#include "PMTWfVector.hh"
#include "PMTWaveform.hh"
#include "PMTPulse.hh"
#include "DmtpcRootTools.hh" 
#include "WaveformTools.hh"
#include "FirFilter.hh"
#include "IirFilter.hh"
#include "TClass.h"
#include "TH1.h"
#include "TH2.h"
#include "TObjArray.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMath.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <complex> 
#include <stdlib.h>
#ifdef HAVE_LIBSAMPLERATE
#include <samplerate.h>
#endif
using namespace std;
#define PI 3.14159026535898
#define get(h,i) h->GetArray()[i]
#define set(h,i,x) h->GetArray()[i]=x
//#define get2D(h,i,j) h->GetArray()[i*(h->GetNbinsX()+2)+j]
//#define set2D(h,i,j,x) h->GetArray()[i*(h->GetNbinsX()+2)+j]=x
using namespace waveform;
/*
=============================================
Basic convolution algorithm
Used by most functions in this namespace

=============================================
*/


TH1F*
tools::convolution(const TH1F* hist, int N, const double* arr, bool removeEdges,const char* name)
{
  int nbins = hist->GetNbinsX();
  double xmax = hist->GetXaxis()->GetXmax();
  double xmin = hist->GetXaxis()->GetXmin();
  
  TH1F* filter = 0;
  if (removeEdges)
  {
    double xsize = hist->GetXaxis()->GetBinWidth(1); 
    filter = new TH1F(name,name,nbins-2*N-2, xmin+xsize*(N+1),xmax-xsize*(N-1));
    double sum = 0;

    for (int i = 0; i<2*N+1; i++){
      sum+= arr[i];
    }

    double value;
    for (int i = N+2; i<=nbins-N-1; i++)
    {
      value = 0;
      for (int j = i-N; j<=i+N; j++)
        value+= arr[j-i+N] * get(hist,j);

      set(filter,i-N-1,value/sum);
    }

  }else{
   
    filter = new TH1F(name,name,nbins, xmin,xmax);

    double value, sum;
    int minBin, maxBin;
    for (int i = 1; i<= nbins; i++)
    {
      value = 0;
      minBin = max(2,i-N);
      maxBin = min(nbins-1,i+N);
      sum = 0;
      for (int j = minBin; j<=maxBin; j++){
        value += arr[j-i+N]*get(hist,j);
        sum += arr[j-i+N];
      }
      if (sum==0) set(filter,i,0);
      else set(filter,i,value/sum);
    }



  }
  return filter;
}

TH1F*
tools::convNoNorm(const TH1F* hist, int N, const double* arr, bool removeEdges,const char* name)
{
  int nbins = hist->GetNbinsX();
  double xmax = hist->GetXaxis()->GetXmax();
  double xmin = hist->GetXaxis()->GetXmin();
  
  TH1F* filter = 0;
  if (removeEdges)
  {
    double xsize = hist->GetXaxis()->GetBinWidth(1); 
    filter = new TH1F(name,name,nbins-2*N-2, xmin+xsize*(N+1),xmax-xsize*(N-1));

    double value;
    for (int i = N+2; i<=nbins-N-1; i++)
    {
      value = 0;
      for (int j = i-N; j<=i+N; j++)
        value+= arr[j-i+N] * get(hist,j);

      set(filter,i-N-1,value);
    }

  }else{

    double value;
    int minBin, maxBin;
    for (int i = 1; i<= nbins; i++)
    {
      value = 0;
      minBin = max(2,i-N);
      maxBin = min(nbins-1,i+N);
      if (minBin==2){
        for (int j = 2; j<=2+2*N; j++)
          value+=arr[j-1+N]*get(hist,j);

      }else if (maxBin==nbins-1){
        for (int j = maxBin-2*N; j<=maxBin; j++)
	  value+=arr[j-i+N]*get(hist,j);

      }else{
        for (int j = minBin; j<=maxBin; j++)
          value += arr[j-i+N]*get(hist,j);
      }
      set(filter,i,value);
    }



  }
  return filter;
}

/*
==================================
Basic window convolutions

===================================
*/
TH1F*
tools::raisedCosWin(const TH1F* hist, int N,double a0, double a1, double a2, double a3, double a4,bool rmEdge, const char* name)
{
  //Set the convolution kernel
  //Default: a = 0.16
  double kernel[1+2*N];
  for (int i = -N; i<=N; i++){
      kernel[i+N] =
       a0 + a1 * cos(PI*i/N) + a2 * cos(2*PI*i/N)
       + a3*cos(3*PI*i/N) + a4*cos(4*PI*i/N);
  }

  return convolution(hist,N,kernel,rmEdge,name);
}

TH1F*
tools::hannWin(const TH1F* hist, int N,bool rmEdge,const char* name)
{
  return raisedCosWin(hist,N,0.5,0.5,0,0,0,rmEdge,name);
}

TH1F*
tools::hammingWin(const TH1F* hist, int N,bool rmEdge,const char* name)
{
  return raisedCosWin(hist,N,0.54,0.46,0,0,0,name);
}

TH1F*
tools::blackmanWin(const TH1F* hist, int N, double alpha,bool rmEdge,const char* name)
{
  return raisedCosWin(hist,N,0.5*(1-alpha),0.5,0.5*alpha,0,0,rmEdge,name);
}


TH1F*
tools::blackmanHarrisWin(const TH1F* hist, int N, bool rmEdge,const char* name)
{
  return raisedCosWin(hist,N,0.35875,0.48829,0.14128,0.01168,0,name);
}

TH1F*
tools::blackmanNutallWin(const TH1F* hist, int N,bool rmEdge,const char* name)
{
  return raisedCosWin(hist,N,0.3635819,0.4891775,0.1365995,0.0106511,0,rmEdge,name);
}

TH1F*
tools::flatTopWin(const TH1F* hist, int N, bool rmEdge,const char* name)
{
  return raisedCosWin(hist,N,1,1.93,1.29,0.388,0.032,rmEdge,name);
}

TH1F*
tools::nutallWin(const TH1F* hist, int N,bool rmEdge, const char* name)
{
  return raisedCosWin(hist,N,0.355768,0.487396,0.144232,0.012604,0,rmEdge,name);
}

TH1F*
tools::bohmanWin(const TH1F* hist, int N,bool rmEdge, const char* name)
{


  //Set the convolution kernel
  double kernel[1+2*N];
  for (int i = -N; i<=N; i++){
      kernel[i+N] =
            (1-abs(i)/N) * cos(PI*abs(i)/N)+ 1./PI * sin(PI*abs(i)/N);
  }

  return convolution(hist,N,kernel,rmEdge,name);
}

TH1F*
tools::kaiserWin(const TH1F* hist, int N,double a, bool rmEdge, const char* name)
{
  //Default alpha: 3
  //Set the convolution kernel
  double kernel[1+2*N];
  for (int i = -N; i<=N; i++){
      kernel[i+N] =
	TMath::BesselI0(PI*a*sqrt(1-i*i*1.0/(N*N) ));
  }

  return convolution(hist,N,kernel,rmEdge,name);

}

TH1F*
tools::lanczosWin(const TH1F* hist, int N,bool rmEdge,const char* name)
{

  //Set the convolution kernel
  double kernel[1+2*N];
  for (int i = -N; i<=N; i++){
      kernel[i+N] =
           i==0?1 :
           sin(i*PI / N) *N / (i*PI);
  }
  return convolution(hist,N,kernel,rmEdge, name);
}

TH1F*
tools::lanczosFilter(const TH1F* hist, int N,bool rmEdge,const char* name)
{

  //Set the convolution kernel
  double kernel[1+2*N];
  for (int i = -N; i<=N; i++){
      kernel[i+N] =
        i==0?1 :
	sin(i*PI / N)*sin(i*PI) *N / (i*i*PI*PI);
  }
  return convolution(hist,N,kernel,rmEdge, name);
}

TH1F*
tools::cosineWin(const TH1F* hist, int N,bool rmEdge,const char* name)
{
  //Set the convolution kernel
  double kernel[1+2*N];
  for (int i = -N; i<=N; i++){
      kernel[i+N] =
	cos(PI*i/(2*N));
  }

  return convolution(hist, N, kernel, rmEdge,name);
}

TH1F* 
tools::rectWin(const TH1F* hist, int N, bool rmEdge,const char* name)
{

  double kernel[2*N+1];
  for (int i = 0; i<2*N+1; i++) kernel[i] = 1./(2*N+1);

  return convolution(hist,N,kernel,rmEdge,name);
}

TH1F*
tools::triangleWin(const TH1F* hist, int N, bool rmEdge,const char* name)
{
  double kernel[2*N+1];
  for (int i = -N; i<=N; i++) kernel[i+N] = 1 - abs(i)/N;

  return convolution(hist,N,kernel,rmEdge,name);
}

TH1F*
tools::gausWin(const TH1F* hist, int N,double sigma, bool rmEdge,const char* name)
{
  double kernel[2*N+1];
  for (int i = -N; i<=N; i++) kernel[i+N] = exp(-i*i / (2*sigma*sigma) );

  return convolution(hist,N,kernel,rmEdge,name);
}

TH1F*
tools::expoWin(const TH1F* hist, int N,double x0, bool rmEdge,const char* name)
{
  double kernel[2*N+1];
  for (int i = -N; i<=N; i++) kernel[i+N] = exp(-abs(i)/x0 );

  return convolution(hist,N,kernel,rmEdge,name);
}


/*
==========================================
True convolution based smoothing
(use integrals of the given functions to determine weights)
==========================================
*/
TH1F* 
tools::gausConv(const TH1F* hist, double sigma, bool rmEdge,double cut,const char* name)
{

  //Set the convolution kernel
  int N = (int) (sigma*sqrt(-2*log(cut)));
  double kernel[1+2*N];
  for (int i = -N; i<=N; i++){
      kernel[i+N] = 
	erf((i+0.5)/(sqrt(2)*sigma)) - erf((i-0.5)/(sqrt(2)*sigma));

  }  
  
  return convolution(hist,N,kernel,rmEdge,name);
}

TH1F* 
tools::triangleConv(const TH1F* hist, int N,bool rmEdge, const char* name)
{


  //Set the convolution kernel
  double kernel[1+2*N];
  for (int i = -N; i<=N; i++) kernel[i+N] =
        i==0?1-0.25/(N+0.5) : 1 - fabs(i)/(N+0.5);

  return convolution(hist,N,kernel,rmEdge,name);
}

TH1F* 
tools::expoConv(const TH1F* hist, double x0, bool rmEdge, double cut,const char* name)
{
  //Set the convolution kernel
  int N = (int) (-x0*log(cut));
  double kernel[1+2*N];
  for (int i = -N; i<=N; i++) kernel[i+N] =
        i==0? 2 * (1-exp(-0.5 / x0)) :
        exp(-(abs(i)-0.5)/x0) - exp(-(abs(i)+0.5)/x0);
  
  
  return convolution(hist,N,kernel,rmEdge,name);

}

TH1F* 
tools::fermiConv(const TH1F* hist, double x0, double t,bool rmEdge, double cut, const char* name)
{
  //Set the convolution kernel
  int N = (int) ceil(x0+t*log((1+exp(-x0/t)) / cut - 1));
  double kernel[1+2*N];
  for (int i = -N; i<=N; i++)
  {
    double y1=i==0?exp(-x0/t)+1:exp((abs(i)-0.5-x0)/t)+1;
    double y2=exp((abs(i)+0.5-x0)/t)+1;
    kernel[i+N] = i==0? 2 * log((y2-1)*y1 /( (y1-1)*y2 )):
           log( (y2-1)*y1/ ((y1-1)*y2));

  }
  

  return convolution(hist,N,kernel,rmEdge,name);

}


/*
==========================================

Other smoothing algorithms
==========================================
*/


TH1F* 
tools::savitzkyGolay(const TH1F* hist, int order, int size,bool rmEdge, const char* name)
{

  TMatrixD A(2*size+1,order+1);
  for (int i = -size; i <= size ; i++)
  {
    for (int j = 0; j<= order; j++)
    {
      A[i+size][j] = j==0 ? 1:pow(1.0*i,1.0*j);
    }
  }
  TMatrixD At(order+1,2*size+1);
  At.Transpose(A);
  TMatrixD AinvAt = (At*A).Invert()*At;

  double kernel[2*size+1];
  for (int i = -size; i<=size; i++)
    kernel[i+size] = AinvAt[0][i+size];

  return convNoNorm(hist,size,kernel,rmEdge,name);
}

TH1F* 
tools::medianFilter(const TH1F* hist, int size, const char* name)
{

  //Create the filtered histogram
  int nbins = hist->GetNbinsX();
  double xmax = hist->GetXaxis()->GetXmax();
  double xmin = hist->GetXaxis()->GetXmin();
  double xwidth = hist->GetXaxis()->GetBinWidth(1);
  TH1F* filter = new TH1F(name,"Median Filter",nbins-2*(size+1),xmin+(size+1)*xwidth,xmax-(size+1)*xwidth);

  vector<double> v(2*size+1);


  for (int i = 1; i<=2*size+1; i++){
    v[i-1] = get(hist,i);
  }

 
  sort(v.begin(),v.end());
  double value;
  vector<double>::iterator iter;
  for (int i = size+2; i<= nbins-size-1; i++)
  {


    set(filter,i-size,v[size]);

    //Find and erase first element

    iter = lower_bound(v.begin(),v.end(),get(hist,i-size));

    v.erase(iter);
    //Add new element in its place
    if (i==nbins-size)break;
 
    value = get(hist,i+size+1);

    iter = lower_bound(v.begin(),v.end(),value);

    if (iter-v.end()==0) v.push_back(value);
    else v.insert(iter,value);

  }//loop over all non-edge bins 

  return filter;

}

/*
======================================================
Other useful functions and transforms
(Integrals, derivatives, etc)
Plan to add FFTW3 support for FFT, reverse FFT and autocorrelation
====================================================
*/
TH1F* 
tools::integ(const TH1F* hist, bool useXwidth, const char* name)
  {


    int nbins = hist->GetNbinsX();
    double xmax = hist->GetXaxis()->GetXmax();
    double xmin = hist->GetXaxis()->GetXmin();

    TH1F* integral = new TH1F(name,"Integral",nbins, xmin, xmax);
 
    if (useXwidth){
      double xwidth = hist->GetXaxis()->GetBinWidth(1);

      set(integral,1,0.5*get(hist,1))*xwidth;
      for (int i = 2; i<=nbins; i++)
        set(integral,i, get(integral,i-1)+0.5 * xwidth *  
           (get(hist,i-1)+get(hist,i))  
           );

    } else{

      set(integral,1,0.5*get(hist,1));
      for (int i = 2; i<=nbins; i++)
        set(integral,i, get(integral,i-1)+0.5 * 
           (get(hist,i-1)+get(hist,i))  
           );

    }
   return integral;
  }

double factorial(int n)
{
  return n<2?1:factorial(n-1)*n;

}

TH1F* 
tools::sgDeriv(const TH1F* hist, int order, int size, int nder, bool rmEdge, const char* name)
{
  double derFact = factorial(nder);

  TMatrixD A(2*size+1,order+1);
  for (int i = -size; i <= size ; i++)
  {
    for (int j = 0; j<= order; j++)
    {
      A[i+size][j] = j==0 ? 1:pow(1.0*i,1.0*j);
    }
  }
  TMatrixD At(order+1,2*size+1);
  At.Transpose(A);
  TMatrixD AinvAt = (At*A).Invert()*At;
  double kernel[2*size+1];

  for (int i = -size; i<=size; i++)
    kernel[i+size] = AinvAt[nder][i+size]*derFact;

  return convNoNorm(hist,size,kernel,rmEdge,name);
}


//================================= FIR and IIR filters ===================================//


TH1 *
tools::FIRFilter(const TH1*hist, unsigned int order, const double * bcoeffs)
{
  if (order < 1) return 0; 

  TH1 * out = (TH1*) hist->Clone(hist->GetName() + TString("_FIR")); 
  out->SetTitle(out->GetName()); 
  for (int i = 0; i < out->GetNbinsX(); i++)
  {
    double ynew = 0; 
    for (unsigned int j =0; j <= order; j++)
    {
      if (i <(int) j) break; 
      ynew +=  bcoeffs[j] * hist->GetBinContent(i+1-j); 
    }

    out->SetBinContent(i+1,ynew); 
  }

  return out; 

}

TH1 *
tools::IIRFilter(const TH1* hist, unsigned int order, const double * acoeffs, const double * bcoeffs)
{
  if (order < 1) return 0; 

  TH1 * out = (TH1*) hist->Clone(hist->GetName() + TString("_IIR")); 
  out->SetTitle(out->GetName()); 
  for (int i = 0; i < out->GetNbinsX(); i++)
  {
    double ynew = 0; 

    for (unsigned int j =0; j <= order; j++)
    {
      if (i < (int)j) break; 
      ynew +=  bcoeffs[j] * hist->GetBinContent(i+1-j); 
      if (j>0 && acoeffs[j] != 0) ynew-= acoeffs[j] * out->GetBinContent(i+1-j); 
    }

    ynew /= acoeffs[0]; 

    out->SetBinContent(i+1,ynew); 
  }

  return out; 
}


//Sorta equivalent to matlab poly. (except indices are reversed so that the index matches the order)
static std::complex<double> * poly(unsigned int n, std::complex<double> * zeros)
{

  for (unsigned int i = 0; i < n; i++)
  {
 //   std::cout << zeros[i] << " " ; 
  }

  std::complex<double> * pol = new std::complex<double>[n+1]; 

  //Start out with (x - zero_0)
  pol[0] = -zeros[0]; 
  pol[1] = std::complex<double>(1,0); 


  for (unsigned int i = 1; i < n; i++)
  {

    //Multiply  by each (x - zero_i)
    for (unsigned int j = i+1;;  j--)
    {
      pol[j] = (j > 0 ? pol[j-1] : 0) - (j < i+1 ? pol[j] * zeros[i] : 0); 
      if (j == 0) break; 
    }

  }

  
  return pol; 
}

void 
tools::butterworthBandStop(double w0, double delta, unsigned int order, double ** acoeffs, double ** bcoeffs)
{
  *(acoeffs) = new double[2*order+1] ;
  *(bcoeffs) = new double[2*order+1] ;

  std::complex<double> * poles = new std::complex<double>[order]; 
  std::complex<double> * tmp = new std::complex<double>[order]; 
  std::complex<double> pole_product(1,0); 

  //std::cout << "Here" << std::endl; 
  //This code is based on the octave-forge implementation of butter.m
  // There are certainly simplifications possible, but I wanted something that worked...  


  //Compute edges
  double wlow = w0 - delta; 
  double whigh = w0 + delta; 

  //Convert frequency to s space 
  double Wlow = TMath::Tan(TMath::Pi() * wlow / 2.); 
  double Whigh = TMath::Tan(TMath::Pi() * whigh / 2.); 

  //std::cout << " Wlow: " << Wlow << std::endl; 
  //std::cout << " Whigh: " << Whigh << std::endl; 

  //Compute Butterworth poles for given order
  for (unsigned int i = 0; i < order; i++)
  {
    poles[i] = std::exp(std::complex<double>(0, TMath::Pi()* (2*(i+1)+order - 1)/(2*order))); 
    tmp[i] = (Whigh - Wlow) / poles[i] / 2.; 
    pole_product*= -poles[i]; 
  }

  //std::cout << "initial poles: "; 
//  for (unsigned int i = 0; i <order; i++) std::cout << poles[i] << " "; 
  //std::cout << std::endl; 

  //std::cout << "tmp: "; 
//  for (unsigned int i = 0; i <order; i++) std::cout << tmp[i] << " "; 
  //std::cout << std::endl; 
  //Transform to band stop filter 
  
  double gain = 1./std::real(pole_product); 
 
  std::complex<double> * stop_poles = new std::complex<double>[order*2]; 
  std::complex<double> * stop_zeros = new std::complex<double>[order*2]; 
  for (unsigned int i = 0; i < order; i++)
  {
    stop_poles[i] = tmp[i] + std::sqrt(std::pow(tmp[i],2) - Wlow*Whigh); 
    stop_poles[order+i] = tmp[i] - std::sqrt(std::pow(tmp[i],2) - Wlow*Whigh); 
  }  

  for (unsigned int i = 0; i < order*2; i++)
  {
    stop_zeros[i] =  -std::sqrt(std::complex<double>(-Whigh*Wlow,0)); 
    if (i &1) stop_zeros[i] *= -1; 
  }

  //std::cout << "sz: "; 
  //for (unsigned int i = 0; i < order*2; i++) std::cout << stop_zeros[i] << " "; 
  //std::cout << std::endl; 
  //std::cout << "sp: "; 
//  for (unsigned int i = 0; i < order*2; i++) std::cout << stop_poles[i] << " "; 
  //std::cout << std::endl; 
  //std::cout << "sg: " << gain << std::endl; 



  //Perform bilinear transform to convert to digital 
  
  std::complex<double> * digi_poles = new std::complex<double>[order*2]; 
  std::complex<double> * digi_zeros = new std::complex<double>[order*2]; 

  std::complex<double> digi_gain = std::complex<double>(gain,0);  

  std::complex<double> one(1,0); 
  for (unsigned int i = 0; i < order*2; i++)
  {
    digi_gain *=  (one -stop_zeros[i]) / (one - stop_poles[i]); 

    digi_poles[i] = (one + stop_poles[i]) / (one - stop_poles[i]); 
    digi_zeros[i] = (one + stop_zeros[i]) / (one - stop_zeros[i]); 

  }

  digi_gain = std::real(digi_gain); 

 // std::cout << "z: "; 
//  for (unsigned int i = 0; i < order*2; i++) std::cout << digi_zeros[i] << " "; 
 // std::cout << std::endl; 
 // std::cout << "p: "; 
//  for (unsigned int i = 0; i < order*2; i++) std::cout << digi_poles[i] << " "; 
 // std::cout << std::endl; 
 // std::cout << "g: " << digi_gain << std::endl; 

  //Now, get coefficients from poles and zeroes
  std::complex<double> * bpoly = poly(2*order, digi_zeros); 
  std::complex<double> * apoly = poly(2*order, digi_poles); 

  for (unsigned int i = 0; i < 2*order+ 1; i++)
  {
       (*bcoeffs)[i] = (double) std::real(digi_gain * bpoly[2*order - i]);             
       (*acoeffs)[i] = (double) std::real(apoly[2*order - i]);             
  }

  delete poles; 
  delete digi_poles; 
  delete digi_zeros; 
  delete stop_poles; 
  delete stop_zeros; 
  delete apoly; 
  delete bpoly; 
  delete tmp; 
}

void 
tools::zTransform(const TH1F* hist, double x, double y,double& a, double& b, bool isReIm)
{
    int nbins = hist->GetNbinsX();
    double mag, arg, magn;
    a=0, b=0;
    if (isReIm){
      mag = sqrt(x*x+y*y);
      arg = atan2(y,x);

   } else{
      mag = x;
      arg = y;
   }
   for (int i = 1; i <= nbins; i++)
   {
     magn = pow(mag,i-1);
     a += magn * cos((i-1)*arg);
     b += magn * sin((i-1)*arg);
   }

}

void 
tools::addWaveformVectorToTObjArray(TObjArray * wfvlist, int index, const char * name, const char * type){

  //  cout << "inside addWaveformVectorToTObjArray( "
  //       << wfvlist
  //       << " , "
  //       << index 
  //       << " , "
  //       << name 
  //       << " , "
  //       << type
  //       << " ) ..." 
  //       << endl;

  TString tsType(type);
  
  if( !( tsType.CompareTo("fast",TString::kIgnoreCase) ) ||
      !( tsType.CompareTo("FastWfVector",TString::kIgnoreCase) ) ){

    wfvlist->AddAt( new FastWfVector(0,0,0,name), index);
    
  } else if ( !( tsType.CompareTo("csp",TString::kIgnoreCase) ) ||
	      !( tsType.CompareTo("CspWfVector",TString::kIgnoreCase) ) ){
    
    wfvlist->AddAt( new CspWfVector(0,0,0,name), index);
    
  } else if ( !( tsType.CompareTo("pmt",TString::kIgnoreCase) ) ||
	      !( tsType.CompareTo("PMTWfVector",TString::kIgnoreCase) ) ){
    
    wfvlist->AddAt( new PMTWfVector(0,0,0,name), index);
  
  } else {
    
    wfvlist->AddAt( new WaveformVector(0,0,0,name), index );

  }
}

TH1*
tools::createFirOutputHist(TH1* in, int nkernel)
{
  //Look for class:
  char c = 'F';
  if (in->IsA()->GetBaseClass("TH1D"))
    c = 'D';

  TH1* h = createFirOutputHist(in->GetNbinsX(),nkernel,c);
  if (!h) return 0;
  h->SetXTitle(in->GetXaxis()->GetTitle());
  h->SetYTitle(in->GetYaxis()->GetTitle());

  double binW = in->GetXaxis()->GetBinWidth(1);
  double xmin, xmax;
  xmin = in->GetXaxis()->GetXmin() + binW * (nkernel/2);
  xmax = in->GetXaxis()->GetXmax() - binW * (nkernel/2);
  if (nkernel%2==0)
    xmax += binW;
  h->GetXaxis()->SetLimits(xmin,xmax);
  return h;
}

TH1*
tools::createFirOutputHist(int nin, int nkernel, char type)
{
  if (nin-nkernel<=0) return 0;
  TH1* h;

  switch (type){
    case 'C':
    case 'c':
      h = new TH1C("FirC","",nin-nkernel+1,0,nin-nkernel+1);
      break;
    case 'S':
    case 's':
      h = new TH1S("FirS","",nin-nkernel+1,0,nin-nkernel+1);
      break;
    case 'I':
    case 'i':
      h = new TH1I("FirI","",nin-nkernel+1,0,nin-nkernel+1);
      break;
    case 'F':
    case 'f':
      h = new TH1F("FirF","",nin-nkernel+1,0,nin-nkernel+1);
      break;
    case 'D':
    case 'd':
      h = new TH1D("FirD","",nin-nkernel+1,0,nin-nkernel+1);
      break;
    default:
      return 0;
  }

  return h;

}

TH1*
tools::FirHistFilter(TH1* in, TH1* out, waveform::FirFilter<double>* fir)
{
  if (out==0) out = createFirOutputHist(in,fir->kernel_size());
  if (   in->IsA()->GetBaseClass("TH1F"))
    fir->run(in->GetNbinsX(),( (TH1F*)in )->GetArray()+1);
  else if ( in->IsA()->GetBaseClass("TH1D"))
    fir->run(in->GetNbinsX(),( (TH1D*)in )->GetArray()+1);
  else {
    float x[in->GetNbinsX()];
    for (int i = 1; i <= in->GetNbinsX(); i++)
      x[i-1] = in->GetBinContent(i);
    fir->run(in->GetNbinsX(),x);
  }
  for (int i = 0; i < fir->output_size(); i++)
    out->SetBinContent(i+1,fir->output()[i]);
  return out;
}

TH1*
tools::createIirOutputHist(TH1* in, int nkernel)
{
  //Look for class:
  char c = 'F';
  if (in->IsA()->GetBaseClass("TH1D"))
    c = 'D';

  TH1* h = createFirOutputHist(in->GetNbinsX(),nkernel,c);
  if (!h) return 0;
  h->SetXTitle(in->GetXaxis()->GetTitle());
  h->SetYTitle(in->GetYaxis()->GetTitle());

  double binW = in->GetXaxis()->GetBinWidth(1);
  double xmin, xmax;
  xmin = in->GetXaxis()->GetXmin() + binW * (nkernel/2);
  xmax = in->GetXaxis()->GetXmax() - binW * (nkernel/2);
  if (nkernel%2==0)
    xmax += binW;
  h->GetXaxis()->SetLimits(xmin,xmax);
  return h;
}

TH1*
tools::createIirOutputHist(int nin, int nkernel, char type)
{
  if (nin-nkernel<=0) return 0;
  TH1* h;

  switch (type){
    case 'C':
    case 'c':
      h = new TH1C("IirC","",nin-nkernel+1,0,nin-nkernel+1);
      break;
    case 'S':
    case 's':
      h = new TH1S("IirS","",nin-nkernel+1,0,nin-nkernel+1);
      break;
    case 'I':
    case 'i':
      h = new TH1I("IirI","",nin-nkernel+1,0,nin-nkernel+1);
      break;
    case 'F':
    case 'f':
      h = new TH1F("IirF","",nin-nkernel+1,0,nin-nkernel+1);
      break;
    case 'D':
    case 'd':
      h = new TH1D("IirD","",nin-nkernel+1,0,nin-nkernel+1);
      break;
    default:
      return 0;
  }

  return h;

}

TH1*
tools::IirHistFilter(TH1* in, TH1* out, waveform::IirFilter<double>* iir)
{
  if (!out) out = createIirOutputHist(in,iir->fir_size());
  if (   in->IsA()->GetBaseClass("TH1F"))
    iir->run(in->GetNbinsX(),( (TH1F*)in )->GetArray()+1);
  else if ( in->IsA()->GetBaseClass("TH1D"))
    iir->run(in->GetNbinsX(),( (TH1D*)in )->GetArray()+1);
  else {
    float x[in->GetNbinsX()];
    for (int i = 1; i <= in->GetNbinsX(); i++)
      x[i-1] = in->GetBinContent(i);
    iir->run(in->GetNbinsX(),x);
  }

  for (int i = 0; i < iir->output_size(); i++)
    out->SetBinContent(i+1,iir->output()[i]);
  return out;
}

unsigned
tools::zeroSuppress(TH1* in, double zero_level,  double noise_ceiling, unsigned min_samples, unsigned window_size, double * mean, double * rms)
{
  vector<int> kill_list; 
  
  //this is very slow way of doing this... probably using a buffer would be way faster... 
/*  for (unsigned i = 1; i <= in->GetNbinsX(); i++)
  {
    int nabove = 0; 
    for (unsigned j = i-window_size; j <= i+window_size; j++)
    {
      if (j < 1) continue; 
      if (j > in->GetNbinsX()) break; 
      if (fabs(in->GetBinContent(j) - zero_level) >  noise_ceiling)
      {
        nabove++; 
      }
    }
    if (nabove < min_samples)
    {
      kill_list.push_back(i); 
    }
  }
 */

   int nabove = 0; 
   for (unsigned i = 1; i <= in->GetNbinsX() + window_size; i++)
   {
     if (i <= in->GetNbinsX() && fabs(in->GetBinContent(i) - zero_level) > noise_ceiling) 
     {
        nabove++; 
     }

     if (i > window_size)
     {
       if (nabove < min_samples)
       {
          kill_list.push_back(i-window_size-1); 
       }

       if (i > 2 * window_size + 1 && fabs(in->GetBinContent(i-2*window_size-1) - zero_level) > noise_ceiling)
       {
         nabove--; 
       }
     }
   }


  double sum = 0; 
  double sum2 = 0; 

  for (unsigned i = 0; i < kill_list.size(); i++)
  {
    if (mean || rms)
    {
      sum+= in->GetBinContent(kill_list[i]);  
    }
    if (rms)
    {
      sum2+= in->GetBinContent(kill_list[i]) * in->GetBinContent(kill_list[i]);  
    }

    in->SetBinContent(kill_list[i],zero_level); 
  }

  if (mean)
  {
    *mean = sum / kill_list.size(); 
  }

  if (rms)
  {
    *rms = sqrt(sum2/kill_list.size() - sum*sum/kill_list.size()/kill_list.size());  
  }

  return kill_list.size(); 
}

TH1 * 
tools::cropZeros(const TH1 * in, TH1 * outptr)
{
  int start_bin = 1, end_bin = in->GetNbinsX();  

  for (int i = 1; i <= in->GetNbinsX(); i++)
  {
    if (in->GetBinContent(i) == 0)
      start_bin++; 
    else
      break; 
  }

  if (start_bin == in->GetNbinsX())
  {
    //empty hist! 
    return 0; 
  }

  for (int i = in->GetNbinsX(); i >=1; i--)
  {
    if (in->GetBinContent(i) == 0)
      end_bin--; 
    else
      break; 
  }

  TH1 * out = DmtpcRootTools::newTH1StealType(in, TString(in->GetName())+TString("_cropped"), TString(in->GetTitle()) +TString( " (cropped)"), 
                                              end_bin - start_bin + 1, in->GetXaxis()->GetBinLowEdge(start_bin), in->GetXaxis()->GetBinLowEdge(end_bin+1),outptr); 
  int j = 1; 
  for (int i = start_bin; i <= end_bin; i++)
  {
    out->SetBinContent(j++, in->GetBinContent(i)); 
  }

  return out; 
}


TH1 * 
tools::resample(const TH1* in, double factor, TH1 * outptr)
{

#ifndef HAVE_LIBSAMPLERATE
  std::cerr<< "Error in waveform::tools::resample():   Not compiled with libsamplerate support. " << std::endl; 
  return NULL; 
#else

  unsigned nin = in->GetNbinsX(); 
  unsigned nout = nin / factor; 
  TH1 * out = DmtpcRootTools::newTH1StealType(in, TString(in->GetName())+TString("_resampled"), TString(in->GetTitle()) +TString( " (resampled)"), 
              nout, in->GetXaxis()->GetXmin(), in->GetXaxis()->GetXmax(),outptr); 

  if (factor == 1)
  {
    
    for (int i = 0; i <nin; i++)
    {
      out->SetBinContent(i+1, in->GetBinContent(i+1)); 
    }
    return out; 
  }

  float *float_data_in = new float[nin];
  float *float_data_out = new float[nout];

  //preserve under/overflow
  out->SetBinContent(0, in->GetBinContent(0)); 
  out->SetBinContent(nout+1, in->GetBinContent(nin+1)); 

  //what we do here depends on the type

  float max_size = out->IsA()->InheritsFrom("TArrayC") ? 1<<8 : 
                   out->IsA()->InheritsFrom("TArrayS") ? 1<<16 : 
                   out->IsA()->InheritsFrom("TArrayI") ? INT_MAX : 
                   FLT_MAX; 
      

  

  for (int i = 0; i <nin; i++)
  {
    float_data_in[i] = (float) in->GetBinContent(i+1) / max_size; 
  }

  SRC_DATA data; 
  data.data_in = float_data_in; 
  data.data_out = float_data_out; 
  data.input_frames = nin; 
  data.output_frames = nout; 
  data.src_ratio = 1./factor; 

  if (int err = src_simple(&data,0,1)) 
  {
    std::cerr <<  "resample error: " << src_strerror(err) << std::endl; 
  }

  for (int i = 0; i < nout; i++)
  {
    out->SetBinContent(i+1, float_data_out[i] * max_size);  
  }


  delete float_data_out;
  delete float_data_in;

  return out; 
#endif


}
