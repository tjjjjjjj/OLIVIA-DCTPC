#ifndef KERNELS_CC
#define KERNELS_CC
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMath.h"
#include <vector>
#include <cmath>
using namespace std;
#define PI 3.14159026535898

/*
==================================
Basic window convolutions

===================================
*/
vector<double>
raisedCosWin(int N,double a0=1, double a1=0, double a2=0, double a3=0, double a4=0)
{
  //Set the convolution kernel
  vector<double> kernel(1+2*N);
  for (int i = -N; i<=N; i++){
      kernel[i+N] =
       a0 + a1 * cos(PI*i/N) + a2 * cos(2*PI*i/N)
       + a3*cos(3*PI*i/N) + a4*cos(4*PI*i/N);
  }

  return kernel;
}

vector<double>
hannWin(int N)
{
  return raisedCosWin(N,0.5,0.5,0,0,0);
}

vector<double>
hammingWin(int N)
{
  return raisedCosWin(N,0.54,0.46,0,0,0);
}

vector<double>
blackmanWin(int N, double alpha)
{
  return raisedCosWin(N,0.5*(1-alpha),0.5,0.5*alpha,0,0);
}


vector<double>
blackmanHarrisWin(int N)
{
  return raisedCosWin(N,0.35875,0.48829,0.14128,0.01168,0);
}

vector<double>
blackmanNutallWin(int N)
{
  return raisedCosWin(N,0.3635819,0.4891775,0.1365995,0.0106511,0);
}

vector<double>
flatTopWin(int N)
{
  return raisedCosWin(N,1,1.93,1.29,0.388,0.032);
}

vector<double>
nutallWin(int N)
{
  return raisedCosWin(N,0.355768,0.487396,0.144232,0.012604,0);
}

vector<double>
bohmanWin(int N)
{

  //Set the convolution kernel
  vector<double> kernel(1+2*N);
  for (int i = -N; i<=N; i++){
      kernel[i+N] =
            (1-fabs(i)/N) * cos(PI*fabs(i)/N)+ 1./PI * sin(PI*fabs(i)/N);
  }

  return kernel;
}

vector<double>
kaiserWin(int N,double a=3)
{
  //Default alpha: 3
  //Set the convolution kernel
  vector<double> kernel(1+2*N);
  for (int i = -N; i<=N; i++){
      kernel[i+N] =
	TMath::BesselI0(PI*a*sqrt(1-i*i*1.0/(N*N) ));
  }

  return kernel;

}

vector<double>
lanczosWin(int N)
{

  //Set the convolution kernel
  vector<double> kernel(1+2*N);
  for (int i = -N; i<=N; i++){
      kernel[i+N] =
           i==0?1 :
           sin(i*PI / N) *N / (i*PI);
  }
  return kernel;
}


vector<double>
cosineWin(int N)
{
  //Set the convolution kernel
  vector<double> kernel(1+2*N);
  for (int i = -N; i<=N; i++){
      kernel[i+N] =
	cos(PI*i/(2*N));
  }

  return kernel;
}

vector<double>
rectWin(int N)
{

  vector<double> kernel(2*N+1);
  for (int i = 0; i<2*N+1; i++) kernel[i] = 1./(2*N+1);

  return kernel;
}

vector<double>
triangleWin(int N)
{
  vector<double> kernel(2*N+1);
  for (int i = -N; i<=N; i++) kernel[i+N] = 1 - fabs(i)/N;

  return kernel;
}

vector<double>
gausWin(int N,double sigma)
{
  vector<double> kernel(2*N+1);
  for (int i = -N; i<=N; i++) kernel[i+N] = exp(-i*i / (2*sigma*sigma) );

  return kernel;
}

vector<double>
expoWin(int N,double x0)
{
  vector<double> kernel(2*N+1);
  for (int i = -N; i<=N; i++) kernel[i+N] = exp(-abs(i)/x0 );

  return kernel;
}


/*
==========================================
True convolution based smoothing
(use integrals of the given functions to determine weights)
==========================================
*/
vector<double>
gausConv(double sigma,double cut)
{

  //Set the convolution kernel
  int N = (int) (sigma*sqrt(-2*log(cut)));
  vector<double> kernel(1+2*N);
  for (int i = -N; i<=N; i++){
      kernel[i+N] = 
	erf((i+0.5)/(sqrt(2)*sigma)) - erf((i-0.5)/(sqrt(2)*sigma));

  }  
  
  return kernel;
}

vector<double>
triangleConv(int N)
{


  //Set the convolution kernel
  vector<double> kernel(1+2*N);
  for (int i = -N; i<=N; i++) kernel[i+N] =
        i==0?1-0.25/(N+0.5) : 1 - fabs(i)/(N+0.5);

  return kernel;
}

vector<double>
expoConv(double x0, double cut)
{
  //Set the convolution kernel
  int N = (int) (-x0*log(cut));
  vector<double> kernel(1+2*N);
  for (int i = -N; i<=N; i++) kernel[i+N] =
        i==0? 2 * (1-exp(-0.5 / x0)) :
        exp(-(abs(i)-0.5)/x0) - exp(-(abs(i)+0.5)/x0);
  
  
  return kernel;

}

vector<double>
fermiConv(double x0, double t,double cut)
{
  //Set the convolution kernel
  int N = (int) ceil(x0+t*log((1+exp(-x0/t)) / cut - 1));
  vector<double> kernel(1+2*N);
  for (int i = -N; i<=N; i++)
  {
    double y1=i==0?exp(-x0/t)+1:exp((abs(i)-0.5-x0)/t)+1;
    double y2=exp((abs(i)+0.5-x0)/t)+1;
    kernel[i+N] = i==0? 2 * log((y2-1)*y1 /( (y1-1)*y2 )):
           log( (y2-1)*y1/ ((y1-1)*y2));

  }
  

  return kernel;

}

vector<double>
savitzkyGolay(int order, int size,int nder)
{
  double derFact = TMath::Factorial(nder);
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

  vector<double> kernel(2*size+1);
  for (int i = -size; i<=size; i++)
    kernel[i+size] = AinvAt[nder][i+size]*derFact;

  return kernel;
}
#endif
