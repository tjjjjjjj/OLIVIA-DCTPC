#include "DmtpcMath.hh"
#include <cmath>
#include "TRandom.h"
#include <iostream>
#include "TMath.h" 
#include "TMatrixD.h" 
#include "TVectorD.h" 
#include "TDecompSVD.h" 
#include "TDecompChol.h" 
#include <gsl/gsl_poly.h>
#include <set>
#include <gsl/gsl_math.h>


double DmtpcMath::Beta(double x,double y)
{
  double lB = lgamma(x)+lgamma(y)-lgamma(x+y);
  return exp(lB);
}
double DmtpcMath::Gamma(double x){return tgamma(x);}
double DmtpcMath::logGamma(double x){return lgamma(x);}

double
DmtpcMath::betaPdf(double x,double par[])
{
  //Bounds check
  if (x<=0||x>=1) return 0;
  //Parameters
  double& alpha = par[0];
  double& beta = par[1];

  //Log of output value:
  double lval = lgamma(alpha+beta)-lgamma(alpha)-lgamma(beta)+(alpha-1)*log(x) + (beta-1)*log(1-x);
  return exp(lval);

}

void
DmtpcMath::betaPdfStats(double x[],double par[])
{
  double& alpha = par[0];
  double& beta = par[1];

  x[MEAN] = alpha/(alpha+beta);
  x[MED] = -NAN;
  if (alpha>1&&beta>1)
    x[MODE] = (alpha-1)/(alpha+beta-2);
  else x[MODE] = NAN;
  x[VAR] =  (alpha*beta)/(alpha+beta)/(alpha+beta)/(alpha+beta+1);
  x[SKEW] = 2*(beta-alpha)*sqrt(alpha+beta+1)/sqrt(alpha*beta)/(alpha+beta+2);
  x[KURT] = 6*( (alpha-beta)*(alpha-beta)*(alpha+beta+1)-alpha*beta*(alpha+beta+2) ) /
            (alpha*beta*(alpha+beta+2)*(alpha+beta+3));
}


double 
DmtpcMath::vonMisesDist(double *x, double *par)
{
  double A = par[0]; 
  double u = par[1]; 
  double k = fabs(par[2]); 

  return A*exp(k * cos(*x-u)) / (2*M_PI * TMath::BesselI0(k));  
}


double 
DmtpcMath::vonMisesDistHT(double *x, double *par)
{
  double f = par[0]; 
  double par2[3]; 
  par2[0] = par[1];
  par2[1] = par[2] + M_PI;
  par2[2] = par[3]; 
  return vonMisesDist(x,par+1) * f + (1-f) * vonMisesDist(x,par2); 
}


double 
DmtpcMath::wrappedNormalDist(double *x, double *par)
{

  double phi = normalizeAngle(*x); 
  double A = par[0]; 
  double u = par[1]; 
  double w = par[2]; 

  int max_k = 3+int(ceil(2 * w / (2*M_PI))); 

  double ans = 0; 
  for (int k = -max_k; k <= max_k; k++)
  {
    ans += exp(-pow(phi - u + 2 * M_PI * k,2)/(2*w*w)); 
  }

  ans *= A/(w * sqrt(2*M_PI)); 

  return ans; 
}


double 
DmtpcMath::wrappedNormalDistHT(double *x, double *par)
{
  double f = par[0]; 
  double par2[3]; 
  par2[0] = par[1];
  par2[1] = par[2] + M_PI;
  par2[2] = par[2]; 
  return wrappedNormalDist(x, par+1) * f + (1-f) * wrappedNormalDist(x,par2);
}


double
DmtpcMath::chiPdf(double x,double par[])
{
  if (x<0) return 0;
  double& ndf = par[0];

  double lval = (1-0.5*ndf)*LOG2+(ndf-1)*log(x)-0.5*x*x-lgamma(0.5*ndf);
  return exp(lval);
}

void
DmtpcMath::chiPdfStats(double x[],double par[])
{
  double& k = par[0];

  x[MEAN] = SQRT2*exp( lgamma(0.5*k+0.5)-lgamma(0.5*k) );
  x[MED] = -NAN;
  x[MODE] = k>=1? sqrt(k-1) : 0;
  x[VAR] = k - x[MEAN]*x[MEAN];
  x[SKEW] = x[MEAN]/pow(x[VAR],1.5)*(1-2*x[VAR]);
  x[KURT] = 2./x[VAR] * (1-x[MEAN]*sqrt(x[VAR])*x[SKEW] - x[VAR]);
}

double
DmtpcMath::chiSquarePdf(double x,double par[])
{
  if (x<=0) return 0;
  double& ndf = par[0];
 
  double lval = - 0.5*ndf*LOG2-lgamma(0.5*ndf)+(0.5*ndf-1)*log(x)-0.5*x;

  return exp(lval);
}

void
DmtpcMath::chiSquarePdfStats(double x[], double par[])
{
  double& ndf = par[0];

  x[MEAN] = ndf;
  x[MED] = ndf*pow(1-2./(9*ndf),3);
  x[MODE] = ndf>2 ? ndf-2 : 0;
  x[VAR] = 2*ndf;
  x[SKEW] = sqrt(8/ndf);
  x[KURT] = 12./ndf; 

}

double
DmtpcMath::crystalBallPdf(double *xx,double par[])
{
  double x = *xx; 
  double y = x - par[0];
  
  double sigma = par[1];
  double n = par[2];
  double alpha = par[3];
  double norm = par[4]; 

  double A = pow(n/fabs(alpha),n)*exp(-0.5*alpha*alpha);
  double B = n/fabs(alpha) - fabs(alpha); 

  if (y > -sigma*alpha)
    return norm * exp( -y*y/(2*sigma*sigma) );
  return norm *A*pow(B-y/sigma,-n);

}


double
DmtpcMath::crystalBall(double y[], double par[])
{
  double x = y[0];
  double xave = par[0];
  double sigma = par[1];
  double n = par[2];
  double alpha = par[3];
  double norm = par[4];
  double A = pow(n/fabs(alpha),n)*exp(-0.5*alpha*alpha);
  double B = n/fabs(alpha) - fabs(alpha); 

  if (x-xave > -sigma*alpha)
    return norm * exp(-(x-xave)*(x-xave)/(2*sigma*sigma));
  return norm*A*pow(B-(x-xave)/sigma,-n);

}


double
DmtpcMath::gammaPdf(double x,double par[])
{
  if (x<=0) return 0;

  double& k = par[0]; 
  double& theta = par[1];

  double lval = -lgamma(k)-k*log(theta)+(k-1)*log(x)-x/theta;
  return exp(lval);
}

void
DmtpcMath::gammaPdfStats(double x[],double par[])
{
  double& k = par[0];
  double& theta = par[1];

  x[MEAN] = k*theta;
  x[MED] = -NAN;
  x[MODE] = k>=1?(k-1)*theta : 0;
  x[VAR] = k*theta*theta;
  x[SKEW] = 2./sqrt(k);
  x[KURT] = 6./k;
}

double
DmtpcMath::inverseGausPdf(double x,double par[])
{
  if (x<=0) return 0;
  double& lambda = par[0];
  double& mu = par[1];

  return sqrt( lambda/(2*PI*x*x*x) )  
         * exp( -lambda*(x-mu)*(x-mu)/(2*mu*mu*x)  );

}

double 
DmtpcMath::inverseGausCdf(double x,double par[])
{
  if (x<=0) return 0;
  double p0[2] = {0,1};
  double& l = par[0];
  double& mu = par[1];
  return normalCdf(sqrt(l/x)*(x/mu-1),p0)
       + exp(2*l/mu) * 
         normalCdf(-sqrt(l/x)*(x/mu+1),p0);

}

void
DmtpcMath::inverseGausPdfStats(double x[],double par[])
{
  double& lambda = par[0];
  double& mu = par[1];
  x[MEAN] = mu;
  x[MED] = -NAN;
  x[MODE] = mu*(sqrt(1+2.25 * mu*mu/(lambda*lambda))-1.5*mu/lambda );
  x[VAR] = mu*mu*mu/lambda;
  x[SKEW] = 3*sqrt(mu/lambda);
  x[KURT] = 15*mu/lambda;
}

double
DmtpcMath::logisticPdf(double x,double par[])
{
  double& mu = par[0];
  double& s = par[1];

  return exp( -(x-mu)/s  ) / s * pow(1+exp(-(x-mu)/s),-2);

}

void
DmtpcMath::logisticPdfStats(double x[],double par[])
{
  double& mu = par[0];
  double& s = par[1];
  x[MEAN] = mu;
  x[MED] = mu; 
  x[MODE] = mu;
  x[VAR] = PI*PI/3 * s*s;
  x[SKEW] = 0;
  x[KURT] = 6./5;

}

double
DmtpcMath::logNormalPdf(double x,double par[])
{
  if (x<=0) return 0;

  double& mean = par[0];
  double& sigma = par[1];

  return 1. / (x*sqrt(2*PI*sigma*sigma))*exp(-(log(x)-mean)*(log(x)-mean)/(2*sigma*sigma));
}

double
DmtpcMath::logNormalCdf(double x,double par[])
{
  if (x<=0) return 0;
  return normalCdf(log(x),par);
}

void
DmtpcMath::logNormalPdfStats(double x[],double par[])
{
  double& mean = par[0];
  double& sigma = par[1];
  double ess = exp(sigma*sigma);
  x[MEAN] = exp(mean) * sqrt(ess);
  x[MED] = exp(mean);
  x[MODE] = exp(mean)/ess;
  x[VAR] = ( ess-1 )*exp(2*mean)*ess;
  x[SKEW] = (ess+2)*sqrt(ess-1);
  x[KURT] = ess*ess*ess*ess+2*ess*ess*ess+3*ess*ess-6;

}

double
DmtpcMath::maxwellPdf(double x,double par[])
{
  return SQRT2/sqrt(PI) * x*x/(par[0]*par[0]*par[0])*exp( -x*x/(2*par[0]*par[0]));
}

void
DmtpcMath::maxwellPdfStats(double x[],double par[])
{
  double a = par[0];
  x[MEAN] = 2*a*SQRT2/sqrt(PI);
  x[MED] = -NAN;
  x[MODE] = sqrt(2)*a;
  x[VAR] = a*a/PI * (3*PI-8);
  x[SKEW] = 2*SQRT2*(16-5*PI)/pow(3*PI-8,1.5);
  x[KURT] = 4*(-96+40*PI-3*PI*PI)/pow(3*PI-8,2);
}

double 
DmtpcMath::normalPdf(double x,double par[])
{
  double& mean = par[0];
  double& sigma = par[1];

  return 1./sqrt(2*PI*sigma*sigma) * exp(-(x-mean)*(x-mean)/(2*sigma*sigma)   );
}

double
DmtpcMath::normalCdf(double x,double par[])
{
  return 0.5*(1+erf((x-par[0])/par[1]/SQRT2) );
}

double
DmtpcMath::gausPdf(double x,double par[]){return normalPdf(x,par);}

void
DmtpcMath::normalPdfStats(double x[],double par[])
{
  x[MEAN] = par[0];
  x[MED] = x[MEAN];
  x[MODE] = x[MEAN];
  x[VAR] = par[1]*par[1];
  x[SKEW] = 0;
  x[KURT] = 0;

}

double
DmtpcMath::rayleighPdf(double x,double par[])
{
  return x/(par[0]*par[0])*exp( -x*x/(2*par[0]*par[0])  );
}

void
DmtpcMath::rayleighPdfStats(double x[],double par[])
{
  double& s = par[0];
  x[MEAN] = s*sqrt(PI)*SQRT1_2;
  x[MED] = s*SQRT2*sqrt(LOG2);
  x[MODE] = s;
  x[VAR] = (2-PI_2)*s*s;
  x[SKEW] = 2*sqrt(PI)*(PI-3)/pow(4-PI,1.5);
  x[KURT] = - (6*PI*PI-24*PI+16)/pow(4-PI,2);
}  

double
DmtpcMath::studentPdf(double x,double par[])
{
  
  double& ndf = par[0];

  double lval = -0.5*log(ndf*PI)+lgamma(0.5*ndf+0.5)-lgamma(0.5*ndf) -(0.5*ndf+0.5)*log(1+x*x/ndf);
  return exp(lval);
}

void
DmtpcMath::studentPdfStats(double x[],double par[])
{
  double& ndf = par[0];
  x[MEAN] = 0;
  x[MODE] = 0;
  x[MED] = 0;
  if (ndf>2) x[VAR] = ndf/(ndf-2);
  else if (ndf>1) x[VAR] = INFINITY;
  else x[VAR] = NAN;
  
  if (ndf>3) x[SKEW] = 0;
  else x[SKEW] = NAN;

  if (ndf>4) x[KURT] = 6./(ndf-4);
  else if (ndf>2) x[KURT] = INFINITY;
  else x[KURT] = NAN;

}


int DmtpcMath::lineLineIntersection(double x1,double y1, double x2, double y2, double x3, double y3, double x4, double y4,double * xint, double * yint)
{
    double denom = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4); 

    if (denom == 0) return 1; 

    *xint = ((x1*y2-y1*x2)*(x3-x4) - (x1-x2)*(x3*y4-y3*x4))/denom; 
    *yint = ((x1*y2-y1*x2)*(y3-y4) - (y1-y2)*(x3*y4-y3*x4))/denom; 

    return 0; 
}



bool  DmtpcMath::lineRectIntersectionPoints(double line_th, double line_r, double rect_xmin, double rect_xmax, double rect_ymin, double rect_ymax, double * x0, double * y0, double * x1, double * y1)
{
  double sinth = sin(line_th); 

  //Check case that m = infinity
  if (sinth == 0)
  {
    if (line_r > rect_xmin && line_r < rect_xmax)
    {
     *x0 = line_r; 
     *x1 = line_r; 
     *y0 = rect_ymin; 
     *y1 = rect_ymax; 
     return true; 
    }
    else
    {
      return false; 
    }
  }

  double costh = cos(line_th); 

  //case that m = 0
  if (costh == 0)
  {
    if (line_r > rect_ymin && line_r < rect_ymin) 
    {
       *y0 = line_r; 
       *y1 = line_r; 
       *x0 = rect_xmin; 
       *x1 = rect_xmax; 
       return true; 
    }
  }
  
 //evaluate line 
  double m =  -1/tan(line_th); 
  double b = line_r/sin(line_th); 
  

  double y_at_xmin = m * rect_xmin + b; 
  double y_at_xmax = m * rect_xmax + b; 
  double x_at_ymax = (rect_ymax - b ) / m; 
  double x_at_ymin = (rect_ymin - b ) / m;

  
  bool hits_W = (y_at_xmin > rect_ymin && y_at_xmin < rect_ymax); 
  bool hits_E = (y_at_xmax > rect_ymin && y_at_xmax < rect_ymax); 
  bool hits_S = (x_at_ymin > rect_xmin && x_at_ymin < rect_xmax); 
  bool hits_N = (x_at_ymax > rect_xmin && x_at_ymax < rect_xmax); 



  if (hits_W && hits_E)
  {
    *x0 = rect_xmin; 
    *y0 = y_at_xmin; 
    *x1 = rect_xmax; 
    *y1 = y_at_xmax; 
  }
  else if (hits_W && hits_N)
  {
    *x0 = rect_xmin; 
    *y0 = y_at_xmin; 
    *x1 = x_at_ymax; 
    *y1 = rect_ymax;
  }
  else if (hits_W && hits_S)
  {
    *x0 = rect_xmin; 
    *y0 = y_at_xmin; 
    *x1 = x_at_ymin; 
    *y1 = rect_ymin;
  }
  else if (hits_S && hits_N)
  {
    if (x_at_ymin < x_at_ymax)
    {
      *x0 = x_at_ymin;  
      *y0 = rect_ymin; 
      *x1 = x_at_ymax; 
      *y1 = rect_ymax;
    }
    else
    {
      *x1 = x_at_ymin;  
      *y1 = rect_ymin; 
      *x0 = x_at_ymax; 
      *y0 = rect_ymax;
    }
  }


  else if (hits_S && hits_E)
  {
    *x0 = x_at_ymin;  
    *y0 = rect_ymin; 
    *x1 = rect_xmax; 
    *y1 = y_at_xmax; 
  }
  else if (hits_N && hits_E)
  {
    *x0 = x_at_ymax;  
    *y0 = rect_ymax; 
    *x1 = rect_xmax; 
    *y1 = y_at_xmax; 
  }
  else return false; 

  return true; 


}

double DmtpcMath::round(double in, double roundto)
{

  if (roundto == 1.)
  {
    int intpart = int(in); 
    return in - intpart >= 0.5 ? intpart + 1: intpart; 
  }

  return round(in/roundto,1.) * roundto; 
}


#ifdef HAVE_MATHMORE
DmtpcMath::Rational::Rational(unsigned m, const double * mParams, unsigned n, const double * nParams) :
   ParFunc(m+n+1), fTop (m), fBottom (n), fBottomParams(n+1) 
{
  fM=m;
  fN=n;
  fBottomParams[0] = 1; 

  SetMParams(mParams); 
  SetNParams(nParams); 
}


DmtpcMath::Rational::Rational(unsigned m, unsigned n, const double * params) :
   ParFunc(m+n+1), fTop (m), fBottom (n), fBottomParams(n+1)  
{
  fM=m;
  fN=n;
  fBottomParams[0] = 1; 
  SetParameters(params); 
}


double DmtpcMath::Rational::DoEvalPar(double x, const double * p) const
{
  return gsl_poly_eval(p,fM+1, x)/ (gsl_poly_eval(p+fM+1, fN,x) *x + 1); 
}

double DmtpcMath::Rational::DoDerivative(double x) const
{
  
  std::cout << "DoDerivative(" << x << "); " << std::endl; 
  double deriv_top = fTop.Derivative(x); 
  double deriv_bottom = fBottom.Derivative(x); 
  double top = fTop(x); 
  double bottom = fBottom(x); 

  return (deriv_top * bottom - deriv_bottom * top) / (bottom*bottom); 

}

double DmtpcMath::Rational::DoParameterDerivative(double x, const double * p , unsigned int ipar) const
{
//  std::cout << "DoParameterDerivative(" << x << ",,"<<ipar<<"); "<<  std::endl; 
  if (ipar <= fM) return gsl_pow_int(x,ipar) / (gsl_poly_eval(p+fM+1, fN,x) *x + 1) ;  // p'/q
  else return -gsl_pow_int(x,ipar - fM) * gsl_poly_eval (p, fM+1,x)  / gsl_pow_int(gsl_poly_eval(p+fM+1, fN,x) *x + 1,2); 
}




const std::vector< std::complex <double> > &  DmtpcMath::Rational::FindPoles() 
{
  return fBottom.FindRoots(); 

}

const std::vector< std::complex <double> > &  DmtpcMath::Rational::FindRoots() 
{

  return fTop.FindRoots(); 
}

ROOT::Math::IGenFunction * DmtpcMath::Rational::Clone() const 
{
  return new Rational(fM,fN,Parameters()); 
}

void DmtpcMath::Rational::SetParameters(const double * params)
{
  if (params) 
  {
    for (unsigned i = 0; i <= fM; i++)
    {
      fParams[i] = params[i]; 
    }

    for (unsigned i = 0; i <fN; i++)
    {
      fParams[i+fM+1] = params[i+fM+1]; 
      fBottomParams[i+1] = params[i+fM+1]; 
    }

    std::cout << "Setting Parameters!" << std::endl; 
    for (unsigned i = 0; i < fM+fN+1; i++)
    {
      std::cout << params[i] << " "; 
    }
    std::cout << std::endl; 

    fTop.SetParameters(&fParams.front()); 
    fBottom.SetParameters(&fBottomParams.front()); 

  }
}
 

void DmtpcMath::Rational::SetMParams(const double * m)
{
  if (!m) return;
  for (unsigned i = 0; i <= fM; i++)
  {
    fParams[i] = m[i]; 
  }
  fTop.SetParameters(&fParams.front()); 
}

void DmtpcMath::Rational::SetNParams(const double * n)
{
   if (!n) return; 
   for (unsigned i = 0; i < fN; i++)
   {
      fParams[i+fM+1] = n[i]; 
      fBottomParams[i+1] = n[i]; 

   }

   fBottom.SetParameters(&fBottomParams.front()); 
}

void DmtpcMath::Rational::InitFitParameters(unsigned N, const double * x, const double *y, bool forcenopoles, unsigned ntries)
{
  unsigned size = fM+fN+1; 
  bool overconstrained = false; 
  if (N <= size)
  {
    size = N;
    overconstrained = true; 
    ntries = 1; 
  }


  TRandom * rand = new TRandom; 
  rand->SetSeed(); 


  std::vector<double> best_params(size); 
  double best_chisq = -1; 

  for (unsigned t = 0; t < ntries; t++)
  {
    bool done = false; 

    std::cout << t << std::endl; 
    while (!done)
    {

      //randomly choose points, unless overconstrained
      std::set<int> chosen; 
      if (overconstrained)
      {
        for (unsigned i = 0; i <N; i++)
        {
          chosen.insert(i); 
        }
      }
      else
      {

        while(chosen.size() < size) chosen.insert(rand->Integer(N)); 

      }


      //Find a solution with those points
      TMatrixD A(size,size); 
      TVectorD ys(size); 

      int r = 0; 
      for (std::set<int>::iterator it = chosen.begin(); it!=chosen.end(); it++)
      {
        int i = *it; 
//        std::cout <<"using point " << i << ": " << x[i] << "," << y[i] << std::endl; 

        for (unsigned j = 0; j < size; j++)
        {
          A[r][j] =  j <= fM ?  gsl_pow_int(x[i],j) : -gsl_pow_int(x[i], j - fM) * y[i]; 
        }

        ys[r] = y[i]; 

        r++; 
      }

//      A.Print(); 
//      ys.Print(); 

      TDecompSVD svd(A); 
      svd.Solve(ys); 
//      ys.Print(); 
    
      SetParameters(ys.GetMatrixArray()); 

      done = true; 
      if (!overconstrained && forcenopoles) 
      {
        const std::vector<std::complex<double> > poles  = FindPoles(); 

        for (unsigned p = 0; p < poles.size(); p++)
        {
          if (poles[p].imag() == 0 && poles[p].real() >= x[0] && poles[p].real() <= x[N-1])
          {
              done = false; 
              break;
          }
        }
      }
    }

    if (!overconstrained)
    {
      double chisq = 0; 
      for (unsigned i = 0; i < N; i++)
      {
        chisq += pow(DoEvalPar(x[i], Parameters()) - y[i],2); 
      }

      std::cout << chisq << std::endl; 

      if (t == 0 || chisq < best_chisq)
      {
        for (unsigned j = 0; j < size; j++)
        {
          best_params[j] = Parameters()[j]; 
        }
        best_chisq = chisq; 

      }
    }

    std::cout << std::endl; 
  }

  if (!overconstrained)
  {
    SetParameters(&best_params.front());   
  }

//  delete rand; 

}

//ClassImp(DmtpcMath::Rational); 

#endif 

double DmtpcMath::normalizeAngle(double phi, double center) 
{
 return phi - (2*M_PI) * TMath::Floor((phi + M_PI - center)/(2*M_PI)); 
}

double DmtpcMath::deltaAngle(double phi0,double phi1)
{
  return normalizeAngle(phi0,phi1) - phi0; 
} 

double DmtpcMath::integralOfLineSegmentConvolvedWithGaussian(double y0, double y1, double * par) 
{

 double m = par[0]; 
 double b = par[1]; 
 double x0 = par[2]; 
 double x1 = par[3]; 
 double w = par[4]; 

 //thank Mathematica CForm for this 
 return (sqrt(2/M_PI)*w*(2*b*(-exp(-pow(x0 - y0,2)/(2.*pow(w,2))) + exp(-pow(x1 - y0,2)/(2.*pow(w,2))) + 
         exp(-pow(x0 - y1,2)/(2.*pow(w,2))) - exp(-pow(x1 - y1,2)/(2.*pow(w,2)))) + 
         m*(-((x0 + y0)/exp(pow(x0 - y0,2)/(2.*pow(w,2)))) + (x1 + y0)/exp(pow(x1 - y0,2)/(2.*pow(w,2))) + 
        (x0 + y1)/exp(pow(x0 - y1,2)/(2.*pow(w,2))) - (x1 + y1)/exp(pow(x1 - y1,2)/(2.*pow(w,2))))) + 
        (2*b*(-x0 + y0) + m*(pow(w,2) - pow(x0,2) + pow(y0,2)))*TMath::Erf((x0 - y0)/(sqrt(2)*w)) - 
        (2*b*(-x1 + y0) + m*(pow(w,2) - pow(x1,2) + pow(y0,2)))*TMath::Erf((x1 - y0)/(sqrt(2)*w)) - 
        (2*b*(-x0 + y1) + m*(pow(w,2) - pow(x0,2) + pow(y1,2)))*TMath::Erf((x0 - y1)/(sqrt(2)*w)) + 
        (2*b*(-x1 + y1) + m*(pow(w,2) - pow(x1,2) + pow(y1,2)))*TMath::Erf((x1 - y1)/(sqrt(2)*w)))/4. ;
}


double DmtpcMath::lineSegmentConvolvedWithGaussian2D(double *xx, double * par)
{
  double x = xx[0]; 
  double y = xx[1]; 

  double m = par[0]; 
  double z0 = par[1]; 
  double x0 = par[2]; 
  double y0 = par[3]; 
  double x1 = par[4]; 
  double y1 = par[5]; 
  double w = fabs(par[6]); 


  //calculate intersection onto line: 
  double xint, yint; 
  
  double slope_proj = y1==y0 ? 1 : (x0-x1)/(y1-y0); 
  double inc = y1==y0? 0 : 1; 

  lineLineIntersection(x0,y0,x1,y1,x,y, x + inc, y + slope_proj,&xint,&yint); 

  double distance_to_line = sqrt(pow(x-xint,2) + pow(y-yint,2)); 
  double distance_from_start =  sqrt(pow(xint-x0,2) + pow(yint-y0,2)); 
  double segment_length = sqrt(pow(x1-x0,2) + pow(y1-y0,2)); 


  //value along line: 
  double par1d[5]; 
  par1d[0] = m; 
  par1d[1] = z0; 
  par1d[2] = 0; 
  par1d[3] = segment_length; 
  par1d[4] = w; 

  double x1d = ( xint > x0 == x1 > x0  || (xint == x0 && yint > y0)) ? distance_from_start : -distance_from_start; 
  double val1d = lineSegmentConvolvedWithGaussian(&x1d, par1d); 

  //value perpendicular to line
  
  double val2d = val1d * TMath::Gaus(distance_to_line,0,w,true); 
  //printf("(x,y,xint,yint,d2l,d2smval1d,val2d): (%f,%f,%f,%f,%f,%f,%f,%f)\n", x,y,xint,yint,distance_to_line,distance_from_start,val1d,val2d); 


  return val2d; 

}

double DmtpcMath::lineSegmentConvolvedWithGaussian(double * xx, double * par)
{
   double x = *xx; 
   double m = par[0]; 
   double b = par[1]; 
   double x0 = par[2]; 
   double x1 = par[3]; 
   double w = par[4]; 

//   std::cout << x  << " , " << m  << " , " << b << " , " << x0 << " , " << x1 << " , " << w << std::endl; 
   double val =   (exp((2.*x0*x - x0*x0 - x*x)/(w*w*2.)) - exp((2.*x1*x - x1*x1 - x*x) / (2.*w*w))) * m * w / sqrt(2.*M_PI)  
              + 0.5 * (b+m*x) * ( -TMath::Erf((x0-x)/(sqrt(2) * w)) + TMath::Erf((x1-x)/(sqrt(2) * w)) ); 

//   std:: cout << "val: " << val << std::endl; 
      
   return val; 
}

double DmtpcMath::sumVar(double dA, double dB, double cA, double cB) 
{ 
  return pow(dA*cA,2) + pow(dB*cB,2);
}

double DmtpcMath::sumError(double dA, double dB, double cA, double cB)
{
  return sqrt(sumVar(dA,dB,cA,cB)) ;
}

double DmtpcMath::productVar(double A, double B, double dA, double dB)
{
  return pow(A*B,2) * (pow(dA/A,2) + pow(dB/B,2));
}

double DmtpcMath::productError(double A, double B, double dA, double dB)
{ 
  return  A*B*sqrt(pow(dA/A,2) + pow(dB/B,2));
} 

double DmtpcMath::quotientVar(double A, double B, double dA, double dB)
{
  return pow(A/B,2) * (pow(dA/A,2) + pow(dB/B,2));
}

double DmtpcMath::quotientError(double A, double B, double dA, double dB) 
{
  return  A/B*sqrt(pow(dA/A,2) + pow(dB/B,2));
} 




double DmtpcMath::MVNpdf(size_t ndim, const double * x, const double * mu, const double * covar)
{
  TMatrixD covarM(ndim,ndim,covar); 
  TDecompChol chol(covarM); 
  TMatrixDSym inv; 
  chol.Invert(inv); 
  return MVNpdf(ndim,x,mu,&inv, covarM.Determinant()); 
}

double DmtpcMath::MVNpdf(size_t ndim, const double * x, const double * mu, const TMatrixDSym * invcovar, double det)
{

  TMatrixD v(ndim,1,x); 
  v -= TMatrixD(ndim,1,mu); 

  TMatrixD vT(1,ndim); 
  vT.Transpose(v); 

  TMatrixDSym work(*invcovar); 



  return 1./sqrt(pow(2*M_PI,ndim) * det) * exp(-(vT*work*v)(0,0)/2); 
}


