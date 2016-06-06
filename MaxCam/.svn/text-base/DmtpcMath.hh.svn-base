#ifndef __DMTPCMATH_HH__
#define __DMTPCMATH_HH__


#ifdef HAVE_MATHMORE
#include "Math/Polynomial.h"
#endif

#include <complex>
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include <cmath>

namespace DmtpcMath
{
  /*
    Constants (most defined in math.h but not on all systems
  */

  /** e */
  static const double EXP = 2.7182818284590452353602874713526625L;
  /** log(2) */
  static const double LOG2 = 0.6931471805599453094172321214581766L;
  /** log(10) */
  static const double LOG10 = 2.3025850929940456840179914546843642L;
  /** log_2(e) */
  static const double LOG2E = 0.6931471805599453094172321214581766L;
  /** log_10(e) */
  static const double LOG10E = 2.3025850929940456840179914546843642L;
  /** pi */
  static const double PI = 3.1415926535897932384626433832795029L;
  /** pi/2 */
  static const double PI_2 = 1.5707963267948966192313216916397514L;
  /** pi/4 */
  static const double PI_4 = 0.7853981633974483096156608458198757L;
  /** 1/pi */
  static const double ONE_PI = 0.3183098861837906715377675267450287L;
  /** 2/pi */
  static const double TWO_PI = 0.6366197723675813430755350534900574L;
  /** 2/sqrt(pi) */
  static const double TWO_SQRTPI = 1.1283791670955125738961589031215452L;
  /** sqrt(2) */
  static const double SQRT2 = 1.4142135623730950488016887242096981L;
  /** 1/sqrt(2) */
  static const double SQRT1_2 = 0.7071067811865475244008443621048490L;
  /** Euler-Mascheroni constant*/
  static const double EULER_GAMMA = 0.5772156659015328606065120900824023L;

  /** 1/sqrt(2*pi) **/ 
  static const double ONE_SQRT2PI = 0.3989422804; 

  /*
  Special functions
  */
  double Beta(double alpha,double beta);
  double Gamma(double x);
  double logGamma(double x);


  /*
  Distributions
  */

  double betaPdf(double x,double par[]);
  void betaPdfStats(double x[],double par[]);

  double chiPdf(double x,double par[]);
  void chiPdfStats(double x[],double par[]);

  double chiSquarePdf(double x,double par[]);
  void chiSquarePdfStats(double x[],double par[]);

  double crystalBallPdf(double *x,double par[]);

  double crystalBall(double x[], double par[]);

  double gammaPdf(double x,double par[]);
  void gammaPdfStats(double x[],double par[]);

  double inverseGausPdf(double x,double par[]);
  double inverseGausCdf(double x, double par[]);
  void inverseGausPdfStats(double x[],double par[]);


  double logisticPdf(double x,double par[]);
  void logisticPdfStats(double x[],double par[]);

  double logNormalPdf(double x,double par[]);
  double logNormalCdf(double x,double par[]);
  void logNormalPdfStats(double x[],double par[]); 

  double maxwellPdf(double x,double par[]);
  void maxwellPdfStats(double x[],double par[]);

  double normalPdf(double x,double par[]);
  double gausPdf(double x,double par[]);
  double normalCdf(double x,double par[]);
  void normalPdfStats(double x[],double par[]);

  double rayleighPdf(double x,double par[]);
  void rayleighPdfStats(double x[],double par[]);

  double studentPdf(double x, double par[]);
  void studentPdfStats(double x[],double par[]);

  /** Von mises distribution, 0 = scale, 1 = mean, 2 = dispersion */
  double vonMisesDist(double *x, double *par); 

  /** Von mises distribution w/ head tail fraction, 0 = HT, 1 = scale, 2 = mean, 3 = dispersion  */
  double vonMisesDistHT(double *x, double *par); 


  /** Wrapped normal distribution, 0 = scale, 1 = mean, 2 = sigma **/ 
  double wrappedNormalDist(double *x, double *par) ;
  /** Wrapped normal distribution w/ head tail fraction 0 = HT, 1 = scale, 2 = mean, 3 = sigma  **/ 
  double wrappedNormalDistHT(double *x, double *par) ;


  /** Wrapped cauchy distribution, 0 = scale, 1 = mean, 2 = sigma **/ 
//  double wrappedCauchyDist(double *x, double *par) ;
  /** Wrapped cauchy distribution w/ head tail fraction 0 = HT, 1 = scale, 2 = mean, 3 = sigma **/ 
//  double wrappedCauchyDistHT(double *x, double *par) ;

  //parameters:  m,b,x0,x1,sigma
  double lineSegmentConvolvedWithGaussian(double *x, double *pars);
  double lineSegmentConvolvedWithGaussian2D(double *x, double *pars);
  double integralOfLineSegmentConvolvedWithGaussian(double y0, double y1, double *pars);
  int lineLineIntersection(double x1,double y1, double x2, double y2, double x3, double y3, double x4, double y4,double * xint, double * yint); 





  static const int MEAN = 0;
  static const int MED = 1;
  static const int MODE = 2;
  static const int VAR = 3;
  static const int SKEW = 4;
  static const int KURT = 5;

  //calculate where a line hits a rectangle with sides orthogonal to coordinate system. Returns false if no intersection. Otherwise, x0,y0,x1,y1 are populated with the intersection points. 
  // parameterization: y = -cotan(theta) + r / sin(theta) 
  bool lineRectIntersectionPoints(double line_theta, double line_r, double rect_xmin, double rect_xmax, double rect_ymin, double rect_ymax, double * x0, double *y0, double * x1, double * y1); 


  double round(double in, double RoundToNearest = 1.); 

  double normalizeAngle(double phi, double center=0.); 
  double deltaAngle(double phi0,double phi1); 

  // propagation of error 
  double sumVar(double dA, double dB, double cA=1, double cB=1);
  double sumError(double dA, double dB, double cA=1, double cB=1);
  double productVar(double A, double B, double dA, double dB);
  double productError(double A, double B, double dA, double dB);
  double quotientVar(double A, double B, double dA, double dB);
  double quotientError(double A, double B, double dA, double dB);


  double MVNpdf( size_t ndim, const double * x, const double * mu, const double * covar); 
  double MVNpdf( size_t ndim, const double * x, const double * mu, const TMatrixDSym* invcovar, double det); 
  
#ifdef HAVE_MATHMORE
  /*
   * Function of form  (a0 + a1*x + a2*x^2 + ..) / (1 + b1*x + b2*x^2 + ..) 
   * m is numerator order, n is denominator
   * there are m + n + 1 parameters (since the zeroth order denominator coefficient is always 1)
   *
   *
   * Doesn't work in CINT right now 
   */
  class Rational 
#ifndef CINT
    : public ::ROOT::Math::ParamFunction< ::ROOT::Math::IParamGradFunction> , public ::ROOT::Math::IGradientOneDim
#endif
  {
      public:
        typedef ::ROOT::Math::ParamFunction< ::ROOT::Math::IParamGradFunction> ParFunc; 
        Rational(unsigned m, const double *mparams, unsigned n, const double * nparams); 
        Rational(unsigned m, unsigned n, const double * p = 0); 

        virtual ~Rational() {}; 

        const std::vector<std::complex <double> > & FindRoots(); 
        const std::vector<std::complex <double> > & FindPoles(); 

        void SetMParams(const double * m); 
        void SetNParams(const double * n); 

        unsigned int M() const { return fM; } 
        unsigned int N() const { return fN; } 

        void InitFitParameters(unsigned N, const double *x, const double *y, bool forcenopoles = true, unsigned ntries =10); 

        virtual void SetParameters(const double * params); 
        ROOT::Math::IGenFunction * Clone() const; 

        void FdF( double x, double & f, double & df) const 
        { 
          f = (*this)(x); 
          df = Derivative(x); 
        }

        const ROOT::Math::Polynomial & Top() { return fTop; } 
        const ROOT::Math::Polynomial & Bottom() { return fBottom; } 


      private: 
        unsigned fM; 
        unsigned fN; 
        ROOT::Math::Polynomial fTop; 
        ROOT::Math::Polynomial fBottom; 
        std::vector<double> fBottomParams; 


        double DoEvalPar(double x, const double * p) const; 
        double DoDerivative(double x) const; 
        double DoParameterDerivative(double x, const double * p, unsigned int ipar) const; 

//        ClassDef(Rational,0); 
  };

#endif

}


#endif
