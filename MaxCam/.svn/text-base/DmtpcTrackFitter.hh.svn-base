#ifndef DMTPC_TRACK_FITTER_HH
#define DMTPC_TRACK_FITTER_HH



#include "TF2.h"
#include "TH2.h"
#include "TCanvas.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/WrappedFunction.h"
#include "Math/AdaptiveIntegratorMultiDim.h"
#include "DmtpcSkimEvent.hh"
#include "TCanvas.h"
#include "MaxCamSRIM.hh"
#include "TArrow.h"
#include <vector>

class DmtpcTrackFitter
{
  public: 
    enum PARTICLE
    {
      ALPHA,
      FLUORINE,
      CARBON  
    }; 

    struct Param
    {
      Param(){;} 
      double val; 
      double err; 
      double errUp; 
      double errDn; 

      void zeroErr() {errUp=0;errDn=0;} 
      void fixErr() { if(errUp==0 || isnan(errUp)){errUp = err;}  if(errDn==0 || isnan(errDn)) {errDn = -err;} } 
      ClassDefNV(Param,1); 
    }; 


    struct Result
    {
      Result() {;}
      Param phi; 
      Param htphi; 
      Param fitE; 
      Param E; 
      Param Ephys; 
      Param range; 
      Param sigma; 
      Param x0; 
      Param y0; 
      Param z0; 
      Param z1; 
      Param delta_z; 

      double chisq; 
      double gain; 
      double likelyIntegralSame; 
      double likelyIntegralSameError; 
      double likelyIntegralTotal; 
      double likelyIntegralTotalError; 
      double prob;
      double probError; 
      double probFast; 
      int ndof; 
      bool success; 
      bool flip; 
      bool minos; 
      bool fullProb; 
      ClassDefNV(Result,1); 
    }; 

    DmtpcTrackFitter( MaxCamSRIM * sr ,  
                    double rangecal=0.16 /* mm/px */,
                    double gain = 18 /* keV/adu */,
                    const TH2 * gainmap=0,
                    double min_gain = 0.2,
                    double read_noise = 10 
                    );  

    virtual ~DmtpcTrackFitter(); 

    int fit(const TH2 * img, const vector<int> * pix, double phi_guess = DBL_MAX, double range_guess = DBL_MAX, double sigma_guess = DBL_MAX); 
    int fit(DmtpcSkimEvent * ev, int c, int t); 

    void setMinimizationMaxCalls(int ncalls){ min->SetMaxFunctionCalls(ncalls); }
    void setIntegrationAbsTol(double val) { integral_abstol = val; } 
    void setIntegrationRelTol(double val) { integral_reltol = val; } 
    void setIntegrationNCalls(int val) { integral_ncalls = val; } 

    void setClusterPadding(int val) { padding = val; } 
    void setZeroOutsideCluster(bool val) { zero_out = val; } 
    void setDraw(bool val = true) { draw = val; } 
    void setVerbose(bool val = true) { verbose = val; } 
    void setDoIntegral(bool val = true) { do_integral = val; } 
    void setAlwaysDoMinos(bool val = true) { always_minos = val; } 
    void setMinosThresh(double val = true) { minos_thresh = val; } 
    void setIntegraNSigma(double val) { integral_nsigma = val; } 
    const Result * getResult() const {return result; }  
    void setQuickEval(bool q) { quickeval = q; } 
    TCanvas * getCanvas() {return canvas; } 

    double SRIMLineFit2DFn(const double * xx) const; 


    private: 
      void fillResultParams(); 
      ROOT::Math::Minimizer *min; 
      ROOT::Math::Functor minfn; 
      double gain; 
      double lengthcal; 
      const TH2 * gainMap; 
      double min_gain;
      Result *result; 
      mutable TF2* fitfn; 
      mutable ROOT::Math::WrappedMultiFunction<TF2&> * wrapped_fitfn; 
      ROOT::Math::AdaptiveIntegratorMultiDim * fit_integrator; 
      MaxCamSRIM *srim; 
      TCanvas * canvas; 
      double minheight; 
      double maxheight; 
      double binwidth; 
      double read_noise; 
      TH2* fithist; 
      TH2* local_gainmap; 
      int padding; 
      bool zero_out; 
      bool verbose; 
      bool draw; 
      bool do_integral; 
      double integral_abstol; 
      double integral_reltol; 
      double integral_nsigma; 
      int integral_ncalls; 
      bool always_minos; 
      double minos_thresh; 
      TArrow * arrow; 

      bool quickeval; 

      ClassDef(DmtpcTrackFitter,0); 
};





#endif
