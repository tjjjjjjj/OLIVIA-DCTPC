#ifndef DMTPC_PROJECTION_HH
#define DMTPC_PROJECTION_HH

#include <vector>
#include "DmtpcSRIMProjection.hh"
#include <string>
#include "TObject.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include <math.h>

class DmtpcGainMap; 
class DmtpcSkimEvent; 
class TH2; 
class MaxCamClusterImage; 
class MaxCamSRIM; 
class TGraphErrors; 
class TF1; 

using std::vector;
using std::string;

class DmtpcProjection : public TNamed
{
  public: 

      struct LineFitParams
      {
        double phi; 
        double slope; 
        double alpha; 
        double alpha_error; 
        double b; 
        double E; 
        double chisq; 
        double transverse_chisq; 
        double range; 
        double range_error; 
        double offset; 
        double offset_error; 
        double sigma; 
        double sigma_error; 
        int ndof; 
        int transverse_ndof; 
        bool success; 
      }; 

      struct SRIMLineFitParams
      {
        double phi; 
        double sinThetaMin;
        double sinThetaSlope; 
        double gain; 
        double chisq; 
        double transverseChisq; 
        double lengthcal; 
        double E; 
        double EPhys; 
        double EError; 
        double EPhysError; 
        double range; 
        double rangeError; 
        double sigma; 
        double sigmaError; 
        double x0; 
        double y1; 
        double y0; 
        double x0Error; 
        double x0ErrorLow; 
        double x0ErrorHigh; 
        double y0Error; 
        double y0ErrorLow; 
        double y0ErrorHigh; 
        double y1Error; 
        double y1ErrorLow; 
        double y1ErrorHigh; 
        double likelyIntegralAll; 
        double likelyIntegralSame; 
        double prob; 
        int ndof; 
        int transverseNdof; 
        bool success; 
      }; 


      struct SRIMFitParams
      {
        double phi; 
        double sin_theta; 
        double gain; 
        double chisq; 
        double lengthcal; 
        double pressure; 
        double E; 
        double range; 
        double sigma; 
        double x0; 
        int ndof; 
        DmtpcSRIMProjection::PARTICLE type; 
        bool success; 
      }; 

      DmtpcProjection(DmtpcSkimEvent * e, int c, int i, const DmtpcGainMap * gm=0, bool flip = false, double spacer_sigma=1, int ninterp = 10);
      DmtpcProjection(const MaxCamClusterImage * clust, int track, double phi, double x, double y, const DmtpcGainMap * gm =0, bool flip = false, double spacer_sigma=1, int ninterp = 10);
      
      DmtpcProjection();
      virtual ~DmtpcProjection();

      const TH1 * getLongitudinalProfile() const {return _long_profile;}
      const TH1 * getInterpolatedLongitudinalProfile() const {return _long_profile_interpolated;}
      const TH1 * getTransverseProfile() const {return _tranv_profile;}
      TH1 * getTransverseProfile() {return _tranv_profile;}
//      const TH1 * getTransverseSlice(int i) const {return _tranv_slices[i];}
//      unsigned getNTransverseSlices() const {return _tranv_slices.size();}
      const TGraphErrors * getEInterpolatedProfile() const {return _EProfile; } 
      const TGraphErrors * getEFitProfile() const {return _EFitProfile; } 
      TGraphErrors * getEFitProfile() {return _EFitProfile; } 

      double getWidth() const { return _width; } 
      TF1 * longitudinalFit(TF1 * function,  vector<string> *  minuit_cmds = 0); 
      TF1 * EFit(TF1 * function, TF1 * rangefn_guess = 0, double jump_limits = 0.2, double jump_error = 0.05,  vector<string> *  minuit_cmds = 0); 

      //Needed for TMinuit fitting. 
      TF1 * getUserFunc() { return _userfunc; } 

      int getNJumps() const { return _njumps; } 
      int getJumpLowBin(int i) const;
      int getJumpHighBin(int i) const;

      int getNStartIgnore() const { return  _start_cutoff ? _long_profile->FindBin(_width) : 0 ; }
      int getNEndIgnore() const { return  _end_cutoff ? _long_profile->FindBin(_width ) : 0 ; }

      void setStartCutoff(bool start_cutoff) { _start_cutoff = start_cutoff; } 
      void setEndCutoff(bool end_cutoff) { _end_cutoff = end_cutoff; } 

      double getLongitudinalSpacer(int i) const { return _long_spacers[i]; }
      double getLongitudinalSpacerWidth(int i) const { return _long_spacer_widths[i];}
      double getNLongitudinalSpacers() const { return (int) _long_spacers.size(); }
      double getSpacerSigma() const { return _spacer_sigma; }

      LineFitParams *doLineFit(const char * opt = "") const; 
      SRIMFitParams *doSRIMFit(DmtpcSRIMProjection::PARTICLE type, double lengthcal, double gain, double pressure,  const char * opt = "") const; 
      SRIMLineFitParams *doSRIMLineFit(MaxCamSRIM * sr, double lengthcal, double gain, const char * opt = "") const; 


      double lineFitFn(const double * xx) const; 
      double SRIMfitFn(const double * xx) const; 
      double SRIMLineFitFn(const double * xx) const; 
     
  private: 
      TH1 * _long_profile; 
      TH1 * _long_profile_interpolated; 
      TH1 * _tranv_profile; 
//      vector<TH1 *> _tranv_slices; 
    
      TF1 * _userfunc;   //!
      double _spacer_sigma; 
      vector<double> _long_spacers; 
      vector<double> _tranv_spacers; 
      vector<double> _long_spacer_widths; 
      vector<double> _tranv_spacer_widths; 
      TGraphErrors * _EProfile; 
      TGraphErrors * _EFitProfile; 
      int _njumps; 
      double _width; 
      double _phi;
      bool _start_cutoff; 
      bool _end_cutoff; 
      bool _flip_spacers; 
      double _max_interpolated; 
      int _theta_param;

      mutable MaxCamSRIM * _srim; //!  //for fitting; 
      
      void init(const MaxCamClusterImage * clust, int track, double phi, double x, double y, const DmtpcGainMap * gm, bool flip = false, double spacer_sigma=1, int ninterp = 10);

      ClassDef(DmtpcProjection,4); 
};


#endif 
