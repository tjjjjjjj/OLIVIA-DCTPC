/*~~~~~ \(*^,^)/ Hi! This is where waveforms get analyzed! Unless you're fixing bugs or adding new algorithms, you probably just want WFknobs.temp. ~~~~~*/

#ifndef WFRECON_HH
#define WFRECON_HH


#include <TMath.h>
#include <TObject.h>
#include <TString.h>
#include "TH1.h"
#include <cstdlib>
#include <vector>
using namespace std;

class WfRecon {

public:

  WfRecon();
  //constructor (empty)

  WfRecon(TString knobfile_);
  //constructor + initializer

  //WfRecon(TString rootfile_, TString braggfile_, TString knobfile_);
  //constructor + initializer + override auxiliary files

  virtual ~WfRecon();
  //destructor
  
  void TurnKnobs();
  void TurnKnobs(TString knobfile_)
  {knobfile = knobfile_;
    TurnKnobs();}
  //Automatically sets the settings below

  void analyze();
  //runs analysis


  /*~~~~~ \(*^w^)/ You can set all the settings on your own too! ~~~~~*/

  void SetRootFile(TString rootfile_){rootfile = rootfile_;}
  void SetBraggFile(TString braggfile_){braggfile = braggfile_;}
  void SetKnobFile(TString knobfile_){knobfile = knobfile_;}

  void SetDecayTime(double decaytime){Td = decaytime;}
  void SetTimeStep(double timestep){dt = timestep;}
  void SetBraggPeak(double braggpeak){peakstd=braggpeak;}
  void SetBraggVal(double braggval){braggmaxval=braggval;}
  void SetxScale(double xscale){xscalestd = xscale;}
  void SetEnergyMult(double energymult){energy_mult = energymult;}
  void SetHalfHeight(int halfheight_){halfheight = halfheight_;}

  void SetTermThresh(double termthresh_){termthresh = termthresh_;}
  void SetTermThresh2(double termthresh2_){termthresh2 = termthresh2_;}
  void SetTermThreshC(double termthreshc){termthreshcout=termthreshc;}
  void SetTerm2Pow(double term2pow_){term2pow = term2pow_;}
  void SetTSBCoeff(double TSB_coeff_){TSB_coeff = TSB_coeff_;}
  void SetTSBPow(double TSB_pow_){TSB_pow = TSB_pow_;}

  void SetDriftCor(double drift_cor){drift_correction = drift_cor;}
  void SetHeightMult(double height_mult_){height_mult = height_mult_;}
  void SetDeconvolution(int decon_on_){decon_on = decon_on_;}
  void SetDeconSD(double decon_SD_){decon_SD = decon_SD_;}
  void SetDeconOrder(int decon_ord_){decon_ord = decon_ord_;}

  void SetDerivOrder(int deriv_order_){deriv_order=deriv_order_;}
  void SetDerivStep(int deriv_step_){deriv_step=deriv_step_;}
  void SetDer2Step(int der2step){d2_step = der2step;}
  void SetDer4Step(int der4step){d4_step = der4step;}
  void SetDer6Step(int der6step){d6_step = der6step;}
  void SetDer8Step(int der8step){d8_step = der8step;}

  void SetRecSD(double recSD){rec_SD = recSD;}
  void SetRecSDMult(double recSDmult){rec_SD_mult = recSDmult;}
  void SetRecOverride(int recoverride){rec_SD_override = recoverride;}
  void SetRecCWidth(double reccwidth){rec_cwidth = reccwidth;}
  void SetRecCorrection(double reccor){rec_cor = reccor;}
  void SetMinSDMult(double minSDmult){min_mult = minSDmult;}
  void SetMinCWidth(double mincwidth){min_cwidth = mincwidth;}
  void SetMinCorrection(double mincor){min_cmult = mincor;}
  void SetGaussAcc(double gausacc){gaus_acc = gausacc;}

  void SetFinalSD(double finalSD){finalsmooth = finalSD;}
  void SetMFCfactor(double MFCfactor){mfcfactor = MFCfactor;}
  void SetDiffFactor(double diffFactor){diff_factor = diffFactor;}

  /*~~~~~ \(*^,^)/ Use these to read out results! ~~~~~*/

  //reconstructed energies
  double GetLeftInt(){return leftint;} //downward-traveling alpha
  double GetLeftInt2(){return leftint2;} //cheating (uses actual origin)
  double GetRightInt(){return rightint;} //upward-traveling alpha
  double GetRightInt2(){return rightint2;} //cheating

  //reconstructed and actual vertices
  int GetjMin(){return jmin;}
  int GetjBragg(){return jbragg;}
  int GetjDiff(){return jdiff;}
  int GetjTerm(){return jterm;}
  int GetjTerm2(){return jterm2;}
  int GetjTerm3(){return jterm3;}
  int GetjDev(){return maxdevloc;}
  int GetVertex(){return origin;}
  int GetMatchVectorSize(){return terminator_results.size();}
  int GetMatchVector(int i){return terminator_results.at(i);}

  //properties of reconstructed Bragg curves
  double GetSD(){return rec_SD;} //strength of smoothing filter used
  int GetDuration(){return wfd_delta;} //full width at 30% height
  int GetPeak1(){return peak1;} //Location of bragg peaks
  int GetPeak2(){return peak2;}
  double GetPeak1Val(){return peak1val;} //Height of peaks
  double getPeak2Val(){return peak2val;}
  int GetHalf1(){return half1;} //Location of 75% height on Bragg curves
  int GetHalf2(){return half2;}

  //These are measurements of how well the reconstructed energy-loss curve is approximated by a pair of best-fit Bragg curves.
  double GetTermFit(){return termdist;}
  double GetMaxDev(){return maxdev;}
  double GetLeftRMS(){return rms_left;}
  double GetRightRMS(){return rms_right;}
  double GetOuterRMS(){return rms_outer;}
  double GetFullRMS(){return rms_full;}

  //Other values that may be convenient to have in the same place
  double GetTimeStep(){return dt;}

private:

  /*~~~~~ (*o,o) lots of cool waveform analysis functions! ~~~~~*/

  double findRMS(int xmin, int xmax, TH1D* curve, TH1D* hist1, TH1D* hist2);
  //compute RMS error of curve. (Error is based on whichever of hist1 or hist2 is closer.)
  
  double maxdevsearch(int xmin, int xmax, TH1D* bragghist3, TH1D* bragghist4, TH1D* mfc2);
  //finds the highest distance from mfc2 to either of the bragghists.

  int terminatorsearch(TH1D* bragghist3, TH1D* bragghist4, TH1D* mfc2, int xmin, int xmax, double threshold, vector<int>& vec);
  //An algorithm sent from the future. Extremely adept at finding the vertex but has flaws that prevent it from achieving its goal on certain types of waveforms.

  int terminatorsearchbackup(TH1D* bragghist3, TH1D* bragghist4, TH1D* mfc2, int xmin, int xmax);
  //Returns an answer in the event that terminatorsearch does not return one.

  int terminatorsearch2(TH1D* bragghist3, TH1D* bragghist4, TH1D* mfc2, int xmin, int xmax, double threshold);
  //Modified version of the Terminator. Generally gives poor results, but tends to give decent results precisely when Terminator gives bad ones.

  int terminatorsearch3(TH1D* bragghist3, TH1D* bragghist4, int xmin, int xmax);
  //A crude but consistent version of Terminator2. Not recommended on its own. May be useful in combination with other algorithms.

  void plotbragg(TH1D* bragghist3, TH1D* bragghist2, double xscaletrue, double xscale, double yscale, int xoff);
  //Scales and shifts a bragg curve to match the observed event.

  void plotbraggrev(TH1D* bragghist4, TH1D* bragghist2, double xscaletrue, double xscale, double yscale, int xoff);
  //For bragg curves that go left-to-right instead of right-to-left.

  double integrate(TH1D* hist, int imin, int imax);
  //Gets area under hist.

  int maxsearch(TH1D* hist, int imin, int imax);
  //Finds max bin of hist within given limits.

  int halfsearch(TH1D* hist, int imin, int imax, double maxval);
  //Finds where a hist reaches half of its maxval within given limits.
  //Note: Currently programmed to actually look for 75% height.

  int absmin(TH1D* hist, int i, int range);
  //Finds max bin of hist within range bins of i.

  double deriv(TH1D* hist, int n, int k, double h);
  //returns first derivative of step size k at bin n with hist step size h.

  void deconvolve(TH1D* hist, double SD, int order);
  //Un-does Gaussian smearing. May contain traces of sorcery.

  double deriv2(TH1D* hist, int n, int k, double h);
  //Second derivative.

  double deriv4(TH1D* hist, int n, int k, double h);
  //Fourth derivative.

  double deriv6(TH1D* hist, int n, int k, double h);
  //Sixth derivative.

  double deriv8(TH1D* hist, int n, int k, double h);
  //Eighth derivative.
  
  int minfinder(TH1D* wfder);
  int minfinder(TH1D* wf, TH1D* wfder);
  //Finds the minimum of certain processed versions of the waveform. Provides weak estimates of the vertex.

  void minprepare(TH1D* wf, TH1D* minfinder_curve);
  //processes a waveform into a format for minfinder.

  double average(int j, vector<double>& vec, TH1D* hist);
  //Averages points on a hist according to weights in vec.

  double median(vector<double>& vec);
  //Returns median entry of a vector.

  double normaldist(double peak, double x, double SD);
  //evaluates a normal distribution at x with given SD and peak location.

  void normalize(vector<double>&vec);
  //Does what it says.

  void makekernel(vector<double>&vec, double SD, double acc);
  //Fills vec with a Gaussian distribution for use with average().

  void smooth(TH1D* hist, double SD, double acc);
  //Applies Gaussian blur to a histogram. Used as a low-pass filter.

  void correct(TH1D* hist, double SD, double acc, double factor);
  //Un-does some undesirable effects of smooth().

  void shifthist(TH1D* hist, int plot_offset);
  //Shifts the values in a histogram on the x-axis.
  
  double histspline(TH1D*hist, int xmin, int xoff, int x, double xscale, double yscale);
  //An interface to read decimal values out of a histogram as if it were a spline curve.


  /*~~~~~ (o.o*) Whoa. I wonder what all these variables do...? ~~~~~*/

  /*~~~~~\(^.^*)/ I know! I'll go read WFknobs.temp! ~~~~~*/

  TString rootfile;
  TString braggfile;
  TString knobfile;

  int wfd_delta;
  int wfd_delta2;
  double vertoff;
  int plotoff;
  int peak1;
  int peak2;
  double peak1val;
  double peak2val;
  int half1;
  int half2;

  double mfc_int;
  double sec1;
  double sec2;
  double sec1true;
  double sec2true;
  double cos1;
  double cos2;
  double cos1true;
  double cos2true;
  double bestfit;
  double maxdev;
  int maxdevloc;
  vector<int> terminator_results;

  double dt;
  double Td;

  double peakstd;
  double braggmaxval;
  double xscalestd;
  double energy_mult;
  double yscalestd;
  int halfheight;

  double termthresh;
  double termthresh2;
  double termthreshcout;
  double termdist;
  double term2pow;
  double TSB_coeff;
  double TSB_pow;

  double drift_correction;
  double height_mult;

  int decon_on;
  double decon_SD;
  int decon_ord;

  int deriv_order;
  int deriv_step;
  int d2_step;
  int d4_step;
  int d6_step;
  int d8_step;

  double rec_SD;
  double rec_SD_mult;
  int rec_SD_override;
  double rec_cwidth;
  double rec_cor;
  double min_mult;
  double min_cwidth;
  double min_cmult;
  double gaus_acc;

  double finalsmooth;
  double mfcfactor;
  double diff_factor;

  TH1D* bragghist;
  TH1D* bragghist2;
  TH1D* bragghist3;
  TH1D* bragghist4;
  TH1D* rough_wf;
  TH1D* histA3;
  TH1D* histB3;
  TH1I* overkill_histogram;
  TH1D* overkill_histogram2;
  TH1D* histsum3;
  TH1D* test_rec_wf;
  TH1D* test_minfinder_curve;
  TH1D* rec_wf;
  TH1D* minfinder_curve;
  TH1D* mfc2;
  TH1D* mfcdiff;

  vector<double> minima;
  int jpeaksearch;
  int jcent;
  int zeroA;
  int zeroB;
  int origin;
  double rms_left;
  double rms_right;
  double rms_outer;
  double rms_full;

  double leftint;
  double leftint2;
  double rightint;
  double rightint2;
  int jmin;
  int jbragg;
  int jdiff;
  int jterm;
  int jterm2;
  int jterm3;

  ClassDef(WfRecon,1);
};

#endif
