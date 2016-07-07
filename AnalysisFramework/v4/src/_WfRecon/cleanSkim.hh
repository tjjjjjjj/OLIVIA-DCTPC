/*~~~~~ \(*^,^)/ Hi! This is where waveforms get analyzed! Unless you're fixing bugs or adding new algorithms, you probably just want WFknobs.temp. ~~~~~*/


#include <TMath.h>
#include <TObject.h>
#include <TString.h>
#include "TH1.h"
#include <cstdlib>
#include <vector>
using namespace std;
  
  void TurnKnobs();
void TurnKnobs(TString knobfile_);
  //Automatically sets the settings below

  void analyze();
  //runs analysis

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
