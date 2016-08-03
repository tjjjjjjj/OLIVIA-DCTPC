#include "TCanvas.h"
#include "WaveformTools.hh"
#include "WaveformAnalysis.hh"
#include "CspWaveform.hh"
#include "CspPulse.hh"
#include "PMTWaveform.hh"
#include "PMTPulse.hh"
#include "FastWaveform.hh"
#include "FastPulse.hh"
#include <TFile.h>
#include <TH1.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
using namespace std;
using namespace waveform;
#define max(x,y) x>y?x:y
#define min(x,y) x<y?x:y
#define hget(h,i) h->GetBinContent(i)
#define hcenter(h,i) h->GetXaxis()->GetBinCenter(i)
#define verase(v,i) v.erase(v.begin()+i)

// Begin Big-Ugly-Double-Alpha-Analysis-Definitions
// (It's long. Do a text-search for "End Big-" to find the end.)

/* Waveform reconstruction variables/functions */

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
vector<double> terminator_values;

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

double termderivfactor;
double dermfcfac;
double dermfcfac2;
double dermfcfac3;

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

void TurnKnobs(int runnum);
void analyze();
double findRMS(int, int, TH1D*, TH1D*, TH1D*);
double maxdevsearch(int, int, TH1D*, TH1D*, TH1D*);
int terminatorsearch(TH1D*, TH1D*, TH1D*, int, int, double, vector<int>&, vector<double>&);
int terminatorsearchbackup(TH1D*,TH1D*,TH1D*,int,int);
int terminatorsearch2(TH1D*,TH1D*,TH1D*,int,int,double);
int terminatorsearch3(TH1D*,TH1D*,int,int);
void plotbragg(TH1D*,TH1D*,double,double,double,int);
void plotbraggrev(TH1D*,TH1D*,double,double,double,int);
double integrate(TH1D*,int,int);
int maxsearch(TH1D*,int,int);
int halfsearch(TH1D*,int,int,double);
int absmin(TH1D*,int,int);
double deriv(TH1D*,int,int,double);
void deconvolve(TH1D*,double,int);
double deriv2(TH1D*,int,int,double);
double deriv4(TH1D*,int,int,double);
double deriv6(TH1D*,int,int,double);
double deriv8(TH1D*,int,int,double);
int minfinder(TH1D*);
int minfinder(TH1D*,TH1D*);
void minprepare(TH1D*,TH1D*);
double average(int,vector<double>&,TH1D*);
double median(vector<double>&);
double normaldist(double,double,double);
void normalize(vector<double>&);
void makekernel(vector<double>&,double,double);

void smooth(TH1D*,double,double);
void correct(TH1D*,double,double,double);
void shifthist(TH1D*,int);
double histspline(TH1D*,int,int,int,double,double);

void TurnKnobs(int runnum)
{
  ifstream ifstr(knobfile);
  if (!ifstr.is_open()){
    cout << knobfile << " not found." << endl;
  }
  string line;
  TString name;
  TString s;
  Int_t n;
  Double_t d;
  while (!ifstr.eof()){
    getline(ifstr,line);
    if (line[0]=='#' || line=="") continue;
    istringstream linestr(line);
    linestr >> name;

    if (name == "RootFile"){

      std::stringstream tmpstr;
      tmpstr.str("");
      linestr >> s;
      tmpstr << s << runnum << ".root";
      rootfile = tmpstr.str();

      cout << endl << endl << endl << rootfile << endl << endl << endl;
    }
    else if (name == "TermDerFac"){
      linestr >> d;
      termderivfactor = d;
    }
    else if (name == "DMFCFac"){
      linestr >> d;
      dermfcfac = d;
    }
    else if (name == "DMFCFac2"){
      linestr >> d;
      dermfcfac2 = d;
    }
    else if (name == "DMFCFac3"){
      linestr >> d;
      dermfcfac3 = d;
    }
    else if (name == "BraggFile"){
      linestr >> s;
      braggfile = s;
    }
    else if (name == "DecayTime"){
      linestr >> d;
      Td = d;
    }
    else if (name == "TimeStep"){
      linestr >> d;
      dt = d;
    }
    else if (name == "BraggPeak"){
      linestr >> d;
      peakstd = d;
    }
    else if (name == "BraggMax"){
      linestr >> d;
      braggmaxval = d;
    }
    else if (name == "xScale"){
      linestr >> d;
      xscalestd = d;
    }
    else if (name == "EnergyMult"){
      linestr >> d;
      energy_mult = d;
    }
    else if (name == "HalfHeight"){
      linestr >> n;
      halfheight = n;
    }
    else if (name == "TermThresh"){
      linestr >> d;
      termthresh = d;
    }
    else if (name == "TermThresh2"){
      linestr >> d;
      termthresh2 = d;
    }
    else if (name == "TermThreshC"){
      linestr >> d;
      termthreshcout = d;
    }
    else if (name == "TSBCoeff"){
      linestr >> d;
      TSB_coeff = d;
    }
    else if (name == "TSBPow"){
      linestr >> d;
      TSB_pow = d;
    }
    else if (name == "DriftCor"){
      linestr >> d;
      drift_correction = d;
    }
    else if (name == "HeightMult"){
      linestr >> d;
      height_mult = d;
    }
    else if (name == "Deconvolution"){
      linestr >> n;
      decon_on = n;
    }
    else if (name == "DeconSD"){
      linestr >> d;
      decon_SD = d;
    }
    else if (name == "DeconOrder"){
      linestr >> n;
      decon_ord = n;
    }
    else if (name == "DerivOrder"){
      linestr >> n;
      deriv_order = n;
    }
    else if (name == "DerivStep"){
      linestr >> n;
      deriv_step = n;
    }
    else if (name == "Der2Step"){
      linestr >> n;
      d2_step = n;
    }
    else if (name == "Der4Step"){
      linestr >> n;
      d4_step = n;
    }
    else if (name == "Der6Step"){
      linestr >> n;
      d6_step = n;
    }
    else if (name == "Der8Step"){
      linestr >> n;
      d8_step = n;
    }
    else if (name == "RecSD"){
      linestr >> d;
      rec_SD = d;
    }
    else if (name == "RecSDMult"){
      linestr >> d;
      rec_SD_mult = d;
    }
    else if (name == "RecOverride"){
      linestr >> n;
      rec_SD_override = n;
    }
    else if (name == "RecCWidth"){
      linestr >> d;
      rec_cwidth = d;
    }
    else if (name == "RecCorrection"){
      linestr >> d;
      rec_cor = d;
    }
    else if (name == "MinSDMult"){
      linestr >> d;
      min_mult = d;
    }
    else if (name == "MinCWidth"){
      linestr >> d;
      min_cwidth = d;
    }
    else if (name == "MinCorrection"){
      linestr >> d;
      min_cmult = d;
    }
    else if (name == "GaussAcc"){
      linestr >> d;
      gaus_acc = d;
    }
    else if (name == "FinalSD"){
      linestr >> d;
      finalsmooth = d;
    }
    else if (name == "MFCfactor"){
      linestr >> d;
      mfcfactor = d;
    }
    else if (name == "DiffFactor"){
      linestr >> d;
      diff_factor = d;
    }
  }
}

void analyze()
//runs analysis                                                                                                           
{
  term2pow = .5;
  yscalestd = xscalestd * peakstd / braggmaxval;

  TFile fbragg(braggfile);
  bragghist = (TH1D*)fbragg.Get("bragghist");
  bragghist2 = (TH1D*)fbragg.Get("bragghist2");
  bragghist3 = new TH1D("asdf1","jkl1",16384,0,16384);
  bragghist4 = new TH1D("asdf2","jkl2",16384,0,16384);

  TFile f(rootfile);
  rough_wf = (TH1D*)f.Get("hist1");
  histA3 = (TH1D*)f.Get("first_hist2");
  histB3 = (TH1D*)f.Get("first_hist3");

  overkill_histogram = (TH1I*)f.Get("plot_offset");
  overkill_histogram2 = (TH1D*)f.Get("vertical_offset");
  plotoff = overkill_histogram->GetBinContent(0);
  vertoff = overkill_histogram2->GetBinContent(0);

  shifthist(histA3,plotoff);
  shifthist(histB3,plotoff);

  histsum3 = new TH1D("asdf3","jkl3",16384,0,16384);
  for (int i = 0; i < 16384; i++)
    histsum3->SetBinContent(i,histA3->GetBinContent(i)+histB3->GetBinContent(i));

  test_rec_wf = new TH1D("test_rec_wf","jkl4",16384,0,16384);
  for (int j = 0; j < 16384; j++) test_rec_wf->SetBinContent(j,rough_wf->GetBinContent(j));
  for (int j = 0; j < 1; j++) smooth(test_rec_wf,rec_SD,gaus_acc);
  for (int j = 0; j < 1; j++) correct(test_rec_wf,rec_cwidth * rec_SD,gaus_acc,rec_cor);

  test_minfinder_curve = new TH1D("asdf5","jkl5",16384,0,16384);
  minprepare(test_rec_wf,test_minfinder_curve);
  wfd_delta2 = int(wfd_delta * dt * 1e9);
  if (rec_SD_override == 0)
    {
      if (wfd_delta2 < 500)
	rec_SD = 10.  + wfd_delta2 * 1./50.;
      else if (wfd_delta2 < 1500)
	rec_SD = 20. + (wfd_delta2 -  500) / 100.;
      else if (wfd_delta2 < 2500)
	rec_SD = 30. + (wfd_delta2 - 1500) / 100.;
      else
	rec_SD = 40. + (wfd_delta2 - 2500) / 200.;
      rec_SD *= rec_SD_mult;
    }

  rec_wf = new TH1D("rec_wf","jkl6",16384,0,16384);
  for (int j = 0; j < 16384; j++) rec_wf->SetBinContent(j,rough_wf->GetBinContent(j));
  for (int j = 0; j < 1; j++) smooth(rec_wf,rec_SD,gaus_acc);
  for (int j = 0; j < 1; j++) correct(rec_wf,rec_cwidth * rec_SD,gaus_acc,rec_cor);

  jmin = 0;
  jbragg = 0;
  jdiff = 0;
  minfinder_curve = new TH1D("asdf7","jkl7",16384,0,16384);
  minprepare(rec_wf, minfinder_curve);
  jmin = minfinder(rec_wf,minfinder_curve);
  mfc2 = new TH1D("asdf8","jkl8",16384,0,16384);
  for (int i = 0; i < 16384; i++)
    mfc2->SetBinContent(i,minfinder_curve->GetBinContent(i));
  if (decon_on == 1)
    {
      deconvolve(mfc2,decon_SD,decon_ord);
      smooth(mfc2,finalsmooth,gaus_acc);
      correct(mfc2,finalsmooth*1.05,gaus_acc,.96);
    }
  mfcdiff = new TH1D("asdf9","jkl9",16384,0,16384);
  for (int i = 0; i < 16384; i++)
    mfcdiff->SetBinContent(i,mfc2->GetBinContent(i)- mfcfactor*minfinder_curve->GetBinContent(i));

  smooth(mfcdiff,6,gaus_acc);
  correct(mfcdiff,6,gaus_acc,1);
  jdiff = absmin(mfcdiff,jmin,60);

  jbragg = minfinder(mfc2);

  minima.push_back(jbragg);
  minima.push_back(jmin);
  minima.push_back(jdiff);
  jpeaksearch = jbragg;
  //jpeaksearch = int(median(minima) + .5);   
  peak1 = maxsearch (mfc2,2000,jpeaksearch);
  peak2 = maxsearch (mfc2,jpeaksearch,7000);
  peak1val = mfc2->GetBinContent(peak1);
  peak2val = mfc2->GetBinContent(peak2);

  
  half1 = halfsearch(mfc2,2000,peak1,peak1val);
  half2 = halfsearch(mfc2,peak2,10000,peak2val);
  sec1 = peak1val / peakstd * height_mult;
  sec1true = sec1;
  sec2 = peak2val / peakstd * height_mult;
  sec2true = sec2;

  cos1 = 1./sec1;
  cos1true = 1./sec1 + drift_correction;
  cos2 = 1./sec2;
  cos2true = 1./sec2 - drift_correction;

  plotbragg(bragghist3,bragghist2,cos1true*xscalestd,cos1*xscalestd,yscalestd,int(half1 - (halfheight * (cos1true*xscalestd))+.5));
  plotbraggrev(bragghist4,bragghist2,cos2true*xscalestd,cos2*xscalestd,yscalestd,int(half2 + (halfheight * (cos2true*xscalestd))+.5));

  jterm  =  terminatorsearch(bragghist3,bragghist4,mfc2,peak1,peak2,termthresh,terminator_results,terminator_values);
  jterm2 = terminatorsearch2(bragghist3,bragghist4,mfc2,peak1,peak2,termthresh2);
  jterm3 = terminatorsearch3(bragghist3,bragghist4,peak1,peak2);

  maxdev = maxdevsearch(peak1,peak2,bragghist3,bragghist4,mfc2);

  zeroA = 0;
  zeroB = 0;
  origin = 0;
  for (int i = 0; i < 16384; i++)
    {
      if (histA3->GetBinContent(i) > 0) zeroA = 1;
      if (histB3->GetBinContent(i) > 0) zeroB = 1;
      if (zeroA == 0 | zeroB == 0) origin = i;
    }

  jcent = jterm;

  leftint  = histspline(bragghist,half1,halfheight,jcent,1/(cos1 * xscalestd),1);
  rightint = histspline(bragghist,half2,halfheight,jcent,-1/(cos2 * xscalestd),1);
  leftint *= energy_mult;
  rightint *= energy_mult;

  jcent = origin;
  leftint2  = histspline(bragghist,half1,halfheight,jcent,1/(cos1true * xscalestd),1);
  rightint2 = histspline(bragghist,half2,halfheight,jcent,-1/(cos2true * xscalestd),1);
  leftint2 *= energy_mult;
  rightint2 *= energy_mult;

  rms_left = findRMS(half1, half1 + (half2 - half1)/4, mfc2, bragghist3, bragghist4);
  rms_right = findRMS(half1 + 3*(half2 - half1)/4, half2, mfc2, bragghist3, bragghist4);
  rms_outer = sqrt((rms_left*rms_left + rms_right*rms_right)/2);
  rms_full = findRMS(half1,half2,mfc2,bragghist3,bragghist4);

  delete bragghist;
  delete bragghist2;
  delete bragghist3;
  delete bragghist4;
  delete rough_wf;
  delete histA3;
  delete histB3;
  delete overkill_histogram;
  delete overkill_histogram2;
  delete histsum3;
  delete test_rec_wf;
  delete test_minfinder_curve;
  delete rec_wf;
  delete minfinder_curve;
  delete mfc2;
  delete mfcdiff;

  fbragg.Close();
  f.Close();
}

/*~~~~~ (*o,o) lots of cool waveform analysis functions! ~~~~~*/

double findRMS(int xmin, int xmax, TH1D* curve, TH1D* hist1, TH1D* hist2)
//compute RMS error of curve. (Error is based on whichever of hist1 or hist2 is closer.)                                  
{
  if (xmax <= xmin) xmax = xmin + 1;
  double sum = 0;
  double diff1;
  double diff2;
  double diff;
  for(int x = xmin; x < xmax; x++)
    {
      diff1 = fabs(pow(curve->GetBinContent(x) - hist1->GetBinContent(x),2));
      diff2 = fabs(pow(curve->GetBinContent(x) - hist2->GetBinContent(x),2));
      diff = (diff1 < diff2) ? diff1 : diff2;
      sum += diff;
    }
  sum /= (xmax - xmin);
  return sqrt(sum);
}

double maxdevsearch(int xmin, int xmax, TH1D* bragghist3, TH1D* bragghist4, TH1D* mfc2)
//finds the highest distance from mfc2 to either of the bragghists.                                                       
{
  double maxdev = 0;
  double valmfc;
  double val3;
  double val4;
  double dist3;
  double dist4;
  double dist;
  double val;
  for (int x = xmin; x < xmax; x++)
    {
      valmfc = mfc2->GetBinContent(x);
      val3 = bragghist3->GetBinContent(x);
      val4 = bragghist4->GetBinContent(x);
      dist3 = fabs(valmfc - val3);
      dist4 = fabs(valmfc - val4);
      dist = (dist3 < dist4) ? dist3 : dist4;
      val =  (dist3 < dist4) ? val3  : val4;
      if (maxdev < dist / val)
	{
	  maxdev = dist / val;
	  maxdevloc = x;
	}
    }
  return maxdev;
}

int terminatorsearch(TH1D* bragghist3, TH1D* bragghist4, TH1D* mfc2, int xmin, int xmax, double threshold, vector<int>& vec, vector<double>& vec2)
//An algorithm sent from the future. Extremely adept at finding the vertex but has flaws that prevent it from achieving its goal on certain types of waveforms.                                                                                      
{
  int thisval = 0;
  int lastval = 0;
  double val3;
  double val4;
  double valmfc;
  double a3 = 0.5;
  double a4 = 1 - a3;
  double thisfit = -2.;
  double der3;
  double der4;
  double dermfc;
  bestfit = -1e100;
  int bestfitloc = 0;
  for (int x = xmin; x < xmax; x++)
    {
      der3 = deriv(bragghist3,x,2,dt);
      der4 = deriv(bragghist4,x,2,dt);
      dermfc = deriv(mfc2,x,2,dt);
      val3 = bragghist3->GetBinContent(x);
      val4 = bragghist4->GetBinContent(x);
      valmfc = mfc2->GetBinContent(x);
      if (val4 == 0) val4 = 100*val3;
      lastval = thisval;

      if (fabs(1 - val3/val4) < threshold)
	{
	  thisval = 0;
	  lastval = 0;
	}

  else
    {
      thisval = (valmfc > (a3 * val3 + a4 * val4 - termderivfactor * (der4 - der3))) ? 1 : -1;
      if (thisval * lastval == -1)
	{
	  thisfit = fabs(val3 - val4) - dermfcfac2 * pow(fabs((val4- val3) + dermfcfac3 * (der4 + der3) - dermfcfac * dermfc),1.5);
	  vec.push_back(x);
	  vec2.push_back(thisfit * 2 / (val3 + val4));
	  if (thisfit > bestfit)
	    {
	      bestfit = thisfit;
	      bestfitloc = x;
	      termdist = 2 * bestfit / (val3 + val4);
	    }
	}
    }
}
if (bestfitloc == 0) bestfitloc = terminatorsearchbackup(bragghist3,bragghist4,mfc2,xmin,xmax);
return bestfitloc;

}

int terminatorsearchbackup(TH1D* bragghist3, TH1D* bragghist4, TH1D* mfc2, int xmin, int xmax)
//Returns an answer in the event that terminatorsearch does not return one.                                               
{
  int bestvalloc = 0;
  double bestval = 1e100;
  double val3;
  double val4;
  double valmfc;
  double valav;
  double valdiff;
  double diff;
  for (int x = xmin; x < xmax; x++)
    {
      val3 = bragghist3->GetBinContent(x);
      val4 = bragghist4->GetBinContent(x);
      valmfc = mfc2->GetBinContent(x);
      valav = (val3 + val4) / 2;
      valdiff = fabs(val3 - val4);
      if (valdiff == 0.) valdiff = 1e-100;
      diff = (fabs(valav - valmfc)+TSB_coeff*pow(valdiff,TSB_pow/2)) / pow(valdiff,TSB_pow);
      if (diff < bestval)
	{
	  bestval = diff;
	  bestvalloc = x;
	}
    }
  return bestvalloc;
}

int terminatorsearch2(TH1D* bragghist3, TH1D* bragghist4, TH1D* mfc2, int xmin, int xmax, double threshold)
//Modified version of the Terminator. Generally gives poor results, but tends to give decent results precisely when Terminator gives bad ones.
{
  double bestval = 1e100;
  int bestvalloc = 0;
  double val3;
  double val4;
  double valav;
  double thisval;
  double valmfc;
  double difference;
  for (int x = xmin; x < xmax; x++)
    {
      val3 = bragghist3->GetBinContent(x);
      val4 = bragghist4->GetBinContent(x);
      valav = (val3 + val4) / 2.;
      valmfc = mfc2->GetBinContent(x);
      if (val4 == 0) val4 = 100*val3;
      difference = fabs(val3 - val4) + .1*threshold*val3;
      if (difference == 0) difference = 1e-100;
      if (fabs(1 - val3/val4) < threshold)
	{
	  thisval = fabs(valmfc - valav) / pow(difference,term2pow);
	  if (thisval < bestval)
	    {
	      bestval = thisval;
	      bestvalloc = x;
	    }}
    }
  return bestvalloc;
}

int terminatorsearch3(TH1D* bragghist3, TH1D* bragghist4, int xmin, int xmax)
//A crude but consistent version of Terminator2. Not recommended on its own. May be useful in combination with other algorithms
{
  int crosspoint = 0;
  int thisval = 0;
  int lastval = 0;
  for (int x = xmin; x < xmax; x++)
    {
      lastval = thisval;
      thisval = (bragghist3->GetBinContent(x) > bragghist4->GetBinContent(x)) ? 1 : -1;
      if (thisval * lastval == -1)
	crosspoint = x;
    }
  return crosspoint;
}

void plotbragg(TH1D* bragghist3, TH1D* bragghist2, double xscaletrue, double xscale, double yscale, int xoff)
//Scales and shifts a bragg curve to match the observed event.
{
  for (int x = 1; x < 16384 - xoff; x++)
    {
      bragghist3->SetBinContent(x-1 + xoff,histspline(bragghist2,0,0,x,1/xscaletrue,yscale / xscale));
      //if (x%10 == 0) cout << x << '|' << bragghist3->GetBinContent(x) << endl;
    }
}


void plotbraggrev(TH1D* bragghist4, TH1D* bragghist2, double xscaletrue, double xscale, double yscale, int xoff)
//For bragg curves that go left-to-right instead of right-to-left.                                                        
{
  for (int x = xoff; x > 0; x--)
    bragghist4->SetBinContent(x-1,histspline(bragghist2,0,0,xoff - x + 1,1/xscaletrue,yscale / xscale));
}
double integrate(TH1D* hist, int imin, int imax)
//Gets area under hist.                                                                                                   
{
  double sum = 0;
  for (int i = imin+1; i < imax; i++)
    sum += hist->GetBinContent(i);
  return sum;
}

int maxsearch(TH1D* hist, int imin, int imax)
//Finds max bin of hist within given limits.                                                                              
{
  int maxloc = 0;
  double maxval = -1e100;
  for (int i = imin; i < imax; i++)
    if (hist->GetBinContent(i) > maxval)
      {
	maxloc = i;
	maxval = hist->GetBinContent(i);
      }
  return maxloc;
}

int halfsearch(TH1D* hist, int imin, int imax, double maxval)
//Finds where a hist reaches half of its maxval within given limits.                                                      
//Note: Currently programmed to actually look for 75% height.                                                             
{
  int lastcheck = 0;
  int thischeck = 0;
  for (int i = imin; i < imax; i++)
    {
      lastcheck = thischeck;
      thischeck = (hist->GetBinContent(i) < .75 * maxval) ? -1 : 1;
      if (thischeck * lastcheck == -1)
	{
	  return i;
	}
    }
  cout << "error in function 'halfsearch'. no value found. " << imin << endl;

  return 0;
}


int absmin(TH1D* hist, int i, int range)
//Finds max bin of hist within range bins of i.                                                                           
{
  int minloc = 0;
  double minval = 1e100;
  for (int j = i - range; j < i + range; j++)
    if (j>0 && j<16384)
      if (hist->GetBinContent(j) < minval)
	{
	  minval = hist->GetBinContent(j);
	  minloc = j;
	}
  return minloc;
}
double deriv(TH1D* hist, int n, int k, double h)
//returns first derivative of step size k at bin n with hist step size h.                                                 
{
  if (deriv_order == 2) return (hist->GetBinContent(n+k) - hist->GetBinContent(n-k)) / (2 * h * k);
  else if (deriv_order == 4) return (8*hist->GetBinContent(n+k)-8*hist->GetBinContent(n-k)-hist->GetBinContent(n+2*k)+hist->GetBinContent(n-2*k))/(12*h*k);
  else
    {
      deriv_order = 2;
      return deriv(hist,n,k,h);
    }
}void deconvolve(TH1D* hist, double SD, int order)
 //Un-does Gaussian smearing. May contain traces of sorcery.
{
  TH1D* hist2 = new TH1D("asdf10","jkl10",16384,0,16384);
  for (int i = 2000; i < 10000; i++)
    {
      if (order == 2)
	hist2->SetBinContent(i,hist->GetBinContent(i) - 0.25 * pow(SD*dt,2) * deriv2(hist,i,d2_step,dt));
      else if (order == 4)
	hist2->SetBinContent(i,hist->GetBinContent(i) - 0.25 * pow(SD*dt,2) * deriv2(hist,i,d2_step,dt) + .125 * pow(SD*dt,4) * deriv4(hist,i,d4_step,dt));
      else if (order == 6)
	  hist2->SetBinContent(i,hist->GetBinContent(i) - 0.25 * pow(SD*dt,2) * deriv2(hist,i,d2_step,dt) + .125 * pow(SD*dt,4) * deriv4(hist,i,d4_step,dt) - 0.09375 * pow(SD*dt,6) * deriv6(hist,i,d6_step,dt));
      else if (order == 8)
	hist2->SetBinContent(i,hist->GetBinContent(i) - 0.25 * pow(SD*dt,2) * deriv2(hist,i,d2_step,dt) + .125 * pow(SD*dt,4) * deriv4(hist,i,d4_step,dt) - .09375 * pow(SD*dt,6)* deriv6(hist,i,d6_step,dt) + .09375 * pow(SD*dt,8) * deriv8(hist,i,d8_step,dt));
      else
	cout << "invalid value for variable 'order' (must be 2, 4, 6, or 8.)" << endl;
    }
      
  for (int i = 0; i < 16384; i++)
    hist->SetBinContent(i,0);
  for (int i = 2100; i < 10000; i++)
    hist->SetBinContent(i,hist2->GetBinContent(i));
}

double deriv2(TH1D* hist, int n, int k, double h)
//Second derivative.
{
  return (hist->GetBinContent(n+k) + hist->GetBinContent(n-k) - 2*hist->GetBinContent(n))/(h*k*h*k);
}
double deriv4(TH1D* hist, int n, int k, double h)
//Fourth derivative.   
{
  return (hist->GetBinContent(n+2*k) - 4* hist->GetBinContent(n+k) + hist->GetBinContent(n-2*k)-4* hist->GetBinContent(n-k) +6*hist->GetBinContent(n))/(pow(h*k,4));
}

double deriv6(TH1D* hist, int n, int k, double h)
{
  return (hist->GetBinContent(n+3*k)-6*hist->GetBinContent(n+2*k)+15*hist->GetBinContent(n+k)-20*hist->GetBinContent(n)+15*hist->GetBinContent(n-k)-6*hist->GetBinContent(n-2*k)+hist->GetBinContent(n-3*k))/(pow(h*k,6));
}
double deriv8(TH1D* hist, int n, int k, double h)
{
  return (hist->GetBinContent(n+4*k) - 8 * hist->GetBinContent(n + 3*k) + 28 * hist->GetBinContent(n + 2*k) - 56 * hist->GetBinContent(n + k) + 70 * hist->GetBinContent(n) - 56 * hist->GetBinContent(n - k) + 28 * hist->GetBinContent(n-2*k) - 8 * hist->GetBinContent(n-3*k) + hist->GetBinContent(n-4*k))/(pow(h*k,8));
}

int minfinder(TH1D* wfder)
{
  double thresh = .2 * wfder->GetBinContent(wfder->GetMaximumBin());
  TH1D* dummy = new TH1D("asdf11","jkl11",16384,0,16384);
  for (int i = 0; i < 16384; i++)
    if (wfder->GetBinContent(i) > thresh)
      dummy->SetBinContent(i,1e6);
  return minfinder(dummy,wfder);
}
int minfinder(TH1D* wf, TH1D* wfder)
//Finds the minimum of certain processed versions of the waveform. Provides weak estimates of the vertex.                 
{   
  int fallen = 0;
  int risen = 0;
  int jmin = 0;
  int cooldownmax = 16;
  double coolfactor = .7;
  double threshold = 1.0016 + .11/rec_SD;

  for(int j = 2000;j < 10000; j++)
    {
      if (threshold * wfder->GetBinContent(j) < wfder->GetBinContent(j-cooldownmax) && wfder->GetBinContent(j) > 1e7 && wf->GetBinContent(j) > 30){fallen = 1;}
      if (wfder->GetBinContent(j) > threshold * wfder->GetBinContent(j-cooldownmax) && fallen == 1 && risen == 0)
	{risen = 1;
	  jmin = int(j - coolfactor * cooldownmax + .5);}
    }
  return jmin;
}
void minprepare(TH1D* wf, TH1D* minfinder_curve)
//processes a waveform into a format for minfinder.                                                                       
{
  double k = 1/Td;

  TH1D* wfder = new TH1D("asdf12","jkl12",16384,0,16384);

  for(int j = 100;j < 16300; j++)
    wfder->SetBinContent(j-1,k*(wf->GetBinContent(j)-vertoff)+deriv(wf,j,deriv_step,dt));
  //small offset accounts for rise time                                                                                   

  double wfd_max = wfder->GetBinContent(wfder->GetMaximumBin());
  double wfd_thresh = 0.3 * wfd_max;
  int wfd_start = 0;
  int wfd_end = 0;
  int stop = 0;
  int countdownmax = 50;
  int countdown = countdownmax;
  for(int j=2200;j<11800 && stop == 0;j++)
    {
      if (wfder->GetBinContent(j)>wfd_thresh)
	{
	  wfd_end = j;
	  countdown = countdownmax;
	}
      else
	if (wfd_start !=0)
	  {
	    countdown -=1;
	    if (countdown == 0) stop = 1;
	  }
      if (wfd_start == 0 && wfd_end > 0) wfd_start = wfd_end;
    }

  wfd_delta = wfd_end - wfd_start;
  termderivfactor = double(wfd_delta) / 1000. * -1.e-8;
  double smoothsize = double(wfd_delta)/26 + 4;
  double correction = 1;
  smoothsize *= min_mult;
  correction *= min_cmult;
  if (min_mult > 0)
    {
      smooth(wfder,smoothsize,gaus_acc);
      correct(wfder,min_cwidth*smoothsize,gaus_acc,correction);
    }
  for(int j=0;j<16384;j++)
    minfinder_curve->SetBinContent(j,wfder->GetBinContent(j));
  delete wfder;
}

double average(int j, vector<double>& vec, TH1D* hist)
//Averages points on a hist according to weights in vec.                                                                  
{
  int sz = (vec.size() - 1) / 2;
  double sum = 0;
  double coeffsum = 0;
  for (int i = -sz; i < sz+1; i++)
    if (i+sz>=0 && i+sz<vec.size()&&i+j>0&&i+j<16384)
      {
	coeffsum += vec.at(i+sz);
	sum += vec.at(i+sz) * hist->GetBinContent(i+j);
      }
  if (coeffsum == 0) coeffsum = 1;
  sum /= coeffsum;
  return sum;
}

double median(vector<double>& vec)
//Returns median entry of a vector. 
{
  if (vec.size()==0){return 0;}
  sort(vec.begin(),vec.end());
  return vec.size()%2 == 1 ? vec.at((vec.size()-1) / 2) : .5*(vec.at(vec.size()/2)+vec.at(vec.size()/2-1));

}

double normaldist(double peak, double x, double SD)
//evaluates a normal distribution at x with given SD and peak location.                                                   
{return exp((x-peak)*(peak-x)/(2*SD*SD)) / (SD * pow(3.141592,0.5));}

void normalize(vector<double>&vec)
//Does what it says.
{
  double sum = 0;
  for (int i = 0; i < vec.size(); i++)
    sum += vec.at(i);
  if (sum == 0) sum = 1;
  for (int i = 0; i < vec.size(); i++)
    vec.at(i) /= sum;
}

void makekernel(vector<double>&vec, double SD, double acc)
//Fills vec with a Gaussian distribution for use with average(). 
{
  if (SD == 0) SD = 1.e-69;
  int size = int(acc*SD+1);
  for (int i = -size; i < size + 1; i++)
    vec.push_back((4*normaldist(0.,double(i),SD)+normaldist(0.,double(i)+.5,SD)+normaldist(0.,double(i)-.5,SD))/6.);
  normalize(vec);
}
void smooth(TH1D* hist, double SD, double acc)
//Applies Gaussian blur to a histogram. Used as a low-pass filter. 
{
  TH1D* hist2 = new TH1D("asdf13","jkl13",16384,0,16384);
  vector<double> vec;
  makekernel(vec,SD,acc);

  for(int j=2000; j < 10200; j++)
    hist2->SetBinContent(j,average(j,vec,hist));
  for(int j=1; j < 16384; j++)
    hist->SetBinContent(j,hist2->GetBinContent(j));
  delete hist2;
}

void correct(TH1D* hist, double SD, double acc, double factor)
//Un-does some of the undesirable effects of smooth(). 
{
  TH1D* hist2 = new TH1D("asdf14","jkl14",16384,0,16384);
  vector<double> vec;
  makekernel(vec,SD,acc);
  for(int j=2000; j < 10200; j++)
    hist2->SetBinContent(j,(1+factor)*hist->GetBinContent(j)-factor*average(j,vec,hist));
  for(int j=1; j < 16384; j++)
    hist->SetBinContent(j,hist2->GetBinContent(j));
  delete hist2;
}

void shifthist(TH1D* hist, int plot_offset)
//Shifts the values in a histogram on the x-axis.                                                                         
{
  for (int i = 16383; i > 0; i--)
    if (i > plot_offset)
      hist->SetBinContent(i,hist->GetBinContent(i+plot_offset));
    else
      hist->SetBinContent(i,0);
}

double histspline(TH1D*hist, int xmin, int xoff, int x, double xscale, double yscale)
//An interface to read decimal values out of a histogram as if it were a spline curve.                                    
{
  double true_x = xscale * (x - xmin) + xoff;
  int int_x = int(true_x);
  double weight1 = fabs(true_x - int_x);
  double weight0 = 1 - weight1;
  if (int_x < 0) int_x = -int_x;
  double value =  weight0 * hist->GetBinContent(int_x) + weight1 * hist->GetBinContent(int_x + 1);
  return value * yscale;
}

// End Big-Ugly-Double-Alpha-Analysis-Definitions

double
analysis::baseline(const TH1* hist, double& rms, int binMin, int binMax)
{

  double base=0;
  rms=0;
  int n = binMax-binMin+1;
  for (int i = binMin; i<=binMax; i++)
  {
    double v = hget(hist,i);
    base+=v;
    rms+=v*v;
  }

  base/=n;
  rms/=n;
  rms -= base*base;
  rms = sqrt(rms);
  return base;

}

int
analysis::minbin(const TH1F* hist)
{
  int minbin=0;
  // neglects the underflow, the overflow, and the last bin
  minbin=TMath::LocMin(hist->GetNbinsX()-1,hist->GetArray()+1)+1;
  
  return minbin;
} 

int
analysis::maxbin(const TH1F* hist)
{
  int maxbin=0;
  // neglects the underflow, the overflow, and the last bin
  maxbin=TMath::LocMax(hist->GetNbinsX()-1,hist->GetArray()+1)+1;

  return maxbin;
}


bool 
 analysis::isPeak(const TH1F* hist, int bin, int nbins)
{
  double binval =hget(hist,bin);
  int minBin = max(1,bin-nbins);
  int maxBin = min(hist->GetNbinsX(),bin+nbins);
  for (int i = minBin; i<= maxBin; i++)
    if (hget(hist,i)>binval) return false;

  return true;
}

bool 
analysis::isTrough(const TH1F* hist, int bin, int nbins)
{
  double binval =hget(hist,bin);
  int minBin = max(1,bin-nbins);
  int maxBin = min(hist->GetNbinsX(),bin+nbins);
  for (int i = minBin; i<= maxBin; i++)
    if (hget(hist,i)<binval) return false;

  return true;
}


double
analysis::integral(const TH1F* hist, double start, double end)
{


  start = max(hist->GetXaxis()->GetXmin(),start);
  end = min(hist->GetXaxis()->GetXmax(),end);
  int startBin = hist->GetXaxis()->FindBin(start);
  int endBin = hist->GetXaxis()->FindBin(end);

  double startCenter = hist->GetXaxis()->GetBinCenter(startBin);
  double endCenter = hist->GetXaxis()->GetBinCenter(endBin);
  
  double width = hist->GetXaxis()->GetBinWidth(1);

  double sum = 0;
  sum += hget(hist,startBin) * (start-startCenter+0.5*width)/width;
  sum += hget(hist,endBin) * (endCenter - 0.5*width-end)/width;
  for(int i = startBin+1; i<endBin; i++) sum += hget(hist,i);
  return sum;
}



void 
analysis::riseTime(const TH1F* hist, const vector<double>& list,
                  vector<double>& values, double startTime, double endTime,
                  bool fromStart)
{

  double binWidth = hist->GetXaxis()->GetBinWidth(1);
  int n = list.size();

  values.clear();
  values.resize(n,-1);
  
  int startBin = hist->GetXaxis()->FindBin(startTime);
  int endBin = hist->GetXaxis()->FindBin(endTime);
  startBin = startBin==1?startBin:startBin-1;
  double v;
  if (fromStart){
    double lastV=hget(hist,startBin);
    for (int i = startBin; i<=endBin ;i++)
    {
      v = hget(hist, i);
    
      for (int j = 0; j<n ; j++){
            
        if (values[j]==-1&&lastV<list[j]&&v>=list[j]){
          double binFrac = (list[j] - lastV) / (v-lastV);
          values[j] = binFrac*binWidth + hcenter(hist,i-1);

          if (j==n-1) break;

        }//if reaches value
 

      }//values
      lastV=v;
    }//bins
  }else{
    double lastV=hget(hist,endBin);
    for (int i = endBin; i>=startBin; i--)
    {
       v = hget(hist,i);
       for (int j=n-1; j>=0; j++){

         if (values[j]==-1&&lastV>list[j]&&v<=list[j]){

           double binFrac = (list[j]-lastV)/(v-lastV);
           values[j] = hcenter(hist,i+1)-binFrac*binWidth;

           if (j==0) break;
         }


       }//values
       lastV=v;
    }
  }//which side to start at
  

}

void 
analysis::fallTime(const TH1F* hist,const vector<double>& list,
                  vector<double>& values, double startTime, double endTime,
                  bool fromStart)
{

  double binWidth = hist->GetXaxis()->GetBinWidth(1);
  int n = list.size();
  values.clear();
  values.resize(n,-1);
  int startBin = hist->GetXaxis()->FindBin(startTime);
  int endBin = hist->GetXaxis()->FindBin(endTime);
  endBin = endBin==hist->GetNbinsX()?endBin:endBin+1;
  double v;
  if (fromStart){
    double lastV=hget(hist,startBin);
    for (int i = startBin; i<=endBin ;i++)
    {
      v = hget(hist, i);
      for (int j = 0; j<n ; j++){

        if (values[j]==-1&&lastV>list[j]&&v<=list[j]){
          double binFrac = (list[j] - lastV) / (v-lastV);
          values[j] = binFrac*binWidth + hcenter(hist,i-1);

          if (j==n-1) break;

        }//if reaches value


      }//values
      lastV=v;
    }//bins
  }else{
    double lastV=hget(hist,endBin);
    for (int i = endBin; i>=startBin; i--)
    {
       v = hget(hist,i);
       for (int j=n-1; j>=0; j--){

         if (values[j]==-1&&lastV<list[j]&&v>=list[j]){

           double binFrac = (list[j]-lastV)/(v-lastV);
           values[j] = hcenter(hist,i+1)-binFrac*binWidth;

           if (j==0) break;
         }


       }//values
       lastV=v;
    }
  }//which side to start at


}

double 
analysis::startTime(const TH1F* hist, int bin,
                     double threshold,int& startBin,int minBin)
{

  if (minBin==-1) minBin=1;

  double binWidth = hist->GetXaxis()->GetBinWidth(1);
  startBin = minBin;
  double time = hcenter(hist,minBin);
  double v, lastV = hget(hist,bin);
  for (int i = bin; i>= minBin; i--)
  {
    v = hget(hist,i);
    if (lastV > threshold&&v<=threshold){
      double binFrac = (threshold-lastV)/(v-lastV);
      time = hcenter(hist,i+1)-binWidth*binFrac;
      //      time = hist->GetBinLowEdge(i+1);
      startBin = binFrac>0.5?i:i+1;
      //      startBin = i+1;
      break;
    }
    lastV = v;
  }
  return time;
}

double 
analysis::endTime(const TH1F* hist, int bin,
               double threshold,int& endBin,int maxBin)
{

  if (maxBin==-1) maxBin=hist->GetNbinsX();

  double binWidth = hist->GetXaxis()->GetBinWidth(1);
  endBin = maxBin;
  double time = hcenter(hist,maxBin);
  double v, lastV = hget(hist,bin);
  for (int i = bin; i<= maxBin; i++)
  {
    v = hget(hist,i);
    if (lastV > threshold&&v<=threshold){
      double binFrac = (threshold-lastV)/(v-lastV);
      time = hcenter(hist,i-1)+binWidth*binFrac;
      //      time = hist->GetBinLowEdge(i);
      endBin = binFrac>0.5?i:i-1;
      break;
    }
    lastV = v;
  }
  return time;

}

vector<int> 
analysis::peaks(const TH1F* hist, double threshold,
                          int minBin,int maxBin)
{

  minBin = minBin<0?1:minBin;
  maxBin = maxBin<0?hist->GetNbinsX():maxBin;
  vector<int> p;
  double v, lastV=hget(hist,minBin);
  bool isInPeak = false;
  for (int i = minBin; i<= maxBin; i++)
  {

    v = hget(hist,i);
    if (lastV<threshold&&v>=threshold){
      isInPeak=true;
      p.push_back(i);
    }else if (lastV>= threshold&&v<threshold) isInPeak=false;

    if (isInPeak){
      if (v>hget(hist,p[p.size()-1])) p[p.size()-1] = i;
    }

    lastV=v;
  }
  return p;

}

vector<int> 
analysis::valleys(const TH1F* hist, double threshold,
                            int minBin,int maxBin)
{

  minBin = minBin<0?1:minBin;
  maxBin = maxBin<0?hist->GetNbinsX():maxBin;
  vector<int> p;
  double v, lastV=hget(hist,minBin);
  bool isInPeak = false;
  for (int i = minBin; i<= maxBin; i++)
  {

    v = hget(hist,i);
    if (lastV>threshold&&v<=threshold){
      isInPeak=true;
      p.push_back(i);
    }else if (lastV<= threshold&&v>threshold) isInPeak=false;
    
    if (isInPeak){
      if (v<hget(hist,p[p.size()-1])) p[p.size()-1] = i;
    }
    
    lastV=v;
  }
  return p;

}

void 
analysis::peaksAndValleys(const TH1F* hist,vector<int>& pkBin,
                         vector<int>& valBin, 
                         int minBin,int maxBin,int nbins)
{

  minBin = minBin<0?1:minBin;
  maxBin = maxBin<0?hist->GetNbinsX():maxBin;

  double v;
  double currentMin=0;
  int bin = 0;
  int npeak = 0;
  for (int i = minBin; i<=maxBin; i++)
  {
    v = hget(hist,i);
    if (isPeak(hist,i,nbins))
    {

      pkBin.push_back(i);
      if (npeak>0) valBin.push_back(bin);
      currentMin = v;
      bin = i;
      npeak++;
    }else{

      bin = v<currentMin?i:bin;
      currentMin = v<currentMin?v:currentMin;
      
    }

  }


}

void 
analysis::mergePeaksByDistance(const TH1F* hist,vector<int>& pkBin, 
                               vector<int>& valBin,int minDist)
{

  if(pkBin.size()<=1) return;

  //First remove all peaks within a certain distance of one another
  vector<int> rmPk(pkBin.size(),0);
  //Find all peaks within this distance from one another
  for (unsigned int i = 1; i<pkBin.size(); i++)
    if (pkBin[i]-pkBin[i-1] <minDist) rmPk[i] = 1;

  //Loop through everything and assign peaks new bins
  for (unsigned int i = 1; i<pkBin.size(); i++)
  {
    if (rmPk[i]==0) continue;
    int bin = hget(hist,pkBin[i])>hget(hist,pkBin[i-1]) ?
              pkBin[i] : pkBin[i-1];

    unsigned int j = i-1;
    while(rmPk[j+1]==1)
    {
      pkBin[j] = bin;
      j--;
    }
  }

  //Eliminate repeated peaks

  for (unsigned int i = 1; i<pkBin.size(); i++)
  {
    if (pkBin[i]==pkBin[i-1]){
      verase(pkBin,i);
      i--;
    }
  }

  //Assign each valley to a possible position

  int minP = pkBin[0], maxP = pkBin[pkBin.size()-1];
  vector<int> valPos;
  for (unsigned int i = 0; i <valBin.size(); i++)
  {
    if (valBin[i]>=maxP||valBin[i]<minP){
       verase(valBin,i);
       i--;
       continue;
    }

    for (unsigned int j = 0; j < pkBin.size()-1; j++)
    {
      if (valBin[i] >= pkBin[j]&&valBin[i]<pkBin[j+1])
      {
        valPos.push_back(j);
        break;
      }

    }
  }

  //cout <<valBin.size()<<endl;                                                            
  if(valBin.size()>1){
    for (unsigned int i = 0; i< valBin.size()-1; i++)
      {
	
	if (valPos[i]!=valPos[i+1]) continue;
	if (hget(hist,valBin[i]) <= hget(hist,valBin[i+1])){
	  verase(valPos,i+1);
	  verase(valBin,i+1);
	}else{
	  verase(valPos,i);
	  verase(valBin,i);
	}
	i--;
      }
  }

}


/*

void 
analysis::mergePeaksByDepth(const TH1F* hist,vector<int>& pkBin, 
                               vector<int>& valBin,double minDepth)
{

  if (pkBin.size()<=1) return;

  //First bin:
  if (hget(hist,pkBin[0])-hget(hist,valBin[0])<minDepth ||
      hget(hist,pkBin[1])-hget(hist,valBin[0])<minDepth )
  {

    if (hget(hist,pkBin[0])>hget(hist,pkBin[1]))
    {//Second peak is spurious
      verase(pkBin,1);
      if (valBin.size()==1) verase(valBin,0);
      else{
        if (hget(hist,valBin[0])>hget(hist,valBin[1]))
          verase(valBin,0);
        else verase(valBin,1);
      }
    }else{//First peak is spurious
      verase(pkBin,0);
      verase(valBin,0);
    }
  


  }




  //Remove all very shallow peaks (likely noise)
  vector<int> rmPk(pkBin.size(),0);
  //Find all peaks to remove 
  rmPk[0] = hget(hist,pkBin[i]-hget(hist,valBin[i])) < minDepth;
  rmPk[pkBin.size()-1] = hget(hist,pkBin[i]) - hget(hist,valBin[i-1]) < minDepth;

  for (unsigned int i = 1; i<pkBin.size()-1; i++)
    if (hget(hist,pkBin[i])-hget(hist,valBin[i-1]) <minDepth ||
        hget(hist,pkBin[i])-hget(hist,valBin[i]) <minDepth)
    rmPk[i] = 1;

  //Loop through everything and assign peaks new bins
  for (unsigned int i = 1; i<pkBin.size(); i++)
  {
    if (rmPk[i]==0) continue;
    int bin = hget(hist,pkBin[i])>hget(hist,pkBin[i-1]) ?
              pkBin[i] : pkBin[i-1];

    unsigned int j = i-1;
    while(rmPk[j+1]==1)
    {
      pkBin[j] = bin;
      j--;
    }
  }
  //Eliminate repeated peaks
  for (unsigned int i = 1; i<pkBin.size(); i++)
  {
    if (pkBin[i]==pkBin[i-1]){
      verase(pkBin,i);
      i--;
    }
  }

  //Assign each valley to a possible position
  int minP = pkBin[0], maxP = pkBin[pkBin.size()-1];
  vector<int> valPos;
  for (unsigned int i = 0; i <valBin.size(); i++)
  {
    if (valBin[i]>=maxP||valBin[i]<minP){
       verase(valBin,i);
       i--;
       continue;
    }

    for (unsigned int j = 0; j < pkBin.size()-1; j++)
    {
      if (valBin[i] >= pkBin[j]&&valBin[i]<pkBin[j+1])
      {
        valPos.push_back(j);
        break;
      }
    }
  }
  //Find best match for each position
  for (unsigned int i = 0; i< valBin.size()-1; i++)
  {
    if (valPos[i]!=valPos[i+1]) continue;
    if (hget(hist,valBin[i]) <= hget(hist,valBin[i+1])){
      verase(valPos,i+1);
      verase(valBin,i+1);
    }else{
      verase(valPos,i);
      verase(valBin,i);
    }
    i--;
  }

      
}
*/

void 
analysis::analyzeCSP(const TH1F* h, CspWaveform& wf, Double_t gausConvSigma)
{
  TH1F*  htemp;

  double base,
         rms,
         hmax,
	 hmaxTime,
	 hmin,
	 hminTime;
  int    hmaxBin,
         hminBin;

  vector<double> rise(6),
                 fall(6);

  double peak,
	 peakTime,
	 startTime,
	 endTime,
	 integral;
  int    peakBin,
         startBin,
	 endBin;

  double threshold;
  vector<double> riseV(6),
                 fallV(6);

  double rf[6] = {0,0.1,0.25,0.5,0.75,0.9};
  double binWidth = h->GetXaxis()->GetBinWidth(1);

//Clear waveform
  wf.clear();

//Smoothing
  htemp = tools::gausConv(h,gausConvSigma);
  
//General waveform parameters
  base = analysis::baseline(htemp,rms,1,500); 
  hmax = htemp->GetMaximum();
  hmaxBin = htemp->GetMaximumBin();
  hmin = htemp->GetMinimum();
  hminBin = htemp->GetMinimumBin();
  hmaxTime = htemp->GetXaxis()->GetBinCenter(hmaxBin);
  hminTime = htemp->GetXaxis()->GetBinCenter(hminBin);

  wf.setBase(base);
  wf.setRMS(rms);
  wf.setWfMax(hmax);
  wf.setWfMaxTime(hmaxTime);
  wf.setWfMaxBin(hmaxBin);
  wf.setWfMin(hmin);
  wf.setWfMinTime(hminTime);
  wf.setWfMinBin(hminBin);

//Pulse parameters: only allow 1 pulse per waveform for CSP

  peak = hmax-base;
  peakBin = hmaxBin;
  peakTime = hmaxTime;

  threshold = base;
  startTime = analysis::startTime(htemp,peakBin,threshold,startBin);
  endTime = analysis::endTime(htemp,peakBin,threshold,endBin);

  for (int i = 0; i < 6; i++)
  {
    riseV[i] = rf[i] * peak + base;
    fallV[i] = rf[5-i] * peak + base;
  }
  
  analysis::riseTime(htemp,riseV,rise,max(startTime-binWidth,htemp->GetXaxis()->GetXmin()),peakTime,true);
  analysis::fallTime(htemp,fallV,fall,peakTime,min(endTime+binWidth,htemp->GetXaxis()->GetXmax()-binWidth*0.001),false);

  for (int i =0; i<6; i++)
  {
    rise[i] = peakTime - rise[i];
    fall[i] = fall[i]-peakTime;
  }

  //rise[0] = peakTime-startTime;
  integral = analysis::integral(htemp,startTime,endTime);
  CspPulse p(peakBin);
  p.setRise( &(rise[0]) );
  p.setFall( &(fall[0]) );
  p.setPeak(peak);
  p.setPeakTime(peakTime);
  p.setStartTime(startTime);
  p.setEndTime(endTime);
  p.setStartBin(startBin);
  p.setEndBin(endBin);
  p.setIntegral(integral);
  wf.add(p);

  delete htemp;
}

void 
analysis::analyzePMT(const TH1F* h, PMTWaveform& wf, Double_t gausConvSigma)
{
  TH1F*  htemp;  TH1F*  htempinverted;

  double base,
         rms,
         hmax,
	 hmaxTime,
	 hmin,
	 hminTime;
  int    hmaxBin,
         hminBin;

  vector<double> rise(6),
                 fall(6);

  double peak,
	 peakTime,
	 startTime,
	 endTime,
	 integral;
  int    peakBin,
         startBin,
	 endBin;

  double threshold;
  vector<double> riseV(6),
                 fallV(6);

  double rf[6] = {0,0.1,0.25,0.5,0.75,0.9};
  double binWidth = h->GetXaxis()->GetBinWidth(1);

//Clear waveform
  wf.clear();

//Smoothing
  htemp = ((TH1F*)h->Clone()); 
//  htemp = tools::gausConv(h,gausConvSigma);

//Inversion
  htempinverted = ((TH1F*)htemp->Clone());
  for (int i=1; i<=htempinverted->GetNbinsX(); i++) htempinverted->SetBinContent( i, -htempinverted->GetBinContent(i) );

//General waveform parameters
  base = analysis::baseline(htemp,rms,1,500); 

  hmax = htemp->GetMaximum();
  hmaxBin = htemp->GetMaximumBin();
  hmin = htemp->GetMinimum();
  hminBin = htemp->GetMinimumBin();
  hmaxTime = htemp->GetXaxis()->GetBinCenter(hmaxBin);
  hminTime = htemp->GetXaxis()->GetBinCenter(hminBin);

  wf.setBase(base);
  wf.setRMS(rms);
  wf.setWfMax(hmax);
  wf.setWfMaxTime(hmaxTime);
  wf.setWfMaxBin(hmaxBin);
  wf.setWfMin(hmin);
  wf.setWfMinTime(hminTime);
  wf.setWfMinBin(hminBin);

//Pulse parameters: only allow 1 pulse per waveform for CSP

  peak = hmin-base;
  peakBin = hminBin;
  peakTime = hminTime;

  threshold = base;

  startTime = analysis::startTime(htempinverted,peakBin,-threshold,startBin);
  endTime = analysis::endTime(htempinverted,peakBin,-threshold,endBin);

  for (int i = 0; i < 6; i++)
  {
    riseV[i] = rf[i] * peak + base;
    fallV[i] = rf[5-i] * peak + base;
  }
  
  analysis::riseTime(htemp,riseV,rise,max(startTime-binWidth,htemp->GetXaxis()->GetXmin()),peakTime,true);
  analysis::fallTime(htemp,fallV,fall,peakTime,min(endTime+binWidth,htemp->GetXaxis()->GetXmax()-binWidth*0.001),false);

  for (int i =0; i<6; i++)
  {
    rise[i] = peakTime - rise[i];
    fall[i] = fall[i]-peakTime;
  }

  //rise[0] = peakTime-startTime;
  integral = analysis::integral(htemp,startTime,endTime);
  PMTPulse p(peakBin);
  p.setRise( &(rise[0]) );
  p.setFall( &(fall[0]) );
  p.setPeak(peak);
  p.setPeakTime(peakTime);
  p.setStartTime(startTime);
  p.setEndTime(endTime);
  p.setStartBin(startBin);
  p.setEndBin(endBin);
  p.setIntegral(integral);
  wf.add(p);

  delete htemp;
  delete htempinverted;
}

void
analysis::analyzeFast(const TH1F* h, FastWaveform& wf, int runnum)
{

  knobfile = "runParameters/LittleDCTPC_far/WfKnobs.temp";
  TurnKnobs(runnum);
  analyze();

  

  double SMOOTH_WIDTH=1.5;
  double MIN_FAST_PEAK=0.65;
  int PEAK_FINDER_MIN_DIST=5;
  int MIN_PEAK_DIST=10;

  TH1F*  htemp;
  double base,
         rms,
         hmax,
	 hmaxTime,
	 hmin,
	 hminTime;
  int    hmaxBin,
         hminBin;

  vector<double> rise(6),
                 fall(6);

  double peak,
	 peakTime,
	 startTime,
	 endTime,
	 integral;
  int    peakBin,
         startBin,
	 endBin;

  double fastPeak=0,
	 fastPeakTime=0,
         trough=0,
	 troughTime=0,
	 slowPeak=0,
	 slowPeakTime=0;

  int	 fastPeakBin=0,
	 troughBin=0,
	 slowPeakBin=0;

  double threshold;
  vector<double> riseV(6),
                 fallV(6);

  vector<int> pkBin,
              valBin;

  int minBin, maxBin;
  double rf[6] = {0,0.1,0.25,0.5,0.75,0.9};
  double binWidth = h->GetXaxis()->GetBinWidth(1);

//Clear waveform
  wf.clear();

//Smoothing
  htemp = tools::gausConv(h,SMOOTH_WIDTH);

//General waveform parameters
  base = analysis::baseline(htemp,rms,1,500); 
  hmax = htemp->GetMaximum();
  hmaxBin = htemp->GetMaximumBin();
  hmin = htemp->GetMinimum();
  hminBin = htemp->GetMinimumBin();
  hmaxTime = htemp->GetXaxis()->GetBinCenter(hmaxBin);
  hminTime = htemp->GetXaxis()->GetBinCenter(hminBin);

  wf.setBase(base);
  wf.setRMS(rms);
  wf.setWfMax(hmax);
  wf.setWfMaxTime(hmaxTime);
  wf.setWfMaxBin(hmaxBin);
  wf.setWfMin(hmin);
  wf.setWfMinTime(hminTime);
  wf.setWfMinBin(hminBin);
  wf.setLeftInt(leftint);
  wf.setLeftInt2(leftint2);
  wf.setRightInt(rightint);
  wf.setRightInt2(rightint2);
  wf.setjMin(jmin);
  wf.setjBragg(jbragg);
  wf.setjTerm(jterm);
  wf.setjTerm2(jterm2);
  wf.setjTerm3(jterm3);
  wf.setOrigin(origin);
  wf.setRecSD(rec_SD);
  wf.setWfdDelta(wfd_delta);
  wf.setPeak1(peak1);
  wf.setPeak2(peak2);
  wf.setPeak1Val(peak1val);
  wf.setPeak2Val(peak2val);
  wf.setHalf1(half1);
  wf.setHalf2(half2);
  wf.setTermDist(termdist);
  wf.setMaxDevLoc(maxdevloc);
  wf.setRMSLeft(rms_left);
  wf.setRMSRight(rms_right);
  wf.setRMSOuter(rms_outer);
  wf.setRMSFull(rms_full);
  wf.setdt(dt);

  //Pulse parameters: only allow 1 pulse per waveform for CSP

  peak = hmax-base;
  peakBin = hmaxBin;
  peakTime = hmaxTime;

  threshold = base;
  startTime = analysis::startTime(htemp,peakBin,threshold,startBin);
  endTime = analysis::endTime(htemp,peakBin,threshold,endBin);

  for (int i = 0; i < 6; i++)
  {
    riseV[i] = rf[i] * peak + base;
    fallV[i] = rf[5-i] * peak + base;
  }
  
  analysis::riseTime(htemp,riseV,rise,max(startTime-binWidth,htemp->GetXaxis()->GetXmin()),peakTime,true);
  analysis::fallTime(htemp,fallV,fall,peakTime,min(endTime+binWidth,htemp->GetXaxis()->GetXmax()-binWidth*0.001),false);
  for (int i = 0; i<6; i++)
  {
    rise[i] = peakTime - rise[i];
    fall[i] = fall[i] - peakTime;
  }

//  rise[0] = peakTime-startTime;
  integral = (startTime-endTime) * base / binWidth +analysis::integral(htemp,startTime,endTime);
 // cout <<"Peak Bin: "<<peakBin<<endl;
 // cout <<"Peak position "<<peakTime<<endl;
 // cout <<"Peak "<<peak-base<<endl;
  FastPulse p(peakBin);
  p.setRise( &(rise[0]) );
  p.setFall( &(fall[0]) );
  p.setPeak(peak);
  p.setPeakTime(peakTime);
  p.setStartTime(startTime);
  p.setEndTime(endTime);
  p.setStartBin(startBin);
  p.setEndBin(endBin);
  p.setIntegral(integral);

  //Fast preamp parameters
  
  minBin = rise[3]<-0.5?1
	   :htemp->GetXaxis()->FindBin(peakTime-rise[3]);
  maxBin = fall[3]<-0.5?htemp->GetNbinsX()
	   :htemp->GetXaxis()->FindBin(fall[3]+peakTime);
  analysis::peaksAndValleys(htemp,pkBin,valBin,
		  minBin, maxBin, PEAK_FINDER_MIN_DIST);

//  cout <<"Number of peaks "<< pkBin.size() << endl;
 //// cout <<"Number of valleys "<< valBin.size()<<endl;

  analysis::mergePeaksByDistance(htemp,pkBin,valBin,MIN_PEAK_DIST);

  if(pkBin.size()>0){
    for (unsigned int i = 0; i < pkBin.size(); i++)
      {
	if (hget(htemp,pkBin[i])-base> MIN_FAST_PEAK*peak)
	  {
	    fastPeakBin = pkBin[i];
	    fastPeakTime = htemp->GetXaxis()->GetBinCenter(pkBin[i]);
	    fastPeak = hget(htemp,pkBin[i])-base;
	    slowPeakBin = pkBin[i];
	    slowPeakTime = htemp->GetXaxis()->GetBinCenter(pkBin[i]);
	    slowPeak = hget(htemp,pkBin[i])-base;
	    troughBin = pkBin[i];
	    troughTime = htemp->GetXaxis()->GetBinCenter(pkBin[i]);
	    trough = hget(htemp,pkBin[i])-base;
	    break;
	  }else{
	  // only erase if there's elements to erase
	  if(pkBin.size()>0)	  pkBin.erase(pkBin.begin());
	  if(valBin.size()>0)     valBin.erase(valBin.begin());
	}
      }
  }

  if( valBin.size()>0 ){
    int curMin = valBin[0];
    // if this is too big, then pkBin[i] overreaches its size
    for (unsigned int i = 1; i<pkBin.size();i++)
      {
	curMin = min(curMin,valBin[i-1]);
	if (pkBin[i]-pkBin[0]<MIN_PEAK_DIST) continue;
	if (hget(htemp,pkBin[i])-base <= slowPeak && slowPeakBin!=fastPeakBin) continue;
	slowPeak = hget(htemp,pkBin[i])-base;
	slowPeakBin = pkBin[i];
	slowPeakTime = htemp->GetXaxis()->GetBinCenter(pkBin[i]);
	trough = hget(htemp,curMin)-base;
	troughBin = curMin;
	troughTime = htemp->GetXaxis()->GetBinCenter(curMin);
      }
  }

  for (int i = 0; i < 6; i++)
  {
    riseV[i] = rf[i] * fastPeak + base;
    fallV[i] = rf[5-i] * slowPeak + base;
  }

  analysis::riseTime(htemp,riseV,rise,max(startTime-binWidth,htemp->GetXaxis()->GetXmin()),fastPeakTime,true);
  analysis::fallTime(htemp,fallV,fall,slowPeakTime,min(endTime+binWidth,htemp->GetXaxis()->GetXmax()-binWidth*0.001),false);


  for (int i = 0; i<6; i++)
  {
    
    rise[i] = fastPeakTime - rise[i];
    fall[i] = fall[i] - slowPeakTime;
  }

  p.setFastRise(&rise[0]);
  p.setSlowFall(&fall[0]);
  p.setFastPeak(fastPeak);
  p.setSlowPeak(slowPeak);
  p.setTroughHeight(trough);
  p.setFastPeakTime(fastPeakTime);
  p.setSlowPeakTime(slowPeakTime);
  p.setTroughTime(troughTime);
  p.setFastPeakBin(fastPeakBin);
  p.setSlowPeakBin(slowPeakBin);
  p.setTroughBin(troughBin);
  wf.add(p);

  delete htemp;
}
