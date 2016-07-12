// waveformtools includes
//  for analysis of generic waveforms
#include "../../../MaxCam/waveformtools/include/WaveformVector.hh"
#include "../../../MaxCam/waveformtools/include/SkimWaveform.hh"
#include "../../../MaxCam/waveformtools/include/DmtpcPulse.hh"
//  for analysis of voltage sensitive amplitifier waveforms
#include "../../../MaxCam/waveformtools/include/FastWfVector.hh"
#include "../../../MaxCam/waveformtools/include/FastWaveform.hh"
#include "../../../MaxCam/waveformtools/include/FastPulse.hh"
//  for analysis of charge sensitive pre-amp waveforms
#include "../../../MaxCam/waveformtools/include/CspWfVector.hh"
#include "../../../MaxCam/waveformtools/include/CspWaveform.hh"
#include "../../../MaxCam/waveformtools/include/CspPulse.hh"
//  for analysis of charge sensitive pmt waveforms
#include "../../../MaxCam/waveformtools/include/PMTWfVector.hh"
#include "../../../MaxCam/waveformtools/include/PMTWaveform.hh"
#include "../../../MaxCam/waveformtools/include/PMTPulse.hh"
//  general functions for waveform analysis
#include "../../../MaxCam/waveformtools/include/WaveformAnalysis.hh"
#include "../../../MaxCam/waveformtools/include/WaveformTools.hh"

#include "../../../MaxCam/MaxCamMC.hh"
#include "../../../MaxCam/MaxCamImageTools.hh"
#include "../../../MaxCam/MaxCamCluster.hh"
#include "../../../MaxCam/DmtpcEvent.hh"
#include "../../../MaxCam/MaxCamConfig.hh"
#include "../../../MaxCam/DmtpcDataset.hh"
#include "../../../MaxCam/DmtpcKeys.hh"
#include "../../../MaxCam/DmtpcSkimDataset.hh"
#include "../../../MaxCam/DmtpcSkimEvent.hh"
#include "../../../MaxCam/DmtpcSkimRunSummary.hh"
#include "../../../MaxCam/MaxCamTriggerGroup.hh"
#include "../../../MaxCam/ScopeWaveformData.hh"
#include "../../../MaxCam/ScopeDataInfo.hh"
#include "TBenchmark.h"
#include "TF1.h"
#include "TSystem.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TProfile.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TList.h"
#include "TCutG.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <list>
#include <stdlib.h>
#include <sstream>
#include <stdio.h>
#include "cleanSkimConfig.hh"
#include "cleanSkimProcessor.hh"

using namespace std;

/********* Global Variables *********/ 

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

void TurnKnobs();
void TurnKnobs(TString);
void analyze();
double findRMS(int, int, TH1D*, TH1D*, TH1D*);
double maxdevsearch(int, int, TH1D*, TH1D*, TH1D*);
int terminatorsearch(TH1D*, TH1D*, TH1D*, int, int, double, vector<int>&);
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

void TurnKnobs()
{
  ifstream ifstr(knobfile);
  if (!ifstr.is_open()){
    cout << knobfile << " not found." << endl;
  }
  string line;
  TString name;
  TString s;
  Int_t n;
  double_t d;
  while (!ifstr.eof()){
    getline(ifstr,line);
    if (line[0]=='#' || line=="") continue;
    istringstream linestr(line);
    linestr >> name;
    if (name == "RootFile"){
      linestr >> s;
      rootfile = s;
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
void TurnKnobs(TString knobfile_)
{knobfile = knobfile_;
  TurnKnobs();}
//Automatically sets the settings below                                                                                   

void analyze()
//runs analysis                                                                                                           
{
  term2pow = .5;
  yscalestd = xscalestd * peakstd / braggmaxval;

  TFile fbragg(braggfile);
  bragghist = (TH1D*)fbragg.Get("bragghist");
  bragghist2 = (TH1D*)fbragg.Get("bragghist2");
  bragghist3 = new TH1D("","",16384,0,16384);
  bragghist4 = new TH1D("","",16384,0,16384);

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

  histsum3 = new TH1D("","",16384,0,16384);
  for (int i = 0; i < 16384; i++)
    histsum3->SetBinContent(i,histA3->GetBinContent(i)+histB3->GetBinContent(i));

  test_rec_wf = new TH1D("test_rec_wf","",16384,0,16384);
  for (int j = 0; j < 16384; j++) test_rec_wf->SetBinContent(j,rough_wf->GetBinContent(j));
  for (int j = 0; j < 1; j++) smooth(test_rec_wf,rec_SD,gaus_acc);
  for (int j = 0; j < 1; j++) correct(test_rec_wf,rec_cwidth * rec_SD,gaus_acc,rec_cor);

  test_minfinder_curve = new TH1D("","",16384,0,16384);
  minprepare(test_rec_wf,test_minfinder_curve);
  wfd_delta2 = int(wfd_delta * dt * 1e9);
  cout << endl << "OISDJGOISJDOGIJSDOGIJSDG" << wfd_delta2 << endl;
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

  rec_wf = new TH1D("rec_wf","",16384,0,16384);
  for (int j = 0; j < 16384; j++) rec_wf->SetBinContent(j,rough_wf->GetBinContent(j));
  for (int j = 0; j < 1; j++) smooth(rec_wf,rec_SD,gaus_acc);
  for (int j = 0; j < 1; j++) correct(rec_wf,rec_cwidth * rec_SD,gaus_acc,rec_cor);

  jmin = 0;
  jbragg = 0;
  jdiff = 0;
  minfinder_curve = new TH1D("","",16384,0,16384);
  minprepare(rec_wf, minfinder_curve);
  jmin = minfinder(rec_wf,minfinder_curve);

  mfc2 = new TH1D("","",16384,0,16384);
  for (int i = 0; i < 16384; i++)
    mfc2->SetBinContent(i,minfinder_curve->GetBinContent(i));
  if (decon_on == 1)
    {
      deconvolve(mfc2,decon_SD,decon_ord);
      smooth(mfc2,finalsmooth,gaus_acc);
      correct(mfc2,finalsmooth*1.05,gaus_acc,.96);
    }
  mfcdiff = new TH1D("","",16384,0,16384);
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

  jterm  =  terminatorsearch(bragghist3,bragghist4,mfc2,peak1,peak2,termthresh,terminator_results);
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

  leftint  = histspline(bragghist,half1,halfheight,jcent,1/(cos1true * xscalestd),1);
  rightint = histspline(bragghist,half2,halfheight,jcent,-1/(cos2true * xscalestd),1);
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

int terminatorsearch(TH1D* bragghist3, TH1D* bragghist4, TH1D* mfc2, int xmin, int xmax, double threshold, vector<int>& vec)
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
  bestfit = -1.;
  int bestfitloc = 0;
  for (int x = xmin; x < xmax; x++)
    {
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
	  thisval = (valmfc > (a3 * val3 + a4 * val4)) ? 1 : -1;
	  if (thisval * lastval == -1)
	    {
	      vec.push_back(x);
	      thisfit = fabs(val3 - val4);
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
	    }
	}
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
      bragghist3->SetBinContent(x-1 + xoff,histspline(bragghist2,0,0,x,1/xscaletrue,yscale / xscale));
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
  }
  void deconvolve(TH1D* hist, double SD, int order)
    //Un-does Gaussian smearing. May contain traces of sorcery.
  {
    TH1D* hist2 = new TH1D("","",16384,0,16384);
    for (int i = 2000; i < 10000; i++)
      {
        if (order == 2)
          hist2->SetBinContent(i,hist->GetBinContent(i) - 0.25 * pow(SD*dt,2) * deriv2(hist,i,d2_step,dt));
        else if (order == 4)
	  hist2->SetBinContent(i,hist->GetBinContent(i) - 0.25 * pow(SD*dt,2) * deriv2(hist,i,d2_step,dt) + .125 * pow(SD*dt,4) * deriv4(hist,i,d4_step,dt));
        else
          cout << "invalid value for variable 'order' (must be 2 or 4.)" <<endl;
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

  int minfinder(TH1D* wfder)
  {
    double thresh = .2 * wfder->GetBinContent(wfder->GetMaximumBin());
    TH1D* dummy = new TH1D("","",16384,0,16384);
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

    TH1D* wfder = new TH1D("","",16384,0,16384);

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
    TH1D* hist2 = new TH1D("","",16384,0,16384);
    vector<double> vec;
    makekernel(vec,SD,acc);

    for(int j=2000; j < 10200; j++)
      hist2->SetBinContent(j,average(j,vec,hist));
    for(int j=1; j < 16384; j++)
      hist->SetBinContent(j,hist2->GetBinContent(j));
  }

  void correct(TH1D* hist, double SD, double acc, double factor)
    //Un-does some of the undesirable effects of smooth(). 
  {
    TH1D* hist2 = new TH1D("","",16384,0,16384);
    vector<double> vec;
    makekernel(vec,SD,acc);
    for(int j=2000; j < 10200; j++)
      hist2->SetBinContent(j,(1+factor)*hist->GetBinContent(j)-factor*average(j,vec,hist));
    for(int j=1; j < 16384; j++)
      hist->SetBinContent(j,hist2->GetBinContent(j));
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

/* End Waveform Reconstruction vars/fns */

//Printout levels: | together 

static const u_int32_t dbg_error     =  1;     //Fatal unexpected condition printouts (these go to cerr)
static const u_int32_t dbg_warning   =  2;     //Non-fatal unexpected condition printouts
static const u_int32_t dbg_progress  =  4;     //Progress printouts
static const u_int32_t dbg_bias      =  8;     //Bias algorithm printouts
static const u_int32_t dbg_cluster   =  16;    //Cluster finding printouts
static const u_int32_t dbg_burnin    =  32;    //Burnin printouts
static const u_int32_t dbg_config    =  64;    //Config printout

static u_int32_t dbg_level = dbg_progress | dbg_error | dbg_warning; 




/******** Prototypes **********/
//static void samepos_copy_next(int * burnin, int * nburnin, list<int *> * same_positions, int ncamera, int );a

static int parse_args( int nargs, char ** args, TString * files, TString * keys, TString * out, char ** ); 
//static void handle_legacy_events(int,int,DmtpcDataset &, TH2F *, double *, double *, TH2F **, int *, double *, int, CleanSkimConfig *);
static bool burnin_test( double xdelta, double ydelta, CleanSkimConfig * conf);
/******** end prototypes ******/

/******** Waveform analysis prototypes **********/
int fillWaveformVectorsInTObjArray( DmtpcDataset & d, TObjArray * wfvlist );
int fillCspWfVector(DmtpcDataset & d, CspWfVector *cspwv, int ich, int nch);
int fillPMTWfVector(DmtpcDataset & d, PMTWfVector *cspwv, int ich, int nch);
int fillFastWfVector(DmtpcDataset & d, FastWfVector *fastwv, int ich, int nch);
/******** end waveform analysis prototypes ******/

/******** cleanSkim *********/
int cleanSkim(DmtpcDataset & d, TString & key, TTree * rawtree, int runnum, 
              DmtpcSkimDataset & sd, bool delete_intermediate, CleanSkimConfig * conf, unsigned algo_index)
{

  std::cout <<"hasOverscan(): " << conf->hasOverscan() << std::endl; 

  TBenchmark* myBench=new TBenchmark();

  TString tmpoutfilename = TString(sd.getFileName()).ReplaceAll(".root","_tmp.root");

  /* Define Temporary out files and trees */
  TFile* tmpoutfile = new TFile(tmpoutfilename,"RECREATE");
  TTree* tmpskimtree;
  TH2F* tempimg=0;
  TH2F* os_tempimg=0;


  /* Get some eventwide information */
  d.getEvent(0);
  const int ncamera = d.event()->ccdData()->GetEntries();
  const int nbinsx =  d.event()->ccdData(0)->GetNbinsX();
  const int nbinsy =  d.event()->ccdData(0)->GetNbinsY();
  const int eventVersion = d.event()->IsA()->GetClassVersion();
  int ncam = ncamera;

  cout << "ncamera " << ncamera << endl;

  /* Support processing on datasets with no actual images */ 
  bool fake_image = false; 
  if (ncamera ==0 || d.event()->ccdData(0)==NULL || d.event()->ccdData(0)->Integral() ==0)
  {
    fake_image = true; 
    tempimg = new TH2F ("tempimg","tempimg",256,0,1024,256,0,1024); 
  }

  cout << "Getting camera id" << endl;
  	
 //Get camera id's 
 string camera_ids[ncamera];  
 for (int u = 0; u < ncamera; u++)
 {
    camera_ids[u] = string( ((MaxCamConfig*)d.event()->ccdConfig(u))->serialNumber.Data()); 

    cout << "camera_ids " << camera_ids << endl;

    if (camera_ids[u] == "")
    {
      camera_ids[u] = std::string(conf->getFallBackCameraID(u)); 
      if (camera_ids[u]!="")
      {
        std::cout << "Fell back to camera id: " << camera_ids[u] << std::endl; 
      }
    }
 }

 cout << "Loading complete" << endl;

  /* Create the temporary skim tree without burnin data.  */
  tmpskimtree = new TTree(key,"Skimmed Events");
  double theta[ncamera][15], phi[ncamera][15], E[ncamera][15];
  double EGainMap[ncamera][15],diffusedrange[ncamera][15];
  double range[ncamera][15], x[ncamera][15], y[ncamera][15];
  int xbegin[ncamera][15], ybegin[ncamera][15], xend[ncamera][15], yend[ncamera][15];
  double skewness[ncamera][15], integral[ncamera];
  bool edge[ncamera][15], spark[ncamera];
  int ntracks[ncamera], lastspark[ncamera],date, time, eventnum;
  double timenow[ncamera];
  int pixels_killed[ncamera], npixel_red[ncamera][15];
  double image_mean[ncamera],image_rms[ncamera];
  double os_mean[ncamera],os_rms[ncamera];
  MaxCamClusterImage* clusti=0;
  double cluster_rms[ncamera][15], cluster_mean[ncamera][15], energy_density[ncamera][15];
  int neighbors [ncamera][15], npixel[ncamera][15]; 
  double maxpixel[ncamera][15], cygnus_angle[ncamera][15]; 
  double moments[ncamera][4][15]; 
  double transverse_moments[ncamera][4][15]; 
  double ra[ncamera][15], dec[ncamera][15]; 
  double glat[ncamera][15],  glon[ncamera][15]; 
  DmtpcEvent* skimevent = new DmtpcEvent();
  TObjArray* clust = new TObjArray(ncamera);
  clust->SetOwner(kTRUE);  
  TObjArray* trigger_groups = new TObjArray;
  vector< string > * cameraSerialNumber = new vector< string >;
  double majoraxis[ncamera][15],minoraxis[ncamera][15];
  // make a fixed length array of waveform vectors, with only as many entries as there are scope channels
  TObjArray* waveform_vectors = new TObjArray( conf->getNChannelsPerTrigger() );
  for(Int_t ich=0; ich<conf->getNChannelsPerTrigger(); ++ich){
    waveform::tools::addWaveformVectorToTObjArray( waveform_vectors,
						   ich,
						   conf->getChannelId(ich),
						   conf->getChannelType(ich) );
  }

  waveform_vectors->SetOwner(kTRUE);
  //  std::cout << "conf->getNChannelsPerTrigger()=" <<  conf->getNChannelsPerTrigger() << std::endl;
  //  std::cout << "waveform_vectors->GetEntries()=" << waveform_vectors->GetEntries() << std::endl;
  //  std::cout << "waveform_vectors->GetSize()=" << waveform_vectors->GetSize() << std::endl;


  /*^^^

  cout << "CHECKPOINT 1" << endl;

  knobfile = "runParameters/LittleDCTPC_far/WfKnobs.temp";
  TurnKnobs();
  analyze();

  cout << "CHECKPOINT 1.1" << endl;
  tmpskimtree2 = new TTree("wf"+key,"waveform analyses");
  cout << "CHECKPOINT 2" << endl;

  //tmpskimtree->Branch("nwaveform",&ncam);
  cout << "CHECKPOINT 3" << endl;
  tmpskimtree->Branch("leftint",&leftint);
  cout << "CHECKPOINT 3.1" << endl;
  tmpskimtree->Branch("leftint2",&leftint2);
  tmpskimtree->Branch("rightint",&rightint);
  tmpskimtree->Branch("rightint2",&rightint2);
  tmpskimtree->Branch("jmin",&jmin);
  tmpskimtree->Branch("jbragg",&jbragg);
  tmpskimtree->Branch("jterm",&jterm);
  tmpskimtree->Branch("jterm2",&jterm2);
  tmpskimtree->Branch("jterm3",&jterm3);
  tmpskimtree->Branch("origin",&origin);
  cout << "CHECKPOINT 4" << endl;

  if (terminator_results.size() > 0)
  {
    for (int i = 0; i < terminator_results.size(); i++)
      tmpskimtree->Branch("termmatch" + i,&terminator_results.at(i));
  }
 
 tmpskimtree->Branch("rec_SD",&rec_SD);
 tmpskimtree->Branch("wfd_delta",&wfd_delta);
 tmpskimtree->Branch("peak1",&peak1);
 tmpskimtree->Branch("peak2",&peak2);
 tmpskimtree->Branch("peak1val",&peak1val);
 tmpskimtree->Branch("peak2val",&peak2val);
 tmpskimtree->Branch("half1",&half1);
 tmpskimtree->Branch("half2",&half2);
 tmpskimtree->Branch("termdist",&termdist);
 tmpskimtree->Branch("maxdevloc",&maxdevloc);
 tmpskimtree->Branch("rms_left",&rms_left);
 tmpskimtree->Branch("rms_right",&rms_right);
 tmpskimtree->Branch("rms_outer",&rms_outer);
 tmpskimtree->Branch("rms_full",&rms_full);
 tmpskimtree->Branch("dt",&dt);

 cout << endl << endl << endl << "HEY LOOK OVER HERE @@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl << endl << endl;

 ^^^*/


 tmpskimtree->Branch("ncamera",&ncam,"ncamera/I");
  tmpskimtree->Branch("spark",&spark,"spark[ncamera]/O");
  tmpskimtree->Branch("lastspark",&lastspark,"lastspark[ncamera]/I");
  tmpskimtree->Branch("integral",&integral,"integral[ncamera]/D");
  tmpskimtree->Branch("eventnum",&eventnum,"eventnum/I");
  tmpskimtree->Branch("theta",theta,"theta[ncamera][15]/D");
  tmpskimtree->Branch("phi",phi,"phi[ncamera][15]/D");
  tmpskimtree->Branch("E",E,"E[ncamera][15]/D");
  tmpskimtree->Branch("EGainMap",EGainMap,"EGainMap[ncamera][15]/D");
  tmpskimtree->Branch("energy_density", energy_density, "energy_density[ncamera][15]/D");
  tmpskimtree->Branch("range",range,"range[ncamera][15]/D");
  tmpskimtree->Branch("majoraxis",majoraxis,"majoraxis[ncamera][15]/D");
  tmpskimtree->Branch("minoraxis",minoraxis,"minoraxis[ncamera][15]/D");
  tmpskimtree->Branch("diffusedrange",diffusedrange,"diffusedrange[ncamera][15]/D");
  tmpskimtree->Branch("x",x,"x[ncamera][15]/D");
  tmpskimtree->Branch("y",y,"y[ncamera][15]/D");  
  tmpskimtree->Branch("xbegin",xbegin,"xbegin[ncamera][15]/I");
  tmpskimtree->Branch("xend",xend,"xend[ncamera][15]/I");
  tmpskimtree->Branch("ybegin",ybegin,"ybegin[ncamera][15]/I");
  tmpskimtree->Branch("yend",yend,"yend[ncamera][15]/I");  
  tmpskimtree->Branch("ntracks",ntracks,"ntracks[ncamera]/I");
  tmpskimtree->Branch("edge",edge,"edge[ncamera][15]/O");
  tmpskimtree->Branch("skewness",skewness,"skewness[ncamera][15]/D");
  tmpskimtree->Branch("ra",ra,"ra[ncamera][15]/D");
  tmpskimtree->Branch("dec",dec,"dec[ncamera][15]/D");
  tmpskimtree->Branch("glat",glat,"glat[ncamera][15]/D");
  tmpskimtree->Branch("glon",glon,"glon[ncamera][15]/D");
  tmpskimtree->Branch("moments",moments,"moments[ncamera][4][15]/D");
  tmpskimtree->Branch("transverse_moments",transverse_moments,"transverse_moments[ncamera][4][15]/D");
  tmpskimtree->Branch("date",&date,"date/I");
  tmpskimtree->Branch("timenow",&timenow,"timenow[ncamera]/D");
  tmpskimtree->Branch("clusters","TObjArray",&clust,128000,0);
  tmpskimtree->Branch("cluster_rms",cluster_rms,"cluster_rms[ncamera][15]/D");
  tmpskimtree->Branch("cluster_mean",cluster_mean,"cluster_mean[ncamera][15]/D");
  tmpskimtree->Branch("maxpixel",maxpixel,"maxpixel[ncamera][15]/D");
  tmpskimtree->Branch("neighbors",neighbors,"neighbors[ncamera][15]/I");
  tmpskimtree->Branch("cygnus_angle", cygnus_angle, "cygnus_angle[ncamera][15]/D"); 
  tmpskimtree->Branch("npixel", npixel, "npixel[ncamera][15]/I"); 
  tmpskimtree->Branch("npixel_red",npixel_red,"npixel_red[ncamera][15]/I");
  tmpskimtree->Branch("pixels_killed",pixels_killed,"pixels_killed[ncamera]/I");
  tmpskimtree->Branch("image_mean",image_mean,"image_mean[ncamera]/D");
  tmpskimtree->Branch("image_rms",image_rms,"image_rms[ncamera]/D");
  //tmpskimtree->Branch("trigger_groups","TObjArray",&trigger_groups,128000,0);
  tmpskimtree->Branch("cameraSerialNumber",&cameraSerialNumber);
  if (algo_index == 0)
    tmpskimtree->Branch("waveform_vectors","TObjArray",&waveform_vectors,128000,0);

  gROOT->cd();

  /* Open up file for bias frame saving */
  TFile * biasoutfile;
  if (algo_index == 0)
  {
    TString biasoutfilename = TString(sd.getFileName()).ReplaceAll(key+".root","bias.root");
    biasoutfilename.ReplaceAll(key+".root","bias.root");
    biasoutfile = new TFile(biasoutfilename,"RECREATE");
  }

  TH2F* biasframe;

  if(conf->Bias())
    {
  biasframe = d.event()->ccdData(0);
    }
  else
    {
  biasframe=new TH2F("biasframe","biasframe",256,-512,512,256,-512,512);
  std::cout<<"NO BIAS"<<std::endl;
    }

  TTree* biastree = new TTree("bias","Bias Information");
  biastree->Branch("biasframe","TH2F",&biasframe,128000,0);
  gROOT->cd();

  double sigmathr[ncamera];
  int nframes[ncamera];
  for(int u=0; u<ncamera; u++){nframes[u]=0;}
  TH2F* secondarybias[ncamera];
  TH2F* os_secondarybias[ncamera];
  double bias_mean[ncamera], bias_rms[ncamera],os_bias_mean[ncamera],os_bias_rms[ncamera];
  double cleaned_os_bias_mean[ncamera];
  double last_mean[ncamera]; // for the spark cut

  // handle the bias frame
  //      define the bias threshold for runs<500 and MC
  for(int u=0; u<ncamera; u++)
    {
      if (fake_image)
      {
         break; 
      }
      if (dbg_level & dbg_bias) std::cout << "Camera " << u << " image cleaning:" << std::endl;

      if(conf->Bias())
      secondarybias[u] = (TH2F*)d.getBiasFrame(u+1)->Clone();
      else
      secondarybias[u] = (TH2F*)biasframe->Clone();	
      // the cluster finding threshold is very sensitive to this number!
      // using 1.5 sigma reproduces cluster finding efficiency from AF v1
      sigmathr[u]=1.5*MaxCamImageTools::getRMS(secondarybias[u]);
      
      
      if (dbg_level & dbg_bias)
      {  
        std::cout << "     camera bias MEAN: " << 
             MaxCamImageTools::getMean(secondarybias[u]) << std::endl;
        std::cout << "     camera bias RMS: " << 
             MaxCamImageTools::getRMS(secondarybias[u]) << std::endl;
      }
    }

  
  const int nev = rawtree->GetEntries();
  
  /* Bias frame stuff */
  for(int u=0; u<ncamera; u++)
    {

      if (fake_image) 
      {
        secondarybias[u] = new TH2F("bias","bias",256,0,1024,256,0,1024); 
      }
      else
      {
        if(nframes[u] !=0 && runnum < 500) {
          if (dbg_level & dbg_bias) std::cout << "run number " << runnum << " using averaged bias " << std::endl;
          secondarybias[u]->Scale(1/double(nframes[u]));
        } else {
	  if(conf->Bias())
          secondarybias[u] = (TH2F*)d.getBiasFrame(u+1)->Clone();
	  else
	  secondarybias[u] = (TH2F*)biasframe->Clone();  
  
          if (conf->hasOverscan())
          {
            os_secondarybias[u]= (TH2F*)d.getBiasFrameOverscan(u+1)->Clone();
          }
          if (dbg_level & dbg_bias) std::cout << "run number " << runnum << " using camera bias " << std::endl;
        }
        // kill 5 sigma outliers in the bias frame, iterate 3 times
        for (int it=0; it<3; it++) {
           double mean, rms,osmean,osrms;
          mean=MaxCamImageTools::getMean(secondarybias[u]);
          rms=MaxCamImageTools::getRMS(secondarybias[u]);
          if (conf->hasOverscan())
          {
            osmean=MaxCamImageTools::getMean(os_secondarybias[u]);
            osrms=MaxCamImageTools::getRMS(os_secondarybias[u]);
          }
          if (conf->isMC()) { // MC
            MaxCamImageTools::killPixels(secondarybias[u],mean+10.*rms);
          } else {
            MaxCamImageTools::killPixels(secondarybias[u],mean*conf->getOutlierFactor());
            if (conf->hasOverscan())
            {
              MaxCamImageTools::killPixels(os_secondarybias[u],osmean*conf->getOutlierFactor());
            }
          }
          if (dbg_level & dbg_bias) 
          {
            std::cout << "     iteration " << it << " averaged bias image mean (no outliers): " << mean << std::endl;
            std::cout << "     iteration " << it << " averaged bias image rms (no outliers): " << rms << std::endl;
          }
          bias_mean[u]=mean;
          bias_rms[u]=rms;
          if (conf->hasOverscan())
          {
	    os_bias_mean[u]=osmean;
            os_bias_rms[u]=osrms;
          }
        }
	
	if (conf->hasOverscan())
          {
	    // we use this later to correct how much images have drifted in overall mean 
	    // since the bias frame was taken.  It's important that we use the cleaned image
	    // for this mean because we might otherwise pick up a hot pixel
	    cleaned_os_bias_mean[u]=MaxCamImageTools::getMean(os_secondarybias[u]);
          }
      }
      biasframe = secondarybias[u];
      if (algo_index==0)
      {
        biasoutfile->cd();
        biastree->Fill();
      }
      gROOT->cd();   
    }
    
  if (dbg_level & dbg_progress) std::cout << "All preclean activities done" << std::endl;
  int nempty=0;

  /* Running list of positions so we don't have to read back through skim file
     all the time. It is a List of a Vector[2] of a double[2].*/
  list<vector<double *>**> positions; 
  

  /* List of spark ref excluded pixels */
  vector<pair<int,int> > * sparkref_running = new vector<pair<int,int> >[ncamera]; 
  list<vector<pair<int,int> >*>  sparkref; 

  list<vector<vector<vector<BurninEncoded_t> > >*> burnin_temp; //Temporary list that is updated a lot
  list<vector<vector<vector<BurninEncoded_t> > >*> burnin; //permanent list 
  
  int position_offset  = 0; 

  //define cut for sparks with overscan
  TCutG* os_cut;
  if (conf->hasOverscan())
  {
    os_cut = new TCutG("os_cut",5);
    os_cut->SetPoint(0,0,1028);
    os_cut->SetPoint(1,1024,1028);
    os_cut->SetPoint(2,1024,1032);
    os_cut->SetPoint(3,0,1032);
    os_cut->SetPoint(4,0,1028);
  }

  //Load gain maps if we are using them; 

  DmtpcGainMap * gainmaps[ncamera]; 

  if (conf->useGainMap())
  {
    TFile gfile(conf->getGainMapFile()); 
    gfile.ls(); 
    for (int u = 0; u < ncamera; u++)
    {
      std::string ser =  camera_ids[u]; 
      gainmaps[u]=  ser == "" ? 0 : (DmtpcGainMap*) gfile.Get(ser.c_str())->Clone();
      if (conf->normalizeGainMap()) 
      {
        gainmaps[u]->getGainMap()->Scale(gainmaps[u]->getGainMap()->GetNbinsX()*gainmaps[u]->getGainMap()->GetNbinsY()/gainmaps[u]->getGainMap()->Integral());
      }
      sd.addGainMap((DmtpcGainMap*)gainmaps[u]->Clone());
    }
    sd.writeGainMaps(); 
    gfile.Close(); 
  }
  else
  {
    for (int u = 0; u < ncamera; u++)  
    {
      gainmaps[u] = 0; 
    }
  }

  // loop over events

  std::cout << "nev " << nev << std::endl;

  for(int i = 0; i<nev; i++)
  {
    //cout<<i<<endl;
    myBench->Start("Event");

    if (dbg_level & dbg_progress) std::cout << "Processing event " << i << std::endl;
    d.getEvent(i);
    int sumtracks=0;
    int sparks=0;
    clust->Clear();
    trigger_groups->Clear(); 
    
    vector<vector<vector<BurninEncoded_t> > > * this_event_burnin= new vector<vector<vector<BurninEncoded_t> > >;
    myBench->Start("Saving scope data");
        
    

    if (algo_index == 0 )
      fillWaveformVectorsInTObjArray( d , waveform_vectors );
    

// if(d.event()->scopeDataInfo(0)!=NULL)
// {    
//      FastWfVector* anode = (FastWfVector*) waveform_vectors->At(0);
//      int size=anode->size();
// 
// 
// 
// for(int j=0;j<size;j++)
//  	{
//  		if(j>=15)
//  		continue;
//  	
//     ScopeWaveformData *hwf2=d.event()->rawScopeData(j);
//     triggertime[j]=hwf2->getTimeStamp();
// 
//  	}
// 
// }

    //FastWfVector* anode = (FastWfVector*) waveform_vectors->At(0);
    //cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!anode size "<<anode->size()<<endl;
    
    myBench->Show("Saving scope data");
    
    myBench->Start("Loop over cameras");
    
    for(int u=0; u<ncamera; u++)
      { 
	timenow[u]=d.event()->timeStamp()->Convert();
	// cout<<"TIMENOW !!!!!!!!!!!!!!!!!!!!!!!!"<<timenow[0]<<endl;
	//save the camera serial numbers
        cameraSerialNumber->push_back(camera_ids[u]); 
	
        vector< vector<BurninEncoded_t> >  this_camera_burnin ;
        MaxCamClusterImage::CAMERA_ORIENTATION orientation = conf->getCameraOrientation(camera_ids[u].c_str());         
	
        if (dbg_level & dbg_progress) std::cout << "\t" << u << ":" << std::endl;
        // keep track of the last frame for the spark cut
        if (i==0) {
          last_mean[u] = bias_mean[u];
          lastspark[u]=0;
        }
	
        if(!fake_image)
	  {
	    tempimg = (TH2F*)d.event()->ccdData(u)->Clone("tempimg");
	    os_tempimg = conf->hasOverscan() ? (TH2F*)d.event()->overscan(u)->Clone("tempimg") : 0;
	  }
	
        double mean, rms, osmean, osrms;
        mean= fake_image ? 0 : MaxCamImageTools::getMean(tempimg);
        rms= fake_image ? 0 : MaxCamImageTools::getRMS(tempimg);
        osmean= fake_image  || !conf->hasOverscan() ? 0 : MaxCamImageTools::getMean(os_tempimg);
        osrms= fake_image || !conf->hasOverscan() ? 0 : MaxCamImageTools::getRMS(os_tempimg);
	
        double mean_nooutliers = 0; 
        double rms_nooutliers = 0;
	
        double os_mean_nooutliers = 0; 
        double os_rms_nooutliers = 0;
	myBench->Start("meanRMSNoOutliers");
        if (!fake_image)
	  {
	    MaxCamImageTools::meanRMSNoOutliers(tempimg, mean_nooutliers, rms_nooutliers);          
	    if (conf->hasOverscan())
	      {
		MaxCamImageTools::meanRMSNoOutliers(os_tempimg, os_mean_nooutliers, os_rms_nooutliers);          
	      }
	  }
	myBench->Show("meanRMSNoOutliers");

        if (dbg_level & dbg_cluster)
	  {
	    std::cout << "     event: " << i << std::endl;
	    std::cout << "     image raw mean: " << mean << std::endl;
	    std::cout << "     image raw rms: " << rms << std::endl;
	    std::cout << "     image n.o. mean: " << mean_nooutliers << std::endl;
	    std::cout << "     image n.o. rms: " << rms_nooutliers << std::endl;
	    std::cout << "     bias mean: " << bias_mean[u] << std::endl;
	    std::cout << "     bias rms: " << bias_rms[u] << std::endl;
	    if (runnum > 0) 
	      std::cout << "     spark test: " << (mean/last_mean[u]) << std::endl;
	    
	  }
        image_mean[u] = mean_nooutliers;
	image_rms[u] = rms_nooutliers;
        os_mean[u] = os_mean_nooutliers;
        os_rms[u] = os_rms_nooutliers;
        // check for sparks
	
        double os_spark_mean = conf->hasOverscan() ? MaxCamImageTools::getMean(os_tempimg,os_cut): 0 ;
	
	if( fake_image || 
	    conf->isMC() ||
	    ( mean/last_mean[u]<conf->getSparkCut(const_cast<char*>(camera_ids[u].c_str())) && !conf->hasOverscan() ) || 
	    ( (mean/last_mean[u])<conf->getSparkCut(const_cast<char*>(camera_ids[u].c_str())) && 
	      ( conf->hasOverscan() && mean-os_spark_mean< conf->getOverscanSparkCut(const_cast<char*>(camera_ids[u].c_str())) ) ) )
	  {
	    spark[u]=0;
	    last_mean[u]=mean;
	  } else {
          spark[u]=1;
          //std::cout << "overscan spark cut for " << ser << ": " << conf->getOverscanSparkCut(const_cast<char*>(ser.Data())) << std::endl; 
	  
          //Populate sparkref 
          TH2F* sparkrefimg = (TH2F*) tempimg->Clone("sparkrefimg"); 
          sparkrefimg->Add(secondarybias[u],-1); 
	  
          int nsat = 0;  
	  myBench->Start("Loading sparkrefimg");
          for (int m = 1; m <= sparkrefimg->GetNbinsX(); m++)
	    {
	      for (int n = 1; n <= sparkrefimg->GetNbinsY(); n++)
		{
		  if (sparkrefimg->GetBinContent(m,n) > conf->getSatThresh())
		    {
		      nsat++; 
		    }
		}
	    }
	  myBench->Show("Loading sparkrefimg");
	  
	  myBench->Start("Thresholding");
          if (nsat > conf->getNSatThresh())
	    {
	      for (int m = 1; m <= sparkrefimg->GetNbinsX(); m++)
		{
		  for (int n = 1; n <= sparkrefimg->GetNbinsY(); n++)
		    {
		      if (sparkrefimg->GetBinContent(m,n) > conf->getSatThresh())
			{
			  std::cout << "Adding sparkref pixels " << m << "," << n << std::endl; 
			  sparkref_running[u].push_back(pair<int,int>(m,n)); 
			}
		    }
		}
	    }
	  myBench->Show("Thresholding");
	  
          sparkrefimg->Delete(); 
        }
        
        if (dbg_level & dbg_cluster) 
	  {
	    if(spark[u] == 1) std::cout << "     SPARK!" << std::endl;
	  }
	
        // clean up image before cluster finding, same way as bias frame
        if (conf->isMC()) { // MC
          if (!conf->isNoKill())
            pixels_killed[u] = MaxCamImageTools::killLonePixels2(tempimg, mean + 10*rms);
        } else if (fake_image)
	  {
            pixels_killed[u] = 0; 
	  } else {
          pixels_killed[u] = MaxCamImageTools::killLonePixels2(tempimg,conf->getOutlierFactor()*mean);
          if (conf->hasOverscan())
	    {
	      MaxCamImageTools::killLonePixels2(os_tempimg,conf->getOutlierFactor()*osmean);
	    }
        }
	
	// we use this later to correct how much images have drifted in overall mean 
	// since the bias frame was taken.  It's important that we use the cleaned image
	// for this mean because we might otherwise pick up a hot pixel
	double cleaned_osmean=conf->hasOverscan() ? MaxCamImageTools::getMean(os_tempimg) : 0;;
	
        if (!fake_image)
	  {
	    tempimg->Add(secondarybias[u],-1);
	    if (conf->hasOverscan())
	      {
		os_tempimg->Add(os_secondarybias[u],-1);
	      }
	  }
        //      subtract off remaining pedestal
        //
	//        double perpx = fake_image ? 0 : tempimg->Integral()/(nbinsx*nbinsy);
	
	
        std::cout << i << "," << u << ": " << bias_mean[u] << "," << mean << ";" 
		  << os_bias_mean[u] << "," << osmean << std::endl;
	
	// this is the correction meant to account for a possible drift from the mean since the bias frame was taken
	// it's extremely important to use the hot pixel *cleaned* overscan means for this drift correction, 
	double perpx = fake_image ? 0 : conf->hasOverscan() ? cleaned_osmean-cleaned_os_bias_mean[u] : tempimg->Integral()/(nbinsx*nbinsy) ;
	
	std::cout << perpx << std::endl;
	if (dbg_level & dbg_cluster) std::cout << "     image pedestal correction: " << perpx << std::endl;
	//look for partials here

	myBench->Start("precluster");

	if(spark[u]==0 && !fake_image)
	  {
	    MaxCamImageTools::subtractPedestal(tempimg,perpx);
	    
	    //look for partials here
	    TH2C* edges = MaxCamImageTools::edgeDetect(tempimg,
						       conf->getPartialBlur(),
						       conf->getPartialLowThresh(),
						       conf->getPartialHighThresh());
	    
	    TH1D* edges_p = edges->ProjectionX("edges_p",2,edges->GetNbinsX()-1);
	    
	    const int validbins=tempimg->GetNbinsX();
	    double edgesval[validbins];
	    int edgesind[validbins];
	    for(int j=2; j<edges_p->GetNbinsX(); j++)
	      {
		edgesval[j-1]=edges_p->GetBinContent(j);
	      }
	    edgesval[0]=0; edgesval[255]=0;
	    TMath::Sort(validbins,edgesval,edgesind);
	    
	    if(edgesval[edgesind[0]]>conf->getPartialNPrimaryThresh() )
	      {
		int p=1;
		bool secedges=false;
		while(edgesval[edgesind[p]] > conf->getPartialNSecondaryThresh())
		  {
		    if(abs(edgesind[p]-edgesind[0])>=conf->getPartialDistanceLow()&& 
		       abs(edgesind[p]-edgesind[0])<=conf->getPartialDistanceHigh())
		      {
			secedges=true;
		      }
		    p++;
		  }
		spark[u]=!secedges;
	      }
	    
	    delete edges; delete edges_p;
	  }
	
	integral[u]=rms; // this variable is written to the output tree
	                 // useful to use rms here since the spark cut
	// now cuts on rms (instead of integral) 
	
	
        // save spark images            
        if(spark[u]==1) 
	  {
	    sparks++;
	    ntracks[u] = 0; 
	    lastspark[u]=0;
	    TH2F * sparkimg = (TH2F*) tempimg->Clone("sparkimg");
	    if(eventVersion==1)
	      {
		clusti= new MaxCamClusterImage(sparkimg,d.event()->timeStamp());
	      }
	    else
	      {
		clusti= new MaxCamClusterImage(sparkimg,d.event()->UTCtimeStamp());
	      }
	    
	    /* Find oversaturated pixels */           
	    
	    
	  }

	myBench->Show("precluster");

	myBench->Start("Cluster finding");
        if (fake_image)
	  {
	    clusti = new MaxCamClusterImage((TH2F*)tempimg->Clone("baseimg"),d.event()->UTCtimeStamp());  
	    ntracks[u] = 0; 
	  }
        // find clusters
        else if(spark[u]==0)
          {
            if(i>0)
              lastspark[u]++;
            std::cout << u << "," << lastspark[u] << std::endl;
	    
            // look for clusters
	    
            TH2F* clustimg = (TH2F*)tempimg->Clone("clustimg");
	    
            //old cluster finding algorithm 
            if (!strcmp(conf->getClusterFindingAlgorithm(algo_index),"ci")) 
	      {
		TH2F* baseimage = (TH2F*)tempimg->Clone("baseimage");
		if(eventVersion==1)
		  {
		    clusti = new MaxCamClusterImage(baseimage,d.event()->timeStamp());
		  }
		else
		  {
		    clusti = new MaxCamClusterImage(baseimage,d.event()->UTCtimeStamp());
		  }
		// rebin
		baseimage->Rebin2D(2,2);
		// blur
		baseimage = (TH2F*)MaxCamImageTools::blur(baseimage,1,conf->getBlurAmount());
		
		if (conf->useGainMap() && conf->getClusterFindCIUseGainMapToMerge())
		  {
		    ntracks[u]=  MaxCamImageTools::findClustersGM(baseimage,
								  clusti,
								  conf->getClusterMinSigma(),
								  conf->getClusterMaxSigma(),
								  conf->getClusterMinSize(),
								  conf->getClusterMinDist(),
								  conf->getClusterJoinMinRxyGlobal(),
								  conf->getClusterJoinMinRxyCluster(),
								  conf->getClusterJoinMaxJoinResidual(),
								  conf->getClusterJoinLeastSquaresWeight(),
								  gainmaps[u],
								  conf->getClusterJoinSpacerWidth()
								  );
		  }
		else
		  {
		    ntracks[u]=  MaxCamImageTools::findClustersCI(baseimage,
								  clusti,
                                                          conf->getClusterMinSigma(),
								  conf->getClusterMaxSigma(),
								  conf->getClusterMinSize(),
								  conf->getClusterMinDist());
		  }
		
		if (dbg_level & dbg_cluster) std::cout << "     unsmeared cluster pixel threshold: " << conf->getClusterMinSigma()*rms_nooutliers << std::endl;
		clusti->changeImageWithThreshold(clustimg,conf->getClusterReducedThreshold()*rms_nooutliers);  
		// now we must delete baseimage
		gROOT->Delete("baseimage"); 
		
            }
            
            //Otherwise use a new one
            else
	      {
		clusti = new MaxCamClusterImage(clustimg,d.event()->UTCtimeStamp());
		if(!strcmp(conf->getClusterFindingAlgorithm(algo_index),"seed"))
		  {
		    ntracks[u] = MaxCamImageTools::findClustersGMSeed(clustimg, 
								      clusti,
								      conf->getClusterSeedThreshold(), 
								      conf->getClusterThreshRingRatio(),
								      conf->getClusterMaxWrongProb(),
								      conf->getClusterMinThreshold(),
								      conf->getClusterBlurRadius(),
								      conf->getGaussianBlurAmount(),
								      conf->getClusterNeighborsThresholdForFilling(),
								      conf->getClusterMinNeighborsToKeepPixel(),
								      conf->getClusterMinSizeUnbinned(),
								      conf->getClusterMinDist(),
								      conf->getClusterJoinMinRxyCluster(),
								      conf->getClusterJoinMinRxyGlobal(),
								      conf->getClusterJoinMaxJoinResidual(),
								      conf->getClusterJoinLeastSquaresWeight(),
								      gainmaps[u],
								      conf->getClusterJoinSpacerWidth()); 
		    
		  }
		else if (!strcmp(conf->getClusterFindingAlgorithm(algo_index),"ad"))
		  {
		    ntracks[u] = MaxCamImageTools::findClustersADHysteresisGM(clustimg, clusti, 
									      conf->getClusterFindADK(), 
									      conf->getClusterFindADLambda(), 
									      MaxCamImageTools::TUKEY, 
									      conf->getClusterFindADNIter(), 
									      conf->getClusterFindADGradientBlurAmount(),
									      MaxCamImageTools::SOBEL, 
									      conf->getClusterFindADHighThresh(),
									      conf->getClusterFindADLowThresh(),
									      conf->getClusterNeighborsThresholdForFilling(),
									      conf->getClusterMinSizeUnbinned(),
									      conf->getClusterMinDist(),
									      conf->getClusterJoinMinRxyCluster(),
									      conf->getClusterJoinMinRxyGlobal(),
									      conf->getClusterJoinMaxJoinResidual(),
									      conf->getClusterJoinLeastSquaresWeight(),
									      gainmaps[u],
									      conf->getClusterJoinSpacerWidth());
		  }
		else
		  {
		    ntracks[u] = 0; 
		    std::cout << "NO CLUSTER FINDING ALGORITHM!!!!" << std::endl; 
		  }
		
		clusti->applyRedThreshold(conf->getClusterReducedThreshold() * rms_nooutliers); 
	      }
	    myBench->Show("Cluster finding");
	    
            sumtracks+=ntracks[u];
	    
            std::cout << "Ntracks: " << ntracks[u] << std::endl;
            //Loop over the tracks here

	    myBench->Start("Loop over tracks");
            for(int v=0; v<ntracks[u]; v++)
              {
                if (v>=15) break; 
                vector<BurninEncoded_t> this_track_burnin;
                
                vector<int> pxs_smear = clusti->getCluster(v);
                if (dbg_level & dbg_cluster) std::cout << "     smeared cluster size: " << int(pxs_smear.size()) << std::endl;
                vector<int> pxs_nosmear = clusti->getClusterRed(v);
                if (dbg_level & dbg_cluster) std::cout << "     unsmeared cluster size: " <<  int(pxs_nosmear.size()) << std::endl;
		
                if(!conf->useGainMap())
		  {
		    E[u][v] = clusti->getIntegral(v);
		    EGainMap[u][v] = 0;
		  }
                else
		  {
		    E[u][v] = clusti->getIntegral(v);
		    EGainMap[u][v] = clusti->getIntegralWithGainMap(v,gainmaps[u]);
		  }
		edge[u][v] = clusti->hitsEdge(v);
                skewness[u][v] = clusti->getSkewness(v,phi[u][v]);
                theta[u][v] = 0; //TODO: When we figure out how to get theta, fix this
                switch(conf->getPhiAlgorithm())
		  {
                  case 1: phi[u][v] = clusti->getPhi(v); break;
                  case 2: phi[u][v] = clusti->getPhi2(v); break;
                  case 3: phi[u][v] = clusti->getPhi3(v); break; 
                  case 4: phi[u][v] = clusti->getPhi4(v); break;
                  default: phi[u][v] = 0; 
		    if (dbg_level & dbg_warning) std::cout << "Warning: invalid phi algorithm selected" << std::endl; 
		  }
                double xb,yb,xe,ye;
                switch(conf->getRangeAlgorithm())
		  {                  
                  case 1: range[u][v] = clusti->getLength(v,xb,yb,xe,ye); break;
                  case 2: range[u][v] = clusti->getLength2(v,phi[u][v],1); break;
                  default: range[u][v] = 0;
		    if (dbg_level & dbg_warning) std::cout << "Warning: invalid range algorithm selected" << std::endl; 
		  }
		diffusedrange[u][v]=clusti->getDiffusedLength(v,xb,yb,xe,ye);
		clusti->getEllipseAxes(v,majoraxis[u][v],minoraxis[u][v]);
		
                clusti->getXY(v,x[u][v],y[u][v]);
                
                cout<<"x y "<<x[u][v]<<" "<<y[u][v]<<endl;
                
                cout<<"MINOR AXIS: "<<endl;
                cout<<minoraxis[u][v]<<endl;
                
				xbegin[u][v]=xb;
				xend[u][v]=xe;
                ybegin[u][v]=yb;
				yend[u][v]=ye;
				
                //Get cluster mean, rms, maxpixel
                
                cluster_mean[u][v] = clusti->getMean(v);
                cluster_rms[u][v] = clusti->getRMS(v,cluster_mean[u][v]);
                energy_density[u][v] = clusti->getEnergyDensity(v); 
                npixel[u][v] = clusti->getCluster(v).size(); 
                npixel_red[u][v] = clusti->getClusterRed(v).size(); 
                int maxBin; 
                maxpixel[u][v] = clusti->getMax(v,&maxBin);
                neighbors[u][v] = clusti->getNumNeighbors(v,conf->getClusterMinSigma(),maxBin);
                cygnus_angle[u][v] = clusti->getCygnusAngle(v, conf->getNorthAngle(),
                                                            orientation,
                                                            conf->getLatitude(), conf->getLongitude(),
                                                            phi[u][v],
                                                            theta[u][v]+ TMath::Pi()/2. //TODO: when we figure out theta, remove this argument
                                                            );
		
		
                clusti->getRADec( phi[u][v], theta[u][v] + TMath::Pi()/2., //TODO: ditto as above
				  d.event()->timeStamp(), 
				  conf->getLatitude(), conf->getLongitude(), 
				  conf->getNorthAngle(), orientation, 
				  ra[u][v],dec[u][v],glat[u][v],glon[u][v]); 
		
                for (int m = 0; m <4; m++)
		  { 
		    moments[u][m][v] = clusti->getMoment(v,m+1,phi[u][v], 4,"pixelPerBin"); 
		    transverse_moments[u][m][v] = clusti->getMoment(v,m+1,phi[u][v]+TMath::Pi()/2., 4,"pixelPerBin"); 
		  }
		
                //Loop over last up to N_LOOPBACK events in search of things at the same position
                int pos_i = 0; //position index (since we're using iterators)

		
                list<vector<double *>**>::iterator pos_iter; //iterator for previous positions
                list<vector<vector<vector<BurninEncoded_t> > >* >::iterator burnin_iter; 
                
                for (pos_iter = positions.begin(), burnin_iter = burnin_temp.begin(); 
                     pos_iter != positions.end() ; 
                     pos_iter++, burnin_iter++)
		  {
		    vector<double *> ** test_positions = *pos_iter; 
		    vector<vector<vector<BurninEncoded_t> > > * test_burnin = *burnin_iter;
                    
                    
		    //Loop over the tracks in these events
		    for(unsigned int z = 0; z<test_positions[u]->size(); z++)
		      {
			double * position = (double *)(test_positions[u]->at(z));
			if (position[0]==0.0 && position[1] == 0.0) continue; //ignore events with no tracks
			double x_delta = x[u][v] - position[0];
			double y_delta = y[u][v] - position[1];
			
			if(burnin_test(x_delta, y_delta, conf))
			  {
			    int test_event_index = position_offset + pos_i; 
			    this_track_burnin.push_back(encodeBurnin(z,test_event_index));
			    
			    //Mark this event on the other events record...
			    int our_index = position_offset + positions.size();
			    assert(test_event_index != our_index);
			    
			    test_burnin->at(u)[z].push_back(encodeBurnin(v,our_index));
			  }
		      }    
		    pos_i++; 
		  }//End Loop over last up to N_LOOKBACK events
		
                this_camera_burnin.push_back(this_track_burnin); 
              }//End Loop over tracks
 
              
	    myBench->Show("Loop over tracks");

            /* initialize the rest of the event array */
            for(int v=ntracks[u]; v<15; v++)
	      {
	      
                theta[u][v] = 0;
                phi[u][v]=0;
                E[u][v]=0;
                range[u][v]=0;
                x[u][v]=0;
                y[u][v]=0;
                xbegin[u][v]=0;
                xend[u][v]=0;
                ybegin[u][v]=0;
                yend[u][v]=0;
                edge[u][v]=0;
                skewness[u][v] = 0;
                cluster_rms[u][v]=0;
                cluster_mean[u][v]=0;
                maxpixel[u][v]=0;
                cygnus_angle[u][v]=0;
                neighbors[u][v]=-1;
                energy_density[u][v]=0;
                npixel[u][v]=0;
                dec[u][v] = 0; 
                ra[u][v] = 0; 
                glat[u][v] = 0;  
                glon[u][v] = 0; 
                npixel_red[u][v] = 0;
                for (int m = 0; m< 4; m++) 
		  {
		    moments[u][m][v] = 0; 
		    transverse_moments[u][m][v] = 0; 
		  }
		
	      }
	  } // End if spark[u] == 0
        
        if (algo_index > 0)
	  {

	  	 //spitz added this because it doesnt work with MC files for some reason
	  	if(conf->Bias())    
	    const_cast<TH2F*>((TH2F*)clusti->getImage())->Delete(); 

	    clusti->forgetImage(); 

	  }

        clust->Add(clusti);
        this_event_burnin->push_back(this_camera_burnin);

	
      } //End loop over cameras
    myBench->Show("Loop over cameras");
    
    /* Copy sparkref to list */
    vector<pair<int,int> > * this_sparkref = new vector<pair<int,int> >[ncamera]; 
    for (int u = 0; u < ncamera; u++)
      {
	for (unsigned int v = 0; v < sparkref_running[u].size(); v++)
	  {
	    this_sparkref[u].push_back(sparkref_running[u][v]); 
	  }
	
	assert (this_sparkref[u].size()==sparkref_running[u].size()); 
      } 
    sparkref.push_back(this_sparkref); 
    //Clear temporary images
    if (!fake_image)
      {
	gROOT->Delete("tempimg"); 
	gROOT->Delete("copy");
      }
    
    // fill tree
    myBench->Start("Filling tree");
    
    if(sumtracks == 0 && sparks==0) nempty++;
    
    // if(sumtracks>0 || nempty%100 ==0 || sparks>0)
    // Now just write out every track... 
    { 
      for(unsigned int csn=0; csn<cameraSerialNumber->size(); csn++)
        {
          //std::cout << (*cameraSerialNumber)[csn] << std::endl;
        }
      eventnum=i;
      skimevent = d.event();
      tmpoutfile->cd();
      //tmpskimtree2->Fill();
      tmpskimtree->Fill();
      gROOT->cd();
      cameraSerialNumber->clear();
      
      trigger_groups->SetOwner(kTRUE); 
      trigger_groups->Clear(); 
      clust->SetOwner(kTRUE); 
      clust->Clear(); 
      //Add self to position list, removing first entry if list is already N_LOOKBACK
      vector<double *> ** these_positions = new vector<double *>*[ncamera];
      for (int u = 0; u < ncamera; u++)
        {
          these_positions[u] = new vector<double *>; 
          for (int v = 0; v < ntracks[u]; v++)
	    {
	      if (v >=15) continue; 
	      double * this_pos = new double[2];
	      this_pos[0] = x[u][v];
	      this_pos[1] = y[u][v];
	      these_positions[u]->push_back(this_pos); 
	    }
        } 
      
      myBench->Show("Filling tree");
      
      positions.push_back(these_positions); 
      
      //Store the burnin_info in the buffer until we go past it enough times
      burnin_temp.push_back(this_event_burnin); 
      
      
      //After we start filling the buffer, 
      // it is time to clear the first entry
      // and to write out the first event
      
      if (positions.size() >(unsigned int) conf->getBurninNumEvents() && !fake_image)
	{
	  
	  //Delete the stored positions. 
	  for (int u = 0; u < ncamera; u++)                 
	    {
	      for (unsigned int z = 0; z < positions.front()[u]->size(); z++)
		{
		  delete positions.front()[u]->at(z);
		}
	      delete positions.front()[u];
	    }
	  delete positions.front();
	  positions.pop_front(); 
	  
	  
	  burnin.push_back(burnin_temp.front());
	  burnin_temp.pop_front(); 
	  
	  position_offset++;
	  gROOT->cd(); 
	}
      myBench->Show("Event");
      myBench->Reset();
    }
  }//End loop over events
  
  
  std::cout << "End of loop over events" << std::endl;
  
  //Empty out remaining event buffer. 
  while (burnin_temp.size() > 0)
    {
      burnin.push_back(burnin_temp.front());
      burnin_temp.pop_front(); 
    }
  
  
  //Write out and clean up
  tmpoutfile->cd();
  tmpskimtree->Write();
  //tmpskimtree2->Write();
  if (algo_index == 0)
    {
      biasoutfile->cd();
      biastree->Write();
    }
  // std::cout << sd.getGainMaps()->GetEntries() << std::endl;
  //  std::cout << sd.getGainMap(0)->GetName() << std::endl;
  
  
  //Save configuration
  if (algo_index == 0)
    {
      std::stringstream str; 
      conf->print(str); 
      sd.setConfig(str.str().c_str()); 
      sd.writeConfig(); 
    }
  
  sd.mergeTempTrees(tmpskimtree,&burnin, &sparkref, runnum);    
  
  if (conf->useGainMap())
    {
      for (int i = 0; i < ncamera; i++) 
	{
	  if (gainmaps[i]) gainmaps[i]->Delete(); 
	}
    }
  
  
  
  while(burnin.size() > 0)
    {
      delete burnin.front(); 
      burnin.pop_front(); 
    }
  
  while(sparkref.size() > 0)
    {
      delete sparkref.front(); 
      sparkref.pop_front(); 
    }
  
  delete tmpskimtree;
  //delete tmpskimtree2;
  tmpoutfile->Close();
  
  delete biastree;
  if (algo_index == 0) 
    {
      biasoutfile->Close();
    }
  
  std::cout << "delete_intermediate" << std::endl;
  
  if (delete_intermediate)
    {
      if (dbg_level & dbg_progress) std::cout << "Deleting " << tmpoutfilename << std::endl;
      unlink(tmpoutfilename.Data());
    }
  
  return 0;
}

static const int nreqmod = 0; 

/* Main function */ 

int main(int argn, char ** argv)
{
   int return_val = 0; 

  /* Argument Parsing */ 
   TString rawdatafiles = "files.txt";
   TString keyfilename = "keys.txt";
   TString outdir = "./skim/"; 
   //TString sumdir = "./sum/";
   char * config = 0;
   if(parse_args(argn, argv, &rawdatafiles, &keyfilename, &outdir,&config)) return -1; 
   TString sumdir = outdir;//summary file directory
   sumdir+="sum/";
  // if (dbg_level & dbg_config) config->print(); 
   
   
  /* Handle keys */
   TString key = "skim"; 
   TString reqmod[nreqmod];
   bool pass; 
   DmtpcKeys k(keyfilename, rawdatafiles, key, outdir, nreqmod, reqmod, pass); 
   
   if (!pass) return -1; 
    
   /* Loop through input files */

   std::cout << "k.getNFiles " << k.getNFiles() << std::endl;

   for (int f = 0; f < k.getNFiles(); f++)
   {

     std::cout << "Processing file " << f << std::endl;

      if (dbg_level & dbg_progress) std::cout << k.getFile(f) << std::endl;
                 
      /* Create DMTPC Database and draw out tree */
      DmtpcDataset d; 
      d.openRootFile(k.getRootDirName() + k.getFile(f)); 
        
      TTree * rawtree = k.getBaseTree(f); 
        
      //Extract run number from file 
      d.getEvent(0); 
      int runnum = d.event()->runNumber(); 
      
      //Add friends
      k.addFriends(rawtree, f); 
      
      key =  "skim"; 
      TString outfilename = k.getFile(f).ReplaceAll(".root",key+".root"); 
      DmtpcSkimDataset sd; 
      
      sd.newRootFile(outdir + outfilename); 

      /* Now ready to process this file */
      CleanSkimConfig * conf = skim::preprocess(&d,config); 
      
      for (unsigned int nprocess = 0; nprocess < conf->getNClusterFindingAlgorithms(); nprocess++)
      {
        if(nprocess > 0) key = TString(conf->getClusterFindingAlgorithm(nprocess)); 

        return_val += cleanSkim(d,key,rawtree,runnum,sd,true,conf,nprocess); 
      }
   }
   return return_val;
}

static bool burnin_test( double x_delta, double y_delta, CleanSkimConfig * conf)
{
  double xd = TMath::Abs(x_delta);
  double yd = TMath::Abs(y_delta); 

  double th = conf->getBurninDistanceThresh(); 
  burnin_method_t method= conf->getBurninMethod();

  switch(method)
  {
    case SQUARE:
      return (xd <= th && yd <= th);
    case CIRCLE:
      return (xd*xd <= th*th && yd*yd <= th*th); 
    case CROSS:
      return (xd <= th || yd <= th);
  }
  
  //if we got here... something's wrong
  if (dbg_level & dbg_warning)
  {
    cerr << "Invalid Burnin Test: " << method << std::endl; 
  }
  return false; 
}

static int parse_args( int nargs, char ** args, TString * files, TString * keys, TString * out, char** cfg)
{
  int c = 0; 
  for (int i = 1; i < nargs; i++)
  {
    /* Handle switches */
    if (args[i][0]=='-')
    {
      //Check for debug level
      if(strcmp(args[i],"-d")==0) 
      {
        char * endptr; 
        dbg_level = (u_int32_t) strtol(args[++i], &endptr, 10);
        if (strlen(endptr) > 0) // This means a digit wasn't passed
        {
           if (dbg_level & dbg_error) cerr << "Non number was passed to debug level " << std::endl; 
           return 1; 
        }
      }

      if(strcmp(args[i],"-c")==0) 
      {
        //free(*cfg); 
        *cfg = args[++i]; 
      }
                  
    }
    else
    {
      switch(c++)
      {
        case 0: *files = TString(args[i]); break; 
        case 1: *keys = TString(args[i]); break;
        case 2: *out = TString(args[i]); break; 
        default: 
          if (dbg_level & dbg_error) cerr << "Unexpected Extra Argument. Aborting!" << std::endl;
          return 1; 
      }
    }
  }
  
  return 0; 
}

int fillWaveformVectorsInTObjArray( DmtpcDataset & d, TObjArray * wfvlist ){

    std::cout << "inside fillWaveformVectorsInTObjArray( "
  	    << &d
  	    << " , " 
  	    << wfvlist
  	    << " ) ..." 
  	    << std::endl;

  


  int nch=wfvlist->GetEntries();
  //  std::cout << "nch=" << nch << std::endl;
  for (int ich=0; ich<nch; ++ich){
    

    // point to the first trigger on this channel
    int nwf=d.event()->rawScopeData()->GetEntries();
        std::cout << "nwf=" << nwf << std::endl;

    // if there's no waveforms, don't do anything
    if(!nwf) 
      {
	//wfvlist->Clear();
break;
      }

    // data quality check
    int ntriggers=d.event()->scopeDataInfo(0)->getNTriggers();
    // the total number of waveforms had better be evenly divided by the number of global triggers
    if(nwf%ntriggers) cout << "!!! ntriggers=" 
			   << ntriggers 
			   << " does not divide nwf=" 
			   << nwf 
			   << " without a remainder.  Asserting!" 
			   << endl;
    assert( !(nwf%ntriggers) );
    // the number of triggers in the config file better match the number of channels
    // implied by the information in the raw data file
    if((nwf/ntriggers)!=nch) cout << "!!! nwf=" 
	 		   << nwf 
			   << " divided by ntriggers=" 
			   << ntriggers 
			   << " is " 
			   << (nwf/ntriggers) 
			   << " which does not equal nch="
			   << nch
			   << ".  Asserting!"
			   << endl;
    assert( (nwf/ntriggers)==nch );
    

    //    std::cout << "nwf=" << nwf << std::endl;
    //    std::cout << "nwf/nch*ich=" << (nwf/nch)*ich << std::endl;
    ScopeWaveformData *hwf=d.event()->rawScopeData(nwf/nch*ich);
    //    std::cout << "hwf=" << hwf << std::endl;
    const char *hname = hwf->GetName();
    //    std::cout << "hname=" << hname << std::endl;

    int boardID, triggerID;
    char channelStr;
    sscanf(hname, "scope_%d_%c_%d", &boardID, &channelStr, &triggerID);

    int channelID = (channelStr=='A') ? 0 : 1;
    
    //    std::cout << boardID << '\t' << triggerID << '\t' << channelStr << std::endl;
    
    TString tsTempWaveformVectorClassType(((TObject*)(*wfvlist)[ich])->IsA()->GetName());
    //    std::cout << "tsTempWaveformVectorClassType=" << tsTempWaveformVectorClassType << std::endl;

    if( tsTempWaveformVectorClassType == (TString)"CspWfVector" ){
      
      //      std::cout << "it's a CspWfVector!" << std::endl;
      CspWfVector *wv=static_cast<CspWfVector*> ((*wfvlist)[ich]);
      wv->clear();

      wv->setBoard(boardID);
      wv->setChan(channelID);
      fillCspWfVector(d, wv, ich, nch);
      
    } else if ( tsTempWaveformVectorClassType == (TString)"FastWfVector" ) {
      
      //      std::cout << "it's a FastWfVector!" << std::endl;
      FastWfVector *wv=static_cast<FastWfVector*> ((*wfvlist)[ich]);
      wv->clear();

      wv->setBoard(boardID);
      wv->setChan(channelID);
      fillFastWfVector(d, wv, ich, nch);
      
    } else if ( tsTempWaveformVectorClassType == (TString)"PMTWfVector" ) {
      
      //      std::cout << "it's a PMTWfVector!" << std::endl;
      PMTWfVector *wv=static_cast<PMTWfVector*> ((*wfvlist)[ich]);
      wv->clear();

      wv->setBoard(boardID);
      wv->setChan(channelID);
      fillPMTWfVector(d, wv, ich, nch);
 
    } else {
      
      //      std::cout << "it's an else!" << std::endl;
      WaveformVector *wv=static_cast<WaveformVector*> ((*wfvlist)[ich]);
      wv->clear();

      //      fillWaveformVector(d, wv, ich, conf);
      
    }

  }
  
  //  std::cout << "leaving fillWaveformVectorsInTObjArray( "
  //	    << &d
  //	    << " , " 
  //	    << wfvlist
  //	    << " ) ..." 
  //	    << std::endl;

  return 0;
}

int fillCspWfVector(DmtpcDataset & d, CspWfVector *cspwv, int ich, int nch)
{
  //  std::cout << "inside fillCspWfVector( ... ) ..." << std::endl;

  // general scope info for this event
  int nent = d.event()->scopeData()->GetEntries(); 
  int ntr = nent/nch; 
  //  std::cout << "Nchannels: " << nch << std::endl;
  //  std::cout << "Ntriggers: " << ntr << std::endl; 
  //  std::cout << "Nentries: " << nent << std::endl; 

  // loop over the triggers for this specific channel, and pull the waveform and use it to 
  // fill a SkimWaveform object for each trigger
  for(int itr=0; itr<ntr; ++itr){ 

    TH1F *hwf=d.event()->scopeData(ich,itr);

    CspWaveform cspwf;

    waveform::analysis::analyzeCSP( hwf , cspwf , 0);
    
    // push this SkimWaveform object onto the WaveformVector for this channel
    cspwv->add(cspwf);
  }
  
  //  std::cout << "cspwv->size()=" << cspwv->size() << std::endl;
  //  std::cout << "exiting fillCspWfVector( ... ) ..." << std::endl;
  
  return 0;
}

int fillPMTWfVector(DmtpcDataset & d, PMTWfVector *pmtwv, int ich, int nch)
{
  //  std::cout << "inside fillPMTWfVector( ... ) ..." << std::endl;

  // general scope info for this event
  int nent = d.event()->scopeData()->GetEntries(); 
  int ntr = nent/nch; 
  //  std::cout << "Nchannels: " << nch << std::endl;
  //  std::cout << "Ntriggers: " << ntr << std::endl; 
  //  std::cout << "Nentries: " << nent << std::endl; 

  // loop over the triggers for this specific channel, and pull the waveform and use it to 
  // fill a SkimWaveform object for each trigger
  for(int itr=0; itr<ntr; ++itr){ 

    TH1F *hwf=d.event()->scopeData(ich,itr);
    
    PMTWaveform pmtwf;
    waveform::analysis::analyzePMT( hwf , pmtwf , 1.0);
    
    // push this SkimWaveform object onto the WaveformVector for this channel
    pmtwv->add(pmtwf);
  }
  
  //  std::cout << "pmtwv->size()=" << pmtwv->size() << std::endl;
  //  std::cout << "exiting fillPMTWfVector( ... ) ..." << std::endl;
  
  return 0;
}

int fillFastWfVector(DmtpcDataset & d, FastWfVector *fastwv, int ich, int nch)
{
  //  std::cout << "inside fillFastWfVector( ... ) ..." << std::endl;
  
  // general scope info for this event
  int nent = d.event()->scopeData()->GetEntries(); 
  int ntr = nent/nch; 
  //  std::cout << "Nchannels: " << nch << std::endl;
  //  std::cout << "Ntriggers: " << ntr << std::endl; 
  //  std::cout << "Nentries: " << nent << std::endl; 
  
  // loop over the triggers for this specific channel, and pull the waveform and use it to 
  // fill a SkimWaveform object for each trigger
  for(int itr=0; itr<ntr; ++itr){ 
    
    TH1F *hwf=d.event()->scopeData(ich,itr);
    
    FastWaveform fastwf;
    //    std::cout << "hwf=" << hwf << std::endl;
    //    std::cout << "ich=" << ich << std::endl;
    //    std::cout << "itr=" << itr << std::endl;
    waveform::analysis::analyzeFast( hwf , fastwf );
    //    std::cout << "just got done with analyzeFast( ... ) ... " << std::endl;
    
    // push this SkimWaveform object onto the WaveformVector for this channel
    fastwv->add(fastwf);
  } 
  
  //  std::cout << "fastwv->size()=" << fastwv->size() << std::endl;
  //  std::cout << "exiting fillFastWfVector( ... ) ..." << std::endl;
  
  return 0;
}
