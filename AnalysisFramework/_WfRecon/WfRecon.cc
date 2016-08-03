#include "TFile.h"
#include <TObject.h>
#include <TMath.h>
#include <cmath>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <TString.h>
#include <string.h>
#include "WfRecon.hh"
#include <stdlib.h>

ClassImp(WfRecon);

using namespace std;

WfRecon::WfRecon(){}

WfRecon::WfRecon(TString knobfile_)
{
  knobfile = knobfile_;
  TurnKnobs();
}

/*
WfRecon::WfRecon(TString rootfile_,TString braggfile_, TString knobfile_)
{
  knobfile = knobfile_;
  TurnKnobs();
  rootfile = rootfile_;
  braggfile = braggfile_;
}
*/

WfRecon::~WfRecon(){}


//There has to be a better way...
void WfRecon::TurnKnobs(){
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

double WfRecon::findRMS(int xmin, int xmax, TH1D* curve, TH1D* hist1, TH1D* hist2)
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

double WfRecon::maxdevsearch(int xmin, int xmax, TH1D* bragghist3, TH1D* bragghist4, TH1D* mfc2)
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

int WfRecon::terminatorsearch(TH1D* bragghist3, TH1D* bragghist4, TH1D* mfc2, int xmin, int xmax, double threshold, vector<int>& vec)
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

int WfRecon::terminatorsearchbackup(TH1D* bragghist3, TH1D* bragghist4, TH1D* mfc2, int xmin, int xmax)
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


int WfRecon::terminatorsearch2(TH1D* bragghist3, TH1D* bragghist4, TH1D* mfc2, int xmin, int xmax, double threshold)
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


int WfRecon::terminatorsearch3(TH1D* bragghist3, TH1D* bragghist4, int xmin, int xmax)
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


void WfRecon::plotbragg(TH1D* bragghist3, TH1D* bragghist2, double xscaletrue, double xscale, double yscale, int xoff)
{
  for (int x = 1; x < 16384 - xoff; x++)
    bragghist3->SetBinContent(x-1 + xoff,histspline(bragghist2,0,0,x,1/xscaletrue,yscale / xscale));
}

void WfRecon::plotbraggrev(TH1D* bragghist4, TH1D* bragghist2, double xscaletrue, double xscale, double yscale, int xoff)
{
  for (int x = xoff; x > 0; x--)
    bragghist4->SetBinContent(x-1,histspline(bragghist2,0,0,xoff - x + 1,1/xscaletrue,yscale / xscale));
}
//dunno why these offsets work but they do
  
double WfRecon::integrate(TH1D* hist, int imin, int imax)
{
  double sum = 0;
  for (int i = imin+1; i < imax; i++)
    sum += hist->GetBinContent(i);
  return sum;
}


int WfRecon::maxsearch(TH1D* hist, int imin, int imax)
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
int WfRecon::halfsearch(TH1D* hist, int imin, int imax, double maxval)
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


int WfRecon::absmin(TH1D* hist, int i, int range)
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
  
double WfRecon::deriv(TH1D* hist, int n, int k, double h)
{
  if (deriv_order == 2) return (hist->GetBinContent(n+k) - hist->GetBinContent(n-k)) / (2 * h * k);
  else if (deriv_order == 4) return (8*hist->GetBinContent(n+k)-8*hist->GetBinContent(n-k)-hist->GetBinContent(n+2*k)+hist->GetBinContent(n-2*k))/(12*h*k);
  else
    {
      deriv_order = 2;
      return deriv(hist,n,k,h);
    }
}

void WfRecon::deconvolve(TH1D* hist, double SD, int order)
{
  TH1D* hist2 = new TH1D("","",16384,0,16384);
  for (int i = 2000; i < 10000; i++)
    {
      if (order == 2)
	hist2->SetBinContent(i,hist->GetBinContent(i) - 0.25 * pow(SD*dt,2) * deriv2(hist,i,d2_step,dt));
      else if (order == 4)
	hist2->SetBinContent(i,hist->GetBinContent(i) - 0.25 * pow(SD*dt,2) * deriv2(hist,i,d2_step,dt) + .125 * pow(SD*dt,4) * deriv4(hist,i,d4_step,dt));
      /*else if (order == 6)
	hist2->SetBinContent(i,hist->GetBinContent(i) - 0.25 * pow(SD*dt,2) * deriv2(hist,i,d2_step,dt) + .125 * pow(SD*dt,4) * deriv4(hist,i,d4_step,dt)) - 64./6 * pow(SD*dt,6)* deriv6(hist,i,d6_step,dt);
      else if (order == 8)
      hist2->SetBinContent(i,hist->GetBinContent(i) - 0.25 * pow(SD*dt,2) * deriv2(hist,i,d2_step,dt) + .125 * pow(SD*dt,4) * deriv4(hist,i,d4_step,dt)) - 64./6 * pow(SD*dt,6)* deriv6(hist,i,d6_step,dt) + 64./6 * pow(SD*dt,8) * deriv8(hist,i,d8_step,dt);*/ 

/* Math is likely wrong. Double-check the coefficients before un-commenting these options. */
      else
	cout << "invalid value for variable 'order' (must be 2 or 4.)" << endl;
    }
  for (int i = 0; i < 16384; i++)
    hist->SetBinContent(i,0);
  for (int i = 2100; i < 10000; i++)
    hist->SetBinContent(i,hist2->GetBinContent(i));
}

double WfRecon::deriv2(TH1D* hist, int n, int k, double h)
{
  return (hist->GetBinContent(n+k) + hist->GetBinContent(n-k) - 2*hist->GetBinContent(n))/(h*k*h*k);
}
double WfRecon::deriv4(TH1D* hist, int n, int k, double h)
{
  return (hist->GetBinContent(n+2*k) - 4* hist->GetBinContent(n+k) + hist->GetBinContent(n-2*k)-4* hist->GetBinContent(n-k) +6*hist->GetBinContent(n))/(pow(h*k,4));
}
double WfRecon::deriv6(TH1D* hist, int n, int k, double h)
{
  return (hist->GetBinContent(n+3*k)-6*hist->GetBinContent(n+2*k)+15*hist->GetBinContent(n+k)-20*hist->GetBinContent(n)+15*hist->GetBinContent(n-k)-6*hist->GetBinContent(n-2*k)+hist->GetBinContent(n-3*k))/(pow(h*k,6));
}
double WfRecon::deriv8(TH1D* hist, int n, int k, double h)
{
  return (hist->GetBinContent(n+4*k) - 8 * hist->GetBinContent(n + 3*k) + 28 * hist->GetBinContent(n + 2*k) - 56 * hist->GetBinContent(n + k) + 70 * hist->GetBinContent(n) - 56 * hist->GetBinContent(n - k) + 28 * hist->GetBinContent(n-2*k) - 8 * hist->GetBinContent(n-3*k) + hist->GetBinContent(n-4*k))/(pow(h*k,8));
}


int WfRecon::minfinder(TH1D* wfder)
{
  double thresh = .2 * wfder->GetBinContent(wfder->GetMaximumBin());
  TH1D* dummy = new TH1D("","",16384,0,16384);
  for (int i = 0; i < 16384; i++)
    if (wfder->GetBinContent(i) > thresh)
      dummy->SetBinContent(i,1e6);
  return minfinder(dummy,wfder);
}
int WfRecon::minfinder(TH1D* wf, TH1D* wfder)
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

void WfRecon::minprepare(TH1D* wf, TH1D* minfinder_curve)
{
  double k = 1/Td;
    
  TH1D* wfder = new TH1D("","",16384,0,16384);
    
  for(int j = 100;j < 16300; j++)
    wfder->SetBinContent(j-1,k*(wf->GetBinContent(j)-vertoff)+deriv(wf,j,deriv_step,dt));
  //small offset accounts for 8ns rise time

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

double WfRecon::average(int j, vector<double>& vec, TH1D* hist)
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

double WfRecon::median(vector<double>& vec){
  if (vec.size()==0){return 0;}
  sort(vec.begin(),vec.end());
  return vec.size()%2 == 1 ? vec.at((vec.size()-1) / 2) : .5*(vec.at(vec.size()/2)+vec.at(vec.size()/2-1));
}

 double WfRecon::normaldist(double peak, double x, double SD)
{return exp((x-peak)*(peak-x)/(2*SD*SD)) / (SD * pow(3.141592,0.5));}

void WfRecon::normalize(vector<double>& vec)
{
  double sum = 0;
  for (int i = 0; i < vec.size(); i++)
    sum += vec.at(i);
  if (sum == 0) sum = 1;
  for (int i = 0; i < vec.size(); i++)
    vec.at(i) /= sum;
}
 
void WfRecon::makekernel(vector<double>& vec, double SD, double acc)
{
  if (SD == 0) SD = 1.e-69;
  int size = int(acc*SD+1);
  for (int i = -size; i < size + 1; i++)
    vec.push_back((4*normaldist(0.,double(i),SD)+normaldist(0.,double(i)+.5,SD)+normaldist(0.,double(i)-.5,SD))/6.);
  normalize(vec);
}

void WfRecon::smooth(TH1D* hist, double SD, double acc)
{
  TH1D* hist2 = new TH1D("","",16384,0,16384);
  vector<double> vec;
  makekernel(vec,SD,acc);

  for(int j=2000; j < 10200; j++)
    hist2->SetBinContent(j,average(j,vec,hist));
  for(int j=1; j < 16384; j++)
    hist->SetBinContent(j,hist2->GetBinContent(j));
}

void WfRecon::correct(TH1D* hist, double SD, double acc, double factor)
{
  TH1D* hist2 = new TH1D("","",16384,0,16384);
  vector<double> vec;
  makekernel(vec,SD,acc);

  for(int j=2000; j < 10200; j++)
    hist2->SetBinContent(j,(1+factor)*hist->GetBinContent(j)-factor*average(j,vec,hist));
  for(int j=1; j < 16384; j++)
    hist->SetBinContent(j,hist2->GetBinContent(j));
}

void WfRecon::shifthist(TH1D* hist, int plot_offset)
{
  for (int i = 16383; i > 0; i--)
    if (i > plot_offset)
      hist->SetBinContent(i,hist->GetBinContent(i+plot_offset));
    else
      hist->SetBinContent(i,0);
}

double WfRecon::histspline(TH1D* hist, int xmin, int xoff, int x, double xscale, double yscale)
{
  double true_x = xscale * (x - xmin) + xoff;
  int int_x = int(true_x);
  double weight1 = fabs(true_x - int_x);
  double weight0 = 1 - weight1;
  if (int_x < 0) int_x = -int_x;
  double value =  weight0 * hist->GetBinContent(int_x) + weight1 * hist->GetBinContent(int_x + 1);
  return value * yscale;
}

void WfRecon::analyze()
{
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
  //int jpeaksearch = int(median(minima) + .5);

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
