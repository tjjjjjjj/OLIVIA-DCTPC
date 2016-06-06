#include "MaxCamWaveformTools.hh"
#include "TCanvas.h"
#include "TSystem.h"
#include "TH1.h"
#include "TF1.h"
#include "TH1.h"
#include "TString.h"
#include "TMath.h"

#include <iostream>
#include <fstream>

using std::cout;
using std::endl;
using std::cerr;
using std::flush;

#include <vector>
using std::vector;
#include <math.h>
#include "TROOT.h"

#include <stdarg.h> // for variable argument list
#include <stdlib.h> // for variable argument list

// non-class member predicates used for STL vector sorting...
// find a better way to do this, but it is hard to make predicates class
// members and there doesn't seem to be an easy way to sort based on the 
// second element of a list of pairs inline...
bool SortOn2nd(const pair<Double_t, Double_t>& lhs, const pair<Double_t, Double_t>& rhs);
void print_stl_list( list< pair<Double_t,Double_t> >& a );//utility function print lists of ints

ClassImp(MaxCamWaveformTools)


MaxCamWaveformTools::MaxCamWaveformTools() {;}
MaxCamWaveformTools::~MaxCamWaveformTools() {;}


MaxCamWaveformTools::MaxCamWaveformTools(TH1F *wf, bool doInvert,
					 int nbaselinebins, int nrebin,
					 float threshold) {

    _threshold=threshold;

    _wf=wf;
    if (doInvert) invertPolarity();
    evaluate(nbaselinebins,nrebin);
}


void
MaxCamWaveformTools::evaluate(int nbaselinebins, int nrebin) {

   TH1F* wfi = (TH1F*)_wf->Clone("wfinternal");
   //Rebin to smooth out waveform
   wfi->Rebin(nrebin);
   wfi->Scale(1.0/double(nrebin));

   // what is this call for?!?!?
   wfi->GetXaxis()->SetRange(1,wfi->GetNbinsX()-1);

   /////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////
   //
   // General waveform information  
   //
   /////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////
   // the waveform's baseline
   calculateBaseline(wfi,nbaselinebins);
   // the RMS of the waveform's baseline (???)
   calculateBaselineRMS(wfi,nbaselinebins);

   // the waveform's minimum
   int minbin = wfi->GetMinimumBin();
   _wfMinimum = wfi->GetMinimum()-_baseline;
   _wfMinimumTime = wfi->GetBinCenter(wfi->GetMinimumBin());

   // the waveform's maximum 
   int maxbin = wfi->GetMaximumBin();
   _peakHeight = wfi->GetMaximum()-_baseline;
   _peakHeightTime = wfi->GetBinCenter(wfi->GetMaximumBin());

   // the average of the waveform in a 5 bin
   // window centered on the waveform's
   // maximum (???)
   _peakHeightAvg=0;
   for(int ii=maxbin-2; ii<=maxbin+2; ii++)
     {
       _peakHeightAvg+=wfi->GetBinContent(ii);
     } 
   _peakHeightAvg /= 5.0;

   // the bin in which the waveform maximum occurs
   int peakHeightBin =  wfi->GetMaximumBin();

   // another name for the waveform's minimum (???)
   _troughDepth = wfi->GetMinimum()-_baseline;

   /////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////
   //
   // IS THIS EVEN RIGHT?
   //
   /////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////
   // the average time the waveform spends negative
   int zeroBin=wfi->FindBin(0);
   wfi->GetXaxis()->SetRange(1,zeroBin);
   _peakAverageNegTime = wfi->Integral() / zeroBin;

   // the average time the waveform spends positive
   wfi->GetXaxis()->SetRange(zeroBin, wfi->GetNbinsX() );
   _peakAveragePosTime = wfi->Integral() / ( wfi->GetNbinsX() - zeroBin );
   /////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////
   //
   // IS THIS EVEN RIGHT? 
   //
   /////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////

   /////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////
   //
   // Pulse finding...
   //
   /////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////
   // no matter what, treat the maximum and minimum of 
   // the waveform as a negative and positive pulse
   // (ie, ignore thresholds).  Then do some simple pulse
   // finding with a threshold, and keep a list of 
   // what we find, plus useful calculated quantities,
   // like rise time, integral, etc.
   
   {
     /////////////////////////////////////////////////////////
     /////////////////////////////////////////////////////////
     //
     // 0) Simplest Pulse finding... (waveform min/max)
     //
     /////////////////////////////////////////////////////////
     /////////////////////////////////////////////////////////

     // fill a MaxCamPulse object for the waveform minimum
     _minPulse = new MaxCamPulse(minbin);
     fillPulse(wfi,_minPulse); 
   }

   {
     /////////////////////////////////////////////////////////
     /////////////////////////////////////////////////////////
     //
     // 1) More involved Pulse finding... (peak-finding)
     //
     /////////////////////////////////////////////////////////
     /////////////////////////////////////////////////////////
     // construct the STL list of peaks and shunt them
     // into our class' list of MaxCamPulses...
     findAllPeaks(wfi);
   }

   /////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////
   //
   // THIS NEEDS TO GET MOVED INTO THE PULSE CODE
   //
   /////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////
   float ch90=_peakHeight*0.9;
   float ch10=_peakHeight*0.1;
   int   ich90=0, ich10=0, nch=_wf->GetNbinsX()-1;
   for (int ii=peakHeightBin; ii>0; ii--) {
     if (!ich10 && (wfi->GetBinContent(ii)-_baseline)<ch10 && wfi->GetBinContent(ii+10)>wfi->GetBinContent(ii) ) { ich10=ii; }
     if (!ich90 && (wfi->GetBinContent(ii)-_baseline)<ch90 ) { ich90=ii; }
     if (ich10 && ich90) break;
   }
   _riseTime=(ich90-ich10)*wfi->GetBinWidth(1);
   
   _peakHeightTime10=wfi->GetBinCenter(ich10);
   _peakHeightTime90=wfi->GetBinCenter(ich90);
   
   int ich90f=0, ich10f=0;
   for(int ii=maxbin; ii<nch; ii++)
     {
       if(!ich90f &&  (wfi->GetBinContent(ii)-_baseline)<ch90){ich90f=ii;}
       if(!ich10f &&  (wfi->GetBinContent(ii)-_baseline)<ch10){ich10f=ii;}
       if(ich10f && ich90f) break;
     }
   
   if(ich10f==0) {ich10f=nch-1;}
   if(ich90f==0) {ich90f=nch-1;}
   
   int ich90fr=0, ich10fr=0;
   for(int ii=nch-1; ii>maxbin; ii--)
     {
       if(!ich90fr && (wfi->GetBinContent(ii)-_baseline)>ch90){ich90fr=ii;}
       if(!ich10fr && (wfi->GetBinContent(ii)-_baseline)>ch10){ich10fr=ii;}
       if(ich10fr && ich90fr) break;
     }
   
   _fallTime=((ich10f+ich10fr)/2-(ich90f+ich90fr)/2)*wfi->GetBinWidth(1);
   _timeAbove90=((ich90f+ich90fr)/2-ich90)*wfi->GetBinWidth(1);
   _timeAbove10=((ich10f+ich10fr)/2-ich10)*wfi->GetBinWidth(1);
   /////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////
   //
   // THIS NEEDS TO GET MOVED INTO THE PULSE CODE
   //
   /////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////
   
   // what is this call for?!?!?
   wfi->GetXaxis()->SetRange( 1, wfi->GetNbinsX() );
   wfi->Delete();
   // delete wfi;
}


void
MaxCamWaveformTools::invertPolarity() {
    int n=_wf->GetNbinsX();
    for (int i=1; i<=n; i++) _wf->SetBinContent( i, -_wf->GetBinContent(i) );
}




void
MaxCamWaveformTools::print() {

  
  cout << "********************************************************" << endl;

  cout << "** void MaxCamWaveformTools::print(): " << endl;
  
  cout << "*** BASELINE="<<_baseline
       << " +/- " << _baselineRMS
       <<  endl
       << "*** PEAKHEIGHT="<<_peakHeight
       <<  endl
       << "*** time10,90,100="<< _peakHeightTime10
       << ", " << _peakHeightTime90
       << ", " << _peakHeightTime
       <<  endl
       << "*** RISE="<<_riseTime
       << endl
       << "*** FALL="<<_fallTime
       << endl
       << "*** MINIMUM="<<_wfMinimum
       << endl
       << "*** AVERAGE T+/T-=" <<_peakAverageNegTime<<"/"<<_peakAveragePosTime
       << endl;
  
  cout << "*** NUMBER OF PULSES FOUND=" << _pulseList.size() << endl;
  if(_threshold) cout << "*** THRESHOLD=" << _threshold << endl;
  else cout << "*** THRESHOLD (=20.*BASELINE RMS)=" << 20.*_baselineRMS << endl; 

  cout << "********************************************************" << endl;
  
}

// with a capital 'P' because pyroot can't handle a class
// member function whose name is print with a lower case 'p'
void
MaxCamWaveformTools::Print() {
  print();
}

TH1F*
MaxCamWaveformTools::smooth(TH1F* hdata, TH1F* hkernel) {

  int ndata   = hdata->GetNbinsX();
  int nkernel = hkernel->GetNbinsX();
  const int ndatadim = ndata;
  const int nkerneldim = nkernel;
  const int nconvdim = ndata + nkernel - 1;
  double f[ndatadim];
  double g[nkerneldim];
  double h[nconvdim];
  for (int ii=0; ii<ndatadim; ii++) {
    f[ii] = hdata->GetBinContent(ii+1);
  }
  double kernel_integral = hkernel->Integral();
  for (int ii=0; ii<nkerneldim; ii++) {
    g[ii] = hkernel->GetBinContent(ii+1)/kernel_integral;
  }

  //MaxCamWaveformTools::convolve(f, ndata, g, nkernel, h);
  convolve(f, ndata, g, nkernel, h);

  // Construct the output histogram, restricting the range to
  // match that of the input dataset
  TH1F* hconv = new TH1F(*hdata);
  for (int ibin=1; ibin<=hconv->GetNbinsX(); ibin++) {
    int hindex = ibin + (nkerneldim+1)/2-2;
    hconv->SetBinContent(ibin, h[hindex]);
  }
  return hconv;
}

TH1F*
MaxCamWaveformTools::makeGaussKernel(float binWidth, float sigma, float nsigma) {

  float mean  = 0.0;

  float gausXmin = -nsigma*sigma;
  float gausXmax = -gausXmin;
  
  TF1 *fgaus = new TF1("fgaus", "gaus", 
		       2*gausXmin, 2*gausXmax);
		       //-2*nsigma*sigma, 2*nsigma*sigma);
  fgaus->SetParameter(0, 1/TMath::Sqrt(2*TMath::Pi()));  // amplitude
  fgaus->SetParameter(1, mean);  // mu
  fgaus->SetParameter(2, sigma);  //sigma

  int nbins = TMath::Ceil((gausXmax-gausXmin)/binWidth);
  //cout << "gausXmax-gausXmin = " << (gausXmax-gausXmin) << endl;
  //cout << "binWidth = " << binWidth << endl;
  //cout << "nbins = " << nbins << endl;
  
  TH1F *hkernel = new TH1F("hkernel", "hkernel", nbins, gausXmin, gausXmax);
  for (int ii=1; ii<=hkernel->GetNbinsX(); ii++) {
    hkernel->SetBinContent(ii, fgaus->Eval(hkernel->GetBinCenter(ii)));
  }

  return hkernel;
}

void 
MaxCamWaveformTools::convolve(double ff[], int nf, double gg[], int ng, double hh[]) {
  
  int nh = nf + ng - 1;
  for (int ii=0; ii<nh; ii++) {
    hh[ii] = 0;
    int start = TMath::Max(ii-ng+1, 0);
    int stop  = TMath::Min(nf, ii);
    for (int jj=start; jj<=stop; jj++) {
      int kk= ii-jj;
      hh[ii] += ff[jj]*gg[kk];
    }
  }
}



TH1F*
MaxCamWaveformTools::convolve(TH1F* hdata, TH1F* hkernel, bool flip, bool draw) {

  TCanvas *cc=0;
  if (draw) {
    cc = new TCanvas("cc", "cc", 600, 600);
    cc->Divide(2,1); cc->cd(1); hdata->Draw(); hkernel->Draw("same"); 
  }

  TH1F* hkernel_shifted = new TH1F(*hkernel);
  hkernel_shifted->SetName("hkernel_shifted");

  int nbins = hdata->GetNbinsX();
  float xmindata = hdata->GetBinLowEdge(1);
  float xmaxdata = hdata->GetBinLowEdge(nbins)+hdata->GetBinWidth(nbins);

  TH1F *hconv = new TH1F("hconv", "hconv",   // is this correct???
			 2*nbins, 2*xmindata, 2*xmaxdata);

  // loop over all lags (shifts) and compute the overlap integral
  int nbins_kernel = hkernel->GetNbinsX();
  int nsteps = 2*nbins_kernel-1;
  for (int istep=1; istep<=nsteps; istep++) {

    // generate the shifted and flipped kernel
    hkernel_shifted->Reset();
    int oldbin;
    for (int ibin=1; ibin<=istep; ibin++) {
      oldbin = flip ? istep+1-ibin : (nbins_kernel+1) - (istep+1-ibin);
      float val = hkernel->GetBinContent(oldbin);
      hkernel_shifted->SetBinContent(ibin, val);
    }

    // compute the overlap integral
    hconv->SetBinContent(istep,((*hkernel_shifted)*(*hdata)).Integral("width"));
    //cout << ". " << flush;
    if (draw) {
      cc->cd(1); hkernel_shifted->Draw("same");
      cc->cd(2); hconv->Draw();
      cc->Update(); gSystem->Sleep(100);
    }
  }
  //cout << endl;
  return hconv;
}

void 
MaxCamWaveformTools::histOfFunc(TH1F* hin, TF1* fxn) {

  int nbins = hin->GetNbinsX();
  for (int ibin=1; ibin<=nbins; ibin++) {
    hin->SetBinContent(ibin, fxn->Eval(hin->GetBinCenter(ibin)));
  }

  //for (int ibin=1; ibin<=nbins; ibin++) {
  //  cout << "ibin, hin = " << ibin << ", " << hin->GetBinContent(ibin) << endl;
  //}
}

// could not get the variable argument list to work in cint
// but it should work.   see:
// http://root.cern.ch/root/html/cint/limitati.html
//void
//MaxCamWaveformTools::histOfFuncName(TH1F* hin, TString fxnname, ...) {
//
//  fxnname.ToUpper();
//  cout << "fxnname = " << fxnname << endl;
//  va_list listPointer;
//  va_start(listPointer, fxnname);
//  
//  if (fxnname.Contains("GAUS")) {
//    TF1 f1("f1", "gaus", -1, 1);
//    for (int iparam=0; iparam<3; iparam++) {
//      f1.SetParameter(iparam, va_arg(listPointer, double));
//    }
//    histOfFunc(hin, &f1);
//  } else {
//    cout << "Unrecognized function:  " << fxnname << endl;
//    cout << " returning without modifying the input histogram" << endl;
//  }
//
//  return;
//}

TGraph*
MaxCamWaveformTools::graphOfHist(TH1* h1, TString grName, bool flip) {

  TGraph *gr = new TGraph(); int ipt=0;  gr->SetName(grName); 
  //gr->SetPoint(ipt++, 0, 0);
  int nbins = h1->GetNbinsX();
  for (int ii=1; ii<=nbins; ii++) {
    float histy = h1->GetBinContent(ii);
    float histxlow = h1->GetBinLowEdge(ii);
    float histxhi  = histxlow+h1->GetBinWidth(ii);
    float x1, x2, y1, y2;
    if (flip) {
      x1 = histy; x2 = histy;
      y1 = histxlow; y2=histxhi;
    } else {
      x1 = histxlow; x2 = histxhi;
      y1 = histy; y2 = histy;
    }
    gr->SetPoint(ipt++, x1, y1); gr->SetPoint(ipt++, x2, y2);
  }
  //gr->SetPoint(ipt++, 0, h1->GetBinLowEdge(nbins)+h1->GetBinWidth(nbins));

  // preserve line style, color
  gr->SetLineWidth(h1->GetLineWidth());
  gr->SetLineColor(h1->GetLineColor());

  return gr;
}


Double_t
MaxCamWaveformTools::calculateBaseline(TH1F* wfi, int nbaselinebins) {
  //Calculate baseline statistics
  //int nbaselinebins = 25; 
  _baseline=0;
  for (int ii=1; ii<nbaselinebins+1; ii++) {
    _baseline += wfi->GetBinContent(ii);
  }
  _baseline /= double(nbaselinebins);
  
  return _baseline;
} 

Double_t
MaxCamWaveformTools::calculateBaselineRMS(TH1F* wfi, int nbaselinebins) {
  // requires that the baseline has already
  // been calculated
  _baselineRMS=0;
  for (int ii=1; ii<nbaselinebins+1; ii++) {
    _baselineRMS += wfi->GetBinContent(ii)*wfi->GetBinContent(ii);
  }
  _baselineRMS /= ((Double_t)nbaselinebins);
  //  cout << "_baselineRMS=" << _baselineRMS << endl;
  //  cout << "_baseline*_baseline=" << _baseline*_baseline << endl;
  // what is the point of this?
  //  _baselineRMS = sqrt(_baselineRMS - _baseline*_baseline );
  _baselineRMS=sqrt(_baselineRMS);

  return _baselineRMS;
}

 Double_t
 MaxCamWaveformTools::calculateIntegral(TH1F* wfi, Int_t startBin, Int_t endBin) {

   Double_t integral=0.;
   //   cout << "startBin=" << startBin << endl;
   for (int ii=startBin; ii<(endBin+1); ii++) {
     //     cout << "ii=" 
     //	  << ii 
     //	  << " (wfi->GetBinContent(ii)-_baseline)="
     //	  << (wfi->GetBinContent(ii)-_baseline)
     //	  << " wfi->GetBinWidth(ii)=" 
     //	  << wfi->GetBinWidth(ii)
     //	  << endl;

     integral+=(wfi->GetBinContent(ii)-_baseline)*wfi->GetBinWidth(ii);
   }
   //   cout << "endBin=" << endBin << endl;

   return integral;
} 

void MaxCamWaveformTools::findAllPeaks( TH1F* wfi ){
  // print out play-by-play for debugging?
  Bool_t debug=false;

  if(debug) cout << "Inside findAllPeaks(...) ..." << endl;
  
  // going to do most of our computations on these things
  // as STL lists because they're easy to sort
  // migrate to TClonesArrays someday?
  list< pair<Double_t,Double_t> > wf_vordered;
  
  // loop over waveform and get (t,V) points
  for (Int_t ii=1; ii<wfi->GetNbinsX()+1; ++ii){

    if(debug) cout << "("
		   << wfi->GetBinCenter(ii)
		   << ","
		   << wfi->GetBinContent(ii)
		   << ") , ";
    
    wf_vordered.push_back( make_pair( wfi->GetBinCenter(ii),
				      wfi->GetBinContent(ii) ) );
  }

  // sort the waveform by voltage, low to high
  wf_vordered.sort(SortOn2nd);
  //  print_stl_list(wf_vordered);

  // a sortable list to hold this channel's pulses' start and end times
  list< pair<Double_t,Double_t> > foundPulseList;
  list< pair<Double_t,Double_t> > foundPulseStartAndEndTimesList;
  
  // the peak threshold (anything below this value, both positive and
  // negative is just ignored...
  Double_t threshold_including_baseline=_baseline-20.*_baselineRMS;
  if(debug) cout << "threshold_including_baseline=" 
		 << threshold_including_baseline
		 << " _baseline=" 
		 << _baseline
		 << " -20.*_baselineRMS="
		 << -20.*_baselineRMS
		 << endl;
  
  if(_threshold) threshold_including_baseline=_baseline-_threshold;
  if(debug) cout << "threshold_including_baseline=" << threshold_including_baseline << endl;
  if(debug) cout << "_threshold=" << _threshold << endl;

  // loop through the voltage sorted list popping off maxima until we go below threshold
  for( list< pair<Double_t,Double_t> >::iterator iter = wf_vordered.begin();
       iter != wf_vordered.end(); ++iter ){
    
    Double_t this_time=(*iter).first;
    Double_t this_voltage=(*iter).second;

    // if the voltage is greater than the threshold
    if( this_voltage<threshold_including_baseline ){

      if(debug) cout << "this_time=" << this_time << " this_voltage=" << this_voltage << endl;

      // we want to make sure that this is the only contribution coming from this pulse
      bool inPreviousPulse=false;
      
      if(debug) cout << "foundPulseStartAndEndTimesList.size()=" << foundPulseStartAndEndTimesList.size() << endl;
      
      // loop over found pulse time widths and see if this time falls into 
      // any previously found pulses
      for( list< pair<Double_t,Double_t> >::iterator foundPulseIter = foundPulseStartAndEndTimesList.begin();
	   foundPulseIter != foundPulseStartAndEndTimesList.end(); ++foundPulseIter ){
	Double_t thisFoundPulseStartTime=(*foundPulseIter).first;
	Double_t thisFoundPulseEndTime=(*foundPulseIter).second;
	
	if(debug) cout << "> (" << thisFoundPulseStartTime << "," << thisFoundPulseEndTime << ")" << endl;

	if( ( this_time<=thisFoundPulseEndTime) &&
	    ( this_time>=thisFoundPulseStartTime) ){
	  inPreviousPulse=true;
	}
      } // end loop to see if this point falls into any previous found pulses
      
      if(debug) cout << "inPreviousPulse=" << inPreviousPulse << endl;

      // not found in previous pulses; store
      if(!inPreviousPulse){

	// so we found a pulse above our desired threshold and not in the time window of
	// any previous pulses that we recorded a maximum from.  We want to make sure 
	// this is a local maximum and put the time window of this pulse into the 
	// list of times to exclude looking in from now on
	bool isLocalMinimum=true;
	
	// loop over the full waveform, and find the peak in the data
	for (Int_t ii=1; ii<wfi->GetNbinsX()+1; ++ii){

	  if(debug) cout << "wfi->GetBinCenter(ii)=" 
			 << wfi->GetBinCenter(ii)
			 << " this_time=" 
			 << this_time
			 << endl;

	  if(debug) cout << "( wfi->GetBinCenter(ii)==this_time )=" << ( wfi->GetBinCenter(ii)==this_time ) << endl;
	  if(debug) cout << "TMath::Abs( wfi->GetBinCenter(ii)-this_time )=" << TMath::Abs( wfi->GetBinCenter(ii)-this_time ) << endl;
	  // a less than only because for some reason, in python
	  // these numbers can be VERY slightly different.
	  // the right thing to do is to keep things in terms of 
	  // an index but that will take some restructuring...
	  if( TMath::Abs(wfi->GetBinCenter(ii)-this_time)<1.e-20 ){
	    // this is the point in the TH1F waveform
	    // that corresponds to the point we're
	    // considering as a peak minimum 
	    // candidate

	    //	    if(debug) cout << "(q+1)*TIME_PER_PT=" << (q+1)*TIME_PER_PT << " this_amplitude=" << this_amplitude << endl;
	    //	    if(debug) cout << "-(raw_data[q]-baseline)=" << -(raw_data[q]-baseline) << endl;
	    // the start time and end time of the pulse that this voltage point is in
	    Double_t startOfPulse, endOfPulse;
	    
	    if(debug) cout << "in this loop?" << endl;
	    if(debug) cout << "wfi->GetBinCenter("
			   << ii 
			   << ")="  
			   << wfi->GetBinCenter(ii) 
			   << " wfi->GetBinContent(" 
			   << ii 
			   << ")=" 
			   << wfi->GetBinContent(ii)
			   << " | this_time="
			   << this_time
			   << " this_voltage=" 
			   << this_voltage 
			   << endl;
	    
	    // find the pulse start
	    startOfPulse=findPulseStart(wfi,ii);
	    if(startOfPulse==-1) isLocalMinimum=false;
	    
	    // find the pulse end
	    endOfPulse=findPulseEnd(wfi,ii);
	    if(endOfPulse==-1) isLocalMinimum=false;
	    
	    if(debug) cout << "isLocalMinimum=" 
			   << isLocalMinimum 
			   << endl;
	    
	    if(isLocalMinimum){
	      if(debug) cout << "startOfPulse=" << startOfPulse << " endOfPulse=" << endOfPulse << endl;
	      // SHOULD PUT SOME CODE IN HERE TO FIGURE OUT WHAT THE PULSE *SHOULD* BE HERE
	      // SO WE ONLY HAVE TO DO THIS PART ONCE...
	      foundPulseStartAndEndTimesList.push_back( make_pair( startOfPulse,
								   endOfPulse ) );
	    }
	  }
	}
	
	// so this minimum is below threshold, isn't in any previous pulses, and is a local maximum
	// record it
	if(isLocalMinimum){
	  foundPulseList.push_back( make_pair( this_voltage, this_time ) );
	}
      }
    }
  }
  
  if (debug) cout << "MaxCamWaveformTools::findAllPeaks(): Found " << foundPulseList.size() << " negative pulses..." << endl;
  
  // loop over the found peaks and shunt them into MaxCamPulse objects...
  for(list< pair<Double_t,Double_t> >::iterator foundPulseIter = foundPulseList.begin(); foundPulseIter!=foundPulseList.end(); foundPulseIter++){
    //    cout << '(' << (*foundPulseIter).first << "," << (*foundPulseIter).second << ')';
    MaxCamPulse* this_pulse=new MaxCamPulse(wfi->FindBin((*foundPulseIter).second));
    //    cout << "this_pulse->getBin()=" << this_pulse->getBin() << endl;
    fillPulse(wfi,this_pulse);
    if(debug) this_pulse->print();
    _pulseList.push_back(*this_pulse);
  }

  if(debug) cout << "_pulseList.size()=" << _pulseList.size() << endl;

  // erase the lists just in case
  wf_vordered.erase(wf_vordered.begin(),wf_vordered.end());
  foundPulseStartAndEndTimesList.erase(foundPulseStartAndEndTimesList.begin(),foundPulseStartAndEndTimesList.end());
  foundPulseList.erase(foundPulseList.begin(),foundPulseList.end());

  if(debug) cout << "Leaving findAllPeaks(...) ..." << endl;
}

void MaxCamWaveformTools::fillPulse(TH1F* wfi, MaxCamPulse* pulse){
  
  //  cout << "in void MaxCamWaveformTools::FillPulse(MaxCamPulse* pulse)..." << endl;
  //  cout << "pulse->getBin()=" << pulse->getBin() << endl;

  Int_t nbin=((Int_t)pulse->getBin());
  
  // get the baseline subtracted pulse height
  pulse->setPulseHeight( ( wfi->GetBinContent( nbin ) -
			   _baseline ) );

  // get the time at which the pulse maximum occurs
  pulse->setPulseHeightTime(wfi->GetBinCenter( nbin ));

  pulse->setPulseStartTime(findPulseStart(wfi,nbin));
  Int_t startBin=wfi->FindBin(pulse->getPulseStartTime());
  pulse->setPulseStartBin(startBin);
  
  Double_t pulseEndTime=findPulseEnd(wfi,nbin);
  pulse->setPulseEndTime(pulseEndTime);
  Int_t endBin=wfi->FindBin(pulse->getPulseEndTime());
  pulse->setPulseEndBin(endBin);

  Double_t pulseIntegral=calculateIntegral(wfi,startBin,endBin);
  pulse->setPulseIntegral(pulseIntegral);

  //  cout << "leaving void MaxCamWaveformTools::FillPulse(MaxCamPulse* pulse)..." << endl;
}
   
// returns -1 if isLocalMinimum is false
Double_t MaxCamWaveformTools::findPulseStart(TH1F* wfi, Int_t pulseBin, Bool_t debug ){
  
  Bool_t isLocalMinimum=true;
  Double_t startOfPulse=-1;
  
  Double_t pulse_max_voltage=wfi->GetBinContent(pulseBin);
  Double_t pulse_max_time=wfi->GetBinCenter(pulseBin);
  
  // loop, inside the pulse, behind our candidate
  if(debug) cout << "-----------------------------------------------------" << endl;
  if(debug) cout << "| looping behind the candidate..." << endl;
  for(Int_t s=(pulseBin-1); s>0; --s){
    
    if(debug) cout << "wfi->GetBinCenter("
	 	   << s 
		   << ")=" 
		   << wfi->GetBinCenter(s) 
		   << " wfi->GetBinContent(" 
		   << s 
		   << ")=" 
		   << wfi->GetBinContent(s)
		   << endl;
    
    Double_t thisPointVoltage=wfi->GetBinContent(s);
    
    if((s==(pulseBin-1))&&(thisPointVoltage<pulse_max_voltage)){
      // the point right before this minimum candidate is 
      // higher.  So it's not a local minimum
      isLocalMinimum=false;
      
      Double_t delta=fabs(thisPointVoltage-pulse_max_voltage);
      // a failure mode I don't understand that occurs if two points close to 
      // each other are the same
      // THIS IS NOT THE RIGHT WAY TO HANDLE THIS; NEED TO KEEP GOING LOWER 
      // IF THE ONE BELOW THIS ONE IS THE SAME UNTIL YOU FIND ONE THAT'S
      // DIFFERENT AND USE THAT INFORMATION TO DECIDE IF THIS IS A LOCAL
      // MINIMUM OR NOT
      // this number is chosen because it is much smaller than the voltage
      // digitization number
      if(delta<1e-6) isLocalMinimum=true;
    }  
    
    // stop looping if we cross into the baseline
    if( thisPointVoltage > _baseline ){

      break;
    } 

    // only pass a non-negative value back if this is  local minimum
    startOfPulse=wfi->GetBinCenter(s);
    if(debug) cout << "startOfPulse=" 
		   << startOfPulse 
		   << endl;
  }
  if(debug) cout << "| done looping behind of the candidate..." << endl;
  if(debug) cout << "-----------------------------------------------------" << endl;
  if(debug) cout << endl;

  if(!isLocalMinimum) return -1;
  return startOfPulse;
}

// returns -1 if isLocalMinimum is false
Double_t MaxCamWaveformTools::findPulseEnd(TH1F* wfi, Int_t pulseBin, Bool_t debug ){

  Bool_t isLocalMinimum=true;
  Double_t endOfPulse=-1;

  Double_t pulse_max_voltage=wfi->GetBinContent(pulseBin);
  Double_t pulse_max_time=wfi->GetBinCenter(pulseBin);
  
  // loop, inside the pulse, ahead our candidate
  if(debug) cout << "-----------------------------------------------------" << endl;
  if(debug) cout << "| looping ahead the candidate..." << endl;
  for(Int_t s=(pulseBin+1); s<wfi->GetNbinsX(); ++s){
    
    if(debug) cout << "wfi->GetBinCenter("
		   << s 
		   << ")=" 
		   << wfi->GetBinCenter(s) 
		   << " wfi->GetBinContent(" 
		   << s 
		   << ")=" 
		   << wfi->GetBinContent(s)
		   << " wfi->GetNbinsX()="
		   << wfi->GetNbinsX()
		   << endl;
    
    Double_t thisPointVoltage=wfi->GetBinContent(s);
    
    if((s==(pulseBin+1))&&(thisPointVoltage<pulse_max_voltage)){
      // the point right before this minimum candidate is 
      // higher.  So it's not a local minimum
      isLocalMinimum=false;
      
      Double_t delta=fabs(thisPointVoltage-pulse_max_voltage);
      // a failure mode I don't understand that occurs if two points close to 
      // each other are the same
      // THIS IS NOT THE RIGHT WAY TO HANDLE THIS; NEED TO KEEP GOING LOWER 
      // IF THE ONE BELOW THIS ONE IS THE SAME UNTIL YOU FIND ONE THAT'S
      // DIFFERENT AND USE THAT INFORMATION TO DECIDE IF THIS IS A LOCAL
      // MINIMUM OR NOT
      // this number is chosen because it is much smaller than the voltage
      // digitization number
      if(delta<1e-6) isLocalMinimum=true;
    }  
    
    // stop looping if we cross into the baseline
    if( thisPointVoltage > _baseline ){
      break;
    }
    
    endOfPulse=wfi->GetBinCenter(s);
    if(debug) cout << "endOfPulse=" 
		   << endOfPulse 
		   << endl;
  }
  if(debug) cout << "| done looping ahead of the candidate..." << endl;
  if(debug) cout << "-----------------------------------------------------" << endl;
  if(debug) cout << endl;

  if(!isLocalMinimum) return -1;
  return endOfPulse;
}  

bool SortOn2nd(const pair<Double_t, Double_t>& lhs, const pair<Double_t, Double_t>& rhs) 
{ 
  return lhs.second < rhs.second; 
} 

void print_stl_list( list< pair<Double_t,Double_t> >& a)
{
  cout << "a.size()=" << a.size() << endl;

  list< pair<Double_t,Double_t> >::iterator iter = a.begin();

  for(;iter!=a.end();iter++){
    cout << '(' << (*iter).first << "," << (*iter).second << ')';
  }
  cout << endl;
  cout << "----------------"<<endl;
}

