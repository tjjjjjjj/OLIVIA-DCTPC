#ifndef DMTPC_Waveform_TOOLS_HH
#define DMTPC_Waveform_TOOLS_HH

#include "TH1F.h"
#include "TGraph.h"
#include "DmtpcPulse.hh"
#include <list>
#include <algorithm>
using std::pair;
using std::make_pair;
using std::list;
using std::copy;

#include <iostream>
using std::cout;
using std::endl;

#include <vector>
using std::vector;

/** Determines reconstructed quantities for a waveform stored in a TH1 */

class DmtpcWaveformTools :public TObject {

public:
   
   /** constructor; runs evaluate
       \param wf TH1 containing the waveform
       \param doInvert true to invert the waveform, false to leave it as is
       \param nbaselinebins number of bins to calculate the baseline from
       \param nrebin rebinning factor for waveform; higher numbers may reduce error, but also reduce precision
   */
    DmtpcWaveformTools(TH1F *wf, bool doInvert=false,
			int nbaselinebins=25, int nrebin=1,
			float threshold=0);

    DmtpcWaveformTools();
    ~DmtpcWaveformTools();
    
   /** method to evalute and populate the reconstructed quantities. Descriptions of the reconstructed quantities are given at their get functions
       \param nbaselinebins number of bins to calculate the baseline from
       \param nrebin rebinning factor for waveform; higher numbers may reduce error, but also reduce precision
    */
    void evaluate(int nbaselinebins, int nrebin);

   /** calculated as the average of the first nbaselinebins */
    float getBaseline() { return _baseline; }

   /** calculated as the rms of the first nbaselinebins */
    float getBaselineRMS() { return _baseline; }

   /** calculated as the maximum of TH1, minus the baseline*/
    float getPeakHeight() { return _peakHeight; }

   /** calculated as the center of the time bin in which the peak height occurs*/
    float getPeakHeightTime() { return _peakHeightTime; }

   /** calculated as the minimum of the TH1, minus the baseline */
    float getMinumum() { return _wfMinimum; }

   /** calculated as the time to go from 10% to 90% of the peak height */
    float getRiseTime() { return _riseTime; }

   /** same as minimum */
    float getTroughDepth() { return _troughDepth; }

   /** calculated as the time to fall from 90% to 10% of peak. If the waveform does not get down to 10% of peak, the fall time is calculated as 90% of peak to the end of readout */
    float getFallTime() { return _fallTime; }

   /** same as Minimum */
    float getWfMinimum() { return _wfMinimum; }

   
    float getPeakAveragePosTime() { return _peakAveragePosTime; }
    float getPeakAverageNegTime() { return _peakAverageNegTime; }
    
   /** prints quantities for the waveform */
    void print();
    void Print();
   /** takes y->-y for all bins for the TH1 */
    void invertPolarity();

    // routines for computing pulse characteristics
    Double_t calculateBaseline(TH1F* wfi, int nbaselinebins);
    Double_t calculateBaselineRMS(TH1F* wfi, int nbaselinebins);
    Double_t calculateIntegral(TH1F* wfi, Int_t startBin, Int_t endBin);
  
    // routines for peak-finding
    void findAllPeaks(TH1F* wfi);
    void fillPulse(TH1F* wfi, DmtpcPulse* pulse);
    void getCopyOfPulseList(list<DmtpcPulse>& list){ list.resize(_pulseList.size()); copy(_pulseList.begin(),_pulseList.end(),list.begin()); }
    Int_t getSizeOfPulseList(void){ return _pulseList.size(); }

    DmtpcPulse returnPulseInPulseList(Int_t i){
      cout << "i=" << i << endl;
      cout << "getSizeOfPulseList()=" << getSizeOfPulseList() << endl;

      if(i<=getSizeOfPulseList()){
	list<DmtpcPulse>::iterator iter=_pulseList.begin();
	for(Int_t j=1; j<i; ++j){
	  cout << "going up one on (*iter)" << endl;
	  ++iter;
	}
	cout << "in returnPulseInPulseList(Int_t i=" << i << "): (*iter).getPulseHeightTime()=" << (*iter).getPulseHeightTime() << endl;
	return (*iter);
      }

      return 0;
    }

    /** threshold voltage used in pulse finding */
    float getThreshold() { 
      if( _threshold ){
	return _threshold;
      } else {
	return 20.*_baselineRMS;
      }
    }

    Double_t findPulseStart(TH1F* wfi, Int_t pulseBin , Bool_t debug=false );
    Double_t findPulseEnd(TH1F* wfi, Int_t pulseBin , Bool_t debug=false );
  
  /**
   *  Smooth one histogram by another using convolution
   *
   *  Assumes homogeneous (same width) bins.
   *  
   *  WARNING:  Lightly tested... (James Battat -- 2010 Nov 3)
   *  
   *  The output convolution may be off by 1 bin or so...
   *  Need to check the edge effects here.
   *  
   * @param[in] hdata   input hist to be convolved
   * @param[in] hkernel convolution kernel
   * @return    hconv   The resulting convolution histogram
  */
  static TH1F* smooth(TH1F* hdata, TH1F* hkernel);


  /**
   *  Generate a gaussian smoothing kernel as a TH1F
   *
   * @param[in] binWidth  Width of the bin of the data to be smoothed
   * @param[in] sigma     RMS of gaussian in same units as binWidth
   * @param[in] nsigma    How many sigma to extend the gaussian 
   * @return    hkernel   A gaussian kernel with mean=0 and sigma as specified
   *
   */
  static TH1F* makeGaussKernel(float binWidth, float sigma=1.0, float nsigma=3.0);


  /**
   *  Compute the convolution of two data arrays
   *
   *  This is a completely generic algorithm.  Beware of "edge effects"
   *  when the convolution kernel only partially overlaps with the
   *  data to be smoothed
   *
   *  Note, that given two functions:
   *     f(t) and g(t)
   * 
   *  the convolution (represented by an asterisk) is defined as:
   *     f * g (t) = INTEGRAL(  f(tau) g(t-tau) dtau )
   *
   *  And it is "symmetric":
   *     f * g = g * f
   * 
   *  @param[in]  ff  array containing data to be convolved
   *  @param[in]  nf  length of the ff array
   *  @param[in]  gg  array containing the convolution kernel
   *  @param[in]  ng  length of kernel array
   *  @param[out] hh  resulting convolution (length is nf+ng-1)
   */
  static void convolve(double ff[], int nf, double gg[], int ng, double hh[]);

  /**
   *  YOU PROBABLY WANT TO BE USING smooth() INSTEAD OF THIS FUNCTION
   * 
   *  Compute the convolution between two histograms.
   *
   *  Histograms must have the same number of bins (this could be relaxed)
   *  Assumes homogeneous (same width) bins.
   *  
   *  WARNING:  Lightly tested... (James Battat -- 2010 Oct 1)
   *  
   *  The output convolution may be off by 1 bin or so...
   *  Need to check the edge effects here.
   *  
   *  
   * @param[in] hdata   input hist to be convolved
   * @param[in] hkernel convolution kernel
   * @param[in] flip    Convolution requires one function to be mirrored.  
   *                    To disable the mirroring (which is equiv to doing 
   *                    a cross-correlation), then set this to false. 
   * @param[in] draw    Animate the convolution as it is calculated
   * @return    hconv   The resulting convolution histogram
  */
  static TH1F* convolve(TH1F* hdata, TH1F* hkernel, bool flip=true, bool draw=false);

  static TGraph* graphOfHist(TH1* h1, TString grName="_gr", bool flip=false);


private:
   
   TH1F *_wf; 
   
   float _baseline;
   float _baselineRMS;
   float _peakHeight;
   float _peakHeightAvg;
   float _peakHeightTime;
   float _peakHeightTime10;
   float _peakHeightTime90;
   float _riseTime;
   float _fallTime;
   float _troughDepth;
   
   float _wfMinimum;
   float _wfMinimumTime;

   float _threshold;
   DmtpcPulse* _maxPulse;
   DmtpcPulse* _minPulse;
   list< DmtpcPulse > _pulseList;

   float _timeAbove10;
   float _timeAbove90;
   
   float _peakAveragePosTime;
   float _peakAverageNegTime;
   
   ClassDef(DmtpcWaveformTools,0)
};

#endif
