// ScopeWaveform.hh
//
#ifndef SCOPE_WAVEFORM_DATA_HH
#define SCOPE_WAVEFORM_DATA_HH

#include "TH1.h"

class TString;

/** 
 Class to contain waveforms from an oscilloscope for data storage. Data is 
 stored in a TH1C in digital units, which stores signed chars, but the data is unsigned, 
 so plotting it directly will not produce expected results. 

 The scopeExpand in DmtpcDataConverter may be used to convert one of these
 waveforms into a TH1F in physical units for plotting and analysis. 

 \author Asher Kaboth (code) / Cosmin Deaconu (documentation) 
 
*/ 
class ScopeWaveformData : public TH1D { 

public:

  // Destructor
  ~ScopeWaveformData();

  // Constructors
  ScopeWaveformData();
  //ScopeWaveformData(const ScopeWaveformData &other);  // copy constructor
  
  /** 
   * Create an empty ScopeWaveformData 
   *
   * @param histoName the name to give to the waveform
   * @param nSamples The number of samples (bins) 
   * @param timeStamp the timestamp of the trace, understood to be seconds
   *                  since the beginning of trigger taking
   * @param leveltoVolt the conversion between one digital unit and one analog unit in volts
   * @param zero the zero level of the data (generally 128) 
   **/
  ScopeWaveformData(TString histoName, int nSamples, 
		    double timeStamp, float leveltoVolt, short zero);
   /** 
   * Create an empty ScopeWaveformData with time units on the x axis
   *
   * @param histoName the name to give to the waveform
   * @param nSamples The number of samples (bins) 
   * @param xlow The minimum time value
   * @param xup The maximum time value
   * @param timeStamp the timestamp of the trace, understood to be seconds
   *                  since the beginning of trigger taking
   * @param leveltoVolt the conversion between one digital unit and one analog unit in volts
   * @param zero the zero level of the data (generally 128) 
   **/
  ScopeWaveformData(TString histoName, int nSamples, float xlow, float xup,
		    double timeStamp, float leveltoVolt, short zero);

  /** 
   * Returns the timestamp of the trace, which is the number of
   * seconds since the beginning of trigger taking 
   *
   * @return timestamp
   */
  double getTimeStamp() const {return _timeStamp;}

  /** 
   * Returns the conversion factor between digital units and analog units (V)
   * 
   * @return unit conversion 
   */
  float getleveltoVolt() const {return _leveltoVolt;}

  /** Returns the baseline level in digital units (generally 128)
   *
   * @return zero level 
   */

  short getZero() const{return _zero;}

  int zeroSuppress(double nsigma = 3, unsigned nrequired = 2, unsigned window_size = 10, int nsamples_baseline=1000, double max_deviation_sigma = 1); 
  void cropZeros();
  void resample(double factor = 2);

  TH1S * getAsShortHist() const; 
  TH1D * getAsSignedHist() const; 
  void convertToSignedHist(); 
  bool isSigned() const { return _signed; } 

private:  
  double _timeStamp;    
  double _nresample; 
  int _nSamples; 
  unsigned _nsuppressed; 
  float _leveltoVolt;
  float _suppressed_rms; 
  short _zero;    
  bool _signed; 

 
  ClassDef(ScopeWaveformData,3)
};

#endif
