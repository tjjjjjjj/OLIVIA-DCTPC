#ifndef CLEAN_SKIM_FUNCTIONS_HH
#define CLEAN_SKIM_FUNCTIONS_HH


/* Most of the functionality of cleanSkim is now implemented here, which should make it a lot easier
 * to maintain. 
 */ 



// waveformtools includes
//  for analysis of generic waveforms
#include "WaveformVector.hh"
#include "SkimWaveform.hh"
#include "DmtpcPulse.hh"
//  for analysis of voltage sensitive amplitifier waveforms
#include "FastWfVector.hh"
#include "FastWaveform.hh"
#include "FastPulse.hh"
//  for analysis of charge sensitive pre-amp waveforms
#include "CspWfVector.hh"
#include "CspWaveform.hh"
#include "CspPulse.hh"
//  for analysis of charge sensitive pmt waveforms
#include "PMTWfVector.hh"
#include "PMTWaveform.hh"
#include "PMTPulse.hh"



#include "DmtpcGainMap.hh"
#include "DmtpcSkimDataset.hh"
#include "DmtpcSkimEvent.hh"
#include "cleanSkimConfig.hh"
#include "Dmtpc4ShooterStitcher.hh"



//Class so that can use friend functionality... 
class cleanSkimFunctions
{

  public: 
  

 static bool checkSpark(const char * camid, double mean, double last_mean, double os_spark_mean, const CleanSkimConfig * conf); 

 static int loadGainMaps(DmtpcGainMap ** gainmaps, DmtpcSkimDataset *d, int n, const CleanSkimConfig * conf, const Dmtpc4ShooterStitcher * stitch, string * camera_ids); //done

 static TH2 * cleanImage(const TH2 * raw_img, const TH2 * bias, int & nkilled, double raw_rms, double raw_mean, const TH2 * raw_overscan, double raw_os_mean, double os_bias_mean, const CleanSkimConfig * conf); 


 static TH2 * cleanBiasFrame(const TH2 * in, const CleanSkimConfig * conf); //done

 //TH2 * cleanBiasFrameStack(TTree * in, CleanSkimConfig * conf); 
 
 static Dmtpc4ShooterStitcher * loadStitch(int run, const CleanSkimConfig * conf); 


 static bool checkPartialSpark(const TH2 * img, const CleanSkimConfig * conf); 


 static void populateSparkRef(TH2 * raw_img, TH2 * biases, const CleanSkimConfig * conf, vector<pair<int,int> > * sparkref_running) ; 
                   
 static int fillWaveformVectorsInTObjArray( DmtpcDataset & d, TObjArray * wvflist) ; 

 static int fillCspWfVector(DmtpcDataset & d, CspWfVector *cspwv, int ich, int nch); 
 static int fillPMTWfVector(DmtpcDataset & d, PMTWfVector *pmtwv, int ich, int nch); 
 static int fillFastWfVector(DmtpcDataset & d, FastWfVector *fastwv, int ich, int nch); 
 static bool burnin_test( double x_delta, double y_delta, const CleanSkimConfig * conf); 
 static int parse_args( int nargs, char ** args, TString * files, TString * keys, TString * out, char** cfg); 

 static vector<BurninEncoded_t> computeTrackBurnin(list<vector<double*> **> * positions,
                                                               list<vector<vector<vector<BurninEncoded_t> > >* > * tmpburnin, 
                                                               int position_offset, double x, double y, int cam, int track, const CleanSkimConfig * ); 


 // THIS IS THE PLACE TO ADD MOST NEW VARIABLES //
 static void fillTrackInfo(DmtpcSkimEvent * ev, int u, int v, MaxCamClusterImage * clusti, const CleanSkimConfig * conf, const DmtpcGainMap * gm, const Dmtpc4ShooterStitcher * stitch, const char * cam_id, DmtpcEvent * orig_ev); 



 static int clusterFind(const TH2F * img,const  DmtpcEvent * ev, MaxCamClusterImage ** clusti, const DmtpcGainMap * maps, 
                        const CleanSkimConfig * cfg, const DmtpcStitchInfo * stitch, const double * image_rms,
                        const double * blurred_image_rms,  const double * image_mean, const double * blurred_image_means,  int indx, int cam); 


 static void updateGlobalSparkRef(vector<pair<int,int> > **sparkref_running, list<vector<pair<int,int> > * > * sparkef, int ncamera); 

 static void updateGlobalBurnin( list<vector<double*>**> * positions, list<vector<vector<vector<BurninEncoded_t> > > *> * burnin_temp, 
                                            list<vector<vector<vector<BurninEncoded_t> > > *> * burnin, DmtpcSkimEvent * tmp_event, const CleanSkimConfig * conf, int ncamera, vector<vector<vector<BurninEncoded_t> > > *  this_event_burnin) ; 



 static TH2 * medianBiasFrameStack(TTree * tree, int frame, const CleanSkimConfig * conf, bool os = false); 
};



#endif
