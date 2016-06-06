#ifndef CLEAN_SKIM_CONFIG_HH
#define CLEAN_SKIM_CONFIG_HH

/* cleanSkimConfig.hh 
  
   author: Cosmin Deaconu <cozzyd@mit.edu> 
   
   This file (and the associated source file) are responsible
   for parsing the configuration file for cleanSkim and providing
   cleanSkim with the configuration information. 
*/

#include <string>
#include <strings.h>
#include <ostream>
#include <iostream>
#include "../../../MaxCam/SimpleConfig.hh"
#include "../../../MaxCam/MaxCamClusterImage.hh"


typedef enum 
{
  SQUARE,
  CIRCLE,
  CROSS
}burnin_method_t;


class CleanSkimConfig :  public SimpleConfig
{
  public:
    CleanSkimConfig();
    CleanSkimConfig(char *f);
    
    /* Prints out the current config to the given output 
       stream. */
    void print(std::ostream & out = std::cout); 

    /* Setters and getters for config */
    burnin_method_t getBurninMethod(){ return burnin_method;}
    void setBurninMethod(burnin_method_t method){ burnin_method = method; }
    
    int getBurninNumEvents(){return burnin_num_events;}
    void setBurninNumEvents(int n){burnin_num_events = n;}

    double getBurninDistanceThresh(){ return burnin_distance_thresh;}
    void setBurninDistanceThresh(double t) {burnin_distance_thresh = t;}

    double getOutlierFactor(){ return outlier_factor; }
    void setOutlierFactor(double f){ outlier_factor = f; }

    double getBlurAmount() { return blur_amount; }
    void setBlurAmount(double b) { blur_amount = b; }

    double getClusterMinSigma() { return cluster_min_sigma; }
    void setClusterMinSigma(double s) { cluster_min_sigma = s;}

    double getClusterMaxSigma() { return cluster_max_sigma; }
    void setClusterMaxSigma(double s) { cluster_max_sigma = s;}

    int getClusterMinSize(){ return cluster_min_size; } 
    void setClusterMinSize(int s) { cluster_min_size = s;}

    double getClusterMinDist() { return cluster_min_dist; } 
    void setClusterMinDist(double d) { cluster_min_dist = d;}

    int getRangeAlgorithm() { return range_algorithm; } 
    void setRangeAlgorithm(int a) {range_algorithm = a; }

    int getPhiAlgorithm() { return phi_algorithm; }
    void setPhiAlgorithm(int a) {phi_algorithm = a; }

    double getNorthAngle() { return north_angle; }
    void setNorthAngle(double a) { north_angle = a; } 

    double getLongitude() { return longitude; }
    void setLongitude(double l) {longitude = l;}

    MaxCamClusterImage::CAMERA_ORIENTATION getCameraOrientation(int camera_num); 
    //void setCameraOrientation(int camera_num, MaxCamClusterImage::CAMERA_ORIENTATION orient); 
  
    bool isMC()const {return mc;} 
    void setMC(bool ismc){mc=ismc;} 
    
    bool isNoKill()const {return nokill;}
    void setNoKill(bool isnokill){nokill = isnokill;}

    double getLatitude() const { return latitude; } 
    void setLatitude(double l ) { latitude = l; } 

    char * getChannelId (int i)const  { return channel_ids[i];}
    int getNChannelsPerTrigger()const  { return n_channels_per_trigger;} 

    double getSatThresh() const { return satThresh; } 
    int getNSatThresh() const { return nSatThresh; } 

    //these will all return the highest value if i exceeds the array bounds 
    int getWaveformNBaselineLength() const { return nbaseline_l;}
    int getWaveformNBaseline(int i = 0) const {  return (i >= nbaseline_l ? nbaseline[nbaseline_l-1] : nbaseline[i]); } 
    int getWaveformNBlurLength() const { return nblur_l;}
    int getWaveformNBlur(int i = 0) const { return (i >= nblur_l ? nblur[nblur_l-1] : nblur[i]); } 
    int getInvertWaveformLength() const { return invert_waveform_l;}
    bool getInvertWaveform(int i = 0) const { return (bool) ( i >= invert_waveform_l ? invert_waveform[invert_waveform_l-1] : invert_waveform[i]); } 
    

  private:
    void parsePair(std::string key, std::string value, int n);
    int burnin_num_events;
    double burnin_distance_thresh;
    double outlier_factor;
    double blur_amount;
    double cluster_min_sigma;
    double cluster_max_sigma;
    int cluster_min_size;
    double cluster_min_dist;
  
    int range_algorithm;
    int phi_algorithm;
     
    double north_angle;
    double longitude;
    double latitude;
    burnin_method_t burnin_method; 
  
    int * top_cameras;
    int n_top_cameras; 
    int * bottom_cameras; 
    int n_bottom_cameras; 
    
    int nbaseline_l; 
    int nblur_l; 
    int invert_waveform_l; 

    int * nbaseline; 
    int * nblur; 
    int * invert_waveform; 

    void setDefaults();
    void clearMemory();  

    bool mc; 
    bool nokill; 

    int n_channels_per_trigger; 
    char ** channel_ids; 

    double satThresh;
    int nSatThresh; 
};

#endif

