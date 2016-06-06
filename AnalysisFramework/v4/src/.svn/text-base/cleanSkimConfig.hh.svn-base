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
#include <map>
#include <vector>
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
    CleanSkimConfig(const char *f);
    
    /* Prints out the current config to the given output 
       stream. */
    void print(std::ostream & out = std::cout); 

    /* Setters and getters for config */
    burnin_method_t getBurninMethod(){ return burnin_method;}
    int getBurninNumEvents(){return burnin_num_events;}
    double getBurninDistanceThresh(){ return burnin_distance_thresh;}
    double getOutlierFactor(){ return outlier_factor; }
    double getBlurAmount() { return blur_amount; }
    double getGaussianBlurAmount() { return gaussian_blur_amount; }
    double getClusterMinSigma() { return cluster_min_sigma; }

    double getClusterMaxSigma() { return cluster_max_sigma; }
    int getClusterMinSize(){ return cluster_min_size; } 
    int getClusterMinSizeUnbinned(){ return cluster_min_size_unbinned; } 
    double getClusterMinDist() { return cluster_min_dist; } 
    int getRangeAlgorithm() { return range_algorithm; } 
    int getPhiAlgorithm() { return phi_algorithm; }
    double getNorthAngle() { return north_angle; }
    double getLongitude() { return longitude; }
    MaxCamClusterImage::CAMERA_ORIENTATION getCameraOrientation(const char * camera_id); 
    bool isMC()const {return mc;} 
    bool normalizeGainMap() const {return normalize_gainmap; } 
    bool Bias() const {return bias; }
    bool isNoKill()const {return nokill;}
    double getLatitude() const { return latitude; } 
    char * getChannelId (int i)const  { return channel_ids[i];}
    char * getChannelType (int i)const  { return channel_types[i];}
    int getNChannelsPerTrigger()const  { return n_channels_per_trigger;} 
    double getSatThresh() const { return satThresh; } 
    int getNSatThresh() const { return nSatThresh; } 
    //these will all return the highest value if i exceeds the array bounds 
    bool useGainMap() const { return usegainMap; } 
    const char * getGainMapFile() const { return gainMapFile; } 
    bool hasOverscan() const{ return overscan; } 
    void setHasOverscan(bool has) { overscan = has;} 
    void setNScopeChannels(int n) { n_channels_per_trigger = n;}
    int getClusterBlurRadius() const{ return cluster_blur_radius; } 
    double getClusterThreshRingRatio() const{ return cluster_thresh_ring_ratio; } 
    double getClusterMaxWrongProb() const{ return cluster_max_wrong_prob; } 
    double getClusterMinThreshold() const{ return cluster_min_threshold; } 
    double getClusterSeedThreshold() const{ return cluster_seed_threshold; } 
    int getClusterNeighborsThresholdForFilling()  const{ return cluster_neighbors_threshold_for_filling; } 
    int getClusterMinNeighborsToKeepPixel() const { return cluster_min_neighbors_to_keep_pixel; } 
    double getClusterJoinMinRxyCluster() const { return cluster_join_min_rxy_cluster;} 
    double getClusterJoinMinRxyGlobal() const { return cluster_join_min_rxy_global;} 
    double getClusterJoinMaxJoinResidual() const { return cluster_join_max_join_residual;} 
    double getClusterJoinLeastSquaresWeight() const { return cluster_join_least_squares_weight;} 
    double getClusterJoinSpacerWidth() const { return cluster_join_spacer_width;}
    char * getClusterFindingAlgorithm(int i = 0) const { return cf_algos[i]; } 
    unsigned getNClusterFindingAlgorithms() const { return cf_algos.size(); } 
    double getOverscanSparkCut(char * cam_id) { return overscan_spark_cut.count(std::string(cam_id)) ? overscan_spark_cut[std::string(cam_id)] : DBL_MAX; } 
    double getSparkCut(char * cam_id) { return spark_cut.count(std::string(cam_id)) ? spark_cut[std::string(cam_id)] : 1.01; } 
    double getClusterFindADGradientBlurAmount() const { return cf_ad_gradient_blur; }
    double getClusterFindADZeroNegatives() const { return cf_ad_zero; }
    double getClusterFindADHighThresh() const { return cf_ad_high_thresh; }
    double getClusterFindADLowThresh() const { return cf_ad_low_thresh; }
    double getClusterFindADK() const { return cf_ad_k; }
    double getClusterFindADLambda() const { return cf_ad_lambda; }
    double getClusterFindCIUseGainMapToMerge() const { return cluster_ci_use_gain_map_to_merge; } 
    double getClusterReducedThreshold() const { return cluster_red_threshold; } 
    const char * getFallBackCameraID(unsigned int i = 0) const { return fallback_camera_ids.size() > i ? fallback_camera_ids[i].c_str() : ""; }

    int getClusterFindADNIter() const { return cf_ad_ntries;} 
    double getPartialBlur() const {return partial_blur;}
    double getPartialLowThresh() const {return partial_low_thresh;}
    double getPartialHighThresh() const {return partial_high_thresh;} 
    double getPartialNPrimaryThresh() const {return partial_n_primary_thresh;}
    double getPartialNSecondaryThresh() const {return partial_n_secondary_thresh;}
    double getPartialDistanceLow() const {return partial_dist_low;}
    double getPartialDistanceHigh() const {return partial_dist_high;}


  private:
    void parsePair(std::string key, std::string value, int n);
    int burnin_num_events;
    double burnin_distance_thresh;
    double cluster_red_threshold; 
    double outlier_factor;
    double blur_amount;
    double cluster_min_sigma;
    double cluster_max_sigma;
    int cluster_min_size;
    int cluster_min_size_unbinned;
    double cluster_min_dist;

    double gaussian_blur_amount; 
    char * cluster_algo; 

    int cf_ad_ntries; 
    bool cf_ad_zero; 
    double cf_ad_gradient_blur; 
    double cf_ad_high_thresh; 
    double cf_ad_low_thresh; 
    double cf_ad_k; 
    double cf_ad_lambda; 

    int cluster_blur_radius; 
    double cluster_thresh_ring_ratio; 
    double cluster_max_wrong_prob; 
    double cluster_min_threshold; 
    double cluster_seed_threshold;; 
    int cluster_neighbors_threshold_for_filling; 
    int cluster_min_neighbors_to_keep_pixel; 

    double cluster_join_min_rxy_cluster;  
    double cluster_join_min_rxy_global; 
    double cluster_join_max_join_residual; 
    double cluster_join_least_squares_weight; 
    double cluster_join_spacer_width; 

    double partial_blur;
    double partial_low_thresh;
    double partial_high_thresh;
    double partial_n_primary_thresh;
    double partial_n_secondary_thresh;
    double partial_dist_low;
    double partial_dist_high;

    int range_algorithm;
    int phi_algorithm;
     
    double north_angle;
    double longitude;
    double latitude;
    burnin_method_t burnin_method; 
  
    char ** top_cameras;
    int n_top_cameras; 
    char ** bottom_cameras; 
    int n_bottom_cameras; 

    bool overscan; 
    bool charge;

    void setDefaults();
    void clearMemory();  

    bool mc; 
    bool nokill; 
    bool normalize_gainmap; 
    bool bias;
    int n_channels_per_trigger; 
    char ** channel_ids; 
    vector<char * > cf_algos; 
    char ** channel_types; 

    bool usegainMap; 
    char * gainMapFile; 
    
    bool cluster_ci_use_gain_map_to_merge; 

    double satThresh;
    int nSatThresh; 

    std::map<std::string, double> overscan_spark_cut; 
    std::map<std::string, double> spark_cut; 
    std::vector<std::string> fallback_camera_ids; 
};

#endif

