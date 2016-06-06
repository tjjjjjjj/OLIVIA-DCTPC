#include "cleanSkimConfig.hh" 
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>

CleanSkimConfig::CleanSkimConfig()
{
  setDefaults(); 
}

CleanSkimConfig::CleanSkimConfig(const char * file)
{
  setDefaults();
  parseFile(file); 
}

void CleanSkimConfig::setDefaults()
{
 burnin_method = SQUARE; 


 compression_algo = ROOT::kZLIB; 
 char ** compression_keys = new char * [5]; 
 compression_keys[0]="global"; 
 compression_keys[1]="zlib"; 
 compression_keys[2]="lzma"; 
 compression_keys[3]="old"; 
 compression_keys[4]="future"; 

 compression_level = 5;
 registerNewInt("compression_level",&compression_level); 
 registerNewEnum("compression_algo",5,compression_keys,0,&compression_algo); 


 overscan = true; 
 registerNewBool("overscan", &overscan); 

 normalize_gainmap = true; 
 registerNewBool("normalize_gainmap",&normalize_gainmap); 

 interpolation_method = strdup("bicubic");
 registerNewString("interpolation_method",&interpolation_method); 

 burnin_num_events = 1000;
 registerNewInt("burnin_num_events",&burnin_num_events); 

 burnin_distance_thresh = 12;
 registerNewDouble("burnin_distance_thresh",&burnin_distance_thresh); 

 outlier_factor = 1.25;
 registerNewDouble("outlier_factor",&outlier_factor); 

 blur_amount = 0.8;
 registerNewDouble("blur_amount",&blur_amount); 

 gaussian_blur_amount = 2; 
 registerNewDouble("gaussian_blur_amount",&gaussian_blur_amount); 
 
 cluster_min_sigma = 3.7;
 registerNewDouble("cluster_min_sigma",&cluster_min_sigma); 

 cluster_max_sigma = 300;
 registerNewDouble("cluster_max_sigma", &cluster_max_sigma); 

 cluster_min_size = 5;
 registerNewInt("cluster_min_size", &cluster_min_size); 

 cluster_min_size_unbinned = 16; 
 registerNewInt("cluster_min_size_unbinned", &cluster_min_size_unbinned); 

 cluster_min_dist = 3600;
 registerNewDouble("cluster_min_dist", &cluster_min_dist); 

 range_algorithm = 1;
 registerNewInt("range_algorithm", &range_algorithm); 

 phi_algorithm = 2;
 registerNewInt("phi_algorithm", &phi_algorithm); 

 north_angle = 240;
 registerNewDouble("north_angle", &north_angle); 

 longitude = -71.110556;
 registerNewDouble("longitude",&longitude); 

 latitude = 42.373611;
 registerNewDouble("latitude",&latitude); 

 partial_blur = 2;
 registerNewDouble("partial_blur",&partial_blur);

partial_low_thresh = 5.5;
 registerNewDouble("partial_low_thresh",&partial_low_thresh);
 
 partial_high_thresh = 18;
 registerNewDouble("partial_high_thresh",&partial_high_thresh);

 partial_n_primary_thresh = 85;
 registerNewDouble("partial_n_primary_thresh",&partial_n_primary_thresh);

 partial_n_secondary_thresh = 70;
 registerNewDouble("partial_n_secondary_thresh",&partial_n_secondary_thresh);

 partial_dist_low = 2;
 registerNewDouble("partial_dist_low",&partial_dist_low);

 partial_dist_high = 8;
 registerNewDouble("partial_dist_high",&partial_dist_high);


 mc = false; 
 registerNewBool("mc",&mc); 

 nokill=false; 
 registerNewBool("nokill",&nokill); 

 cluster_red_threshold = 3.7; 
 registerNewDouble("cluster_red_threshold",&cluster_red_threshold); 

 _stitch = false; 
 registerNewBool("stitch", &_stitch); 

 //Initialize camera orientations
 n_top_cameras = 1; 
 n_bottom_cameras = 1; 
 top_cameras = new char * [1];  
 bottom_cameras = new char * [1]; 
 top_cameras[0] = strdup("081264"); 
 bottom_cameras[0] = strdup("100439"); 

 //Initialize scope information 
 n_channels_per_trigger = 4; 
 channel_ids = new char*[4]; 
 channel_ids[0] = "top";
 channel_ids[1] = "bottom"; 
 channel_ids[2] = "mesh_bottom";
 channel_ids[3] = "mesh_top"; 
 channel_types = new char*[4]; 
 channel_types[0] = "csp";
 channel_types[1] = "csp"; 
 channel_types[2] = "fast";
 channel_types[3] = "fast"; 

 // Initialize sparkref info
 satThresh = 55000; 
 registerNewDouble("sat_thresh", &satThresh); 

 nSatThresh = 4; 
 registerNewInt("n_sat_thresh", &nSatThresh); 

  usegainMap = false; 
  gainMapFile = NULL;
  stitch_dir = NULL; 
  stitch_file = strdup("auto"); 

  cluster_blur_radius = 4; 
  registerNewInt("cluster_find_seed_blur_radius",&cluster_blur_radius); 

  cluster_thresh_ring_ratio = 0.75;
  registerNewDouble("cluster_find_seed_thresh_ring_ratio", &cluster_thresh_ring_ratio); 

  cluster_max_wrong_prob = 1e-33; 
  registerNewDouble("cluster_find_seed_max_wrong_prob", &cluster_max_wrong_prob); 

  cluster_min_threshold = 0.75; 
  registerNewDouble("cluster_find_seed_min_threshold", &cluster_min_threshold); 

  cluster_seed_threshold = 5.2;
  registerNewDouble("cluster_find_seed_seed_threshold",&cluster_seed_threshold); 

  cluster_neighbors_threshold_for_filling = 4; 
  registerNewInt("cluster_neighbors_threshold_for_filling", &cluster_neighbors_threshold_for_filling);

  cluster_min_neighbors_to_keep_pixel = 2; 
  registerNewInt("cluster_min_neighbors_to_keep_pixel", &cluster_min_neighbors_to_keep_pixel);

  cluster_join_min_rxy_cluster = 0.8; 
  registerNewDouble("cluster_join_min_rxy_cluster",&cluster_join_min_rxy_cluster); 

  cluster_join_min_rxy_global = 0.65; 
  registerNewDouble("cluster_join_min_rxy_global",&cluster_join_min_rxy_global); 

  cluster_join_max_join_residual = 5; 
  registerNewDouble("cluster_join_max_join_residual", &cluster_join_max_join_residual); 

  cluster_join_least_squares_weight = 2; 
  registerNewDouble("cluster_join_least_squares_weight", &cluster_join_least_squares_weight); 
     
  cf_ad_ntries = 30; 
  registerNewInt("cluster_find_ad_niter",&cf_ad_ntries); 

  cf_ad_zero = true; 
  registerNewBool("cluster_find_ad_zero_neg", &cf_ad_zero); 

  cf_ad_gradient_blur = 2.5; 
  registerNewDouble("cluster_find_ad_gradient_blur",&cf_ad_gradient_blur); 

  cf_ad_high_thresh = 8.5; 
  registerNewDouble("cluster_find_ad_high_thresh",&cf_ad_high_thresh); 

  cf_ad_low_thresh = 4; 
  registerNewDouble("cluster_find_ad_low_thresh",&cf_ad_low_thresh); 

  cf_algos.push_back(strdup("ci")); 

  cluster_join_spacer_width = 2; 
  registerNewDouble("cluster_join_spacer_width",&cluster_join_spacer_width); 

  cf_ad_k = 25; 
  registerNewDouble("cluster_find_ad_k",&cf_ad_k); 

  cf_ad_lambda = 1; 
  registerNewDouble("cluster_find_ad_lambda",&cf_ad_lambda); 

  cf_ring_space_sigma = 4; 
  registerNewDouble("cluster_find_ring_space_sigma",&cf_ring_space_sigma); 
  cf_ring_rms_sigma = 1; 
  registerNewDouble("cluster_find_ring_rms_sigma",&cf_ring_rms_sigma); 
  cf_ring_core_thresh_low = 0.75;
  registerNewDouble("cluster_find_ring_core_thresh_low",&cf_ring_core_thresh_low); 
  cf_ring_core_thresh_high = 1.4;
  registerNewDouble("cluster_find_ring_core_thresh_high",&cf_ring_core_thresh_high); 
  cf_ring_ring_thresh = 0.5;
  registerNewDouble("cluster_find_ring_ring_thresh",&cf_ring_ring_thresh); 
  cf_ring_nsigma = 3; 
  registerNewDouble("cluster_find_ring_ring_nsigma",&cf_ring_nsigma); 
  cf_ring_ncleanup = 0; 
  registerNewInt("cluster_find_ring_ncleanup",&cf_ring_ncleanup); 


  overscan_spark_cut["081264"] = -1; 
  overscan_spark_cut["100439"] = -7.5; 

  cluster_ci_use_gain_map_to_merge = false; 
  registerNewBool("cluster_find_ci_use_gain_map_for_merging", &cluster_ci_use_gain_map_to_merge); 

  spacer_join_r = 30; 
  spacer_join_theta = 0.2; 
  registerNewDouble("spacer_join_r",&spacer_join_r); 
  registerNewDouble("spacer_join_theta",&spacer_join_theta); 


  char ** bias_clean_keys = new char * [2]; 
  bias_clean_keys[0] = "traditional"; 
  bias_clean_keys[1] = "median_stack"; 
  registerNewEnum("bias_clean_method",2, bias_clean_keys,0, &bias_clean_method); 

  char ** image_clean_keys = new char * [3]; 
  image_clean_keys[0] = "traditional"; 
  image_clean_keys[1] = "bias_subtract_first"; 
  image_clean_keys[2] = "median_difference"; 
  registerNewEnum("image_clean_method",3, image_clean_keys,0, &image_clean_method); 

  char ** image_ratio_keys = new char * [2]; 
  image_ratio_keys[0] = "mean"; 
  image_ratio_keys[1] = "median"; 
  registerNewEnum("image_ratio_method",2, image_ratio_keys,0, &image_ratio_method); 

  char ** save_type_keys = new char * [3]; 
  save_type_keys[0] = "float"; 
  save_type_keys[1] = "short"; 
  save_type_keys[2] = "int"; 
  registerNewEnum("hist_save_type",3,save_type_keys,0, &hist_save_type); 


  round_amt = 0; 
  registerNewDouble("round_amount",&round_amt); 

  round_in_clusters = false; 
  registerNewBool("round_in_clusters",&round_in_clusters); 

  seed_cluster_find_reproduce_v4_bug = true; 
  registerNewBool("seed_cluster_find_reproduce_v4_bug",&seed_cluster_find_reproduce_v4_bug); 

}

void parse_int_list(const char * str, int * length, int ** values, bool del = true)
{
  if (del) delete (*values); 
  //count number of tokens  
  int ntokens = 1; 
  for (unsigned int i = 0; i < strlen(str); i++)
   {
       if (str[i] == ',') ntokens++;  
   }
   *length = ntokens;   
   *values = new int[ntokens]; 
   int n = 0; 
   char * ch;  
   char * strcp= strdup(str); //Need non-const string
   ch = strtok(strcp,",");
   while (ch!=NULL)
   {
     (*values)[n++] = atoi(ch); 
     ch = strtok(NULL,",");   
   }
   free(strcp); 
   return;
}

void parse_channel_ids (const char * str, char *** channel_ids, int * arr_l)
{
  
  //Check for none
  if (strcasecmp(str,"none")==0)
  {
    *arr_l = 0; 
    *channel_ids = NULL; 
    return; 
  }
  
  int ntokens = 1; 
  //Figure out how many tokens there are
  for (unsigned int i = 0; i < strlen(str); i++)
  {
    if (str[i] == ',') ntokens++;  
  }

  *arr_l = ntokens; 
  *channel_ids = (char **) malloc(sizeof(char*) * ntokens); 

  if (ntokens==1)
  {
    (*channel_ids)[0] = strdup(str); 
    return; 
  }
  
  int n = 0; 
  char * ch;  
  char * strcp= strdup(str); //Need non-const string
  ch = strtok(strcp,",");
  while (ch!=NULL)
  {
    std::cout << ch << std::endl; 
    (*channel_ids)[n++] = strdup(ch);
    ch = strtok(NULL,",");   
  }
  free(strcp); 
}
  

void parse_camera_orientations(const char * str, char *** camera_arr, int * arr_l)
{
  //Check for none
  if (strcasecmp(str,"none")==0)
  {
    *arr_l = 0; 
    *camera_arr = NULL; 
    return; 
  }
  
  int ntokens = 1; 
  //Figure out how many tokens there are
  for (unsigned int i = 0; i < strlen(str); i++)
  {
    if (str[i] == ',') ntokens++;  
  }

  *arr_l = ntokens; 
  *camera_arr = (char **) malloc(sizeof(char *) * ntokens); 

  if (ntokens==1)
  {
    (*camera_arr)[0] = strdup(str); 
    return; 
  }
  
  int n = 0; 
  char * ch;  
  char * strcp= strdup(str); //Need non-const string
  ch = strtok(strcp,",");
  while (ch!=NULL)
  {
    (*camera_arr)[n++] = strdup(ch);
    ch = strtok(NULL,",");   
  }
  free(strcp); 
}

void CleanSkimConfig::parsePair(std::string key, std::string value, int n) 
{
  if (true) 
  {
    std::cout << "key: \"" << key << "\"\tvalue: \"" << value << "\"" << std::endl; 
  }

  if (key == "cluster_find_algo")
  {
    for (unsigned i = 0; i < cf_algos.size(); i++)
    {
      free(cf_algos[i]); 
    }
    cf_algos.clear(); 

    char * ch; 
    char * strcp = strdup(value.c_str()); 
    ch = strtok(strcp,","); 
    while (ch)
    {
      cf_algos.push_back(strdup(ch)); 
      std::cout << "Adding algorithm: " << ch << std::endl; 
      ch = strtok(0,","); 
    }
    free (strcp); 
    return;
  }

  if (key == "overscan_spark_cut") 
  {
    overscan_spark_cut.clear(); 

    if (!strcasecmp(value.c_str(),"none") || !strcasecmp(value.c_str(),""))
    {
       return;   
    }

    int ntokens = 1; 

    for (unsigned int i = 0; i < value.length() ; i++)
    {
      if (value[i] == ',') ntokens++; 
    }

    char buf[80]; 
    double f; 
    if (ntokens == 1)
    {
      sscanf(value.c_str(),"%[^:]:%lf",buf,&f);  
      overscan_spark_cut[buf] = f; 
      return;
    }
   
   char * ch; 
   char * strcp = strdup(value.c_str()); 
   ch = strtok(strcp,","); 
   while (ch!=NULL)
   {
     sscanf(ch,"%[^:]:%lf",buf,&f);  
     overscan_spark_cut[buf] = f; 
     std::cout << buf << ": " << f << std::endl; 
     ch = strtok(NULL,","); 
   }

   free (strcp); 
   return;
  }

  if (key == "gain_map_file")
  {
    usegainMap = true; 
    if (gainMapFile) free(gainMapFile); 
    gainMapFile = strdup(value.c_str()); 
    return; 

  }

  if (key == "stitch_file")
  {
    if (stitch_file) free (stitch_file) ; 
    stitch_file = strdup(value.c_str()); 
    return;
  }
  if (key == "stitch_dir")
  {
    if (stitch_dir) free (stitch_dir) ; 
    stitch_dir = strdup(value.c_str()); 
    return;
  }

  if (key == "spark_cut") 
  {
    spark_cut.clear(); 

    if (!strcasecmp(value.c_str(),"none") || !strcasecmp(value.c_str(),""))
    {
       return;   
    }

    int ntokens = 1; 

    for (unsigned int i = 0; i < value.length() ; i++)
    {
      if (value[i] == ',') ntokens++; 
    }

    char buf[80]; 
    double f; 
    if (ntokens == 1)
    {
      sscanf(value.c_str(),"%[^:]:%lf",buf,&f);  
      spark_cut[buf] = f; 
      return;
    }
   
   char * ch; 
   char * strcp = strdup(value.c_str()); 
   ch = strtok(strcp,","); 
   while (ch!=NULL)
   {
     sscanf(ch,"%[^:]:%lf",buf,&f);  
     spark_cut[buf] = f; 
     std::cout << buf << ": " << f << std::endl; 
     ch = strtok(NULL,","); 
   }

   free (strcp); 
   return;
  }

  if (key == "gain_map_file")
  {    usegainMap = true; 
    if (gainMapFile) free(gainMapFile); 
    gainMapFile = strdup(value.c_str()); 
    return; 

  }

  if (key == "burnin_method")
  {
    //std::string has no case insensitive compare... grumble grumble
    if (strcasecmp(value.c_str(),"square")==0)
    {
      burnin_method = SQUARE;
      return;
    }
    if (strcasecmp(value.c_str(),"circle")==0)
    {
      burnin_method = CIRCLE;
      return;
    }
    if (strcasecmp(value.c_str(),"cross")==0)
    {
      burnin_method = CROSS;
      return;
    }

    burnin_method = SQUARE; 
    return; 
  }

  if (key == "cameras_orientation_top")
  {
    for (int i = 0; i < n_top_cameras; i++) free(top_cameras[i]); 
    delete top_cameras;
    parse_camera_orientations(value.c_str(), &top_cameras, &n_top_cameras);
    return; 
  }
  
  if (key == "cameras_orientation_bottom")
  {
    for (int i = 0; i < n_bottom_cameras; i++) free(bottom_cameras[i]); 
    delete bottom_cameras; 
    parse_camera_orientations(value.c_str(), &bottom_cameras, &n_bottom_cameras);
    return; 
  }


  if (key == "fallback_camera_ids")
  {
    fallback_camera_ids.clear(); 

    char * ch;  
    char * strcp= strdup(value.c_str()); //Need non-const string
    ch = strtok(strcp,",");
    while (ch!=NULL)
    {
      fallback_camera_ids.push_back(std::string(ch)); 
      ch = strtok(NULL,",");   
    }
    free(strcp); 
    return;
  }

  if (key == "channel_ids")
  {
    delete channel_ids; 
    std::cout << "parsing channel_ids" << std::endl; 
    parse_channel_ids(value.c_str(), &channel_ids, &n_channels_per_trigger); 
    return; 
  }

  if (key == "channel_types")
  {
    delete channel_types; 
    std::cout << "parsing channel_types" << std::endl; 
    parse_channel_ids(value.c_str(), &channel_types, &n_channels_per_trigger); 
    return; 
  }

  if (parseRegistered(key,value,n)) 
  {
    return; 
  }
      
//If we get here, that means there was an unrecognized key
  std::cerr << "Unrecognized key \""<< key <<"\" on line " 
            << n << ". Skipping." << std::endl;

}


MaxCamClusterImage::CAMERA_ORIENTATION 
CleanSkimConfig::getCameraOrientation(const char * cam) const
{
  for (int i = 0; i < n_top_cameras; i++)
  {
    if (!strcmp(top_cameras[i],cam))
      return MaxCamClusterImage::TOP; 
  }

  for (int i = 0; i < n_bottom_cameras; i++)
  {
    if (!strcmp(bottom_cameras[i],cam))
       return MaxCamClusterImage::BOTTOM; 
  }

  //If we got this far, that means the camera wasn't specified. 
  //Warn and return an up camera

  std::cerr << "WARNING: No orientation specified for camera " << cam << "." << std::endl;
  std::cerr << "  Assuming orientation is TOP" << std::endl; 
  return MaxCamClusterImage::TOP; 
}


void CleanSkimConfig::print(std::ostream & out) const
{

  printRegistered(out); 

  out << "burnin_method = ";
  switch (burnin_method)
  {
    case CROSS:
      out << "CROSS"  << std::endl; break;
    case SQUARE:
      out << "SQUARE" << std::endl; break;
    case CIRCLE:
      out << "CIRCLE" << std::endl; break;
    default: out << "INVALID" << std::endl; 
  }

  out << "cameras_orientation_top = "; 
  if (n_top_cameras == 0) out << "none" << std::endl; 
  else
  {
    out << top_cameras[0]; 
    for (int i = 1; i < n_top_cameras; i++)
    {
      out << "," << top_cameras[i]; 
    } out << std::endl; }

  out << "cameras_orientation_bottom = "; 
  if (n_bottom_cameras == 0) out << "none" << std::endl; 
  else
  {
    out << bottom_cameras[0]; 
    for (int i = 1; i < n_bottom_cameras; i++)
    {
      out << "," << bottom_cameras[i]; 
    }
    out << std::endl;
  }

  out << "channel_ids = "; 
  if (n_channels_per_trigger == 0) out << "none" << std::endl; 
  else
  {
    out << channel_ids[0]; 
    for (int i = 1; i < n_channels_per_trigger; i++)
    {
      out << "," << channel_ids[i]; 
    }
    out << std::endl; 
  }

  out << "channel_types = "; 
  if (n_channels_per_trigger == 0) out << "none" << std::endl; 
  else
  {
    out << channel_types[0]; 
    for (int i = 1; i < n_channels_per_trigger; i++)
    {
      out << "," << channel_types[i]; 
    }
    out << std::endl; 
  }

  if (gainMapFile!=NULL)
  out << "gain_map_file = " << gainMapFile << std::endl; 
 
  out << "overscan_spark_cut = "; 

  int i = 0; 
  for (std::map<std::string,double>::const_iterator it = overscan_spark_cut.begin(); it!= overscan_spark_cut.end(); it++)
  {
    if (i++>0) out << " , " << std::endl; 
    out << it->first << ":" << it->second; 
  }
  out << std::endl; 

  out << "spark_cut = "; 

  i = 0; 
  for (std::map<std::string,double>::const_iterator it = spark_cut.begin(); it!= spark_cut.end(); it++)
  {
    if (i++>0) out << " , " << std::endl; 
    out << it->first << ":" << it->second; 
  }
  out << std::endl; 
  
}
