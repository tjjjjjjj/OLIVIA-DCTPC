#include "cleanSkimConfig.hh" 
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>

CleanSkimConfig::CleanSkimConfig()
{
  setDefaults(); 
}

CleanSkimConfig::CleanSkimConfig(char * file)
{
  setDefaults();
  parseFile(file); 
}

void CleanSkimConfig::setDefaults()
{
 burnin_method = SQUARE; 
 burnin_num_events = 1000;
 burnin_distance_thresh = 12;
 outlier_factor = 1.25;
 blur_amount = 0.8;
 cluster_min_sigma = 3.7;
 cluster_max_sigma = 300;
 cluster_min_size = 5;
 cluster_min_dist = 3600;
 range_algorithm = 1;
 phi_algorithm = 2;
 north_angle = 240;
 longitude = -71.110556;
 latitude = 42.373611;
 mc = false; 
 nokill=false; 

 //Initialize camera orientations
 n_top_cameras = 1; 
 n_bottom_cameras = 1; 
 top_cameras = new int[1];  
 bottom_cameras = new int[1]; 
 top_cameras[0] = 0; 
 bottom_cameras[0] = 1; 

 //Initialize scope information 

 n_channels_per_trigger = 2; 
 channel_ids = new char*[2]; 
 channel_ids[0] = "top";
 channel_ids[1] = "bottom"; 

 // Initialize sparkref info
 satThresh = 55000; 
 nSatThresh = 4; 

//Waveform tools defaults

nblur_l = 1; 
nbaseline_l = 1; 
invert_waveform_l = 1; 

nblur = new int[1]; 
nblur[0] = 1; 

nbaseline = new int[1];
nbaseline[0] = 25; 

invert_waveform = new int[1];
invert_waveform[0] = 0; 

}

void parse_int_list(const char * str, int * length, int ** values, bool del = true)
{
  if (del) delete (*values); 
  //count number of tokens  
  int ntokens = 1; 
  for (int i = 0; i < strlen(str); i++)
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
  for (int i = 0; i < strlen(str); i++)
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
  

void parse_camera_orientations(const char * str, int ** camera_arr, int * arr_l)
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
  for (int i = 0; i < strlen(str); i++)
  {
    if (str[i] == ',') ntokens++;  
  }

  *arr_l = ntokens; 
  *camera_arr = (int *) malloc(sizeof(int) * ntokens); 

  if (ntokens==1)
  {
    (*camera_arr)[0] = atoi(str); 
    return; 
  }
  
  int n = 0; 
  char * ch;  
  char * strcp= strdup(str); //Need non-const string
  ch = strtok(strcp,",");
  while (ch!=NULL)
  {
    (*camera_arr)[n++] = atoi(ch);
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

  if (key == "burnin_num_events")
  {
    burnin_num_events = atoi(value.c_str()); 
    return; 
  }

  if (key == "burnin_distance_thresh")
  {
    burnin_distance_thresh = atof(value.c_str()); 
    return;
  }

  if (key == "outlier_factor")
  {
    outlier_factor = atof(value.c_str()); 
    return;
  }

  if (key == "blur_amount")
  {
    blur_amount = atof(value.c_str()); 
    return;
  }
  
  if (key == "cluster_min_sigma")
  {
    cluster_min_sigma = atof(value.c_str()); 
    return;
  }

  if (key == "cluster_max_sigma")
  {
    cluster_max_sigma = atof(value.c_str()); 
    return;
  }

  if (key == "cluster_min_size")
  {
    cluster_min_size = atoi(value.c_str()); 
    return;
  }

  if (key == "cluster_min_dist")
  {
    cluster_min_dist = atof(value.c_str()); 
    return;
  }

  if (key == "range_algorithm")
  {
    range_algorithm = atoi(value.c_str()); 
    return;
  }

  if (key == "phi_algorithm")
  {
    phi_algorithm = atoi(value.c_str()); 
    return;
  }

  if (key == "north_angle")
  {
    north_angle = atof(value.c_str()); 
    return; 
  }

  if (key == "longitude")
  {
    longitude = atof(value.c_str()); 
    return; 
  }

  if (key == "latitude")
  {
    latitude = atof(value.c_str()); 
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
    delete top_cameras;
    parse_camera_orientations(value.c_str(), &top_cameras, &n_top_cameras);
    return; 
  }
  
  if (key == "cameras_orientation_bottom")
  {
    delete bottom_cameras; 
    parse_camera_orientations(value.c_str(), &bottom_cameras, &n_bottom_cameras);
    return; 
  }

  if (key=="mc")
  {
    mc = atoi(value.c_str())==1;
    return; 
  }
  if (key=="nokill")
  {
    nokill = atoi(value.c_str())==1;
    return; 
  }


  if (key == "channel_ids")
  {
    delete channel_ids; 
    std::cout << "parsing channel_ids" << std::endl; 
    parse_channel_ids(value.c_str(), &channel_ids, &n_channels_per_trigger); 
    return; 
  }

  if (key == "sat_thresh")
  {
    satThresh = atof(value.c_str());  
    return;
  }

  if (key == "n_sat_thresh")
  {
     nSatThresh = atoi(value.c_str()); 
     return;
  }
  if (key == "nbaseline")
  {
    parse_int_list(value.c_str(),&nbaseline_l, &nbaseline);   
    return;
   }
  if (key == "nblur")
  {
    parse_int_list(value.c_str(),&nblur_l,&nblur);  
    return; 
  }
  if (key == "invert_waveform")
  {
    parse_int_list(value.c_str(),&invert_waveform_l,&invert_waveform); 
    return; 
  }

  //If we get here, that means there was an unrecognized key
  std::cerr << "Unrecognized key \""<< key <<"\" on line " 
            << n << ". Skipping." << std::endl;

}


MaxCamClusterImage::CAMERA_ORIENTATION 
CleanSkimConfig::getCameraOrientation(int cam)
{
  for (int i = 0; i < n_top_cameras; i++)
  {
    if (top_cameras[i]==cam)
      return MaxCamClusterImage::TOP; 
  }

  for (int i = 0; i < n_bottom_cameras; i++)
  {
    if (bottom_cameras[i]==cam)
       return MaxCamClusterImage::BOTTOM; 
  }

  //If we got this far, that means the camera wasn't specified. 
  //Warn and return an up camera

  std::cerr << "WARNING: No orientation specified for camera " << cam << "." << std::endl;
  std::cerr << "  Assuming orientation is TOP" << std::endl; 
  return MaxCamClusterImage::TOP; 
}


void CleanSkimConfig::print(std::ostream & out)
{
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

  out << "burnin_num_events = " << burnin_num_events << std::endl;
  out << "burnin_distance_thresh = " << burnin_distance_thresh << std::endl;
  out << "outlier_factor = " << outlier_factor << std::endl;
  out << "blur_amount = " << blur_amount << std::endl;
  out << "cluster_min_sigma = " << cluster_min_sigma << std::endl;
  out << "cluster_max_sigma = " << cluster_max_sigma << std::endl;
  out << "cluster_min_size = " << cluster_min_size << std::endl;
  out << "cluster_min_dist = " << cluster_min_dist << std::endl;
  out << "range_algorithm = " << range_algorithm << std::endl;
  out << "phi_algorithm = " << phi_algorithm << std::endl;
  out << "north_angle = " << north_angle << std::endl;
  out << "latitude = " << latitude << std::endl;
  out << "longitude = " << longitude << std::endl;
  
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

  out << "mc = " << ((int)mc) << std::endl; 
  out << "nokill = " << ((int)nokill) << std::endl; 
 
  out << "sat_thresh = " << satThresh << std::endl; 
  out << "n_sat_thresh = " << nSatThresh << std::endl; 

  out << "nbaseline = ";
  for (int i = 0; i < nbaseline_l; i++) 
  {
    if (i > 0) out << ",";     
    out << nbaseline[i]; 
  }
  out << std::endl; 
  out << "nblur = ";
  for (int i = 0; i < nblur_l; i++) 
  {
    if (i > 0) out << ",";     
    out << nblur[i]; 
  }
  out << std::endl; 
  out << "invert_waveform = ";
  for (int i = 0; i < invert_waveform_l; i++) 
  {
    if (i > 0) out << ",";     
    out << invert_waveform[i]; 
  }
  out << std::endl; 
}

