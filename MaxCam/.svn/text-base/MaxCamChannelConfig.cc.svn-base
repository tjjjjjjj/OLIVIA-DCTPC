#include "MaxCamChannelConfig.hh"
#include <strings.h>
#include <string>
#include <istream>
#include <cstring>
#include <cstdlib>

static char * default_dev = "/dev/comedi0"; 
static char * default_serial = "/dev/ttyS4"; //For backwards compatibility
static char * default_dbaccess = "/etc/dbaccess.txt";



void MaxCamChannelConfig::init()
{
  _dev = default_dev; 
  _serial = default_serial; 
  _aref = DIFF;
  _readChannel = -1; 
  _writeChannel = -1; 
  _multiply = 1; 
  _readRange = 0; 
  _writeRange = 0; 
  _readSubdev = 0; 
  _writeSubdev = 1; 
  _calMethod = NONE; 
  _readOffset = 0; 
  _writeOffset = 0; 
  _dir = NOT_DIGITAL; 
  _calFile = NULL; 
  _alg = DIRECT; 
  _dv = 0.001; 
  _eps = 0.001; 
  _eps_zero = 0.005; 
  _dt = 200; 
  _invert = false; 
  _dbaccess = default_dbaccess;
  _hvtype = COMEDI; 
  _hardLimit = 1e99;
  _offset_wait = 3000; 
  _maxstep = -1; 
  
  registerNewString("device",&_dev); 

  char ** ref_keys = new char* [3];
  ref_keys[0] = "ground"; 
  ref_keys[1] = "common";
  ref_keys[2] = "diff"; 
  registerNewEnum("ref",3,ref_keys,NULL, &_aref); 

  char ** hv_keys = new char*[3]; 
  hv_keys[0] = "comedi"; 
  hv_keys[1] = "caen"; 
  hv_keys[2] = "mcc"; 
  registerNewEnum("hvType",3,hv_keys,NULL, &_hvtype); 

  registerNewInt("readChannel",&_readChannel); 
  registerNewInt("writeChannel",&_writeChannel); 
  registerNewUInt("readRange",&_readRange); 
  registerNewUInt("writeRange",&_writeRange); 
  registerNewUInt("readSubdev",&_readSubdev); 
  registerNewUInt("writeSubdev",&_writeSubdev); 

  char ** calib_keys = new char* [4];
  calib_keys[0]="none";
  calib_keys[1]="hardcal";
  calib_keys[2]="softcal",
  calib_keys[3]="offset"; 
  registerNewEnum("calibrationMethod",4,calib_keys,NULL,&_calMethod); 

  registerNewDouble("readOffset",&_readOffset); 
  registerNewDouble("writeOffset",&_writeOffset); 
  registerNewDouble("hardLimit",&_hardLimit); 
  registerNewString("calibrationFile",&_calFile); 
  
  char ** ramp_keys = new char * [2];
  ramp_keys[0] = "direct"; 
  ramp_keys[1] = "follow"; 
  registerNewEnum("rampAlgorithm",2,ramp_keys,NULL,&_alg); 

  char ** dio_keys = new char * [3]; 
  dio_keys[0] = "input"; 
  dio_keys[1] = "output"; 
  dio_keys[2] = "not_digital"; 
  registerNewEnum("dioDirection",3,dio_keys,NULL,&_dir); 

  registerNewDouble("dv",&_dv); 
  registerNewDouble("dt",&_dt); 
  registerNewBool("invert",&_invert); 
  registerNewDouble("epsilon",&_eps); 
  registerNewDouble("epsilonZero",&_eps_zero); 
  registerNewString("title",&_title); 
  registerNewString("serial",&_serial); 
  registerNewString("dbaccess",&_dbaccess); 
  registerNewDouble("writeMultiply",&_multiply); 
  registerNewDouble("offsetLimit",&_offset_limit); 
  registerNewDouble("offsetWait",&_offset_wait); 
  registerNewDouble("maxStep",&_maxstep); 
}

MaxCamChannelConfig::MaxCamChannelConfig(char * name, char * title)
{
  _name = name; 
  _title = title; 
  init(); 
}

MaxCamChannelConfig::~MaxCamChannelConfig()
{
 //TODO
}

MaxCamChannelConfig * MaxCamChannelConfig::loadConfig(const char * name, const char * dir)
{
  if (dir == NULL) dir = getenv("MCCHANNELCFGS"); 
  if (dir == NULL) std::cerr << "WARNING: MCCHANNELCFGS not defined" << std::endl; 
 
  char path[512]; 
  unsigned int dirlength = strlen(dir); 
  strcpy(path,dir); 
  path[dirlength] = '/'; 
  strcpy(path+1+dirlength, name); 

  MaxCamChannelConfig * ret = new MaxCamChannelConfig(strdup(name),strdup(name)); 
  ret->parseFile(path); 
  return ret; 

}

