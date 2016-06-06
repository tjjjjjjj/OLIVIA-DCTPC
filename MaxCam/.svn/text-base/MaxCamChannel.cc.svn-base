#include "MaxCamChannel.hh"
//#include "MaxCam.hh"
#include "MaxCamSerial.hh"

// ROOT includes
#include "TStopwatch.h"
#include "TSystem.h"
#include "TRegexp.h"
#include "TString.h"


//#include <string> //for ni6229info()

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::setw;
#include <termios.h>
//#include <fcntl.h>
//#include <unistd.h>
#include <assert.h>

#include <stdint.h>
#include <math.h>

#ifdef DM_DAQ
#include <linux/kernel.h>
#include <comedilib.h>
#include <mysql/mysql.h>
#endif


/** 
 * Used to read crededentials from
 * a dbaccess file
 */
namespace DBAccess {
  char server[256];///< The host that the database is on
  char user[256];///< username for the database
  char pass[256];///< password for the database
  char database[256];///< database name
  /** 
   * Set the database info fields from the given file.
   *
   * @param fname the file to read from
   */
  void set(char *fname) { 
    ifstream infile(fname);
    infile >> server >> user >> pass >> database;
  }
};

static uint16_t mcc_volts_to_adu(double v, bool bipolar)
{
   if (bipolar) 
    return  (uint16_t)  ((v/10. + 1.) * ((1 << 15)-1)); 

   return  (uint16_t)  ((v/10.) * ((1 << 16) -1)); 
}

#ifdef DM_DAQ
/* Returns true on success */ 
static bool connectToDB(MYSQL * mysql, const char * action, const char * name= "", int ntries = 10)
{
  int nloop =0; 
  while (!mysql_real_connect( mysql, DBAccess::server, DBAccess::user, 
            DBAccess::pass, DBAccess::database,0, NULL, 0 )) 
    {
      nloop++;
      cout <<name << ": Cannot connect after trials " << nloop << endl;
      gSystem->Sleep(100);
      if (nloop>ntries) {
  cout << name << ": Cannot connect to DB to "<<action << endl;
  return false;
      }
    }
  return true; 
}
#endif 




ClassImp(MaxCamChannel)
/* Class to interface between various hardware devices. Analog DAC's use comedilib, other devices
 * may use serial. */

void MaxCamChannel::init()
{

  if (cfg->getDBAccess()) DBAccess::set(cfg->getDBAccess()); 
#ifdef DM_DAQ

  mutex_db = 0;
  mutex_comedi = 0;
  mutex_usb = 0;
  mutex_caen = 0; 
  if (writeID >= 0 || readID >=0)
  {

    /* Apply hardware calibration if necessary */ 
    if (cfg->getCalibrationMethod()==HARDCAL)
    {
      char * calfile; 
      comedi_t * hw = comedi_open(cfg->getDevice()); 
      if (cfg->getCalibrationFile()==NULL)
      {
        calfile = comedi_get_default_calibration_path(hw); 
      }
      else calfile = cfg->getCalibrationFile(); 

      comedi_apply_calibration(hw,cfg->getReadSubdev(),cfg->getReadChannel(), cfg->getReadRange(),cfg->getRef(),calfile); 
      comedi_apply_calibration(hw,cfg->getWriteSubdev(),cfg->getWriteChannel(), cfg->getWriteRange(),cfg->getRef(),calfile); 
      comedi_close(hw); 
    } 

    /* Initialize software calibrator if necessary */ 
    if (cfg->getCalibrationMethod() == SOFTCAL)
    {
      calibrator = new MaxCamSoftCal((comedi_t *) comedi_open(cfg->getDevice()), cfg->getCalibrationFile());     
      cout << calibrator << endl; 
    }

    
    offset_guess_init = false; 
    offset_guess = 0; 
          
    //Try to load offset guess from database if following
    if (cfg->getRampAlgorithm() == FOLLOW)
    {
      MYSQL mysql; 
      mysql_init(&mysql); 
      bool success = connectToDB(&mysql, "load offset guesses", GetName()); 
      if (success)
      {
        char sql[200];
        sprintf(sql, "SELECT value FROM offset_guesses WHERE channel_name = \"%s\" ORDER BY timestamp DESC LIMIT 1", GetName()); 
        cout << sql << endl; 
        mysql_real_query(&mysql, sql, strlen(sql)); 
        MYSQL_RES * res = mysql_store_result(&mysql); 
        if (res == NULL) cout << "res is null" << endl; 
        if (res == NULL || !mysql_num_rows(res))
        {
          cout << "No offset guesses for " << GetName() << " in DB" <<endl; 
        }
        else
        {
          MYSQL_ROW row = mysql_fetch_row(res);
          offset_guess = atof(row[0]); 
          offset_guess_init = true; 
        }

        mysql_free_result(res); 
      }

      mysql_close(&mysql); 
    }

    if (cfg->getDigitalDirection()!=NOT_DIGITAL)        
    {
      comedi_t * device = comedi_open(cfg->getDevice()); 
      unsigned int dir = cfg->getDigitalDirection(); 
      comedi_dio_config(device, dir == INPUT ? cfg->getReadSubdev() : cfg->getWriteSubdev(),
                                dir == INPUT ? cfg->getReadChannel() : cfg->getWriteChannel(),
                                dir); 
      comedi_close(device); 

    }

    if (cfg->getHVType() == MCC3101)
    {
#ifdef DM_MCC
      if (!hid_is_initialised())
      {
        int ret = hid_init(); 
        printf("Ret: %d\n",ret); 
        if (ret != HID_RET_SUCCESS)
        {
          fprintf(stderr, "hid_init failed with return code %d\n", ret);
        }
      }

      if (PMD_Find_Interface(&(hid.hid), 0, USB3101_PID))
      {
        fprintf(stderr, "No interface found.\n");
        hid.hid = 0; 
      }
#endif
    }
  }

#endif
}


/* Legacy constructor that does not use MaxCamChannelConfig */ 
MaxCamChannel::MaxCamChannel(const char *name, const char *title, int idRead, int idWrite, char *dbparam, bool invert, MAXCAM_HV_TYPE_T type) : 
  TNamed(name, title),
  currentValue(0),
  currentValueRMS(0),
  setValue(0),
  readID(idRead),
  writeID(idWrite),
  invertPolarity(invert) {
  cfg = new MaxCamChannelConfig(strdup(name),strdup(title)); 
  cfg->setWriteChannel(writeID); 
  cfg->setReadChannel(readID); 
  cfg->setInvert(invert);
  cfg->setDBAccess(dbparam); 
  cfg->setHVType(type); 

  init();  
}



/* Copy constructor TODO: Implement copying of config */ 
MaxCamChannel::MaxCamChannel(const MaxCamChannel &other) :
  TNamed(other),
  currentValue(other.currentValue),
  currentValueRMS(other.currentValueRMS),
  setValue(other.setValue),
  readID( other.readID),
  writeID( other.writeID),
  invertPolarity( other.invertPolarity)
  
{}


/* Create a max cam channel based on a key. Depends on existence of 
   MCCHANNELCFGS environmental variable */ 
MaxCamChannel::MaxCamChannel(const char * key)
{
  cfg = MaxCamChannelConfig::loadConfig(key); 
  readID = cfg->getReadChannel();

  if(key=="temp0")
    cout << "TEMPERATURE readID: " << readID << endl;

  writeID = cfg->getWriteChannel();   
  invertPolarity = cfg->getInvert(); 
  fTitle = strdup(cfg->getTitle());   
  fName = key;
  init();  
}

MaxCamChannel::MaxCamChannel(MaxCamChannelConfig * c) 
{
  cfg = c;
  readID = cfg->getReadChannel(); 
  writeID = cfg->getWriteChannel();   
  invertPolarity = cfg->getInvert(); 
  fTitle = strdup(cfg->getTitle());   
  fName = strdup(cfg->getName()); 
  init(); 
}


void
MaxCamChannel::print() {
  cout << setw(12) << GetName() << ": " 
       << setw(12) << currentValue   << "  +/- " << setw(12) << currentValueRMS   
       << setw(12) << " set=" << setw(12) << setValue << endl;
}


int
MaxCamChannel::setDio(int val) {
  cout << "val = " << val << endl;


#ifdef DM_DAQ


  if (writeID<0) {
    cout << "MaxCamChannel: error writing digital output for channel " << writeID << endl;
    return -1;
  }

  /* FIXME:  ensure val is reasonable */
  if (val >= 1) val = 1;
  if (val <= 0) val = 0;
  if (val < 1) val = 0; // hmm...
  
  int ret;
  //int direction;

  char* deviceName = cfg->getDevice();
  int subdev = cfg->getWriteSubdev();
  cout << "device name = " << deviceName << endl;
  cout << "subdev = " << subdev << endl;

  comedi_t *dev = comedi_open(deviceName);
  if (!dev) {
    cout << "MaxCamChannel: error opening dio device " << deviceName << endl;
    return -1;    
  }

  /* FIXME:  add check that this channel is configured to be digital out */
  //ret = comedi_dio_get_config(dev, subdev, writeID, &direction);
  //if (direction != COMEDI_OUTPUT) {
  //  cout << "MaxCamChannel: error, need to configure digital channel as output." << writeID << endl;
  //  return -1;
  //}

  // and finally set the value
  ret = comedi_dio_write(dev, subdev, writeID, val);
  return ret;

#else
  cout << GetName() << ": DB not used " << endl;
  return 0;
#endif
  

  }


int
MaxCamChannel::rampFromDB(const char *table) {
  // Set voltage using setval field from database.
  // See also readValueFromDB and ramp. This function can be used
  // on slow-control machine with ADC connections.
  //
    
#ifdef DM_DAQ
  if (!table) table=GetName();

  // get value from DB
  double newsetval = readValueFromDB(table, 2);


// make lock
  MYSQL mysql;
  bool lock = mutex_db != 0; 
  if (lock) mutex_db->Lock(); 
  mysql_init(&mysql);
  bool success = connectToDB(&mysql,"ramp HV", GetName(), 1); 
  if (!success) 
  {
    if (lock) mutex_db->UnLock(); 
    return -1; 
  }

  TString sqlLock=TString("LOCK TABLES") + table + TString(" READ");
  mysql_real_query(&mysql, (const char*)sqlLock, sqlLock.Length());

  if (lock) mutex_db->UnLock(); 

  // talk to ADC
  ramp(newsetval);


  // save
  if ( saveValueToDB(table) ) {
  cout << "MaxCamChannel " << GetName() <<": cannot write setValue into DB" << endl;   
  return -1;  
  }

  // wait  
  // WHY?
//  gSystem->Sleep(1500);

  // release lock
  TString sqlUnlock=TString("UNLOCK TABLES");
  if (lock) mutex_db->Lock(); 
  mysql_real_query(&mysql, (const char*)sqlUnlock, sqlUnlock.Length());
  mysql_close(&mysql);
  if (lock) mutex_db->UnLock(); 
#else
    cout << GetName() << ": DB not used " << table << endl;
#endif

  return 0;
}


int
MaxCamChannel::ramp(double newVal, bool checkValue, bool * interrupt) {
  // Set voltage.

#ifdef DM_DAQ
 
  if (writeID<0) {
    cout << "MaxCamChannel " << GetName() 
   << ": cannot write to channel with id = " << writeID << endl;
    return -1;
  }

  if (newVal > cfg->getHardLimit())
  {
    newVal = cfg->getHardLimit(); 
    cout << "Warning: newVal larger than hard limit. newVal clipped to " << cfg->getHardLimit() << endl; 
  }

  readADC();


  double delta = newVal - currentValue; 

  //If already at desired value, do not continue to ramp
  if (newVal == 0 
        ? fabs(delta) < cfg->getEpsilonZero()
        : fabs(delta) < cfg->getEpsilon()
     )
  {
    return 0; 
  }



  if (cfg->getMaxStep() > 0 && cfg->getMaxStep() < fabs(delta))
  {
    int nsteps = int(delta / cfg->getMaxStep() + 1); 
    double step_val = delta/nsteps; 
    double initValue = currentValue; 
    for (int step = 0; step < nsteps-1; step++)
    {
      ramp((step+1) * step_val + initValue, checkValue, interrupt); 
      gSystem->Sleep((int)cfg->getDT()); 
      saveValueToDB(GetName(), 0, false); 
      if (interrupt && *interrupt) return 1; 
    }
  }
  
  //cout << "Value for " << GetName() << " new=" << newVal << "   set="  << setValue << "  val=" << currentValue<<endl;
  setBusy(true); 
  if (cfg->getRampAlgorithm() == DIRECT)
    {
      setValue = newVal;
      setDAC();
    
      TStopwatch sw;
      sw.Start();
      while(checkValue) {
        if (fabs(readADC()-setValue)<0.01) break;
        gSystem->Sleep((int)cfg->getDT()); 
        //cout << setValue << "  " << readADC() << endl;
        if (sw.RealTime()>10) {
            cout << "MaxCamChannel " << GetName() 
            << ": setting voltage: " << currentValue 
            << " -> " <<  setValue << "  " << endl;
            sw.Start();
        }
        else sw.Start(0);


        if (interrupt && *interrupt)
        {
          cout << GetName() << " Interrupted!" << endl; 
          break; 
        }
      }
    }

    else if (cfg->getRampAlgorithm() == FOLLOW)
    {
      if (!offset_guess_init)
        setValue= 0;  
      else setValue= newVal + offset_guess; 

      
      if (mutex_db) mutex_db->Lock(); 
      setDAC();
      if (mutex_db) mutex_db->UnLock(); 

      if (offset_guess_init) 
      {
        double waitAmt = cfg->getOffsetWait(); 
        int nintervals = int(waitAmt/500) + 1; 
        for (int interval = 0; interval < nintervals; interval++)
        {
          gSystem->Sleep((unsigned)waitAmt/nintervals); 
          saveValueToDB(GetName(), 0, false); 
        }
      }


      readADC(); 

      bool save = true; 

      while (newVal == 0 
                ? fabs(currentValue - newVal) > cfg->getEpsilonZero() 
                : fabs(currentValue - newVal) > cfg->getEpsilon()
           )
      {



        cout << currentValue << " : " << newVal << " : " << setValue<<  endl;
        if (currentValue > newVal) 
        { 
          setValue-= cfg->getDV();    
        }
        else
        {
            setValue+= cfg->getDV(); 
        }

        if (setValue< -1.) setValue= -1; 

        if (setValue> cfg->getHardLimit()) setValue= cfg->getHardLimit(); 

        if (mutex_db) mutex_db->Lock(); 
        setDAC();
        if (mutex_db) mutex_db->UnLock(); 
        saveValueToDB(GetName(),0,false);  
        gSystem->Sleep((int)cfg->getDT()); 
        readADC(); 

        if (cfg->getOffsetLimit() > 0 && fabs(setValue - newVal) > cfg->getOffsetLimit() && newVal != 0)
        {
          cout << GetName() << ": Offset limit exceeded. Giving up! " <<endl; 
          save = false; 
          break; 
        }

        if (interrupt && *interrupt)
        {
          cout << GetName() << " Interrupted!" << endl; 
          save = false;
          break; 
        }

      }

      if(save && newVal != 0 && setValue!= cfg->getHardLimit())
      {
        double last_offset_guess = 0; 
        if (offset_guess_init == true) last_offset_guess = offset_guess; 
        offset_guess_init = true; 
        offset_guess = setValue-currentValue; 

        //Store new offset in db if different by more than half volt from last offset guess
        if (fabs(offset_guess - last_offset_guess) > 0.0005)
        {
            saveOffsets(); 
        }
      }
      setValue = newVal; 
    }
    setBusy(false); 
#else 
    cout << GetName() << ": Cannot set new value " << newVal << endl;
    
#endif

    return 0;
  }



int MaxCamChannel::saveOffsets()
{
  cout << "Saving offsets" << endl; 
#ifdef DM_DAQ
  MYSQL mysql;  
  bool lock = mutex_db!=0; 
  if (lock) mutex_db->Lock(); 
  mysql_init(&mysql); 
  bool success = connectToDB(&mysql,"save offsets", GetName()); 

  if (success)
    {
      char sql[200]; 
      sprintf(sql,"INSERT INTO offset_guesses (value, channel_name, timestamp) VALUES ( %f , \"%s\" , NOW());", offset_guess, GetName());
      cout << sql << endl; 
      int result = mysql_real_query(&mysql,sql,strlen(sql)); 
      cout << "Result of query is " << result << endl; 
    }
  mysql_close(&mysql); 
  if(lock) mutex_db->UnLock(); 

  if (success) return 0; 
#endif
  return 1; 
}

int
MaxCamChannel::setDAC(char *deviceName) {
  // Set DAC channel. DAS1200 has two output channels 
  // channel       pins
  //       0      35-36
  //       1      37-38

  // NI6229 Has 4 output channels 

  int ret=0;

  if (setValue > cfg->getHardLimit()) setValue = cfg->getHardLimit(); 


#ifdef DM_DAQ
  if (cfg->getHVType() == COMEDI)
  {
    int range = cfg->getWriteRange();    
    int aref = cfg->getRef();  
    int subDevice = cfg->getWriteSubdev(); // 1 for DAC with two output channels

    //TODO: Fix this check 
    if (writeID<0 || writeID>3) {
      cout << "MaxCamChannel: error setting voltage for channel " << writeID << endl;
      return -1;
    }

    lsampl_t data;
    comedi_t *cf;
    bool lock = mutex_comedi != 0; 
    if (lock) mutex_comedi->Lock(); 
    if (deviceName==NULL) cf = comedi_open(cfg->getDevice());
    else cf  = comedi_open(deviceName);

    if (!cf) {
      cout << "MaxCamChannel: error opening DAC device " << deviceName << endl;
      if (lock) mutex_comedi->UnLock(); 
      return -1;    
    }

    int maxdata = comedi_get_maxdata(cf, subDevice, writeID);
    comedi_range *cr = comedi_get_range(cf, subDevice, writeID, range);
    if (!cr) {
      cout << "MaxCamChannel: error getting range for channel " << writeID << endl;
      if (lock) mutex_comedi->UnLock(); 
      return -1;    
    }

    if (cfg->getCalibrationMethod() == SOFTCAL)
    {
      data = calibrator->fromPhysical(setValue * cfg->getWriteMultiply(), cfg->getWriteSubdev(), cfg->getWriteChannel(), cfg->getWriteRange());
    } 
    else if (cfg->getCalibrationMethod() == OFFSET)
    {
      data = comedi_from_phys(setValue * cfg->getWriteMultiply() + cfg->getWriteOffset(), cr, maxdata);
    }
    else
    {
      data = comedi_from_phys(setValue * cfg->getWriteMultiply(), cr,maxdata); 
    }

    ret = comedi_data_write(cf, subDevice, writeID, range, aref, data);

    comedi_close(cf);
    if (lock) mutex_comedi->UnLock(); 
  }
  else if (cfg->getHVType() == CAEN)
  {
    char cmd[128]; 
    sprintf(cmd,"$BD:0,CMD:SET,CH:%d,PAR:VSET,VAL:%06.1f", writeID, setValue * cfg->getWriteMultiply()); 
    rawCAENcmd(cmd,NULL,0.1,deviceName);    
  }
  else if (cfg->getHVType() == MCC3101)
  {
#ifdef DM_MCC                
    bool lock = mutex_usb != 0; 
    if (lock) mutex_usb->Lock(); 
    std::cout << "Multiply factor: " << cfg->getWriteMultiply() << std::endl; 
    uint32_t write = mcc_volts_to_adu( cfg->getCalibrationMethod() == OFFSET ? setValue  * cfg->getWriteMultiply()+ cfg->getWriteOffset() : setValue * cfg->getWriteMultiply(),
                                       cfg->getWriteRange()); 
    usbAOutConfig_USB31XX(hid.hid,cfg->getWriteChannel(), cfg->getWriteRange()); 
    usbAOut_USB31XX(hid.hid, cfg->getWriteChannel(), write, 0); 
                  
    if (lock) mutex_usb->UnLock(); 

#else
     cerr << "NOT COMPILED WITH DM_MCC!!!!" << endl; 
#endif

  }
  else cout << GetName() << ": device not used " << deviceName << endl;
#endif 

  return ret;
}


double
MaxCamChannel::readADC(int nread, char *deviceName) {
  // Read ADC channel. DAS1200 has 16 analog input channels specified 
  // by the first argument (0-15) for connections to pins 2-17. 
  // 'nread' is the number of readout that are averaged.
  // 'deviceName' is the linux device (default /dev/comedi0).

#ifdef DM_DAQ
  if (cfg->getHVType() == COMEDI || cfg->getHVType() == MCC3101)
  {
    if (readID<0 || readID>32) {
      cout << "MaxCamChannel: error reading voltage for channel " << readID << endl;
      return -1;
    }

    int subDevice = cfg->getReadSubdev();
    int range = cfg->getReadRange(); 
    // not default, but a discrepancy between what Cosmin's new addition has and
    // what James has, which seems to work
    //  int aref = cfg->getRef();  
    //  int aref=AREF_GROUND;
    int aref=cfg->getRef();

    lsampl_t data;
    if (deviceName == NULL) deviceName = cfg->getDevice(); 
    bool lock = mutex_comedi !=0; 
    if (lock) mutex_comedi->Lock(); 
    comedi_t *cf = comedi_open(deviceName);
    int maxdata = comedi_get_maxdata(cf, subDevice, readID);
    comedi_range *cr = comedi_get_range(cf, subDevice, readID, range);

    currentValue=0, currentValueRMS=0;
    for (int i=0; i<nread; i++) {
      comedi_data_read(cf, subDevice, readID, range, aref, &data);
      double volts; 
      if (cfg->getCalibrationMethod() == SOFTCAL)
      {
        volts = calibrator->toPhysical(data, subDevice, readID, range);    
      }
      else
      {
        volts = comedi_to_phys(data, cr, maxdata);
        if (cfg->getCalibrationMethod() == OFFSET)
          volts += cfg->getReadOffset(); 
      }
      
      currentValue += volts;
      currentValueRMS += volts*volts;
    }
    currentValue/=nread;

    currentValueRMS/=nread;
    currentValueRMS = sqrt( currentValueRMS-currentValue*currentValue); 
    currentValueRMS/=(nread-1);
   
    comedi_close(cf);
    if (lock) mutex_comedi->UnLock(); 
  }
  else if(cfg->getHVType() == CAEN)
  {
    // get the value of the hv
    char cmd[128]; 
    sprintf(cmd,"$BD:0,CMD:MON,CH:%d,PAR:VMON", readID); 
    char * result; 
    rawCAENcmd(cmd,&result,0.1,deviceName); 
    char * val = strstr(result,"VAL:"); 
    if (val == NULL)
     {
      //Handle error so we don't seg fault
      return -1; 
     }
    val+=4; //4 is strlen("VAL:"); 

    // get the polarity
    // sleep a little bit so we get a good read
    gSystem->Sleep(100);
    sprintf(cmd,"$BD:0,CMD:MON,CH:%d,PAR:POL", readID);
    rawCAENcmd(cmd,&result,0.1,deviceName);
    char * pol = strstr(result,"VAL:");
    if (pol == NULL)
      {
  //Handle error so we don't seg fault                                                                                                                   
  return -1;
      }
    pol+=4; //4 is strlen("POL:");   
    cout << "pol=" << pol << endl;
    if(pol[0]=='+') currentValue = atof(val);
    else if(pol[0]=='-') currentValue = -atof(val);
    else return -1;
    
    currentValueRMS = 0; // Is there some good way to estimate this?
  }
#else
  cout << GetName() << ": device not used " << deviceName << "  " << nread << endl;

#endif

  if (isPolarityInverted()) currentValue=-currentValue;

  return currentValue;
}


bool
MaxCamChannel::subdevIsAnalogIn(int subdev, char *deviceName){
  bool retval = false;
#ifdef DM_DAQ
  comedi_t *dev=comedi_open(deviceName);
  int subdevtype=comedi_get_subdevice_type(dev, subdev);
  if (subdevtype==COMEDI_SUBD_AI){
    retval = true;
    return retval;
  }
  else {
    retval = false;
    return retval;
  }
#endif
  return retval;
}

bool
MaxCamChannel::subdevIsAnalogOut(int subdev, char *deviceName){
  bool retval = false;
#ifdef DM_DAQ
  comedi_t *dev=comedi_open(deviceName);
  int subdevtype=comedi_get_subdevice_type(dev, subdev);
  if (subdevtype==COMEDI_SUBD_AO){
    retval = true;
    return retval;
  }
  else {
    retval = false;
    return retval;
  }
#endif
  return retval;
}


double
MaxCamChannel::readConvectron(char *port_name, bool verbose) {
  if (port_name == NULL)
  {
      port_name = cfg->getSerial(); 
  }
  long baudRate = 19200;
  int dataBits  = 8;
  char parity   = 'n';
  int stopBits  = 1;
  int xonxoff   = 0;
  int rtscts    = 0;

  SerialHandle rs232 = MaxCamSerial::Open(port_name);
  MaxCamSerial::Init(rs232, baudRate, dataBits, parity, stopBits, 
         xonxoff, rtscts);
  MaxCamSerial::Send(rs232, "RD\r", 0.1);
  int readlength = 255;
  bool trimWhitespace = true;
  string readstr = MaxCamSerial::ReadPort(rs232, readlength, trimWhitespace);
  MaxCamSerial::Close(rs232);

  if (verbose) {
    cout << "readstr = [" << readstr << "]" << endl;
    cout << "readstr.length() = " << readstr.length() << endl;
    cout << "readstr.size()   = " << readstr.size() << endl;
    cout << endl;
  }

  if (!readstr.compare("SNSR UNP")) {
    if (verbose)
      cout << "no sensor" << endl;
    return -10;
  } else if (!readstr.compare("SNSR OVP")) {
    if (verbose)
      cout << "sensor overpressure" << endl;
    return -20;
  } else if (!readstr.compare("OPN SNSR")) {
    if (verbose)
      cout << "defective transducer" << endl;
    return -30;
  }

  // otherwise, convert the string to a double
  if (verbose)
    cout << "valid pressure reading" << endl;
  // convert c-style string to float
  double pressure = atof(readstr.c_str());

  // find a better fix for this that lets the user
  // know that there's been a problem
  //  if (readstr.length()>0) currentValue = pressure; 
  return pressure;
}

double
MaxCamChannel::readInficonController(char *port_name, int channel, long baudRate) {
 if (port_name == NULL)
    {
      port_name = cfg->getSerial(); 
    }

  int dataBits  = 8;
  char parity   = 'n';
  int stopBits  = 1;
  int xonxoff   = 0;
  int rtscts    = 0;

  char ENQ = 0x05;  // "enquire" for the data after an ACK

  SerialHandle rs232 = MaxCamSerial::Open(port_name);
  MaxCamSerial::Init(rs232, baudRate, dataBits, parity, stopBits, 
         xonxoff, rtscts);  
  MaxCamSerial::Flush(rs232);

  string cmdStr;
  if (channel == 1) {
    cmdStr = "PR1\r";
  } else if (channel == 2) {
    cmdStr = "PR2\r";
  } else if (channel == 3) {
    cmdStr = "PR3\r";
  } else {
    cout << "error in MaxCamChannel::readInficonController()" << endl;
    cout << "   unrecognized channel number " << channel << endl;
    cout << "   returning -1" << endl;
    return -1;
  }
  char cmd[1024];
  strcpy(cmd, cmdStr.c_str());
  //MaxCamSerial::Send(rs232, cmd, 0.1);
  int bytes_sent;
  bytes_sent = MaxCamSerial::SerialWrite(rs232, cmd, 4);
  //cout << "bytes_sent = " << bytes_sent << endl;
  usleep(100000); // this is very necessary!

  inficonGetAck(rs232);

  // send ENQ
  bytes_sent = MaxCamSerial::SerialWrite(rs232, &ENQ, 1);
  //MaxCamSerial::Send(rs232, &ENQ, 0.1);

  // Read response
  // for some reason, MaxCamSerial::ReadPort(rs232,64) occasionally fails
  // so I am doing the read "manually" using MaxCamSerial::SerialRead()
  //
  //string sensorString = MaxCamSerial::ReadPort(rs232, 64);
  char c[1024];
  int ncount = 0;
  string sensorString;
  int totaltries=0;
  sensorString.clear();
  while(1) {
    int nread = MaxCamSerial::SerialRead(rs232, &c[ncount], 1);
    ++totaltries;
    if(totaltries>1e6){
      break; // give up if it fails 100x or more
    }
    if (nread<1) continue;
    sensorString += c[ncount];
    if (++ncount==13) break;  // the response from the VGC402 is 13 chars long
  }
  //  cout << "sensorString.size(), val = [" << sensorString.size() << "], [" 
  //       << sensorString << "]" << endl;
  
  // wrap up
  MaxCamSerial::Close(rs232);

  // For some reason, the response is sometimes truncated.
  // So, instead of getting e.g. 2.047E+02 (torr), you just get 2.047
  // One way around this is to check the length of the response
  // The string should be 13 characters long
  //    <PM>B.BBBBE<PM>BB 
  // 
  // But the following code could easily lead to an infinite loop!
  // need a better way to deal with this...
  //size_t expectedSize = 13;
  //if (sensorString.size() != expectedSize) {
  //  cout << "failed!!! retrying" << endl;
  //  return readInficonController(port_name, channel, baudRate);
  //}

  // should probably verify that A=0
  string statusStr = sensorString.substr(0, sensorString.find_first_of(","));
  //  cout << "sensorString.find_first_of(',') = " << sensorString.find_first_of(",") << endl;
  //  cout << "statusStr = " << statusStr << endl;
  int status = atoi(statusStr.c_str());
  enum { MSRMT_OK, UNDERRANGE, OVERRANGE, 
   SENSOR_ERROR, SENSOR_OFF, NO_SENSOR, 
   ID_ERROR, BPG_HPG_ERROR };
  if (status == UNDERRANGE)
    return -10*UNDERRANGE;
  else if (status == OVERRANGE)
    return -10*OVERRANGE;
  else if (status == SENSOR_ERROR)
    return -10*SENSOR_ERROR;
  else if (status == SENSOR_OFF)
    return -10*SENSOR_OFF;
  else if (status == NO_SENSOR)
    return -10*NO_SENSOR;
  else if (status == ID_ERROR)
    return -10*ID_ERROR;
  else if (status == BPG_HPG_ERROR)
    return -10*BPG_HPG_ERROR;

  // then extract the entire substring following the comma
  TString valueStr = sensorString.substr(sensorString.find_first_of(",")+1);
  //  cout << "valueStr = " << valueStr << endl;

  // need to see if the value that we read matches the pattern of what
  // we expect to read; if it doesn't, it's a bad read, and we need to throw
  // it away.  Use regular expressions.
  TRegexp expected_inficon_value_regexp_pattern("[0-9]\\.[0-9][0-9][0-9][0-9]E[+-][0-9][0-9]");

  // and return the value (<PM>B.BBBBE<PM>BB) converted to a double
  //  cout << "valueStr.Data() = " << valueStr.Data() << endl;
  //  cout << "valueStr.Contains(expected_inficon_value_regexp_pattern)="
  //       << valueStr.Contains(expected_inficon_value_regexp_pattern)
  //       << endl;
  //  cout << "valueStr(expected_inficon_value_regexp_pattern)="
  //       << valueStr(expected_inficon_value_regexp_pattern)
  //       << endl;
  if(totaltries<1e6 && valueStr.Contains(expected_inficon_value_regexp_pattern)){
    currentValue=atof(valueStr(expected_inficon_value_regexp_pattern).Data());
  } else {
    // don't change the current value, so that the value won't change in 
    // the database; for now this is probably the best option ...
  }
  currentValueRMS=0;
  return currentValue;
}

int
MaxCamChannel::inficonGetAck(SerialHandle rs232) {
  // reads the <ACK> byte from the inficon controller
  // and then reads two more bytes (<CR><LF>)
  // in preparation for data transfer

  char ACK = 0x06;  // "acknowledge" sent back by the controller
  //char NAK = 0x15;  // "negative acknowledge" sent by controller

  char readbuf[256];
  int  readlength = 1;
  int  nloop = 0;
  int  bytes_read = 0;
  while (1) {
    nloop++;
    //printf("readport\n");
    //MaxCamSerial::ReadPort(rs232, readbuf, readlength);
    bytes_read = MaxCamSerial::SerialRead(rs232, readbuf, readlength);
    //if (bytes_read > 0) {
    //  printf("bytes_read [%d] ", bytes_read);
    //  printf("%03d, [0x%x]\n", nloop, (unsigned int)readbuf[0]);
    //}
    if ( (bytes_read && readbuf[0] == ACK) || nloop>100) {
      //printf("readbuf done\n");
      // read the <CR> and <LF> and return
      MaxCamSerial::ReadPort(rs232, readbuf, 2);
      break;
    }
  }

  return 0;
}

double
MaxCamChannel::readSerialInficon(char *port_name) {
  // Function reads serial port that receives stream of bytes with
  // measurements of pressure from Inficon pressure gauge. 


  if (port_name == NULL)
    {
      port_name = cfg->getSerial(); 
    }
  int fd;

  //----------------------------------------------------------
  fd = open(port_name, O_RDWR | O_NOCTTY | O_NDELAY);
  if (fd < 0) 
    cout<<GetName() << ": Could not connect"<<endl;
  //else 
  //cout<<GetName() << ": Connected to "<<port_name<<endl;
  //----------------------------------------------------------

  struct termios options;
  if (tcgetattr(fd, &options)) {
    cout << GetName() << ": tcgetattr failed " << endl;
    return -1;
  }
  options.c_cflag |= (CLOCAL | CREAD);
  options.c_cflag &= ~PARENB; /* Clear parity enable */
  options.c_cflag &= ~CSTOPB;
  options.c_cflag &= ~CSIZE;
  options.c_cflag |= CS8;
  options.c_cflag &= ~CRTSCTS;
  options.c_cflag |= (IXON | IXOFF); /* Software flow control */
  options.c_cflag |= HUPCL; /* Drop DTR on close */

  options.c_iflag &= ~INPCK;  /* Enable parity checking */

  options.c_cc[VMIN] = 0;
  options.c_cc[VTIME] = 10;

  options.c_lflag &= ~(ICANON | ECHO | ECHOE | ISIG);
  options.c_lflag = 0; /* no local flags */

  options.c_oflag &= ~(IXON | IXOFF | IXANY); /* no flow control */
  options.c_oflag &= ~OPOST; /* No output processing */
  options.c_oflag &= ~ONLCR; /* Don't convert linefeeds */
  
  cfsetospeed(&options, B9600);
  
  /* Clear the line */
  tcflush(fd,TCIFLUSH);



  // Set the new options for the port...
  if (tcsetattr(fd, TCSANOW, &options)) {
    cout << "tcsetattr failed\n" << endl;
    close(fd);
    return -1;
  }


  char PR1[]="PR1\r";
  int bytes_send = write(fd, PR1, 4 );
  assert (bytes_send==4);

  char c[1024];
  int bytes_read=0;
  while (1) {
    bytes_read = read(fd, c, 1 );
    if (bytes_read && c[0]==6) break;
  }
  
  std::string buffer;
  buffer.clear();
  char ENQ=5;
  bytes_send = write(fd, &ENQ, 1 );
  assert (bytes_send==1);

  int ncount=0;
  while (1) {
    int nread = read(fd, &c[ncount], 1 );
    if (nread<1) continue;
    buffer += c[ncount];
    if (++ncount==15) break;
  }

  int flg;
  sscanf( buffer.c_str(), "%d, %le", &flg, &currentValue);
  switch(flg) {
    //case 0: cout << GetName() << ": Current value is "     << currentValue << endl; break;
  case 1: cout << GetName() << ": Underrange, P set to " << currentValue << endl; break;
  case 2: cout << GetName() << ": Overrange, P set to "  << currentValue << endl; break;
  }
  currentValueRMS=0;
  close(fd);

  return currentValue;
}



double
MaxCamChannel::readSerialLVG200TC(char *port_name) {
  // Function reads serial port that receives stream of bytes with
  // measurements of pressure from LVG-200TC pressure gauge. 

  int fd;

  //----------------------------------------------------------
  fd = open(port_name, O_RDWR | O_NOCTTY | O_NDELAY);
  if (fd < 0) 
    cout<<GetName() << ": Could not connect"<<endl;
  else 
    cout<<GetName() << ": Connected to "<<port_name<<endl;
  //----------------------------------------------------------

  struct termios options;
  if (tcgetattr(fd, &options)) {
    cout << GetName() << ": tcgetattr failed " << endl;
    return -1;
  }
  options.c_cflag |= (CLOCAL | CREAD);
  options.c_cflag &= ~PARENB; /* Clear parity enable */
  options.c_cflag &= ~CSTOPB;
  options.c_cflag &= ~CSIZE;
  options.c_cflag |= CS8;
  options.c_cflag &= ~CRTSCTS;
  options.c_cflag |= (IXON | IXOFF); /* Software flow control */
  options.c_cflag |= HUPCL; /* Drop DTR on close */

  options.c_iflag &= ~INPCK;  /* Enable parity checking */

  options.c_cc[VMIN] = 0;
  options.c_cc[VTIME] = 10;

  options.c_lflag &= ~(ICANON | ECHO | ECHOE | ISIG);
  options.c_lflag = 0; /* no local flags */

  options.c_oflag &= ~(IXON | IXOFF | IXANY); /* no flow control */
  options.c_oflag &= ~OPOST; /* No output processing */
  options.c_oflag &= ~ONLCR; /* Don't convert linefeeds */
  
  cfsetospeed(&options, B9600);
  
  /* Clear the line */
  tcflush(fd,TCIFLUSH);



  // Set the new options for the port...
  if (tcsetattr(fd, TCSANOW, &options)) {
    cout << "tcsetattr failed\n" << endl;
    close(fd);
    return -1;
  }
  else {
    cout << "Port initialized" << endl;
  }
  
  char c[32];
  bool hasStarted=false;
  std::string buffer;
  int nloop=0;
  while (1) {
    int nread = read(fd, c, 1 );
    if (nloop++>10000) {
      cout << GetName() <<": No data on serial port!!" << endl;
      break;
    }
    if (nread<1) continue;
    if ( c[0]=='\n' ) {
      hasStarted=true;
      buffer.clear();
    }
    if (hasStarted && c[0]=='a') break;

    buffer += c[0];
  }  
   
  currentValue=atof( buffer.c_str() );
  currentValueRMS=0;
  close(fd);

  return currentValue;
}

void
MaxCamChannel::readFromDB(const char *table, bool doPrint, MY_MYSQL*  handle) {
  // Update all values from DB filled by slow control.
  // Call this function from a process that is taking data.
  //

  if (!table) table=GetName();

  currentValue    = readValueFromDB(table,0,doPrint,handle);
  //currentValueRMS = readValueFromDB(table,1,doPrint);
  //setValue        = readValueFromDB(table,2,doPrint);
}



float
MaxCamChannel::readValueFromDB(const char *table, unsigned int ifield, bool doPrint, MY_MYSQL* handle) {
  // A primitive access to a DB containing slow control data.
  // Only look for the last row which has the following fields:
  // 0=value
  // 1=rms
  // 2=preset value
  //

  float ret=0;

#ifdef DM_DAQ

  MYSQL *  mysql;
  MYSQL mysql_struct; 
  
  bool lock = mutex_db !=0; 
  if (lock) mutex_db->Lock(); 

  if (!handle)
  {
    mysql = &mysql_struct;
    mysql_init(mysql);
    bool success = connectToDB(mysql,"save value", GetName()); 
    if (!success) 
    {
      if (lock) mutex_db->UnLock(); 
      return -1; 
    }
  }
  else
  {
    mysql =(MYSQL*) handle;
    mysql_ping(mysql);
  }

  TString sql=TString("SELECT * FROM ") + table + TString(" ORDER BY timestamp DESC LIMIT 1");
  mysql_real_query(mysql, (const char*)sql, sql.Length());

  //cout << sql << endl;

  MYSQL_RES *res=mysql_store_result(mysql);

  //cout << "table = " << table;
  if (!mysql_num_rows(res)) {
    cout << GetName() << ": No entries in table " << table << endl;
    if (!handle)
      mysql_close(mysql);
    if (lock) mutex_db->UnLock(); 
    return 0;
  }

  MYSQL_ROW  row=mysql_fetch_row(res);

  if (ifield<0 || ifield>mysql_num_fields(res)-1) { 
    cout << "MaxCamChannel: from field = " << ifield << endl;
    if (!handle)
      mysql_close(mysql);
    if (lock)
      mutex_db->UnLock(); 
    return -1; 
  }

  unsigned int  num_fields = mysql_num_fields(res);
  MYSQL_FIELD *fields = mysql_fetch_fields(res);

  currentValue = atof( row[0] );
  currentValueRMS = atof( row[1] );
  setValue = atof( row[2] );
  ret = atof( row[ifield] );
  
  if (doPrint) {
    cout << setw(12) << table << ":";
    for(unsigned int i = 0; i < num_fields; i++) {
      cout <<  setw(12) << fields[i].name << "=" << setw(12) << row[i] ; 
    }
    cout << endl;
  }

  mysql_free_result(res);
  if (!handle)
    mysql_close(mysql);
  if (lock) mutex_db->UnLock(); 
#else
  cout << GetName() << ": DB not used " << table << "  " << ifield << endl;

#endif

  return ret;
}


int
MaxCamChannel::saveValueToDB(const char *table, MY_MYSQL* handle, bool save_set) {
  // A primitive access to a DB containing slow control data.

  if (!table) table=GetName();

#ifdef DM_DAQ

  MYSQL *  mysql;
  bool lock = (mutex_db!=0); 
  if (lock) mutex_db->Lock(); 
  MYSQL mysql_struct; 
  
  if (!handle)
  {
    mysql = &mysql_struct; 
    mysql_init(mysql);
    bool success = connectToDB(mysql,"save value",GetName());
    if (!success)
    {
      mutex_db->UnLock(); 
      return 1; 
    }
  }
  else
  {
    mysql = (MYSQL*)handle; 
    mysql_ping(mysql);
  }

  TString sql = ""; 
  double set = setValue; 
  if (!save_set)
  {
    sql = TString("SELECT setval FROM ") + table + TString(" ORDER BY timestamp DESC LIMIT 1"); 
    mysql_real_query(mysql, sql.Data(), sql.Length()); 
    MYSQL_RES * res = mysql_store_result(mysql); 
    MYSQL_ROW row = mysql_fetch_row(res); 
    set = atof(row[0]); 
    mysql_free_result(res); 
  }

  sql=TString("INSERT INTO ") + table 
    + TString(" ( value, rms, setval, timestamp) VALUES( ");

  sql += currentValue;
  sql += TString(", ");
  sql += currentValueRMS;
  sql += TString(", ");
  sql += set;
  sql += TString(", NOW() )");
       
  //cout << sql << endl;

  mysql_real_query(mysql, (const char*)sql, sql.Length());

  if (!handle)
    mysql_close(mysql);
  if (lock) mutex_db->UnLock(); 
#else
  cout << GetName() << ": DB not used " << table << endl;

#endif

  return 0;
}





int
MaxCamChannel::saveValueToDBAndCheck(const char *table,MY_MYSQL* handle) {

  if (!table) table=GetName();

  double currentValue0=currentValue;
  double currentValueRMS0=currentValueRMS; 
  double setValue0=setValue;
 
  int ret=0;
  double delta=1e-5;
  while (1) {

    ret=saveValueToDB(table,handle);
    
    readFromDB(table,handle);
    
    if ( fabs(currentValue-currentValue0)<delta &&
   fabs(currentValueRMS-currentValueRMS0)<delta &&
   fabs(setValue-setValue0)<delta  ) break;

    cout << GetName() << ": Setting values again (delta=" << delta << ")" << endl;
    if (fabs(currentValue-currentValue0)>delta)        
      cout << currentValue << "!=" << currentValue0 << "  delta=" << fabs(currentValue-currentValue0) << endl;
    if (fabs(currentValueRMS-currentValueRMS0)>delta)  
      cout << currentValueRMS << "!=" << currentValueRMS0 << endl;
    if (fabs(setValue-setValue0)>delta)                
      cout << setValue << "!=" << setValue0 << endl;

 
    currentValue=currentValue0;
    currentValueRMS=currentValueRMS0; 
    setValue=setValue0;

  }

  return ret;
}

int
MaxCamChannel::cleanOldValues(const char *table, unsigned int maxRows) {

#ifdef DM_DAQ

  if (!table) table=GetName();

  MYSQL mysql;
  bool lock = mutex_db !=0;  
  if (lock) mutex_db->Lock(); 
  mysql_init(&mysql);
  bool success = connectToDB(&mysql,"clean old values",GetName()); 
  if(!success)
  {
    if (lock) mutex_db->UnLock(); 
    return 0; 
  }

  TString sql=TString("SELECT count(*) FROM ") + table;
  //cout << sql << endl;
  mysql_real_query(&mysql, (const char*)sql, sql.Length() );
  
  MYSQL_RES *res=mysql_store_result(&mysql);
  MYSQL_ROW row = mysql_fetch_row(res); 
  unsigned long nrow = atoi(row[0]); 
  if (nrow>maxRows) {
    // delete oldest entry
    sql = TString("DELETE FROM ") + table + TString(" ORDER BY timestamp ASC LIMIT "); 
    sql += nrow-maxRows;
    //cout << sql << "  nrow="<<nrow<<"  maxrow="<<maxRows<<endl;
    mysql_real_query(&mysql, (const char*)sql, sql.Length() );
    nrow =  mysql_num_rows(res);
  }
  
  mysql_free_result(res);
  mysql_close(&mysql);

  if (lock) mutex_db->UnLock(); 
  return nrow>maxRows;
#else
  cout << GetName() << ": DB not used " << table << endl;

  return 0;
#endif
}


bool MaxCamChannel::setBusy(bool busy)
{
#ifdef DM_DAQ
  MYSQL mysql;
  bool lock = mutex_db != 0; 
  if (lock) mutex_db->Lock(); 
  mysql_init(&mysql);
  bool success = connectToDB(&mysql,"set busy",GetName());
  if(!success) 
  {
    if (lock) mutex_db->UnLock(); 
    return false; 
  }

  char sql[200]; 
  const char * val = busy ? "1" : "0";
  sprintf(sql,"UPDATE busy SET %s = %s;",GetName(),val); 
  int ret =  mysql_real_query(&mysql,sql,strlen(sql)); 
  cout << "query: " << sql << " returned " << ret << endl;
  mysql_close(&mysql); 
  if (lock) mutex_db->UnLock(); 
#endif
  return false; 
}

bool MaxCamChannel::isBusy()
{
  float busy=0;
  
#ifdef DM_DAQ
  MYSQL mysql;
  bool lock = mutex_db !=0; 
  if (lock) mutex_db->Lock(); 
  mysql_init(&mysql);
  bool success = connectToDB(&mysql,"get busy",GetName());
  if(!success)
  {
    mutex_db->UnLock(); 
    return false; 
  }
  
  char sql[500]; 
  sprintf(sql,"SELECT %s  FROM  busy LIMIT 1",GetName());
  int ret =  mysql_real_query(&mysql,sql,strlen(sql)); 
  cout << "query: " << sql << " returned " << ret << endl;
  MYSQL_RES * res = mysql_store_result(&mysql); 

  if (res == NULL || !mysql_num_rows(res))
    {
      cout << "No busy entry for " << GetName() << " in DB" <<endl; 
    }
  
  else
    {
      MYSQL_ROW row = mysql_fetch_row(res);
      busy = atof(row[0]); 
    }
  mysql_free_result(res);
  mysql_close(&mysql); 
  if (lock) mutex_db->UnLock(); 
#endif
  if(busy==1)
    return true;
  else
    return false; 

}

double MaxCamChannel::getAverage(int n)
{

  double avg=0;

#ifdef DM_DAQ
  
  MYSQL mysql;
  
  bool lock = mutex_db != NULL; 
  if (lock) mutex_db->Lock(); 

  mysql_init(&mysql);
  bool success = connectToDB(&mysql,"get values", GetName()); 
  if (!success) 
  {
    if (lock) mutex_db->UnLock(); 
    return -1; 
  }

  TString sql=TString("SELECT value FROM ") + GetName() + TString(" ORDER BY timestamp DESC LIMIT ");
  sql+=n;
  int ret = mysql_real_query(&mysql, (const char*)sql, sql.Length());
  cout << "query: " << sql << " returned " << ret << endl;

  MYSQL_RES *res=mysql_store_result(&mysql);

  if ((int)mysql_num_rows(res)<n) {
    cout << GetName() << ": Not enough data  " << endl;
    mysql_close(&mysql);
    if (lock) mutex_db->UnLock(); 
    return -1;
  }

  for(int i=0; i<n; i++)
    {
      MYSQL_ROW  row=mysql_fetch_row(res);
      avg+=atof(row[0]);
    }

  avg=avg/double(n);

  mysql_free_result(res);
  mysql_close(&mysql);
  if (lock) mutex_db->UnLock(); 
#endif

  return avg;

}

int MaxCamChannel::rawCAENcmd(const char * cmd, char ** result, float pause , char * port_name)
{

#ifdef DM_DAQ
  /* Get config port if NULL is passed */
  if (port_name == NULL)
  {
    port_name = cfg->getSerial(); 
  }  
 
 /* Initialize options correctly */ 
  bool lock = mutex_caen !=0; 
  if (lock) mutex_caen->Lock(); 
  SerialHandle caen = MaxCamSerial::Open(port_name); 
  int baudRate = 9600; 
  int dataBits = 8;
  int stopBits = 1; 
  int xonxoff=1; 
  int rtscts = 0;
  char parity = 'n';

  size_t l = strlen(cmd); 
  char * cmdbuf = (char*)malloc(l+3); 
  char *res = (char*) calloc(64,1); 

  int ret = MaxCamSerial::Init(caen,baudRate,dataBits,parity,stopBits,xonxoff,rtscts);
  if (!ret) 
  {
    cout << "MaxCamSerial init failed" << endl; 
    goto serial_close;
  }
  MaxCamSerial::Flush(caen); 

  /* add CR to the end  and make non-const copy*/ 
  memcpy(cmdbuf,cmd,l); 
  cmdbuf[l] = '\r'; 
  cmdbuf[l+1] = '\n'; 
  cmdbuf[l+2] = '\0'; 

  ret = MaxCamSerial::Send(caen,cmdbuf,pause); 
  if (!ret) goto serial_close; 
 
  //64 bytes ought to be enough for everyone 
  ret = read(caen,res,63);  
  cout << ret << " bytes read: " << res << endl; 

  if (result == NULL) free(res); 
  else *result = res; 

  serial_close:
    ret+=MaxCamSerial::Close(caen); 

  if (lock) mutex_caen->UnLock(); 
  return ret; 
#else
  return -1; 
#endif
}

MaxCamChannel::~MaxCamChannel()
{
  if (hid.t_int)
  {
    #ifdef DM_MCC
    hid_delete_HIDInterface(&(hid.hid));
    #endif
  }
}
