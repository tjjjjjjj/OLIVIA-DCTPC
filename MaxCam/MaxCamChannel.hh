#ifndef MAXCAM_CHANNEL_HH
#define MAXCAM_CHANNEL_HH

#include "TNamed.h"
#include "MaxCamSerial.hh"
#include "MaxCamSoftCal.hh"
#include "MaxCamChannelConfig.hh"
#include "TMutex.h" 

#ifdef DM_MCC
#include <stdbool.h>
#include <fcntl.h>
#include <asm/types.h>
#include <sys/types.h>
#include <libhid/pmd.h>
#include <libhid/usb-3100.h>
#endif

//Typedef because  cint hates mysql

typedef unsigned int MY_MYSQL;


using namespace std;

/** 
 *
 * MaxCamChannel is a class for reading and setting values from hardware or a slow control database. 
 * It supports many hardware types although not always in a polymorphic fashion. 
 *
 * It can be heavily configured using MaxCamChannelConfig or lightly configured in the constructor. 
 *
 * @see MaxCamChannelConfig
 * @author Cosmin Deaconu <cozzyd@mit.edu>
 */

class MaxCamChannel : public TNamed {

public:

  // Ctors
  /** Creates an en empty MaxCamChannel */ 
  MaxCamChannel() { hid.t_int=0; };


  /** 
   * Creates a MaxCamChannel with the given parameters
   *
   * \attention This is here for legacy or interactive purposes. Normally, one should use a MaxCamChannelConfig. 
   *
   * \param name The ROOT name of the channel
   * \param title The ROOT title of the channel
   * \param idRead the channel id of the read channel or -1 for no read channel
   * \param idWrite the channel id of the write channel or -1 for no write channel
   * \param dbaccess the name of the file containing database access information
   * \param invert true to invert the polarity
   * \param type the MAXCAM_HV_TYPE of the channel. Determines how readADC and writeADC work 
   */
  MaxCamChannel(const char *name, const char *title, int idRead=-1, int idWrite=-1, char *dbaccess=0, bool invert=false, MAXCAM_HV_TYPE_T type = COMEDI);

  /**
   * Creates a MaxCamChannel based on a MaxCamChannelConfig 
   * @sa MaxCamChannelConfig
   * @param cfg The configuration file to use.
   */
  MaxCamChannel(MaxCamChannelConfig * cfg);

  /** 
   * Creates a MaxCamChannel based on a MaxCamChannelConfig name. The name
   * must correspond to a file in the path defined by the environmental variable
   * MCCHANNELCFGS
   * @param the name of the file to parse
   */
  MaxCamChannel(const char * config_name); 

  /** 
   * Copy Constructor
   * \warning Not yet fully implemented
   * \param other another MaxCamChannel to clone
   */
  MaxCamChannel(const MaxCamChannel &other);

  /** Destructor */
  virtual ~MaxCamChannel();

  /** 
   * Reads a pressure value from a Convectron controller
   * \param port_name the serial port to use (e.g. /dev/ttyS0) or NULL to use the value from the config
   * \param verbose if true, more output is dumped
   */
  double readConvectron(char *port_name, bool verbose=false);

  /** 
   Get a pressure reading via serial from an Inficon Vacuum 
   Gauge Controller (VGC).

   Works for models VGC402 or VGC403 (and maybe VGC401, too).
  
   <p>Communication with the VGC40X goes like this:
   <pre>
   Commands are a 3 character ASCII mnemonic, optionally followed by
   associated parameters, then followed by a carriage return or a 
   carriage return and line feed
     e.g. 

      <MNE>[params]<CR>[<LF>]
       where <MNE> stands for mnemonic.
  
       The controller will respond with an acknowledgement:
       e.g.
          <ACK><CR><LF>  or <NAK><CR><LF>
    
     Sending (computer to controller)
      S: <MNE>[params]<CR>[<LF>]
      R: <ACK><CR><LF>
       
    
     Receiving 
       First, send a command to the controller.
       Then look for <ACK> or <NAK>
       If <ACK>, the controller has placed information for retrieval
       in its output buffer, ready for you to retrieve.
       Then ask to retrieve data by sending <ENQ> ("enquire")
       Finally controller will send contents of output buffer to computer.

       e.g.
       S: PR1<CR>
       R: <ACK><CR><LF>
       S: <ENQ>
       R: 0, 7.1000E+02

     <ENQ> = 0x05
     <ACK> = 0x06
     <NAK> = 0x15

     For pressure readings (PR1\r, PR2\r, PR3\r), the return value 
     will look like this:
     A,<PM>B.BBBBE<PM>BB<CR><LF>
     where <PM> = plus or minus
           <CR> = carriage return (0x0D)
           <LF> = line feed (0x0A)
     and
                          A = 0  -- Measurement data ok
                              1  -- Underrange
                              2  -- Overrange
                              3  -- Sensor error
                              4  -- Sensor switched off
                              5  -- No sensor
                              6  -- Identification error
                              7  -- BPG/HPG error
    
     and
        <PM>B.BBBBE<PM>BB  = the pressure reading in the current 
                             unit of measurement

   </pre>

   @param port_name  the serial port e.g. /dev/ttyS0 or NULL to use the value from the config
   @param channel    the VGC402 and VGC403 have 2 and 3 sensor channels,
               respectively.  This specifies which channel to read 
               (1,2, or 3)
   @param baudRate   Serial comm. baud rate.  Defaults to 9600 which is
               the default baud rate for the VGC402.

   @return the read value

 */
 double readInficonController(char *port_name, int channel, long baudRate=9600);

 /** Reads the ACK bye from the infinicon controller (01 I believe) and then two more bytes
  *  (CR LR) in preparation for data transfer
  *  @param rs232 Serial handle to device
  */
 int inficonGetAck(SerialHandle rs232);


  /** 
   * Reads the currentValue from the CAEN or COMEDI device,
   * stores result in currentValue as well as returning it
   * the RMS is stored in currentValueRMS
   *
   * \param nread the number of reads to average (works only on COMEDI devices)
   * \param deviceName the name of the device (/dev/something) or NULL to use from configuration
   *
   * \return the average of the value
   */
  double readADC(int nread=100, char *deviceName=NULL);

  /** 
  * Function reads serial port that receives a stream of bytes with measurement of 
  * pressure from LVG-200TC pressure gauge. Populates currentValue and currentValueRMS. 
  * @param port_name the serial port (e.g. /dev/ttyS0) to use or NULL to use value from config
  * @return the pressure value
  */
  double readSerialLVG200TC(char *port_name=NULL);

  /** 
  * Function reads serial port that receives a stream of bytes with measurement of 
  * pressure from inficon pressure gauge. Populates currentValue and currentValueRMS.
  * @param port_name the serial port (e.g. /dev/ttyS0) to use or NULL to use value from config
  * @return the pressure value
  */
  double readSerialInficon(char *port_name=NULL);


  /**
   * Reads channel value from the slow control database. This is useful when you're taking data 
   * to fill in the current pressure et al.
   * @param table the name of the table to read or NULL to use the channel name. 
   * @param ifield the field number in the table to read. 
   * @param doPrint determines whether or not to print the row or not
   * @param db_handle if not NULL, the _MYSQL connection associated with that pointer is used instead of making a new connection
   * @return the read value of the field
   */
  float  readValueFromDB(const char *table, unsigned int ifield, bool doPrint=false, MY_MYSQL* db_handle = NULL);

  /** Updates the currentValue from the database
   *
   * @param table the name of the table. If NULL, defaults to the channel name. 
   * @param doPrint whether or no to print out the row to stdout
   * @param db_handle if not NULL, the MYSQL connection associated with that pointer is used instead of making a new connection
   **/
  void   readFromDB(const char *table=0, bool doPrint=false, MY_MYSQL* db_handle = NULL);


  /** saves the currentValue, currentValueRMS, setValue and timestamp to 
   * the specified database. 
   * \sa saveValueToDBAndCheck()
   * @param table the database table to save to. If NULL, defaults to the channel name.
   * @param db_handle if not NULL, the _MYSQL connection associated with that pointer is used instead of making a new connection
   * @param save_set if false, the set value is not saved 
   * @return 0 on success, 1 on failure
   */
  int    saveValueToDB(const char *table=0, MY_MYSQL* db_handle = NULL, bool save_set = true);

 /** saves the currentValue, currentValueRMS and setValue to the specified
  * and checks that they were actually written. If they weren't, it tries again indefinitely. 
  * @param table  the dtabase table to save to or NULL to use the channel name
   * @param db_handle if not NULL, the MYSQL connection associated with that pointer is used instead of making a new connection
  * @return the return value of saveValueToDB
  *
  */
  int    saveValueToDBAndCheck(const char *table=0, MY_MYSQL* db_handle = NULL);

  /** Sets the DAC value to setValue. 
   * \warning Note that CAEN voltages should be in  volts and COMEDI voltages should be in kv (typically)
   * @param deviceName the name of the device to use (e.g. /dev/comedi0 for comedi things and /dev/ttyUSB0 for CAEN)
   * @return 0 on success
   */
  int    setDAC(char *deviceName=NULL);

  /** * \deprecated
   * Purges all but maxROWS value from the given table in the slow control database. Note that this should probably no longer be used as using SQL indexes has sped up
   * operation on loarge databases 
   *
   * @param the name of the table (or null to use the channel name)
   * @param maxRows the maximum number of surviving rows
   * returns 0 if successful or could not connect to database. returns 1 if there are still more rows than maxRows
   */
  int    cleanOldValues(const char *table=0, unsigned int maxRows=100);

  /**
   * returns the currentValue
   * @return currentValue
   * @see currentValue
   */
  double val() { return currentValue; }
  /** 
   * The currentValue is populated by most read calls. 
   */
  double currentValue;

  /**
   * returns the currentValueRMS
   * @return currentValueRMS
   * @see currentValueRMS
   */
  double rms() { return currentValueRMS; }
  /** 
   * The currentValueRMS is populated by most read calls. 
   */
  double currentValueRMS; 

  
  /**
   * returns the setValue
   * @return setValue
   * @see setValue
   */
  double set() { return setValue; }    
  /** 
   * The setValue is determines what will be passed to setDAC.
   * This is a public method that is writable from outside the class. 
   */
  double setValue; 

  // set a digital output line
  // val = 0 --> low
  //     = 1 --> high
  int setDio(int val);


  /**
   * ramp sets the DAC to the appropriate value (if it's not alraedy there)
   * the same value is read.
   *
   * There are two ramp algorithms, specifiable in MaxCamChannelConfig (see documentation there for details). 
   *
   *  @param newVal value to set to 
   *  @param checkValue check on read channel to make sure value is set. This parameter only makes a difference for the DIRECT ramp algorithm as ramp needs to always check. 
   *  @param interrupt If nonzero, the bool this points to will be checked periodically and the code will return if interrupted. 
   *  @return 0 on success
   */
  int ramp(double newVal, bool checkValue=true, bool * interrupt = 0);

  /** 
   * ramp according to the setValue in the database for this table. 
   * @param table name to read from (or NULL to use channel name)
   * @return 0 on success 
   */
  int rampFromDB(const char *table=0);

  /** 
   * Toggles inverting the polarity
   * @param invert whether or not to invert
   */
  void setInvertPolarity(bool invert) { invertPolarity=invert; }
  /**
   * Checks if the polarity is inverted. 
   * @return whether or not the polarity is inverted.
   */
  bool isPolarityInverted() { return invertPolarity; }
 
  /** 
   * Print out to stdout the name, currentValue, currentValueRMS and setValue of the channel 
   */ 
  void print();

  /** 
   * Checks if a COMEDI subdevice number corresponds to analog in
   * @param subdev the subdev number of the subdevice to check
   * @param deviceName the device name of the comedi device (usually /dev/comedi0)
   * @return true if the given subdevice is an analog in device
   */
  bool subdevIsAnalogIn(int subdev, char *deviceName); 

  /** 
   * Checks if a COMEDI subdevice number corresponds to analog out
   * @param subdev the subdev number of the subdevice to check
   * @param deviceName the device name of the comedi device (usually /dev/comedi0)
   * @return true if the given subdevice is an analog out device
   */
  bool subdevIsAnalogOut(int subdev, char *deviceName); 
 
  /** 
   * Sets the channel busy status in the slow control database. The busy status is used to indicate that
   * the device is doing something (ramping or refilling or whatever it may be).  
   *
   * For this to work properly, there needs to be a columnn in the busy table with the same name as this channel. 
   *
   * @see isBusy()
   * @param busy whether or not the device should be set to busy
   * @return true on success
   */
  bool setBusy(bool busy); 


  /** Checks if the channel is busy or not in the slow control database.
   *
   * For this to work properly, there needs to be a columnn in the busy table with the same name as this channel. 
   * @return true if the channel is busy. 
   * @see setBusy
   */
  bool isBusy();

  /** 
   * returns the average of the last n entries in the database.
   *
   * Uses the channel name to figure out which table to use. 
   *
   * @param n the number of entries to average
   * @return the average value of the entries
   */
  double getAverage(int n);

  /** 
   * Returns the MaxCamChannelConfig associated with this MaxCamChannel 
   */
  MaxCamChannelConfig * getConfig() { return cfg; } 

  /**
   * Perform a raw CAEN command cmd. The result is stored in the 
   pointer pointed to by result unless result is NULL; 
   <p>CR is auto appended to cmd
   <p>In general, the command should have the form (although this
   code does not sanity check): 
   <pre>$BD:**,CMD:***,CH*,PAR:***,VAL:***.**</pre>
   \param cmd command to execute
   \param result pointer to pointer to store result. Auto allocated. If NULL not stored. 
   \param pause seconds to pause after sending command
   \param port name or NULL to get from config
  
   \return Return status of 0 indicates success. 

  */
  int rawCAENcmd(const char * cmd, char ** result, float pause = 0.1, char * port_name = NULL); 

  void synchronizeComedi(TMutex * m) {mutex_comedi = m; } 
  void synchronizeDB(TMutex * m) {mutex_db = m; } 
  void synchronizeUSB(TMutex * m) {mutex_usb = m; } 
  void synchronizeCAEN(TMutex * m) {mutex_caen = m; } 

  unsigned * getMCCPtr() { return hid.t_int;}; 
  void setMCCPtr(unsigned * ptr) { hid.t_int = ptr;} 


private:

  int readID; ///< channel id for reading
  int writeID; ///< channel id for writing
  bool invertPolarity; //</ invert polarity

  /**
   * saves the offset guesses into the database for the channel. The database
   * needs to have an offset_guesses table. 
   *
   * This is used by the follow ramp method to get an initial value for
   * the write offset so that ramping can occur faster. 
   *
   * @returns 0 on success
   */
  int saveOffsets();  

  /** 
   * Initializes the device calibrations and dbaccess for the channel. Also loads 
   * offset guesses if the ramp algorithm is follow.
   */
  void init(); 

  /** The config file for this 
   * channel
   * Not persisted. 
   */ 
  MaxCamChannelConfig * cfg; //!

  /** The calibrator for this channel.
   * Not persisted.
   */
  MaxCamSoftCal * calibrator; //!

  bool offset_guess_init; ///< true if the offset guess has been initialized
  double offset_guess; ///< the current offset gguess

  union 
  {
    unsigned * t_int;
#ifdef DM_MCC
    HIDInterface * hid;
#endif
  } hid; //!

  TMutex *mutex_comedi; //!
  TMutex *mutex_db; //!
  TMutex *mutex_usb; //!
  TMutex *mutex_caen; //!

  ClassDef(MaxCamChannel,6)
};

#endif

