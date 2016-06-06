# to run this:
#
# nohup python ~/projects/DarkMatter/MaxCam/logging/DetectorLogger.py > /var/www/html/tmp/logger.nohup.out &
# need to add output redirect, etc, and run from a common location, not from jbattat's directory

##########
# To do
###########
#
# * log data file is not written to after each fill
#   so a viewer would not be able to see recent data...
#   Figure out what triggers a write and do that after each fill
# * currently no timeout for data reads (serial comm, telnet instances, etc)
#   this can block and lead to data loss (I've seen it happen)
# * run serial reads as separate threads?
# * verify automatic log-file roll-over at midnight GMT
# * log-file safety.  if code crashes/restarts, ensure log file is preserved
# * embed the DmtpcLoggableParam within DmtpcLogEntry so that DmtpcLoggableParam doesn't
#   have to be handled/managed by this code
# * communication with database so that DAQ can access certain vars
# * get power status from C++ code, not from external python script...
# * add reality checks on parameter reads. (especially the pressure gauge readings).
# * proper time axes and appropriate labels on TGraphs
# * "AutoSave" option to allow log file to be viewed as it is written
# * need logfile viewing code, both for interactive mode and for automated plot generation
# * could have a single TTimeStamp instance, no?
#   (if the entire program is single-threaded)
# * will time_since_start overflow?  what about XXX_due (e.g. pressure_due)?
#
#
#* Timestamp/value for several parameters
#* Data comes from multiple computers
#  (e.g. CCD temperatures from DAQ, the rest from Slow Control)
#* Permanence
#* Access to data for real-time status and recent-history plots
#  Plotting in ROOT perhaps?
#* Logging at different cadences for each parameter
#* Data load will be something like 50-75 parameters (and timestamps)
#  logged approximately every 20 seconds.
#  That comes to ~78-120 million entries per year
#
#################################
# Brief description of the code #
#################################
#
#  Code to log the environmental parameters of the 4shooter detector
#
#  data is "logged" to the database
#  we could have an ascii log as well for the user to check
#  to see if there are any problems (e.g. error w/ serial comm, etc)
#
#  For information on how to stop the logging program, see the
#  comments in the DetectorLogger() function below.
#
#  Author:  James Battat  jbattat@mit.edu
#  2010 February 4

#  What needs to be logged (differential vs. single-ended)
#  0  ADC Temperature0
#     ADC Temperature1
#     ADC Temperature2
#     ADC Temperature3
#     ADC Humidity
#  5  ADC TiltmeterX
#     ADC TiltmeterY
#     ADC TiltmeterZ
#     ADC +HV (Anode)
#     ADC -HV (Cathode)
# 10  ADC PMT_HV_0
#     ADC PMT_HV_1
#     ADC PMT_HV_2
#     ADC Big gate valve position indicator
#     ADC Small gate valve position indicator
# 15  ADC 
#     ADC 
#     ADC 
#     ADC 
#     ADC 
# 20  ADC 
#     ADC 
#     ADC 
#     ADC 
#     ADC 
# 25  ADC 
#     ADC 
#     ADC 
#     ADC 
#     ADC 
# 30  ADC 
# 31  ADC 
#
#  
#  CCD Parameters --> but how do you get them?
#    CCD0_temperature   
#    CCD1_temperature
#    CCD2_temperature
#    CCD3_temperature
#
#  Serial Varian Turbo controller state (rpm, temperature)
#  Serial Inficon (CDG, BPG)
#  Serial Convectron
#  Serial MFC (is not logged, but is used)
#
#  telnet power1
#  telnet power2
#  telnet power3
#  telnet power4
#  telnet power5
#  telnet power6
#  telnet power7
#  telnet power8  Scroll pump
#  telnet power9  Big gate
#  telnet power10 Small gate
#  telnet power11 Rough valve
#  telnet power12 MFC valve
#  telnet power13
#  telnet power14 
#  telnet power15
#  telnet power16 
#
#

import os.path  # check if file exists
import time
import sys
import DmtpcLoggingUtils as dlu
import DmtpcSQLUtils as dsu

EXIT_SUCCESS = 0

dbName = "DMTPC_TEST"

############################################################################
# perhaps these constants could be moved out to an ASCII configuration file
DMTPC_TPS   = 1000 # ticks per second  (for DMTPC_TPS=1000, 1 tick is 1 ms)
MAX_ENTRIES = 100  # maximum number of log values between fills
SLEEP_SEC   = 0.1  # sleep time in seconds
#
# The Inficon controller (VGC402) reads out two inficon pressure gauges
# The VGC402 has 2 channels.  
# Channel 1 is the BPG (dual gauge, atmosphere to 1e-11 torr)
# Channel 2 is the CDG (capacitance diaphragm gauge, ~1 to 200 torr)
INFICON_SERIAL_PORT = "/dev/ttyS0"
INFICON_BPG_CHAN = 1
INFICON_CDG_CHAN = 2 
############################################################################

verbose = 0

time_prog_started = 0 # start time in ticks 
time_current      = 0 # current epoch time in ticks
time_since_start  = 0 # number of ticks since program started

############################
#  Define logging cadences #
############################
power_interval    = 10*DMTPC_TPS 
pres_interval     = 10*DMTPC_TPS

# keep track of when to log things
# these vals are incremented after each read
power_due  = 0   # read power status of IP power strip
pres_due   = 0   # pressures

#############################################################################
#############################################################################
# set_time(): updates time_since_start & time_current for current execution time
#
# Static In/Out time_current current epoch time, in DMTPC_TPS units
# Out time_since_start  --  time since execution started
#
# All scheduled processing uses 'time_since_start' to decide when to process.
#############################################################################
def set_time():
    global time_current
    global time_since_start
    global time_prog_started
    
    time_current = time.time()*DMTPC_TPS
    time_since_start = time_current - time_prog_started # ticks since execution started
    return time_since_start

#############################################################################
#############################################################################
# begin data logging routines
#############################################################################
def get_pressures():
    if (verbose > 2):
        print "   reading pressures"
    # CDG reading
    presCDG = dlu.inficonReadPressureCDG()
    # BPG reading
    presBPG = dlu.inficonReadPressureBPG()

    presConvectron = dlu.convectronReadPressure()

    # add pressures to db
    print "time, CDG, BPG, Conv = ", time.strftime("%Y-%m-%dT%H:%M%S", time.gmtime()), ", ", presCDG, ", ", presBPG, ", ", presConvectron

    sqlString = "INSERT INTO pressure (value_cdg, value_bpg, value_convectron, timestamp) VALUES('%s', '%s', '%s', UTC_TIMESTAMP())" % (str(presCDG), str(presBPG), str(presConvectron))
    dsu.executeSql(sqlString, dbName=dbName)

def get_powers():
    dlu.synaccessReadPowerStatus()

def do_system_logging():

    global time_since_start
    global pres_due
    global pres_interval
    global power_due
    global power_interval
    
    # handle opening logfile on startup and roll-over of log files
    if (verbose > 1):
        print "do_system_logging()"
        print "time_since_start = ", time_since_start
        print "pres_due = ", pres_due
    if (verbose > 2):
        print "pres_interval = ", pres_interval
        
    if (time_since_start >= pres_due):
        get_pressures()
        pres_due += pres_interval
    #if (time_since_start >= power_due):
    #    get_powers()
    #    power_due += power_interval

def do_waiting():
  if (verbose > 2):
      print "do_waiting()" 
  time.sleep(SLEEP_SEC);  # argument is seconds
  set_time();

#############################################################################
#############################################################################
# init_time():  Initialize the time base and the current time
#
# global out:  time_prog_started, time_current are set to the current time
#############################################################################
def init_time():
    # initialize time
    # time.time() returns seconds since epoch (float)
    global time_prog_started
    global time_current

    time_prog_started = time.time()*DMTPC_TPS  
    time_current = time_prog_started

#############################################################################
#############################################################################
# init_loggable_params():  Initialize the variables to be logged
#############################################################################
def init_loggable_params():
    pass

#############################################################################
#############################################################################
# DetectorLogger(nloop=0):  "infinite" loop for logging of system parameters
#
#  if nloop = 0 then run indefinitely 
#  set nloop to a positive value to execute that many iterations
#  Each iteration takes at least SLEEP_SEC seconds
#
#  By default, the logging code runs in an infinite loop.  To stop
#  the logging code, simply create a file named .stoplog.  For
#  example:
#
#     > touch .stoplog 
#
#  will stop the logger (as long as the file permissions are ok and you've done
#  it in the correct directory).
#############################################################################
def DetectorLogger(nloop=0):

    global time_prog_started
    global time_current
    global time_since_start

    #init_loggable_params()
    init_time()

    print "time_prog_started, time_current = ", time_prog_started, ", ", time_current
    
    ii=0;
    ok_to_log = True
    while (ok_to_log):
        if (verbose > 0):
            print "ii, time_since_start, ok_to_log = %d, %f, %s" % (ii, time_since_start, ok_to_log)
        do_system_logging()
        do_waiting();
    
        ii += 1

        if (verbose > 2):
            print "ii, nloop = ", ii, ", ", nloop
            print "(ii >= nloop) = ", (ii >= nloop)
        if ( (nloop > 0) and (ii >= nloop) ):
        #if ( (ii >= 20 ) ):
            #print "***** debugging mode, just do two iterations *****"
            ok_to_log = False
        if (stop_requested()):
            ok_to_log = False
        if (verbose > 2):
            print "       ok_to_log = ", ok_to_log

    return EXIT_SUCCESS

def stop_requested():
    # check to see if a file called .logstop exists
    stopfile = "/var/www/cgi-bin/.stoplog"
    if os.path.exists(stopfile):
        return True;
    return False;

    
if __name__ == "__main__":
    if len(sys.argv) > 1:
        nloop = int(sys.argv[1])
    else:
        nloop = 0
    print "nloop = ", nloop
    DetectorLogger(nloop)
