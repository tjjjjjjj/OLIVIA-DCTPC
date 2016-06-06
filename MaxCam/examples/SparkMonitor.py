# ======================================
# DMTPC Spark Monitoring 
#
# James Battat
# jbattat@mit.edu
# 2011 April 8
# =============================
#
#   Monitor, count, and timestamp transient events.
#   Designed to monitor, sense and log with timestamp 
#   spark events by listening to the current monitor port 
#   on a HV supply
#   
#   You can stop the loop at anytime by making a file called .stopenv
#   in this directory.  
#     > touch .stopenv
#   The file will automatically be deleted the next time this program runs
#
#   Run with -h or --help for list of options
#   
# =============================
# Bertan Hardware Description
# =============================
#   The Bertan current monitor ports puts out 10 V per milliamp
#
# =============================
# U3 Hardware Description 
# =============================
# 
#   Highlights:
#   * U3 counters trigger on falling edges
#   * Counters/Timers must be on pin FIO4 or higher
#   * Voltage protection on the EIO lines is better than the FIO lines
#     (see below for details)
# 
#   16 I/O lines (8 FIO and 8 EIO) can be individually configured
#     as digital input, digital output or analog input
#     <= 2 of these lines can be configured as Timers (Timer0, Timer1)
#     <= 2 of these lines can be configured as Counters (Counter0, Counter1)
#
#   8 Flexible I/O (FIO) lines:
#     FIO0, FIO1, FIO2, ..., FIO7
#
#   8 Extended I/O (EIO) lines:
#     EIO0, EIO1, EIO2, ..., EIO7
#  
#   Voltage protection on I/O lines
#     The FIO lines can withstand continuous voltages of +/- 10V
#     The EIO/CIO lines can withstand continuous voltages of +/- 6V
#
#   Digital "Low"  input = -0.3 to 0.8 V
#   Digital "High" input =  2.0 to 5.8 V
#   Output Low Voltage (No Load) 0 V
#     FIO      Sinking  1 mA  0.55 V
#     EIO/CIO  Sinking  1 mA  0.18 V
#     EIO/CIO  Sinking  5 mA  0.9  V
#   Output High Voltage (No Load) 3.3 V
#     FIO      Sourcing 1 mA  2.75 V
#     EIO/CIO  Sourcing 1 mA  3.12 V
#     EIO/CIO  Sourcing 5 mA  2.4  V
#   Counter Input Frequency (Hardware v1.21+) 8 MHz	
#   Input Timer Total Edge Rate (v1.21+) 
#     No Stream       30000 edges/s	
#     While Streaming  7000 edges/s
#
#   Digital I/O lines use 3.3 V logic and are 5 Volt tolerant.
#
#   DB15 Pinout
#     1  Vs     9  CIO0	
#     2  CIO1  10  CIO2	
#     3  CIO3  11  GND	
#     4  EIO0  12  EIO1	
#     5  EIO2  13  EIO3	
#     6  EIO4  14  EIO5	
#     7  EIO6  15  EIO7	
#     8  GND		
#
#   Counters
#     32 bit --   (2**32)-1 = 4,294,967,295
#     2 counters (Counter0, Counter1)
#     Each counter consists of a 32-bit register that accumulates 
#     the number of falling edges detected on the external pin. 
#
#   Setting up the timers and counters
#     Ordering of timers and counters is:
#       Timer0, Timer1, Counter0, Counter1
#     Where the first item is TimerCounterPinOffset lines above FIO0
#     and the rest are sequential above the first.
#     FIO0..7 are the first 8 lines and EIO0..7 are the next 8.
#
#     TimerCounterPinOffset must be >= 4.
#
#     For example, if you wanted to have Counter0 be on EIO0, and
#     you did not enable any Timers, then you would want 
#     TimerCounterPinOffset=8
# 
#     If, however, you also enabled a timer, then 
#     TimerCounterPinOffset=8 would put Timer0 on EIO0
#     and Counter0 on EIO1.
# 
# =============================
# Polling the U3
# =============================
#   Two modes:  Command/response & Streaming
#   Command/response reads take 0.6-4.0 ms
#   We're using Command/response with a sleep
#
#
# =============================
# To do:
# =============================
#   * Write to SQL database (see log_spark())
#   * Can this code also work for the cathode?
#     b/c of its non-zero dc current level
#   * Deal with a wrapped counter
#     if ctr val < old_ctr_val then compute the number of ticks
#     taking into account the counter wrap
#     if it's indeed a 32-bit counter, then this would be
#     (2**32)-1 = 4,294,967,295
#     At a spark per second, this would be 136 years
#     so likely not a critical upgrade...
#

import sys
import os
import time
import u3
from optparse import OptionParser

def init_daq(debug=False):
    d = u3.U3()
    d.debug = debug

    #----------------------------------
    # Configure the counter channel
    #----------------------------------
    # TimerCounterPinOffset specifies the number of lines
    #            above FIO0 for the counter.  MUST BE in [4,8].
    #            Assuming that no timers are in play, then  
    #            4 --> Counter0 on FIO4
    #            5 --> Counter0 on FIO5
    #            6 --> Counter0 on FIO6
    #            7 --> Counter0 on FIO7
    #            8 --> Counter0 on EIO0
    #            9 --> Counter0 on EIO0
    # FIOAnalog  Each bit sets whether that FIO line is analog input (=1)
    #            or digital I/O (=0)
    # EIOAnalog  Same as fioAnalog, but for the EIO lines
    pinOffset = 8  
    fioAnalog = 0
    eioAnalog = 0
    d.configIO(EnableCounter0=True, FIOAnalog=fioAnalog, EIOAnalog=eioAnalog, 
               TimerCounterPinOffset=pinOffset)
    return d

# =============================================================
# ---------- Define a class to hold info on the spark ---------
# ---------- and ensure that class attributes are frozen ------
# ---------- i.e. that attributes cannot be added on the fly --
def frozen(set):
    """Raise an error when trying to set an undeclared name."""
    def set_attr(self,name,value):
        if hasattr(self,name):
            set(self,name,value) 
        else:
            raise AttributeError("You cannot add attributes to %s" % self)
    return set_attr
    
class Frozen(object):
    """Subclasses of Frozen are frozen, i.e. it is impossibile to add
    new attributes to them and their instances."""
    __setattr__=frozen(object.__setattr__)
    class __metaclass__(type):
        __setattr__=frozen(type.__setattr__)
        
class SparkInfo(Frozen):
    timestamp = "" # [string] timestamp in sql format "YYYY-MM-DD HH:MM:SS"
    secfrac   = 0  # [integer] sub-second tagging of spark time (milliseconds)
    counter   = 0  # [integer] value on counter for last spark
    utcsec    = 0  # [float] return value of time.time()
    def __init__(self, timestamp, secfrac, counter, utcsec):
        self.timestamp = timestamp
        self.secfrac   = secfrac
        self.counter   = counter
        self.utcsec    = utcsec
# ------------ end of Spark class definition and freezing -----
# =============================================================

def log_spark(sparkinfo):
    millisec_str =  "%03d" % sparkinfo.secfrac
    counter_str  = "%010d" % sparkinfo.counter
    utcsec_str   = "%.6f" % sparkinfo.utcsec
    print sparkinfo.timestamp+"."+millisec_str+"  "+counter_str+"  "+utcsec_str
    sys.stdout.flush()
    # add code here to connect and write to the db
    #

def handle_spark(ctr_val):
    # determine the time of the spark
    utcsec_since_epoch = time.time()  # returns a float, with fractional seconds

    # get the fractional second, down to 1 ms
    millisec = int(round(1000*(utcsec_since_epoch-int(utcsec_since_epoch))))

    # Generate MySQL format timestamp
    sqltime = time.strftime("%Y-%m-%d %H:%M:%S",time.gmtime(utcsec_since_epoch))
    # and a data nugget to log
    sparkinfo = SparkInfo(sqltime, millisec, ctr_val, utcsec_since_epoch)

    # write to the database
    log_spark(sparkinfo)

def check_for_spark(daq, old_ctr_val):
    out = daq.getFeedback(u3.Counter(counter = 0, Reset=False))
    ctr_val = out[0]

    if (old_ctr_val < ctr_val):  # spark occurred
        handle_spark(ctr_val)
    elif (old_ctr_val > ctr_val): # counter wrap occurred
        print "WARNING:  looks like the counter wrapped"
        print "ctr_val, old_ctr_val = ", \
            ctr_val, ", ", old_ctr_val
        print "but this is highly unlikely..."

    # no matter what happens, update the old_ctr_val
    # most of the time, ctr_val == old_ctr_val already...
    return ctr_val
    
def do_waiting(sleep_seconds):
    time.sleep(sleep_seconds)

def stop_requested(filename):
    if (os.path.exists(filename)):
        return True
    return False

def clear_stop_request(filename, verbose=False):
    progdir = "."
    dirlist = os.listdir(progdir)
    srf_exists = False
    for ff in dirlist:
        if (verbose):
            print "  found file: ", ff
        if ff == filename:
            print "found stop request file:  ", filename
            srf_exists = True
            break

    if srf_exists:  # try to remove it
        print "Trying to remove it"
        os.remove(filename)
        
    return
    
def parse_options():
    usage = "usage: %prog [options] type"
    parser = OptionParser(usage=usage)
    parser.add_option("-n", "--nloop", action="store", type="string", 
                      dest="nloop", metavar="NLOOP",
                      help="Number of iterations (0 = infinite)") 
    parser.add_option("-s", "--sleep", action="store", type="string", 
                      dest="sleep",  metavar="SLEEP",
                      help="Duration of sleep (decimal seconds)") 
    parser.add_option("-v", "--verbose", action="store_true", 
                      dest="verbose",  metavar="VERBOSE",
                      help="provide verbose output") 
    parser.add_option("-d", "--debug_daq", action="store_true", 
                      dest="debug",  metavar="VERBOSE",
                      help="provide debug output for DAQ (LabJack U3)") 

    parser.set_defaults(nloop="0",sleep="0.1", verbose=False)

    (options, args) = parser.parse_args()

    return (options, args)


def SparkMonitor(options):

    # User control inputs 
    # (see init_parser() for parameter descriptions)
    nloop     = int(options.nloop)    
    sleep_sec = float(options.sleep)  
    verbose   = options.verbose
    daq_debug = options.debug

    # Start fresh
    iloop             = 0L
    old_counter_val   = 0L
    counter_wrap_val  = 0L # Counter wrap handling is not yet implemented...
    ok_to_loop        = True
    stop_request_file = ".stopsparkmon"
    clear_stop_request(stop_request_file, verbose)

    # Initialize the device
    d = init_daq(daq_debug)

    # Poll the counter in a loop
    while(ok_to_loop):
        old_counter_val = check_for_spark(d, old_counter_val)
        do_waiting(sleep_sec)
        iloop = iloop+1

        # Should we quit?
        if stop_requested(stop_request_file):
            ok_to_loop = False
        if (nloop > 0) and (iloop >= nloop):
            ok_to_loop = False


if __name__ == "__main__":

    (options, args) = parse_options()

    if options.verbose:
        print "options      = ", options

    # Start looking for sparks
    SparkMonitor(options)
