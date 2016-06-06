# a collection of routines that are useful for communicating 
# with DMTPC hardware and for logging data

import sys
import serial
import telnetlib

# serial connections on rubin.lns.mit.edu
# /dev/ttyS0 Inficon_VGC402
# /dev/ttyS1 Convectron_375
# /dev/ttyS2 MFC
# /dev/ttyS3 free
# /dev/ttyS4 software_not_enabled

def convectronReadPressure(port=1, baudrate=19200, bytesize=serial.EIGHTBITS, parity=serial.PARITY_NONE, stopbits=serial.STOPBITS_ONE, timeout_sec=1.0, xonxoff=False, rtscts=False, verbose=0):
    # port     - serial port number.  e.g. port=1 means /dev/ttyS1 (in linux)
    # baudrate - the serial baud rate
    requestPressure = "RD\r"

    ser = serial.Serial()
    ser.port = port
    ser.timeout = timeout_sec
    ser.baudrate = baudrate
    ser.bytesize = bytesize
    ser.parity = parity
    ser.stopbits = stopbits
    ser.xonxoff = xonxoff
    ser.rtscts = rtscts
    
    ser.open()
    if not ser.isOpen():
        print "error, serial port did not open"
        print ser
        print " returning ...."
        sys.exit()

    # send the request for the pressure reading
    ser.write(requestPressure)
    response = ser.readline(None, '\r')
    # wrap up serial comm
    ser.close()
    if (verbose > 1):
        print "response = ["+response+"]"
    # strip the carriage return
    response = response.strip('\r')
    if (verbose > 1):
        print "response = ["+response+"]"

    if response.find("OPN SNSR") > -1:
        print "error:  defective transducer"
        val = -1
    elif response.find("SNSR UNP") > -1:
        print "error:  sensor is unplugged"
        val = -2
    elif response.find("SNSR OVP") > -1:
        print "error:  pressure is higher than 999 torr"
        val = -3
    else:
        val = float(response)
        
    if (verbose):
        print "val = ", val
    return val

def convectronReadVersion(port=1, baudrate=19200, bytesize=serial.EIGHTBITS, parity=serial.PARITY_NONE, stopbits=serial.STOPBITS_ONE, timeout_sec=1.0, xonxoff=False, rtscts=False, verbose=0):
    # get the version number of the controller

    ser = serial.Serial()
    ser.port = port
    ser.timeout = timeout_sec
    ser.baudrate = baudrate
    ser.bytesize = bytesize
    ser.parity = parity
    ser.stopbits = stopbits
    ser.xonxoff = xonxoff
    ser.rtscts = rtscts
    
    ser.open()
    if not ser.isOpen():
        print "error, serial port did not open"
        print ser
        print " returning ...."
        sys.exit()

    ser.write("VER\r")
    response = ser.readline(None, '\r')
    if (verbose > 1):
        print "response = ["+response+"]"
    response = response.strip("\r")
    if (verbose):
        print "response = ["+response+"]"
    return response



INFICON_BPG_CHAN = 1  # combo bayerd-alpert pirani gauge (atm to ~1e-10 torr)
INFICON_CDG_CHAN = 2  # capacitance diaphragm gauge (0-200 torr, accurate)

def inficonGetPressureReading(channel, port=0, baudrate=9600, timeout_sec=1.0, verbose=0):
    ser = serial.Serial()
    ser.port = port
    ser.timeout = timeout_sec
    ser.baudrate = baudrate

    ser.open()
    if not ser.isOpen():
        print "error, serial port did not open"
        print ser
        print " returning ...."
        return "-1"

    # send the request for the pressure reading
    if ( (channel != INFICON_BPG_CHAN) and (channel != INFICON_CDG_CHAN) ):
        print "error, invalid channel: ", channel
        print "exiting..."
        sys.exit()
    sendString = "PR"
    sendString += str(channel)
    sendString += "\r"
    if (verbose > 1):
        print sendString

    ser.write(sendString)
    response = ser.readline(None, '\n')
    if (verbose > 0):
        print "response = ["+response+"]"

    #ACK_STR = r"'\x06'"
    ACK_HEX = '\x06'
    NAK_HEX = '\x15'
    ENQ_HEX = '\x05'
    
    for char in response:
        if (verbose > 0):
            print "char = ", char, "repr(char) = ["+repr(char)+"]"
            if repr(char) == ACK_STR:
                print "found ACK"
        if char == ACK_HEX:
            if (verbose > 0):
                print "found ACK_HEX"
            break
        elif char == NAK_HEX:
            print "ERROR:  found NAK"

    bytesSent = ser.write(ENQ_HEX)
    if (verbose > 0):
        print "bytesSent = ", bytesSent

    response = ser.readline(None, '\n')

    # wrap up serial comm
    ser.close()

    response = response.strip()
    if (verbose > 0):
        print "response = ["+response+"]"
    
    return response

def inficonProcessPressureReading(response, verbose=0):
    # expects a response like:
    #  a,<PM>b.bbbbE<PM>bb<CR><LF>
    # a = Channel status
    #     0 = measurement data OK
    #     1 = underrange
    #     2 = overrange
    #     3 = Sensor error
    #     4 = Sensor switched off
    #     5 = No sensor
    #     6 = identification error
    #     7 = BPG/HPG error
    # b.bbbbE<PM>bb is the pressure reading
    # <PM> = plus or minus (a single character)
    # <CR> = carriage return
    # <LF> = line feed

    status_str, value_str =  response.split(",")
    if (verbose > 0):
        print "status_str = ["+status_str+"]"
        print "value_str  = ["+value_str+"]"
    value  = float(value_str)
    status = int(status_str)
    if (verbose > 0):
        print "status = ", status
        print "value = ", value

    if status != 0:
        return -10*status
    else:
        return value

def inficonReadPressureBPG():
    return inficonReadPressure(INFICON_BPG_CHAN)

def inficonReadPressureCDG():
    return inficonReadPressure(INFICON_CDG_CHAN)
    
def inficonReadPressure(channel):
    response = inficonGetPressureReading(channel, baudrate=9600,
                                         timeout_sec=1.0)
    return inficonProcessPressureReading(response)

def synaccessReadPowerStatus(host="cannon.lns.mit.edu", verbose=0):
    # python script to get the power status of the Synaccess IP Power strip
    # via telnetting into the IP Power strip
    # No username or password is required
    
    # the command "pshow" requests the status of all outlets
    # and receives a response which is a table of port number, port name and port status
    
    # Here is a typical session:
    
    #$ telnet cannon.lns.mit.edu
    #Trying 198.125.161.235...
    #Connected to cannon.lns.mit.edu.
    #Escape character is '^]'.
    #Tel
    #
    #
    #************************************************************
    #*                                                          *
    #*                                                          *
    #* Synaccess Networks Inc., Carlsbad, CA, USA. Copyright(c) *
    #*                                                          *
    #*           System  NPB-20                                 *
    #*                                                          *
    #*                                                          *
    #************************************************************
    #
    #
    #HW:7.0 SW:6.6
    #>
    #>Type "help" for a list of commands.
    #>Make sure to set Telnet mode to Local Echo Off
    #>pshow
    #
    #
    #************************************************************
    #*                                                          *
    #*                                                          *
    #*     Power Outlet Port Parameters and Status              *
    #*                                                          *
    #*                                                          *
    #************************************************************
    #>
    #>
    #
    #Port |       Name |  status | Reserved By | Timer   | AutoPing
    #-----+------------+---------+-------------|---------|----------
    #   1 |      Port1 |     Off |       Open  |   Off   |     No
    #   2 |      Port2 |     Off |       Open  |   Off   |     No
    #   3 |      Port3 |     Off |       Open  |   Off   |     No
    #   4 |      Port4 |     Off |       Open  |   Off   |     No
    #   5 |      Port5 |     Off |       Open  |   Off   |     No
    #   6 |      Port6 |     Off |       Open  |   Off   |     No
    #   7 |      Port7 |     Off |       Open  |   Off   |     No
    #   8 |      Port8 |     Off |       Open  |   Off   |     No
    #   9 |      Port9 |      On |       Open  |   Off   |     No
    #  10 |     Port10 |      On |       Open  |   Off   |     No
    #  11 |     Port11 |      On |       Open  |   Off   |     No
    #  12 |     Port12 |      On |       Open  |   Off   |     No
    #  13 |     Port13 |      On |       Open  |   Off   |     No
    #  14 |     Port14 |      On |       Open  |   Off   |     No
    #  15 |     Port15 |      On |       Open  |   Off   |     No
    #  16 |     Port16 |      On |       Open  |   Off   |     No
    #
    # Total AC current draw:
    #  2.86 Amps. Max detected   5.25 Amps.
    # Power reboot duration: 5 seconds.
    #***Type "pTmShow" to view outlet Timer parameters.
    #>                                               

    tn = telnetlib.Telnet(host)

    lines = tn.read_until("Echo Off")
    outstring = "".join(lines[:])
    lines = outstring.split("\n\r")
    if (verbose > 0):
        for line in lines:
            print "["+line+"]"
        print " ------------------------ "
        print " sending pshow command"
        print " ------------------------ "
    tn.write("pshow\n")
    lines = tn.read_until("Total AC current draw:")
    tn.close()

    outstring = "".join(lines[:])
    lines = outstring.split("\n\r")

    powerNames = ["Port1",  "Port2",  "Port3",  "Port4",
                  "Port5",  "Port6",  "Port7",  "ScrollPump",
                  "BigGate",  "SmallGate", "RoughValve", "MFCValve",
                  "Port13", "Port14", "Port15", "Port16"]
    powerDict = {}

    for pname in powerNames:
        powerDict[pname] = "unknown"
        for line in lines:
            if line.find(pname) > -1:
                toks = line.split("|")
                status = toks[2].strip()
                powerDict[pname] = status
                break

    if (verbose > 0):        
        for kk in powerNames:
            print "kk, vv:  ", kk, ", ", powerDict[kk]

    # couldn't get the urllib2 approach to work...
    
    #import urllib2
    #urlToRead = "http://cannon.lns.mit.edu/synOpStatus.shtml"
    #urlToRead = "http://dmatter:seedark@cannon.lns.mit.edu/synOpStatus.shtml"
    #
    #realm  = "Admin"
    #uri    = "http://cannon.lns.mit.edu"
    #user   = "dmatter"
    #passwd = "seedark"
    #protocol = 'http'
    #proxyServer = "http://cannon.lns.mit.edu:80"
    #
    #o = urllib2.build_opener()
    #f = o.open(urlToRead)
    
    #auth = urllib2.HTTPBasicAuthHandler()
    #auth.add_password(realm, uri, user, passwd)
    #o = urllib2.build_opener(auth)
    #f = o.open(urlToRead)
    
    #phand = urllib2.ProxyHandler( {protocol : proxyServer} )
    #pauth = urllib2.HTTPBasicAuthHandler()
    #pauth.add_password(realm, uri, user, passwd)
    #o = urllib2.build_opener(phand, pauth)
    #f = o.open(urlToRead)
    
    #for line in f:
    #    print line


