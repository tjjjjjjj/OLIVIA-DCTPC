import sys, time
import serial

INFICON_SERIAL_PORT = 0 # /dev/ttyS0
INFICON_BAUDRATE = 9600
INFICON_TIMEOUT = 1.0 # seconds

BPG_CHAN = 1  # combo bayerd-alpert pirani gauge (atm to ~1e-10 torr)
CDG_CHAN = 2  # capacitance diaphragm gauge (0-200 torr, accurate)

verbose = 0

def getInficonReading(channel, baudrate=9600, timeout_sec=1.0):
    ser = serial.Serial()
    ser.port = INFICON_SERIAL_PORT
    ser.timeout = timeout
    ser.baudrate = baudrate

    ser.open()
    #print "ser.isOpen = ", ser.isOpen()
    if not ser.isOpen():
        print "error, serial port did not open"
        print ser
        print " returning ...."
        return "-1"

    # send the request for the pressure reading
    if ( (channel != BPG_CHAN) and (channel != CDG_CHAN) ):
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
    #print "response = ["+response+"]"

    #ACK_STR = r"'\x06'"
    ACK_HEX = '\x06'
    NAK_HEX = '\x15'
    ENQ_HEX = '\x05'
    
    for char in response:
        #print "char = ", char, "repr(char) = ["+repr(char)+"]"
        #if repr(char) == ACK_STR:
        #    print "found ACK"
        if char == ACK_HEX:
            #print "found ACK_HEX"
            break
        elif char == NAK_HEX:
            print "ERROR:  found NAK"
    #print "end"

    bytesSent = ser.write(ENQ_HEX)
    #print "bytesSent = ", bytesSent

    response = ser.readline(None, '\n')

    # wrap up serial comm
    ser.close()

    #print "response = ["+response+"]"
    response = response.strip()
    if (verbose > 0):
        print "response = ["+response+"]"
    
    return response

def processInficonReading(response):
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
    #print "status_str = ["+status_str+"]"
    #print "value_str  = ["+value_str+"]"
    #
    value  = float(value_str)
    status = int(status_str)
    #print "status = ", status
    #print "value = ", value

    if status != 0:
        return -10*status
    else:
        return value

def getPressureBPG():
    return getPressure(BPG_CHAN)

def getPressureCDG():
    return getPressure(CDG_CHAN)
    
def getPressure(channel):
    response = getInficonReading(channel, baudrate=INFICON_BAUDRATE,
                                 timeout_sec=INFICON_TIMEOUT)
    return processInficonReading(response)

# if there is an error, return a pressure reading that is -10 times
# the error number

if __name__ == "__main__":
    
    baudrate = 9600  # 8 N 1 is the default
    timeout = 1 # 1 second timeout

    for ii in range(10):
        print "getPressureBPG() = ", getPressureBPG()
        print "getPressureCDG() = ", getPressureCDG()
        time.sleep(0.1)
