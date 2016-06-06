import sys
import serial

requestPressure = "RD\r"
port = 1 # /dev/ttyS1
baudrate = 19200  # 8 N 1 is the default
timeout = 1 # 1 second timeout

ser = serial.Serial()
ser.port = port
ser.timeout = timeout
ser.baudrate = baudrate

ser.open()
if not ser.isOpen():
    print "error, serial port did not open"
    print ser
    print " returning ...."
    sys.exit()

# send the request for the pressure reading
ser.write(requestPressure)
#nbytes = 9 # 8 characters and a carriage return
#response = ser.read(nbytes)
response = ser.readline(None, '\r')
print "response = ["+response+"]"
# strip the carriage return
response = response.strip('\r')
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

print "val = ", val


# get the version number of the controller
ser.write("VER\r")
response = ser.readline(None, '\r')
print "response = ["+response+"]"
response = response.strip("\r")
print "response = ["+response+"]"

# is NIST calibration valid?
ser.write("CA\r")
response = ser.readline(None, '\r')
print "response = ["+response+"]"
response = response.strip("\r")
print "response = ["+response+"]"

# wrap up
ser.close()
