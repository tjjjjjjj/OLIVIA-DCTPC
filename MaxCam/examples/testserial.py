#!/usr/local/bin/python

import serial
import string
import binascii

def calcCRC(cmdStr="test"):

    #  This code was taken directly from a Sierra Instruments 
    #  sample program.  No effort was made to understand it.
    #
    #  Multiplying by 0x0100 just shifts the bits left by two slots
    #  in other words, 'a' = 0x61
    #  and so 'a'*0x0100 = 0x6100
    #
    #  ^ is the xor operator.  1^1=0; 1^0=1; 0^0=0.
    #  3^5 = 6
    #  0xffff^0x6100 = 0x93ff
    #
    #  & is the and operator.  1&1=1; 1&0=0; 0&0=0;

    crc = 0xffff
    print crc
    for ii in range(len(cmdStr)):
        crc=crc^(ord(cmdStr[ii])*0x0100)
        print crc
        for jj in range(8):
            if (crc & 0x8000) == 0x8000:
                crc = ((crc<<1)^0x1021) & 0xffff
            else:
                crc = (crc<<1) & 0xffff

    if ( (crc & 0xff00) == 0x0d00):
        crc += 0x0100
    if ( (crc & 0x00ff) == 0x000d):
        crc += 0x0001
    if ( (crc & 0xff00) == 0x0000):
        crc += 0x0100
    if ( (crc & 0x00ff) == 0x0000):
        crc += 0x0001

    #    print "hex:",hex(crc),"decimal:",crc
    # get the two checksum characters this corresponds to
    #    print hex(crc)[0:2]+hex(crc)[2:4]
    # print hex(crc)[0:2]+hex(crc)[4:6]
    # return the hex string
    crcchar1=hex(crc)[2:4]
    crcchar2=hex(crc)[4:6]
    return binascii.unhexlify(crcchar1+crcchar2)

def makeCommand(cmdStr="test"):
    # get the checksum for this command
    checksum=calcCRC(cmdStr)
    return cmdStr+checksum+'\r'


# end of calcCRC()

cmd=makeCommand('?Ver')

#ser = serial.Serial(0)
#x = ser.read()          # read one byte
#s = ser.read(10)        # read up to ten bytes (timeout)
#line = ser.readline()   # read a '\n' terminated line
#ser.close()
#
#print line
#
#print repr(line)

# open
#opencmd.append('V')
#opencmd.append('l')
#opencmd.append('v')
#opencmd.append('i')
#opencmd.append('1')
#opencmd.append(0x22)
#opencmd.append(0x33)
#opencmd.append(0x0d)  
#opencmd.append(0x00)
#
# open
# DOESN'T LIKE 0 TERMINATION (putting \x00 on the end)
#print ser.write('\x56\x6C\x76\x69\x31\x22\x33'+'\r')
