#!/usr/bin/python

import serial
import string
import binascii
import time

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
    #    print crc
    for ii in range(len(cmdStr)):
        crc=crc^(ord(cmdStr[ii])*0x0100)
        #        print crc
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
    goodhexrep=str("0x%04x" % crc)
    #    print goodhexrep
    #    print hex(crc)
    crcchar1=goodhexrep[2:4]
    crcchar2=goodhexrep[4:6]
    #    print crcchar1
    #    print crcchar2

    return binascii.unhexlify(crcchar1+crcchar2)

def makeCommand(cmdStr="test"):
    # get the checksum for this command
    checksum=calcCRC(cmdStr)
    return cmdStr+checksum+'\r'

def queryMFC(cmdStart,nBytesInAnswer,ser):
    cmd=makeCommand(cmdStart)

    ser.flushInput()
    ser.write(cmd)
    time.sleep(0.3)

    buffer = ser.read(ser.inWaiting())
    readableBuffer=repr(buffer)

    # get rid of the '?'
    if '?' in cmdStart:
        cmdStart=cmdStart.replace('?','')
    
    index=readableBuffer.find(cmdStart)
    if index!=-1:
        response=readableBuffer[(index+len(cmdStart)):(index+len(cmdStart))+nBytesInAnswer]
    else:
        response='no suitable answer'
        
    return response


def getVersion(ser):
    #    print '...in getVersion()'
    version=queryMFC('?Ver',5,ser)
    return version

def getSerial(ser):
    serial=queryMFC('?Srn',6,ser)
    return int(serial)

def getFactor(ser):
    factor=queryMFC('?Ffs',5,ser)
    return float(factor)

def getFlowPoint(ser):
    flowpoint=queryMFC('?Sinv',5,ser)
    return flowpoint

def getSource(ser):
    source=queryMFC('?Sini',1,ser)
    return int(source)

def getPassword(ser):
    password=queryMFC('?Xwrd',4,ser)
    return password

def getValve(ser):
    valve=queryMFC('?Vlvi',1,ser)
    return int(valve)

def getOutput(ser):
    output=queryMFC('?Outi',1,ser)
    return int(output)

def getUnits(ser):
    units=queryMFC('?Unti',1,ser)
    return units

def getMax(ser):
    max=queryMFC('?Fscl',5,ser)
    return float(max)

def getGas(ser):
    gas=queryMFC('?Gasi',1,ser)
    return int(gas)

def getFlow(ser):
    flow=queryMFC('Flow',5,ser)
    return float(flow)

# default to close
def setValve(ser=0, ii=2):
    # VlviValueCCR
    # Value can be an integer between 1 and 3
    # 1=Automatic, controls the set point value
    # 2=Closed
    # 3=Purge
    if int(ii)<1 or int(ii)>3:
        print 'invalid input range, doing nothing'
    else:
        # valid input range!
        cmd=makeCommand('Vlvi'+str(ii))
        print cmd
        ser.flushInput()
        for i in range(10):
            ser.write(cmd)

# default to close
def setFlowPoint(ser=0, fp=0):
    # valid input range!
    cmdStart=''
    if(fp<10):
        cmdStart="Sinv%4.3f" % fp
    if ((fp>=10) and (fp<100)):
        cmdStart="Sinv%4.2f" % fp
    if ((fp>=100) and (fp<1000)):
        cmdStart="Sinv%4.1f" % fp

    cmd=makeCommand(cmdStart)
    ser.flushInput()
    for i in range(10):
        ser.write(cmd)

# default to close
def setSource(ser=0, ii=1):
    # SiniValueCCR
    # Value can be an integer between 1 and 5
    # 1=Module/RS232 as the set point source
    # 2=0-5 V analog set point source
    # 3=0-10 V analog set point source
    # 4=1-5 volt analog set point source
    # 5=4-20 ma analog set point source
    if int(ii)<1 or int(ii)>5:
        print 'invalid input range, doing nothing'
    else:
        # valid input range!
        cmd=makeCommand('Sini'+str(ii))
        print cmd
        ser.flushInput()
        for i in range(10):
            ser.write(cmd)

def setUnits(ser=0, ii=1):
     # UntiValueCCR
     # Value can be an integer between 1 and 3
     # 1=scc/s
     # 2=scc/m
     # 3=scc/H
     # 4=Ncc/s
     # 5=Ncc/m
     # 6=Ncc/H
     # 7=SCF/s
     # 8=SCF/m
     # 9=SCF/
     # 10=NM^3/s
     # 11=NM^3/M
     # 12=NM^3/H
     # 13=SM^3/s
     # 14=SM^3/M
     # 15=SM^3/H
     # 16=sl/s
     # 17=sl/M
     # 18=sl/H
     # 19=NL/s
     # 20=NL/M
     # 21=NL/H
     # 22=g/s
     # 23=g/M
     # 24=g/H
     # 25=kg/s
     # 26=kg/M
     # 27=kg/H
     # 28=lb/s
     # 29=lb/M
     # 30=lb/H
    if int(ii)<1 or int(ii)>30:
        print 'invalid input range, doing nothing'
    else:
        # valid input range!
        cmd=makeCommand('Unti'+str(ii))
        print cmd
        ser.flushInput()
        for i in range(10):
            ser.write(cmd)

def setGas(ser=0, ii=1):
    # GasiValueCCR
    # Value can be an integer between 1 and ?
    # 1=Air
    # 2=Argon
    # 3=Carbon Dioxide
    # 4=Carbon Monoxide
    # 5=Helium
    # 6=90% Argon/10% Carbon Dioxide
    # 7=90% Argon/10% Methane
    # 8=Carbon Tetra-Fluoride
    # 9=Helium 3
    # 10=CS2   
    if int(ii)<1 or int(ii)>10:
        print 'invalid input range, doing nothing'
    else:
        # valid input range!
        cmd=makeCommand('Gasi'+str(ii))
        print cmd
        ser.flushInput()
        for i in range(10):
            ser.write(cmd)

def printBuffer(ser=0):
    buffer = ser.read(ser.inWaiting())
    readableBuffer=repr(buffer)
    print readableBuffer
    

# example of how to use the mfc code
#ser = serial.Serial('/dev/ttyS2')
#
## // tell controller to take all commands from serial
#setSource(ser,1);
##  // set units; 
#setUnits(ser,23); #23==g/min (max flow rate is 0.494 g/min)
##  // set the gas;
#setGas(ser,7); # 7==90% Argon/10% Methane
#
## tell controller to take all commands from serial
#setSource(1)
#
## tell controller to take all commands from serial
#setUnits(23)
#
## set the gas
#setGas(8)
#
## set a flow point
## for air, the maximum is 0.286
## for CF4 the maximum is 0.365
#setFlowPoint(0.)
#
## to start flowing
#setValve(ser, 3)
#    # 1=Automatic, controls the set point value
#    # 2=Closed
#    # 3=Purge
##time.sleep(10)
#
## to close
##setValve(ser, 2)
#
## purge
## I think the flow set point has to be zero, to purge
##setFlowPoint(0.0)
##setValve(ser, 3)
##time.sleep(1.5)
##setValve(ser, 2)
#
##print 'version=',getVersion(ser)
##print 'serial=',getSerial(ser)
##print 'factor=',getFactor(ser)
##print 'flowpoint=',getFlowPoint(ser)
##print 'source=',getSource(ser)
##print 'password=',getPassword(ser)
##print 'valve=',getValve(ser)
##print 'output=',getOutput(ser)
#
## something wrong with this one; sort out
###print 'units=',getUnits(ser)
#
##print 'max=',getMax(ser)
##print 'gas=',getGas(ser)
##print 'flow=',getFlow(ser)
#
#ser.close()
