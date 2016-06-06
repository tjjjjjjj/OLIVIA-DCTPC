#!/usr/bin/python

import mfc
import serial
import time
import sys

print "sys.argv[1] = ", sys.argv[1]
nsec = float(sys.argv[1])

ser = serial.Serial('/dev/ttyS2')

# tell controller to take all commands from serial
mfc.setSource(ser, 1)

# flow in units of g/m
mfc.setUnits(ser, 23)

# set the gas
# 8==CF4
# 1==Air
mfc.setGas(ser, 8)

# purge for n seconds
mfc.setFlowPoint(ser, 0.0)
mfc.setValve(ser, 3)
time.sleep(nsec)
mfc.setValve(ser, 2)

ser.close()
