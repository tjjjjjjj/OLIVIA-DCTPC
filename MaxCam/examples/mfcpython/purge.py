#!/usr/bin/python

import mfc
import serial

ser = serial.Serial('/dev/ttyS2')

# tell controller to take all commands from serial
mfc.setSource(ser, 1)

# flow in units of g/m
mfc.setUnits(ser, 23)

# set the gas
# 8==CF4
# 1==Air
mfc.setGas(ser, 8)

# set to purge 
mfc.setFlowPoint(ser, 0.0)
mfc.setValve(ser, 3)

ser.close()