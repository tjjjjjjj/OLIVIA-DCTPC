#!/usr/bin/python

import mfc
import serial

ser = serial.Serial('/dev/ttyS2')

# tell controller to take all commands from serial
mfc.setSource(1)

# flow in units of g/m
mfc.setUnits(23)

# set the gas
# 1==Air
mfc.setGas(1)

mfc.setFlowPoint(0.0)

ser.close()
