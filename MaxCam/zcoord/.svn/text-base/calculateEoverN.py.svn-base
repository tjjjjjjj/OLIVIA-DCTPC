#!/usr/local/bin/python

import sys

voltage=float(sys.argv[1])

pressure=80*101325/760. # Torr

# need the number of molecules per cm^3
temperature=300 # Kelvin
kB=1.38e-23 # J/Kelvin

N=pressure/(1.e6*kB*temperature) # the 1e6 is to go to cm, instead of meters

#voltage=23.1 # volts
#voltage=5.e3 
spacing=2*2.54 # cm
#spacing=20

EoverN=voltage/(spacing*N) # V/(cm*(N/cm^3))
print "E/N="+str(EoverN)+" V*cm^2"
print "************or************"
print "E/N="+str(EoverN*1.e17)+" Townsends"

