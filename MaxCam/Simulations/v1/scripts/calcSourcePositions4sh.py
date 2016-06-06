#! /usr/bin/env python

import sys
import math

lengthcal = 0.179 #mm/px 
coll_length = 1.25 * (25.4)
source_width = 0.1 * (25.4) 

if (len(sys.argv) < 4) : 
  print "calcSourcePositions4sh.py x y track_phi" 
  sys.exit(1)


x = (float(sys.argv[1]) - 512) * lengthcal; 
y = (float(sys.argv[2]) - 512) * lengthcal; 
phi = float(sys.argv[3])  

print "Start Point: " + repr(x) + " , " + repr(y)

source_dist = coll_length - source_width; 
x_prime = x + math.cos(phi) * source_dist
y_prime = y + math.sin(phi) * source_dist


print "Position: "  + repr(x_prime) + " , " + repr(y_prime)
print "Angle: "  + repr(-math.cos(phi)) + " , " + repr(-math.sin(phi))


