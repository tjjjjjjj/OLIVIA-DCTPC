#! /usr/bin/env python

import os
import sys

f = open(sys.argv[1])
i=0

for g in f.readlines(): 
    g = g.strip(); 
    print("---\nNow processing file "+repr(i))
    i+=1
    cmd1 = "echo ====="+g+"=====: >> wr2.log"
    cmd15 = "echo " +g + " >> wr2.done"
    print(cmd1)
    os.system(cmd1)
    print (cmd15)
    os.system(cmd15)
    cmd2 = "bin/wr2analysis " +g  +" wr2pass.root wr2counts >> wr2.log 2>>wr2.err" 
    print(cmd2)
    os.system(cmd2)

