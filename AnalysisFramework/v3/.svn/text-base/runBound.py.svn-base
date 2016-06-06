#! /usr/bin/env python

import os
import sys


if len(sys.argv) < 4:
	print ""
	print "ERROR:  must provide three args"
	print "> python "+sys.argv[0]+" RUNMIN RUNMAX DIR"
	print "where:"
	print "        RUNMIN = number of lowest run"
	print "        RUNMAX = number of highest run"
	print "        DIR    = directory holding run files"
	print "Assumes filenames are:  dmtpc_runXXXXX.root"
	sys.exit()
	
runmin = int(sys.argv[1])
runmax = int(sys.argv[2])
dir    = sys.argv[3]

print "runmin, max = ", runmin, ", ", runmax
print "dir = ",dir

for ii in range(runmax-runmin+1):
	id = ii+runmin
	print("---\nNow processing file "+repr(id))

	runnum_str = "%05d" % (id)
	cmd1 = "./prep_runbasicfromdir "+ runnum_str+" "+dir
	print cmd1
	os.system(cmd1)
	cmd2 = "./bin/cleanSkim > logs/out_"+runnum_str+".log"
	print cmd2
	os.system(cmd2)  

