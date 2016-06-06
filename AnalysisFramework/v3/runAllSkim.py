#! /usr/bin/env python

import os

i = 0


done = []

for g in os.listdir("skim"):
	if g.find("skim") > -1:
		done.append(g[9:14])

	
for f in os.listdir("./data/"):
	
	print("---\nNow processing file "+repr(i))
	i+=1
	n = f[9:14]
	if n in done:
		print "skipping"
		continue 
	cmd1 = "./prep_runbasicfromdir "+ n+" ./data/"
	print cmd1
	os.system(cmd1)
	cmd2 = "./bin/cleanSkim > logs/out"+n+".log"
	print cmd2
	os.system(cmd2)  

