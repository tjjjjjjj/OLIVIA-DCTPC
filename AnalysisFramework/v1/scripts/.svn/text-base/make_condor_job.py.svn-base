#! /usr/bin/env python

#! /usr/bin/env python
import os
import sys
import glob
import re

input1=sys.argv[1]  # first run to process 

workingdir=os.getcwd()+'/'

output=open('job_scripts/condor_spec_'+str(input1), 'w')
output.write("Universe        = vanilla\n")
output.write("Executable = /net/cockroach/data02/jocelyn/data_processing_14April2009/projects/DarkMatter/AnalysisFramework/v1/scripts/runSkimAndSort_Condor\n")
output.write("Arguments       = "+str(input1)+"\n")
output.write("GetEnv          = True\n")
output.write("Initialdir = /net/cockroach/data02/jocelyn/data_processing_14April2009/\n")
output.write("Output          = /net/cockroach/data02/jocelyn/data_processing_14April2009/condor_out/"+str(input1)+".out\n")
output.write("Error           = /net/cockroach/data02/jocelyn/data_processing_14April2009/condor_err/"+str(input1)+".err\n")
output.write("Log             = /net/cockroach/data02/jocelyn/data_processing_14April2009/condor_log/"+str(input1)+".log\n")
output.write("CoreSize        = 0\n")
output.write("Notification        = Never\n")
output.write("Queue\n")
