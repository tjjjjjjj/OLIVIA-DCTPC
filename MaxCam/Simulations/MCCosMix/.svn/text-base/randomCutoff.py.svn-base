#! /usr/bin/env python

import optparse
import MySQLdb
import os
import sys
import string
from random import choice

if __name__ == "__main__":

    parser=optparse.OptionParser()
    parser.add_option('-k','--keyword',help='Run Keyword',dest='keyword',\
                          action='store',type='string',default="")    
    parser.add_option('-l','--location',help='Location',dest='location',\
                          action='store',type='string',default="WIPP")
    parser.add_option('-e','--exposure',help='Exposure in ms',dest='exposure',\
                          action='store',type='float',default=1000)
    parser.add_option('-a','--anodemin',help='Anode Minumum Voltage (kV)',\
                      dest='anoden',action='store',type='float',default=-5.0)
    parser.add_option('-A','--anodemax',help='Anode Maxiumum Voltage (kV)',\
                      dest='anodex',action='store',type='float',default=8)
    parser.add_option('-i','--currentmin',help='Anode Minumum Current (mA)',\
                      dest='ain',action='store',type='float',default=-0.001)
    parser.add_option('-I','--currentmax',help='Anode Maxiumum Current (mA)',\
                      dest='aix',action='store',type='float',default=1)
    parser.add_option('-c','--cathodemin',help='Cathode Minumum Voltage (kV)',\
                      dest='cathoden',action='store',type='float',default=-5.0)
    parser.add_option('-C','--cathodemax',help='Cathode Maximum Voltage (kV)',\
                      dest='cathodex',action='store',type='float',default=10)
    parser.add_option('-p','--pressuremin',help='Minimum Pressure (torr)',\
                      dest='pressuren',action='store',type='float',default=0.0)
    parser.add_option('-P','--pressuremax',help='Maximum Pressure (torr)',\
                      dest='pressurex',action='store',type='float',default=900)
    parser.add_option('-s','--scopetracesmin',help='Minimum Scope Traces Avg',\
                      dest='scopen',action='store',type='float',default=0.0)
    parser.add_option('-S','--scopetracesmax',help='Maximum Scope Traces Avg',\
                      dest='scopex',action='store',type='float',default=25.0)
    parser.add_option('-d','--datemin',\
                      help='Start Date/Time, format \"YYYY-MM-DD HH:MM\"',\
                      dest='timen',action='store',type='string',\
                      default="2010-01-01 0:00")
    parser.add_option('-D','--datemax',\
                      help='End Date/Time, format \"YYYY-MM-DD HH:MM\"',\
                      dest='timex',action='store',type='string',\
                      default="2013-01-01 0:00")
    parser.add_option('-t','--tempmin',help='Minumum Temperature (deg C)',\
                      dest='tempn',action='store',type='float',default=20)
    parser.add_option('-T','--tempmax',help='Maximum Temperature (deg C)',\
                      dest='tempx',action='store',type='float',default=40)
    parser.add_option('-q','--sql',action="store_true",dest="sql",\
                      default=False,help="print the sql command; useful for debug")
    parser.add_option('-b','--blank',help='Fraction of blank frames',\
                      dest='blank',action='store',type='float',default=0.0)
    parser.add_option('-n','--inputpath',help='Path to input data file',\
                      dest='ipath',action='store',type='string',\
                      default="/net/zwicky/dmtpc/data/10L/")
    parser.add_option('-o','--outputpath',help='Path to output MC/data file',\
                      dest='opath',action='store',type='string',\
                      default="../v1/output/data/")


    (opts,args)=parser.parse_args()
    
    conn = MySQLdb.connect (host = "pippin.lns.mit.edu",
                            user = "dmatter",
                            passwd = "seedark",
                            db = "DM_SLOWCONTROL")

    cursor = conn.cursor()

    cmd = "SELECT file FROM run_desc WHERE (keyword=\""+opts.keyword+\
          "\") AND (location=\""+opts.location+"\") AND (exposure="+\
          str(opts.exposure)+") AND (anode_hv_avg>"+str(opts.anoden)+\
          ") AND (anode_hv_avg<"+str(opts.anodex)+") AND (drift_hv_avg>"+\
          str(opts.cathoden)+") AND (drift_hv_avg<"+str(opts.cathodex)+\
          ") AND (timestamp>\""+opts.timen+"\") AND (timestamp<\""+opts.timex+\
          "\") AND (temp0_avg>"+str(opts.tempn)+") AND (temp0_avg<"+\
          str(opts.tempx)+") AND (anode_i_avg>"+str(opts.ain)+\
          ") AND (anode_i_avg<"+str(opts.aix)+") AND (pressure_avg>"+\
          str(opts.pressuren)+") AND (pressure_avg<"+str(opts.pressurex)+\
          ") "#AND (ntriggers_avg >"+str(opts.scopen) +") AND (ntriggers_avg<"+\
#          str(opts.scopex)+")"

    if(opts.sql):
        print cmd
    
    cursor.execute(cmd)

    result = cursor.fetchall()

    #make a list of all the results

    filelist=[]
    for row in result:
        filelist.append(row[0])

    file = choice(filelist)
    file = opts.ipath + file

    if(len(args) != 1):
        print "Wrong number of arguments! Quitting"
        sys.exit(0)

    cmd = "./MCCutoff "+args[0]+" "+file+" "+str(opts.blank)+" "+opts.opath
    print cmd
    os.system(cmd)
