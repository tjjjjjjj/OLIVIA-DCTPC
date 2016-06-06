#! /usr/bin/env python
import os
import commands
import sys
import datetime

now = datetime.datetime.now()

firstrun=8313
lastrun=8314
camera=100534

sqlcmd="sqlite3 -line mccfg.db 'select * from mccfg where \"CameraProperties.SerialNumber\"=\"0 "+str(camera)+"\" and \"RunProperties.MinimumMatchingDataRunNumber\"="+str(firstrun)+" and \"RunProperties.MaximumMatchingDataRunNumber\"="+str(lastrun)+" order by mccfgdb_timestamp limit 1'"

print sqlcmd
ret=commands.getstatusoutput(sqlcmd)

parameters=ret[1].split('\n')

cameracfgfile=open('cameraProperties0.example','w')
chambercfgfile=open('chamberProperties.example','w')
runcfgfile=open('AlphaRunProperties.example','w')

header=( '#', '# '+str(camera)+' ('+str(firstrun)+'-'+str(lastrun)+')', '# retrieved from sqlite/mccfg.db '+str(now),'#' )

for i in range(len(header)):
    cameracfgfile.write(header[i]+'\n')
    chambercfgfile.write(header[i]+'\n')
    runcfgfile.write(header[i]+'\n')

for i in range(len(parameters)):
    
    par=parameters[i].lstrip(' ')

    if par.startswith('CameraProperties.'):
        if par.split()[len(par.split())-1]!='=':
            cameracfgfile.write(par[len('CameraProperties.'):len(par)]+'\n')
    elif par.startswith('RunProperties.'):
        if par.split()[len(par.split())-1]!='=':
            runcfgfile.write(par[len('RunProperties.'):len(par)]+'\n')

    elif par.startswith('DetectorProperties.'):
        if par.split()[len(par.split())-1]!='=':
            chambercfgfile.write(par[len('DetectorProperties.'):len(par)]+'\n')

    else:
        print par

cameracfgfile.close()
chambercfgfile.close()
runcfgfile.close()
