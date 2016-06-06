#!/usr/bin/python
print "Content-type: text/html"
print

import os
os.environ[ 'HOME' ] = '/usr/local/apache2/html/tmp/' 

import sys
from datetime import datetime
from datetime import timedelta
from numpy import zeros
from numpy import arange
import DB 
import random

import cgi
import cgitb
cgitb.enable()
import MySQLdb

def makePlot(table,StartT=None,EndT=None,Log=False,WithSetVal=False,Value="value",Unit="", Filename=""):
  import matplotlib
  matplotlib.use('Agg')
  from matplotlib import pyplot

  cur = DB.openDB()
  if cur is None:
    print "Could not connect to database!"
    sys.exit(1)
  if WithSetVal is not True:
    cmd = "SELECT "+Value+", timestamp FROM "
    cmd+=table 
    cmd+=" WHERE timestamp > %s && timestamp < %s ORDER BY timestamp ASC"
 
  else:
    cmd = "SELECT "+Value+",timestamp,setval FROM " 
    cmd+=table
    cmd+=" WHERE timestamp > %s && timestamp < %s ORDER BY timestamp ASC" 

  deltaT = timedelta(0,0,0,0,0,5)
  if StartT is None and EndT is None:
    EndT = datetime.utcnow()
    StartT = EndT - deltaT
  elif StartT is None:
    StartT = EndT-deltaT
  elif EndT is None:
    EndT = StartT+deltaT

#strip microseconds from values:
  StartT = StartT - timedelta(0,0,StartT.microsecond)
  EndT = EndT - timedelta(0,0,EndT.microsecond)
  deltaT = StartT-EndT
  totalT = (deltaT.microseconds + 10**6*(deltaT.seconds+3600*24*deltaT.days))/10**6 #in seconds
 
  try:
    num =  cur.execute(cmd,(StartT.isoformat(" "),EndT.isoformat(" ")))
  except MySQLdb.Error, e:
    print "Error %d: %s"%(e.args[0],e.args[1])
    sys.exit(2)

  if num < 2:
    print "gfx/norecent.png"
    return
  vals = zeros(num)
  times = zeros(num)
  if WithSetVal is True:
    set = zeros(num)
  else:
    set = None
  for i in range(num):
    res = cur.fetchone()
    vals[i] = res[0]
    dT = res[1]-StartT
    times[i] = (dT.microseconds + 10**6*(dT.seconds+3600*24*dT.days))/10**6
    if WithSetVal is True:
      set[i]=res[2]
  if (Log is True):
    if WithSetVal is False:
      pyplot.semilogy(times,vals,"x",basey=10)
    else:
      pyplot.semilogy(times,vals,"bx",label="true",basey=10)
      pyplot.semilogy(times,set,"r--",linewidth=2.0,label="set",basey=10)
      pyplot.figlegend((pyplot.subplot(111).lines[0],pyplot.subplot(111).lines[1]),("true","set"),"upper right")
  else:
    if WithSetVal is False:
      pyplot.plot(times,vals,"x")
    else:
      pyplot.plot(times,vals,"bx",label="true")
      pyplot.plot(times,set,"r--",label="set",linewidth=2.0)
      pyplot.figlegend((pyplot.subplot(111).lines[0],pyplot.subplot(111).lines[1]),("true","set"),"upper right")

  pyplot.title(table)
  pyplot.xlabel("Time")
  if Unit=="":
    pyplot.ylabel(Value)
  else:
    pyplot.ylabel(Value+" ["+Unit+"]")
  pyplot.grid(True)
  xtick = (arange(5)+1.)/5*(times[-1]-times[0])+times[0]
  tickname = range(5)
  for i in range(5):
    tdelta = timedelta(0,xtick[i])
    tempTime = StartT+tdelta
    tickname[i] = tempTime.strftime("%Y-%m-%d\n%H:%M:%S")

  pyplot.xticks(xtick,tickname)

  if Filename:
    pyplot.savefig("../html/"+Filename)
  else:
    pyplot.show()
  
  pyplot.clf()
  pyplot.cla()
  print Filename


if __name__=="__main__":
 
  fields = cgi.FieldStorage()
  table = ""
  if "Table" in fields:
    table = fields["Table"].value

  value = ""
  if "Value" in fields:
    value = fields["Value"].value

  startDate = ""
  if "StartDate" in fields:
    startDate = fields["StartDate"].value

  if "StartTime" in fields:
    startTime = fields["StartTime"].value

  endDate = ""
  if "EndDate" in fields:
    endDate = fields["EndDate"].value

  endTime = ""
  if "EndTime" in fields:
    endTime = fields["EndTime"].value

  unit = ""
  if "Unit" in fields:
    unit = fields["Unit"].value

  withSet = False
  if "UseSet" in fields:
    withSet = fields["UseSet"].value == "True"

  logPlot = False
  if "Log" in fields:
    logPlot = fields["Log"].value == "True"

  start = None
  end = None
 
  rndm = random.randint(0,1e9)
  file = "tmp/tmp_plot"+str(rndm)+".png"
  if startDate!="" and startTime!="":
    a=startDate.split("-")
    if (len(a)!=3): 
      a=startDate.split("/")
    if (len(a)!=3): 
      a=startDate.split(":")
    b=startTime.split(":")
    if (len(b)!=3): 
      b=startTime.split("-")
    if (len(b)!=3): 
      b=startTime.split("/")
    start=datetime(int(a[0]),int(a[1]),int(a[2]),int(b[0]),int(b[1]),int(b[2]))

  if endDate!="" and endTime!="":
    a=endDate.split("-")
    if (len(a)!=3): 
      a=endDate.split("/")
    if (len(a)!=3): 
      a=endDate.split(":")
    b=endTime.split(":")
    if (len(b)!=3): 
      b=endTime.split("-")
    if (len(b)!=3): 
      b=endTime.split("/")
    end=datetime(int(a[0]),int(a[1]),int(a[2]),int(b[0]),int(b[1]),int(b[2]))

  makePlot(table,StartT=start,EndT=end,Value=value,Unit=unit,WithSetVal=withSet,Log=logPlot,Filename=file)

