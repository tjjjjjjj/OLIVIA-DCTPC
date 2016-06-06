#!/usr/bin/python
## @file MonitoringPlots.py
# Module for plotting time-varying values of parameters in a MySQL database
#@author Jeremy Lopez
#@date 06 October 2010
from datetime import datetime
from datetime import timedelta
from os import remove
import time
import numpy
import glob

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
import DB #Includes function openDB which returns MySQLdb DictCursor 
import sys
import MySQLdb

## Function which will creates plots for monitoring on the DAQ webpage
def monitoringPlots(N_CAM=2,MAX_NUM=40000):
  
  cur = DB.openDB()
  if cur is None:
    print "Could not connect to database!"
    return 1
  tmp_dir = "/usr/local/apache2/html/tmp/"
  deltaT = timedelta(0,0,0,0,0,5)
  EndT = datetime.utcnow()
  StartT = EndT - deltaT
#strip microseconds from values:
  StartT = StartT - timedelta(0,0,StartT.microsecond)
  EndT = EndT - timedelta(0,0,EndT.microsecond)
  deltaT = StartT-EndT
  totalT = (deltaT.microseconds + 10**6*(deltaT.seconds+3600*24*deltaT.days))/10**6 #in seconds
  stT = StartT.isoformat(" ")
  enT = EndT.isoformat(" ")
  cmd = "(SELECT 'pressure' AS col_name,value, setval,timestamp, rms FROM pressure WHERE timestamp>\'"+stT+"\' && timestamp<\'"+enT+"\' ORDER BY timestamp DESC LIMIT "+str(MAX_NUM)+") UNION" 
  cmd += "(SELECT 'wire_i' AS col_name,value, setval,timestamp, rms FROM wire_i WHERE timestamp>\'"+stT+"\' && timestamp<\'"+enT+"\' ORDER BY timestamp DESC LIMIT "+str(MAX_NUM)+")UNION" 
  cmd += "(SELECT 'wire_hv' AS col_name,value, setval,timestamp, rms FROM wire_hv WHERE timestamp>\'"+stT+"\' && timestamp<\'"+enT+"\' ORDER BY timestamp DESC LIMIT "+str(MAX_NUM)+")UNION" 
  cmd += "(SELECT 'mesh_hv' AS col_name,value, setval,timestamp, rms FROM mesh_hv WHERE timestamp>\'"+stT+"\' && timestamp<\'"+enT+"\' ORDER BY timestamp DESC LIMIT "+str(MAX_NUM)+")UNION" 
  cmd += "(SELECT 'temp0' AS col_name,value, setval,timestamp, rms FROM temp0 WHERE timestamp>\'"+stT+"\' && timestamp<\'"+enT+"\' ORDER BY timestamp DESC LIMIT "+str(MAX_NUM)+")UNION" 
  cmd += "(SELECT 'ccd' AS col_name,avgpixel, temperature,timestamp, ccdid FROM ccd WHERE timestamp>\'"+stT+"\' && timestamp<\'"+enT+"\' ORDER BY timestamp DESC LIMIT "+str(N_CAM*MAX_NUM)+")UNION" 
  cmd += "(SELECT 'scope' AS col_name,scopeid,nwf,timestamp,esum FROM scope WHERE timestamp>\'"+stT+"\' && timestamp<\'"+enT+"\' ORDER BY timestamp DESC LIMIT "+str(MAX_NUM)+")UNION"
  cmd += "(SELECT 'humidity' AS col_name,value,setval,timestamp,rms FROM humidity WHERE timestamp>\'"+stT+"\'&&timestamp<\'"+enT+"\' ORDER BY timestamp DESC LIMIT "+str(MAX_NUM)+")"

  COL_NAME=0
  VALUE=1
  SETVAL=2
  TIMESTAMP=3
  AVGPIXEL=1
  TEMPERATURE=2
  RMS = 4
  CCDID = 4
  NWF = 2
  

  now = datetime.now()
 
  try:
    num =  cur.execute(cmd)
  except MySQLdb.Error, e:
    print "Error %d: %s"%(e.args[0],e.args[1])
    return 2
  cmdtime = datetime.now()-now
  print cmdtime
  now = datetime.now()
  print "Command finished"
  print num," values returned"
  presTime = numpy.empty(MAX_NUM)
  wirehvTime = numpy.empty(MAX_NUM)
  wireiTime = numpy.empty(MAX_NUM)
  meshhvTime = numpy.empty(MAX_NUM)
  tempTime = numpy.empty(MAX_NUM)
  ccd0Time = numpy.empty(MAX_NUM)
  ccd1Time = numpy.empty(MAX_NUM)
  pres = numpy.empty(MAX_NUM)
  wirehv = numpy.empty(MAX_NUM)
  wirei = numpy.empty(MAX_NUM)
  meshhv = numpy.empty(MAX_NUM)
  temp = numpy.empty(MAX_NUM)
  ccd0temp = numpy.empty(MAX_NUM)
  avgpixel0 = numpy.empty(MAX_NUM)
  ccd1temp = numpy.empty(MAX_NUM)
  avgpixel1 = numpy.empty(MAX_NUM)

  humTime = numpy.empty(MAX_NUM)
  humidity = numpy.empty(MAX_NUM)

  scopeTime = numpy.empty(MAX_NUM)
  nWaveform = numpy.empty(MAX_NUM)

  wirehvSet = numpy.empty(MAX_NUM)
  meshhvSet = numpy.empty(MAX_NUM)
  

  i=0
  presNum = 0 
  wirehvNum = 0
  meshhvNum = 0
  wireiNum = 0
  tempNum = 0
  ccd1Num = 0
  ccd0Num = 0
  scopeNum = 0
  humNum = 0

#  for res in result:   
  for i in range(num):
    res = cur.fetchone()
    name = res[COL_NAME]
#Get the time of this row
    dT = res[TIMESTAMP]-StartT
    theTime = (dT.microseconds + 10**6*(dT.seconds+3600*24*dT.days))/10**6
#Check the name and save to the correct array
    if name=="pressure":
      presTime[presNum] = theTime
      pres[presNum] = res[VALUE]
      presNum=presNum+1
    elif name=="wire_hv":
      wirehvTime[wirehvNum] = theTime
      wirehv[wirehvNum] = res[VALUE]
      wirehvSet[wirehvNum] = res[SETVAL]
      wirehvNum = wirehvNum+1
    elif name=="wire_i":
      wireiTime[wireiNum] = theTime
      wirei[wireiNum] =res[VALUE]
      wireiNum = wireiNum+1
    elif name=="mesh_hv":
      meshhvTime[meshhvNum] = theTime
      meshhv[meshhvNum] = res[VALUE]
      meshhvSet[meshhvNum] = res[SETVAL]
      meshhvNum = meshhvNum+1
    elif name=="temp0":
      tempTime[tempNum] = theTime
      temp[tempNum] = res[VALUE]
      tempNum = tempNum+1
    elif name=="ccd":

      if res[CCDID][0:6]=="A80333":
        ccd0Time[ccd0Num] = theTime
        ccd0temp[ccd0Num] = res[TEMPERATURE]
        avgpixel0[ccd0Num] = res[AVGPIXEL]
        ccd0Num = ccd0Num + 1
      elif res[CCDID][0:6]=="A80334":
        ccd1Time[ccd1Num] = theTime
        ccd1temp[ccd1Num] = res[TEMPERATURE]
        avgpixel1[ccd1Num] = res[AVGPIXEL]
        ccd1Num = ccd1Num+1
    elif name=="scope":
      scopeTime[scopeNum] = theTime
      nWaveform[scopeNum] = res[NWF]
      scopeNum = scopeNum + 1
    elif name=="humidity":
      humTime[humNum] = theTime
      humidity[humNum] = res[VALUE]
      humNum = humNum + 1
  looptime = datetime.now()-now
  now=datetime.now()
  print looptime

  plotPres=False
  plotHV=False
  plotWirei=False
  plotTemp=False
  plotCCD=False
  plotScope=False
  plotHum=False
  #timestamp
  ts = str(int(time.time()))
  if presNum>1:
    plotPres=True
    xtick = (numpy.arange(5)+1.)/5*(presTime[0]-presTime[presNum-1])+presTime[presNum-1]
    tickname = range(5)
    for i in range(5):
      tdelta = timedelta(0,xtick[i])
      time_tmp = StartT+tdelta
      tickname[i] = time_tmp.strftime("%Y-%m-%d\n%H:%M:%S")

    if min(pres[0:presNum])*25 < max(pres[0:presNum]):
      pyplot.semilogy(presTime[0:presNum],pres[0:presNum],"b.",basey=10)
    else:
      pyplot.plot(presTime[0:presNum],pres[0:presNum],"b.")
    pyplot.xlabel("Time")
    pyplot.ylabel("Pressure [torr]")
    pyplot.grid(True)
    pyplot.title("Chamber Pressure")
    pyplot.xlim(presTime[presNum-1],presTime[0])
    pyplot.xticks(xtick,tickname)
    pyplot.savefig(tmp_dir+"pressure"+ts+".png")
    pyplot.clf()
    pyplot.cla()    
  else:
    print "Not enough pressure data"

 
  plotAnodeHV = False
  plotMeshHV = False
  mintime = 1e9
  maxtime = -1e9
  if wirehvNum>1:
    wp = pyplot.plot(wirehvTime[0:wirehvNum],wirehv[0:wirehvNum],"b.",label = "Anode V",linewidth=2.0)
    wp1 = pyplot.plot(wirehvTime[0:wirehvNum],wirehvSet[0:wirehvNum],"g--",label = "Anode Set",linewidth=2.0)
    plotAnodeHV = True
    mintime = wirehvTime[wirehvNum-1]
    maxtime = wirehvTime[0]
  else:
    print "Not enough wire hv data"

  if meshhvNum>1:
    mp = pyplot.plot(meshhvTime[0:meshhvNum],meshhv[0:meshhvNum],"r.",label="Drift V",linewidth=2.0)
    mp1 = pyplot.plot(meshhvTime[0:meshhvNum],meshhvSet[0:meshhvNum],"--",color="violet",label = "Drift Set",linewidth=2.0)
    plotMeshHV = True
    mintime = min(mintime,meshhvTime[meshhvNum-1])
    maxtime = max(maxtime,meshhvTime[0])
  else:
    print "Not enough mesh hv data"
  
  if plotAnodeHV or plotMeshHV :
    plotHV=True
    pyplot.xlabel("Time")
    pyplot.ylabel("Voltage [kV]")
    pyplot.title("HV Plots")
    pyplot.grid(True)
    pyplot.xlim(mintime,maxtime)
    xtick = (numpy.arange(5)+1.)/5*(maxtime-mintime)+mintime
    tickname=range(5)
    for i in range(5):
      tdelta = timedelta(0,xtick[i])
      time_tmp = StartT+tdelta
      tickname[i] = time_tmp.strftime("%Y-%m-%d\n%H:%M:%S")
    pyplot.xticks(xtick,tickname)
#    if plotAnodeHV and plotMeshHV:
#      pyplot.figlegend((wp,wp1,mp,mp1),("Anode V","Anode Set","Drift V", "Drift Set"),"upper right")
#    elif plotAnodeHV:
#      pyplot.figlegend((wp,wp1),("Anode V","Anode Set"),"upper right")
#    else:
#      pyplot.figlegend((mp,mp1),("Drift V","Drift Set"),"upper right")
    pyplot.legend(loc="upper left")
    pyplot.savefig(tmp_dir+"hv_plot"+ts+".png")
    pyplot.clf()
    pyplot.cla()
  else:
    print "No hv data!"

  if wireiNum>1:
    plotWirei=True
    xtick = (numpy.arange(5)+1.)/5*(wireiTime[0]-wireiTime[wireiNum-1])+wireiTime[wireiNum-1]
    tickname = range(5)
    for i in range(5):
      tdelta = timedelta(0,xtick[i])
      time_tmp = StartT+tdelta
      tickname[i] = time_tmp.strftime("%Y-%m-%d\n%H:%M:%S")

    pyplot.plot(wireiTime[0:wireiNum],wirei[0:wireiNum],"b.")
    pyplot.xlabel("Time")
    pyplot.ylabel("Current [mA]")
    pyplot.title("Anode Current")
    pyplot.grid(True)
    pyplot.xlim(wireiTime[wireiNum-1],wireiTime[0])
    pyplot.xticks(xtick,tickname)
    pyplot.savefig(tmp_dir+"wire_i"+ts+".png")
    pyplot.clf()
    pyplot.cla()
  else:
    print "Not enough wire i data"


  if tempNum > 1:
    plotTemp = True
    xtick = (numpy.arange(5)+1.)/5*(tempTime[0]-tempTime[tempNum-1])+tempTime[tempNum-1]
    tickname = range(5)
    for i in range(5):
      tdelta = timedelta(0,xtick[i])
      time_tmp = StartT+tdelta
      tickname[i] = time_tmp.strftime("%Y-%m-%d\n%H:%M:%S")
    pyplot.plot(tempTime[0:tempNum],temp[0:tempNum],"b.")
    pyplot.xlabel("Time")
    pyplot.ylabel(r'Temperature [$^\circ$C]')
    pyplot.title("Chamber Temperature")
    pyplot.grid(True)
    pyplot.xlim(tempTime[tempNum-1],tempTime[0])
    pyplot.xticks(xtick,tickname)
    pyplot.savefig(tmp_dir+"temperature"+ts+".png")
    pyplot.clf()
    pyplot.cla()
  else:
    print "Not enough temperature data"

  if ccd0Num>1 or ccd1Num>1:
    plotCCD=True
    minTime = 0
    maxTime = 0
    if ccd0Num<=1:
      minTime = ccd1Time[ccd1Num-1]
      maxTime = ccd1Time[0]
    elif ccd1Num<=1:
      minTime = ccd0Time[ccd0Num-1]
      maxTime = ccd0Time[0]
    else:
      maxTime = min(ccd0Time[0],ccd1Time[0])
      minTime = max(ccd0Time[ccd0Num-1],ccd1Time[ccd1Num-1])
    xtick = (numpy.arange(5)+1.)/5*(maxTime - minTime)+minTime
    tickname = range(5)
    for i in range(5):
      tdelta = timedelta(0,xtick[i])
      time_tmp = StartT+tdelta
      tickname[i] = time_tmp.strftime("%Y-%m-%d\n%H:%M:%S")

    if ccd0Num > 1:
      pyplot.plot(ccd0Time[0:ccd0Num],ccd0temp[0:ccd0Num],"r.",label="Top")
    if ccd1Num > 1:
      pyplot.plot(ccd1Time[0:ccd1Num],ccd1temp[0:ccd1Num],"g.",label="Bottom")
    pyplot.xlabel("Time")
    pyplot.ylabel(r'Temperature [$^\circ$C]')
    pyplot.title("CCD Temperature")
    pyplot.grid(True)
    pyplot.xlim(minTime,maxTime)
    pyplot.xticks(xtick,tickname)
    pyplot.legend(loc="upper left")
    pyplot.savefig(tmp_dir+"ccdtemp"+ts+".png")
    pyplot.clf()
    pyplot.cla()
    if ccd0Num>1:
      pyplot.plot(ccd0Time[0:ccd0Num],avgpixel0[0:ccd0Num],"r.",label="Top")
    if ccd1Num>1:
      pyplot.plot(ccd1Time[0:ccd1Num],avgpixel1[0:ccd1Num],"g.",label="Bottom")
    pyplot.xlabel("Time")
    pyplot.ylabel("Average Pixel")
    pyplot.title("CCD Average Pixel Value")
    pyplot.grid(True)
    pyplot.xlim(minTime,maxTime)
    pyplot.xticks(xtick,tickname)
    pyplot.legend(loc="upper left")
    pyplot.savefig(tmp_dir+"avgpixel"+ts+".png")
    pyplot.clf()
    pyplot.cla()
  else:
    print "Not enough ccd data"

  if humNum > 1:
    plotHum = True
    xtick = (numpy.arange(5)+1.)/5*(humTime[0]-humTime[humNum-1])+humTime[humNum-1]
    tickname = range(5)
    for i in range(5):
      tdelta = timedelta(0,xtick[i])
      time_tmp = StartT+tdelta
      tickname[i] = time_tmp.strftime("%Y-%m-%d\n%H:%M:%S")
    pyplot.plot(humTime[0:humNum],humidity[0:humNum],"b.")
    pyplot.xlabel("Time")
    pyplot.ylabel('% Humidity')
    pyplot.title("Humidity")
    pyplot.grid(True)
    pyplot.xlim(humTime[humNum-1],humTime[0])
    pyplot.xticks(xtick,tickname)
    pyplot.savefig(tmp_dir+"humidity"+ts+".png")
    pyplot.clf()
    pyplot.cla()
  else:
    print "Not enough humidity data"

  if scopeNum > 1:
    plotScope = True
    xtick = (numpy.arange(5)+1.)/5*(scopeTime[0]-scopeTime[scopeNum-1])+scopeTime[scopeNum-1]
    tickname = range(5)
    for i in range(5):
      tdelta = timedelta(0,xtick[i])
      time_tmp = StartT+tdelta
      tickname[i] = time_tmp.strftime("%Y-%m-%d\n%H:%M:%S")
    pyplot.plot(scopeTime[0:scopeNum],nWaveform[0:scopeNum],"b.")
    pyplot.xlabel("Time")
    pyplot.ylabel('Number of Waveforms')
    pyplot.title("Scope Triggers")
    pyplot.grid(True)
    pyplot.xlim(scopeTime[scopeNum-1],scopeTime[0])
    pyplot.xticks(xtick,tickname)
    pyplot.savefig(tmp_dir+"scope_nwf"+ts+".png")
    pyplot.clf()
    pyplot.cla()
  else:
    print "Not enough temperature data"


  #path when accessing from internet
  path=""
  pressurePlot=""
  norecent=path+"gfx/norecent.png"
  if plotPres:
    pressurePlot=path+"tmp/pressure"+ts+".png"
  else:
    pressurePlot=norecent
  avgpixelPlot=""
  ccdtempPlot=""
  if plotCCD:
    avgpixelPlot=path+"tmp/avgpixel"+ts+".png"
    ccdtempPlot=path+"tmp/ccdtemp"+ts+".png"
  else:
    avgpixelPlot=path+"gfx/norecent.png"
    ccdtempPlot=path+"gfx/norecent.png"
  wireiPlot=""
  if plotWirei:
    wireiPlot=path+"tmp/wire_i"+ts+".png"
  else:
    wireiPlot=path+"gfx/norecent.png"
  hvPlot=""
  if plotHV:
    hvPlot=path+"tmp/hv_plot"+ts+".png"
  else:
    hvPlot=path+"gfx/norecent.png"
  tempPlot=""
  if plotTemp:
    tempPlot=path+"tmp/temperature"+ts+".png"
  else:
    tempPlot=path+"gfx/norecent.png"
  if plotHum:
    humPlot=path+"tmp/humidity"+ts+".png"
  else:
    humPlot=path+"gfx/norecent.png"
  if plotTemp:
    nwfPlot=path+"tmp/scope_nwf"+ts+".png"
  else:
    nwfPlot=path+"gfx/norecent.png"
 
  file=open(tmp_dir+"monitoring_plots.html","w")
  file.write("<table border='0' width='100%' cellpadding='5'>\n")
  file.write("<tr>\n")
  file.write("<td width='50%' valign='top'>\n")

  file.write("<img src=\""+pressurePlot+"\" style=\"width:100%;height:100%;\" onClick='load_img_in_popup(\""+pressurePlot+"\");'/>\n")

  file.write("<img src= \""+avgpixelPlot+"\" style=\"width:100%;height:100%;\" onClick='load_img_in_popup(\""+avgpixelPlot+"\");'/>\n")

  file.write("<img src=\""+wireiPlot+"\" style=\"width:100%;height:100%;\" onClick='load_img_in_popup(\""+wireiPlot+"\");'/>\n")

  file.write("<img src=\""+humPlot+"\" style=\"width:100%;height:100%;\" onClick='load_img_in_popup(\""+humPlot+"\");'/>\n")

  file.write("</td>\n")
  file.write("<td width='50%' valign='top'>\n")

  file.write("<img src=\""+hvPlot+"\" style=\"width:100%;height:100%;\" onClick='load_img_in_popup(\""+hvPlot+"\");'/>\n")

  file.write("<img src=\""+nwfPlot+"\" style=\"width:100%;height:100%;\" onClick='load_img_in_popup(\""+nwfPlot+"\");'/>\n")

  file.write("<img src=\""+tempPlot+"\" style=\"width:100%;height:100%;\" onClick='load_img_in_popup(\""+tempPlot+"\");'/>\n")

  file.write("<img src=\""+ccdtempPlot+"\" style=\"width:100%;height:100%;\" onClick='load_img_in_popup(\""+ccdtempPlot+"\");'/>\n")
  file.write("</td>\n")
  file.write("</tr>\n")
  file.write("</table>\n")
  file.close()
 
  #delete older images:
  imgs = glob.glob(tmp_dir+"*.png") 
  for im in imgs:
    if im.rfind(ts)==-1:
      remove(im)

  return 0

## Command line function to create all of the plots.
#  Creates plots and sleeps for 45 sec to keep CPU usage low
def main(*args):
  while(1):
    ret = monitoringPlots()
    if ret==1:
      print "Error occurred! No plots generated."
    time.sleep(30)

if __name__=="__main__":
  main(sys.argv[1:])

