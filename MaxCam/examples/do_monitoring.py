#! /usr/bin/env python
import os
import sys
import time

# the current time
timestring = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time()))

# useful on my machine
#export LD_LIBRARY_PATH=/afs/lns.mit.edu/user/shawnh/data01/mysql-5.0.77-linux-i686/lib:$LD_LIBRARY_PATH

debug=0

# run like
# ./do_monitoring.py 816 825
workingdir=os.getcwd()+'/'

startRun=sys.argv[1].split('-')[0]
endRun=sys.argv[1].split('-')[1]
# whether or not we're running on mitdm00 (isLocal=0) or we're trying to
# run this on our machine
isLocal=0

print sys.argv

if len(sys.argv)>2:
    # the user may have specified that this
    # should be done locally
    isLocal=sys.argv[2]

# make the directory for the data quality monitoring plots
outputDirectoryName='/usr/local/apache2/htdocs/tmp/dqm/dqm_runs'+str(startRun)+'-'+str(endRun)
if int(isLocal)>0:
    outputDirectoryName='dqm_runs'+str(startRun)+'-'+str(endRun)
    
outputName='dqm_runs'+str(startRun)+'-'+str(endRun)

if os.path.exists(outputDirectoryName+'/'+outputName+'.txt')==False or debug==1:
    command='mkdir '+outputDirectoryName
    print command
    os.system(command)
    
    # make the data quality monitoring plots
    command='root -l -b -q process_monitoring.C"("'+str(startRun)+'","'+str(endRun)+'","'+str(isLocal)+'")"'
    os.system(command)
    
    # make the data monitoring webpage
    output=open(outputDirectoryName+'/'+'index.html', 'w')
    
    output.write("<html>\n")
    output.write("<body>\n\n")

    output.write("file generated "+timestring+"<br><br>")
    
    if os.path.exists(outputDirectoryName+'/'+outputName+'.txt')==True:
        inputFile=open(outputDirectoryName+'/'+outputName+'.txt','r')
        inputLines=inputFile.readlines()

        # an explanation of the variables in the data quality monitoring plot
        explanationFileName=outputName+'_exp.html'
        
        # make this stuff a link to an explanation of what these variables are
        if os.path.exists(outputDirectoryName+'/'+explanationFileName)==True:
            output.write('<a href="'+explanationFileName+'" style="color:black;text-decoration:none">')

        # write the exposure information that the monitoring plot above should have written into the web page
        for i in range(len(inputLines)):
            if "run=" in inputLines[i]:
                for j in range(10):
                    output.write('&nbsp;')
                
            output.write(inputLines[i])
            output.write('<br>\n')
        
        if os.path.exists(outputDirectoryName+'/'+explanationFileName)==True:
            output.write('</a>')
            
        extension='runs'+str(startRun)+'-'+str(endRun)

        for root,dirs,files, in os.walk(outputDirectoryName):
            for name in files:
                if ".gif" in name and extension in name:
                    modifier=name.replace('.jpeg','')
                    modifier2=modifier
                    output.write("\n")
                    output.write("<img src=\""+name+"\"ALT=\""+name+"\">\n")
                    output.write("<br>")
                    output.write(modifier2+"<br>")
                    output.write("\n")

        # write what versions were used to make this; just the first lines in
        # the scripts
        processMonitoring=open('process_monitoring.C','r')
        monitoring=open('monitoring.C','r')
        
        processMonitoringLines=processMonitoring.readlines()
        monitoringLines=monitoring.readlines()
        
        processMonitoringLines[0].replace("//","")
        monitoringLines[0].replace("//","")

        output.write('<center><h1><font color="purple">.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.</font></h1></center>')
        output.write('These plots and numbers were produced with:<br>')
        output.write(processMonitoringLines[0].replace("//","")+'<br>')
        output.write(monitoringLines[0].replace("//","")+'<br>')

        output.write("\n</html>\n")
        output.write("</body>\n\n")

else:
    print 'DQM plots for this run range already exist!  Don\'t do anything!'







