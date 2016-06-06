import time
import datetime
import pylab, matplotlib.dates
import DmtpcSQLUtils as dsu
import DmtpcPlotUtils as dpu

dbName = "DMTPC_TEST"

def getAll(npoints=0):

    cmd = "SELECT value_bpg, value_cdg, value_convectron, timestamp FROM pressure ORDER BY timestamp DESC"
    if npoints > 0:
        cmd += " limit "+str(npoints)
    db = dsu.openSQL(dbName)
    cur = db.cursor()
    cur.execute(cmd)
    resp = cur.fetchall()
    cur.close()
    db.close()

    bpgs  = []
    cdgs  = []
    convs = []
    dates = []

    for rr in resp:
        bpgs.append(float(rr[0]))
        cdgs.append(float(rr[1]))
        convs.append(float(rr[2]))
        dates.append(rr[3])

    figure = pylab.figure(1)
    subplot = pylab.subplot(111)
    ax = figure.gca()

    timeRange = dates[-1]-dates[0]
    year_range = timeRange.days/365.
    #lines = ax.semilogy(dates, bpgs, 'k.')
    #lines = ax.semilogy(dates, cdgs, 'r.')
    lines = ax.plot(dates, bpgs, 'k.')
    lines = ax.plot(dates, cdgs, 'r.')
    lines = ax.plot(dates, convs, 'g.')
    dpu.format_line_ticks(ax, year_range)
    pylab.show()
    
# get a list of BPG pressure data and timestamps
# and plot in pylab
def getBPG(npoints=0, autoupdate=False):

    #cmd = "SELECT value_bpg, timestamp FROM pressure"
    cmd = "SELECT value_bpg, timestamp FROM pressure ORDER BY timestamp DESC"
    if npoints > 0:
        cmd += " limit "+str(npoints)

    db = dsu.openSQL(dbName)
    cur = db.cursor()
    cur.execute(cmd)
    resp = cur.fetchall()
    cur.close()
    db.close()

    vals  = []
    dates = []

    # interactive mode on
    if (autoupdate):
        pylab.ion()
    
    for rr in resp:
        vals.append(float(rr[0]))
        dates.append(rr[1])

    figure = pylab.figure(1)
    subplot = pylab.subplot(111)
    ax  = figure.gca()

    timeRange  = dates[-1]-dates[0]
    year_range = timeRange.days/365.
    #ax.plot_date(pylab.date2num(dates), vals, fmt='k.')
    #figure.autofmt_xdate()
    lines = ax.semilogy(dates, vals, 'k.')
    dpu.format_line_ticks(ax, year_range)

    if (autoupdate):
        figure.canvas.draw()
    else:
        pylab.show()

    if (autoupdate):
        print "autoupdate"
        time.sleep(5.0)
        #getBPG(npoints=npoints, autoupdate=autoupdate)
        
    #pylab.show()
    #pylab.savefig("junk.png")
