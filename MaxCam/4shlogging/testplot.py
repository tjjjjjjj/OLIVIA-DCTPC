import sys
import DmtpcLogViewingUtils as dlvu

if len(sys.argv) > 1:
    npoints = int(sys.argv[1])
else:
    npoints = 500
    
print "npoints = ", npoints
dlvu.getBPG(npoints=npoints)
#dlvu.getAll(npoints=npoints)