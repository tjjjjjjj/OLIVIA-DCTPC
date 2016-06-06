import pylab
import sys

filename = sys.argv[1]
# open the file
#filename = "pressurefile.dat"

fh = open(filename, "r")

datetimes   = []
inficons    = []
convectrons = []
for line in fh:
    if line.startswith("time") or line.startswith("#"):
        continue
    datetime, inficon, convectron = line.split()
    datetimes.append(datetime)
    inficons.append(float(inficon))
    convectrons.append(float(convectron))

fh.close()
pylab.semilogy(inficons, 'k.')
pylab.ylabel("Pressure (torr)")
pylab.xlabel("sample number")
pylab.savefig(filename+".pdf")

