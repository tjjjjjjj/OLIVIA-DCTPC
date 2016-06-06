# python script to get the power status of the Synaccess IP Power strip
# via telnetting into the IP Power strip
# No username or password is required

# the command "pshow" requests the status of all outlets
# and receives a response which is a table of port number, port name and port status

# Here is a typical session:

#$ telnet cannon.lns.mit.edu
#Trying 198.125.161.235...
#Connected to cannon.lns.mit.edu.
#Escape character is '^]'.
#Tel
#
#
#************************************************************
#*                                                          *
#*                                                          *
#* Synaccess Networks Inc., Carlsbad, CA, USA. Copyright(c) *
#*                                                          *
#*           System  NPB-20                                 *
#*                                                          *
#*                                                          *
#************************************************************
#
#
#HW:7.0 SW:6.6
#>
#>Type "help" for a list of commands.
#>Make sure to set Telnet mode to Local Echo Off
#>pshow
#
#
#************************************************************
#*                                                          *
#*                                                          *
#*     Power Outlet Port Parameters and Status              *
#*                                                          *
#*                                                          *
#************************************************************
#>
#>
#
#Port |       Name |  status | Reserved By | Timer   | AutoPing
#-----+------------+---------+-------------|---------|----------
#   1 |      Port1 |     Off |       Open  |   Off   |     No
#   2 |      Port2 |     Off |       Open  |   Off   |     No
#   3 |      Port3 |     Off |       Open  |   Off   |     No
#   4 |      Port4 |     Off |       Open  |   Off   |     No
#   5 |      Port5 |     Off |       Open  |   Off   |     No
#   6 |      Port6 |     Off |       Open  |   Off   |     No
#   7 |      Port7 |     Off |       Open  |   Off   |     No
#   8 |      Port8 |     Off |       Open  |   Off   |     No
#   9 |      Port9 |      On |       Open  |   Off   |     No
#  10 |     Port10 |      On |       Open  |   Off   |     No
#  11 |     Port11 |      On |       Open  |   Off   |     No
#  12 |     Port12 |      On |       Open  |   Off   |     No
#  13 |     Port13 |      On |       Open  |   Off   |     No
#  14 |     Port14 |      On |       Open  |   Off   |     No
#  15 |     Port15 |      On |       Open  |   Off   |     No
#  16 |     Port16 |      On |       Open  |   Off   |     No
#
# Total AC current draw:
#  2.86 Amps. Max detected   5.25 Amps.
# Power reboot duration: 5 seconds.
#***Type "pTmShow" to view outlet Timer parameters.
#>                                               



import telnetlib

HOST = "cannon.lns.mit.edu"
tn = telnetlib.Telnet(HOST)

#tn.read_until("Echo Off")
lines = tn.read_until("Echo Off")
outstring = "".join(lines[:])
lines = outstring.split("\n\r")
#for line in lines:
#    print "["+line+"]"
#print " ------------------------ "
#print " sending pshow command"
#print " ------------------------ "
tn.write("pshow\n")
lines = tn.read_until("Total AC current draw:")
outstring = "".join(lines[:])
lines = outstring.split("\n\r")

powerNames = ["Port1",  "Port2",  "Port3",  "Port4",
              "Port5",  "Port6",  "Port7",  "Port8",
              "Port9",  "Port10", "Port11", "Port12",
              "Port13", "Port14", "Port15", "Port16"]
powerDict = {}

for pname in powerNames:
    powerDict[pname] = "unknown"
    for line in lines:
        if line.find(pname) > -1:
            toks = line.split("|")
            status = toks[2].strip()
            powerDict[pname] = status
            break
tn.close()

fout = open("power.tmp", 'w')
for kk in powerNames:
    outstr = kk+"    "+powerDict[kk]+"\n";
    fout.write(outstr)
fout.close()

# couldn't get the urllib2 approach to work...

#import urllib2
#urlToRead = "http://cannon.lns.mit.edu/synOpStatus.shtml"
#urlToRead = "http://dmatter:seedark@cannon.lns.mit.edu/synOpStatus.shtml"
#
#realm  = "Admin"
#uri    = "http://cannon.lns.mit.edu"
#user   = "dmatter"
#passwd = "seedark"
#protocol = 'http'
#proxyServer = "http://cannon.lns.mit.edu:80"
#
#o = urllib2.build_opener()
#f = o.open(urlToRead)

#auth = urllib2.HTTPBasicAuthHandler()
#auth.add_password(realm, uri, user, passwd)
#o = urllib2.build_opener(auth)
#f = o.open(urlToRead)

#phand = urllib2.ProxyHandler( {protocol : proxyServer} )
#pauth = urllib2.HTTPBasicAuthHandler()
#pauth.add_password(realm, uri, user, passwd)
#o = urllib2.build_opener(phand, pauth)
#f = o.open(urlToRead)

#for line in f:
#    print line
