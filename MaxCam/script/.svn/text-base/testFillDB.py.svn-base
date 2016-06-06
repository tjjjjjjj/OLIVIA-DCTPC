# sample code to stuff values into the database

import MySQLdb as sql
import time

dbName = "DMTPC_TEST"
hostname = "localhost"
user = "dmatter"
passwd = "seedark"

# connect to the MySQL server
db = sql.connect(host=hostname, user=user, passwd=passwd)
cur = db.cursor()

# activate a database
cmd = "USE "+dbName
cur.execute(cmd)

table = "pressure"
cdgbase = 1.0
bpgbase = 3.0
convbase = 5.0
for ii in range(5):
    cdgval  = cdgbase*(ii+1)
    bpgval  = bpgbase*(ii+1)
    convval = convbase*(ii+1)
    cmd = "INSERT INTO " + table + " (value_cdg, value_bpg, value_convectron, timestamp) VALUES( "
    cmd += str(cdgval) + ", "
    cmd += str(bpgval) + ", "
    cmd += str(convval)+ ", "
    cmd += "NOW() )"

    print cmd
    cur.execute(cmd)
    time.sleep(1)

cur.close()
db.close()
