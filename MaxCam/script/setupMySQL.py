import MySQLdb as sql

db = sql.connect(host="localhost", user="dmatter", passwd="seedark")
cur = db.cursor()

dbName = "DMTPC_TEST"
cmd = "SELECT IF(EXISTS(SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA WHERE SCHEMA_NAME = '" + dbName + "'), 'yes', 'no')"
print cmd
cur.execute(cmd)
response = cur.fetchone()
if response[0] == 'no':
    #then make the database
    print "should make the database"
    #cmd = "CREATE DATABASE " + dbName + ";"
    #cur.execute(cmd)
else:
    print "database already exists."

cmd = "use "+dbName+";"
cur.execute(cmd);

cmd = "create table pressure_CDG (value FLOAT NOT NULL, rms FLOAT DEFAULT 0, setval FLOAT, timestamp TIMESTAMP);"
cur.execute(cmd);

## try opening a database, reading rows from a table
## "SHOW TABLES FROM <DB_NAME>"
## e.g.
## "SHOW TABLES FROM DMTPC_WIPP"
#cmd = "SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA WHERE SCHEMA_NAME='"+dbName+"'"
##cmd = "SELECT TABLE_NAME FROM INFORMATION_SCHEMA.TABLES WHERE INFORMATION_SCHEMA.SCHEMATA='"+dbName+"'"
#cur.execute(cmd)
#response = cur.fetchall()
#print response
#print ''
#
##cmd = "SHOW TABLES FROM "+dbName
#cmd = "SELECT TABLES FROM INFORMATION_SCHEMA.TABLES"
#cur.execute(cmd)
#r = db.store_result()
#print "r.fetch_row() = ", r.fetch_row()
#print ''
#
## get names of columns
#cmd = "SELECT COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS WHERE TABLE_NAME='ccd'"
#cur.execute(cmd)
#response = cur.fetchall()
#print response
#colnames = []
#for resp in response:
#    colnames.append(resp[0])
#print colnames


cur.close()
db.close()
