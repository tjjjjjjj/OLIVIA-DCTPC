##@file DB.py
# This module is designed to contain functions allowing us to manipulate a MySQL database
#@author Jeremy Lopez
#@date 06 October 2010
import MySQLdb

##Function that opens a database and returns a cursor if successful
# @param hostname The database hostname
# @param userid The database username
# @param password The user's password
# @param database The name of the database
def openDB(hostname="pippin.lns.mit.edu",userid="dmatter",password="seedark",database="DM_SLOWCONTROL"):
  try: 
    db = MySQLdb.connect(host=hostname,user=userid,passwd=password,db=database)
  except MySQLdb.Error,e:
    print "Error %d: %s"%(e.args[0],e.args[1])
    return None
  cur = db.cursor()

  return cur
  

