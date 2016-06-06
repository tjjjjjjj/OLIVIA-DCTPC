#! /usr/bin/env python
import os

mccfgdb_userid=os.environ['USER']
mccfgdb_firstrun=8313
mccfgdb_lastrun=8314

sqlcmd="""'insert into mccfg (mccfgdb_userid,mccfgdb_firstrun,mccfgdb_lastrun) values (\""""+mccfgdb_userid+"""\","""+str(mccfgdb_firstrun)+""","""+str(mccfgdb_lastrun)+""");'"""
print sqlcmd
#sqlcmd="""'select * from sqlite_master'"""

#os.system("sqlite3 -line mccfg.db "+sqlcmd)

# print
sqlcmdlist="""'select * from mccfg where (mccfgdb_detId=\"4sh\" and """+str(mccfgdb_firstrun)+""">=mccfgdb_firstrun and """+str(mccfgdb_lastrun)+"""<=mccfgdb_lastrun) order by mccfgdb_timestamp asc'"""
print sqlcmdlist
os.system("sqlite3 -line mccfg.db "+sqlcmdlist)
