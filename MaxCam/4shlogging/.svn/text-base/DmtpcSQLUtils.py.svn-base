import MySQLdb as sql

def executeSql(cmd, dbName="DMTPC_TEST"):
    db = openSQL(dbName)
    cur = db.cursor()
    cur.execute(cmd)
    cur.close()
    db.close()

def openSQL(dbName):
    sql_hostname = "localhost"
    sql_username = "dmatter"
    sql_passwd   = "seedark"

    db = sql.connect(host=sql_hostname, user=sql_username, passwd=sql_passwd)
    cur = db.cursor()
    cur.execute("use "+dbName)
    return db

