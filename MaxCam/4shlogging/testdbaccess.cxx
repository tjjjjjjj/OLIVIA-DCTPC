#include "TMySQLServer.h"
#include "TSQLResult.h"
#include "TSQLRow.h"

void sqltest() {
  TSQLServer *db = TMySQLServer::Connect("mysql://localhost/DMTPC_TEST", "dmatter", "seedark");
  printf("Server info: %s\n", db->ServerInfo());
  
  TSQLRow *row;
  TSQLResult *res;

  const char *sql = "select value_bpg, timestamp from pressure";
  res = db->Query(sql);
  int nrows = res->GetRowCount();
  printf("Got %d rows in result\n", nrows);

  int nfields = res->GetFieldCount();
  for (int ii=0; ii<nfields; ii++)
    printf("%40s", res->GetFieldName(ii));
  printf("\n");
  for (int ii=0; ii<nfields*40; ii++)
    printf("=");
  printf("\n");


  for (int ii=0; ii<nrows; ii++) {
    row = res->Next();
    for (int jj=0; jj<nfields; jj++) {
      printf("%40s", row->GetField(jj));
    }
    printf("\n");
    delete row;
  }

}
