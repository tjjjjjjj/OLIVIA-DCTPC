float x[100000];
int n;
int day0, time0;
TTree *nt;


void read() {
    //gSystem->Load("libMySQL.so");

  nt=new TTree("nt","");
  nt->Branch("pressure",&x,"val/F:rms:set:time");
  
  
  TMySQLServer  mysql("mysql://mitdm00.mit.edu:/DM_TEST_26544", "dmatter","seedark");
  
  TMySQLResult *res= mysql.Query("select * from pressure");
  
  TMySQLRow *row;
  
  time_t time0;
  time_t event_time;
  
  nt = new TTree("nt","");
  nt->Branch("data",&x, "val/F:rms:set:time");
  
  
  int n=0;
  while ( row=(TMySQLRow*)res->Next() ) {
    
    
    TDatime dt( row->GetField(3) ); 
    event_time = dt.Convert();
    
    if (!n) { time0=event_time; }

    x[3] = event_time - time0;
    
    x[3] /=3600;
  
    for (int i=0; i<3; i++) x[i] = atof(row->GetField(i));
    
    nt->Fill();
    n++;

    
  }
  

}
