// $Id: monitoring.C,v 1.5 2009/05/07 23:05:34 shawnh Exp $
void monitoring(TString extension="EXT",float startTime,float endTime, Int_t isLocal=0);
int printDatabase( TMySQLServer & , TString table);
int plotPressure( TMySQLServer & , TCanvas & );
int plotCCD( float startTime, float endTime, TMySQLServer & , TCanvas & , TString extension="EXT", Int_t isLocal=0 );
TGraph* plotGeneric( TString variable ,  float startTime, float endTime, TMySQLServer & , TCanvas & , TString drawOptions="alp" , TString yTitle="default", Int_t shouldDraw=1, TString extension="EXT",Int_t isLocal=0);
float ReturnPlottableTime(TDatime &);
int writeTree( float startTime, float endTime, TMySQLServer & , TString extension="EXT", Int_t isLocal=0 );
TTree* writeGeneric(TString variable, float startTime, float endTime, TMySQLServer &mysql);
TTree* writeCCD(float startTime, float endTime, TMySQLServer &mysql);


void monitoring(TString extension, float startTime, float endTime, Int_t isLocal){

  Float_t min_pressure=74; // torr
  Float_t max_pressure=77; // torr
  Float_t min_anode=0.65; // torr
  Float_t max_anode=0.80; // torr
  Float_t min_drift=4.95; // torr
  Float_t max_drift=5.05; // torr
  
  TString directoryPath="/usr/local/apache2/htdocs/tmp/dqm/";
  if(isLocal) directoryPath="";

  // open a canvas plot
  TCanvas *c1=new TCanvas("c1","c1",0,0,800,800);

  // in order for this line to work, have to do some black magic
  // on the LNS machines
  gSystem->Load("libMySQL.so");

  // open a connection to the mysql server
  TMySQLServer  mysql("mysql://mitdm00.mit.edu:/DM_TEST_26544", "dmatter","seedark");
  
  TGraph *temp;

  cout << "Making the pressure DQM plot...<br>" << endl;
  // make the pressure plot
  temp=plotGeneric("pressure", startTime, endTime, mysql,*c1,"alp","Pressure (Torr)",0,extension,isLocal);
  c1->Clear();
  temp->SetMinimum(min_pressure);
  temp->SetMaximum(max_pressure);
  temp->Draw("alp");
  c1->SaveAs(directoryPath+"dqm_"+extension+"/pressure_"+extension+".gif");

  cout << "Making the wire high voltage DQM plot...<br>" << endl;
  // make the anode plot
  temp=plotGeneric("wire_hv", startTime, endTime, mysql,*c1,"alp","Anode (kV)",0,extension,isLocal);
  c1->Clear();
  temp->SetMinimum(min_anode);
  temp->SetMaximum(max_anode);
  temp->Draw("alp");
  c1->SaveAs(directoryPath+"dqm_"+extension+"/anode_"+extension+".gif");

  cout << "Making the mesh high voltage DQM plot...<br>" << endl;
  // make the drift plot
  temp=plotGeneric("mesh_hv", startTime, endTime, mysql,*c1,"alp","Drift (kV)",0,isLocal);
  c1->Clear();
  temp->SetMinimum(min_drift);
  temp->SetMaximum(max_drift);
  temp->Draw("alp");
  c1->SaveAs(directoryPath+"dqm_"+extension+"/drift_"+extension+".gif");
  
  cout << "Making the temp0 DQM plot...<br>" << endl;
  temp=plotGeneric("temp0", startTime, endTime, mysql,*c1,"ap","temp0",1,extension,isLocal);
  // neither of these is meaningful, really
  //  temp=plotGeneric("temp1", startTime, endTime, mysql,*c1,"ap","temp1",1,extension);
  //  temp=plotGeneric("temp2", startTime, endTime, mysql,*c1,"ap","temp2",1,extension);
  
  cout << "Making the CCD DQM plots...<br>" << endl;
  plotCCD(startTime, endTime, mysql,*c1,extension,isLocal);  
  //  printDatabase(mysql,"pressure");

  // now write the information for further use to a root file
  cout << "Writing data to file for this run range...<br>" << endl;
  writeTree(startTime,endTime,mysql,extension,isLocal);
}

TGraph* plotGeneric(TString variable, float startTime, float endTime, TMySQLServer &mysql, TCanvas &c1, TString drawOptions, TString yTitle, Int_t shouldDraw, TString extension,Int_t isLocal){
  
  TString directoryPath="/usr/local/apache2/htdocs/tmp/dqm/";
  if(isLocal) directoryPath="";

  // pull the pressure rows out of the mysql database
  TString queryString="select * from "+variable;
  TMySQLResult *res= mysql.Query(queryString);
  
  if(yTitle=="default") yTitle=variable;
  
  // put the pressure information in a TGraphErrors
  TGraphErrors *variableGraph=new TGraphErrors();
  
  int n=0;
  while ( row=(TMySQLRow*)res->Next() ) {
    
    float value=atof(row->GetField(0));
    float rms=atof(row->GetField(1));
    float setval=atof(row->GetField(2));

    TDatime dt( row->GetField(3) ); 
    float timestamp=ReturnPlottableTime(dt);
    
    if(timestamp>startTime&&timestamp<endTime){
      variableGraph->SetPoint(n,timestamp,value);
      variableGraph->SetPointError(n,0,rms);
      n++;
    }
    
    //    cout << row->GetField(0) << " " << row->GetField(1) << " " << row->GetField(2) << " " << row->GetField(3) << endl;
    //    cout << value << " " << rms << " " << setval << " " << timestamp << endl;
  }
  
  c1.Clear();
  TString T0="%b/%d %I%p";
  variableGraph->GetXaxis()->SetTimeDisplay(1);
  variableGraph->GetXaxis()->SetTimeFormat((const char*)T0);
  variableGraph->GetXaxis()->SetTitle("Time");
  variableGraph->GetYaxis()->SetTitle(yTitle);
  variableGraph->SetLineColor(kRed);
  variableGraph->SetMarkerColor(kRed);
  variableGraph->SetMarkerStyle(2);
  variableGraph->SetMarkerSize(0.5);

  if(shouldDraw){
    variableGraph->Draw(drawOptions);
    c1.SaveAs(directoryPath+"dqm_"+extension+"/value_in_"+variable+"_"+extension+".gif");
  }

  delete res;

  //  delete variableGraph;
  return variableGraph;
} 

int plotCCD(float startTime, float endTime, TMySQLServer &mysql, TCanvas &c1, TString extension, Int_t isLocal){
  
  TString directoryPath="/usr/local/apache2/htdocs/tmp/dqm/";
  if(isLocal) directoryPath="";

  // pull the ccd rows out of the mysql database
  TMySQLResult *res= mysql.Query("select * from ccd");
  
  // put the pressure information in a TGraphErrors
  //  TGraphErrors *ccdTemperatureGraph=new TGraphErrors();
  //  TGraphErrors *ccdAveragePixelYieldGraph=new TGraphErrors();
  //  TGraphErrors *ccdExposureGraph=new TGraphErrors();

  // the ccd average pixel yield for each camera individually
  TGraphErrors *ccdAveragePixelYieldCamera1Graph=new TGraphErrors();
  TGraphErrors *ccdAveragePixelYieldCamera2Graph=new TGraphErrors();

  // the ccd temperature for each camera individually
  TGraphErrors *ccdTemperatureCamera1Graph=new TGraphErrors();
  TGraphErrors *ccdTemperatureCamera2Graph=new TGraphErrors();

  // the ccd exposure for each camera individually
  TGraphErrors *ccdExposureCamera1Graph=new TGraphErrors();
  TGraphErrors *ccdExposureCamera2Graph=new TGraphErrors();

  int n=0;
  int ncamera1=0;
  int ncamera2=0;
  while ( row=(TMySQLRow*)res->Next() ) {

    // CRASHES WHEN WE TRY TO GET THE CCDID VARIABLE (FIELD 5)
    float temperature=atof(row->GetField(0));
    float exposure=atof(row->GetField(1));
    float daqtime=atof(row->GetField(2));
    float avgpixel=atof(row->GetField(3));

    TDatime dt( row->GetField(4) ); 
    float timestamp=ReturnPlottableTime(dt);

    if(timestamp>startTime&&timestamp<endTime){

      // fill the CCD temperature graph
      //      ccdTemperatureGraph->SetPoint(n,timestamp,temperature);
      //      ccdTemperatureGraph->SetPointError(n,0,0);
      //      
      // fill the CCD average pixel yield graph
      //      ccdAveragePixelYieldGraph->SetPoint(n,timestamp,avgpixel);
      //      ccdAveragePixelYieldGraph->SetPointError(n,0,0);
      //
      // fill the CCD exposure graph
      //      ccdExposureGraph->SetPoint(n,timestamp,exposure);
      //      ccdExposureGraph->SetPointError(n,0,0);
      
      // fill the CCD average pixel yield for the cameras individually
      if(n%2==0){
	// for the first camera
	// average pixel yield
	ccdAveragePixelYieldCamera1Graph->SetPoint(ncamera1,timestamp,avgpixel);
	ccdAveragePixelYieldCamera1Graph->SetPointError(ncamera1,0,0);

	// average temperature
	ccdTemperatureCamera1Graph->SetPoint(ncamera1,timestamp,temperature);
	ccdTemperatureCamera1Graph->SetPointError(ncamera1,0,0);

	// exposure
	ccdExposureCamera1Graph->SetPoint(ncamera1,timestamp,exposure);
	ccdExposureCamera1Graph->SetPointError(ncamera1,0,0);

	++ncamera1;
      } else {
	// for the second camera
	// average pixel yield
	ccdAveragePixelYieldCamera2Graph->SetPoint(ncamera2,timestamp,avgpixel);
	ccdAveragePixelYieldCamera2Graph->SetPointError(ncamera2,0,0);
	
	// average temperature
	ccdTemperatureCamera2Graph->SetPoint(ncamera2,timestamp,temperature);
	ccdTemperatureCamera2Graph->SetPointError(ncamera2,0,0);

	// exposure
	ccdExposureCamera2Graph->SetPoint(ncamera2,timestamp,exposure);
	ccdExposureCamera2Graph->SetPointError(ncamera2,0,0);

	++ncamera2;
      }
      
      n++;
    }

    //    cout << row->GetField(0) << " " 
    //	 << row->GetField(1) << " " 
    //	 << row->GetField(2) << " " 
    //	 << row->GetField(3) << " " 
    //	 << row->GetField(4) << " " 
    //	 << endl;
  }

  // make the ccd temperature graph
  c1.Clear();
  TString T0="%b/%d %I%p";
  /*
  ccdTemperatureGraph->GetXaxis()->SetTimeDisplay(1);
  ccdTemperatureGraph->GetXaxis()->SetTimeFormat((const char*)T0);
  ccdTemperatureGraph->GetXaxis()->SetTitle("Time");
  ccdTemperatureGraph->GetYaxis()->SetTitle("ccd temperature");
  ccdTemperatureGraph->GetYaxis()->SetTitleOffset(1.5);
  ccdTemperatureGraph->SetLineColor(kRed);
  ccdTemperatureGraph->SetMarkerColor(kRed);
  ccdTemperatureGraph->SetMarkerStyle(2);
  ccdTemperatureGraph->SetMarkerSize(0.5);
  ccdTemperatureGraph->Draw("ap");
  */
  ccdTemperatureCamera1Graph->GetXaxis()->SetTimeDisplay(1);
  ccdTemperatureCamera1Graph->GetXaxis()->SetTimeFormat((const char*)T0);
  ccdTemperatureCamera1Graph->GetXaxis()->SetTitle("Time");
  ccdTemperatureCamera1Graph->GetYaxis()->SetTitle("ccd temperature");
  ccdTemperatureCamera1Graph->GetYaxis()->SetTitleOffset(1.5);
  ccdTemperatureCamera1Graph->SetLineColor(kRed);
  ccdTemperatureCamera1Graph->SetMarkerColor(kRed);
  ccdTemperatureCamera1Graph->SetMarkerStyle(2);
  ccdTemperatureCamera1Graph->SetMarkerSize(0.5);
  ccdTemperatureCamera1Graph->Draw("ap");  
  ccdTemperatureCamera2Graph->SetLineColor(kBlue);
  ccdTemperatureCamera2Graph->SetMarkerColor(kBlue);
  ccdTemperatureCamera2Graph->SetMarkerStyle(2);
  ccdTemperatureCamera2Graph->SetMarkerSize(0.5);
  ccdTemperatureCamera2Graph->Draw("p");
  c1.SaveAs(directoryPath+"dqm_"+extension+"/temperature_in_ccd_"+extension+".gif");
  
  // make the ccd exposure graph
  c1.Clear();
  /*  ccdExposureGraph->GetXaxis()->SetTimeDisplay(1);
  ccdExposureGraph->GetXaxis()->SetTimeFormat((const char*)T0);
  ccdExposureGraph->GetXaxis()->SetTitle("Time");
  ccdExposureGraph->GetYaxis()->SetTitle("ccd exposure");
  ccdExposureGraph->GetYaxis()->SetTitleOffset(1.5);
  ccdExposureGraph->SetLineColor(kRed);
  ccdExposureGraph->SetMarkerColor(kRed);
  ccdExposureGraph->SetMarkerStyle(2);
  ccdExposureGraph->SetMarkerSize(0.5);
  ccdExposureGraph->Draw("ap");
  */
  ccdExposureCamera1Graph->GetXaxis()->SetTimeDisplay(1);
  ccdExposureCamera1Graph->GetXaxis()->SetTimeFormat((const char*)T0);
  ccdExposureCamera1Graph->GetXaxis()->SetTitle("Time");
  ccdExposureCamera1Graph->GetYaxis()->SetTitle("ccd exposure");
  ccdExposureCamera1Graph->GetYaxis()->SetTitleOffset(1.5);
  ccdExposureCamera1Graph->SetLineColor(kRed);
  ccdExposureCamera1Graph->SetMarkerColor(kRed);
  ccdExposureCamera1Graph->SetMarkerStyle(2);
  ccdExposureCamera1Graph->SetMarkerSize(0.5);
  ccdExposureCamera1Graph->Draw("ap");  
  ccdExposureCamera2Graph->SetLineColor(kBlue);
  ccdExposureCamera2Graph->SetMarkerColor(kBlue);
  ccdExposureCamera2Graph->SetMarkerStyle(2);
  ccdExposureCamera2Graph->SetMarkerSize(0.5);
  ccdExposureCamera2Graph->Draw("p");
  c1.SaveAs(directoryPath+"dqm_"+extension+"/exposure_in_ccd_"+extension+".gif");

  // make the ccd average pixel yield graph
  c1.Clear();
  /* this is for all the pixel yields together
  ccdAveragePixelYieldGraph->GetXaxis()->SetTimeDisplay(1);
  ccdAveragePixelYieldGraph->GetXaxis()->SetTimeFormat((const char*)T0);
  ccdAveragePixelYieldGraph->GetXaxis()->SetTitle("Time");
  ccdAveragePixelYieldGraph->GetYaxis()->SetTitle("ccd average pixel yield");
  ccdAveragePixelYieldGraph->GetYaxis()->SetTitleOffset(1.5);
  ccdAveragePixelYieldGraph->SetLineColor(kRed);
  ccdAveragePixelYieldGraph->SetMarkerColor(kRed);
  ccdAveragePixelYieldGraph->SetMarkerStyle(2);
  ccdAveragePixelYieldGraph->SetMarkerSize(0.5);
  ccdAveragePixelYieldGraph->Draw("ap");
  */
  ccdAveragePixelYieldCamera1Graph->GetXaxis()->SetTimeDisplay(1);
  ccdAveragePixelYieldCamera1Graph->GetXaxis()->SetTimeFormat((const char*)T0);
  ccdAveragePixelYieldCamera1Graph->GetXaxis()->SetTitle("Time");
  ccdAveragePixelYieldCamera1Graph->GetYaxis()->SetTitle("ccd average pixel yield");
  ccdAveragePixelYieldCamera1Graph->GetYaxis()->SetTitleOffset(1.5);
  ccdAveragePixelYieldCamera1Graph->SetLineColor(kRed);
  ccdAveragePixelYieldCamera1Graph->SetMarkerColor(kRed);
  ccdAveragePixelYieldCamera1Graph->SetMarkerStyle(2);
  ccdAveragePixelYieldCamera1Graph->SetMarkerSize(0.5);
  ccdAveragePixelYieldCamera1Graph->SetMinimum(1000);
  ccdAveragePixelYieldCamera1Graph->SetMaximum(1800);
  ccdAveragePixelYieldCamera1Graph->Draw("ap");  
  ccdAveragePixelYieldCamera2Graph->SetLineColor(kBlue);
  ccdAveragePixelYieldCamera2Graph->SetMarkerColor(kBlue);
  ccdAveragePixelYieldCamera2Graph->SetMarkerStyle(2);
  ccdAveragePixelYieldCamera2Graph->SetMarkerSize(0.5);
  ccdAveragePixelYieldCamera2Graph->Draw("p");
  c1.SaveAs(directoryPath+"dqm_"+extension+"/avg_pixel_in_ccd_"+extension+".gif");

  //  delete ccdTemperatureGraph;
  //  delete ccdExposureGraph;
  //  delete ccdAveragePixelYieldGraph;
  delete ccdAveragePixelYieldCamera1Graph;
  delete ccdAveragePixelYieldCamera2Graph;
  delete ccdTemperatureCamera1Graph;
  delete ccdTemperatureCamera2Graph;
  delete ccdExposureCamera1Graph;
  delete ccdExposureCamera2Graph;
  delete res;
  return 1;
}

int printDatabase(TMySQLServer &mysql, TString table){
  printf("Server info: %s\n", mysql.ServerInfo());
  
  TSQLRow *row;
  TSQLResult *res;
  
  // list databases available on server
  printf("\nList all databases on server %s\n", mysql.GetHost());
  res = mysql.GetDataBases();
  while ((row = res->Next())) {
    printf("%s\n", row->GetField(0));
    delete row;
  }
  delete res;
  
  // list tables in database "DM_TEST_26544" (the permission tables)
   printf("\nList all tables in database \"DM_TEST_26544\" on server %s\n",
	  mysql.GetHost());
   res = mysql.GetTables("DM_TEST_26544");
   while ((row = res->Next())) {
     printf("%s\n", row->GetField(0));
     delete row;
   }
   delete res;

   // list columns in table "temp0" in database "mysql"
   printf("\nList all columns in table \"");
   cout << table;
   printf("\" in database \"DM_TEST_26544\" on server %s\n",mysql.GetHost());
   res = mysql.GetColumns("DM_TEST_26544", table);
   while ((row = res->Next())) {
      printf("%s\n", row->GetField(0));
      delete row;
   }
   delete res;

   // exit successfully
   return 1;
}

float ReturnPlottableTime(TDatime &dt){
  float timestamp=dt.Convert();
  return timestamp;
}


int writeTree(float startTime, float endTime, TMySQLServer &mysql, TString extension, Int_t isLocal){

  TString directoryPath="/usr/local/apache2/htdocs/tmp/dqm/";
  if(isLocal) directoryPath="";

  Float_t x[1000];
  TFile *outputTreeFile=new TFile(directoryPath+"dqm_"+extension+"/dqm_"+extension+".root","RECREATE");

  cout << "Saving the pressure information to " << "dqm_"+extension+".root" << "...<br>" << endl;
  TTree *pressureTree=writeGeneric("pressure",startTime,endTime,mysql);
  pressureTree->Write();
  cout << "Saving the wire high voltage information to " << "dqm_"+extension+".root" << "...<br>" << endl;
  TTree *wireHVTree=writeGeneric("wire_hv",startTime,endTime,mysql);
  wireHVTree->Write();
  cout << "Saving the mesh high voltage information to " << "dqm_"+extension+".root" << "...<br>" << endl;
  TTree *meshHVTree=writeGeneric("mesh_hv",startTime,endTime,mysql);
  meshHVTree->Write();
  cout << "Saving the temp0 information to " << "dqm_"+extension+".root" << "...<br>" << endl;
  TTree *temp0Tree=writeGeneric("temp0",startTime,endTime,mysql);
  temp0Tree->Write();
  cout << "Saving the CCD information to " << "dqm_"+extension+".root" << "...<br>" << endl;
  TTree *ccdTree=writeCCD(startTime,endTime,mysql);
  ccdTree->Write();
  outputTreeFile->Close();

  return 1;
}

TTree* writeGeneric(TString variable, float startTime, float endTime, TMySQLServer &mysql){
  
  Float_t x[1000];
  
  variableTree=new TTree(variable,"");
  variableTree->Branch(variable,&x,"val/F:rms:set:time");
  
  // pull the pressure rows out of the mysql database
  TString queryString="select * from "+variable;
  TMySQLResult *res= mysql.Query(queryString);
  
  while ( row=(TMySQLRow*)res->Next() ) {
    x[0]=atof(row->GetField(0));
    x[1]=atof(row->GetField(1));
    x[2]=atof(row->GetField(2));
    
    TDatime dt( row->GetField(3) ); 
    x[3]=ReturnPlottableTime(dt);
    
    // save data from an hour before start of first run 
    // and from an hour after the end of the last run
    if(x[3]>=(startTime-3600)&&x[3]<(endTime+3600)){
      variableTree->Fill();
    }
  }
  
  delete res;
  return variableTree;
} 

TTree* writeCCD(float startTime, float endTime, TMySQLServer &mysql){
  
  Float_t x[1000];

  ccdTree=new TTree("ccd","");
  ccdTree->Branch("ccd",&x,"temperature/F:exposure/F:daqtime/F:avgpixel/F:timestamp/F:ccdid/F");

  // pull the ccd rows out of the mysql database
  TMySQLResult *res= mysql.Query("select * from ccd");

  Int_t n=0;
  while ( row=(TMySQLRow*)res->Next() ) {

    // ccdid
    x[0]=atof(row->GetField(0));
    x[1]=atof(row->GetField(1));
    x[2]=atof(row->GetField(2));
    x[3]=atof(row->GetField(3));
    
    TDatime dt( row->GetField(4) ); 
    x[4]=ReturnPlottableTime(dt);

    if(x[4]>(startTime-3600)&&x[4]<(endTime+3600)){
      
      if(n%2==0){
	x[5]=0;
      } else {
	x[5]=1;
      }
      
      ccdTree->Fill();
      
      n++;
    }
  }

  delete res;
  return ccdTree;
}
