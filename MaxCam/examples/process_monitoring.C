// $Id: process_monitoring.C,v 1.5 2009/05/21 13:51:24 shawnh Exp $

#include "monitoring.C"
#include <fstream>
#include <iomanip>

void process_monitoring(Int_t minRun=800, Int_t maxRun=801, Int_t isLocal=0);
void addRun(int i);
void GetExposures(TString extension="EXT", Int_t, Int_t, TString, Int_t isLocal=0);
void getRunSparkRate(Double_t runStartTime, Double_t runEndTime, Double_t [], Double_t []);
void getRunSparkRateFromFile(TString extension="EXT",Double_t runStartTime, Double_t runEndTime, Double_t [], Double_t []);
void postMonitoring(TString extension="EXT", Int_t minRun, Int_t maxRun);
void DressGraph(TGraph &, TString yTitle, Int_t color=kRed);
TH1F* myCreateYieldHisto(TH1* h);

void process_monitoring(Int_t minRun, Int_t maxRun, Int_t isLocal)
{
  cout << "<br>" << endl;

  // make sure the plot labels aren't too large 
  // due to the default sizing of the text labels
  // in the MaxCam rootlogon.C file
  gStyle->SetLabelSize(0.02,"x");
  gStyle->SetTitleSize(0.03,"x");
  gStyle->SetLabelSize(0.02,"y");
  gStyle->SetTitleSize(0.03,"y");
  gStyle->SetLabelSize(0.02,"z");
  gStyle->SetTitleSize(0.03,"z");

  if(isLocal>0) gROOT->ProcessLine(".L ../DmtpcDataset.hh");
  else gROOT->ProcessLine(".L DmtpcDataset.hh");

  DmtpcDataset *d=new DmtpcDataset();

  TString minRunName=makeRunName(minRun);
  TString maxRunName=makeRunName(maxRun);
  
  TString extension="runs";
  extension+=minRun;
  extension+='-';
  extension+=maxRun;
  
  // get the first dataset and the last dataset into the
  // chain
  d->openRootFile(minRunName);
  d->chain()->Add(maxRunName);
  
  // get the first time
  d->getEvent(0);
  float firstTime=ReturnPlottableTime((*d->event()->timeStamp()));
  d->getEvent(d->chain()->GetEntries()-1);
  float lastTime=ReturnPlottableTime((*d->event()->timeStamp()));

  cout << "<br>" << endl;

  if(isLocal<2){
    // make the monitoring plots
    cout << "making all monitoring plots...<br>" << endl;
    monitoring(extension,firstTime,lastTime,isLocal);
    
    // now we want to produce the exposures and etc.
    cout << "Getting all exposure information...<br>" << endl;
    GetExposures(extension,minRun,maxRun,extension,isLocal);
  } else {
    // we want to make monitoring plots after the fact;
    // so probably filling in holes.  We'll take the info
    // directly from the dmtpc data files

    // make the monitoring plots
    cout << "making all monitoring plots...<br>" << endl;
    postMonitoring(extension,minRun,maxRun);
    
    // now we want to produce the exposures and etc.
    cout << "Getting all exposure information...<br>" << endl;
    GetExposures(extension,minRun,maxRun,extension,isLocal);
  }
}

TString makeRunName(int i) {
  TString fname="data/dmtpc_run";
  if (i<100) fname+="000";
  else if (i<1000) fname+="00";
  else if (i<10000) fname+="0";
  fname += i;
  fname += ".root";
  return fname;
}

void GetExposures(TString extension, Int_t minRun, Int_t maxRun, TString moniker, Int_t isLocal){

  DmtpcDataset *d=new DmtpcDataset;
  d->openRootFile(makeRunName(minRun));
  for (int i=minRun+1; i<maxRun+1; i++) {
    d->chain()->Add( makeRunName(i) );
  }
  
  ofstream outputFile;
  TString outputFileName="/usr/local/apache2/htdocs/tmp/dqm/dqm_"+moniker+"/dqm_"+moniker+".txt";
  if(isLocal>0) outputFileName="dqm_"+moniker+"/dqm_"+moniker+".txt";

  ofstream variableExplanationFile;
  TString variableExplanationFileName="/usr/local/apache2/htdocs/tmp/dqm/dqm_"+moniker+"/dqm_"+moniker+"_exp.html";
  if(isLocal>0) variableExplanationFileName="dqm_"+moniker+"/dqm_"+moniker+"_exp.html";

  outputFile.open(outputFileName);
  
  d->getEvent(0);
  TString T0="%b/%d%F";
  T0 += d->event()->timeStamp()->GetYear();
  T0 += "-";
  T0 += d->event()->timeStamp()->GetMonth();
  T0 += "-";
  T0 += d->event()->timeStamp()->GetDay();
  T0 += " ";
  T0 += d->event()->timeStamp()->GetHour();
  T0 +=":";
  T0 += d->event()->timeStamp()->GetMinute();
  T0 +=":";
  T0 += d->event()->timeStamp()->GetSecond();
  
  float pressure = d->event()->experimentConfig("pressure")->currentValue;
  float temperature = d->event()->experimentConfig("temp0")->currentValue+273.15;
  float rho=MaxCamSRIM::density(pressure, 88, temperature); // g/cm3
  float volume = 14.5*14.5*19.5 + 16.5*16.5*19.5; // cm3
  float m=volume*rho; // in g

  variableExplanationFile.open(variableExplanationFileName);
  variableExplanationFile << "<html>" << endl;
  variableExplanationFile << "<body>" << endl << endl;
  variableExplanationFile << "<h1><font color=\"red\">Pressure=</font>d->event()->experimentConfig(\"pressure\")->currentValue for the first event in each run (Torr)</h1>" << endl;
  variableExplanationFile << "<h1><font color=\"red\">Volume=</font>14.5*14.5*19.5 + 16.5*16.5*19.5 (cm^3)</h1>" << endl;
  variableExplanationFile << "<h1><font color=\"red\">Density=</font>MaxCamSRIM::density(pressure, 88, temperature) (g/cm^3)</h1>" << endl;
  variableExplanationFile << "<h1><font color=\"red\">Mass=</font>Volume*DensityMax (g)</h1>" << endl;
  variableExplanationFile << "<h1><font color=\"red\">T0=</font>Hours since first event in this set of runs (hr)</h1>" << endl;
  variableExplanationFile << "<h1><font color=\"red\">dT=</font>Total time in hours that the camera spent taking exposures during this run (hr)</h1>" << endl;
  variableExplanationFile << "<h1><font color=\"red\">exposure=</font>the sum for this run and all previous runs of the exposure time, per run in days, times the Mass in grams (g days)</h1>" << endl;
  variableExplanationFile << "<h1><font color=\"red\">eff=</font>The fraction of time the camera has spent being exposed divided by the total time elapsed thus far for the sum of this run and all previous runs</h1>" << endl;
  variableExplanationFile << "<h1><font color=\"red\">SR#=</font>The spark rate for this run for camera # in hours, counting as sparks the events with average pixel yields above threshold (>THRESHOLD) (sparks/hr)</h1>" << endl;

  variableExplanationFile << "</body>" << endl;
  variableExplanationFile << "</html>" << endl;
  variableExplanationFile.close();
  
  outputFile << "Pressure .... " << pressure << " Torr" << endl;
  outputFile << "Volume ...... " << volume << " cm3" << endl;
  outputFile << "Density ..... " << rho << " g/cm3" << endl;
  outputFile << "Mass ........ " << m << "g"<<endl;

  vector<float> x, y;
  int n=0;
  
  time_t t0=0;
  time_t event_time=0;
  float tot=0;
  int eventStep=d->chain()->GetTree()->GetEntries();
  float totExpo=0;

  Int_t whichRun=minRun;

  for (int i=0; i<d->chain()->GetEntries(); i+=eventStep) {
    d->getEvent(i);
    event_time=d->event()->timeStamp()->Convert();
    if (i==0) t0=event_time;
    event_time = event_time - t0;
    float eventExposure  = d->event()->ccdConfig(0)->exposureTime*1e-3; // s
    float exposureTime = eventExposure * eventStep /  3600.; // h
    if (n) totExpo+=exposureTime;
    tot += exposureTime / 24  * m; // g*day
    x.push_back(event_time);
    y.push_back(tot);
    outputFile << " run=" << d->event()->runNumber() 
	       << "  T0=" << event_time/3600. << "h"
	       << "  dT=" << exposureTime << "h"
	       << "  exposure=" << tot<< "g-day";
    if (n) outputFile << " eff="<< totExpo/event_time*3600;
    
    // get the last time in the event
    Double_t runStartTime=d->event()->timeStamp()->Convert();
    d->getEvent(i+eventStep-1);
    Double_t runEndTime=d->event()->timeStamp()->Convert();
    
    // the spark rates on cameras 1 and 2
    Double_t sparkRates[2]={0.};
    // the thresholds used to get the spark rates on cameras 1 and 2
    Double_t thresholds[2]={0.};
    
    // this gets the spark rate if we're pulling it from the database
    if(isLocal<2) getRunSparkRate(runStartTime,runEndTime,sparkRates,thresholds); 
    if(isLocal==2){
      getRunSparkRateFromFile(extension,runStartTime,runEndTime,sparkRates,thresholds);
    }
    
    // output the spark rates
    outputFile << " <font color=\"red\">SR1=" << setprecision(3) << sparkRates[0] 
	       << "sparks/hr (>" << setprecision(4) << thresholds[0] << ")</font>"
	       << " <font color=\"blue\">SR2=" << setprecision(3) << sparkRates[1]
	       << "sparks/hr (>" << setprecision(4) << thresholds[1] << ")</font>";

    outputFile << endl;
    
    n++;
    ++whichRun;
  }
  
  outputFile.close();
}

void getRunSparkRate(Double_t runStartTime, Double_t runEndTime, Double_t sparkRates[], Double_t thresholds[]){
  // in order for this line to work, have to do some black magic
  // on the LNS machines
  gSystem->Load("libMySQL.so");

  // open a connection to the mysql server
  TMySQLServer  mysql("mysql://mitdm00.mit.edu:/DM_TEST_26544", "dmatter","seedark");

 // pull the ccd rows out of the mysql database
  TMySQLResult *res= mysql.Query("select * from ccd");
  
  Float_t x[100];
  // the values from camera 1
  nt1 = new TTree("nt1","");
  nt1->Branch("data",&x, "val/F:time");
  // the values from camera 2
  nt2 = new TTree("nt2","");
  nt2->Branch("data",&x, "val/F:time");
  
  int n=0;
  Double_t camera1PixelYieldMean=0;
  Double_t camera2PixelYieldMean=0;
  int ncamera1=0;
  int ncamera2=0;
  while ( row=(TMySQLRow*)res->Next() ) {
    
    float avgpixel=atof(row->GetField(3));
    TDatime dt( row->GetField(4) ); 
    float timestamp=ReturnPlottableTime(dt);
    
    x[0]=avgpixel;
    x[1]=timestamp;
    
    if(timestamp<=runEndTime&& timestamp>=runStartTime){
      if(n%2==0){
	nt1->Fill();
	camera1PixelYieldMean+=avgpixel;
	++ncamera1;
      } else {
	nt2->Fill();
	camera2PixelYieldMean+=avgpixel;
	++ncamera2;
      }
    }

    n++;
  }
  
  camera1PixelYieldMean/=ncamera1;
  camera2PixelYieldMean/=ncamera2;

  //  cout << "camera1PixelYieldMean=" << camera1PixelYieldMean << endl;
  //  cout << "camera2PixelYieldMean=" << camera2PixelYieldMean << endl;

  // default averages and cuts for camera 1
  Double_t defaultCamera1PixelYieldSparkCut=1077;
  Double_t defaultCamera1PixelYieldMean=1075.95;
  
  // default averages and cuts for camera 2
  Double_t defaultCamera2PixelYieldSparkCut=1048;
  Double_t defaultCamera2PixelYieldMean=1044.08;

  char camera1SparkCut[1000];
  // an overall shift from run to run
  Double_t thisRunCamera1PixelYieldSparkCut=defaultCamera1PixelYieldSparkCut+(camera1PixelYieldMean-defaultCamera1PixelYieldMean);
  //  cout << "thisRunCamera1PixelYieldSparkCut=" << thisRunCamera1PixelYieldSparkCut << endl;
  sprintf(camera1SparkCut,"val>%f",thisRunCamera1PixelYieldSparkCut);
  //  cout << camera1SparkCut << endl;

  char camera2SparkCut[1000];
  // an overall shift from run to run
  Double_t thisRunCamera2PixelYieldSparkCut=defaultCamera2PixelYieldSparkCut+(camera2PixelYieldMean-defaultCamera2PixelYieldMean);
  //  cout << "thisRunCamera2PixelYieldSparkCut=" << thisRunCamera2PixelYieldSparkCut << endl;
  sprintf(camera2SparkCut,"val>%f",thisRunCamera2PixelYieldSparkCut);
  //  cout << camera2SparkCut << endl;
  
  Long64_t nCamera1Sparks=nt1->Project("","",camera1SparkCut);
  Long64_t nCamera2Sparks=nt2->Project("","",camera2SparkCut);

  Double_t runTimeInHours=(runEndTime-runStartTime)/3600.;
  sparkRates[0]=nCamera1Sparks/runTimeInHours;
  sparkRates[1]=nCamera2Sparks/runTimeInHours;

  thresholds[0]=thisRunCamera1PixelYieldSparkCut;
  thresholds[1]=thisRunCamera2PixelYieldSparkCut;
  
  delete res;
  delete nt1;
  delete nt2;
}

void getRunSparkRateFromFile(TString extension,Double_t runStartTime, Double_t runEndTime, Double_t sparkRates[], Double_t thresholds[]){
  //  cout << "getting the spark rate from the root file we made" << endl;
  
  TFile fdqm("dqm_"+extension+"/dqm_"+extension+".root");
  
  TTree *ccdTree=(TTree*)fdqm.Get("ccd");
  
  Float_t minimumCCD1AveragePixel=ccdTree->GetMinimum("avgpixel");
  Float_t maximumCCD1AveragePixel=ccdTree->GetMaximum("avgpixel");
  
  Float_t minimumCCD2AveragePixel=ccdTree->GetMinimum("avgpixel");
  Float_t maximumCCD2AveragePixel=ccdTree->GetMaximum("avgpixel");
  
  Int_t nbins=100;
  TH1F* hccd1Yields=new TH1F("hccd1Yields","",nbins,minimumCCD1AveragePixel,maximumCCD1AveragePixel);
  TH1F* hccd2Yields=new TH1F("hccd2Yields","",nbins,minimumCCD2AveragePixel,maximumCCD2AveragePixel);
  
  char camera1Cut[1000];
  sprintf(camera1Cut,"(ccdid==0)&&(timestamp>%f)&&(timestamp<%f)",runStartTime,runEndTime);
  char camera2Cut[1000];
  sprintf(camera2Cut,"(ccdid==1)&&(timestamp>%f)&&(timestamp<%f)",runStartTime,runEndTime);
  
  ccdTree->Project("hccd1Yields","avgpixel",camera1Cut);
  ccdTree->Project("hccd2Yields","avgpixel",camera2Cut);
  
  Double_t camera1PixelYieldMean=hccd1Yields->GetMean();
  Double_t camera2PixelYieldMean=hccd2Yields->GetMean();
  
  //  cout << "camera1PixelYieldMean=" << camera1PixelYieldMean << endl;
  //  cout << "camera2PixelYieldMean=" << camera2PixelYieldMean << endl;
  
  // default averages and cuts for camera 1
  Double_t defaultCamera1PixelYieldSparkCut=1077;
  Double_t defaultCamera1PixelYieldMean=1075.95;
  
  // default averages and cuts for camera 2
  Double_t defaultCamera2PixelYieldSparkCut=1048;
  Double_t defaultCamera2PixelYieldMean=1044.08;
  
  char camera1SparkCut[1000];
  // an overall shift from run to run
  Double_t thisRunCamera1PixelYieldSparkCut=defaultCamera1PixelYieldSparkCut+(camera1PixelYieldMean-defaultCamera1PixelYieldMean);
  //  cout << "thisRunCamera1PixelYieldSparkCut=" << thisRunCamera1PixelYieldSparkCut << endl;
  sprintf(camera1SparkCut,"(avgpixel>%f)&&(ccdid==0)&&(timestamp>%f)&&(timestamp<%f)",thisRunCamera1PixelYieldSparkCut,runStartTime,runEndTime);
  //  cout << camera1SparkCut << endl;
  
  char camera2SparkCut[1000];
  // an overall shift from run to run
  Double_t thisRunCamera2PixelYieldSparkCut=defaultCamera2PixelYieldSparkCut+(camera2PixelYieldMean-defaultCamera2PixelYieldMean);
  //  cout << "thisRunCamera2PixelYieldSparkCut=" << thisRunCamera2PixelYieldSparkCut << endl;
  sprintf(camera2SparkCut,"(avgpixel>%f)&&(ccdid==1)&&(timestamp>%f)&&(timestamp<%f)",thisRunCamera2PixelYieldSparkCut,runStartTime,runEndTime);
  //  cout << camera2SparkCut << endl;
  
  Long64_t nCamera1Sparks=ccdTree->Project("","",camera1SparkCut);
  Long64_t nCamera2Sparks=ccdTree->Project("","",camera2SparkCut);
  
  //  cout << "nCamera1Sparks=" << nCamera1Sparks << endl;
  //  cout << "nCamera2Sparks=" << nCamera2Sparks << endl;
  
  Double_t runTimeInHours=(runEndTime-runStartTime)/3600.;
  sparkRates[0]=nCamera1Sparks/runTimeInHours;
  sparkRates[1]=nCamera2Sparks/runTimeInHours;

  thresholds[0]=thisRunCamera1PixelYieldSparkCut;
  thresholds[1]=thisRunCamera2PixelYieldSparkCut;

  fdqm.Close();
}

void postMonitoring(TString extension, Int_t minRun, Int_t maxRun){

  // data monitoring plot limits
  Float_t min_pressure=74; // torr
  Float_t max_pressure=77; // torr
  Float_t min_anode=0.65; // torr
  Float_t max_anode=0.80; // torr
  Float_t min_drift=4.95; // torr
  Float_t max_drift=5.05; // torr

  // open a canvas plot
  TCanvas *c1=new TCanvas("c1","c1",0,0,800,800);

  DmtpcDataset *d=new DmtpcDataset;
  d->openRootFile(makeRunName(minRun));
  for (int i=minRun+1; i<maxRun+1; i++) {
    d->chain()->Add( makeRunName(i) ); 
  }

  d->getEvent(0);

  int eventStep=1;

  // dqm information output file
  TFile *outputTreeFile=new TFile("dqm_"+extension+"/dqm_"+extension+".root","RECREATE");
  // the CCD tree
  Float_t xccd[1000];
  TTree *ccdTree=new TTree("ccd","");
  ccdTree->Branch("ccd",&xccd,"temperature/F:exposure/F:daqtime/F:avgpixel/F:timestamp/F:ccdid/F");
  // the pressure tree
  Float_t xpressure[1000];
  TTree *pressureTree=new TTree("pressure","");
  pressureTree->Branch("pressure",&xpressure,"val/F:rms:set:time");
  // the pressure tree
  Float_t xwire_hv[1000];
  TTree *wire_hvTree=new TTree("wire_hv","");
  wire_hvTree->Branch("wire_hv",&xwire_hv,"val/F:rms:set:time");
  // the pressure tree
  Float_t xmesh_hv[1000];
  TTree *mesh_hvTree=new TTree("mesh_hv","");
  mesh_hvTree->Branch("mesh_hv",&xmesh_hv,"val/F:rms:set:time");
  // the temp0 tree
  Float_t xtemp0[1000];
  TTree *temp0Tree=new TTree("temp0","");
  temp0Tree->Branch("temp0",&xtemp0,"val/F:rms:set:time");

  // external slow conditions
  TGraph *pressureGraph=new TGraph();
  TGraph *anodeHVGraph=new TGraph();
  TGraph *driftHVGraph=new TGraph();
  TGraph *temp0Graph=new TGraph();
  
  // ccd info
  // exposure info
  TGraph *ccd1ExposureGraph=new TGraph();
  TGraph *ccd2ExposureGraph=new TGraph();

  // temperature info
  TGraph *ccd1TemperatureGraph=new TGraph();
  TGraph *ccd2TemperatureGraph=new TGraph();

  // temperature info
  TGraph *ccd1YieldGraph=new TGraph();
  TGraph *ccd2YieldGraph=new TGraph();
  
  Int_t counter=0;
  for (int i=0; i<d->chain()->GetEntries(); i+=eventStep) {
    d->getEvent(i);
   
    Float_t time = ReturnPlottableTime((*d->event()->timeStamp()));
    Float_t pressure = d->event()->experimentConfig("pressure")->currentValue;
    Float_t temp0 = d->event()->experimentConfig("temp0")->currentValue; 
    Float_t driftHV = d->event()->experimentConfig("driftHV")->currentValue; 
    Float_t anodeHV = d->event()->experimentConfig("anodeHV")->currentValue; 
    Float_t ccd1Exposure  = d->event()->ccdConfig(0)->exposureTime; 
    Float_t ccd1Temperature  = d->event()->ccdConfig(0)->CCDTemp;
    Float_t ccd1DAQTime  = d->event()->ccdConfig(0)->daqTime;
    Float_t ccd2Exposure  = d->event()->ccdConfig(1)->exposureTime;
    Float_t ccd2Temperature  = d->event()->ccdConfig(1)->CCDTemp;
    Float_t ccd2DAQTime  = d->event()->ccdConfig(1)->daqTime;

    pressureGraph->SetPoint(counter,time,pressure);
    anodeHVGraph->SetPoint(counter,time,anodeHV);
    driftHVGraph->SetPoint(counter,time,driftHV);
    temp0Graph->SetPoint(counter,time,temp0);
    ccd1ExposureGraph->SetPoint(counter,time,ccd1Exposure);
    ccd2ExposureGraph->SetPoint(counter,time,ccd2Exposure);
    ccd1TemperatureGraph->SetPoint(counter,time,ccd1Temperature);
    ccd2TemperatureGraph->SetPoint(counter,time,ccd2Temperature);
    
    TH2F *hccd1Image = d->event()->ccdData(0);
    TH1F* hccd1Yield=myCreateYieldHisto(hccd1Image);
    Double_t ccd1Yield=hccd1Yield->GetMean();
    hccd1Image->Clear();
    hccd1Yield->Clear();
    
    TH2F *hccd2Image = d->event()->ccdData(1);
    TH1F* hccd2Yield=myCreateYieldHisto(hccd2Image);
    Double_t ccd2Yield=hccd2Yield->GetMean();
    hccd2Image->Clear();
    hccd2Yield->Clear();

    ccd1YieldGraph->SetPoint(counter,time,ccd1Yield);
    ccd2YieldGraph->SetPoint(counter,time,ccd2Yield);

    // fill the ccd tree
    // for ccd camera 1
    xccd[0]=ccd1Temperature; xccd[1]=ccd1Exposure; xccd[2]=ccd1DAQTime; xccd[3]=ccd1Yield; xccd[4]=time; xccd[5]=0;
    ccdTree->Fill();
    // for ccd camera 2
    xccd[0]=ccd2Temperature; xccd[1]=ccd2Exposure; xccd[2]=ccd2DAQTime; xccd[3]=ccd2Yield; xccd[4]=time; xccd[5]=1;
    ccdTree->Fill();

    // fill the mesh hv tree
    xmesh_hv[0]=driftHV; xmesh_hv[1]=d->event()->experimentConfig("driftHV")->currentValueRMS; xmesh_hv[2]=d->event()->experimentConfig("driftHV")->setValue; xmesh_hv[3]=time;
    mesh_hvTree->Fill();
    
    // fill the pressure tree
    xpressure[0]=pressure; xpressure[1]=d->event()->experimentConfig("pressure")->currentValueRMS; xpressure[2]=d->event()->experimentConfig("pressure")->setValue; xpressure[3]=time;
    pressureTree->Fill();

    // fill the wire hv tree
    xwire_hv[0]=anodeHV; xwire_hv[1]=d->event()->experimentConfig("anodeHV")->currentValueRMS; xwire_hv[2]=d->event()->experimentConfig("anodeHV")->setValue; xwire_hv[3]=time;
    wire_hvTree->Fill();

    // fill the temp0 tree
    xtemp0[0]=temp0; xtemp0[1]=d->event()->experimentConfig("temp0")->currentValueRMS; xtemp0[2]=d->event()->experimentConfig("temp0")->setValue; xtemp0[3]=time;
    temp0Tree->Fill();


    ++counter;
  }
  
  // make the pressure graph
  DressGraph(*pressureGraph,"Pressure (Torr)");
  c1->Clear();
  pressureGraph->SetMinimum(min_pressure);
  pressureGraph->SetMaximum(max_pressure);
  pressureGraph->Draw("alp");
  c1->SaveAs("dqm_"+extension+"/pressure_"+extension+".gif");

  // make the anode plot
  DressGraph(*anodeHVGraph,"Anode (kV)");
  c1->Clear();
  anodeHVGraph->SetMinimum(min_anode);
  anodeHVGraph->SetMaximum(max_anode);
  anodeHVGraph->Draw("alp");
  c1->SaveAs("dqm_"+extension+"/anode_"+extension+".gif");

  // make the drift plot
  DressGraph(*driftHVGraph,"Drift (kV)");
  c1->Clear();
  driftHVGraph->SetMinimum(min_drift);
  driftHVGraph->SetMaximum(max_drift);
  driftHVGraph->Draw("alp");
  c1->SaveAs("dqm_"+extension+"/drift_"+extension+".gif");

  // make the temp0 plot
  DressGraph(*temp0Graph,"temp0");
  c1->Clear();
  temp0Graph->Draw("ap");
  c1->SaveAs("dqm_"+extension+"/value_in_temp0_"+extension+".gif");
  
  // make the ccd temperature comparison plot
  DressGraph(*ccd1TemperatureGraph,"ccd temperature");
  DressGraph(*ccd2TemperatureGraph,"ccd temperature",kBlue);
  c1->Clear();
  ccd1TemperatureGraph->Draw("ap");
  ccd2TemperatureGraph->Draw("psame");
  c1->SaveAs("dqm_"+extension+"/temperature_in_ccd_"+extension+".gif");

  // make the ccd exposure comparison plot
  DressGraph(*ccd1ExposureGraph,"ccd exposure");
  DressGraph(*ccd2ExposureGraph,"ccd exposure",kBlue);
  c1->Clear();
  ccd1ExposureGraph->Draw("ap");
  ccd2ExposureGraph->Draw("psame");
  c1->SaveAs("dqm_"+extension+"/exposure_in_ccd_"+extension+".gif");
  
  // make the ccd yield comparison plot
  DressGraph(*ccd1YieldGraph,"ccd average pixel yield");
  DressGraph(*ccd2YieldGraph,"ccd average pixel yield",kBlue);
  ccd1YieldGraph->SetMinimum(1000);
  ccd1YieldGraph->SetMaximum(1800);
  c1->Clear();
  ccd1YieldGraph->Draw("ap");
  ccd2YieldGraph->Draw("psame");
  c1->SaveAs("dqm_"+extension+"/avg_pixel_in_ccd_"+extension+".gif");
  
  ccdTree->Write();
  pressureTree->Write();
  mesh_hvTree->Write();
  wire_hvTree->Write();
  temp0Tree->Write();
  outputTreeFile->Close();
}

void DressGraph(TGraph &gr, TString yTitle, Int_t color){
  TString T0="%b/%d %I%p";
  gr.GetXaxis()->SetTimeDisplay(1);
  gr.GetXaxis()->SetTimeFormat((const char*)T0);
  gr.GetXaxis()->SetTitle("Time");
  gr.GetYaxis()->SetTitle(yTitle);
  gr.SetLineColor(color);
  gr.SetMarkerColor(color);
  gr.SetMarkerStyle(2);
  gr.SetMarkerSize(0.5);
}

// a copy, as of 5/19/2009, of the createYieldHisto()
// in MaxCamImageTools.cc v1.37
TH1F* myCreateYieldHisto(TH1* h) {

    int nbin=h->GetNbinsX()*h->GetNbinsY()*h->GetNbinsZ();
    
    Float_t       min = h->GetMinimum();
    Float_t       max = h->GetMaximum();
    
    TH1F *yieldHisto=new TH1F("yieldHisto","",100, min,max);
    for (int i=0; i<nbin; i++) {
        yieldHisto->Fill( h->GetBinContent(i) );
    }
    return yieldHisto;
}
