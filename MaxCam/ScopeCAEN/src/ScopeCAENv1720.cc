//This script is created to utilize the CAEN v1720 Digintazier
//Here we set all classes and function that are need to get data from the digitaizer
//	
//Pietro Giampa, RHUL 2013

#include <iostream>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#include <CAENDigitizer.h>
#include <CAENDigitizerType.h>
#include "../include/ScopeCAENv1720.hh"	

using std::cout;
using std::endl;
using std::ofstream;
using std::ios;
using std::flush;

const int ScopeCAENv1720::MaxNChannels = 4;
const int ScopeCAENv1720::MaxBits = 12;
const int ScopeCAENv1720::MAXNB = 1;

//Set name of the digitizer we are using
char* ScopeCAENv1720::GetName() { return "CAENv1720";}

ScopeCAENv1720::~ScopeCAENv1720(){

}
//Main functions that controll the CAEN digitizer
ScopeCAENv1720::ScopeCAENv1720() {

  //Set Up Name and Directory of the output file
  fDataPath = "output/";
  fFileName =  fDataPath + "dmtpc_CAEN_";
  int RunNumber;
  cout << "Enter Run Number (0-9999) : " << endl;
  std::cin >> RunNumber;

  //Check Input Run Number Is Correct
  while (RunNumber > 9999) {
    cout << "Invalid Run Number !!! Must Best Between 0 and 9999" << endl;
    std::cin >> RunNumber;
  }

  //Set Correct Name with Set Run Number
  if (RunNumber < 10) {
    fFileName += "000";
    fFileName += RunNumber;
  } else if (RunNumber >= 10 && RunNumber < 100) {
    fFileName += "00";
    fFileName += RunNumber;
  } else if (RunNumber >= 100 && RunNumber < 1000) {
    fFileName += "0";
    fFileName += RunNumber;
  } else if (RunNumber >= 1000) {
    fFileName += RunNumber;
  }
  fFileName += ".root";

  //Call Readout function
  readout(fFileName);

}

//This function is dedicated to check the keyboard input
//for the Basic Readout functio, but it can also be utilized
//anywhere esle and is build based on CAENKeyComm.cc
//Stole from CAEN, (fuck those guys)
int
ScopeCAENv1720::KeyCheck(){

  //Define variable
  int c = 0;

  //Check with CAENKeyComm.cc
  //For function kbhit
  if(!kbhit())
    return 0;

  //Check with CAENKeyComm.cc
  //Forfunction getch
  c = getch();
  switch (c) {
  case 's':
    return 9;
    break;
  case 'k':
    return 1;
    break;
  case 'q':
    return 2;
    break;
  }

  return 0;
}

//basic readout function for DPP_CI firmware 
//based on the sample given by CAEN
void
ScopeCAENv1720::readout(TString fFileName) {

  //----------------------------------------------------------------------------------------------+
  //Initializing
  CAEN_DGTZ_ErrorCode ret;
  int handle;
  CAEN_DGTZ_BoardInfo_t BoardInfo;
  CAEN_DGTZ_EventInfo_t eventInfo;
  //void *Evt = NULL;
  CAEN_DGTZ_UINT16_EVENT_t *Evt = NULL;
  char *buffer = NULL;
  int i;
  int count[MAXNB];
  count[0] = 0;
  char * evtptr = NULL;
  uint32_t size,bsize;
  uint32_t numEvents;
  int Quit = 0;
  int acqRun = 0;
  int MajorNumber = 0;
  i = sizeof(CAEN_DGTZ_TriggerMode_t);
  
  //Waveform
  TH1S* WFhistCH0 = new TH1S("WFhistCH0", "Waveform From CAENv1720", 4096, 0, 4096);            //used to stored waveform CH0
  TH1S* WFhistCH1 = new TH1S("WFhistCH1", "Waveform From CAENv1720", 4096, 0, 4096);            //used to stored waveform CH1
  TH1S* WFhistCH2 = new TH1S("WFhistCH2", "Waveform From CAENv1720", 4096, 0, 4096);            //used to stored waveform CH2
  TH1S* WFhistCH3 = new TH1S("WFhistCH3", "Waveform From CAENv1720", 4096, 0, 4096);            //used to stored waveform CH3
  TH1S* WFhistCH4 = new TH1S("WFhistCH4", "Waveform From CAENv1720", 4096, 0, 4096);            //used to stored waveform CH4
  TH1S* WFhistCH5 = new TH1S("WFhistCH5", "Waveform From CAENv1720", 4096, 0, 4096);            //used to stored waveform CH5
  TH1S* WFhistCH6 = new TH1S("WFhistCH6", "Waveform From CAENv1720", 4096, 0, 4096);            //used to stored waveform CH6
  TH1S* WFhistCH7 = new TH1S("WFhistCH7", "Waveform From CAENv1720", 4096, 0, 4096);            //used to stored waveform CH7

  TH1S *Waveform = new TH1S("Waveform", "Waveform From CAENv1720", 4096, 0, 4096);        //used to save just the last waveform of the readout, for a quick test that nothing is wrong
  
  //Trigger Tag Time
  Int_t TttCH0=0;   //save Trigger time tag for CH0 (CAEN units [?])
  Int_t TttCH1=0;   //save Trigger time tag for CH1 (CAEN units [?])
  Int_t TttCH2=0;   //save Trigger time tag for CH2 (CAEN units [?])
  Int_t TttCH3=0;   //save Trigger time tag for CH3 (CAEN units [?])
  Int_t TttCH4=0;   //save Trigger time tag for CH4 (CAEN units [?])
  Int_t TttCH5=0;   //save Trigger time tag for CH5 (CAEN units [?])
  Int_t TttCH6=0;   //save Trigger time tag for CH6 (CAEN units [?])
  Int_t TttCH7=0;   //save Trigger time tag for CH7 (CAEN units [?])

  //----------------------------------------------------------------------------------------------+
  //Set Up output file
  TFile *Data = new TFile(fFileName, "RECREATE");
  TTree *dmtpc = new TTree("dmtpc", "CAEN Digitizer");
  dmtpc->Branch("Waveform_Channel0", "TH1S", &WFhistCH0);
  dmtpc->Branch("Waveform_Channel1", "TH1S", &WFhistCH1);
  dmtpc->Branch("Waveform_Channel2", "TH1S", &WFhistCH2);
  dmtpc->Branch("Waveform_Channel3", "TH1S", &WFhistCH3);
  dmtpc->Branch("Waveform_Channel4", "TH1S", &WFhistCH4);
  dmtpc->Branch("Waveform_Channel5", "TH1S", &WFhistCH5);
  dmtpc->Branch("Waveform_Channel6", "TH1S", &WFhistCH6);
  dmtpc->Branch("Waveform_Channel7", "TH1S", &WFhistCH7);
  dmtpc->Branch("TriggerTime_Channel0", &TttCH0);
  dmtpc->Branch("TriggerTime_Channel1", &TttCH1);
  dmtpc->Branch("TriggerTime_Channel2", &TttCH2);
  dmtpc->Branch("TriggerTime_Channel3", &TttCH3);
  dmtpc->Branch("TriggerTime_Channel4", &TttCH4);
  dmtpc->Branch("TriggerTime_Channel5", &TttCH5);
  dmtpc->Branch("TriggerTime_Channel6", &TttCH6);
  dmtpc->Branch("TriggerTime_Channel7", &TttCH7);

  //----------------------------------------------------------------------------------------------+
  //Open Digitizer, and get handle
  //Loop through Boards
  
  ret = CAEN_DGTZ_OpenDigitizer(CAEN_DGTZ_PCI_OpticalLink, 0, 0, 0x32100000, &handle);  //IMPORTANT : CHECK BOARD ADDRESS 
  if (ret) { printf("Can't open digitizer\n"); }
  ret = CAEN_DGTZ_GetInfo(handle, &BoardInfo);
  if (ret) { printf("Can't read board info\n"); }
  printf("\nConnected to CAEN Digitizer Model %s, recognized as board %d\n", BoardInfo.ModelName, 0);
  printf("ROC FPGA Release is %s\n", BoardInfo.ROC_FirmwareRel);
  printf("AMC FPGA Release is %s\n", BoardInfo.AMC_FirmwareRel);
  printf("Channels %d\n", BoardInfo.Channels);
  printf("Serial Number %d\n", BoardInfo.SerialNumber);
  
  sscanf(BoardInfo.AMC_FirmwareRel, "%d", &MajorNumber);
  if (MajorNumber != 130) {
    printf("This digitizer has not a DPP-CI firmware\n");
  }
  
  //----------------------------------------------------------------------------------------------+
  //Set up digitizer properties
  //loop through boards
  
  ret = CAEN_DGTZ_Reset(handle);                                                  /* Reset Digitizer */
  ret = CAEN_DGTZ_SetRecordLength(handle,4096);                                   /* Set the lenght of each waveform (in samples) */
  ret = CAEN_DGTZ_SetPostTriggerSize(handle, 80);                                    /* Set Post Trigger size */
  ret = CAEN_DGTZ_SetIOLevel(handle, CAEN_DGTZ_IOLevel_NIM);                         /* Set IO level */
  ret = CAEN_DGTZ_SetChannelEnableMask(handle,0xFF);                              /* Enable channel 0 */

  ret = CAEN_DGTZ_SetChannelDCOffset(handle, -1, 0xFFFF);                            /* Set DC offset  */
  ret = CAEN_DGTZ_SetChannelTriggerThreshold(handle,-1,3000);                     /* Set selfTrigger threshold (need to check this)*/
  ret = CAEN_DGTZ_SetChannelSelfTrigger(handle,CAEN_DGTZ_TRGMODE_ACQ_ONLY,-1);    /* Set trigger on channel 0 to be ACQ_ONLY */
  ret = CAEN_DGTZ_SetSWTriggerMode(handle,CAEN_DGTZ_TRGMODE_ACQ_ONLY);            /* Set the behaviour when a SW tirgger arrives */

  ret = CAEN_DGTZ_SetMaxNumEventsBLT(handle,3);                                /* Set the max number of events to transfer in a sigle readout */
  ret = CAEN_DGTZ_SetAcquisitionMode(handle,CAEN_DGTZ_SW_CONTROLLED);          /* Set the acquisition mode */
  
  cout << "Digitizer Programmed \n" << endl;
  
  //----------------------------------------------------------------------------------------------+
  //Allocate memory 
  //Loop through Boards
  
  ret = CAEN_DGTZ_AllocateEvent(handle, &Evt);
  ret = CAEN_DGTZ_MallocReadoutBuffer(handle,&buffer,&size);
  if (ret) { printf("Problem With Malloc Buffer"); }
  
  cout << "Memory Allocated \n" << endl;
  
  //----------------------------------------------------------------------------------------------+
  //Acquisition Data
  printf("Enter Command : \n");
  printf("Press 's' to start the acquisition\n");
  printf("Press 'k' to stop the acquisition\n");
  printf("Press 'r' to restart the acquisition\n");
  printf("Press 't' to send software trigger\n");
  printf("Press 'q' to quit\n");
  printf("\n");
  //Loop to set DAQ
  while (!Quit){
    //Check what has been typed
    if (kbhit()){
      //set up commands
      char c;
      c = getch();
      //Quit Command
      if (c=='q') {
	Quit=1;
	cout << "Quit Acquisition \n" << endl;
      }
      //Send SW Trigger
      if (c=='t') {
	ret = CAEN_DGTZ_SendSWtrigger(handle);
	if (ret) { printf("Error SW trigger %d\n", ret); }
	printf("Send Software Trigger\n");
      }
      //Restart
      if (c=='r'){
	ret = CAEN_DGTZ_SWStopAcquisition(handle);
	if (ret) { printf("Error Stopping Acquisition %d\n", ret); }
	printf("Restarted\n");
	ret = CAEN_DGTZ_ClearData(handle);
	if (ret) { printf("Error Cleaning Data %d\n", ret); }
      }
      //Start Acquisition
      if (c=='s'){
	ret = CAEN_DGTZ_SWStartAcquisition(handle);
	if (ret) { printf("Error Starting Acquisition %d\n", ret); }
	printf("Acquisition Started for Board 0\n");
	acqRun=1;
      }
      //Stop Acquisition
      if (c=='k'){
	ret = CAEN_DGTZ_SWStopAcquisition(handle);
	if (ret) { printf("Error Stopping Acquisition %d\n", ret); }
	printf("Acquisition Stopped for Board 0\n");
	acqRun=0;
      }
    }
    if (!acqRun) {
      Sleep(10);
      continue;
    }
    
    //----------------------------------------------------------------------------------------------+
    //Data Analysis
    
    //get data from board
    ret = CAEN_DGTZ_ReadData(handle,CAEN_DGTZ_SLAVE_TERMINATED_READOUT_MBLT,buffer,&bsize); /* Read the buffer from the digitizer */
    if (ret) { printf("Readout Error %d\n", ret); }
    
    //get num events
    ret = CAEN_DGTZ_GetNumEvents(handle,buffer,bsize,&numEvents);
    if (ret) { printf("Readout Error %d\n", ret); }
    count[0] +=numEvents;
    
    //Loop through events block
    for (i=0;i<numEvents;i++) {
      
      //Clear Histo
      WFhistCH0->Reset();
      WFhistCH1->Reset();
      WFhistCH2->Reset();
      WFhistCH3->Reset();
      WFhistCH4->Reset();
      WFhistCH5->Reset();
      WFhistCH6->Reset();
      WFhistCH7->Reset();      

      //get event info
      ret = CAEN_DGTZ_GetEventInfo(handle,buffer,bsize,i,&eventInfo,&evtptr);
      
      //decode event
      ret = CAEN_DGTZ_DecodeEvent(handle,evtptr,&Evt);
      
      //Set Event into tree
      //loop through channels
      for (int ch=0; ch<8; ch++){
	//Set Event into tree
        uint32_t p=0;
	
	//check if the channel has been initialized
	if (Evt->ChSize[ch] > 0){

          //loop through the channel samples
          while (p<Evt->ChSize[ch]) {
            if (ch == 0) { WFhistCH0->SetBinContent(p+1, (int)Evt->DataChannel[ch][p]);}
            else if (ch == 1) { WFhistCH1->SetBinContent(p+1, (int)Evt->DataChannel[ch][p]); }
            else if (ch == 2) { WFhistCH2->SetBinContent(p+1, (int)Evt->DataChannel[ch][p]); }
            else if (ch == 3) { WFhistCH->SetBinContent(p+1, (int)Evt->DataChannel[ch][p]); }
            else if (ch == 4) { WFhistCH4->SetBinContent(p+1, (int)Evt->DataChannel[ch][p]); }
            else if (ch == 5) { WFhistCH5->SetBinContent(p+1, (int)Evt->DataChannel[ch][p]); }
            else if (ch == 6) { WFhistCH6->SetBinContent(p+1, (int)Evt->DataChannel[ch][p]); }
            else if (ch == 7) { WFhistCH7->SetBinContent(p+1, (int)Evt->DataChannel[ch][p]); }
            p++;
          }

	  if (ch == 0) { TttCH0 = (int)eventInfo.TriggerTimeTag;}
	  else if (ch == 1) { TttCH1 = (int)eventInfo.TriggerTimeTag;}
	  else if (ch == 2) { TttCH2 = (int)eventInfo.TriggerTimeTag; }
	  else if (ch == 3) { TttCH3 = (int)eventInfo.TriggerTimeTag;}
	  else if (ch == 4) { TttCH4 = (int)eventInfo.TriggerTimeTag }
	  else if (ch == 5) { TttCH5 = (int)eventInfo.TriggerTimeTag; }
	  else if (ch == 6) { TttCH6 = (int)eventInfo.TriggerTimeTag;}
	  else if (ch == 7) { TttCH7 = (int)eventInfo.TriggerTimeTag; }

        }
      }
      
      //Store in Data/TTree
      dmtpc->Fill();
      
      //free event
      ret = CAEN_DGTZ_FreeEvent(handle,&Evt);
      
    }//end lopp through events
    
  }//end data taking
  
  printf("Total Number of Retrived Events %d\n", count[0]);
  cout << "Data Analysis Done \n" << endl;
  
  //Open and Write File
  dmtpc->Print();
  Waveform->Write();
  Data->Write();
  Data->Close();
  
  //----------------------------------------------------------------------------------------------+
  //Reset the board
  //loop through boards
  
  //Reset memory
  ret = CAEN_DGTZ_FreeReadoutBuffer(&buffer);
  ret = CAEN_DGTZ_CloseDigitizer(handle);
  if (ret) { printf("Error Closing Digitizer %d\n", ret); }
  
  cout << "Board Reset and Closed" << endl;
  
}
