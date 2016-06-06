//This script is created to utilize the CAEN v1720 Digintazier
//Here we set all classes and function that are need to get data from the digitaizer
//	
//Pietro Giampa, RHUL 
//June 2012

#include <iostream>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#include <CAENDigitizer.h>
#include <CAENDigitizerType.h>
#include <CAENComm.h>
#include "CAENKeyComm.hh"
#include "ScopeCAENv1720.hh"    
#include "ScopeHandlerCAEN.hh"	
#include "TStopwatch.h"		

using std::cout;
using std::endl;
using std::ofstream;
using std::ios;
using std::flush;

const int ScopeCAENv1720::MaxNChannels = 8;
const int ScopeCAENv1720::MaxBits = 14;
const int ScopeCAENv1720::MaxNb = 1;

//Set name of the digitizer we are using
char* ScopeCAENv1720::GetName() { return "CAENv1720";}

ScopeCAENv1720::~ScopeCAENv1720(){

}
//Main functions that controll the CAEN digitizer
ScopeCAENv1720::ScopeCAENv1720() {

  /*
  //Call function to set parameters
  DigitizerParams_t *Params;
  Params = SetStandardParameters();

  //Call function to set DPP parameters
  CAEN_DGTZ_DPP_PHA_Params_t *ddpParams;
  ddpParams = SetDPPParameters();

  //Call function to get handle
  int* handle;
  handle = GetHandle(Params);

  cout << "\t" << endl;
  
  //Quick test that the function sets the DPP params 
  //And that we can easily call new variabs 
  cout << "Trigger Threshold value set at : " << ddpParams[0].thr[0] << "\t" << endl;
  cout << "VMEBaseAddress : " << Params[0].VMEBaseAddress << endl;

  cout << "\t" << endl;

  //Basic data acquisition
  BasicReadout(handle, Params, ddpParams);
  */

  //Basic data acquisition
  //No DPP setting used yet
  //Send 10 SW triggers and printf their time tag
  BasicReadout2();
}

//This function is used to return the handle for each board
//in the system using function from the CAEN library
int* ScopeCAENv1720::GetHandle(DigitizerParams_t* Params){

  cout << "Opening Scope : " << GetName() << endl;

  //Define variables
  int *handle = new int[MaxNb];
  int GHN = 0;
  CAEN_DGTZ_BoardInfo_t BoardInfo;

  //Loop through boards
  for (int b=0; b<MaxNb; b++) {

    //Extract handle from boards
    GHN = CAEN_DGTZ_OpenDigitizer(Params[b].LinkType, 0, b, Params[b].VMEBaseAddress, &handle[b]);

    //Check if we connect with digitizer
    if (GHN) {
      printf("Failed to extract handle\n");
      return 0;
    }
    
    //Now that we have the handle that identifies the device
    GHN = CAEN_DGTZ_GetInfo(handle[b],&BoardInfo);
    if (GHN) {
      printf("Can't read board info\n");
      return 0;
    }else {
      printf("\nConnected to CAEN Digitizer Model %s, recognized as board %d\n", BoardInfo.ModelName, b);
      printf("ROC FPGA Release is %s\n", BoardInfo.ROC_FirmwareRel);
      printf("AMC FPGA Release is %s\n", BoardInfo.AMC_FirmwareRel);
    }

  }
  return handle;
}

//This function is used to set up a series of DPP parameters
//That we are going to use to open/read/trigger the CAEN digitizer
CAEN_DGTZ_DPP_PHA_Params_t* ScopeCAENv1720::SetDPPParameters(){

  //define functions CAENDigitizer API
  CAEN_DGTZ_DPP_PHA_Params_t *DPPParams = new CAEN_DGTZ_DPP_PHA_Params_t[MaxNb];
  
  //Set the parameters for each channel
  for (int bb=0; bb<MaxNb; bb++){
    for (int ch=0; ch<MaxNChannels; ch++){
      DPPParams[bb].thr[ch] = 200;
      DPPParams[bb].k[ch] = 1000;
      DPPParams[bb].m[ch] = 500;
      DPPParams[bb].M[ch] = 200;
      DPPParams[bb].ftd[ch] = 30;
      DPPParams[bb].a[ch] = 2;
      DPPParams[bb].b[ch] = 100;
      DPPParams[bb].trgho[ch] = 600;
      DPPParams[bb].nsbl[ch] = 2;
      DPPParams[bb].nspk[ch] = 2;
      DPPParams[bb].pkho[ch] = 770;
      DPPParams[bb].blho[ch] = 100;
      DPPParams[bb].enf[ch] = 1.0;      
    }
  }
  return DPPParams;
}

//This function is used to set up a series of DPP parameters
//That we are going to use to open/read/trigger the CAEN digitizer
//This set a standard value for each parameter that can be change by the user
DigitizerParams_t* ScopeCAENv1720::SetStandardParameters(){

  //Set variables for Parameters
  DigitizerParams_t *Params = new DigitizerParams_t[MaxNb];

  //Set the parameters for each channel
  for (int b=0; b<MaxNb; b++){
    for (int ch=0; ch<MaxNChannels; ch++){
      
      //Basic parameters needed by the board
      Params[b].LinkType = CAEN_DGTZ_PCI_OpticalLink;
      Params[b].VMEBaseAddress = 0;
      Params[b].IOlev = CAEN_DGTZ_IOLevel_TTL;
      
      //Acquisition Parameters
      Params[b].AcqMode = CAEN_DGTZ_DPP_ACQ_MODE_Oscilloscope; //we can also have 
                                                               //CAEN_DGTZ_DPP_ACQ_MODE_List or 
                                                               //CAEN_DGTZ_DPP_ACQ_MODE_Mixed
      Params[b].RecordLength = 4096;           //Num of samples of the waveforms 
                                            //only for Oscilloscope mode
      Params[b].ChannelMask = 0xFF;
      Params[b].EventAggr = 10;
      Params[b].PulsePolarity = CAEN_DGTZ_PulsePolarityNegative;
    }
  }
  return Params;
}

//This function allows the user to reset the digitizer
void
ScopeCAENv1720::ResetDigi(int* handle){
  
  int reset=0;
  cout << "Resetting Digitizier : " << GetName() << endl;
  
  //Loop through Boards
  for (int b=0; b<MaxNb; b++){
    //Calling function from CAEN library to reset Scope
    reset |= CAEN_DGTZ_Reset(handle[b]);
  }
  //make sure we are resetting the Digitizer
  if (reset) { printf("Couldn't Reset Digitizer"); }
}

//Set the Handle recieved by the openScope function
void
ScopeCAENv1720::setScopeHandle(int h) { _h = h; }

//Function allows us to call the scope Handle anytime
int
ScopeCAENv1720::getScopeHandle() { return _h; }

//This function is dedicated to check the keyboard input 
//for the Basic Readout functio, but it can also be utilized 
//anywhere esle and is build based on CAENKeyComm.cc
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

//This function allows the user to quit the data acquisition
//Process by making sure that the Digitizer is properly closed 
//and all the memory are cleared 
int
ScopeCAENv1720::StopReadout(int* handle){

  //Define variable
  char *buffer=NULL;
  int STOP=0, c=0;

  //Calling function from CAEN library to free buffer
  STOP = CAEN_DGTZ_FreeReadoutBuffer(&buffer);

  //Loop through boards
  for (int b=0; b<MaxNb; b++){
    
    //Calling function from CAEN library to close the digitizer
    STOP = CAEN_DGTZ_CloseDigitizer(handle[b]);
    printf("Press 'Enter' key to exit\n");
    c = getchar();
    return 0;
  }
  return 0;
}

//This function makes use of the Digitizer API functions to set up
//The Digitizer's initial configuration and all we need to set for read out
//And the data acquisition phase, including DPP parameters for boards 
//Connected via USB
int 
ScopeCAENv1720::DigiConfig(int* handle, DigitizerParams_t* Params, CAEN_DGTZ_DPP_PHA_Params_t* DPPParams) {

  //Set variables for this function
  int conf = 0;
  int i=0;
  
  //Loop through boards
  for (int b=0; b<MaxNb; b++){

    //Reset digitizer, just in case
    //so that we make sure there is no conflict with the previous setting
    conf |= CAEN_DGTZ_Reset(handle[b]);
    
    //Write out the Register for the Digi 
    //Not sure about the setting for this function need to run it to check
    conf |= CAEN_DGTZ_WriteRegister(handle[b], 0x8000, 0x01000114);
    
    //Set the DPP acquisition mode
    //If one would like to change this always check the manual first
    //The 3rd variables can also be set to one of these 4 variables
    //CAEN_SGTZ_DPP_SAVE_EnergyOnly
    //CAEN_SGTZ_DPP_SAVE_TimeOnly
    //CAEN_SGTZ_DPP_SAVE_EnergyAndTime
    //CAEN_SGTZ_DPP_SAVE_None  
    conf |= CAEN_DGTZ_SetDPPAcquisitionMode(handle[b], Params[0].AcqMode, CAEN_DGTZ_DPP_SAVE_PARAM_EnergyAndTime);
    
    //Set Digitizer acquisition mode via Params
    //on top of the handle thw function needs to be fed with either
    //CAEN_DGTZ_SW_CONTROLLED or CAEN_DGTZ_S_IN_CONTROLLED
    conf |= CAEN_DGTZ_SetAcquisitionMode(handle[b], CAEN_DGTZ_SW_CONTROLLED);
    
    //Set Number of samples for each waveform
    conf |= CAEN_DGTZ_SetRecordLength(handle[b], Params[0].RecordLength);
    
    //Set I/O level for the Digitizer
    //on top of the handle thw function needs to be fed with either
    //CAEN_DGTZ_IOLevel_NIM or CAEN_DGTZ_IOLevel_TTL
    conf |= CAEN_DGTZ_SetIOLevel(handle[b], Params[0].IOlev);
    
    //Set directions for Digitizer if an external trigger arrives
    //Again before change this check the manual please
    //The 2nd variable can also be set to one of the following 4 variables
    //CAEN_DGTZ_TRGMODE_DISABLED or CAEN_DGTZ_TRGMODE_EXTOUT_ONLY
    //CAEN_DGTZ_TRGMODE_ACQ_ONLY or CAEN_DGTZ_TRGMODE_ACQ_AND_EXTOUT
    conf |= CAEN_DGTZ_SetExtTriggerInputMode(handle[b], CAEN_DGTZ_TRGMODE_ACQ_ONLY);

    //Set the mode for a software trigger 
    conf |= CAEN_DGTZ_SetSWTriggerMode(handle[b],CAEN_DGTZ_TRGMODE_ACQ_ONLY);
    
    //Set Self-Trigger Threshold
    conf = CAEN_DGTZ_SetChannelTriggerThreshold(handle[b],0,32768);

    //Set trigger to be ACQ_ONLY
    conf = CAEN_DGTZ_SetChannelSelfTrigger(handle[b],CAEN_DGTZ_TRGMODE_ACQ_ONLY,1);

    //Set up the Digitizer's channels
    conf |= CAEN_DGTZ_SetChannelEnableMask(handle[b], 1);
    
    //Set number of events that the Digitizer needs to accumulate in the board memory
    //before being available for the readout process
    conf |= CAEN_DGTZ_SetDPPEventAggregation(handle[b], Params[0].EventAggr, 0);
    
    //Set up DPP parametes for the specific channelMask
    conf |= CAEN_DGTZ_SetDPPParameters(handle[b], Params[0].ChannelMask, &DPPParams);
    
    //Set Max number of events that can be read per DAQ
    conf = CAEN_DGTZ_SetMaxNumEventsBLT(handle[b],30);

    //Prepare setting for each channel
    for (i=0; i<MaxNChannels; i++){
      if (Params[0].ChannelMask & (1<<i)){
	
	//Set DC offset to the input signal to adapt it 
	//to the digitizer's dynamic range 
	conf |= CAEN_DGTZ_SetChannelDCOffset(handle[b], i, 0x8000);
	
	//Set the pre-trigger size
	conf |= CAEN_DGTZ_SetDPPPreTriggerSize(handle[b], i, 80);
	
	//Set the polarity for the given channel 
	//on top of the handle thw function needs to be fed with either
	//CAEN_DGTZ_PulsePolarityPositive or CAEN_DGTZ_PulsePolarityNegative
	conf |= CAEN_DGTZ_SetChannelPulsePolarity(handle[b], i, Params[0].PulsePolarity);
      }      
    }
  }
  //check whether there has been a problem during the setting up or not
  if (conf) {
    cout << "Warning : errors found during the programming of the Digitizer " << GetName() << " some setting may not be executed" << endl;
    return conf;
  } else {
    cout << "Digitizer "<< GetName() << " has been programed correctly" << endl;
    return 0;
  }
}

//This function makes use of the Digitizer API functions to set up
//The Digitizer's initial configuration for a basic read out process
//And the data acquisition phase, no DPP parameters present
//Connected only via OpticalLink 
int
ScopeCAENv1720::BasicDigiConfig(int* handle, DigitizerParams_t* Params, CAEN_DGTZ_DPP_PHA_Params_t* DPPParams){

  //Define Variables
  int conf=0;

  //Loop through boards
  for (int g=0; g<MaxNb; g++){
    
    conf = CAEN_DGTZ_Reset(handle[g]);
    conf = CAEN_DGTZ_SetRecordLength(handle[g],4096);
    conf = CAEN_DGTZ_SetChannelEnableMask(handle[g],1);
    conf = CAEN_DGTZ_SetChannelTriggerThreshold(handle[g],0,32768);
    conf = CAEN_DGTZ_SetChannelSelfTrigger(handle[g],CAEN_DGTZ_TRGMODE_ACQ_ONLY,1);
    conf = CAEN_DGTZ_SetSWTriggerMode(handle[g],CAEN_DGTZ_TRGMODE_ACQ_ONLY);
    conf = CAEN_DGTZ_SetMaxNumEventsBLT(handle[g],500);
    conf = CAEN_DGTZ_SetAcquisitionMode(handle[g],CAEN_DGTZ_SW_CONTROLLED);

  }

  //check whether there has been a problem during the setting up or not
  if (conf) {
    cout << "Warning : errors found during the programming of the Digitizer " << GetName() << " some setting may not be executed" << endl;
    return conf;
  } else {
    cout << "Digitizer "<< GetName() << " has been programed correctly" << endl;
    return 0;
  }
}

//This function allows the user to allocate the memory amount 
//Before the read out, and it needs to be called only after 
//the function that set all the parameters DigiConfig()
int
ScopeCAENv1720::AllocateMemoryAmount(int* handle){
  
  //Define varible
  int malloc = 0;
  char *buffer = NULL;
  uint32_t Size;

  //Loop through boards
  for (int b=0; b<MaxNb; b++) {
    //Calling function from CAEN library to allocate memory amount
    malloc = CAEN_DGTZ_MallocReadoutBuffer(handle[b],&buffer,&Size);
  }
  //Make sure we run the function properly
  if (malloc){
    cout << "Can't allocate memory buffers" << endl;
  }else{
    cout << "Allocating Memory Amount for : " << GetName() << endl;
  }
  return 0;
}

//This function allows the user to get a basic readout
//of the digitizer data, this doesn't store our save any
//Histograms but prints out on the screen a series of values
void
ScopeCAENv1720::BasicReadout(int* handle, DigitizerParams_t* Params, CAEN_DGTZ_DPP_PHA_Params_t* DPPParams){

  //Define variables
  CAEN_DGTZ_EventInfo_t eventInfo;
  void *Evt=NULL;
  char *buffer = NULL;
  char *evtptr = NULL;
  uint32_t size,bsize;
  uint32_t numEvents;
  int NEVT=0;
  int count=0;
  int READOUT=0;
  int acqRun=0;
  int Quit=0;

  //Define Histograms ???

  //Prepare Digitizer for Data acquisition
  //By setting up all the function by using DigiConfig
  BasicDigiConfig(handle, Params, DPPParams);

  //Loop through Boards
  for (int g=0; g<MaxNb; g++){

    //Allocate Memory buffers
    //for the basic case we don't you the specific function to do this  
    READOUT = CAEN_DGTZ_MallocReadoutBuffer(handle[g],&buffer,&size);
  }

  //-------------------------//
  //Start Acquisition Process//
  //-------------------------//
  
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
      if (c=='q') {Quit=1;}
      
      //Send SW Trigger
      if (c=='t') {
	for (int b=0; b<MaxNb; b++) {
	  CAEN_DGTZ_SetSWTriggerMode(handle[b],CAEN_DGTZ_TRGMODE_ACQ_ONLY);
	  CAEN_DGTZ_SendSWtrigger(handle[b]);
	  printf("Send Software Trigger\n");
	}
      }
      
      //Restart 
      if (c=='r'){
	for (int a=0; a<MaxNb; a++){
	  CAEN_DGTZ_SWStopAcquisition(handle[a]);
	  printf("Restarted\n");
	  CAEN_DGTZ_ClearData(handle[a]);
	  CAEN_DGTZ_SWStartAcquisition(handle[a]);
	}
      }

      //Start Acquisition
      if (c=='s'){
	for (int d=0; d<MaxNb; d++){
	  CAEN_DGTZ_SWStartAcquisition(handle[d]);
	  printf("Acquisition Started for Board %d\n", d);
	}
	acqRun=1;
      }

      //Stop Acquisition
      if (c=='k'){
	for (int e=0; e<MaxNb; e++){
	  CAEN_DGTZ_SWStopAcquisition(handle[e]);
	  printf("Acquisition Stopped for Board %d\n", e);
	}
	acqRun=0;
      }
    }
    if (!acqRun) {
      Sleep(10);
      continue;
    }
  }
  
  //Looping through boards for readout data
  for (int h=0; h<MaxNb; h++){
    
    //Digitizer buffer readout
    READOUT = CAEN_DGTZ_ReadData(handle[h], CAEN_DGTZ_SLAVE_TERMINATED_READOUT_MBLT, buffer, &bsize);
    
    if (READOUT) {
      cout << "Can't readout" << endl;
    }

    //extract the number of events from the read buffer
    READOUT = CAEN_DGTZ_GetNumEvents(handle[h], buffer, bsize, &numEvents);
    
    //Store n event in the array
    count += numEvents;
    NEVT = count;
    
    //Loop through events
    for (int i=0; i<NEVT; i++){
      
      //Extract informations from the event
      READOUT = CAEN_DGTZ_GetEventInfo(handle[i], buffer, bsize, i, &eventInfo, &evtptr);
      
      //Decode the event and read data
      READOUT = CAEN_DGTZ_DecodeEvent(handle[i], evtptr, &Evt);
      READOUT = CAEN_DGTZ_FreeEvent(handle[i], &Evt);
    }
    
    printf("\nBoard %d: Retrieved %d Events\n", h, count);
 
  }
}

//This function allows to perform a very basic data acquisition test
//What it does is sending 10 SW trigger and read the event info
//Printing on the screen the respective trigger time tag
//Before running check the setting, mainly the connection type
void
ScopeCAENv1720::BasicReadout2(){
  
  //Define Same Variables
  CAEN_DGTZ_ErrorCode read;
  int handle[MaxNb];
  CAEN_DGTZ_BoardInfo_t BoardInfo;
  CAEN_DGTZ_EventInfo_t eventInfo;
  void *Evt = NULL;
  char *buffer = NULL;
  int i,b,ev;
  int count[MaxNb];
  char * evtptr = NULL;
  uint32_t size,bsize;
  uint32_t numEvents;
  i = sizeof(CAEN_DGTZ_TriggerMode_t);
  
  //Loop Through Boards
  //Open the Digitizer and print info on the screen
  //Set up the Digitizer for the Acquisition process
  for(b=0; b<MaxNb; b++){
    read = CAEN_DGTZ_OpenDigitizer(CAEN_DGTZ_PCI_OpticalLink,0,0,0x11110000,&handle[b]);
    if(read != CAEN_DGTZ_Success) {
      printf("Can't open digitizer\n");
    }
    read = CAEN_DGTZ_GetInfo(handle[b], &BoardInfo);
    printf("\nConnected to CAEN Digitizer Model %s, recognized as board %d\n", BoardInfo.ModelName, b);
    printf("\tROC FPGA Release is %s\n", BoardInfo.ROC_FirmwareRel);
    printf("\tAMC FPGA Release is %s\n", BoardInfo.AMC_FirmwareRel);
    printf("\n");
    
    read = CAEN_DGTZ_Reset(handle[b]);                                               /* Reset Digitizer */
    read = CAEN_DGTZ_GetInfo(handle[b], &BoardInfo);                                 /* Get Board Info */
    read = CAEN_DGTZ_SetRecordLength(handle[b],4096);                                /* Set the lenght of each waveform (in samples) */
    read = CAEN_DGTZ_SetChannelEnableMask(handle[b],1);                              /* Enable channel 0 */
    read = CAEN_DGTZ_SetChannelTriggerThreshold(handle[b],0,32768);                  /* Set selfTrigger threshold */
    read = CAEN_DGTZ_SetChannelSelfTrigger(handle[b],CAEN_DGTZ_TRGMODE_ACQ_ONLY,1);  /* Set trigger on channel 0 to be ACQ_ONLY */
    read = CAEN_DGTZ_SetSWTriggerMode(handle[b],CAEN_DGTZ_TRGMODE_ACQ_ONLY);         /* Set the behaviour when a SW tirgger arrives */
    read = CAEN_DGTZ_SetMaxNumEventsBLT(handle[b],3);                                /* Set the max number of events to transfer in a sigle readout */
    read = CAEN_DGTZ_SetAcquisitionMode(handle[b],CAEN_DGTZ_SW_CONTROLLED);          /* Set the acquisition mode */
    if(read != CAEN_DGTZ_Success) {
      printf("Errors during Digitizer Configuration.\n");
    }
  }
  
  //Allocate Buffer Memory
  for(b=0; b<MaxNb; b++){
    read = CAEN_DGTZ_MallocReadoutBuffer(handle[b],&buffer,&size);
  }
  
  //Start Basic Acquisition
  for(b=0; b<MaxNb; b++)
    read = CAEN_DGTZ_SWStartAcquisition(handle[b]);
  
  // Start acquisition loop
  // Send 10 Software triggers and print their time tag
  for (ev=0; ev<11; ev++) {
    //Loop through the Boards
    for(b=0; b<MaxNb; b++) {
      //Send Software Trigger
      read = CAEN_DGTZ_SendSWtrigger(handle[b]);
      //Read data collected with the sw trigger
      read = CAEN_DGTZ_ReadData(handle[b],CAEN_DGTZ_SLAVE_TERMINATED_READOUT_MBLT,buffer,&bsize); /* Read the buffer from the digitizer */
      //Extract the number of events
      read = CAEN_DGTZ_GetNumEvents(handle[b],buffer,bsize,&numEvents);
      count[b] +=numEvents;
      //Loop through events
      for (i=0;i<numEvents;i++) {
	//Get event infos
        read = CAEN_DGTZ_GetEventInfo(handle[b],buffer,bsize,i,&eventInfo,&evtptr);
	//Decode the info before extracted
        read = CAEN_DGTZ_DecodeEvent(handle[b],evtptr,&Evt);
	
	//Print Event Info
        cout << "Event Counter : " << eventInfo.EventCounter << endl;
        cout << "Trigger Time Tag : " << eventInfo.TriggerTimeTag << endl;
	
	//Free event
        read = CAEN_DGTZ_FreeEvent(handle[b],&Evt);
      } // end loop over events
    } // end of loop on boards
  } // end of readout loop
  
  //Free Buffer
  read = CAEN_DGTZ_FreeReadoutBuffer(&buffer);
  //Print out Number of Events
  //And close the digitizer
  for (b=0; b<MaxNb; b++){
    read = CAEN_DGTZ_CloseDigitizer(handle[b]);
  }
  
}
