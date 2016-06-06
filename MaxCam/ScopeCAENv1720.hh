// ScopeCAENv1720.hh
//
#ifndef SCOPE_CAEN_V1720_HH
#define SCOPE_CAEN_V1720_HH

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <CAENDigitizer.h>
#include <CAENDigitizerType.h>
#include <CAENComm.h>
#include "CAENKeyComm.hh"
#include "ScopeCAEN.hh"

//
// CAEN v1720 Class (PCI Digitizer card -- an oscilloscope for your computer)
// 
// 

//Define the variable type DigitizerParams_t
//That we are going to use to store some important setting
typedef struct
{
  CAEN_DGTZ_ConnectionType LinkType;
  uint32_t VMEBaseAddress;
  uint32_t RecordLength;
  uint32_t ChannelMask;
  int EventAggr;
  CAEN_DGTZ_PulsePolarity_t PulsePolarity;
  CAEN_DGTZ_DPP_AcqMode_t AcqMode;
  CAEN_DGTZ_IOLevel_t IOlev;
} DigitizerParams_t;

class ScopeCAENv1720 {
	
public:

  //---------//
  //Main Body//
  //---------//

  // Destructor
  ~ScopeCAENv1720();

  // Constructor
  ScopeCAENv1720();

  //---------//
  //Functions//
  //---------//

  //Reset Digitizer
  void ResetDigi(int* handle);

  //Configure Digitizer
  int DigiConfig(int* handle, DigitizerParams_t* Params, CAEN_DGTZ_DPP_PHA_Params_t* DPPParams);

  //Configure Digitizer
  int BasicDigiConfig(int* handle, DigitizerParams_t* Params, CAEN_DGTZ_DPP_PHA_Params_t* DPPParams);

  //Allocate Memory Buffer
  int AllocateMemoryAmount(int* handle);

  //Set acq parameters
  DigitizerParams_t* SetStandardParameters();

  //Set DPP parameters
  CAEN_DGTZ_DPP_PHA_Params_t* SetDPPParameters();

  //Keyboards command check
  int KeyCheck();

  //Close Digitizer
  int StopReadout(int* handle);

  //Extract handle from Digitizer
  int* GetHandle(DigitizerParams_t* Params);
	
  //-----------------//
  //Readout Functions//
  //-----------------//

  //Basic readout function
  void BasicReadout(int* handle, DigitizerParams_t* Params, CAEN_DGTZ_DPP_PHA_Params_t* DPPParams);

  //Basic readout function
  void BasicReadout2();

private:
	
  const static int MaxNChannels;
  const static int MaxBits;
  const static int MaxNb;
  
  int _h;
  void setScopeHandle(int h);
  int getScopeHandle();
  char* GetName();

 };

#endif
