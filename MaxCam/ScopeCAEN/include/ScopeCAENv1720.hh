// ScopeCAENv1720.hh
//
// Pietro Giampa, RHUL 2013

#ifndef SCOPE_CAEN_V1720_HH
#define SCOPE_CAEN_V1720_HH

//Generic stuff
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

//CAEN and MaxCam stuff
#include "CAENDigitizer.h"
#include "CAENDigitizerType.h"
#include "CAENKeyComm.hh"

//ROOT stuff
#include "TH2.h"
#include "TClonesArray.h"
#include "TObject.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "TH1.h"
#include "TTree.h"
#include "TGraph.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TString.h"
#include "TFile.h"
#include "TTimeStamp.h"
#include "TDatime.h"
#include "TObjArray.h"

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

typedef struct
{
  int RecordL;
  int ESize;
  int Ecounter;
  int TrigTime;
  
} EvtInfo;

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

  //Keyboards command check
  int KeyCheck();
	
  //-----------------//
  //Readout Functions//
  //-----------------//

  //Basic readout function
  void readout(TString fFileName);

private:
	
  const static int MaxNChannels;
  const static int MaxBits;
  const static int MAXNB;
  
  int _h;
  void setScopeHandle(int h);
  int getScopeHandle();
  char* GetName();
  Int_t fNEvents;
  Int_t fRunNum;
  TString fDataPath;
  TString fFileName;

 };

#endif
