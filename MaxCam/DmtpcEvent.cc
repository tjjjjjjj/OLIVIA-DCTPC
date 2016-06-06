#include "DmtpcEvent.hh"
#include "MaxCamImage.hh"
#include "TDatime.h"
#include "TTimeStamp.h"
#include "MaxCamConfig.hh"
#include "ScopeDataInfo.hh"
#include "McDarkTrack.hh"
#include "McDarkDigi.hh"
#include "DmtpcDataConverter.hh"
#include "TMath.h"
#include "TH2F.h"
#include "TH1F.h"
#include <strings.h>
#include <assert.h>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;




ClassImp(DmtpcEvent)

// needs to at least accommodate ScopeAlazarATS860::N_RECORDS_PER_BUFFER
const int DmtpcEvent::MAX_N_TRIGGERS      = 250;
const int DmtpcEvent::MAX_N_SCOPE_CONFIGS = 250;  
const int DmtpcEvent::MAX_N_CONFIG_PARAMETERS = 10;
//const int DmtpcEvent::MAX_N_TRIGGERS      = 1;
//const int DmtpcEvent::MAX_N_SCOPE_CONFIGS = 1;  
const int DmtpcEvent::MAX_N_MCDARK_TRACKS = 10000;
const int DmtpcEvent::MAX_N_MCDARK_CCD_DIGIS = 1000000;
//const int DmtpcEvent::MAX_N_MCDARK_CCD_DIGIS = 1;
//const int DmtpcEvent::MAX_N_MCDARK_TRACKS = 1;


//____________________
//
// Class that stores event information.
//
DmtpcEvent::DmtpcEvent()  {
    
    // Event info initialization
    _eventNumber=0;
    _runNumber=0;
    _dtime = new TDatime;
    _tstime = new TTimeStamp;
    isCleared = false; 
    // Configuration parameters and environmental monitoring
    _experConfig = new TClonesArray("MaxCamChannel", MAX_N_CONFIG_PARAMETERS);    
    
    // CCD
    _ccdData   = new TClonesArray("TH2S", 2,1);
    _overscan   = new TClonesArray("TH2S", 2,1);

    _ccdDataClonesCache = new TClonesArray("TH2F",2,1); 
    _overscanClonesCache = new TClonesArray("TH2F",2,1); 

    _ccdConfig = new TClonesArray("MaxCamConfig", 2,1);

    _cache_dirty = false; 
    // Scope
    //_scopeData     = new TClonesArray("ScopeWaveform", MAX_N_TRIGGERS); 
    _scopeData     = new TClonesArray("ScopeWaveformData", MAX_N_TRIGGERS); 
    _scopeDataClonesCache     = new TClonesArray("TH1F", MAX_N_TRIGGERS); 
    _scopeDataInfo = new TClonesArray("ScopeDataInfo", MAX_N_SCOPE_CONFIGS); 

    
    // MC Tracks
    _mcDarkTrackList = new TClonesArray("McDarkTrack", MAX_N_MCDARK_TRACKS);

    // MC Hits
    _mcDarkCcdDigiList = new TClonesArray("McDarkDigi", MAX_N_MCDARK_CCD_DIGIS);

}


bool DmtpcEvent::check_cache(unsigned int index, const char * which)
{


  if (!strcasecmp(which,"all"))
  {
    return check_cache(index,"ccd") && check_cache(index,"scope") && check_cache(index,"overscan");
  }

  if (!strcasecmp(which,"scope"))
  {
    return _scopeDataCache.size() > index && _scopeDataCache[index] != 0; 
  }

  if (!strcasecmp(which,"ccd"))
  {
    return _ccdDataCache.size() > index && _ccdDataCache[index] != 0; 
  }

  if (!strcasecmp(which,"overscan"))
  {
    return _overscanCache.size() > index && _overscanCache[index] != 0; 
  }

  return false;
}

void DmtpcEvent::make_all_cache(const char * which)
{
  if(!strcasecmp(which,"all")) 
  {
    make_all_cache("ccd"); 
    make_all_cache("scope"); 
    make_all_cache("overscan"); 
  }

  TClonesArray * arr = NULL;

  if (!strcasecmp(which,"ccd")) arr = _ccdData; 
  if (!strcasecmp(which,"scope")) arr = _scopeData; 
  if (!strcasecmp(which,"overscan")) arr = _overscan; 

  if (arr==NULL) return; 

  for (int u = 0; u < arr->GetEntries(); u++)
  {
    make_cache(u,which); 
  }
}



void DmtpcEvent::make_cache(unsigned int index, const char * which)
{
  _cache_dirty = true; 
  if (check_cache(index,which)) return; 

  //cout << "Making " << which << " cache for i=" << index  << endl; 

  if(!strcasecmp(which,"all")) 
  {
    make_cache(index,"ccd"); 
    make_cache(index,"scope"); 
    make_cache(index,"overscan"); 
  }


  if (!strcasecmp(which,"scope"))
  {
    //Bounds check
    if ((int)index >= _scopeData->GetEntries()) return; 

    //fill in with 0's if necessary...
    if (_scopeDataCache.size() <= index) 
    {
      for (unsigned int u = _scopeDataCache.size(); u <= index; u++)
      {
        _scopeDataCache.push_back(0);
      }
    }


    _scopeDataCache[index] = DmtpcDataConverter::scopeExpand((ScopeWaveformData*) _scopeData->At(index),
                                                             _scopeData->At(index)->GetName(),
                                                             (TH1F*) (*_scopeDataClonesCache)[index]); 
    _scopeDataCache[index]->SetDirectory(0);
  }

  if (!strcasecmp(which,"ccd"))
  {
    //Bounds check
    if ((int)index >= _ccdData->GetEntries()) return; 

    //fill in with 0's if necessary...
    if (_ccdDataCache.size() <= index) 
    {
      for (unsigned int u = _ccdDataCache.size(); u <= index; u++)
      {
        _ccdDataCache.push_back(0);
      }
    }


    _ccdDataCache[index] = DmtpcDataConverter::ccdExpand((TH2S*)_ccdData->At(index),
                                                         _ccdData->At(index)->GetName(),
                                                         (TH2F*) (*_ccdDataClonesCache)[index]); 
    _ccdDataCache[index]->SetDirectory(0);
    assert(_ccdDataCache[index] == (*_ccdDataClonesCache)[index]); 
  }

  if (!strcasecmp(which,"overscan"))
  {
    //Bounds check
    if ((int)index >= _overscan->GetEntries()) return; 

    //fill in with 0's if necessary...
    if (_overscanCache.size() <= index) 
    {
      for (unsigned int u = _overscanCache.size(); u <= index; u++)
      {
        _overscanCache.push_back(0);
      }
    }

   ////expand clones array if larger size here 
   //if (_overscanClonesCache->GetEntries() < _overscan->GetEntries())
   //  _overscanClonesCache->Expand(_overscan->GetEntries()); 

    _overscanCache[index] = DmtpcDataConverter::ccdExpand((TH2S*)_overscan->At(index),
                                                      _overscan->At(index)->GetName(),
                                                      (TH2F*) (*_overscanClonesCache)[index]); 
    _overscanCache[index]->SetDirectory(0);
  }
}

void DmtpcEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class DmtpcEvent.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      _ccdData->Streamer(R__b);
      _overscan->Streamer(R__b);
      _ccdConfig->Streamer(R__b);
      _scopeData->Streamer(R__b);
      _scopeDataInfo->Streamer(R__b);
      _experConfig->Streamer(R__b);
      R__b >> _dtime;
      R__b >> _tstime;
      R__b >> _runNumber;
      R__b >> _eventNumber;
      _mcDarkTrackList->Streamer(R__b);
      _mcDarkCcdDigiList->Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, DmtpcEvent::IsA());
      clearCache(); 
   } else {
      R__c = R__b.WriteVersion(DmtpcEvent::IsA(), kTRUE);
      TObject::Streamer(R__b);
      _ccdData->Streamer(R__b);
      _overscan->Streamer(R__b);
      _ccdConfig->Streamer(R__b);
      _scopeData->Streamer(R__b);
      _scopeDataInfo->Streamer(R__b);
      _experConfig->Streamer(R__b);
      R__b << _dtime;
      R__b << _tstime;
      R__b << _runNumber;
      R__b << _eventNumber;
      _mcDarkTrackList->Streamer(R__b);
      _mcDarkCcdDigiList->Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

DmtpcEvent::DmtpcEvent(const DmtpcEvent &other, bool copy_data) : TObject() {
  // Copy constructor

  _eventNumber=other._eventNumber;
  _runNumber=other._runNumber;
  _dtime=new TDatime(*other._dtime);
  _tstime=new TTimeStamp(*other._tstime);
  isCleared = false; 
  _cache_dirty = false; 

  _experConfig = new TClonesArray(*other._experConfig);    
  if (copy_data)
  {
    _ccdData   = new TClonesArray(*other._ccdData);
    _overscan   = new TClonesArray(*other._overscan);
    _scopeData     = new TClonesArray(*other._scopeData);
  }
  else
  {
    _scopeData     = new TClonesArray("ScopeWaveformData", MAX_N_TRIGGERS); 
    _ccdData   = new TClonesArray("TH2S", 2,1);
    _overscan   = new TClonesArray("TH2S", 2,1);
  }
  _scopeDataInfo = new TClonesArray(*other._scopeDataInfo);
  _ccdConfig = new TClonesArray(*other._ccdConfig);

  _mcDarkTrackList = new TClonesArray(*other._mcDarkTrackList);
  _mcDarkCcdDigiList = new TClonesArray(*other._mcDarkCcdDigiList);


  //Do not copy cache, just create anew
  _ccdDataClonesCache = new TClonesArray("TH2",2,1); 
  _overscanClonesCache = new TClonesArray("TH2",2,1); 
  _scopeDataClonesCache     = new TClonesArray("TH1", MAX_N_TRIGGERS); 

  
  //  other.Print();;
  //  cout << GetName() << "Copy constructor not done" << endl;
}

DmtpcEvent::~DmtpcEvent() {
  delete _dtime;
  delete _tstime;
//  if (isCleared) return; 
  _experConfig->SetOwner(true);
  _ccdData->SetOwner(true);
  _overscan->SetOwner(true);
  _ccdConfig->SetOwner(true);
  _scopeData->SetOwner(true);
  _scopeDataInfo->SetOwner(true);
  _mcDarkTrackList->SetOwner(true);
  _mcDarkCcdDigiList->SetOwner(true); 

 _ccdDataClonesCache->SetOwner(true);
 _scopeDataClonesCache->SetOwner(true);
 _overscanClonesCache->SetOwner(true);

  delete _ccdDataClonesCache;
  delete _scopeDataClonesCache;
  delete _overscanClonesCache;

  
  delete _experConfig;
  delete _ccdData;
  delete _overscan;
  delete _ccdConfig;
  delete _scopeData;
  delete _scopeDataInfo;
  delete _mcDarkTrackList;
  delete _mcDarkCcdDigiList; 

}

MaxCamChannel*
DmtpcEvent::experimentConfig(TString name) {
    // Experimental configuration
    
    return (MaxCamChannel*)_experConfig->FindObject(name);

}


void
DmtpcEvent::printExperimentConfig() {
    // Print environmental parameters for this event.

    for ( int i=0; i<experimentConfig()->GetEntries(); i++) {
        cout << experimentConfig(i) << endl;
    }
}


void
DmtpcEvent::printMcTracks() {
    // Print MC tracks for this event.

    for ( int i=0; i<mcTrack()->GetEntries(); i++) {
        mcTrack(i)->print();
    }
}


TClonesArray * DmtpcEvent::ccdData()
{
  if (_ccdData == 0) return 0; 

  if (_ccdData->GetEntries() >0 && !strcmp(_ccdData->At(0)->ClassName(),"TH2F")); 
  {
    return _ccdData; 
  }
  make_all_cache("ccd"); 
  return _ccdDataClonesCache; 
}

TClonesArray * DmtpcEvent::overscan()
{
  if (_overscan == 0) return 0; 
  if (_overscan->GetEntries() >0 && !strcmp(_overscan->At(0)->ClassName(),"TH2F")); 
  {
    return _overscan; 
  }
  make_all_cache("overscan"); 
  return _overscanClonesCache; 
}

TClonesArray * DmtpcEvent::scopeData()
{
  if (_scopeData == 0) return 0; 
  if (_scopeData->GetEntries() >0 && !strcmp(_scopeData->At(0)->ClassName(),"TH1F")); 
  {
    return _scopeData; 
  }
  make_all_cache("scope"); 
  return _scopeDataClonesCache; 
}

TH2F * DmtpcEvent::ccdData(int i)
{
  if (i < 0 || i >= _ccdData->GetEntries()) 
    return NULL;

  if (!strcmp(_ccdData->At(i)->ClassName(),"TH2F"))
  {
    return (TH2F*) _ccdData->At(i); 
  }

  make_cache(i,"ccd"); 
  return _ccdDataCache[i]; 
}

TH2F * DmtpcEvent::overscan(int i)
{
  if (i < 0 || i >= _overscan->GetEntries()) 
    return NULL;

  if (!strcmp(_overscan->At(i)->ClassName(),"TH2F"))
  {
    return (TH2F*)_overscan->At(i); 
  }

  make_cache(i,"overscan"); 
  return _overscanCache[i]; 
}

TH1F * DmtpcEvent::scopeData(int i)
{
  if (i < 0 || i >= _scopeData->GetEntries()) 
    return NULL;

  if (!strcmp(_scopeData->At(i)->ClassName(),"TH1F"))
  {
    return (TH1F*)_scopeData->At(i); 
  }

  make_cache(i,"scope"); 
  return _scopeDataCache[i]; 
}


void DmtpcEvent::clearCache()
{

  if (!_cache_dirty) return; 
  _ccdDataClonesCache->Delete(); 
  _ccdDataClonesCache = new TClonesArray("TH2F",2,1); 
  _scopeDataClonesCache->Delete(); 
  _scopeDataClonesCache = new TClonesArray("TH1F",2,1); 
  _overscanClonesCache->Delete(); 
  _overscanClonesCache = new TClonesArray("TH2F",2,1); 

   for (unsigned int i = 0; i < _overscanCache.size(); i++) _overscanCache[i] = 0; 
   for (unsigned int i = 0; i < _ccdDataCache.size(); i++) _ccdDataCache[i] = 0; 
   for (unsigned int i = 0; i < _scopeDataCache.size(); i++) _scopeDataCache[i] = 0; 
  // cout << "Clearing Cache" << endl; 
  _cache_dirty = false; 
}


TH1F* DmtpcEvent::scopeData(int ichan, int itrig) {

  int nwf = _scopeData->GetEntries();
  if (nwf == 0) return NULL;

  int ntrig = scopeDataInfo(0)->getNTriggers();
  //int nchan = nwf/ntrig;
  return scopeData(ichan*ntrig + itrig);
}

TH1F* DmtpcEvent::scopeData(int trigger, int board, int channel) {
  // Get waveform for ith triger, board, channel.
  // This assumes that order of waveforms is not defined so 
  // waveform info is deduced from the waveform name.
  // It requires looping over all waveforms in the event
  // until a right waveform is found.
  //
  
  trigger++; //to match convention trigger=1,2,3...
	
  int nwf=rawScopeData()->GetEntries();
  for (int iwf=0; iwf<nwf; iwf++) {
    ScopeWaveformData *hwf=rawScopeData(iwf);
    const char *hname = hwf->GetName();
    assert(hname);

    int boardID, triggerID;
    char channelStr;
    sscanf(hname, "scope_%d_%c_%d", &boardID, &channelStr, &triggerID);
    int channelID = (channelStr=='A') ? 0 : 1;
	  
    if (trigger==triggerID && board==boardID && channel==channelID) return scopeData(iwf);
  }

  // not found
  return 0;
}
