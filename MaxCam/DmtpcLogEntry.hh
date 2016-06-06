#ifndef DMTPC_LOG_ENTRY_HH
#define DMTPC_LOG_ENTRY_HH

#include "TROOT.h"
#include "TObject.h"
#include "TTimeStamp.h"
#include "TClonesArray.h"
#include "DmtpcLoggableParam.hh"

#include <iostream>
#include <vector>
using namespace std;

class TString;

/** A class to log environmental data (an alternative to MySQL) 
    
    Allows for different parameters to be logged at different intervals.
    Each log entry is a TClonesArray of 
    have an arbitrary number of entries.

    The data is all saved into a tree and so all DmtpcLogEntry
    instances must be logged simultaneously.  To solve this, each 
    DmtpcLogEntry holds a TClonesArray of DmtpcLoggableParam 
    allowing each TTree->Fill() to store an arbitrary number of 
    samples of a particular parameter.

    For example, if you are recording pressure every 1 second and
    temperature every 5 seconds and the TTree is filled every 10
    seconds, then a single fill will contain 10 pressure readings and
    2 temperature readings.
*/
class DmtpcLogEntry : public TObject {
    
public:
    
  DmtpcLogEntry();

  /** 
   * Constructor
   *
   * @param[in] paramName Name of the logged parameter
   * @param[in] paramDesc Description of the logged parameter
   */
  DmtpcLogEntry(TString paramName, TString paramDesc);

  /** 
   * Copy Constructor
   *
   */
  DmtpcLogEntry(const DmtpcLogEntry &other);
  virtual ~DmtpcLogEntry();

  /** 
   * Append a logable parameter to this log entry by using new with
   * placement feature of a TClonesArray
   *
   * @param[in]  dlp Pointer to a sample to be logged
   * @return        void
   */

  virtual void Append(DmtpcLoggableParam *dlp) { 
    Int_t lastid = _samples->GetEntries();
    new ( (*_samples)[lastid] ) DmtpcLoggableParam(*dlp);
  }
  virtual void SetParamName(TString pname) { _paramName=pname; }
  virtual void SetParamDesc(TString pname) { _paramName=pname; }

  virtual Float_t GetValue(Int_t ii) { return ((DmtpcLoggableParam*)_samples->At(ii))->GetValue(); }

  virtual TTimeStamp* GetTimeStamp(Int_t ii) { return ((DmtpcLoggableParam*)_samples->At(ii))->GetTimeStamp(); }
  virtual TString GetParamName() { return _paramName; }
  virtual TString GetParamDesc() { return _paramDesc; }

  virtual Int_t GetNSamples() { return _samples->GetEntries(); }
  TClonesArray* GetSamples() { return _samples; }

  /** 
   * Clear all samples in the DmtpcLogEntry 
   */
  virtual void Reset() {  _samples->Clear();  }

  
  /** 
   * Print the contents of this class in an easy-to-read way
   */
  void print();

protected:
    
private:
  TClonesArray *_samples;
  TString _paramName;     // name of parameter being logged (no spaces)
  TString _paramDesc;     // description of this parameter (spaces ok)

  ClassDef(DmtpcLogEntry,1)

};
#endif

