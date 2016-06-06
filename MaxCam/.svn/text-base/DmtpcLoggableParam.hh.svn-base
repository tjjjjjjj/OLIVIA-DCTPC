#ifndef DMTPC_LOGGABLE_PARAMETER_HH
#define DMTPC_LOGGABLE_PARAMETER_HH

#include "TROOT.h"
#include "TObject.h"
//#include "TString.h"
#include "TTimeStamp.h"

//class TTimeStamp;
class TString;

/** 
 * 
 * A base class to hold a timestamp and a value and other information
 * for "environmental" data that is to be logged.  All logged data
 * should be stored either in this class or a class that inherits from
 * this class.
 *
 * e.g. continual logging of temperature, pressure, turbo status,
 * camera status, power status, valve status, etc, etc
 * 
 * See also DmtpcLogEntry
 *
 */
class DmtpcLoggableParam : public TObject {
    
public:
    
  DmtpcLoggableParam();
  DmtpcLoggableParam(TString paramName, TString paramDesc);
  DmtpcLoggableParam(const DmtpcLoggableParam &other);
  virtual ~DmtpcLoggableParam();

  virtual void SetValue(Float_t newVal) { _value=newVal; }
  virtual void SetTimeStamp(TTimeStamp* newTime) { _timeStamp=newTime; }
  virtual void SetParamName(TString pname) { _paramName=pname; }
  virtual void SetParamDesc(TString pname) { _paramName=pname; }

  virtual Float_t GetValue() { return _value; }
  virtual TTimeStamp* GetTimeStamp() { return _timeStamp; }
  virtual TString GetParamName() { return _paramName; }
  virtual TString GetParamDesc() { return _paramDesc; }

  void print();

protected:
    
private:
  Float_t _value;         // value of parameter
  TTimeStamp* _timeStamp; // time when value was obtained
  TString _paramName;     // name of parameter being logged (no spaces)
  TString _paramDesc;     // description of this parameter (spaces ok)

  ClassDef(DmtpcLoggableParam,1)

};
#endif

