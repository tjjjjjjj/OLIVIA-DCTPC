#ifndef LOGGABLE_PARAMETER_HH
#define LOGGABLE_PARAMETER_HH

#include "TROOT.h"
#include "TObject.h"
#include "TTimeStamp.h"

class TString;

class LoggableParam : public TObject {
    
public:
    
  LoggableParam();
  LoggableParam(TString paramName, TString paramDesc);
  LoggableParam(const LoggableParam &other);
  virtual ~LoggableParam();

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

  ClassDef(LoggableParam,1)

};
#endif

