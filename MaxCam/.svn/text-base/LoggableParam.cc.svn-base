#include "LoggableParam.hh"
#include "TROOT.h"

#include <iostream>
using std::cout;
using std::endl;

ClassImp(LoggableParam)

// ctor
LoggableParam::LoggableParam(TString paramName, TString paramDesc) : TObject(), 
  _value(0.0), _timeStamp(0), _paramName(paramName), _paramDesc(paramDesc) { 
}

// copy ctor
LoggableParam::LoggableParam(const LoggableParam &other) : TObject(other) {
  _value     = other._value;
  _timeStamp = new TTimeStamp;
  *_timeStamp = *other._timeStamp;
  _paramName = other._paramName;
  _paramDesc = other._paramDesc;
}

// Empty default ctor
LoggableParam::LoggableParam() { }
LoggableParam::~LoggableParam() { }

void 
LoggableParam::print() {
  cout << "name, title, _value = [" << GetParamName() << "], [" << GetParamDesc() << "], " << GetValue() << endl;
  cout << _timeStamp->AsString() << endl;
}
