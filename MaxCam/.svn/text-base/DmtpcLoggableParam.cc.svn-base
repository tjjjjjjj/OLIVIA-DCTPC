#include "DmtpcLoggableParam.hh"
#include "TROOT.h"

#include <iostream>
using std::cout;
using std::endl;

ClassImp(DmtpcLoggableParam)

// ctor
DmtpcLoggableParam::DmtpcLoggableParam(TString paramName, TString paramDesc) : TObject() {
    //_value(0.0), _timeStamp(0), _paramName(paramName), _paramDesc(paramDesc) { 
  _value = 0.0;
  _timeStamp = 0;
  _paramName = paramName;
  _paramDesc = paramDesc;
}

// copy ctor
DmtpcLoggableParam::DmtpcLoggableParam(const DmtpcLoggableParam &other) : TObject(other) {
  _value     = other._value;
  if (other._timeStamp) {
    _timeStamp = new TTimeStamp();
    *_timeStamp = *(other._timeStamp);
  } else {
    _timeStamp = 0;
  }
  _paramName = other._paramName;
  _paramDesc = other._paramDesc;
}

// Empty default ctor
DmtpcLoggableParam::DmtpcLoggableParam() { 
  _value = 0.0;
  _timeStamp = 0;
  _paramName = "";
  _paramDesc = "";
}

DmtpcLoggableParam::~DmtpcLoggableParam() { 
  cout << "DmtpcLoggableParam dtor not yet done" << endl;
}

void 
DmtpcLoggableParam::print() {
  cout << "name, title, _value = [" << GetParamName() << "], [" << GetParamDesc() << "], " << GetValue() << endl;
  cout << _timeStamp->AsString() << endl;
}
