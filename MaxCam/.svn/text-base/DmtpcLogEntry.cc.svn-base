#include "DmtpcLogEntry.hh"
#include "DmtpcLoggableParam.hh"
#include "TROOT.h"

#include <iostream>
using std::cout;
using std::endl;

ClassImp(DmtpcLogEntry)

// ctor
DmtpcLogEntry::DmtpcLogEntry(TString paramName, TString paramDesc) : TObject() {
  cout << "ctor" << endl;
  _samples = new TClonesArray("DmtpcLoggableParam");
  _paramName = paramName;
  _paramDesc = paramDesc;
}

// copy ctor
DmtpcLogEntry::DmtpcLogEntry(const DmtpcLogEntry &other) : TObject(other) {
  cout << "copy ctor" << endl;
  if (other._samples) {
    _samples = (TClonesArray*)other._samples->Clone();
  } else {
    _samples = 0;
  }
  _paramName = other._paramName;
  _paramDesc = other._paramDesc;
}

// Empty default ctor
DmtpcLogEntry::DmtpcLogEntry() { 
  cout << "default ctor" << endl;
  _samples = 0;
  _paramName = "";
  _paramDesc = "";
}

DmtpcLogEntry::~DmtpcLogEntry() { 
  cout << "DmtpcLogEntry dtor not yet done" << endl;
}

void 
DmtpcLogEntry::print() {
  cout << "name, title, _value = [" << GetParamName() << "], [" << GetParamDesc() << "], " << endl;
  for (Int_t ii=0; ii<_samples->GetEntries(); ii++) {
    cout << ((DmtpcLoggableParam*)_samples->At(ii))->GetTimeStamp()->AsString("s") << " --> "
	 << ((DmtpcLoggableParam*)_samples->At(ii))->GetValue() << endl;
  }
}
