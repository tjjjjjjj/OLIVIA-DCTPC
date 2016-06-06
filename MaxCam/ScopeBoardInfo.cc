//ScopeBoardInfo.hh
//
#include <iostream>

#include "ScopeBoardInfo.hh"
//#include "Scope.hh"

using std::cout;
using std::endl;;
//using std::string;

ClassImp(ScopeBoardInfo)

ScopeBoardInfo::~ScopeBoardInfo() {}
ScopeBoardInfo::ScopeBoardInfo() {}

int
ScopeBoardInfo::setBoardChannelNumbers(const AU8 chans[], int nchans) {
  for (int ii=0; ii<nchans; ii++) {
    _chanNumbers.push_back(chans[ii]);
    cout << "channel number = " << (int)_chanNumbers.at(ii) << endl;
  }
  return 0;
}

int
//ScopeBoardInfo::setBoardType(const string* bt) {
ScopeBoardInfo::setBoardType(const TString* bt) {
  //if (!Scope::isValidType(bt)) {
  //  assert(!"ERROR: Invalid scope board type");
  //  return -1;
  //}
  _boardType = bt;
  return 0;
}

void ScopeBoardInfo::setRecordPreSize(AU32 rps) { _recordPreSize = rps; }
void ScopeBoardInfo::setRecordPostSize(AU32 rps) { _recordPostSize = rps; }
void ScopeBoardInfo::setBytesPerRecord(int bpr) { _bytesPerRecord = bpr; }
void ScopeBoardInfo::setNRecordsPerBuffer(int nrec) { _nRecordsPerBuffer = nrec; }

AU8  ScopeBoardInfo::getChanNum(int ii) { return _chanNumbers.at(ii); }
AU32 ScopeBoardInfo::getRecordSize() { return (_recordPreSize + _recordPostSize); }
AU32 ScopeBoardInfo::getRecordPreSize() { return _recordPreSize; }
AU32 ScopeBoardInfo::getRecordPostSize() { return _recordPostSize; }
int  ScopeBoardInfo::getBytesPerRecord() { return _bytesPerRecord; }
int  ScopeBoardInfo::getNRecordsPerBuffer() { return _nRecordsPerBuffer; }

//
// PRIVATE
//
