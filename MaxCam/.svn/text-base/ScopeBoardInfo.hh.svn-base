// ScopeBoardInfo.hh
//
#ifndef SCOPE_BOARD_INFO_HH
#define SCOPE_BOARD_INFO_HH

//#include <string>
#include <vector>

#include "TObject.h"
#include "TString.h"

//#include "AlazarApi.h"
//#include "AlazarCmd.h"

#include "ScopeTypes.hh"

//
// Class to contain information about the scope board
//
class ScopeBoardInfo : public TObject {

public:

  // Destructor
  ~ScopeBoardInfo();
  // Constructors
  ScopeBoardInfo();

  //int  setBoardType(const std::string* bt);
  int  setBoardType(const TString* bt);
  int  setBoardChannelNumbers(const AU8 chans[], int nchans);
  void setRecordPreSize(AU32 rps);
  void setRecordPostSize(AU32 rps);
  void setBytesPerRecord(int bpr);
  void setNRecordsPerBuffer(int nrec);

  AU8  getChanNum(int ii);
  int  getNRecordsPerBuffer();
  AU32 getRecordSize();
  AU32 getRecordPreSize();
  AU32 getRecordPostSize();
  int  getBytesPerRecord();

private:
  int _bitsPerSample;
  //const std::string* _boardType;
  const TString* _boardType;
  std::vector<AU8> _chanNumbers;
  AU32 _recordPreSize, _recordPostSize;  // pre/post trigger
  int _bytesPerRecord, _nRecordsPerBuffer;

ClassDef(ScopeBoardInfo, 1)
};

#endif
