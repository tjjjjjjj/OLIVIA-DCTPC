//ScopeData.cc
//
#include "ScopeData.hh"

#include <iostream>
#include <vector>
//class ScopeWaveformData;
//class ScopeDataChannel;

using std::cout;
using std::endl;
using std::vector;

ScopeData::~ScopeData() {
  cout << GetName() << ": destructor" << endl;

  /*
  // this may cause a problem if multiple ScopeData instances
  // are created with the same ScopeConfig.  If so, then 
  // the same memory location will be "delete"d multiple times
  // which is likely a fatal error.
  delete _scopeConfig;

  // delete the data channels
  for (vector<ScopeDataChannel*>::iterator ii=_channels.begin();
       ii != _channels.end(); ++ii) {
    delete *ii;
  }
  */
}
ScopeData::ScopeData() {
  cout << GetName() << ": constructor" << endl;
  setNValidTriggers(0);
}

ScopeData::ScopeData(ScopeConfig* sc) {
  // create a ScopeData set and initialize the tree
  // based on the ScopeConfig class
  cout << GetName() << ": constructor" << endl;
  setNValidTriggers(0);
  setScopeConfig(sc);
  initDataset();
}



int
ScopeData::addDataChannel(ScopeDataChannel* dataChan) {
  cout << GetName() << ".addDataChannel()" << endl;
  _channels.push_back( dataChan );
  // return the index into the vector that can be used
  // to access the newly added data channel
  return _channels.size()-1;
}

ScopeDataChannel* 
ScopeData::dataChan(int id) { return _channels.at(id); }

ScopeConfig*
ScopeData::getScopeConfig() { return _scopeConfig; }

void 
ScopeData::setScopeConfig(ScopeConfig* sc) { _scopeConfig = sc; }

ScopeConfig* 
ScopeData::config() { return _scopeConfig; }

int
ScopeData::getNValidTriggers() { return _nValidTriggers; }
void 
ScopeData::setNValidTriggers(int nval) { _nValidTriggers = nval; }

vector<ScopeWaveformData*>
ScopeData::getWaveforms() {

  // Loop over all ScopeDataChannels to generate
  // a 1-D vector of all ScopeWaveformData* from the last event

  vector<ScopeWaveformData*> swfs;

  int nwaveforms = 0;
  for (vector<ScopeDataChannel*>::iterator ii=_channels.begin();
       ii != _channels.end(); ii++) {
    for (vector<ScopeWaveformData*>::iterator jj=(*((*ii)->wfs())).begin();
	 jj != (*((*ii)->wfs())).end(); jj++) {
    
      swfs.push_back( *jj );
      nwaveforms++;
    }
  }
  //cout << "nwaveforms = " << nwaveforms << endl;
  //cout << "nvalidtriggers = " << getNValidTriggers() << endl;

  return swfs;
}

void
ScopeData::clearData() {
  // loop over scope data channels
  for (vector<ScopeDataChannel*>::iterator ii=_channels.begin(); 
       ii != _channels.end(); ii++) {
    // loop over scope waveforms in each channel
    for (vector<ScopeWaveformData*>::iterator jj=(*ii)->wfs()->begin(); jj != (*ii)->wfs()->end(); jj++) {
      // delete the ScopeWaveformData objects
      //gROOT->Delete((*jj)->GetName());
      delete *jj;
    }
    // then clear out the entries from the vector 
    // (the entries are the pointers)
    (*ii)->wfs()->clear();
    (*ii)->wftimes()->clear();
  }
  // and reset the number of triggers
  setNValidTriggers(0);
}

//
// PRIVATE METHODS
//
void
ScopeData::initDataset() {
  cout << GetName() << ".initDataset()" << endl;

  ScopeConfig* sc = getScopeConfig();

  // loop over each channel and add a channel to the data set
  vector<ScopeChannelConfig*> channels = sc->chans();
  int jj = 0;
  int scopeChannelNumber;
  for (vector<ScopeChannelConfig*>::iterator ii=channels.begin();
       ii != channels.end(); ++ii) {
    scopeChannelNumber = (int)sc->board()->getChanNum(jj);
    cout << "sc->board()->getNRecordsPerBuffer() = " 
	 <<  sc->board()->getNRecordsPerBuffer() << endl;
    cout << "sc->board()->getBytesPerRecord()    = " 
	 <<  sc->board()->getBytesPerRecord() << endl;
    addDataChannel( new ScopeDataChannel(scopeChannelNumber, 
					 (*ii)->label(),
                                         sc->board()->getBytesPerRecord()) );
    jj++;
  }
}

char* 
ScopeData::GetName() { return "ScopeData"; }

