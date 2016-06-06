//ScopeConfig.cc
//
#include <iostream>
#include "ScopeConfig.hh"

#include "ScopeBoardInfo.hh"
#include "ScopeClockConfig.hh"
#include "ScopeTriggerConfig.hh"
#include "ScopeChannelConfig.hh"
//#include <vector>

ClassImp(ScopeConfig)

using std::cout;
using std::endl;
using std::vector;


ScopeConfig::~ScopeConfig() {
  cout << GetName() << ": destructor" << endl;

  //delete _clockCfg;
  //delete _trgCfg;

  // delete the channel configurations
  for (vector<ScopeChannelConfig*>::iterator ii=_chanCfg.begin();
       ii != _chanCfg.end(); ++ii) {
    delete *ii;
  }
}

ScopeConfig::ScopeConfig() {}

ScopeConfig::ScopeConfig(int junk) {
  cout << GetName() << " junk constructor" << endl;
}

  // set up records, clock and trigger, but not channels 
  // (use addChannelConfig() to add channels)
ScopeConfig::ScopeConfig(ScopeBoardInfo* boardcfg,
			 ScopeClockConfig* clockcfg,
			 ScopeTriggerConfig* trgcfg) {
  setBoardConfig(boardcfg);
  setClockConfig(clockcfg);
  setTriggerConfig(trgcfg);
}

// set up a 2 channel scope
ScopeConfig::ScopeConfig(ScopeBoardInfo* boardcfg,
			 ScopeClockConfig* clockcfg,
			 ScopeTriggerConfig* trgcfg, 
			 ScopeChannelConfig* chancfg1, 
			 ScopeChannelConfig* chancfg2) {
  ScopeConfig(boardcfg, clockcfg, trgcfg);
  addChannelConfig(chancfg1);
  addChannelConfig(chancfg2);
}
// copy constructor
ScopeConfig::ScopeConfig(ScopeConfig& sc) {
  cout << GetName() << " copy constructor not done" << endl;
  //ScopeConfig(sc->board(), sc->clock(), sc->trig());
  //// also, copy any channels that are present
  //for (vector<ScopeChannelConfig*>::iterator ii=sc->chans().begin();
  //     ii != sc->chans().end(); ++ii) {
  //  addChannelConfig(*ii);
  //}
}

void ScopeConfig::setBoardConfig(ScopeBoardInfo* bi) { _boardCfg = bi; }
void ScopeConfig::setClockConfig(ScopeClockConfig* cfg) { _clockCfg = cfg; }
void ScopeConfig::setTriggerConfig(ScopeTriggerConfig* cfg) { _trgCfg = cfg; }
void ScopeConfig::addChannelConfig(ScopeChannelConfig* cfg) { _chanCfg.push_back(cfg); }

ScopeBoardInfo*     ScopeConfig::board() { return _boardCfg; }
ScopeClockConfig*   ScopeConfig::clock() { return _clockCfg; }
ScopeTriggerConfig* ScopeConfig::trig() { return _trgCfg; }
ScopeChannelConfig* ScopeConfig::chan(int ii) { return _chanCfg.at(ii); }

vector<ScopeChannelConfig*> ScopeConfig::chans() { return _chanCfg; }

void ScopeConfig::print() {  cout << GetName() << ".print()" << endl; }

void ScopeConfig::printall() {
  cout << "ScopeConfig:" << endl;
  cout << "  BoardInfo"  << endl;
  cout << "    to do" << endl;
  cout << "  ClockConfig" << endl;
  cout << "    rate       " << std::hex << clock()->getClockRate() << std::dec << endl;
  cout << "    source     " << clock()->getClockSource() << endl;
  cout << "    edge       " << clock()->getClockEdge() << endl;
  cout << "    decimation " << clock()->getClockDecimation() << endl;
  cout << "  TriggerConfig" << endl;
  cout << "    trigOp     " << trig()->trigOp() << endl;
  cout << "    trigEn1    " << trig()->trigEn1() << endl;
  cout << "    trigEn2    " << trig()->trigEn2() << endl;
  cout << "    trigSrc1   " << trig()->trigSrc1() << endl;
  cout << "    trigSrc2   " << trig()->trigSrc2() << endl;
  cout << "    trigSlope1 " << trig()->trigSlope1() << endl;
  cout << "    trigSlope2 " << trig()->trigSlope2() << endl;
  cout << "    trigLevel1 " << trig()->trigLevel1() << endl;
  cout << "    trigLevel2 " << trig()->trigLevel2() << endl;
  for (int ichan=0; ichan<chans().size(); ichan++) {
    cout << "  ChannelConfig (" << ichan << ")" << endl;
    cout << "    chanNum    " << chan(ichan)->getChannelNumber() << endl;
    cout << "    voltRange  " << std::hex << chan(ichan)->getVoltageRange() << std::dec << endl;
    cout << "    impedance  " << chan(ichan)->getInputImpedance() << endl;
    cout << "    coupling   " << chan(ichan)->getCoupling() << endl;
  }

}

// 
// PRIVATE
//
const char* ScopeConfig::GetName() const { return "ScopeConfig"; }

