#ifndef DMTPC_PMT_QE_HH
#define DMTPC_PMT_QE_HH

#include "TGraph.h"
#include "TRandom3.h"

class DmtpcPMTQE : public TGraph {

public:

  // Ctors
  //DmtpcPMTQE() {};

  DmtpcPMTQE(TString filename="", double val=-1.0, TRandom3 *rand=0);

  DmtpcPMTQE(const DmtpcPMTQE &other);

  virtual ~DmtpcPMTQE();

  bool isDetected(float lambda);

private:
  TRandom3 *rnd;

  ClassDef(DmtpcPMTQE,1)
};

#endif

