#ifndef MAXCAM_QE_HH
#define MAXCAM_QE_HH

#include "TGraph.h"
#include "TRandom.h"

class MaxCamQE : public TGraph {

public:

  // Ctors
  MaxCamQE() {};

  MaxCamQE(const char *fileName, TRandom *rand=0);

  MaxCamQE(const MaxCamQE &other);

  virtual ~MaxCamQE();

  bool isDetected(float lambda);

private:
  TRandom *rnd;

  ClassDef(MaxCamQE,1)
};

#endif

