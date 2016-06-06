#ifndef SIMTOOLS_HH
#define SIMTOOLS_HH

#include "TDatime.h"
#include "TTimeStamp.h"
#include "TH2.h"
#include "TString.h"

namespace SimTools{

  double getTrackLength(TH2F* trackImage,double phi);
  TTimeStamp* makeTimeStamp(double t);
  TDatime* makeTDatime(double t);
  TDatime* convertStamptoDatime(TTimeStamp *stamp);
  double getTimeDouble(TTimeStamp *stamp);
  double getTimeDouble(TDatime *datime);
  double findParticleMass(TString name);
  double findParticleA(TString name);

};

#endif
