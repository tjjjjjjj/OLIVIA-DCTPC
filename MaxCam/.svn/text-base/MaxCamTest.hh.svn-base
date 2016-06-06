#ifndef MAXCAM_TEST_HH
#define MAXCAM_TEST_HH
#include "TROOT.h"

#include "MaxCamTestBase.hh"

class MaxCam;
class MaxCamConfig;
class MaxCamChannel;
class TH2F;
class TH1F;
class TDatime;
class TTree;
class TString;
class TFile;
class MaxCamTrack;
#include <vector>
using std::vector;

class MaxCamTest : public MaxCamTestBase {

public:

  // Constructors
  MaxCamTest(int debugLevel=0, TString filePath="test.root", TString cameraType="apogee");

  virtual ~MaxCamTest() {};

  //virtual int begin();
  virtual int event();
  virtual int end();

  int saveRun(char *keyword, char *description, char *location="26-544", 
	      char *scpdest="mitdm00.mit.edu:/data/",
	      char *fname="dbaccess.txt");

  void writeValueToDB(TString table, double val);

  void writeCCDConfigToDB(const char *fname="dbaccess.txt");

  void recoverFromDischarge(UInt_t sleepTime=900000);

  void scanVoltage(float min=2.6, float max=2.86, float step=0.05, bool wire=true);

private:

  int *_eventCounter; // total exposures for this run, accepted events

  ClassDef(MaxCamTest,0)
};

#endif

