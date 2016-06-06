#ifndef RGASCAN_HH
#define RGASCAN_HH

#include "TVectorD.h"
using namespace std;

class RGAview {
private:
  char datafile[50];
  char rootfile[50];
  Int_t scan_number;
  Int_t nmasses;
  Int_t mass_number;
  void readheader(Double_t * massarray);
  void readdata(Double_t * massarray);
  void convert_to_TDatime(TDatime *ttime, TString timestring);
public:
  RGAview(char RGAdata[50], char filename[50]);
  RGAview(const RGAview & object);
  ~RGAview();
  virtual int readRGA();
  virtual void view_one_mass(int n);
  virtual void next_mass();
  virtual void view_one_scan(int n);
  virtual void next_scan();
  virtual void view_total_pressure();

ClassDef(RGAview, 1)
};
  
#endif  
