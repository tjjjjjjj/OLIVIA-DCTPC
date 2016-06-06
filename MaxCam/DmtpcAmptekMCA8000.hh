#ifndef DMTPC_Amptek_MCA8000_HH
#define DMTPC_Amptek_MCA8000_HH

#include "TH1F.h"

/** 
 * Class to ....
 */

class DmtpcAmptekMCA8000 : public TObject {

public:
    
  //
  // Ctors
  DmtpcAmptekMCA8000();
  DmtpcAmptekMCA8000(const char* filename);
  
  void show() { _spec->Draw(); }
  TH1F* spec() { return _spec; }

private:  
  
  TString filename;
  TH1F* _spec;

  ClassDef(DmtpcAmptekMCA8000,1)
        
};

#endif
