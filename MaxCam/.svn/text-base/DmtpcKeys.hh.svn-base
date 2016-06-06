#ifndef DMTPC_KEYS_HH
#define DMTPC_KEYS_HH

#include "TClonesArray.h"
#include "TTree.h"
#include "TObject.h"
#include "TFile.h"
#include <vector>

using namespace std;

/** DEPRECATED! 
    \deprecated Former class for managing data processing.
*/

class DmtpcKeys : public TObject {

public:
      
      DmtpcKeys(TString keyfilename, 
		TString rawdatafiles,
		TString key, TString newdir, 
		int nreqmod, TString* reqmod, bool& pass);
      ~DmtpcKeys();  

      int getNKeys() {return _keys.size();}
      int getNFiles() {return _files.size();}
      TString getKey(int i) {return _keys[i];};
      TString getDir(int i) {return _keydir[i];}
      TString getFile(int i) {return _files[i];}
      vector<TString> getListOfKeys() {return _keys;}
      vector<TString> getListOfDirs() {return _keydir;}
      vector<TString> getListOfFiles() {return _files;}
      TString getRootDirName() {return _rootdirname;}
      bool contains(TString testkey);
      bool contains(int n, TString* testkeys);
      void writeKey();
      void addFriends(TTree* tree, int f);
      bool openOutputFile(TFile* routfile, int i);
      bool isSkim() {return _isSkim;}
      TTree* getBaseTree(int f);
      
      

private:    

      vector< TString > _keys;
      vector< TString > _keydir;
      vector< TString > _files;
      TString _rootdirname;
      TString _keyfilename;
      TString _key;
      TString _keynewdir;
      bool _isSkim;
      int _skimindex;
      
      ClassDef(DmtpcKeys,1)
        
};

#endif
