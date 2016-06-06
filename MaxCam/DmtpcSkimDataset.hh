#ifndef DMTPC_SKIM_DATASET_HH
#define DMTPC_SKIM_DATASET_HH

#include "TROOT.h"
#include "TFile.h"
#include "TObject.h"
#include "TH2F.h"
#include "TTree.h"
#include <vector>
#include <map>
#include <list>
#include "DmtpcSkimEvent.hh"
#include "DmtpcDataset.hh"
#include "DmtpcGainMap.hh"
#include "TMap.h" 
#include "TIterator.h" 
#include "waveformtools/include/WaveformVector.hh"
#include "waveformtools/include/WaveformTools.hh"

class TChain;
class TTree;
class TH2F;
class TObjString;


#define DEFAULT_DATA_DIR  "/net/zwicky/dmtpc/data/" 

using namespace std; 

/**
Class to hold reconstructed data trees
\author Jeremy Lopez (documentation)
*/
class DmtpcSkimDataset : public TObject {

public:
/**
Constructor
*/
    DmtpcSkimDataset(unsigned maxtrees = 3);
    
/**
Copy constructor
\param other DmtpcSkimDataset to copy
*/
    DmtpcSkimDataset(const DmtpcSkimDataset &other);
    
/**
Destructor
*/
    virtual ~DmtpcSkimDataset();
    
/**
Always returns "DmtpcSkimDataset"
*/
    virtual const char* GetName() const { return "DmtpcSkimDataset"; }
    
/**
\return tree holding reconstructed event info
*/
    TTree* tree(const char * key = "skim")  { return (TTree*)chain(key); }
    
/**
\return tree holding reconstructed event info as a TChain object
*/
    TChain* chain(const char * key = "skim") { return _indices.count(key) ? _trees[_indices[key]] : 0 ; }
    
/**
\return pointer to the file containing the trees
*/
    TFile* file() { return _file; }
    
/**
Create or open a root file for writing/reading
\param fname the name of the file
\param foption a file opening option.  If "create" will create a new file and open for writing. If "" will open for reading only.  Others will report an error and quit.  This will not overwrite an existing file.
*/
    void createRootFile( const char *fname, TString foption);
    
/**
Opens a new file for reading.  Just call createRootFile(fname,"")
\param fname the name of the file
*/
    void openRootFile( const char *fname); 
    
/**
If a file is open, deletes the event tree if it exists and closes the file.
*/
    void closeRootFile();
 
/**
Creates a new root file for writing.  Just calls createRootFile(fname,"create")
\param fname the name of the file
*/
    void newRootFile(const char *fname);
    
    const char* getFileName(); 
   
/**
 * \deprecated
  Does nothing right now. Here for legacy purposes.
*/
    void clearEventMemory();
    
/**
\return the reconstructed event
*/  
    DmtpcSkimEvent* event(const char * key = "skim"){ return _indices.count(key) ? _events[_indices[key]] : 0 ;} 

/** 
\return a pointer to the indices 
*/
    const std::map<std::string,  unsigned> * getTreeIndices() const { return &_indices; } 

/**
\return the original event or NULL if the DmtpcEvents are not loaded
*/
    DmtpcEvent* orig_event() { return orig_loaded ? d->event() : NULL;}

/**
Get an entry from the data tree
\param i the number of the entry to get
*/
    void getEvent(int i);
  
/**
 *\deprecated
  Find the entry corresponding to a particular event in the original data tree. This should 
  no longer be needed. 
  \param n the event to find (number from the original dataset)
  \param startingGuess an initial value to check in order to speed up the search
  \return the index of this event or the nearest one in the skim dataset
*/
    int getEventNumber(int n, int startingGuess=-1); 
    
/**
  \return the total number of events/entries in the tree
*/
    int nevents() { return _nevents; }

/**
  Load the bias frames from this run into memory
  \param load true to load false to not load
  \param fname the name of the file where the bias frames are held
*/
    void loadBiasFrames(bool load, const char * fname=""); 
    
    TH2F** getBiasFrames(){ return _biasFrames; };
    TH2F* getBiasFrame(int i){ return _biasFrames[i]; }
    DmtpcGainMap* getGainMap(int i){return (_gainMaps && _gainMaps->GetEntries() > i ) ? (DmtpcGainMap*)_gainMaps->At(i) : 0;}
    DmtpcGainMap* getGainMap(TString serialNumber); 
    TObjArray* getGainMaps() {return _gainMaps;}
    void addGainMap(DmtpcGainMap* gainMap) {_gainMaps->Add(gainMap);}
    void writeGainMaps(); 

/**
Merges the given trees into an a skim tree using DmtpcSkimEvent 
\param tmpskim a tree holding track information (except RBI info)
\param burnin a structure holding information about potential RBI events
\param sparkref a structure holding information about saturated pixels from sparks
\param run_n the run number
*/
   void mergeTempTrees(TTree * tmpskim,  list<vector<vector<vector<BurninEncoded_t> > > * > * burnin, 
                                         list<vector<pair<int,int> > *> * sparkref, int run_n); 


  void mergeTempTreesv5(TTree * tmpskim,  list<vector<vector<vector<BurninEncoded_t> > > * > * burnin, 
                                         list<vector<pair<int,int> > *> * sparkref, int run_n); 
            
/**
Loads the original events. Note that now instead of being stored in the skim file, this method will
transparently load the original file. . 
/param load whether or not to load the events. if they are loaded and load is false, the file will be closed. 
/param orig_dataset the filename of the original run file or NULL to try to figure it out
*/
   void loadDmtpcEvent(bool load = true, const char * orig_dataset = NULL);


void setConfig(const char * config) { if (_config){_config->Delete(); } _config = new TObjString(config); } 
void writeConfig() { _file->cd(); _config->Write("config"); gROOT->cd(); }

const char * getConfig(){ return _config->String().Data(); } 


void writeStitch(Dmtpc4ShooterStitcher * stitch); 
Dmtpc4ShooterStitcher * getStitch() { return (Dmtpc4ShooterStitcher *) _file->Get("stitch"); } 


/**
Loads MaxCamClusterImage objects containing lists of clusters and processed images.  Default is true, but set to false if not needed to speed up readout.
\param l true to load, false to not load
*/
   void loadClusters(bool l) { tree()->SetBranchStatus("_clusters",l); }
private:
    

   TChain ** _trees; 
   DmtpcSkimEvent ** _events; 
   unsigned ntrees; 
   unsigned _max;
   std::map<std::string, unsigned> _indices;  

   TFile* _file; ///< event file
   TH2F** _biasFrames;  
   int _nevents; ///<number of events

   TObjArray* _gainMaps;///<TObjArray of gainMaps   

   TObjString * _config;
   bool bias_loaded;    
   int nbias; 

   int current_index; 

   bool orig_loaded; 
   DmtpcDataset * d; 
  ClassDef(DmtpcSkimDataset,0)
};

#endif

