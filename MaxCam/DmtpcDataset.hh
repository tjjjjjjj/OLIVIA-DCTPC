
#ifndef DMTPC_DATASET_HH
#define DMTPC_DATASET_HH

#include <assert.h>

#include "TROOT.h"
#include "TFile.h"
#include "TObject.h"
#include "TTree.h"
#include "TObjString.h"
#include "TH2F.h"
#include "TTimeStamp.h"
class TChain;
class TTree;
class TH2F;
class TObjString;

#include "DmtpcEvent.hh"
#include "DmtpcDetectorPart.hh"

/**
A class for reading and writing raw data files.
*/
class DmtpcDataset : public TObject {

public:
    //
    // Ctors

/**
Constructor
*/
    DmtpcDataset();

/**
Copy constructor
\param other DmtpcDataset to be copied
*/    
    DmtpcDataset(const DmtpcDataset &other);
    
/**
Destructor
*/
    virtual ~DmtpcDataset();
    
/**
Get the name of the dataset
\return Always returns "DmtpcDataset"
*/
    virtual const char* GetName() const { return "DmtpcDataset"; }

/**
Write the tree, keyword, location and comment to the current file.
*/    
    void write() { _file->Write(); _keyword->Write("keyword"); _location->Write("location"); _comment->Write("comment");}
    
/**
\return the tree of DmtpcEvents
*/
    TTree* tree() { return _imageTree; }
    
/**
\return the tree of DmtpcEvents as a TChain object
*/
    TChain* chain() { return (TChain*)_imageTree; }
    
/**
\return a pointer to the file
*/
    TFile* file() { return _file; }
    
/**
Creates and opens a root file for writing/reading.
\param fname the file name (ends in .root)
\param foption file opening option. "recreate" to create a new file, "" for reading only
*/
    void createRootFile( const char *fname, TString foption);
    
/**
Opens a preexisting root file.  Calls createRootFile(fname,"").
\sa createRootFile()
\param fname the file name
*/
    void openRootFile( const char *fname); 

/**
Deletes any open trees and closes the file, if open.
*/
    void closeRootFile();   
 
/**
Fill tree entry with current DmtpcEvent.
*/
    void fill(); 

/**
Sets tree to "SaveSelf" autosave mode.  See root TTree documentation.
*/
    void autoSave();

/**
\return the name of the current file
*/
    const char* getFileName(); 
    
/**
Clear the memory of the pointers in the current DmtpcEvent.  Necessary for instances where CINT won't properly call the DmtpcEvent destructor.
*/
    void clearEventMemory();
    
/**
\return pointer to the current DmtpcEvent.
*/
    DmtpcEvent* event() { return _event; }
    
/**
DAQ only. Saves run to a directory and updates a run database.
\param scpdest the file destination
\param fname a text file holding MySQL database info
\return 0 upon completion
*/
    int saveRun(char *scpdest="mitdm10.mit.edu:/data/", 
                char *fname="dbaccess.txt");
    
/**
Loads a tree entry (DmtpcEvent) into memory
\param i the tree entry to load
*/
    void getEvent(int i);
    
/**
Retrieves a bias frame for this run
\param iccd the array number of the ccd
\return bias frame histogram
*/
    TH2F *getBiasFrame(int iccd);

/**
Retrieves the overscan bias frame for this run
\param iccd the array number of the ccd
\return bias frame histogram for overscan
*/
  TH2F *getBiasFrameOverscan(int iccd);

/**
Returns the number of bias frames present in this run
\return bias frame histogram

Actually, wherever you see "bias", you should replace 
it with "dark" since they are finite-length exposures
*/
    Int_t getNbiasFrames();

/**
Returns the number of cameras present in the run without
requiring a call to getEvent()
\return number of cameras in the dataset
*/
  Int_t getNcameras();


/**
Subtract bias frame from image held in current DmtpcEvent
\param iccd the array number of the CCD we want to use
*/
    void darksub(int iccd) { event()->ccdData(iccd)->Add(getBiasFrame(iccd+1), -1); }

/**
Draw an image from the current DmtpcEvent
\param iccd the array number of the CCD we want to use
*/
    void showimg(int iccd) { event()->ccdData(iccd)->Draw("colz"); }
    
/**
DAQ only. Find an unused run number from a MySQL run database and use it to set the file name
\param dbAccessFile a text file holding information to read from a database
\return 0 upon completion
*/
    int setRunNumberFromDB(const char *dbAccessFile);
    
    void setComment(TString c);
    TObjString *comment() { return _comment; }

    void setKeyword(TString c);
    TObjString *keyword() { return _keyword; }
    
    void setLocation(TString c);
    TObjString *location() { return _location; }

    void setDetId(TString c);
    TObjString *detId() {return _detId;}
    
    TClonesArray* listOfDetectorParts() { return _listOfDetectorParts; }
/**
\param i the number of the part to retrieve
\return a part from the list of detector parts
*/
    DmtpcDetectorPart* getDetectorPart(int i) { return (DmtpcDetectorPart*)listOfDetectorParts()->At(i); }
    
  /**
     Only valid during data collection
     \returns where the file was saved to in saveRun.
  */
  TString getSavedLocation() {return _savedLocation;}

private:
    
    TTree* _imageTree; ///< tree holding many events
    
    TFile* _file; ///< file containing the tree
    
    DmtpcEvent *_event;///< Data from one event
    
    TObjString *_comment; ///< brief comment about run
    TObjString *_location; ///< run location
    TObjString *_keyword; ///< run keyword
    TObjString  * _detId; ///< detector ID tag. used for setting file names 
 
    TString _savedLocation;
    
    TClonesArray *_listOfDetectorParts;///< Array of detector parts

    TTimeStamp* _inittime;
    

  ClassDef(DmtpcDataset,2)
};

#endif

