#ifndef DMTPC_STITCH_GUI_FRAME_HH
#define DMTPC_STITCH_GUI_FRAME_HH

#define DEFAULT_DATA_DIR "/net/zwicky/esata01/dmtpc/data/4sh/:/data/4sh/"

#include "TGFrame.h"
#include "TApplication.h"
#include "../Dmtpc4ShooterStitcher.hh"
#include "TRootEmbeddedCanvas.h"
#include "TGNumberEntry.h"
#include "TGLabel.h"
#include "TGStatusBar.h"
#include <vector>


#define LENS_CORRECTION_ORDER 4


class DmtpcStitchGUIFrame : public TGMainFrame
{

  public: 

    DmtpcStitchGUIFrame(const TGWindow * p, UInt_t w, UInt_t h); 
    void CloseWindow() { gApplication->Terminate(0); } 
    void Load(); 
    void Train(); 
    void LoadStitch(); 
    void SaveStitch(); 
    void SaveStitched(); 
    void SaveOverlays(); 
    void DrawRaw(); 
    void DrawStitched(); 
    void DrawPolar(); 
    void DrawEdges(); 
    void DrawMedians(); 
    void SetPreTrained(); 
    void SetPreLoad(); 
    void Map(); 
    void SetTrained(); 
    void SetBusy(const char * text = "Busy!"); 
    void SetUnbusy(); 
    void Adjust(); 
    void HandleRadio(Int_t id); 

    virtual ~DmtpcStitchGUIFrame() {}; 


  private: 
    
    static void * doTraining(void * arg); 
    static void * doLoading(void * arg); 
    static void * doAdjusting(void * arg); 
    int ncam; 
    void drawSpacersMulti(); 
    void drawRingsMulti(); 
    TGNumberEntry * weights[4]; 
    TGNumberEntry * xcenter[4]; 
    TGNumberEntry * ycenter[4]; 
    TGNumberEntry * scale[4]; 
    TGNumberEntry * rot[4]; 
    TGLabel * ccdnames[4]; 


    TGNumberEntry * lensCorrection[LENS_CORRECTION_ORDER + 1];
    int radioChoice; 
    
    TGNumberEntry * blurLevel; 
    TGNumberEntry * edgeLowThreshold; 
    TGNumberEntry * edgeHighThreshold; 
    TGNumberEntry * minEdgeNeigbhors; 
    TGNumberEntry * linHoughRbins; 
    TGNumberEntry * linHoughThetabins; 
    TGNumberEntry * linHoughMinVotes; 
    TGNumberEntry * spacerJoinRthresh; 
    TGNumberEntry * spacerJoinThetathresh; 
    TGNumberEntry * circHoughFirstPassNbinsX; 
    TGNumberEntry * circHoughFirstPassNbinsY; 
    TGNumberEntry * circHoughFirstPassNbinsZ; 
    TGNumberEntry * circHoughSecondPassNbinsX; 
    TGNumberEntry * circHoughSecondPassNbinsY; 
    TGNumberEntry * circHoughSecondPassNbinsZ; 
    TGNumberEntry * circHoughMinX; 
    TGNumberEntry * circHoughMinY; 
    TGNumberEntry * circHoughMinZ; 
    TGNumberEntry * circHoughMaxX; 
    TGNumberEntry * circHoughMaxY; 
    TGNumberEntry * circHoughMaxZ; 
    TGNumberEntry * nwidthsSecondPass; 
    TGNumberEntry * nSpectrPeaksR; 
    TGNumberEntry * medianNBins; 
    TGNumberEntry * medianNiter; 
    TGNumberEntry * ccdWithLed; 
    TGNumberEntry * LEDBorderWidth; 
    TGNumberEntry * LEDThresh; 
//    TGNumberEntry * ImageHighThresh; 

    Dmtpc4ShooterStitcher * stitch; 
    TH2 * stitched; 
    TH2 * summed[4]; 
    TString cameraNames[4]; 
    int runNumber; 

    TRootEmbeddedCanvas * canvas; 
    TGLabel * fileLabel; 
    TGLabel * mapLabel; 
    TGStatusBar * busy; 
    TGGroupFrame * trainFrame; 
    TGGroupFrame * adjustframe; 
    char * filename; 
    char * mapname; 

    TGTextButton * trainButton; 
    TGTextButton * adjustButton; 
    TGTextButton * loadButton; 
    TGTextButton * mapButton; 
    TGTextButton * loadStitchButton; 
    TGTextButton * saveStitchButton; 
    TGTextButton * saveStitchedButton; 
    TGTextButton * saveOverlayButton; 
    TGCheckButton * drawSpacers; 
    TGCheckButton * drawRings; 

    void doGUI(); 


    ClassDef(DmtpcStitchGUIFrame,0); 
}; 

#endif
