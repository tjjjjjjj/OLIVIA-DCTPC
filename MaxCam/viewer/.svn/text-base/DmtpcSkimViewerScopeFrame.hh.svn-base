#ifndef DMTPC_SKIM_VIEWER_SCOPE_FRAME
#define DMTPC_SKIM_VIEWER_SCOPE_FRAME

#include "TGFrame.h"
#include "../DmtpcSkimDataset.hh"
#include "../DmtpcSkimEvent.hh"
#include "../DmtpcEvent.hh"
#include "../ScopeDataInfo.hh"

// waveformtools includes
//  for analysis of generic waveforms
#include "../waveformtools/include/WaveformVector.hh"
#include "../waveformtools/include/SkimWaveform.hh"
#include "../waveformtools/include/DmtpcPulse.hh"
//  for analysis of voltage sensitive amplitifier waveforms
#include "../waveformtools/include/FastWfVector.hh"
#include "../waveformtools/include/FastWaveform.hh"
#include "../waveformtools/include/FastPulse.hh"
//  for analysis of charge sensitive pre-amp waveforms
#include "../waveformtools/include/CspWfVector.hh"
#include "../waveformtools/include/CspWaveform.hh"
#include "../waveformtools/include/CspPulse.hh"
//  for analysis of pmt waveforms
#include "../waveformtools/include/PMTWfVector.hh"
#include "../waveformtools/include/PMTWaveform.hh"
#include "../waveformtools/include/PMTPulse.hh"
//  general functions for waveform analysis
#include "../waveformtools/include/WaveformAnalysis.hh"
#include "../waveformtools/include/WaveformTools.hh"

#include "TArrow.h"
#include "TLine.h"
#include "TRootEmbeddedCanvas.h"
#include "TGComboBox.h"
#include "TGNumberEntry.h"
#include "TGLabel.h"
#include "TObjArray.h"
#include "TCanvas.h"

class DmtpcSkimViewerScopeFrame : public TGMainFrame
{

  public:
    DmtpcSkimViewerScopeFrame(const TGWindow *p, UInt_t w, UInt_t h, bool * show); 
    virtual ~DmtpcSkimViewerScopeFrame(){canvas->Delete(); *show_traces = false; }
    void Display(DmtpcEvent * evt, TObjArray * wfvecs); 
    void GroupChanged(Long_t id); 
    void ChannelChanged(Int_t id); 
    void DrawTrace(); 
    void DrawAll(); 
    void UpdateAll(); 
    void CloseAll(){show_all=false;} 

    void ShowAlone(); 
    void UpdateAlone(); 
    void CloseAlone(){show_alone=false;} 
  
    TString GenerateWfVectorText(TObject * wfvec);
    void DrawWfAnalysis(TObject * wfvec);

  private: 
    TRootEmbeddedCanvas * canvas; 
    bool * show_traces; //!
    TGNumberEntry * group_select; 
    TGComboBox * channel_select; 
    TGLabel * trace_text; 
    TObjArray * waveform_vectors; 
    DmtpcEvent * original_event;
    int current_group;
    int current_chan;
    TCanvas * all; 
    TCanvas * alone_c; 
    bool show_all;
    bool show_alone; 

    TObjArray * algorithm_draw_objects;

  ClassDef(DmtpcSkimViewerScopeFrame,0);
};

#endif

