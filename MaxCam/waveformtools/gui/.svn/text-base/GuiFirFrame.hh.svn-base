#ifndef __GUI_FIR_HH__
#define __GUI_FIR_HH__

#include "TGFrame.h"
/*
class TGCompositeFrame;
class TRootEmbeddedCanvas;
class TGTab;
class TGLabel;
class TGNumberEntry;
class TGTextButton;
class TGCheckButton;
class TGComboBox;
class TGWindow;
class TGMenuBar;
class TGPopupMenu;

class TH1;
class TF1;
*/

#include "TH1.h"
#include "TF1.h"
#include "TGMenu.h"
#include "TGWindow.h"
#include "TGComboBox.h"
#include "TGButton.h"
#include "TGNumberEntry.h"
#include "TGLabel.h"
#include "TGTab.h"
#include "TRootEmbeddedCanvas.h"


#include <vector>
class GuiFirFrame : public TGMainFrame
{
  public:
    GuiFirFrame(const TGWindow* win,unsigned int w=547,unsigned int h=437);
    virtual ~GuiFirFrame();

    void SetLogPlot();
    void SetCoefficient();
    void SetParameter();
    void SetFilter(bool update);
    void SetFilter(){SetFilter(false);}
    void UpdateFilter();
    void SetAutoUpdate();
    void SetAutoFreq();
    void UpdateCoefficients(bool updateFil);
    void UpdateCoefficients();
    void UpdateFrequency();
    void HandleFileMenu(int x);
    void HandleViewMenu(int x);
    void HandleOptionMenu(int x);
    void Terminate();
 private:
 
    void SetGUI(unsigned int w, unsigned int h);

    TGVerticalFrame* fLeftFrame;
    TGVerticalFrame* fRightFrame;
    TGHorizontalFrame* fTabFrame;
    TGHorizontalFrame* fBottomButtons;
    TGHorizontalFrame* fCoeffsFrame;
    TGHorizontalFrame* fCoeffsFrame2;
    TGHorizontalFrame* fParFrame;
    TGCompositeFrame* fGuiCompFrame;
    TGTab* fHistoTab;
    //
    TGCompositeFrame* fCoeffFrame;
    TRootEmbeddedCanvas* fCoeffECan;
    //
    TGCompositeFrame* fMagFrame;
    TRootEmbeddedCanvas* fMagECan;
    //
    TGCompositeFrame* fPhaseFrame;
    TRootEmbeddedCanvas* fPhaseECan;
    //
    TGCompositeFrame* fRealFrame;
    TRootEmbeddedCanvas* fRealECan;
    //
    TGCompositeFrame* fImagFrame;
    TRootEmbeddedCanvas* fImagECan;
    //
    TGLabel* fHistoBinLabel;
    TGNumberEntry* fHistoBinEntry;
    //
    TGTextButton* fParSetButton;
    TGNumberEntry* fParNumberEntry;
    TGNumberEntry* fParEntry;
    //
    TGTextButton* fCoeffSetButton;
    TGNumberEntry* fCoeffNumberEntry;
    TGNumberEntry* fCoeffEntry;
    //
    TGCheckButton* fAutoUpdateButton;
    TGTextButton* fCalcFreqButton;
    TGTextButton* fCoeffUpdateButton;
    TGComboBox* fFilterComboBox;
    TGLabel* fFilterLabel;
    TGNumberEntry* fNumCoeffEntry;
    TGLabel* fCoeffsLabel;
    TGLabel* fFilterProps;
    //
    TGCheckButton* fLogPlotCheck;

    TGMenuBar* fMenuBar;
    TGPopupMenu* fFileMenu;
    TGPopupMenu* fViewMenu;
    TGPopupMenu* fOptionMenu;


    //Histograms and coefficient properties
    TH1* coeffH; 
    TH1* magH;
    TH1* phH;
    TH1* reH;
    TH1* imH;
    bool fAutoUpdate;
    bool fAutoFreq;
    int fNCoeff;
    int fNPar;
    std::vector<double> fCoeffs; 
    std::vector<double> fPars;
    TF1* fFilterFun;
    int fNbin;

    ClassDef(GuiFirFrame,1)

};

#endif
