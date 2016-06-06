#ifndef MAXCAM_READ_HH
#define MAXCAM_READ_HH
#include "TROOT.h"
#include "TChain.h"
#include "TDatime.h"

#include "MaxCamDataset.hh"

class MaxCam;
class MaxCamConfig;
class MaxCamChannel;
class MaxCamMC;
class MaxCamSegment;
class TH2F;
class TH1F;
class TH1D;
class TDatime;
class TString;
class TFile;
class TF1;
class TCutG;
class TVector2;
class TLorentzVector;
class TH1D;


#include <vector>
using std::vector;

class MaxCamRead : public MaxCamDataset {

public:
    
    // Constructors
    MaxCamRead(TString filePath="test.root", TString outfile="");
    
    virtual ~MaxCamRead() {};
    
    virtual TString GetName() { return TString("MaxCamRead"); }
    
    void getEvent(int i);
    
    long NEvents() { return _imageTree->GetEntries(); }
    
    void setROI( int x0, int x1, int y0, int y1);
    
    TH1F* makeYieldHisto(TH2F* img=0, float min=-1, float max=-1);
    
    int convertIntoFits();
    
    void readBiasFrame();
    
    void deleteBiasFrame();
    
    void findHotPixels(TString hotFile);

    void findHotPixels(int nevents=100, float th=3, float rate=0.05);

    void deleteHotPixels();

    void addHotPixel(int binx, int biny);
    
    void addRun(TString fname);
    
    void addWire(double xwire);
    
    void removeWire(int iwire);
    
    vector<double> wireList; // list of wires (coordinates)
    
    vector<int> wireBinList; // list of wires (bins)
 
    vector<int> wireIdList; // list of wires to be merged; 

    void findWires(float sigma, float threshold, int nev=-1, TString opt="x", float maxAvgSum=-1);

    TH1D* wireYield(int iwire, int dwire, int type=0, TString opt="x"); 
    
    TH1F* findSegmentLength(int iwire, int dwire, float threshold, int pixelsBelowThreshold, TString opt="x");
    TH1F* mergeWireStrips(int dwire0, int dwire1, TString opt="x", TString hname="project" );
    TH2F* getImageTh() { return _imageth; }

    
    void fillWireBasic(int id, float y, float yerr, float ampl, float mean, float width, float chi2, float widthL) {
        Wire[id].id=id;
        Wire[id].y=y;  Wire[id].yerr=yerr;  Wire[id].ampl=ampl;
        Wire[id].mean=mean;  Wire[id].width=width;  Wire[id].chi2=chi2; Wire[id].widthL=widthL;
    }
        
        
    struct wireData {
        // wire
        float y, yerr, ampl, mean, width, chi2, widthL, widthR, skewness, ycount, shapeChi2, antiShapeChi2;
        int ndof;
        int id;
        float esrim; // energy from length based on SRIM 
    } *Wire; // extracted information for wire
    
    float getWidthL(int iw) { return Wire[iw].widthL; }
    float getWidthR(int iw) { return Wire[iw].widthR; }
    void  copyShapeToAntishape(int iw) { Wire[iw].antiShapeChi2=Wire[iw].shapeChi2; }
    
    struct wireChamberData {
        // experiment
        int time;
        float wirehv;
        float meshhv;
        float pressure, setpress;
        int   expotime;
        int event;
        int run;
        int fitnumber;
        float emc;
    } Wires; // extracted information for wire chamber
    
    TTree *nt; // ntuple with fit information
    
    virtual void makeNtuple(int minwire, int maxwire, int width, int type=0, TString opt="x");
    
    TH2F* accumulatePressure(float pressure, float avgmin=-1e10, float avgmax=1e10, int fst=-1, int nev=-1);
    
    void setRectangularCut(float x0, float y0, float x1, float y1, float d, TString opt="sg");
    void select(TH2F* image, TString opt="sg");
    TH2F* rotate(TH2F* image, int nrot=20);
    
    int  getRunNumber();
    void setEventInfo(int iev);
    void setFitNumber(int ifit) { Wires.fitnumber=ifit; }
    
    TF1* getProfileFunction() { return _fun; }
    TF1* getStoppingFunction() { return _fStopping; }
    
    int wireMinROI, wireMaxROI;

    static MaxCamMC *srim;
    
    void setVerbose(int v) { _verbose=v; }

    int  getVerbose() { return _verbose; }


    void resetDischargeCounter() {_ndischarge=0;}
    int  getDischargeCounter() { return _ndischarge; }
    bool isDischargeEvent(float threshold=20);

    bool hasSegments(int &imax, int &jmax, float &threshold);
    
    bool hasTracks();

    MaxCamSegment* segment() { return _segment; }


    void closeOutputFile() { _outputFile->Write();  _outputFile->Close(); }

    bool isMC() { return _recoil ? true : false; }
    
private:
    
    TCutG *_cutSG; // rectangular cut in image for sg
    TCutG *_cutBG; // rectangular cut for bg
    TVector2 *_xaxis, *_yaxis; // axes of cut
    
    TH2F  *_imageth; // analysis image.
    
    TH2F  *_biasFrame; // bias frame.
    
    int _x0,_x1,_y0,_y1; // ROI

    vector<int> _hotPixels;
     
    TF1 *_fun; // Profile perpendicular to wire: a gaussian + flat bg
    
    TF1 *_fStopping; // Profile along wire: dEdx from SRIM + flat background
    
    TDatime _T0;

    TLorentzVector *_recoil, *_projectile;

    void createWireBinList();

    TH1D * _xproj;

    int _verbose;

    int _ndischarge;

    MaxCamSegment *_segment;
    
    TFile *_outputFile;

    
    ClassDef(MaxCamRead,0)


};

#endif

