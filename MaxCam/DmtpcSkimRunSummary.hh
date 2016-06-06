#ifndef DMTPCSKIMRUNSUMMARY
#define DMTPCSKIMRUNSUMMARY

#include "TROOT.h"
#include "TSystem.h"
#include "TObject.h"
#include "TVectorD.h"
#include "TGraph.h"
#include "TH2.h"
#include "TH1.h"
#include "TString.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TAxis.h"
#include <vector>
using namespace std;

class DmtpcSkimRunSummary : public TObject{

public:

   DmtpcSkimRunSummary();
   DmtpcSkimRunSummary(int cam, int nev);

//   DmtpcSkimRunSummary(const DmtpcSkimRunSummary &other);

   virtual ~DmtpcSkimRunSummary();
  
   virtual const char* GetName() const {return "DmtpcSkimRunSummary";}

   void setNEvents(int ncamera, int nev);

   void setOutDir(TString dir){outdir = dir;}
   TString getOutDir(){return outdir;}

   void setOutName(TString file){outname = file;}
   TString getOutName(){return outname;}

   void setImageKey(TString k){key = k;}
 
   int openFile();
   int openFile(TString dir, TString file);
   void closeFile();
   void setStyle(); 

   int setImageRMS(int cam, int n, double rms);
   int setImageMean(int cam, int n, double mean);
   int setNTracks(int cam, int n, int nt);
   int setNburnin(int cam, int n, int nb);
   int setMaxPixel(int cam, int n, double mp);


   void fillSpark(int cam){totalSparks[cam]++;}
   void fillBurnin(int cam, int rbi){totalBurnin[cam]+= rbi;}

   void createImages(int cam);
   void outputCamera(int cam);
   void beginFile();
   void endFile();

   void outputAll();
   int makeTeXFile();
   int makeTeXFile(TString dir, TString file);
   int pdfLatex();
   int tar();
   void cleanup(bool deletePics = false, bool deleteTeX = false, bool deletePDF = false);
   int runAll(bool deletePics = false, bool deleteTeX = false, bool deletePDF = false);

private:
  
  FILE * outfile;	//!
  int nevents;
  int ncam;
  TString key;
  TString outdir;
  TString outname;
  vector<double> imageRMS;
  vector<double> imageMean;
  vector<double> maxPixel;
  vector<int> totalSparks;
  vector<int> totalBurnin;
  vector<int> totalTracks;
  vector<int> ntracks;
  vector<int> nburnin;
  vector<int> nkilled;
  vector<double> aveMean;
  vector<double> aveRMS;
  vector<double> aveMaxPixel;
  

  ClassDef(DmtpcSkimRunSummary,1);
};



#endif
