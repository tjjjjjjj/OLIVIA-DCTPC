//=================================================
// class Varia
//=================================================
#ifndef Tools_hh
#define Tools_hh

#include <math.h>
#include "TH1F.h"
#include "TString.h"
#include "TGraph.h"
#include "TROOT.h"

class Tools {

public:
  Tools();
  virtual ~Tools() {}

  TGraph* RejVsEffUnit(const TH1F* h1, const TH1F* h2); 
  void RejVsEff(int nsets, const TH1F* h1, const TH1F* h2,const TH1F* h3,const TH1F* h4, float lumi1,float lumi2, float lumi3,float lumi4,TString myCan, bool draw ); 
  void RejVsEffComp(const TH1F* h1, const TH1F* h2, const TH1F* h3, const TH1F* h4,TString myCan, bool draw ); 
  void OptimizeCut2(const TH1F* h1, const TH1F* h2, float lumi1,float lumi2,TString myCan ); 
  void OptimizeCut3(const TH1F* h1, const TH1F* h2, const TH1F* h3, float lumi1 ,float lumi2,float lumi3,
                   TString myCan, bool draw ); 

  void Draw4Hist(TH1F* H_Signal,TH1F* H_Backgd,TH1F* H_Signal_old,TH1F* H_Backgd_old,TString myCanName); 
  void DrawTwoHistosNorm(const TH1F* h1, const TH1F* h2, TString canvasName);  

  void DrawNHist(Int_t count, ...);

  void OptimizeCutN(bool draw, Int_t numHists, TString myCan, Int_t sigColor, Int_t bkgColor, ... );
  void CompareData_MC(bool draw, Int_t numHists, TString myCan, Int_t sigColor, Int_t bkgColor, ... );
  
  void OptimizeCutN2(TH1F* histoArray[], double lums[], Int_t sigColor, Int_t bkgColor, bool draw, Int_t numHists, TString myCan); 
  void OptimizeAndCompareWData(const TH1F* h1, const TH1F* h2, const TH1F* h3, const TH1F* h4,TString myCan="testCan");
  
  void DrawNHist2(TH1F* histArray[], Int_t colorArray[], Int_t arraySize);

  void PrintSignificance(TH1F* hsig1 ,TH1F* hsig2 ); 

  //ClassDef(Tools,0)

};

#endif
