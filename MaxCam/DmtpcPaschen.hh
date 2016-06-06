#ifndef DMTPC_PASCHEN_HH
#define DMTPC_PASCHEN_HH

#include "TString.h"
#include "TF1.h"
#include <iostream>
using std::cout;
using std::endl;

class TGraph;

class DmtpcPaschen  {

public:

    // Ctors
  DmtpcPaschen(const char* fileName="Paschen_CF4.dat", TString delim="", TString opt="" );
  DmtpcPaschen(const DmtpcPaschen &other);
  virtual ~DmtpcPaschen();
  
  virtual TString GetName() { return TString("DmtpcPaschen"); }
    
  TGraph* getPaschenCurve() { return _pas; }
  void readPaschenData();

  void makeFitFunction();
  TF1* getFitFunction() { return _fit; }
  void testme() { cout << "_fit= " << _fit << endl; }
  //void testme() { cout << "_fit->GetParameter(0) = " << _fit->GetParameter(0) << endl; }
private:
  TGraph    *_pas;
  TString* _fname;
  TString* _delim;
  TF1* _fit;

  ClassDef(DmtpcPaschen,0)
};

#endif

