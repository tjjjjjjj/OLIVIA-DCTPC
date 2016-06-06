#ifndef ANALYSIS_RESULT_HH
#define ANALYSIS_RESULT_HH

#include "analysis_constants.hh"
#include "TNamed.h"

class AnalysisResult : public TObject
{

  public:
    AnalysisResult(); 
    AnalysisResult(int ncam, int nev, UInt_t * data = NULL); 
    AnalysisResult(const AnalysisResult & other); 

    void set_result(int i, int c, int t, UInt_t result) {_data[_ncam*15*i + 15*c + t] = result;} 
    UInt_t get_result(int i, int c, int t){ return _data[_ncam*15*i+15*c+t]; }
    Int_t get_nev() { return _nev;}
    Int_t get_ncam() { return _ncam;}
    UInt_t * get_data() { return _data;}
    Int_t get_arr_size() { return _arr_size;}

    bool pass_rbi(int i, int c, int t) {return (bool)(get_result(i,c,t) & PASS_RBI);}
    bool pass_alpha(int i, int c, int t) {return (bool)(get_result(i,c,t)& PASS_ALPHA);}
    bool pass_worm(int i, int c, int t) {return (bool)(get_result(i,c,t) & PASS_WORM);}
    bool pass_considered(int i, int c, int t) {return (bool)(get_result(i,c,t) & CONSIDERED);}

    ~AnalysisResult(); 

  private:
    Int_t _ncam;
    Int_t _nev; 
    Int_t _arr_size; 
    UInt_t * _data; //[_arr_size]
    
  ClassDef(AnalysisResult,1);

};
#endif
