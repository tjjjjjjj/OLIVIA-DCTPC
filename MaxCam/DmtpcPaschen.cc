#include "DmtpcPaschen.hh"
#include "TGraph.h"
#include "TF1.h"

#include <iostream>
using std::cout;
using std::endl;

ClassImp(DmtpcPaschen)

//____________________
//
DmtpcPaschen::DmtpcPaschen(const char *fileName, TString delim, TString opt)  {
  // Constructor uses name of file with Paschen curve data
  //
  _fname = new TString(fileName);
  _delim = new TString(delim);
  readPaschenData();
  // nothing done with options yet...
  //opt.ToLower(); 
  //if (opt.Contains("scat"))      fillScatteringTable( fileName );
  //else if (opt.Contains("fis"))  fillFissionTable( fileName );
  //else if (opt.Contains("cs"))  fillCSTable( fileName );
  //else assert(0);
}


DmtpcPaschen::DmtpcPaschen(const DmtpcPaschen &other) {
        _pas=other._pas;
        cout << GetName() << "Copy Constructor not done" << endl;
}


DmtpcPaschen::~DmtpcPaschen() {
  cout << "DmtpcPaschen() destructor not done..." << endl;
}

void DmtpcPaschen::readPaschenData() {
  _pas = new TGraph(_fname->Data());//"", _delim->Data());
}

void DmtpcPaschen::makeFitFunction() {
  cout << "=====" << endl;
  cout << "makeFitFunction(): at present, assume that the data is CF4 paschen curve data" << endl;
  cout << "=====" << endl;
  _fit = new TF1("_fit", "[0]*x/(log(x-[2])+[1])", 0.05, 2);
  _fit->SetParameter(0, 521.);
  _fit->SetParameter(1, 0.5);
  _fit->SetParameter(2, -0.6);
  _fit->SetLineColor(kRed);

  _pas->Fit(_fit);

}


  //// below is from TGraph() constructor
  //Double_t x, y;
  //
  ////gSystem->ExpandPathName(_fname);
  //ifstream infile(_fname.Data());
  //if (!infile.good()) {
  //  Error("DmtpcPaschen:", "Cannot open file: %s", _fname->Data());
  //  return;
  //}
  //
  ////if (!CtorAllocate()) return;
  //
  //std::string line;
  //Int_t np = 0;
  //const char* format = "%lg %lg";
  //
  //// No delimiters specified (standard constructor).
  //if (delim=="") {
  //  
  //  while (std::getline(infile, line, '\n')) {
  //    if (2 != sscanf(line.c_str(), format, &x, &y)) {
  //	continue; //skip empty and ill-formed lines
  //    }
  //    _pas.SetPoint(np, x, y);
  //    np++;
  //  }
  //  
  //  // A delimiter has been specified in "option"
  //} else {
  //  cout << "delimiter code not yet implemented" << endl;
  //  //// Checking format and creating its boolean counterpart
  //  //TString format_ = TString(format) ;
  //  //format_.ReplaceAll(" ", "") ;
  //  //format_.ReplaceAll("\t", "") ;
  //  //format_.ReplaceAll("lg", "") ;
  //  //format_.ReplaceAll("s", "") ;
  //  //format_.ReplaceAll("%*", "0") ;
  //  //format_.ReplaceAll("%", "1") ;
  //  //if (!format_.IsDigit()) {
  //  //  Error("TGraph", "Incorrect input format! Allowed formats are {\"%%lg\",\"%%*lg\" or \"%%*s\"}");
  //  //  return;
  //  //}
  //  //Int_t ntokens = format_.Length() ;
  //  //if (ntokens < 2) {
  //  //  Error("TGraph", "Incorrect input format! Only %d tag(s) in format whereas 2 \"%%lg\" tags are expected!", ntokens);
  //  //  return;
  //  //}
  //  //Int_t ntokensToBeSaved = 0 ;
  //  //Bool_t * isTokenToBeSaved = new Bool_t [ntokens] ;
  //  //for (Int_t idx = 0; idx < ntokens; idx++) {
  //  //  isTokenToBeSaved[idx] = TString::Format("%c", format_[idx]).Atoi() ; //atoi(&format_[idx]) does not work for some reason...
  //  //  if (isTokenToBeSaved[idx] == 1) {
  //  //	ntokensToBeSaved++ ;
  //  //  }
  //  //}
  //  //if (ntokens >= 2 && ntokensToBeSaved != 2) { //first condition not to repeat the previous error message
  //  //  Error("TGraph", "Incorrect input format! There are %d \"%%lg\" tag(s) in format whereas 2 and only 2 are expected!", ntokensToBeSaved);
  //  //  return;
  //  //}
  //  //
  //  //// Initializing loop variables
  //  //Bool_t isLineToBeSkipped = kFALSE ; //empty and ill-formed lines
  //  //char * token = NULL ;
  //  //TString token_str = "" ;
  //  //Int_t token_idx = 0 ;
  //  //Double_t * value = new Double_t [2] ; //x,y buffers
  //  //Int_t value_idx = 0 ;
  //  //
  //  //// Looping
  //  //while (std::getline(infile, line, '\n')) {
  //  //  if (line[line.size() - 1] == char(13)) {  // removing DOS CR character 
  //  //	line.erase(line.end() - 1, line.end()) ;
  //  //  }
  //  //  if (line != "") {
  //  //	token = strtok(const_cast<char*>(line.c_str()), option) ;
  //  //	while (token != NULL && value_idx < 2) {
  //  //	  if (isTokenToBeSaved[token_idx]) {
  //  //	    token_str = TString(token) ;
  //  //	    token_str.ReplaceAll("\t", "") ;
  //  //	    if (!token_str.IsFloat()) {
  //  //	      isLineToBeSkipped = kTRUE ;
  //  //	      break ;
  //  //	    } else {
  //  //	      value[value_idx] = token_str.Atof() ;
  //  //	      value_idx++ ;
  //  //	    }
  //  //	  }
  //  //	  token = strtok(NULL, option) ; //next token
  //  //	  token_idx++ ;
  //  //	}
  //  //	if (!isLineToBeSkipped && value_idx == 2) {
  //  //	  x = value[0] ;
  //  //	  y = value[1] ;
  //  //	  SetPoint(np, x, y) ;
  //  //	  np++ ;
  //  //	}
  //  //  }
  //  //  isLineToBeSkipped = kFALSE ;
  //  //  token = NULL ;
  //  //  token_idx = 0 ;
  //  //  value_idx = 0 ;
  //  //}
  //
  //// Cleaning
  ////delete [] isTokenToBeSaved ;
  ////delete [] value ;
  ////delete token ;
  //}
  //
  //infile.close();

