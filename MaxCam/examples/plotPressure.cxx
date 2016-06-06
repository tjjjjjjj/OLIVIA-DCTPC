// does not work....
// jbattat August 28, 2009

#include "TFile.h"
//#include "TString.h"
#include <iostream>
#include "Riostream.h"
#include <string>

using namespace std;

void plot() {

  ifstream in;
  in.open("pressurefile.dat");
  Float_t inficon, convectron;
  Int_t nlines = 0;
  string datetime;

  //TFile *f = new TFile("pressure.root", "RECREATE");
  while (1) {
    cout << "nlines = " << nlines << endl;
    if (nlines > 0) {
      in >> datetime >> inficon >> convectron;
      if (nlines < 5) 
	printf("inficon, convectron = %f, %f\n", 
	       inficon, convectron);
      if (!in.good()) break;
    }
    nlines++;
  }


}
