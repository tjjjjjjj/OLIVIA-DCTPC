#include "DmtpcAmptekMCA8000.hh"
#include "TH1F.h"

#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <iostream>
using std::string;
using std::ifstream;
using std::cout;
using std::endl;

ClassImp(DmtpcAmptekMCA8000)

DmtpcAmptekMCA8000::DmtpcAmptekMCA8000() {
}

DmtpcAmptekMCA8000::DmtpcAmptekMCA8000(const char* filename) {
  string line;
  ifstream specfile(filename);
  bool parseData = false;

  int nbins = 2048; // need to find this on the fly...
  _spec = new TH1F("_spec","_spec",nbins,0,nbins-1);

  int ndata = 0;
  while (!specfile.eof() ) {
    getline(specfile, line);
    // trim leading and trailing whitespace

    //cout << "line = " << line << endl;
    //cout << "line.find('<<DATA>>') = " << line.find("<<DATA>>") << endl;
    // ignore empty lines
    if (line.size() == 0) continue;
    // ignore comment lines (lines whose first non-whitespace char is #)
    //if (line.find("#") == 0) continue;
    if (line.find("<<DATA>>")==0) {
      cout << "found <<DATA>>" << endl;
      parseData = true;
    }
    if (line.find("<<END>>")==0) {
      cout << "found <<END>>" << endl;
      parseData  = false;
    }
    if (parseData) {
      _spec->SetBinContent(ndata, atoi(line.c_str()));
    }
    ndata++;
  }
  _spec->Draw();
}
