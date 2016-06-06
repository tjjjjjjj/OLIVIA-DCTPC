// call python code

#include <fstream>

#include <iostream>
using namespace std;
#include "TString.h"

void getpower() {

  ifstream fin("power.tmp");
  TString portname;
  TString portstatus;
  while (!fin.eof()) {
    fin >> portname >> portstatus;
    cout << "portname, status = " << portname << ", " << portstatus << endl;
  }
  fin.close();
  
}
