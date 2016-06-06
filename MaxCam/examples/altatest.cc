#include "../MaxCamAlta.hh"
#include "ApnCamera.h"
#include <iostream>
using std::cout;
using std::endl;

int main() {
  MaxCamAlta ccd = MaxCamAlta(1);

  char sernumber[256];
  long buflength=256;
  ccd.getAlta()->GetCameraSerialNumber( sernumber, &buflength);
  cout << "sernumber = " << sernumber << endl;
}
