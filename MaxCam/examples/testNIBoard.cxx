#include <iostream>
#include "../MaxCamChannel.hh"

using namespace std;

void test() {

  MaxCamChannel* mcc = new MaxCamChannel("ni6229chan0", "ni6229chan0", 0, 0);
  
  double volts = mcc->ni6229ReadADC();
  cout << mcc->val() << "  " << mcc->rms() << endl;
}
