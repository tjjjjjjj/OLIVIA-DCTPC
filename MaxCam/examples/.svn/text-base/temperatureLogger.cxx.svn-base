#include <iostream>
#include <fstream>
#include "TStopwatch.h"
#include "TTimeStamp.h"

#include "../MaxCamChannel.hh"

using namespace std;

void test() {

  MaxCamChannel* mcc0 = new MaxCamChannel("ni6229chan0", "ni6229chan0", 16, 0);
  MaxCamChannel* mcc1 = new MaxCamChannel("ni6229chan0", "ni6229chan0", 17, 0);
  MaxCamChannel* mcc2 = new MaxCamChannel("ni6229chan0", "ni6229chan0", 18, 0);
  MaxCamChannel* mcc3 = new MaxCamChannel("ni6229chan0", "ni6229chan0", 19, 0);

  
  ofstream myfile;
  myfile.open("temperatures.dat");
  myfile<<"# time(yyyymmdd-hhmmss) T0 T1 T2 T3" << endl;

  float dt=0.0;
  TStopwatch sw;
  TTimeStamp time;
  for (int ii=0; ii<3600; ii++) {
    mcc0->ni6229ReadADC();
    mcc1->ni6229ReadADC();
    mcc2->ni6229ReadADC();
    mcc3->ni6229ReadADC();
    
    myfile << time.GetDate() << "-" << time.GetTime() << " " 
	 << mcc0->val() << " " << mcc1->val() << " " 
	 << mcc2->val() << " " <<mcc3->val() << endl;

    // "sleep" for a second
    sw.Start();
    dt = 0.0;
    while (dt<1.0) {
      dt = sw.RealTime();
      sw.Start(kFALSE);
    }
    //printf("dt = %f\n", dt);

  }

}
