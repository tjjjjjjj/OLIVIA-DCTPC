{
  gSystem->Load("/usr/local/lib/dmtpcDAQ_cxx.so");
 
  ccdType="apogee";
  //scopeType="dummy";
  globalExpoTime=10;
 	
  init();

  int n=exper->nCCD();
  for (int i=0; i<n; i++) {
    while(1) {
      cout << "CCD " << i << " Temp=" << exper->ccd(i)->getTemperature()<<endl;
      gSystem->Sleep(5000);
      if  (exper->ccd(i)->getTemperature()<-19) break;
    }
  }

  cout << "CAMERA INITIALIZATION COMPLETE" << endl;
  end();

  gSystem->Exit(0);
}
