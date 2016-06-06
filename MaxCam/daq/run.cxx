void run(int nevt=1000, int exposureTime=1000, 
	 TString comment="default comment", 
	 TString keyword="default keyword", 
 TString location="default location",
 TString detid="")
{

  TStopwatch gt; 
  gt.Start(); 
  gSystem->Load("/usr/local/lib/dmtpcDAQ_cxx.so");
 
  ccdType="apogee";
  scopeType="dummy";
  globalExpoTime=exposureTime;
	
  init(0,"working.root",detid);
	

  // ------ EDIT IF NEEDED ------
  // Comment: entered into a database and root file:
  exper->data()->setComment(comment);

  // Keyword: database + root file
  exper->data()->setKeyword(keyword);

  // Location
  exper->data()->setLocation(location);

  // -------- END EDIT ----------


  // take data
  for (int iccd=0; iccd<exper->nCCD(); iccd++) { exper->ccd(iccd)->openShutter(); }
  event(nevt);
  for (int iccd=0; iccd<exper->nCCD(); iccd++) { exper->ccd(iccd)->closeShutter(); }

  // end run
  end();

  // save data
  saveRun("/ccddata/"); 

  cout << "Complete Runtime: " << gt.RealTime() << endl; 
  gSystem->Exit(0);

  return 0;

}
