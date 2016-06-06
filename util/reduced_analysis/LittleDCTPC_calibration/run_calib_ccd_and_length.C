{
  gROOT->Reset();
  gROOT->LoadMacro("DCTPCTree_calib_ccd_and_length.C");
  DCTPCTree t;
  t->Loop();
}
