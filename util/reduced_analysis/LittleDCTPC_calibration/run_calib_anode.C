{
  gROOT->Reset();
  gROOT->LoadMacro("DCTPCTree_calib_anode.C");
  DCTPCTree t;
  t->Loop();
}
