{
  gROOT->Reset();
  gROOT->LoadMacro("DCTPCTree_position_calib.C");
  DCTPCTree t;
  t->Loop();
}
