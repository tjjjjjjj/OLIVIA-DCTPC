{
  gROOT->Reset();
  gROOT->LoadMacro("DCTPCTree_efficiency_ccd.C");
  DCTPCTree t;
  t->Loop();
}
