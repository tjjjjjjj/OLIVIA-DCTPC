{
  gROOT->Reset();
  gROOT->LoadMacro("DCTPCTree_rate_weighted.C");
  DCTPCTree t;
  t->Loop();
}
