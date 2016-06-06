{
  gROOT->Reset();
  gROOT->LoadMacro("DCTPCTree_rate.C");
  DCTPCTree t;
  t->Loop();
}
