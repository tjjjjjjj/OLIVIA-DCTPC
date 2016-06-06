{
  gROOT->Reset();
  gROOT->LoadMacro("DCTPCTree_misc.C");
  DCTPCTree t;
  t->Loop();
}
