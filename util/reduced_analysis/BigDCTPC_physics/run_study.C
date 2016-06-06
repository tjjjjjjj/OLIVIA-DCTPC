{
  gROOT->Reset();
  gROOT->LoadMacro("DCTPCTree_study.C");
  DCTPCTree t;
  t->Loop();
}
