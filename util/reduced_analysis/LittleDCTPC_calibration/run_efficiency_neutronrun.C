{
  gROOT->Reset();
  gROOT->LoadMacro("DCTPCTree_efficiency_neutronrun.C");
  DCTPCTree t;
  t->Loop();
}
