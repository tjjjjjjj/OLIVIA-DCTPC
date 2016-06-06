{
  gROOT->Reset();
  gROOT->LoadMacro("DCTPCTree_energy_comparisons.C");
  DCTPCTree t;
  t->Loop();
}
