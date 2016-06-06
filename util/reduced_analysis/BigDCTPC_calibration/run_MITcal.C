{
  gROOT->Reset();
  gROOT->LoadMacro("MITcal.C");
  DCTPCTree t;
  t->Loop();
}
