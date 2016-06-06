{
  gROOT->Reset();
  gROOT->LoadMacro("DCTPCTree_sequence_energy_calib.C");
  DCTPCTree t;
  t->Loop();
}
