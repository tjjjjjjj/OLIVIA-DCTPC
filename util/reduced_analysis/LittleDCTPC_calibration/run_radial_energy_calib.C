{
  gROOT->Reset();
  gROOT->LoadMacro("DCTPCTree_radial_energy_calib.C");
  DCTPCTree t;
  t->Loop();
}
