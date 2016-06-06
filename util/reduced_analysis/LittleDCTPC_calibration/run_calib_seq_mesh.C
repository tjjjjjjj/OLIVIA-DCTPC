{
  gROOT->Reset();
  gROOT->LoadMacro("DCTPCTree_calib_seq_mesh.C");
  DCTPCTree t;
  t->Loop();
}
