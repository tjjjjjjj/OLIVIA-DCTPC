{
  gROOT->Reset();
  gROOT->LoadMacro("DCTPCTree_stability.C");
  DCTPC_runtree t;
  t->Loop();
}
