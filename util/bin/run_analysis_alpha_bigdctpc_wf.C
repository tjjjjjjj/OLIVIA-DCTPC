{
  gROOT->Reset();
  gROOT->LoadMacro("analysis_125CF4_bigdctpc_alpha_wf.C");
  analysis_125CF4_bigdctpc_alpha_wf(889,954);
  gApplication->Terminate();
}
