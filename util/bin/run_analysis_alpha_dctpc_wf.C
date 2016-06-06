{
  gROOT->Reset();
  gROOT->LoadMacro("analysis_125CF4_dctpc_alpha_wf.C");
  analysis_125CF4_dctpc_alpha_wf(10672,10680);
  gApplication->Terminate();
}
