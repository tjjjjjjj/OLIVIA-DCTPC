{
  gROOT->Reset();
  gROOT->LoadMacro("analysis_125CF4_bigdctpc_alpha_sim.C");
  analysis_125CF4_bigdctpc_alpha_sim(10000,10050);
  gApplication->Terminate();
}
