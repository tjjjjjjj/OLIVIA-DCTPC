{
  gROOT->Reset();
  gROOT->LoadMacro("analysis_125CF4_bigdctpc_neutron.C");
  analysis(980,985);//37-98 is 700 V, 114-148 is 715 V
  //151- is alpha run
  //  gApplication->Terminate();
}
