{
  gROOT->Reset();
  gROOT->LoadMacro("analysis_125CF4_dctpc_alpha_neon.C");
  analysis_125CF4_dctpc_alpha(229,293);//Am241 run at DC far
  //analysis_125CF4_dctpc_alpha(1041,1158);//1158//Cf252 runs at MIT
  gApplication->Terminate();
}
