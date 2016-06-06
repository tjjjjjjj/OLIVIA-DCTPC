{
  gROOT->Reset();
  gROOT->LoadMacro("analysis_125CF4_dctpc_alpha.C");
  analysis_125CF4_dctpc_alpha(9532,9532);//Am241 run at DC far
  //analysis_125CF4_dctpc_alpha(1041,1158);//1158//Cf252 runs at MIT
  gApplication->Terminate();
}
