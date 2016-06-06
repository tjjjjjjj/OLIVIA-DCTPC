// -*- mode: C++ -*-

void html() {

  //gSystem->Load("fli-sdk-1.50/libfli/libfli.so");
  gSystem->Load("MaxCam.so");

   THtml html;
   html.MakeAll(1, "MaxCam*");
}
