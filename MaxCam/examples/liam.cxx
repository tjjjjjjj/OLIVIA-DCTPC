//
// This file is based on Denis' Mesh_NIM_plots.cxx code
//

// data is stored in /net/zwicky/dmtpc/data
TString dataDir("/net/zwicky/dmtpc/data/");

doliam() {
  
  //////////////////////////////////////////
  /////  Don't delete this block  //////////
  //////////////////////////////////////////
  // useful data files
  // ccdrun_00383.root 
  //   events:  418, 948, 762, 1248
  // ccdrun_00384.root 
  //   events:  238, 1202
  // ccdrun_00388.root 
  //   events:  1462
  //////////////////////////////////////////
  //////////////////////////////////////////

  MaxCamRead* ana= new MaxCamRead(dataDir+"ccdrun_00383.root");
  // stores bias frame to be subtracted from raw image on getEvent() call
  ana->readBiasFrame();
  // reads a list of hot pixels and sets them to zero 
  // (suppresses them from images)
  ana->findHotPixels("hotpixels.dat");

  // total number of images
  Int_t nentries = ana->tree()->GetEntries();
  cout << "n events = " << nentries << endl;

  // get and plot a particular event
  ana->getEvent(1248);
  TH2F *image = ana->ccdImage();
  image->Draw("colz");

}
