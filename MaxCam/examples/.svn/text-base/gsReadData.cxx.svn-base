// 
//   Access data event-by-event
// 

void loadLibs() { 

   gSystem->Load("/usr/local/lib/libcfitsio.so");
   gSystem->Load("../MaxCam/MaxCam_linux.so");

}

void readEvent(int eventNo) {

  loadLibs() ; 
  // Access a run file on xroot server:
  // MaxCamRead *ana=new MaxCamRead("root://mitbbr11.lns.mit.edu//data/ccdrun_00039.root");
  // or as a local file
  MaxCamRead *ana=new MaxCamRead("../Runs/ccdrun_00039.root");

  // Read bias image if saved into the file
  ana->readBiasFrame();

  // Read specific event
  ana->getEvent(eventNo); 

  // Display image
  gStyle->SetPalette(1);
  ana->ccdImage()->Draw("colz");

  // fit and display tracks
  MaxCamTrack trfit(ana->ccdImage());
  trfit.makeTracks();
  cout << "Tracks ..... " << trfit.nTracks() << endl;
  for (int i=0; i<trfit.nTracks(); i++) trfit.getTrack(i)->Draw("same");

  // Find HV, pressure and CCD info
  cout << "Wire HV  .... " << ana->wire()->currentValue << endl;
  cout << "Mesh HV  .... " << ana->mesh()->currentValue << endl;
  cout << "Pressure  ... " << ana->pressure()->currentValue << endl;
  cout << "CCD expo .... " << ana->ccdConfig()->exposureTime << " ms" << endl; 
  cout << "CCD temp  ... " << ana->ccdConfig()->CCDTemp << endl;
  cout << "Time Stamp .. " << ana->timeStamp()->Print();
}




void gsReadData(int numberOfEventsToRead, int firstEvent) {

  loadLibs() ; 

  // Access a run file on xroot server:
  // MaxCamRead *ana=new MaxCamRead("root://mitbbr11.lns.mit.edu//data/ccdrun_00039.root");
  // or as a local file
  MaxCamRead *ana=new MaxCamRead("../Runs/ccdrun_00039.root");

  MaxCamTrack *trfit = new MaxCamTrack(); 

  // Read bias image if saved into the file
  ana->readBiasFrame();

  for (int i=firstEvent; i<firstEvent+numberOfEventsToRead; i++) { 
  
  // Read specific event
    cout << "Reading event number " << i << endl; 
    ana->getEvent(i); 

    // ----- 
    // Display image
    gStyle->SetPalette(1);  ana->ccdImage()->Draw("colz");

    // fit and display tracks
    trfit->setImage ( ana->ccdImage() );  
    trfit->makeTracks();
    cout << "Tracks ..... " << trfit->nTracks() << endl;
    for (int ii=0; ii<trfit->nTracks(); ii++) trfit->getTrack(ii)->Draw("same");

    // -----
    cout << "Wire HV  .... " << ana->wire()->currentValue << endl;
    cout << "Mesh HV  .... " << ana->mesh()->currentValue << endl;
    cout << "Pressure  ... " << ana->pressure()->currentValue << endl;
    cout << "CCD expo .... " << ana->ccdConfig()->exposureTime << " ms" << endl; 
    cout << "CCD temp  ... " << ana->ccdConfig()->CCDTemp << endl;
    cout << "Time Stamp .. " << ana->timeStamp()->Print();
    
   } 

  delete ana; 

}


// 
//   Access all data 
// 
void displayTree() {

  //  gSystem->Load("/usr/local/lib/MaxCam.so");
  loadLibs() ; 

  // Access a run file on xroot server:
  //MaxCamRead *ana=new MaxCamRead("root://mitbbr11.lns.mit.edu//data/ccdrun_00039.root");
  // or as a local file
  MaxCamRead *ana=new MaxCamRead("../Runs/ccdrun_00039.root");

  // Display pixel integral vs. time stamp
  ana->tree()->Draw("ccdImage.Integral():timeStamp.Get()");

}
