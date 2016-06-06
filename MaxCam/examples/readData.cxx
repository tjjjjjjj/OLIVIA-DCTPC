// 
//   Access data event-by-event
//
//   readEvent(12);
//   readEvent(5);
// 
//
//
MaxCamRead *ana=0;
void init() {
  cout << "Initialization..."<<endl;

   gSystem->Load("/usr/local/lib/MaxCam.so");

  // Access a run file on xroot server:
  //ana=new MaxCamRead("root://mitbbr11.lns.mit.edu//data/ccdrun_00039.root");
  // or as a local file
  ana=new MaxCamRead("data/ccdrun_00139.root");

  // Read bias image if saved into the file
  ana->readBiasFrame(); 

  // Read hot pixels from file
  ana->findHotPixels("hotpixels.dat");
  // or make a search for channels that are 5 sigma 
  // away from mean 10% of time using 100 events of the
  // curent data sample
  //ana->findHotPixels(100, 5, 0.1);
  // hot pixel can be added with MaxCamRead::addHotPixel(binx,biny)

  cout << "Total number of events="<< ana->tree()->GetEntries()<<endl;
}


void readEvent(int eventNo) {

  if (!ana) init();

  // Read specific event
  ana->getEvent(eventNo); 

  // Display image
  ana->ccdImage()->DrawCopy("colz");

  // Find HV, pressure and CCD info
  cout << "Wire HV  .... " << ana->wire()->currentValue << endl;
  cout << "Mesh HV  .... " << ana->mesh()->currentValue << endl;
  cout << "Pressure  ... " << ana->pressure()->currentValue << endl;
  cout << "CCD expo .... " << ana->ccdConfig()->exposureTime << " ms" << endl; 
  cout << "CCD temp  ... " << ana->ccdConfig()->CCDTemp << endl;
  cout << "Time Stamp .. " << ana->timeStamp()->Print();
}


// 
//   Access all data 
// 
void displayTree() {

  if (!ana) init();

  // Display pixel integral vs. time stamp
  ana->tree()->Draw("ccdConfig.CCDTemp:timeStamp.Get()");

}
