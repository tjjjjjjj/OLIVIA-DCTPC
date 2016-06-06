//____________________________________________________________
//
//  An example for taking ccd images, saving data and reading
//  data files.
//  Script should be run on DAQ machine only.
// 
//____________________________________________________________

MaxCamTest *exper;
TH2F *hAcc=0;
TH1F *hxAcc=0, *hta=0, *htb=0;
int totalSaved=0;
TCanvas *cy=0, *cw0=0, *cw1=0;
TF1 *fun=0;

//________________________________________
//
//    TAKE CCD DATA
//________________________________________

void init(int debug=0, TString fname="test.root") {

  exper= new MaxCamTest(debug,fname);

  // Setup camera
  //exper->ccd()->setNumberOfFlushes(1); // cleaning of CCD before exposure
  //exper->ccd()->setVBin(8); exper->ccd()->setHBin(8); 

  exper->begin();

  //exper->makeBiasFrame(10);
  exper->makeBiasFrame("biasframe.root");
  //exper->findHotPixels(100, 4, 0.1);
  //exper->findHotPixels("hotpixels.dat");
  
  long expotime = 400;
  exper->ccd()->setNormalFrame();
  exper->ccd()->setExposureTime(expotime);       

  hta = new TH1F("hta","",1000,-0.5,0.5); 
  htb = new TH1F("htb","",1000,0,512);

  exper->createCanvas();
  cy= new TCanvas("cy", "", 850, 0, 400,250);
  //cw0= new TCanvas("cw0", "", 400,250);
  //cw1= new TCanvas("cw1", "", 400,250);
  exper->getImagePad()->cd();

  fun = new TF1("fun","[0]+[1]*exp(-0.5*(x-[2])**2/[3]**2)"); 

  setPressureAndSave(0);
}




void setMeshAndSave(double val) {
  exper->writeValueToDB("mesh_hv",val);
}

void setWireAndSave(double val) {
  exper->writeValueToDB("wire_hv",val);
}

void setPressureAndSave(double val) {
  exper->writeValueToDB("pressure",val);
}



int ntrigger=0;
int xmin0=12, xmax0=17;
int xmin1=50, xmax1=55; 
//
float yieldMin=-100, yieldMax=500;
// save=0 none
// save=1 triggered only
// save=-1 every 
int nrecovery=0;
int event(int n=1, int save=0) {    
  

  for (int i=0; i<n; i++) {

    bool isGood=false;

    exper->setTriggerTrials(ntrigger);
    exper->setSaveFlag(false);
    exper->event();

    cout << "EVENT HAS FINISHED"<<endl;

    // display ccd image
    exper->getImagePad()->cd(); TH2F* tmphist=exper->drawImage("colz",yieldMin,yieldMax);  


    // display projections
    //MaxCamImageTools::killLonePixels( tmphist, 5); tmphist->Draw("colz");
    cy->cd(); TH1F* hYY=exper->drawYields(tmphist, "gaus"); hYY->SetStats(11111); hYY->DrawCopy(); hYY->SetStats(); float mean=hYY->GetMean(); 
    if ( hYY->GetMean() < 100 ) { if (!hxAcc) hxAcc=(TH1F*)hYY->Clone("hxAcc"); else hxAcc->Add( hYY ); cout << "USE EVENT"<<endl; }
    delete hYY;
    cout << "Mean="<< mean<<endl;
    
    exper->getYProjPad()->cd(); TH1D* hY=tmphist->ProjectionY("_py"); hY->SetFillColor(1); hY->SetStats(); hY->DrawCopy("hbar"); 

    exper->getXProjPad()->cd(); TH1D* hX=tmphist->ProjectionX("_px");  hX->SetFillColor(1); hX->DrawCopy("bar"); 
    delete tmphist;   

    exper->getImageCanvas()->Update(); cy->Update(); //cw0->Update(); cw1->Update();



    // save event?
    if (!save) continue;

    else if (save>0) {
      
      if (mean>100&&0) { // discharge!!
	exper->recoverFromDischarge(9000);
	nrecovery++;
	if (nrecovery>5) return -1;
	continue;
      }
      
      //if (exper->trfit()->nTracks()!=1) continue; // only 1 track?
      if (exper->trfit()->nTracks()<1) continue; // at least 1 track?
      
      TF1* track=exper->trfit()->getTrack(0);
      if (!track) continue;
      float slope = track->GetParameter(1);
      //hta->Fill(slope);
      float intercept = track->GetParameter(0);
      //htb->Fill( intercept );
      delete track;
      cout << "Slope="<< slope << "   Intercept="<<intercept<<endl;
      //if ( fabs(slope-0.0)>0.15 ) continue;
      //if ( fabs(slope+0.7)>0.15 ) continue;
      //if ( fabs(intercept-340)>100) continue;// && fabs(intercept-420)>100 ) continue;
    }

    totalSaved++;
    exper->saveEvent();
    
    const char* e0="\033[44;37m";
    const char* en = "\033[0m";    
    cout <<e0 << "***** EVENT SAVED ****** (n=" << totalSaved << ")" << en<<endl;
    
  }

  return 0;

}
  

void end() {
  //hta->Write();
  //htb->Write();

  exper->end();
}

void saveRun(char *key, char *desc) {
  exper->saveRun(key, desc, "NW13-039");
}




void scanVoltage(float min=2.6, float max=2.86, float step=0.05, bool wire=true) {


    
  for (float v=min; v<max; v+=step) {
      cout << "Setting wire voltage to " << v << "kV" << endl;
      if (wire) setWireAndSave(v); 
      else setMeshAndSave(v);
      gSystem->Sleep(5000);
      event(200,-1);
  }


}


void timeRun(int n=24, long sleep=3600000) {

  for (int i=0; i<n; i++) {
    event(100,-1);
    gSystem->Sleep(sleep);
  }


}



void scanExposure() {
    
  for (int i=50; i<500; i+=25) {
    cout << "Setting exposure  to " << i << "msec" << endl;
    exper->ccd()->setExposureTime(i);
    event(50,1);
  }


}


void Fe55(int n=100) {
	if (!exper) init();
	yieldMin=10; yieldMax=110;
	event(n,0);
}
