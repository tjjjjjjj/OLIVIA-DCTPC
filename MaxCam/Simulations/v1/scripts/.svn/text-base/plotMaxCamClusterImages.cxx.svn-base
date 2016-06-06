TString skimNameOfNumber(int runnum);
TString skimNameOfRun(TString runname);
TString runNameOfNumber(int runnum,TString detid="4sh");
TString compareNameOfRun(int runnum);

// function definitions
TString skimNameOfNumber(int runnum) {
  return skimNameOfRun(runNameOfNumber(runnum));
}

TString skimNameOfRun(TString runname) {
  return runname.ReplaceAll(".root", "skim.root");
}
TString runNameOfNumber(int runnum, TString detid) {
  return TString::Format("dmtpc_%s_mc_%05d.root", detid.Data(), runnum);
}

TString compareNameOfRun(int runnum) {
  return runNameOfNumber(runnum).ReplaceAll(".root", "_compare.root");
}

void plotMaxCamClusterImages(int runnum=247, 
			     TString path="/net/zwicky/dmtpc/data/4sh/mc/");


void plotMaxCamClusterImages(int runnum, 
			     TString path){

  TString runname = path+runNameOfNumber(runnum);
  cout << "runname = " << runname << endl;
  //  TString skimname = skimNameOfRun(runname);
  //  cout << "skimname = " << skimname << endl;

  c1 = new TCanvas("c1", "c1",1000,400);
  c1->Divide(3,1);
  
  TFile fMC(runname);

  DmtpcDataset d;
  d.openRootFile(runname);

  TObjArray *array=new TObjArray();
  Simulation->SetBranchAddress("trueClusterArray",&array);


  for(Int_t j=0; j<Simulation->GetEntries(); ++j){

    Simulation->GetEvent(j);

    c1->cd(1);
    d.getEvent(j);
    d.event()->ccdData(0)->Draw("COLZ");
    
    c1->cd(2);
    cout << "array->GetEntries()=" << array->GetEntries() << endl;
    MaxCamClusterImage *obj = (MaxCamClusterImage*)array->At(0);

    // draw the true image
    obj->getImage()->Draw("COLZ");

    // draw the true cluster
    c1->cd(3);
    obj->drawCluster(0,1);

    c1->Update();
    getchar();
  }
}
