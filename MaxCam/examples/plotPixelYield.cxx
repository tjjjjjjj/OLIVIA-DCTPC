
TH1F *hwall=0;
DmtpcDataset *d=0;
TF1 *fun;
TString fname="none";


vector<float> pTime, pNoise, pHot, pBias; 



void plot(TString what="", TString fname="test.root", TString opt="") {

  fun = new TF1("fun","gaus", 0, 20);

  if (!d) {
    d= new DmtpcDataset;
    d->openRootFile(fname);
  }

  TH2F *bias=d->getBiasFrame(1);

  what.ToLower();

  hwall= new TH1F("hwall","",200, 0, 20);
  
  int n=d->tree()->GetEntries();
  for (int i=0; i<n; i++) {
    d->tree()->GetEvent(i);
    int nwf=d->event()->ccdData()->GetEntries();

    for (int j=0; j<nwf; j++) {
      pTime.push_back( d->event()->ccdConfig(0)->exposureTime*1e-3 );
      TH2F *image=d->event()->ccdData(j);
      image->Add( bias, -1);
      TH1F* yield=MaxCamImageTools::createYieldHisto( image, -200, 200 );
      fun->SetParameters( yield->GetMaximum(), yield->GetMean(), yield->GetRMS() );
      yield->Fit("fun","LQ");
      //yield->Draw(); c1->Update(); getchar();
      hwall->Fill( yield->GetFunction("fun")->GetParameter(2) );
      pNoise.push_back(  yield->GetFunction("fun")->GetParameter(2) );
      pBias.push_back(  yield->GetFunction("fun")->GetParameter(1) );
      float th=yield->GetFunction("fun")->GetParameter(2)*5 + yield->GetFunction("fun")->GetParameter(1);
      delete yield;
      yield=MaxCamImageTools::createYieldHisto( image, th, 65536 );
      pHot.push_back(  yield->Integral() );
      delete yield;
    }

  }
  
  //hwall->DrawCopy(opt);

  cout << "exposure="<<d->event()->ccdConfig(0)->exposureTime 
       << "msec  sigma=" << hwall->GetMean() << "+/-" << hwall->GetRMS() 
       << endl;
  delete hwall;
  delete d;
  d=0;
}


void ALTA_080133() {
  plot("","data/dmtpc_run00530.root");
  plot("","data/dmtpc_run00531.root");
  plot("","data/dmtpc_run00532.root");
  plot("","data/dmtpc_run00533.root");
  plot("","data/dmtpc_run00534.root");
  plot("","data/dmtpc_run00535.root");
  plot("","data/dmtpc_run00536.root");
  plot("","data/dmtpc_run00537.root");
  plot("","data/dmtpc_run00538.root");
  plot("","data/dmtpc_run00539.root");
  plot("","data/dmtpc_run00540.root");
}

void ALTA_A80333() {
  plot("","data/dmtpc_run00550.root");
  plot("","data/dmtpc_run00551.root");
  plot("","data/dmtpc_run00552.root");
  plot("","data/dmtpc_run00553.root");
  plot("","data/dmtpc_run00554.root");
  plot("","data/dmtpc_run00555.root");
  plot("","data/dmtpc_run00556.root");
  plot("","data/dmtpc_run00557.root");
  plot("","data/dmtpc_run00558.root");
  plot("","data/dmtpc_run00559.root");
  plot("","data/dmtpc_run00560.root");
}

void ALTA_A80334() {
  plot("","data/dmtpc_run00563.root");
  plot("","data/dmtpc_run00564.root");
  plot("","data/dmtpc_run00565.root");
  plot("","data/dmtpc_run00566.root");
  plot("","data/dmtpc_run00567.root");
  plot("","data/dmtpc_run00568.root");
  plot("","data/dmtpc_run00569.root");
  plot("","data/dmtpc_run00570.root");
  plot("","data/dmtpc_run00571.root");
  plot("","data/dmtpc_run00572.root");
  plot("","data/dmtpc_run00573.root");
}


void allNoise() {
  //gStyle->SetOptFit(1);
  ALTA_080133();
  ALTA_A80333();
  ALTA_A80334();
  g1=new TGraph(110, &pTime[0], &pNoise[0]); g1->SetMarkerColor(2); 
  g2=new TGraph(110, &pTime[110], &pNoise[110]); g2->SetMarkerColor(4);
  g3=new TGraph(110, &pTime[220], &pNoise[220]);
  c1->Clear();
  g1->SetMinimum(5); g1->SetMaximum(15); g1->Draw("AP"); g2->Draw("P"); g3->Draw("P");
  fd=new TF1("fd","sqrt([0]**2+[1]*x)",0,100);
  g1->Fit("fd"); g2->Fit("fd"); g3->Fit("fd");
}


void gain(TString fname="data/dmtpc_run00575.root") {

  d= new DmtpcDataset;
  d->openRootFile(fname);
  
  for (int i=0; i<d->tree()->GetEntries(); i++) {
    d->getEvent(i);
    pTime.push_back( d->event()->ccdConfig(0)->exposureTime*1e-3 );
    TH1F* y=MaxCamImageTools::createYieldHisto(d->event()->ccdData(0));
    y->Fit("gaus","QL");
    pNoise.push_back(pow(y->GetFunction("gaus")->GetParameter(2),2));
    pBias.push_back(y->GetFunction("gaus")->GetParameter(1));
    delete y;
  }

  delete d; d=0;
}


void allGain() {
  gain("data/dmtpc_run00575.root");
  gain("data/dmtpc_run00576.root");
  g1=new TGraph(4, &pBias[0], &pNoise[0]); g1->Fit("pol1"); 
  g2=new TGraph(4, &pBias[9], &pNoise[9]); g2->Fit("pol1");
  fr=new TH2F("fr","",100,900, 1200, 100, 0, 300); fr->Draw(); g1->Draw("P"); g2->Draw("P");
}
