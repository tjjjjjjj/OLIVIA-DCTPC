{ 

  cout << "overscan.cxx" << endl;

  MaxCamAlta *cam = new MaxCamAlta(1);
  MaxCamAlta::discoverCameras();
  cam->openCamera(1);
  //cam->reset();
  cam->resetWithFlush();

  cam->setExposureTime(500); // ms
  cam->setDigitizeOverscan(true);
  //cam->setDigitizeOverscan(false);

  Int_t pixPerBin = 4;
  cam->setHBin(pixPerBin);
  cam->setVBin(pixPerBin);

  int nimg = 1;
  TCanvas *c1 = new TCanvas("c1", "", 1000, 450);
  c1->Divide(2);
  for (int ii=0; ii<nimg; ii++) {
    cam->flushCamera(); 
    cam->openShutter();
    cam->expose();
    cam->closeShutter();
    cam->grabImage();
  
    printf("\n\n\n createHisto()\n");
    c1->cd(1);
    TH2F* mainImg = cam->createHisto("myhist");
    mainImg->DrawCopy("colz");
    //cout << "   nbins x, y = " << mainImg->GetNbinsX() << ", " << mainImg->GetNbinsY() << endl;

    c1->cd(1);
    TH2F* fullImg = cam->createFullHisto("myfullhist");
    fullImg->DrawCopy("colz");

    c1->cd(2);
    TH2F* overImg = cam->createOverscanHisto("junk");
    //overImg->SetMinimum(1190);
    //overImg->SetMaximum(1260);
    overImg->DrawCopy("colz");
    cout << "overImg->GetMinimum(), Maximum() = " << overImg->GetMinimum() << ", " << overImg->GetMaximum() << endl;

    c1->Update();
    //gSystem->Sleep(500);
  }

//  cout << "cam->getConfiguration()->meanForExposedPixels = " << cam->getConfiguration()->meanForExposedPixels << endl;
//  cout << "cam->getConfiguration()->rmsForExposedPixels = " << cam->getConfiguration()->rmsForExposedPixels << endl;
//  cout << "MaxCamImageTools::getMean(mainImg) = " << MaxCamImageTools::getMean(mainImg) << endl;
//  cout << "MaxCamImageTools::getRMS(mainImg)  = " << MaxCamImageTools::getRMS(mainImg) << endl;
//
//  cout << "cam->getConfiguration()->meanForOverscanPixels = " << cam->getConfiguration()->meanForOverscanPixels << endl;
//  cout << "cam->getConfiguration()->rmsForOverscanPixels = " << cam->getConfiguration()->rmsForOverscanPixels << endl;
//  cout << "MaxCamImageTools::getMean(overImg) = " << MaxCamImageTools::getMean(overImg) << endl;
//  cout << "MaxCamImageTools::getRMS(overImg)  = " << MaxCamImageTools::getRMS(overImg) << endl;
//
//
//  cout << "cam->getConfiguration()->ul_x, ul_y = " << cam->getConfiguration()->ul_x << ", " << cam->getConfiguration()->ul_y << endl;
//  cout << "cam->getConfiguration()->lr_x, lr_y = " << cam->getConfiguration()->lr_x << ", " << cam->getConfiguration()->lr_y << endl;
//  cout << "cam->getConfiguration()->hbin, vbin = " <<   cam->getConfiguration()->hbin << ", " << cam->getConfiguration()->vbin << endl;
//  cout << "cam->getConfiguration()->row_width, img_rows = " <<   cam->getConfiguration()->row_width << ", " << cam->getConfiguration()->img_rows << endl;
//  cout << "cam->getConfiguration()->serialNumber = " <<   cam->getConfiguration()->serialNumber << endl;
//  cout << "cam->getConfiguration()->overscanColumns = " <<   cam->getConfiguration()->overscanColumns << endl;
//  cout << "cam->getConfiguration()->overscanRows = " <<   cam->getConfiguration()->overscanRows << endl;

  cam->closeCamera();
  delete cam;
  delete overImg;
}

