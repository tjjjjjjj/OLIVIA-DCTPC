{ 

  int nbins = 100;
  TH2F *h1 = new TH2F("name", "title", nbins, 0, 100, nbins, 0, 100);
  for (int ii=1; ii<=nbins; ii++) {
    for (int jj=1; jj<=nbins; jj++) {
      h1->SetBinContent(ii, jj, ii*jj);
    }
  }

  TCanvas *cc = new TCanvas("cc");
  cc->Divide(2,2);
  cc->cd(1);
  h1->DrawCopy("colz");

  TString dir;

  TH2F *htemp;

  cc->cd(2);
  htemp = MaxCamImageTools::rotateRight(h1);
  htemp->DrawCopy("colz");

  cc->cd(3);
  htemp = MaxCamImageTools::rotateLeft(h1);
  htemp->DrawCopy("colz");

  cc->cd(4);
  htemp = MaxCamImageTools::rotate180(h1);
  htemp->DrawCopy("colz");

}
