{
  TFile *f=new TFile("test.root");
  TH1D *hist=(TH1D*)f->Get(";1");

  hist->Draw();
}
