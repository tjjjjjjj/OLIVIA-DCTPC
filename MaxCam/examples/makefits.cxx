#include <iostream>
//#include "MaxCamImageTools.hh"

using std::cout;
using std::endl;

DmtpcDataset d;

int camNum;

int viewAll(int startNum=0) {

  TCanvas* c = new TCanvas();
  // get number of images
  int nImages = d.tree()->GetEntries();
  TString junk;
  for (int ii=startNum; ii<nImages; ii++) {
    cout << "Image Number " << ii << endl;
    view(ii);
    c->Update();
    cout << "continue?  (y|n)" << endl;
    cin >> junk;
    if (junk == "n") break;
  }

}

int init(TString fname="dmtpc_run00018.root", int ccdNum=0) {

  camNum = ccdNum;

  TString dir = "/export/data03/ddujmic/data/";

  // Denis says that in run 18 there are signals in events 46, 172 and 205
  // on ccd Number 1 (bottom)
  //TString fname = "dmtpc_run00018.root";   

  fname = dir+fname;

  cout << "opening " << fname << endl;
  d.openRootFile(fname);

  return 0;
}

float mean(int eventNumber) {
  d.tree()->GetEvent(eventNumber);
  float mean = d.event()->ccdData(camNum)->Integral()/(d.event()->ccdData(camNum)->GetNbinsX()*d.event()->ccdData(camNum)->GetNbinsY());
  cout << "mean = " << mean << endl;
  return mean;
}

float biasMean() {
  float mean = d.getBiasFrame(camNum)->Integral()/(d.getBiasFrame(camNum)->GetNbinsX()*d.getBiasFrame(camNum)->GetNbinsY());
  cout << "mean = " << mean << endl;
  return mean;
}


int view(int eventNumber, int useBias=0) {
  float mean = mean(eventNumber);
  d.tree()->GetEvent(eventNumber);

  float min = mean-60;
  float max = mean+60;
  cout << "min = " << min << endl;
  cout << "max = " << max << endl;
  d.event()->ccdData(camNum)->SetMinimum(min);
  d.event()->ccdData(camNum)->SetMaximum(max);

  d.event()->ccdData(camNum)->Draw("colz");
  return 0;
}

int fitsOf(int eventNumber, TString outfile) {

  // add extension to filename
  outfile = outfile+".fits";

  d.tree()->GetEvent(eventNumber);
  int status = MaxCamImageTools::convertIntoFits(d.event()->ccdData(camNum), outfile);
  cout << "status = " << status << endl;
}


void setCamNum(int cnum) {
  camNum = cnum;
}

int viewBias() {
  cout << "camNum = " << camNum << endl;

  float mean = biasMean();
  float min = mean-60;
  float max = mean+60;
  d.getBiasFrame(camNum)->SetMinimum(min);
  d.getBiasFrame(camNum)->SetMaximum(max);

  d.getBiasFrame(camNum)->Draw("colz");
}

int biasToFits(TString outfile="bias") {
  outfile = outfile + ".fits";
    
  int status = MaxCamImageTools::convertIntoFits(d.getBiasFrame(camNum), outfile);
  cout << "status = " << status << endl;
}
