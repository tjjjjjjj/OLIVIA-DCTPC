#include "DmtpcDataset.hh" 
#include <vector>
#include "MaxCamImageTools_Template.hh" 
#include "DmtpcRootTools.hh" 
#include "TCanvas.h"
#include "DmtpcIterators.hh" 
#include "TApplication.h"


int main(int nargs, char ** args)
{

  DmtpcRootTools::setColorStandard1(); 
  TApplication app("app",0,0); 
  DmtpcDataset d; 
  d.openRootFile(args[1]); 
  int cam = atoi(args[2]); 

  int n = d.tree()->GetEntries(); 
  DmtpcDatasetImageIterator beg(&d,cam); 
  DmtpcDatasetImageIterator end(&d,cam,n); 

  TF1 * gaus = new TF1("fitfn","gaus",0,65535); 
  vector<TH2*> *fits = MaxCamImageTools::valueFit<float> (beg,end,n,gaus,"N",1<<20,true); 

  TCanvas c; 
  c.Divide(2,2); 
  c.cd(1); 
  (*fits)[0]->Draw("colz"); 
  c.cd(2); 
  (*fits)[1]->Draw("colz"); 
  c.cd(3); 
  (*fits)[2]->Draw("colz"); 
  c.cd(4); 
  (*fits)[3]->Draw("colz"); 
  app.Run(); 

}
