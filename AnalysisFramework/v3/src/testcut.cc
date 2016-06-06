#include "AnalysisCut.hh"
#include "../../../MaxCam/DmtpcSkimEvent.hh"
#include "../../../MaxCam/DmtpcSkimDataset.hh"
#include "TFile.h"
#include "TTree.h"
#include <ostream>
#include <iostream>
#include <vector>


bool cutfunction(DmtpcSkimEvent* ev, int c, int t, double* param)
{
   if((c+t) < 4*param[0])
      return 1;
   else
      return 0;

}

bool cuttfunction(DmtpcSkimEvent* ev, int c, int t, vector<void*> param)
{
   double par0 = *((double*)param[0]);
   if((c+t) < 4*par0)
      return 1;
   else
      return 0;
}


int main()
{
   AnalysisCut* cut = new AnalysisCut("testcut",cutfunction,1);
   cut->setParameter(0,2);
   AnalysisCut* cutt = new AnalysisCut("testcutt",cuttfunction);
   vector<void*> cuttparam;
   double par0 = 2;
   cuttparam.push_back(&par0);
   cutt->setParameters(cuttparam);

   TFile* outfile = new TFile("testcut.root","RECREATE");
   TTree* outtree = new TTree("passes","passes");
   outtree->Branch("testcut","AnalysisCut",&cut,128000,1);
   outtree->Branch("testcutt","AnalysisCut",&cutt,128000,1);

   DmtpcSkimDataset d;
   d.openRootFile("/net/zwicky/dmtpc/analysis/10L/skim/dmtpc_10L_02000skim.root");

   for(int i=0; i<d.nevents(); i++)
   {
      if(i%50==0) cout << "Event: " << i << endl;
      d.getEvent(i);
      cut->evaluateEvent(d.event());
      cutt->evaluateEvent(d.event());
      outfile->cd();
      outtree->Fill();
   }

   outtree->Write();
   outfile->Close();

   return 0;
}
