#define DCTPCTree_cxx
#include "../input_files/DCTPCTree.h"
#include "../input_files/Calibration_constants_LittleDCTPC.h";
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

void DCTPCTree::Loop()

{

TStopwatch timer;
  gROOT->SetStyle("Plain");
  gStyle->SetEndErrorSize(3);
  gStyle->SetPalette(1,0);
  gStyle->SetLineWidth(2);
  gStyle->SetHistLineWidth(2);
  gStyle->SetOptStat(kTRUE);
  gStyle->SetOptFit(kTRUE);
  TH1::AddDirectory(false);
  
  TFile *outtree = new TFile("$LittleDCTPC_physics_input");//("../input_files/allie.root");
  TTree *dctreepc = (TTree*)outtree->Get("dctpc_eventinfo");
  DCTPCTree aStep(dctreepc);

  ofstream seqcalibconstants;
  seqcalibconstants.open("Sequence_Calibration_Constants.txt");

  TH1D * hist_emesh=new TH1D("Emesh", "Emesh", 100, 2000, 10000);
  TH1D * hist_emeshseq_list[19];
  for(int seqnum=2;seqnum<=18;seqnum++)
  {
    hist_emeshseq_list[seqnum]=new TH1D("Emesh_seq", "Emesh_seq", 100, 2000, 10000);
  }

 
  //create event loop
  Long64_t nentries=dctreepc->GetEntries();
  cout<<"Number of Entries: "<<nentries<<endl;

  
  string filename;
  string titlename;
  string titlenum;
 for (int event = 0; event<nentries; event++)
    {
      aStep.GetEntry(event);
      if(event%1000000==0)
	cout<<((double)event/(double)nentries)*100.<<"%"<<endl;

      /////////////////////////////////////////////////////////////

      if(aStep.Emesh_kev<1000)
	{continue;}

      if(aStep.Mesh_max_v>0.38)
	{continue;}

      ////////////////////////////////////////////////////

      hist_emesh->Fill(aStep.Emesh_kev);
      hist_emeshseq_list[aStep.SequenceNum]->Fill(aStep.Emesh_kev);
    
    }

 hist_emesh->SetTitle("Mesh Energy for All Sequences");
 hist_emesh->SetXTitle("Mesh Energy kev");
 hist_emesh->SetYTitle("Counts");
 hist_emesh->Draw();
 c1->SaveAs("Mesh_Mean.gif");
 cout<<"Set mean: "<<hist_emesh->GetMean()<<endl;
 cout<<"Set mode: "<<hist_emesh->GetMaximumBin()<<endl;

 for(seqnum=2;seqnum<=18;seqnum++)
   {
     ostringstream convert;
     convert<<seqnum;
     filename = "mesh_";
     filename+ = convert.str();
     filename+=".gif";
     titlenum=convert.str();
     titlename="Sequence Mesh Energy Calibration - Sequence: ";
     titlename+=titlenum;
     hist_emeshseq_list[seqnum]->SetTitle(titlename.c_str());
     hist_emeshseq_list[seqnum]->SetXTitle("Mesh Energy/kev");
     hist_emeshseq_list[seqnum]->SetYTitle("Counts");
     hist_emeshseq_list[seqnum]->Draw();
     c1->SaveAs(filename.c_str());
     seqcalibconstants<<(hist_emesh->GetMean())/(hist_emeshseq_list[seqnum]->GetMean())<<" , ";
     cout<< "Sequence "<<seqnum<<" mean: "<<hist_emeshseq_list[seqnum]->GetMean()<<endl;
     cout<<"Sequence "<<seqnum<<" mode: "<<hist_emeshseq_list[seqnum]->GetMaximumBin()<<endl;

     cout<<"Set mean / Sequence "<<seqnum<< " mean = "<<(hist_emesh->GetMean())/(hist_emeshseq_list[seqnum]->GetMean())<<endl;
     cout<<"Sequence "<<seqnum<<" calibration constant = "<<(hist_emesh->GetMean())/(hist_emeshseq_list[seqnum]->GetMean())<<endl; //this is the number you mulitply all of the data in a sequence by to get the actual calibrated data (seq data/actual data = seq mean/actual mean)

   }

 seqcalibconstants.close();

}
