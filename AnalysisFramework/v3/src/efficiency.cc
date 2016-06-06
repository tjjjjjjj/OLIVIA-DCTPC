#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TChain.h"
#include "../../../MaxCam/DmtpcSkimEvent.hh"
#include "recoilEnergy.hh"
#include "TApplication.h"

int main (int nargs, char ** args)
{
  TApplication * app = new TApplication("blah",0,0);  

  TChain * cam0pass = new TChain("skim"); 
  TChain * cam1pass = new TChain("skim"); 
  cam0pass->Add("pass/100421021.root"); 
  cam0pass->Add("pass/100421022.root"); 
  cam1pass->Add("pass/100421023.root"); 
  cam1pass->Add("pass/100421024.root"); 
  
  
  TChain * cam0all = new TChain("skim"); 
  TChain * cam1all = new TChain("skim"); 

  cam0all->SetBranchStatus("*cluster*",0);
  cam0all->SetBranchStatus("*orig_event*",0);
  cam1all->SetBranchStatus("*cluster*",0);
  cam1all->SetBranchStatus("*orig_event*",0);
  cam0pass->SetBranchStatus("*cluster*",0);
  cam0pass->SetBranchStatus("*orig_event*",0);
  cam1pass->SetBranchStatus("*cluster*",0);
  cam1pass->SetBranchStatus("*orig_event*",0);

  cam0all->Add("skim/dmtpc_run100421021skim.root"); 
  cam0all->Add("skim/dmtpc_run100421022skim.root"); 
  cam1all->Add("skim/dmtpc_run100421023skim.root"); 
  cam1all->Add("skim/dmtpc_run100421024skim.root"); 

  DmtpcSkimEvent * ev = new DmtpcSkimEvent; 
  int track;
  cam0pass->SetBranchAddress("event",&ev); 
  cam1pass->SetBranchAddress("event",&ev); 
  cam0pass->SetBranchAddress("track",&track); 
  cam1pass->SetBranchAddress("track",&track); 
  cam0all->SetBranchAddress("event",&ev); 
  cam1all->SetBranchAddress("event",&ev); 

  TH1F * pass0E = new TH1F("pass0E","pass0E", 40, 0, 200); 
  TH1F * pass1E = new TH1F("pass1E","pass1E", 40, 0, 200); 
  TH1F * all0E = new TH1F("all0E","all0E", 40, 0, 200); 
  TH1F * all1E = new TH1F("all1E","all1E", 40, 0, 200); 

  
  std::cout << "1" <<std::endl;
  for (int i = 0; i < cam0pass->GetEntries(); i++)
  {
    cam0pass->GetEntry(i); 
    pass0E->Fill(RecoilEnergy::getRecoilEnergy(ev->E(0,track),0)); 
  }

  std::cout << "2" <<std::endl;

  for (int i = 0; i < cam1pass->GetEntries(); i++)
  {
    cam1pass->GetEntry(i); 
    pass1E->Fill(RecoilEnergy::getRecoilEnergy(ev->E(0,track),1)); 
  }

  std::cout << "3" <<std::endl;
  for (int i = 0; i < cam0all->GetEntries(); i++)
  {
    cam0all->GetEntry(i); 
    for (int t = 0; t < ev->ntracks(0); t++)
      all0E->Fill(RecoilEnergy::getRecoilEnergy(ev->E(0,t),0)); 
  }

  std::cout << "4" <<std::endl;
  for (int i = 0; i < cam1all->GetEntries(); i++)
  {
    cam1all->GetEntry(i); 
    for (int t = 0; t < ev->ntracks(0); t++)
      all1E->Fill(RecoilEnergy::getRecoilEnergy(ev->E(0,t),1));
  }

  all0E->Sumw2(); 
  all1E->Sumw2(); 
  pass0E->Sumw2(); 
  pass1E->Sumw2(); 

  TH1F * eff0 = (TH1F*) pass0E->Clone("eff0");
  eff0->SetTitle("Camera 0 Efficiency");
  eff0->SetXTitle("Recoil Energy (keV)");
  TH1F * eff1 = (TH1F*) pass1E->Clone("eff1");
  eff1->SetTitle("Camera 1 Efficiency");
  eff1->SetXTitle("Recoil Energy (keV)");
  
  eff0->Divide(all0E); 
  eff1->Divide(all1E); 


  TH1F * effAvg = (TH1F*)eff0->Clone("eff");
  effAvg->SetTitle("Average Efficiency");
  effAvg->SetXTitle("Recoil Energy (keV)");

  effAvg->Add(eff1);
  effAvg->Scale(0.5); 

  TCanvas * c = new TCanvas("c","c",1500,500);
  c->Divide(3,1); 

  c->cd(1);
  eff0->Draw();

  c->cd(2);
  eff1->Draw();

  c->cd(3); 
  effAvg->Draw(); 
  c->Update();
  c->SaveAs("mc_eff.pdf"); 
 
  TFile f("efficiency.root","RECREATE"); 

  eff0->Clone()->Write();
  eff1->Clone()->Write();
  effAvg->Clone()->Write();

  f.Close(); 

  app->Run(); 

}
