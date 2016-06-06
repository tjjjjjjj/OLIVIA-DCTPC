#include "MultiVariate.hh"
#include <vector>
#include <string>
#include <iostream>
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TApplication.h"
#include <stdlib.h>
#include "style.hh"

#include "cosmic_runs.h"

/* Example program for using multivariate framework.  */ 

int main (int nargs, char ** args)
{

  int wr = 3; 
  int cam = 0 ; 
  int use_new = 0; 
  if (nargs > 1) cam = atoi(args[1]); 
  if (nargs > 2) wr = atoi(args[2]); 
  if (nargs > 3) use_new = atoi(args[3]); 
  
  TApplication app("app",NULL,NULL); 
  AnalysisStyle::setStyle(); 
  std::vector<std::string> signal; 
  std::vector<std::string> worms; 

  if (wr == 3)
  {
   worms.push_back("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01402skim.root"); 
   worms.push_back("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01403skim.root"); 
   worms.push_back("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01404skim.root"); 
   worms.push_back("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01405skim.root"); 
   worms.push_back("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01406skim.root"); 
   worms.push_back("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01407skim.root"); 
   worms.push_back("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01408skim.root"); 
   worms.push_back("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01409skim.root"); 
   worms.push_back("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01410skim.root"); 
  }

  else if (wr == 4)
  {
    worms.push_back("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_02613skim.root");
    worms.push_back("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_02614skim.root");
    worms.push_back("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_02615skim.root");
    worms.push_back("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_02616skim.root");
    worms.push_back("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_02617skim.root");
    worms.push_back("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_02618skim.root");
    worms.push_back("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_02619skim.root");
    worms.push_back("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_02620skim.root");
    worms.push_back("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_02621skim.root");
  }
  else if (wr == 5)
  {
    if (use_new)
    {
      WR5_TRAIN_NEW
    }
    else
    {
      WR5_TRAIN_OLD
    }
  }


  TTree * worms_ch = MultiVariate::buildMVTree("worms",&worms,2,1 << cam);


  std::cout << "worm tree built" << std::endl; 

  if (wr == 3)
  {
    TString base = "/net/zwicky/dmtpc/cozzyd/projects/DarkMatter/AnalysisFramework/v3/skim/dmtpc_run"; 
    for (int i = 101202001; i <=  101202018; i++)
    {
      TString f = base;
      f+=i; 
      f+="skim.root"; 
      signal.push_back(f.Data()); 
    }
    for (int i = 101203001; i <=  101203010; i++)
    {
      TString f = base;
      f+=i; 
      f+="skim.root"; 
      signal.push_back(f.Data()); 
    }
  }

  else if (wr == 4)
  {

    TString base = "/net/zwicky/dmtpc/cozzyd/projects/DarkMatter/AnalysisFramework/v3/skim/dmtpc_run"; 
    for (int i = 110110001; i <=  110110030; i++)
    {
      TString f = base;
      f+=i; 
      f+="skim.root"; 
      signal.push_back(f.Data()); 
    }


  }

  else if (wr == 5)
  {
    if (cam ==0)
    {
      signal.push_back("/net/mitdm01/data01/cozzyd/projects/DarkMatter/AnalysisFramework/v3/skim/dmtpc_mc_00005skim.root");
      signal.push_back("/net/mitdm01/data01/cozzyd/projects/DarkMatter/AnalysisFramework/v3/skim/dmtpc_mc_00006skim.root");
      signal.push_back("/net/mitdm01/data01/cozzyd/projects/DarkMatter/AnalysisFramework/v3/skim/dmtpc_mc_00007skim.root");
      signal.push_back("/net/mitdm01/data01/cozzyd/projects/DarkMatter/AnalysisFramework/v3/skim/dmtpc_mc_00008skim.root");
    }

    else if (cam==1)
    {
      signal.push_back("/net/mitdm01/data01/cozzyd/projects/DarkMatter/AnalysisFramework/v3/skim/dmtpc_mc_00013skim.root");
      signal.push_back("/net/mitdm01/data01/cozzyd/projects/DarkMatter/AnalysisFramework/v3/skim/dmtpc_mc_00014skim.root");
      signal.push_back("/net/mitdm01/data01/cozzyd/projects/DarkMatter/AnalysisFramework/v3/skim/dmtpc_mc_00015skim.root");
      signal.push_back("/net/mitdm01/data01/cozzyd/projects/DarkMatter/AnalysisFramework/v3/skim/dmtpc_mc_00016skim.root");
    }

  }

  //TTree * signal_ch = MultiVariate::buildMVTree("signal",&signal,0,1 << cam); 
  TTree * signal_ch = MultiVariate::buildMVTree("signal",&signal,0,1); 
  std::cout << "signal tree built" << std::endl; 
  std::cout << signal_ch->GetEntries() << std::endl; 

  TString name = "wr";
  name+= wr;
  name+="_worms_cam";
  name+=cam;
  if (use_new)
  {
    name+="_beta"; 
  }

  MultiVariate::runMVSimple(name.Data(),"/net/zwicky/dmtpc/cozzyd/wr_worms_mv.root",signal_ch,worms_ch);

  MultiVariate::MultiVariateResult * res = new MultiVariate::MultiVariateResult(); 

  res->setResult(name.Data()); 


  TChain worms_test("skim"); 
 
  if (wr==3) 
  {
    worms_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01411skim.root");
    worms_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01412skim.root");
    worms_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01413skim.root");
    worms_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01414skim.root");
    worms_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01415skim.root");
    worms_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01416skim.root");
    worms_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01417skim.root");
    worms_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01418skim.root");
    worms_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01419skim.root");
    worms_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01420skim.root");
  }

  if(wr==4)
  {
    worms_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_02622skim.root");
    worms_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_02623skim.root");
    worms_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_02624skim.root");
    worms_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_02625skim.root");
    worms_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_02626skim.root");
    worms_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_02627skim.root");
    worms_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_02628skim.root");
    worms_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_02629skim.root");
  }
  else if (wr==5)
  {
    if (use_new)
    {
      WR5_TEST_NEW
    }
    else
    {
      WR5_TEST_OLD
    }

  }

 

  TChain sig_test("skim"); 

  if (wr ==3)
  {
    TString base = "/net/zwicky/dmtpc/cozzyd/projects/DarkMatter/AnalysisFramework/v3/skim/dmtpc_run"; 
    for (int i = 101203011; i <=  101203032; i++)
      {
        TString f = base;
        f+=i; 
        f+="skim.root"; 
        sig_test.Add(f.Data()); 
     }
  }

  if (wr==4)
  {

    TString base = "/net/zwicky/dmtpc/cozzyd/projects/DarkMatter/AnalysisFramework/v3/skim/dmtpc_run"; 
    for (int i = 110110031; i <=  110110050; i++)
    {
      TString f = base;
      f+=i; 
      f+="skim.root"; 
      sig_test.Add(f.Data()); 
    }


  }

  else if (wr==5)
  {
    if (cam == 0)
    {
      sig_test.Add("/net/mitdm01/data01/cozzyd/projects/DarkMatter/AnalysisFramework/v3/skim/dmtpc_mc_00009skim.root");
      sig_test.Add("/net/mitdm01/data01/cozzyd/projects/DarkMatter/AnalysisFramework/v3/skim/dmtpc_mc_00010skim.root");
      sig_test.Add("/net/mitdm01/data01/cozzyd/projects/DarkMatter/AnalysisFramework/v3/skim/dmtpc_mc_00011skim.root");
      sig_test.Add("/net/mitdm01/data01/cozzyd/projects/DarkMatter/AnalysisFramework/v3/skim/dmtpc_mc_00012skim.root");
    }
    else if (cam == 1)
    {
      sig_test.Add("/net/mitdm01/data01/cozzyd/projects/DarkMatter/AnalysisFramework/v3/skim/dmtpc_mc_00017skim.root");
      sig_test.Add("/net/mitdm01/data01/cozzyd/projects/DarkMatter/AnalysisFramework/v3/skim/dmtpc_mc_00018skim.root");
      sig_test.Add("/net/mitdm01/data01/cozzyd/projects/DarkMatter/AnalysisFramework/v3/skim/dmtpc_mc_00019skim.root");
      sig_test.Add("/net/mitdm01/data01/cozzyd/projects/DarkMatter/AnalysisFramework/v3/skim/dmtpc_mc_00020skim.root");
    }
  }


  MultiVariate::MultiVariateEvaluation eval = MultiVariate::evaluate(res,&sig_test,&worms_test,0,cam,-1);
                                              

  TCanvas c3("eff","eff",600,600); 
  eval.efficiency()->SetLineColor(kBlue); 
  eval.efficiency()->Draw("al"); 

  eval.background()->SetLineColor(kRed); 
  eval.background()->Draw("same l"); 

  eval.purity()->SetLineColor(kGreen); 
  eval.purity()->Draw("same l"); 
  double cutoff = eval.cutoff(); 
  cout << "Cutoff: " << cutoff << endl;

/*

  TChain cf_test("skim"); 
  cf_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01277skim.root"); 
  cf_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01278skim.root"); 
  cf_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01279skim.root"); 
  cf_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01280skim.root"); 
  cf_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01281skim.root"); 
  cf_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01282skim.root"); 
  cf_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01283skim.root"); 
  cf_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01284skim.root"); 
  cf_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01285skim.root"); 
  cf_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01286skim.root"); 
  cf_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01287skim.root"); 
  cf_test.Add("/net/zwicky/dmtpc/production/skimout/skim/dmtpc_10L_01288skim.root"); 

  TH2F * cfre = new TH2F("cfre","cf252 r v e",100,0,10000,100,0,500);
  TH2F * cfre_cut = new TH2F("cfre_cut","cf252_cut r v e",100,0,10000,100,0,500);

  cf_test.SetBranchStatus("*clusters*",0); 
  cf_test.SetBranchStatus("*trigger*",0); 
  DmtpcSkimEvent *e = 0; 
  cf_test.SetBranchAddress("event",&e); 
  for (int i = 0; i < cf_test.GetEntries(); i++)
  {
    cf_test.GetEntry(i); 
    for (int t = 0; t < e->ntracks(cam); t++)
    {
      if (e->edge(cam,t) || e->nburnin(cam,t) > 2 || e->range(cam,t) <= 0) continue; 
      cfre->Fill(e->E(cam,t), e->range(cam,t));
      double c = res->getClassifier(e,cam,t); 
      if (c > cutoff) 
        cfre_cut->Fill(e->E(cam,t), e->range(cam,t));
    }
  }

  TCanvas cf("cf","cf",400,400); 
  cfre->Draw("colz"); 
  TCanvas cf_cut("cfcut","cfcut",400,400); 
  cfre_cut->Draw("colz"); 
  */
  app.Run(); 
}
