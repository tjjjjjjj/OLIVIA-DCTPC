{
//  char * names[5] = {"v5", "v5bilinear", "v5rot1", "v4", "neighbor" }; 
//  char * files[5] = {"../cmp/stitched_neutrons.root", "../cmp/stitched_neutrons_bilinear.root", "../cmp/stitched_neutrons_rot1.root", "../cmp/stitched_neutrons_v4.root", "../cmp/stitched_neutrons_neighbor.root"}; 

  char * names[3] = {"v5", "v5new","v4"}; 
  char * files[3] = {"../cmp/stitched_neutrons.root", "../cmp/stitched_neutrons_new.root", "../cmp/stitched_neutrons_v4.root"}; 
   
  int ncanv = 18; 

  char * canvas_names[ncanv] = {"energy","range","phi","phienergy", "deltaE","deltaRange","npixel","width","nfound","cluster_rms","cluster_mean", "deltaEInteg", "deltaEInteg1cam", "npixel_ncam", "npixel_which_cam",  "npixel_mcrange", "npixel_deltaE_mcInteg", "npixel_mcInteg"}; 
  
  TCanvas * canvases[ncanv]; 

  for (int c = 0; c < ncanv; c++)
  {
    canvases[c] = new TCanvas(canvas_names[c], canvas_names[c], 800, 800); 
    canvases[c]->Divide(2,2); 
  }
   

  
  double * E[3]; 
  double * nfound[3]; 


  int n; 

  

  for (int i = 0; i < 3; i++)
  {
    TFile f(files[i]); 
    TTree * t = (TTree*) f.Get("cmp"); 
    t->SetLineColor(i+1); 
    //energy 
    canvases[0]->cd(i+1); 
    t->Draw("recoE:nfound","","NQ"); 
    E[i] = new double[t->GetEntries()]; 
    memcpy(E[i],t->GetV1(), t->GetEntries() * sizeof(double)); 
    nfound[i] = new double[t->GetEntries()]; 
    memcpy(nfound[i],t->GetV2(), t->GetEntries() * sizeof(double)); 
    n = t->GetEntries(); 
    t->Draw("recoE", "nfound == 1",""); 
    ((TH1*)gPad->GetPrimitive("htemp"))->SetTitle(names[i]); 
    canvases[0]->cd(6); 
    t->Draw("recoE", "nfound == 1",i==0 ? "" : "same"); 

    //range
    canvases[1]->cd(i+1); 
    t->Draw("mcRange:recoRange", "nfound == 1","colz"); 
    ((TH1*)gPad->GetPrimitive("htemp"))->SetTitle(names[i]); 

    //phi
    canvases[2]->cd(i+1); 
    t->Draw("mcPhi-recoPhi ", "nfound == 1",""); 
    ((TH1*)gPad->GetPrimitive("htemp"))->SetTitle(names[i]); 
    canvases[2]->cd(6); 
    t->Draw("mcPhi-recoPhi", "nfound == 1",i==0 ? "" : "same"); 
    
    //phienergy
    canvases[3]->cd(i+1); 
    t->Draw("mcPhi-recoPhi:mcE", "nfound == 1","colz"); 
    ((TH1*)gPad->GetPrimitive("htemp"))->SetTitle(names[i]); 

    //delta Range
    canvases[5]->cd(i+1); 
    t->Draw("mcRange-recoRange", "nfound == 1"); 
    ((TH1*)gPad->GetPrimitive("htemp"))->SetTitle(names[i]); 
    canvases[5]->cd(6); 
    t->Draw("mcRange-recoRange", "nfound == 1", i==0 ? "" : "same"); 
    //r 
    canvases[6]->cd(i+1); 
    t->Draw("npixel", "nfound == 1",""); 
    ((TH1*)gPad->GetPrimitive("htemp"))->SetTitle(names[i]); 
    canvases[6]->cd(6); 
    t->Draw("npixel", "nfound == 1", i==0 ? "" : "same"); 

    //width 
    canvases[7]->cd(i+1); 
    t->Draw("recoWidth", "nfound == 1",""); 
    ((TH1*)gPad->GetPrimitive("htemp"))->SetTitle(names[i]); 
    canvases[7]->cd(6); 
    t->Draw("recoWidth", "nfound == 1",i==0 ? "" : "same"); 

    //nfound 
    canvases[8]->cd(i+1); 
    t->Draw("nfound"); 
    ((TH1*)gPad->GetPrimitive("htemp"))->SetTitle(names[i]); 
    canvases[8]->cd(6); 
    t->Draw("nfound", "",i==0 ? "" : "same"); 
 
    //cluster_rms 
    canvases[9]->cd(i+1); 
    t->Draw("cluster_rms","nfound==1"); 
    ((TH1*)gPad->GetPrimitive("htemp"))->SetTitle(names[i]); 
    canvases[9]->cd(6); 
    t->Draw("cluster_rms", "nfound==1",i==0 ? "" : "same"); 
 
    //cluster_mean 
    canvases[10]->cd(i+1); 
    t->Draw("cluster_mean","nfound==1"); 
    ((TH1*)gPad->GetPrimitive("htemp"))->SetTitle(names[i]); 
    canvases[10]->cd(6); 
    t->Draw("cluster_mean", "nfound==1",i==0 ? "" : "same"); 

    //mcInteg 
    canvases[11]->cd(i+1); 
    t->Draw("recoE - mcInteg","nfound==1"); 
    ((TH1*)gPad->GetPrimitive("htemp"))->SetTitle(names[i]); 
    canvases[11]->cd(6); 
    t->Draw("recoE - mcInteg", "nfound==1",i==0 ? "" : "same"); 

    //mcInteg 
    canvases[12]->cd(i+1); 
    t->Draw("recoE - mcInteg","nfound==1 && ncam == 1"); 
    ((TH1*)gPad->GetPrimitive("htemp"))->SetTitle(names[i]); 
    canvases[12]->cd(6); 
    t->Draw("recoE - mcInteg", "nfound==1 && ncam == 1",i==0 ? "" : "same"); 

    //npixel_ncam 
    canvases[13]->cd(i+1); 
    t->Draw("npixel:ncam","nfound==1","colz"); 
    ((TH1*)gPad->GetPrimitive("htemp"))->SetTitle(names[i]); 

    //npixel_which_cam 
    canvases[14]->cd(i+1); 
    t->Draw("npixel:which_cam","nfound==1","colz"); 
    ((TH1*)gPad->GetPrimitive("htemp"))->SetTitle(names[i]); 

    //npixel range
    canvases[15]->cd(i+1); 
    t->Draw("npixel:mcRange","nfound==1 && ncam==1","colz"); 
    ((TH1*)gPad->GetPrimitive("htemp"))->SetTitle(names[i]); 

    //npixel deltaE
    canvases[16]->cd(i+1); 
    t->Draw("npixel:recoE-mcInteg","nfound==1 && ncam==1","colz"); 
    ((TH1*)gPad->GetPrimitive("htemp"))->SetTitle(names[i]); 
    
    //npixel E
    canvases[16]->cd(i+1); 
    t->Draw("npixel:mcInteg","nfound==1 && ncam==1","colz"); 
    ((TH1*)gPad->GetPrimitive("htemp"))->SetTitle(names[i]); 



  }


  for (int i = 0; i <5; i++)
  {
    
    canvases[4]->cd(i+1); 
    TLegend * leg = new TLegend(0.6,0.65,0.88,0.85); 
    for (int j = 0; j < 5; j++)
    {
      if (i==j) continue; 

      cout <<"here" <<endl; 
      TString thisname = TString(names[i]) +TString("-") +  TString(names[j]); 
//      sprintf(thisname,"delta%s_%s",names[i],names[j]); 
      TH1I * delta = new TH1I(thisname,thisname,100,-400,400); 
      for (int k = 0; k < n; k++)
      {
        if (nfound[i][k] != 1. || nfound[j][k] !=1.)
        {
          continue; 
        }

        delta->Fill(E[i][k] - E[j][k]); 
      }

      delta->SetLineColor(j+1); 
      if (j==0 || (i == 0 && j==1) )
      {
        delta->Draw(); 
      }
      else
      {
        delta->Draw("same"); 
      }
      leg->AddEntry(delta,thisname,"l"); 
    }

    leg->Draw(); 
  }


  for (int c = 0; c < ncanv; c++)
  {
    canvases[c]->Update(); 
  }


}
