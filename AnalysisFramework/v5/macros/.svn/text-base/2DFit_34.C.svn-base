{

  gSystem->Load("/net/fritz/data02/cozzyd/slide_maker/libSlideMaker.so"); 
  gSystem->Load("/net/zwicky/dmtpc/cozzyd/projects/DarkMatter/MaxCam/libMaxCam.so"); 
  gROOT->LoadMacro("phiVsX.C"); 


  DmtpcRootTools::setColorStandard1(); 

  gStyle->SetNumberContours(99); 


  char * out_prefix = "2DFit_34/"; 

  BeamerSlideMakerSettings slide_settings; 
  slide_settings->author="Cosmin Deaconu" ; 
  slide_settings->title="2D Fitting"; 
  slide_settings->short_title="2D Fitting"; 
  slide_settings->style="Boadilla"; 
  slide_settings->color_style="beaver"; 
  slide_settings->navigation_symbols=false; 
  slide_settings->institute="MIT"; 
  slide_settings->title_page = false; 

  BeamerSlideMaker slides(out_prefix,"2DFitting", &slide_settings); 


  TChain * fit = new TChain("fit"); 
  fit->Add("../2dfit_34_new.root"); 

  TF1 line("line","pol1",0,1e4); 
  line.SetParameters(0,1); 



  /// chisq /// 
  TCanvas * chisqc = new TCanvas("chisqc","chisqc",800,600); 
  TH1I *chisqh = new TH1I ("chisqh","chisqh",50,0, 5); 
  fit->Draw("chisq/ndof >> chisqh"); 
  chisqh->GetXaxis()->SetTitle("chi^2/ndof"); 
  chisqh->SetLineColor(2); 
  chisqh->DrawCopy(); 
  slides.AddPlotSlide("$\\chi^2/ndof$","chisq",chisqc,"pdf","width=4.5in"); 
 

  /// Energy comparison ///
  TCanvas * energy_comparison = new TCanvas("energycomparison","energycomparison",800,600); 
  TH2I * energyh = new TH2I("energyh","energyh",100,0,5000,100,0,5000); 
  fit->Draw("mcE: E.val >> energyh","","colz"); 
  energyh->GetYaxis()->SetTitle("MC Evis (adu)"); 
  energyh->GetXaxis()->SetTitle("Fit Evis (adu)"); 
  energyh->SetTitle("fit E reconstruction"); 
  line.Draw("same"); 
  slides.AddPlotSlide("Fit Energy Comparison","energyc",energy_comparison,"pdf","width=4.5in"); 

  /// Energy comparison2 ///
  TCanvas * energy_comparison2 = new TCanvas("energycomparison2","energycomparison2",800,600); 
  TH2I * energyh2 = new TH2I("energyh2","energyh2",30,0,300,60,-0.5,0.5); 
  fit->Draw(" (E.val - mcE)/mcE : mcEr >> energyh2","","colz"); 
  energyh2->GetXaxis()->SetTitle("MC E Recoil (keV)"); 
  energyh2->SetTitle("E Resolution"); 
  energyh2->GetYaxis()->SetTitle("(Fit Evis - MC E vis) / (MC E vis)"); 
  slides.AddPlotSlide("Fit Energy Resolution","energyc2",energy_comparison2,"pdf","width=4.5in"); 


  /// Reco Energy comparison ///
  TCanvas * reco_energy_comparison = new TCanvas("recoenergycomparison","recoenergycomparison",800,600); 
  TH2I * recoenergyh = new TH2I("recoenergyh","recoenergyh",100,0,5000,100,0,5000); 
  fit->Draw("mcE : recoE >> recoenergyh","","colz"); 
  recoenergyh->GetYaxis()->SetTitle("MC Evis (adu)"); 
  recoenergyh->GetXaxis()->SetTitle("Reco Evis (adu)"); 
  recoenergyh->SetTitle("E reco"); 
  line.Draw("same"); 
  slides.AddPlotSlide("Reco Energy Comparison","recoenergyc",reco_energy_comparison,"pdf","width=4.5in"); 

  /// range comparison /// 
  TCanvas * range_comparison = new TCanvas("rangecomparison","rangecomparison",800,600); 
  TH2I * rangeh = new TH2I("rangeh","rangeh",80,0,80,80,0,80); 
  fit->Draw("mcRange * sin(mcTheta) : range.val >> rangeh","","colz"); 
  rangeh->GetYaxis()->SetTitle("MC Range (pix)"); 
  rangeh->GetXaxis()->SetTitle("Fit Range (pix)"); 
  rangeh->SetTitle("fit Range reconstruction"); 
  line.Draw("same"); 
  slides.AddPlotSlide("Fit range Comparison","rangec",range_comparison,"pdf","width=4.5in"); 

  /// range comparison2 ///
  TCanvas * range_comparison2 = new TCanvas("rangecomparison2","rangecomparison2",800,600); 
  TH2I * rangeh2 = new TH2I("rangeh2","rangeh2",60,-0.5,0.5,80,0,80); 
  fit->Draw(" mcRange* sin(mcTheta) : (range.val - mcRange * sin(mcTheta))/(mcRange * sin(mcTheta)) >> rangeh2","","colz"); 
  rangeh2->GetYaxis()->SetTitle("MC Range (pix)"); 
  rangeh2->SetTitle("Range Resolution"); 
  rangeh2->GetXaxis()->SetTitle("(Fit Range - MC Range) / (MC Range)"); 
  slides.AddPlotSlide("Fit range Resolution","rangec2",range_comparison2,"pdf","width=4.5in"); 

  /// Reco range comparison ///
  TCanvas * reco_range_comparison = new TCanvas("recorangecomparison","recorangecomparison",800,600); 
  TH2I * recorangeh = new TH2I("recorangeh","recorangeh",80,0,80,80,0,80); 
  fit->Draw(" mcRange * sin(mcTheta) : recoRange >> recorangeh","","colz"); 
  recorangeh->GetYaxis()->SetTitle("MC Range (pix)"); 
  recorangeh->GetXaxis()->SetTitle("Reco Range (pix)"); 
  recorangeh->SetTitle("Range reco"); 
  line.Draw("same"); 
  slides.AddPlotSlide("Reco Range Comparison","recorangec",reco_range_comparison,"pdf","width=4.5in"); 

  /// sigma 
  TCanvas * sigma_comparison = new TCanvas("sigmacomparison","sigmacomparison",800,600); 
  TH2I * sigmah = new TH2I("sigmah","sigmah",25,0,100,75,0,300); 
  fit->Draw("mcZ : pow(sigma.val,2)  >> sigmah","","colz"); 
  sigmah->GetYaxis()->SetTitle("MC Z (mm)"); 
  sigmah->GetXaxis()->SetTitle("Fit #sigma^2 (pix^2)"); 
  sigmah->SetTitle("fit sigma reconstruction"); 
  slides.AddPlotSlide("Fit sigma Comparison","sigmac",sigma_comparison,"pdf","width=4.5in"); 

  /// sigmarange 
  TCanvas * sigmarange_comparison = new TCanvas("sigmarangecomparison","sigmarangecomparison",800,600); 
  TH2I * sigmarangeh = new TH2I("sigmarangeh","sigmarangeh",25,0,12.5,30,0,60); 
  fit->Draw("range.val : sigma.val  >> sigmarangeh","","colz"); 
  sigmarangeh->GetYaxis()->SetTitle("Fit Range (pix)"); 
  sigmarangeh->GetXaxis()->SetTitle("Fit #sigma (pix)"); 
  sigmarangeh->SetTitle("fit sigmarange reconstruction"); 
  double corr = sigmarangeh->GetCorrelationFactor(); 
  slides.AddPlotSlide(TString::Format("Fit Range vs. $\\sigma$ (corr =%f)",corr),"sigmarangec",sigmarange_comparison,"pdf","width=4.5in"); 

  /// sigmachisq 
  TCanvas * sigmachisq_comparison = new TCanvas("sigmachisqcomparison","sigmachisqcomparison",800,600); 
  TH2I * sigmachisqh = new TH2I("sigmachisqh","sigmachisqh",25,0,12.5,30,0,6); 
  fit->Draw("chisq/ndof : sigma.val  >> sigmachisqh","","colz"); 
  sigmachisqh->GetYaxis()->SetTitle("#chi^{2}/ndof "); 
  sigmachisqh->GetXaxis()->SetTitle("Fit #sigma (pix)"); 
  sigmachisqh->SetTitle("fit sigmachisq reconstruction"); 
  slides.AddPlotSlide("Fit $\\chi^2$ vs. $\\sigma$ ","sigmachisqc",sigmachisq_comparison,"pdf","width=4.5in"); 



  /// slope  vs sin theta
//  TCanvas * slope = new TCanvas("slope","slope",800,600); 
 // TH2I * slopeh = new TH2I("slopeh","slopeh",30,0,3,30,0,60); 
  //fit->Draw(" abs(delta_z.val) / range.val : 1/sin(mcTheta) >> slopeh","","colz"); 
 // slopeh->GetYaxis()->SetTitle("slope (adu/pix^2)"); 
 // slopeh->GetXaxis()->SetTitle("sin (MC #theta)"); 
 // slopeh->SetTitle("slope"); 
 // slides.AddPlotSlide("Slope","slope",slope); 

  /// raw phi ///
  TCanvas * rawPhi = new TCanvas("rawphi","rawphi",800,600); 
  TH1I *rawphi_h = new TH1I ("rawphi_h","rawphi_h",60,-TMath::Pi(), TMath::Pi()); 
  fit->Draw("atan2(sin(htphi.val-mcPhi),cos(htphi.val-mcPhi)) >> rawphi_h"); 
  rawphi_h->SetLineColor(4); 
  rawphi_h->SetTitle("#phi difference"); 
  rawphi_h->DrawCopy(); 
  rawphi_h->SetLineColor(2); 
  fit->Draw("atan2(sin(phi2-mcPhi),cos(phi2-mcPhi)) >> rawphi_h","","same"); 
  rawphi_h->DrawCopy("same" ); 
  slides.AddPlotSlide("$\\phi$ , (blue is fit, red is getPhi2)","rawPhi", rawPhi); 


  /// phi Energy ///
  TH2I *phienergy_h = new TH2I ("phienergy_h","#phi vs. Energy", 15,0,300, 36,-TMath::Pi(), TMath::Pi()); 
  fit->Draw("atan2(sin(htphi.val-mcPhi),cos(htphi.val-mcPhi)) : mcEr >> phienergy_h","","goff"); 
  phienergy_h->GetXaxis()->SetTitle("MC Erecoil (keV)"); 
  TCanvas * phiEnergy = phiVsX(phienergy_h); 
  slides.AddPlotSlide("fit $\\phi$ vs. energy","phiEnergy", phiEnergy,"pdf","width=4.5in"); 

  /// phi Range ///
  TH2I *phirange_h = new TH2I ("phirange_h","#phi vs. MC Range", 15,0,60, 36,-TMath::Pi(), TMath::Pi()); 
  fit->Draw("atan2(sin(htphi.val-mcPhi),cos(htphi.val-mcPhi)) : mcRange * sin(mcTheta) >> phirange_h","","goff"); 
  phirange_h->GetXaxis()->SetTitle("MC Range (px)"); 
  TCanvas * phirange = phiVsX(phirange_h); 
  slides.AddPlotSlide("fit $\\phi$ vs. MC Range","phirange", phirange,"pdf","width=4.5in"); 

  /// phi2 Energy ///
  TH2I *phi2energy_h = new TH2I ("phi2energy_h","phi2 vs. Energy", 15,0,300, 36,-TMath::Pi(), TMath::Pi()); 
  fit->Draw("atan2(sin(phi2-mcPhi),cos(phi2.val-mcPhi)) : mcEr >> phi2energy_h","","goff"); 
  phi2energy_h->GetXaxis()->SetTitle("MC Erecoil (keV)"); 
  TCanvas * phi2Energy = phiVsX(phi2energy_h); 
  slides.AddPlotSlide("Reco (GetPhi2) $\\phi$ vs. energy","phi2Energy", phi2Energy,"pdf","width=4.5in"); 

  /// phi2 Range ///
  TH2I *phi2range_h = new TH2I ("phi2range_h","phi2 vs. MC Range", 15,0,60, 36,-TMath::Pi(), TMath::Pi()); 
  fit->Draw("atan2(sin(phi2.val-mcPhi),cos(phi2.val-mcPhi)) : mcRange * sin(mcTheta) >> phi2range_h","","goff"); 
  phi2range_h->GetXaxis()->SetTitle("MC Range (px)"); 
  TCanvas * phi2range = phiVsX(phi2range_h); 
  slides.AddPlotSlide("Reco (GetPhi2) $\\phi$ vs. MC Range","phi2range", phi2range,"pdf","width=4.5in"); 


  /// phi Range ///
  TH2I *phifitrange_h = new TH2I ("phifitrange_h","#phi vs. Fit Range", 15,0,60, 36,-TMath::Pi(), TMath::Pi()); 
  fit->Draw("atan2(sin(htphi.val-mcPhi),cos(htphi.val-mcPhi)) : range.val >> phifitrange_h","","goff"); 
  phifitrange_h->GetXaxis()->SetTitle("Fit Range (px)"); 
  TCanvas * phifitrange = phiVsX(phifitrange_h); 
  slides.AddPlotSlide("fit $\\phi$ vs. Fit Range","phifitrange", phifitrange,"pdf","width=4.5in"); 


  /// phi prob ///
  TH2I *phiprob_h = new TH2I ("phiprob_h","#phi vs. prob", 15,0,1.2, 36,-TMath::Pi(), TMath::Pi()); 
  fit->Draw("atan2(sin(htphi.val-mcPhi),cos(htphi.val-mcPhi)) : prob  >> phiprob_h","","goff"); 
  phiprob_h->GetXaxis()->SetTitle("Prob"); 
  TCanvas * phiprob = phiVsX(phiprob_h); 
  slides.AddPlotSlide("fit $\\phi$ vs. prob","phiprob", phiprob,"pdf","width=4.5in"); 

  /// phi chisq ///
  TH2I *phichisq_h = new TH2I ("phichisq_h","#phi vs. #chi^{2} / ndof", 15,0,3, 36,-TMath::Pi(), TMath::Pi()); 
  fit->Draw("atan2(sin(htphi.val-mcPhi),cos(htphi.val-mcPhi)) : chisq/ndof  >> phichisq_h","","goff"); 
  phichisq_h->GetXaxis()->SetTitle("#chi^{2} /ndof"); 
  TCanvas * phichisq = phiVsX(phichisq_h); 
  slides.AddPlotSlide("fit $\\phi$ vs. $\\chi^2$","phichisq", phichisq,"pdf","width=4.5in"); 

  /// phi range/sigma ///
  TH2I *phiroversigma_h = new TH2I ("phiroversigma_h","#phi vs. range / #sigma", 20,0,10, 36,-TMath::Pi(), TMath::Pi()); 
  fit->Draw("atan2(sin(htphi.val-mcPhi),cos(htphi.val-mcPhi)) : range.val / sigma.val  >> phiroversigma_h","","goff"); 
  phiroversigma_h->GetXaxis()->SetTitle("range / #sigma"); 
  TCanvas * phiroversigma = phiVsX(phiroversigma_h); 
  slides.AddPlotSlide("fit $\\phi$ vs. fit range/ fit $\\sigma$","phiroversigma", phiroversigma,"pdf","width=4.5in"); 

  /// phi theta ///
  TH2I *phitheta_h = new TH2I ("phitheta_h","#phi vs. MC  #theta", 15,0,180, 36,-TMath::Pi(), TMath::Pi()); 
  fit->Draw("atan2(sin(htphi.val-mcPhi),cos(htphi.val-mcPhi)) : mcTheta * 180 / TMath::Pi()  >> phitheta_h","","goff"); 
  phitheta_h->GetXaxis()->SetTitle("MC #theta"); 
  TCanvas * phitheta = phiVsX(phitheta_h); 
  slides.AddPlotSlide("fit $\\phi$ vs. theta","phitheta", phitheta,"pdf","width=4.5in"); 

  /// phi sigma ///
  TH2I *phisigma_h = new TH2I ("phisigma_h","#phi vs. #sigma", 15,0,15, 36,-TMath::Pi(), TMath::Pi()); 
  fit->Draw("atan2(sin(htphi.val-mcPhi),cos(htphi.val-mcPhi)) : sigma.val  >> phisigma_h","","goff"); 
  phisigma_h->GetXaxis()->SetTitle("#sigma (px)"); 
  TCanvas * phisigma = phiVsX(phisigma_h); 
  slides.AddPlotSlide("fit $\\phi$ vs. sigma","phisigma", phisigma,"pdf","width=4.5in"); 

  /// phi slope ///
  TH2I *phidelta_z_h = new TH2I ("phidelta_z_h","#phi vs. #Delta dE/dx  / range", 20,0,40, 36,-TMath::Pi(), TMath::Pi()); 
  fit->Draw("atan2(sin(htphi.val-mcPhi),cos(htphi.val-mcPhi)) : delta_z.val / range.val  >> phidelta_z_h","","goff"); 
  phidelta_z_h->GetXaxis()->SetTitle("#Delta dE/dX / range  (adu/px^{2})"); 
  TCanvas * phidelta_z = phiVsX(phidelta_z_h); 
  slides.AddPlotSlide("fit $\\phi$ vs. $\\Delta dE/dX$","phidelta_z", phidelta_z,"pdf","width=4.5in"); 
 
 
  /// phi err ///
  TH2I *phiphierr_h = new TH2I ("phiphierr_h","#phi vs. #phi err", 20,0,TMath::Pi()/8, 36,-TMath::Pi(), TMath::Pi()); 
  fit->Draw("atan2(sin(htphi.val-mcPhi),cos(htphi.val-mcPhi)) : phi.err  >> phiphierr_h","","goff"); 
  phiphierr_h->GetXaxis()->SetTitle("#phi Error"); 
  TCanvas * phiphierr = phiVsX(phiphierr_h); 
  slides.AddPlotSlide("fit $\\phi$ vs. $\\phi$ error","phiphierr", phiphierr,"pdf","width=4.5in"); 

   /// phi sigma err ///
  TH2I *phisigmaerr_h = new TH2I ("phisigmaerr_h","#phi vs. rel #sigma err", 20,0,0.25, 36,-TMath::Pi(), TMath::Pi()); 
  fit->Draw("atan2(sin(htphi.val-mcPhi),cos(htphi.val-mcPhi)) : sigma.err/sigma.val  >> phisigmaerr_h","","goff"); 
  phisigmaerr_h->GetXaxis()->SetTitle("rel #sigma Error"); 
  TCanvas * phisigmaerr = phiVsX(phisigmaerr_h); 
  slides.AddPlotSlide("fit $\\phi$ vs. $\\sigma$ error","phisigmaerr", phisigmaerr,"pdf","width=4.5in"); 
 






}
