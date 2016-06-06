{
  gSystem->Load("/net/fritz/data02/cozzyd/slide_maker/libSlideMaker.so"); 
  gSystem->Load("/net/zwicky/dmtpc/cozzyd/projects/DarkMatter/MaxCam/libMaxCam.so"); 
  gStyle->SetCanvasPreferGL(true); 

  DmtpcRootTools::setColorJet(); 

  char * out_prefix = "4shDirectionality/"; 

  BeamerSlideMakerSettings slide_settings; 
  slide_settings->author="Cosmin Deaconu" ; 
  slide_settings->title="4sh Directionality With Fitting"; 
  slide_settings->short_title="4sh Directionality"; 
  slide_settings->style="Boadilla"; 
  slide_settings->color_style="beaver"; 
  slide_settings->navigation_symbols=false; 
  slide_settings->institute="MIT"; 
  slide_settings->title_page = false; 


  BeamerSlideMaker slides(out_prefix,"4shDirectionality", &slide_settings); 


  char * cam_names[4] =  {"A80333","110121","100534","A80334"}; 
  double xcenters[4]={23.250000,14.750000,24.750000,14.250000};
  double ycenters[4]={1005.250000,996.250000,1017.250000,1017.250000};
  double right_direction_each[4] = {-TMath::Pi()/2, 0, TMath::Pi(), TMath::Pi()/2}; 
  double right_direction_all = -1.8; 


  TChain * fit = new TChain("fit"); 
  fit->Add("../IntersectionFit.root"); 
  fit->Add("../IntersectionFitHigh.root"); 



  /***************************************************************/
  /***---------------- chisq  -------  --------------------------***/
  /***************************************************************/

  
  TCanvas * chisqc = new TCanvas("chisqc","chisqc",800,600); 
  chisqc->Divide(2,2); 

  for (int icam = 0; icam <4; icam++)
  {
    chisqc->cd(icam+1);
    TH1I *chisqh = new TH1I ("chisqh","chisqh",50,0, 15); 
    fit->Draw("fitChisq/fitNdof >> chisqh",TString::Format("cam==%d",icam)); 
    chisqh->SetLineColor(2); 
    chisqh->SetTitle(cam_names[icam]); 
    chisqh->GetXaxis()->SetTitle("chi^2/ndof"); 
    chisqh->DrawCopy("same" ); 
  }


  /***************************************************************/
  /***---------------- phi  -------  --------------------------***/
  /***************************************************************/

  
  TCanvas * rawPhi = new TCanvas("rawphi","rawphi",800,600); 
  rawPhi->Divide(2,2); 

  for (int icam = 0; icam <4; icam++)
  {
    rawPhi->cd(icam+1);
    TH1I *rawphi_h = new TH1I ("rawphi_h","rawphi_h",20,-TMath::Pi(), TMath::Pi()); 
    fit->Draw("fitPhi >> rawphi_h",TString::Format("cam==%d",icam)); 
    rawphi_h->SetLineColor(4); 
    rawphi_h->SetTitle(cam_names[icam]); 
    rawphi_h->DrawCopy(); 
    rawphi_h->SetLineColor(2); 
    fit->Draw("phi2 >> rawphi_h",TString::Format("cam==%d",icam),"same"); 
    rawphi_h->DrawCopy("same" ); 
  }


  /***************************************************************/
  /***---------------- combined phi  --------------------------***/
  /***************************************************************/

  TCanvas * combined = new TCanvas("combined","combined",800,600); 

  TH1I *combined_h = new TH1I ("combined_h","combined_h",20,-TMath::Pi(), TMath::Pi()); 
  fit.Draw("atan2(sin(fitPhi-offset), cos(fitPhi-offset)) >> combined_h"); 
  combined_h->SetLineColor(4); 
  combined_h->DrawCopy(); 
  fit.Draw("atan2(sin(phi2-offset), cos(phi2-offset)) >> combined_h","","same"); 
  combined_h->SetLineColor(2); 
  combined_h->DrawCopy("same"); 





  /***************************************************************/
  /***----------------split combined phi contrib---------------***/
  /***************************************************************/

  TCanvas * split = new TCanvas("split","split",800,600); 
  fit.Draw("atan2(sin(fitPhi-offset),cos(fitPhi-offset)) >> combined_h"); 
  combined_h->SetLineColor(1); 
  combined_h->SetMinimum(0); 
  combined_h->DrawCopy(); 
  fit.Draw("atan2(sin(fitPhi-offset),cos(fitPhi-offset)) >> combined_h","cam==0","same"); 
  combined_h->SetLineColor(2); 
  combined_h->DrawCopy("same"); 
  fit.Draw("atan2(sin(fitPhi-offset),cos(fitPhi-offset)) >> combined_h","cam==1","same"); 
  combined_h->SetLineColor(3); 
  combined_h->DrawCopy("same"); 
  fit.Draw("atan2(sin(fitPhi-offset),cos(fitPhi-offset)) >> combined_h","cam==2","same"); 
  combined_h->SetLineColor(4); 
  combined_h->DrawCopy("same"); 
  fit.Draw("atan2(sin(fitPhi-offset),cos(fitPhi-offset)) >> combined_h","cam==3","same"); 
  combined_h->SetLineColor(5); 
  combined_h->DrawCopy("same"); 


 



  /***************************************************************/
  /***----------------phi vs. Prob each camera ----------------***/
  /***************************************************************/

  TCanvas * phiVsProb = new TCanvas("phiprobc","phiprobc",800,600); 
  phiVsProb->Divide(2,2); 
  for (int icam = 0; icam < 4; icam++)
  {
    TH2I * phiprob = new TH2I ("phivprob","phivprob",16,0.4,1.2, 20,-TMath::Pi(), TMath::Pi()); 
    phiVsProb->cd(icam+1); 
    fit->Draw("fitPhi:fitProbability >> phivprob",TString::Format("cam==%d && fitChisq/fitNdof< 10",icam)); 
    phiprob->SetTitle(cam_names[icam]); 
    phiprob->GetXaxis()->SetTitle("implied probability"); 
    phiprob->GetYaxis()->SetTitle("#phi"); 
    phiprob->DrawCopy("colz"); 
  }


  /***************************************************************/
  /***----------------r/sigma vs. Prob each camera ----------------***/
  /***************************************************************/

  TCanvas * rsVsProb = new TCanvas("rsprobc","rsprobc",800,600); 
  rsVsProb->Divide(2,2); 
  for (int icam = 0; icam < 4; icam++)
  {
    TH2I * rsprob = new TH2I ("rsprob","rsprob",16,0.4,1.2, 20,0,10); 
    rsVsProb->cd(icam+1); 
    fit->Draw("fitRange/sigma:fitProbability >> rsprob",TString::Format("cam==%d && fitChisq/fitNdof< 10",icam)); 
    rsprob->SetTitle(cam_names[icam]); 
    rsprob->GetXaxis()->SetTitle("implied probability"); 
    rsprob->GetYaxis()->SetTitle("fit range / sigma"); 
    rsprob->DrawCopy("colz"); 
  }



  /***************************************************************/
  /***----------------- visE vs. meshE each camera ----------------***/
  /***************************************************************/

  TCanvas * visEmeshE = new TCanvas("vemec","vemec",800,600); 
  visEmeshE->Divide(2,2); 
  for (int icam = 0; icam < 4; icam++)
  {
    TH2I * veme = new TH2I ("veme","veme",30,0,300, 30,0, 300); 
    visEmeshE->cd(icam+1); 
    fit->Draw("recoE/gain:chargeE >> veme",TString::Format("cam==%d && chargeE > 0",icam)); 
    veme->SetTitle(cam_names[icam]); 
    veme->GetXaxis()->SetTitle("Emesh (keV)"); 
    veme->GetYaxis()->SetTitle("Evis (keV)"); 
    veme->DrawCopy("colz"); 
    veme->ProfileX()->DrawCopy("same"); 
  }

  /***************************************************************/
  /***----------------- fitRange vs. meshE each camera ----------------***/
  /***************************************************************/

  TCanvas * rvec = new TCanvas("rvec","rvec",800,600); 
  rvec->Divide(2,2); 
  for (int icam = 0; icam < 4; icam++)
  {
    TH2I * rve = new TH2I ("rve","rve",30,0, 300,30,0,10); 
    rvec->cd(icam+1); 
    fit->Draw("fitRange*0.16:chargeE >> rve",TString::Format("cam==%d && chargeE > 0",icam)); 
    rve->SetTitle(cam_names[icam]); 
    rve->GetXaxis()->SetTitle("Emesh (keV)"); 
    veme->GetYaxis()->SetTitle("Fit Range (mm)"); 
    rve->DrawCopy("colz"); 
    rve->ProfileX()->DrawCopy("same"); 
  }




  /***************************************************************/
  /***----------------- E vs. Prob each camera ----------------***/
  /***************************************************************/

  TCanvas * EVsProb = new TCanvas("eprobc","eprobc",800,600); 
  EVsProb->Divide(2,2); 
  for (int icam = 0; icam < 4; icam++)
  {
    TH2I * eprob = new TH2I ("eprob","eprob",16,0.4,1.2, 30,0, 300); 
    EVsProb->cd(icam+1); 
    fit->Draw("recoE/gain:fitProbability >> eprob",TString::Format("cam==%d",icam)); 
    eprob->SetTitle(cam_names[icam]); 
    eprob->GetXaxis()->SetTitle("implied probability"); 
    eprob->GetYaxis()->SetTitle("Evis (keV)"); 
    eprob->DrawCopy("colz"); 
  }

  /***************************************************************/
  /***-----------------mesh E vs. Prob each camera ----------------***/
  /***************************************************************/

  TCanvas * mEVsProb = new TCanvas("meprobc","meprobc",800,600); 
  mEVsProb->Divide(2,2); 
  for (int icam = 0; icam < 4; icam++)
  {
    TH2I * meprob = new TH2I ("meprob","meprob",16,0.4,1.2, 30,0, 300); 
    mEVsProb->cd(icam+1); 
    fit->Draw("chargeE:fitProbability >> meprob",TString::Format("cam==%d && chargeE > 0",icam)); 
    meprob->SetTitle(cam_names[icam]); 
    meprob->GetXaxis()->SetTitle("implied probability"); 
    meprob->GetYaxis()->SetTitle("Emesh (keV)"); 
    meprob->DrawCopy("colz"); 
  }


  /***************************************************************/
  /***--------------sigma dist E slices ----------------***/
  /***************************************************************/

  TCanvas * sigeslicec = new TCanvas("sigeslicec","sigeslicec",800,600); 
  TH1I * sigmaeslice = new TH1I ("sigmaeslice","sigmaslice",20,0, 1.6); 

  fit->Draw("sigma * 0.16 >> sigmaeslice","chargeE>0"); 

  sigmaeslice->GetXaxis()->SetTitle("fit sigma (mm)"); 
  sigmaeslice->SetLineColor(1); 
  sigmaeslice->DrawCopy(); 


  fit->Draw("sigma * 0.16 >> sigmaeslice", "chargeE>0 && chargeE <= 100","goff"); 
  sigmaeslice->SetLineColor(2); 
  sigmaeslice->DrawCopy("same"); 

  fit->Draw("sigma * 0.16 >> sigmaeslice", "chargeE>100 && chargeE <= 200","goff"); 
  sigmaeslice->SetLineColor(3); 
  sigmaeslice->DrawCopy("same"); 

  fit->Draw("sigma * 0.16 >> sigmaeslice", "chargeE>200","goff"); 
  sigmaeslice->SetLineColor(4); 
  sigmaeslice->DrawCopy("same"); 



  /***************************************************************/
  /***--------------sigma dist each camera ----------------***/
  /***************************************************************/

  TCanvas * sigc = new TCanvas("sigc","sigc",800,600); 
  sigc->Divide(2,2); 
  for (int icam = 0; icam < 4; icam++)
  {
    sigc->cd(icam+1); 
    TH1I * sigma = new TH1I ("sigma","sigma",20,0, 10); 
    fit->Draw("sigma >> sigma",TString::Format("cam==%d",icam)); 
    sigma->SetTitle(cam_names[icam]); 
    sigma->GetXaxis()->SetTitle("fit sigma (pixels)"); 
    sigma->DrawCopy(); 
  }

  /***************************************************************/
  /***--------------sigma vs. Evis dist each camera ----------------***/
  /***************************************************************/

  TCanvas * sigec = new TCanvas("sigec","sigec",800,600); 
  sigec->Divide(2,2); 
  for (int icam = 0; icam < 4; icam++)
  {
    TH2I * sigmae = new TH2I ("sigmae","sigmae",20,0, 1.6, 30,0,300); 
    sigec->cd(icam+1); 
    fit->Draw("recoE/gain : sigma*0.16 >> sigmae",TString::Format("cam==%d",icam)); 
    sigmae->SetTitle(cam_names[icam]); 
    sigmae->GetXaxis()->SetTitle("fit sigma (mm)"); 
    sigmae->GetYaxis()->SetTitle("EVis (kev)"); 
    sigmae->DrawCopy("colz"); 
  }

  /***************************************************************/
  /***--------------sigma vs. Emesh dist each camera ----------------***/
  /***************************************************************/

  TCanvas * sigmec = new TCanvas("sigmec","sigmec",800,600); 
  sigmec->Divide(2,2); 
  for (int icam = 0; icam < 4; icam++)
  {
    TH2I * sigmame = new TH2I ("sigmame","sigmame",20,0, 1.6, 30,0,300); 
    sigmec->cd(icam+1); 
    fit->Draw("chargeE : sigma * 0.16 >> sigmame",TString::Format("cam==%d && chargeE>0",icam)); 
    sigmame->SetTitle(cam_names[icam]); 
    sigmame->GetXaxis()->SetTitle("fit sigma (mm)"); 
    sigmame->GetYaxis()->SetTitle("EMesh (kev)"); 
    sigmame->DrawCopy("colz"); 
  }

  /***************************************************************/
  /***--------------sigma vs. range dist each camera ----------------***/
  /***************************************************************/

  TCanvas * sigrc = new TCanvas("sigrc","sigrc",800,600); 
  sigrc->Divide(2,2); 
  for (int icam = 0; icam < 4; icam++)
  {
    TH2I * sigmar = new TH2I ("sigmar","sigmar",20,0, 1.6, 30,0,10); 
    sigrc->cd(icam+1); 
    fit->Draw("fitRange * 0.16 : sigma * 0.16 >> sigmar",TString::Format("cam==%d",icam)); 
    sigmar->SetTitle(cam_names[icam]); 
    sigmar->GetXaxis()->SetTitle("fit sigma (mm)"); 
    sigmar->GetYaxis()->SetTitle("Range (mm)"); 
    sigmar->DrawCopy("colz"); 
  }






  /*******************************************************/
  /***---phi spectra with probability cut each camera -***/
  /*******************************************************/

  TCanvas * cutPhi = new TCanvas("cutphi","cutphi",800,600); 
  cutPhi->Divide(2,2); 

  for (int icam = 0; icam <4; icam++)
  {
    cutPhi->cd(icam+1);
    TH1I *cutphi_h = new TH1I ("cutphi_h","cutphi_h",20,-TMath::Pi(), TMath::Pi()); 
    fit->Draw("fitPhi >> cutphi_h",TString::Format("cam==%d",icam)); 
    cutphi_h->SetLineColor(1); 
    cutphi_h->SetTitle(cam_names[icam]); 
    cutphi_h->GetXaxis()->SetTitle("fit #phi"); 
    cutphi_h->SetMinimum(0); 
    cutphi_h->DrawCopy(); 

    fit->Draw("fitPhi >> cutphi_h",TString::Format("cam==%d && fitProbability > 0.8 && fitChisq/fitNdof < 10",icam),"same"); 
    cutphi_h->SetLineColor(2); 
    cutphi_h->DrawCopy("same"); 

    fit->Draw("fitPhi >> cutphi_h",TString::Format("cam==%d && fitProbability > 0.6 && fitProbability <= 0.8 && fitChisq/fitNdof < 10",icam),"same"); 
    cutphi_h->SetLineColor(3); 
    cutphi_h->DrawCopy("same"); 

    fit->Draw("fitPhi >> cutphi_h",TString::Format("cam==%d && fitProbability <= 0.6 && fitChisq/fitNdof < 10",icam),"same"); 
    cutphi_h->SetLineColor(4); 
    cutphi_h->DrawCopy("same"); 
  }




  /*******************************************************/
  /***---probability energy spectra each camera -------***/
  /*******************************************************/
  TCanvas * cutE = new TCanvas("cutE","cutE",800,600); 
  cutE->Divide(2,2); 

  for (int icam = 0; icam <4; icam++)
  {
    cutE->cd(icam+1);
    TH1I *cute_h = new TH1I ("cute_h","cute_h",30,0, 300); 
    fit->Draw("recoE/gain >> cute_h",TString::Format("cam==%d",icam)); 
    cute_h->SetLineColor(1); 
    cute_h->SetTitle(cam_names[icam]); 
    cute_h->GetXaxis()->SetTitle("Evis (keV)"); 
    cute_h->SetMinimum(0); 
    cute_h->DrawCopy(); 

    fit->Draw("recoE/gain >> cute_h",TString::Format("cam==%d && fitProbability > 0.8 && fitChisq/fitNdof < 10",icam),"same"); 
    cute_h->SetLineColor(2); 
    cute_h->DrawCopy("same"); 

    fit->Draw("recoE/gain >> cute_h",TString::Format("cam==%d && fitProbability > 0.6 && fitProbability <= 0.8 && fitChisq/fitNdof < 10",icam),"same"); 
    cute_h->SetLineColor(3); 
    cute_h->DrawCopy("same"); 

    fit->Draw("recoE/gain >> cute_h",TString::Format("cam==%d && fitProbability <= 0.6 && fitChisq/fitNdof < 10",icam),"same"); 
    cute_h->SetLineColor(4); 
    cute_h->DrawCopy("same"); 
  }

  /*******************************************************/
  /***---probability mesh energy spectra each camera -------***/
  /*******************************************************/
  TCanvas * cutmE = new TCanvas("cutmE","cutmE",800,600); 
  cutmE->Divide(2,2); 

  for (int icam = 0; icam <4; icam++)
  {
    cutmE->cd(icam+1);
    TH1I *cutme_h = new TH1I ("cutme_h","cutme_h",30,0, 300); 
    fit->Draw("chargeE >> cutme_h",TString::Format("cam==%d && chargeE>0",icam)); 
    cutme_h->SetLineColor(1); 
    cutme_h->SetTitle(cam_names[icam]); 
    cutme_h->GetXaxis()->SetTitle("Emesh (keV)"); 
    cutme_h->SetMinimum(0); 
    cutme_h->DrawCopy(); 

    fit->Draw("chargeE >> cutme_h",TString::Format("cam==%d && fitProbability > 0.8 && fitChisq/fitNdof < 10 && chargeE > 0",icam),"same"); 
    cutme_h->SetLineColor(2); 
    cutme_h->DrawCopy("same"); 

    fit->Draw("chargeE >> cutme_h",TString::Format("cam==%d && fitProbability > 0.6 && fitProbability <= 0.8 && fitChisq/fitNdof < 10 && chargeE>0",icam),"same"); 
    cutme_h->SetLineColor(3); 
    cutme_h->DrawCopy("same"); 

    fit->Draw("chargeE >> cutme_h",TString::Format("cam==%d && fitProbability <= 0.6 && fitChisq/fitNdof < 10 && chargeE>0",icam),"same"); 
    cutme_h->SetLineColor(4); 
    cutme_h->DrawCopy("same"); 
  }



  /*******************************************************/
  /***---phi spectra with probability cut combined -***/
  /*******************************************************/

  TCanvas * cutPhiAll = new TCanvas("cutphiall","cutphiall",800,600); 

  TH1I *cutphiall_h = new TH1I ("cutphiall_h","cutphiall_h",20,-TMath::Pi(), TMath::Pi()); 
  fit->Draw("atan2(sin(fitPhi-offset),cos(fitPhi-offset)) >> cutphiall_h"); 
  cutphiall_h->SetLineColor(1); 
  cutphiall_h->SetTitle("combined"); 
  cutphiall_h->GetXaxis()->SetTitle("fit #phi"); 
  cutphiall_h->SetMinimum(0); 
  cutphiall_h->DrawCopy(); 

  fit->Draw("atan2(sin(fitPhi-offset),cos(fitPhi-offset)) >> cutphiall_h","fitProbability > 0.8 && fitChisq/fitNdof < 10","goff"); 
  cutphiall_h->SetLineColor(2); 
  cutphiall_h->DrawCopy("same"); 

  fit->Draw("atan2(sin(fitPhi-offset),cos(fitPhi-offset)) >> cutphiall_h","fitProbability > 0.6 && fitProbability <= 0.8 && fitChisq/fitNdof < 10","goff"); 
  cutphiall_h->SetLineColor(3); 
  cutphiall_h->DrawCopy("same"); 

  fit->Draw("atan2(sin(fitPhi-offset),cos(fitPhi-offset)) >> cutphiall_h","fitProbability <= 0.6 && fitChisq/fitNdof < 10","goff"); 
  cutphiall_h->SetLineColor(4); 
  cutphiall_h->DrawCopy("same"); 


  /***************************************************************/
  /***----------------- Head Tail each camera mesh p slices---------***/
  /***************************************************************/

  TCanvas * htpmc = new TCanvas("htpmc","htpmc",800,600); 
  htpmc->Divide(2,2); 
  for (int icam = 0; icam < 4; icam++)
  {
    htpmc->cd(icam+1); 
    TH2* htpm = new TH2F("htpm","htpm", 6,0,300,2,-0.5,1.5); 
    fit->Draw(TString::Format("cos(fitPhi-%f) > 0 : chargeE >> htpm", right_direction_each[icam]), TString::Format("cam==%d && chargeE>0",icam),"goff"); 
    TProfile * prof = htpm->ProfileX(); 
    prof->SetMaximum(1); 
    prof->SetMinimum(0); 

    prof->SetTitle(cam_names[icam]); 
    prof->GetXaxis()->SetTitle("Emesh keV"); 
    prof->GetYaxis()->SetTitle("fraction correct"); 
    prof->SetLineColor(1); 
    prof->SetMarkerColor(1); 
    prof->DrawCopy(); 

    fit->Draw(TString::Format("cos(fitPhi-%f) > 0 : chargeE >> htpm", right_direction_each[icam]), TString::Format("cam==%d && fitProbability > 0.8 && fitChisq/fitNdof < 10 && chargeE>0",icam),"goff"); 
    TProfile * prof = htpm->ProfileX(); 
    prof->SetLineColor(2); 
    prof->SetMarkerColor(2); 
    prof->DrawCopy("same"); 

    fit->Draw(TString::Format("cos(fitPhi-%f) > 0 : chargeE >> htpm", right_direction_each[icam]), TString::Format("cam==%d && fitProbability > 0.6 && fitProbability <= 0.8 && fitChisq/fitNdof < 10 && chargeE>0",icam),"goff"); 
    TProfile * prof = htpm->ProfileX(); 
    prof->SetLineColor(3); 
    prof->SetMarkerColor(3); 
    prof->DrawCopy("same"); 

    fit->Draw(TString::Format("cos(fitPhi-%f) > 0 : chargeE >> htpm", right_direction_each[icam]), TString::Format("cam==%d && fitProbability <= 0.6 && fitChisq/fitNdof < 10 && chargeE>0",icam),"goff"); 
    TProfile * prof = htpm->ProfileX(); 
    prof->SetLineColor(4); 
    prof->SetMarkerColor(4); 
    prof->DrawCopy("same"); 

    TLine l(0,0.5,300,0.5); 
    l.SetLineStyle(3); 
    l.SetLineColor(6); 
    l.Draw(); 
    
  }



  /***************************************************************/
  /***----------------- Head Tail each camera p slices---------***/
  /***************************************************************/

  TCanvas * htpc = new TCanvas("htpc","htpc",800,600); 
  htpc->Divide(2,2); 
  for (int icam = 0; icam < 4; icam++)
  {
    htpc->cd(icam+1); 
    TH2* htp = new TH2F("htp","htp", 6,0,300,2,-0.5,1.5); 
    fit->Draw(TString::Format("cos(fitPhi-%f) > 0 : recoE/gain >> htp", right_direction_each[icam]), TString::Format("cam==%d",icam),"goff"); 
    TProfile * prof = htp->ProfileX(); 
    prof->SetMaximum(1); 
    prof->SetMinimum(0); 

    prof->SetTitle(cam_names[icam]); 
    prof->GetXaxis()->SetTitle("Evis keV"); 
    prof->GetYaxis()->SetTitle("fraction correct"); 
    prof->SetLineColor(1); 
    prof->SetMarkerColor(1); 
    prof->DrawCopy(); 

    fit->Draw(TString::Format("cos(fitPhi-%f) > 0 : recoE/gain >> htp", right_direction_each[icam]), TString::Format("cam==%d && fitProbability > 0.8 && fitChisq/fitNdof < 10",icam),"goff"); 
    TProfile * prof = htp->ProfileX(); 
    prof->SetLineColor(2); 
    prof->SetMarkerColor(2); 
    prof->DrawCopy("same"); 

    fit->Draw(TString::Format("cos(fitPhi-%f) > 0 : recoE/gain >> htp", right_direction_each[icam]), TString::Format("cam==%d && fitProbability > 0.6 && fitProbability <= 0.8 && fitChisq/fitNdof < 10",icam),"goff"); 
    TProfile * prof = htp->ProfileX(); 
    prof->SetLineColor(3); 
    prof->SetMarkerColor(3); 
    prof->DrawCopy("same"); 

    fit->Draw(TString::Format("cos(fitPhi-%f) > 0 : recoE/gain >> htp", right_direction_each[icam]), TString::Format("cam==%d && fitProbability <= 0.6 && fitChisq/fitNdof < 10",icam),"goff"); 
    TProfile * prof = htp->ProfileX(); 
    prof->SetLineColor(4); 
    prof->SetMarkerColor(4); 
    prof->DrawCopy("same"); 

    TLine l(0,0.5,300,0.5); 
    l.SetLineStyle(3); 
    l.SetLineColor(6); 
    l.Draw(); 
  }


  /***************************************************************/
  /***----------------- Head Tail combined    p slices---------***/
  /***************************************************************/

  TCanvas * htpac = new TCanvas("htpac","htpac",800,600); 
  TH2* htap = new TH2F("htap","htap", 12,0,300,2,-0.5,1.5); 
  fit->Draw(TString::Format("cos(fitPhi-offset-%f) > 0 : recoE/gain >> htap", right_direction_all),"","goff"); 
  TProfile * prof = htap->ProfileX(); 
  prof->SetMaximum(1); 
  prof->SetMinimum(0); 
  prof->SetTitle("combined"); 
  prof->GetXaxis()->SetTitle("Evis keV"); 
  prof->GetYaxis()->SetTitle("fraction correct"); 
  prof->SetLineColor(1); 
  prof->SetMarkerColor(1); 
  prof->DrawCopy(); 

  fit->Draw(TString::Format("cos(fitPhi-offset-%f) > 0 : recoE/gain >> htap", right_direction_all),"fitProbability > 0.8 && fitChisq/fitNdof < 10","goff"); 
  TProfile * prof = htap->ProfileX(); 
  prof->SetLineColor(2); 
  prof->SetMarkerColor(2); 
  prof->DrawCopy("same"); 

  fit->Draw(TString::Format("cos(fitPhi-offset-%f) > 0 : recoE/gain >> htap", right_direction_all),"fitProbability > 0.6 && fitProbability <= 0.8 && fitChisq/fitNdof < 10","goff"); 
  TProfile * prof = htap->ProfileX(); 
  prof->SetLineColor(3); 
  prof->SetMarkerColor(3); 
  prof->DrawCopy("same"); 

  fit->Draw(TString::Format("cos(fitPhi-offset-%f) > 0 : recoE/gain >> htap", right_direction_all),"fitProbability <= 0.6 && fitChisq/fitNdof < 10","goff"); 
  TProfile * prof = htap->ProfileX(); 
  prof->SetLineColor(4); 
  prof->SetMarkerColor(4); 
  prof->DrawCopy("same"); 

  TLine l(0,0.5,300,0.5); 
  l.SetLineStyle(3); 
  l.SetLineColor(6); 
  l.Draw(); 
  /***************************************************************/
  /***----------------- Head Tail combined  mesh E  p slices---------***/
  /***************************************************************/

  TCanvas * htpmac = new TCanvas("htpmac","htpmac",800,600); 
  TH2* htmap = new TH2F("htmap","htmap", 12,0,300,2,-0.5,1.5); 
  fit->Draw(TString::Format("cos(fitPhi-offset-%f) > 0 : chargeE >> htmap", right_direction_all),"chargeE>0","goff"); 
  TProfile * prof = htmap->ProfileX(); 
  prof->SetMaximum(1); 
  prof->SetMinimum(0); 
  prof->SetTitle("combined"); 
  prof->GetXaxis()->SetTitle("Emesh keV"); 
  prof->GetYaxis()->SetTitle("fraction correct"); 
  prof->SetLineColor(1); 
  prof->SetMarkerColor(1); 
  prof->DrawCopy(); 

  fit->Draw(TString::Format("cos(fitPhi-offset-%f) > 0 : chargeE >> htmap", right_direction_all),"fitProbability > 0.8 && fitChisq/fitNdof < 10 && chargeE > 0","goff"); 
  TProfile * prof = htmap->ProfileX(); 
  prof->SetLineColor(2); 
  prof->SetMarkerColor(2); 
  prof->DrawCopy("same"); 

  fit->Draw(TString::Format("cos(fitPhi-offset-%f) > 0 : chargeE >> htmap", right_direction_all),"fitProbability > 0.6 && fitProbability <= 0.8 && fitChisq/fitNdof < 10 && chargeE > 0","goff"); 
  TProfile * prof = htmap->ProfileX(); 
  prof->SetLineColor(3); 
  prof->SetMarkerColor(3); 
  prof->DrawCopy("same"); 

  fit->Draw(TString::Format("cos(fitPhi-offset-%f) > 0 : chargeE >> htmap", right_direction_all),"fitProbability <= 0.6 && fitChisq/fitNdof < 10 && chargeE >0","goff"); 
  TProfile * prof = htmap->ProfileX(); 
  prof->SetLineColor(4); 
  prof->SetMarkerColor(4); 
  prof->DrawCopy("same"); 

  l.Draw(); 



  /***************************************************************/
  /***----------------- Head Tail v range each camera p slices---------***/
  /***************************************************************/

  TCanvas * htrpc = new TCanvas("htrpc","htrpc",800,600); 
  htrpc->Divide(2,2); 
  for (int icam = 0; icam < 4; icam++)
  {
    htrpc->cd(icam+1); 
    TH2* htpr = new TH2F("htpr","htpr", 10,0,80,2,-0.5,1.5); 
    fit->Draw(TString::Format("cos(fitPhi-%f) > 0 : fitRange >> htp", right_direction_each[icam]), TString::Format("cam==%d",icam),"goff"); 
    TProfile * prof = htpr->ProfileX(); 
    prof->SetMaximum(1); 
    prof->SetMinimum(0); 

    prof->SetTitle(cam_names[icam]); 
    prof->GetXaxis()->SetTitle("fit range (pixels)"); 
    prof->GetYaxis()->SetTitle("fraction correct"); 
    prof->SetLineColor(1); 
    prof->SetMarkerColor(1); 
    prof->DrawCopy(); 

    fit->Draw(TString::Format("cos(fitPhi-%f) > 0 : fitRange >> htpr", right_direction_each[icam]), TString::Format("cam==%d && fitProbability > 0.8 && fitChisq/fitNdof < 10",icam),"goff"); 
    TProfile * prof = htpr->ProfileX(); 
    prof->SetLineColor(2); 
    prof->SetMarkerColor(2); 
    prof->DrawCopy("same"); 

    fit->Draw(TString::Format("cos(fitPhi-%f) > 0 : fitRange >> htpr", right_direction_each[icam]), TString::Format("cam==%d && fitProbability > 0.6 && fitProbability <= 0.8 && fitChisq/fitNdof < 10",icam),"goff"); 
    TProfile * prof = htpr->ProfileX(); 
    prof->SetLineColor(3); 
    prof->SetMarkerColor(3); 
    prof->DrawCopy("same"); 

    fit->Draw(TString::Format("cos(fitPhi-%f) > 0 : fitRange >> htpr", right_direction_each[icam]), TString::Format("cam==%d && fitProbability <= 0.6 && fitChisq/fitNdof < 10",icam),"goff"); 
    TProfile * prof = htpr->ProfileX(); 
    prof->SetLineColor(4); 
    prof->SetMarkerColor(4); 
    prof->DrawCopy("same"); 
  }


  /***************************************************************/
  /***----------------- Head Tail range combined    p slices---------***/
  /***************************************************************/

  TCanvas * htprac = new TCanvas("htprac","htprac",800,600); 
  TH2* htrap = new TH2F("htrap","htrap", 10,0,80,2,-0.5,1.5); 
  fit->Draw(TString::Format("cos(fitPhi-offset-%f) > 0 : fitRange >> htrap", right_direction_all),"","goff"); 
  TProfile * prof = htrap->ProfileX(); 
  prof->SetMaximum(1); 
  prof->SetMinimum(0); 
  prof->SetTitle("combined"); 
  prof->GetXaxis()->SetTitle("fit range in pixels"); 
  prof->GetYaxis()->SetTitle("fraction correct"); 
  prof->SetLineColor(1); 
  prof->SetMarkerColor(1); 
  prof->DrawCopy(); 

  fit->Draw(TString::Format("cos(fitPhi-offset-%f) > 0 : fitRange >> htrap", right_direction_all),"fitProbability > 0.8 && fitChisq/fitNdof < 10","goff"); 
  TProfile * prof = htrap->ProfileX(); 
  prof->SetLineColor(2); 
  prof->SetMarkerColor(2); 
  prof->DrawCopy("same"); 

  fit->Draw(TString::Format("cos(fitPhi-offset-%f) > 0 : fitRange >> htrap", right_direction_all),"fitProbability > 0.6 && fitProbability <= 0.8 && fitChisq/fitNdof < 10","goff"); 
  TProfile * prof = htrap->ProfileX(); 
  prof->SetLineColor(3); 
  prof->SetMarkerColor(3); 
  prof->DrawCopy("same"); 

  fit->Draw(TString::Format("cos(fitPhi-offset-%f) > 0 : fitRange >> htrap", right_direction_all),"fitProbability <= 0.6 && fitChisq/fitNdof < 10","goff"); 
  TProfile * prof = htrap->ProfileX(); 
  prof->SetLineColor(4); 
  prof->SetMarkerColor(4); 
  prof->DrawCopy("same"); 



  /***************************************************************/
  /***----------------- Head Tail v sigma each camera p slices---------***/
  /***************************************************************/

  TCanvas * htspc = new TCanvas("htspc","htspc",800,600); 
  htspc->Divide(2,2); 
  for (int icam = 0; icam < 4; icam++)
  {
    htspc->cd(icam+1); 
    TH2* htps = new TH2F("htps","htps", 12,2,8,2,-0.5,1.5); 
    fit->Draw(TString::Format("cos(fitPhi-%f) > 0 : sigma >> htps", right_direction_each[icam]), TString::Format("cam==%d",icam),"goff"); 
    TProfile * prof = htps->ProfileX(); 
    prof->SetMaximum(1); 
    prof->SetMinimum(0); 

    prof->SetTitle(cam_names[icam]); 
    prof->GetXaxis()->SetTitle("fit sigma"); 
    prof->GetYaxis()->SetTitle("fraction correct"); 
    prof->SetLineColor(1); 
    prof->SetMarkerColor(1); 
    prof->DrawCopy(); 

    fit->Draw(TString::Format("cos(fitPhi-%f) > 0 : sigma >> htps", right_direction_each[icam]), TString::Format("cam==%d && fitProbability > 0.8 && fitChisq/fitNdof < 10",icam),"goff"); 
    TProfile * prof = htps->ProfileX(); 
    prof->SetLineColor(2); 
    prof->SetMarkerColor(2); 
    prof->DrawCopy("same"); 

    fit->Draw(TString::Format("cos(fitPhi-%f) > 0 : sigma >> htps", right_direction_each[icam]), TString::Format("cam==%d && fitProbability > 0.6 && fitProbability <= 0.8 && fitChisq/fitNdof < 10",icam),"goff"); 
    TProfile * prof = htps->ProfileX(); 
    prof->SetLineColor(3); 
    prof->SetMarkerColor(3); 
    prof->DrawCopy("same"); 

    fit->Draw(TString::Format("cos(fitPhi-%f) > 0 : sigma >> htps", right_direction_each[icam]), TString::Format("cam==%d && fitProbability <= 0.6 && fitChisq/fitNdof < 10 ",icam),"goff"); 
    TProfile * prof = htps->ProfileX(); 
    prof->SetLineColor(4); 
    prof->SetMarkerColor(4); 
    prof->DrawCopy("same"); 
  }


  /***************************************************************/
  /***----------------- Head Tail sigma combined p slices---------***/
  /***************************************************************/

  TCanvas * htpsac = new TCanvas("htpsac","htpsac",800,600); 
  TH2* htsap = new TH2F("htsap","htsap", 12,2,8,2,-0.5,1.5); 
  fit->Draw(TString::Format("cos(fitPhi-offset-%f) > 0 : sigma >> htsap", right_direction_all),"","goff"); 
  TProfile * prof = htsap->ProfileX(); 
  prof->SetMaximum(1); 
  prof->SetMinimum(0); 
  prof->SetTitle("combined"); 
  prof->GetXaxis()->SetTitle("fit sigma"); 
  prof->GetYaxis()->SetTitle("fraction correct"); 
  prof->SetLineColor(1); 
  prof->SetMarkerColor(1); 
  prof->DrawCopy(); 

  fit->Draw(TString::Format("cos(fitPhi-offset-%f) > 0 : sigma >> htsap", right_direction_all),"fitProbability > 0.8 && fitChisq/fitNdof < 10 ","goff"); 
  TProfile * prof = htsap->ProfileX(); 
  prof->SetLineColor(2); 
  prof->SetMarkerColor(2); 
  prof->DrawCopy("same"); 

  fit->Draw(TString::Format("cos(fitPhi-offset-%f) > 0 : sigma >> htsap", right_direction_all),"fitProbability > 0.6 && fitProbability <= 0.8 && fitChisq/fitNdof < 10 ","goff"); 
  TProfile * prof = htsap->ProfileX(); 
  prof->SetLineColor(3); 
  prof->SetMarkerColor(3); 
  prof->DrawCopy("same"); 

  fit->Draw(TString::Format("cos(fitPhi-offset-%f) > 0 : sigma >> htsap", right_direction_all),"fitProbability <= 0.6 && fitChisq/fitNdof < 10 ","goff"); 
  TProfile * prof = htsap->ProfileX(); 
  prof->SetLineColor(4); 
  prof->SetMarkerColor(4); 
  prof->DrawCopy("same"); 


  /***************************************************************/
  /***----------------- Head Tail v range/sigma each camera p slices---------***/
  /***************************************************************/

  TCanvas * htrspc = new TCanvas("htrspc","htrspc",800,600); 
  htrspc->Divide(2,2); 
  for (int icam = 0; icam < 4; icam++)
  {
    htrspc->cd(icam+1); 
    TH2* htprs = new TH2F("htprs","htprs", 12,0,12,2,-0.5,1.5); 
    fit->Draw(TString::Format("cos(fitPhi-%f) > 0 : fitRange/sigma >> htprs", right_direction_each[icam]), TString::Format("cam==%d",icam),"goff"); 
    TProfile * prof = htprs->ProfileX(); 
    prof->SetMaximum(1); 
    prof->SetMinimum(0); 

    prof->SetTitle(cam_names[icam]); 
    prof->GetXaxis()->SetTitle("fit range / sigma"); 
    prof->GetYaxis()->SetTitle("fraction correct"); 
    prof->SetLineColor(1); 
    prof->SetMarkerColor(1); 
    prof->DrawCopy(); 

    fit->Draw(TString::Format("cos(fitPhi-%f) > 0 : fitRange/sigma >> htprs", right_direction_each[icam]), TString::Format("cam==%d && fitProbability > 0.8 && fitChisq/fitNdof < 10",icam),"goff"); 
    TProfile * prof = htprs->ProfileX(); 
    prof->SetLineColor(2); 
    prof->SetMarkerColor(2); 
    prof->DrawCopy("same"); 

    fit->Draw(TString::Format("cos(fitPhi-%f) > 0 : fitRange/sigma >> htprs", right_direction_each[icam]), TString::Format("cam==%d && fitProbability > 0.6 && fitProbability <= 0.8 && fitChisq/fitNdof < 10",icam),"goff"); 
    TProfile * prof = htprs->ProfileX(); 
    prof->SetLineColor(3); 
    prof->SetMarkerColor(3); 
    prof->DrawCopy("same"); 

    fit->Draw(TString::Format("cos(fitPhi-%f) > 0 : fitRange/sigma >> htprs", right_direction_each[icam]), TString::Format("cam==%d && fitProbability <= 0.6 && fitChisq/fitNdof < 10 ",icam),"goff"); 
    TProfile * prof = htprs->ProfileX(); 
    prof->SetLineColor(4); 
    prof->SetMarkerColor(4); 
    prof->DrawCopy("same"); 
  }


  /***************************************************************/
  /***----------------- Head Tail range/sigma combined p slices---------***/
  /***************************************************************/

  TCanvas * htprsac = new TCanvas("htprsac","htprsac",800,600); 
  TH2* htrsap = new TH2F("htrsap","htrsap", 12,0,12,2,-0.5,1.5); 
  fit->Draw(TString::Format("cos(fitPhi-offset-%f) > 0 : fitRange/sigma >> htrsap", right_direction_all),"","goff"); 
  TProfile * prof = htrsap->ProfileX(); 
  prof->SetMaximum(1); 
  prof->SetMinimum(0); 
  prof->SetTitle("combined"); 
  prof->GetXaxis()->SetTitle("fit range/sigma"); 
  prof->GetYaxis()->SetTitle("fraction correct"); 
  prof->SetLineColor(1); 
  prof->SetMarkerColor(1); 
  prof->DrawCopy(); 

  fit->Draw(TString::Format("cos(fitPhi-offset-%f) > 0 : fitRange/sigma >> htrsap", right_direction_all),"fitProbability > 0.8 && fitChisq/fitNdof < 10 ","goff"); 
  TProfile * prof = htrsap->ProfileX(); 
  prof->SetLineColor(2); 
  prof->SetMarkerColor(2); 
  prof->DrawCopy("same"); 

  fit->Draw(TString::Format("cos(fitPhi-offset-%f) > 0 : fitRange/sigma >> htrsap", right_direction_all),"fitProbability > 0.6 && fitProbability <= 0.8 && fitChisq/fitNdof < 10 ","goff"); 
  TProfile * prof = htrsap->ProfileX(); 
  prof->SetLineColor(3); 
  prof->SetMarkerColor(3); 
  prof->DrawCopy("same"); 

  fit->Draw(TString::Format("cos(fitPhi-offset-%f) > 0 : fitRange/sigma >> htrsap", right_direction_all),"fitProbability <= 0.6 && fitChisq/fitNdof < 10 ","goff"); 
  TProfile * prof = htrsap->ProfileX(); 
  prof->SetLineColor(4); 
  prof->SetMarkerColor(4); 
  prof->DrawCopy("same"); 





 


  /***************************************************************/
  /***----------------- Head Tail each camera E slices---------***/
  /***************************************************************/

  TCanvas * htec = new TCanvas("htec","htec",800,600); 
  htec->Divide(2,2); 
  for (int icam = 0; icam < 4; icam++)
  {
    htec->cd(icam+1); 
    TH2* hte = new TH2F("hte","hte", 8,0.4,1.2,2,-0.5,1.5); 
    fit->Draw(TString::Format("cos(fitPhi-%f) > 0 : fitProbability >> hte", right_direction_each[icam]), TString::Format("cam==%d",icam),"goff"); 
    TProfile * prof = hte->ProfileX(); 
    prof->SetMaximum(1); 
    prof->SetMinimum(0); 
    prof->SetTitle(cam_names[icam]); 
    prof->GetXaxis()->SetTitle("implied probability"); 
    prof->GetYaxis()->SetTitle("fraction correct"); 
    prof->SetMarkerSize(1); 
    prof->SetLineColor(1); 
    prof->SetMarkerColor(1); 
    prof->DrawCopy(); 

    fit->Draw(TString::Format("cos(fitPhi-%f) > 0 : fitProbability >> hte", right_direction_each[icam]), TString::Format("cam==%d && recoE/gain > 150",icam),"goff"); 
    TProfile * prof = hte->ProfileX(); 
    prof->SetLineColor(2); 
    prof->SetMarkerColor(2); 
    prof->DrawCopy("same"); 

    fit->Draw(TString::Format("cos(fitPhi-%f) > 0 : fitProbability >> hte", right_direction_each[icam]), TString::Format("cam==%d && recoE/gain > 75 && recoE/gain <= 150",icam),"goff"); 
    TProfile * prof = hte->ProfileX(); 
    prof->SetLineColor(3); 
    prof->SetMarkerColor(3); 
    prof->DrawCopy("same"); 

    fit->Draw(TString::Format("cos(fitPhi-%f) > 0 : fitProbability >> hte", right_direction_each[icam]), TString::Format("cam==%d && recoE/gain <= 75",icam),"goff"); 
    TProfile * prof = hte->ProfileX(); 
    prof->SetLineColor(4); 
    prof->SetMarkerColor(4); 
    prof->DrawCopy("same"); 
    
    //redraw base
    fit->Draw(TString::Format("cos(fitPhi-%f) > 0 : fitProbability >> hte", right_direction_each[icam]), TString::Format("cam==%d",icam),"goff"); 
    TProfile * prof = hte->ProfileX(); 
    prof->SetLineColor(1); 
    prof->SetMarkerColor(1); 
    prof->SetMarkerSize(1.5); 
    prof->DrawCopy("same"); 
 
  }


  /***************************************************************/
  /***----------------- Head Tail combined    E slices---------***/
  /***************************************************************/

  TCanvas * hteac = new TCanvas("hteac","hteac",800,600); 
  TH2* htea = new TH2F("htea","htea", 8,0.4,1.2,2,-0.5,1.5); 
  fit->Draw(TString::Format("cos(fitPhi-offset-%f) > 0 : fitProbability >> htea", right_direction_all), "","goff"); 
  TProfile * prof = htea->ProfileX(); 
  prof->SetMaximum(1); 
  prof->SetMinimum(0); 
  prof->SetTitle("combined"); 
  prof->GetXaxis()->SetTitle("implied probability"); 
  prof->GetYaxis()->SetTitle("fraction correct"); 
  prof->SetMarkerSize(1); 
  prof->SetLineColor(1); 
  prof->SetMarkerColor(1); 
  prof->DrawCopy(); 

  fit->Draw(TString::Format("cos(fitPhi-offset-%f) > 0 : fitProbability >> htea", right_direction_all), "recoE/gain > 150","goff"); 
  TProfile * prof = htea->ProfileX(); 
  prof->SetLineColor(2); 
  prof->SetMarkerColor(2); 
  prof->DrawCopy("same"); 
 
  fit->Draw(TString::Format("cos(fitPhi-offset-%f) > 0 : fitProbability >> htea", right_direction_all), "recoE/gain > 75 && recoE/gain <= 150","goff"); 
  TProfile * prof = htea->ProfileX(); 
  prof->SetLineColor(3); 
  prof->SetMarkerColor(3); 
  prof->DrawCopy("same"); 
 
  fit->Draw(TString::Format("cos(fitPhi-offset-%f) > 0 : fitProbability >> htea", right_direction_all), "recoE/gain <= 75","goff"); 
  TProfile * prof = htea->ProfileX(); 
  prof->SetLineColor(4); 
  prof->SetMarkerColor(4); 
  prof->DrawCopy("same"); 
    
    //redraw base
  fit->Draw(TString::Format("cos(fitPhi-offset-%f) > 0 : fitProbability >> htea", right_direction_all), "","goff"); 
  TProfile * prof = htea->ProfileX(); 
  prof->SetLineColor(1); 
  prof->SetMarkerColor(1); 
  prof->SetMarkerSize(1.5); 
  prof->DrawCopy("same"); 






  /*******************************************************/
  /***-------------------      xy plots       ---------***/
  /*******************************************************/

  const int nxbins = 10; 
  const int nybins = 10; 

  TCanvas * densityc = new TCanvas("densityc","densityc",1000,1000); 
  TCanvas * mdensityc = new TCanvas("mdensityc","mdensityc",1000,1000); 
  TCanvas * ratioc = new TCanvas("ratioc","ratioc",1000,1000); 
  TH2 * density = new TH2I("density","",nxbins,-1000,1000,nybins,-1000,1000); 
  TH2 * mdensity = new TH2I("mdensity","",nxbins,-1000,1000,nybins,-1000,1000); 
  TH2 * ratio = new TH2F("ratio","",nxbins,-1000,1000,nybins,-1000,1000); 

  TCanvas * energyc = new TCanvas("energyc","energyc",1000,1000); 
  TCanvas * menergyc = new TCanvas("menergyc","menergyc",1000,1000); 
  TH2 * energy= new TH2F("energy","",nxbins,-1000,1000,nybins,-1000,1000); 
  TH2 * menergy= new TH2F("menergy","",nxbins,-1000,1000,nybins,-1000,1000); 
  TH2 * htxy= new TH2F("htxy","HTxy",nxbins,-1000,1000,nybins,-1000,1000); 

  TCanvas * rosec = new TCanvas("rosec","rosec",1000,1000); 
  TH2 * rose = new TH2F("rose","Red with p > 0.8",nxbins,-1000,1000,nybins,-1000,1000); 

  TH1* rose_phis_all[nxbins][nybins]; 
  TH1* rose_phis_cut[nxbins][nybins]; 


  int ibin; 
  for (unsigned i = 0; i < nxbins; i++)
  {
    for (unsigned j = 0; j < nybins; j++)
    {
      rose_phis_all[i][j] = new TH1S(TString::Format("phihist_%04d",ibin),"PhiDist!!",12,0,360); 
      rose_phis_cut[i][j] = new TH1S(TString::Format("phihist_%04d_cut",ibin++),"CutPhiDist",12,0,360); 
    }
  }
  

  

  double x, y, gain, recoE, fitPhi, offset, fitProbability, chisq, meshE; 
  int cam,ndof; 

  fit->SetBranchAddress("x",&x);
  fit->SetBranchAddress("y",&y);
  fit->SetBranchAddress("gain",&gain);
  fit->SetBranchAddress("recoE",&recoE);
  fit->SetBranchAddress("fitPhi",&fitPhi);
  fit->SetBranchAddress("offset",&offset);
  fit->SetBranchAddress("fitProbability",&fitProbability);
  fit->SetBranchAddress("fitChisq",&chisq);
  fit->SetBranchAddress("fitNdof",&ndof);
  fit->SetBranchAddress("cam",&cam);
  fit->SetBranchAddress("chargeE",&meshE);

  for (unsigned i = 0; i <fit->GetEntries(); i++)
  {
    fit->GetEntry(i); 
    double newX = (x-xcenters[cam]) * cos(offset) + (y-ycenters[cam]) * sin(-offset); 
    double newY = (x-xcenters[cam]) * sin(offset) + (y-ycenters[cam]) * cos(offset); 
    double phinew = fitPhi - offset; 

    double Evis = recoE / gain; 
    //int color = fitProbability > 0.8 ? 2 : fitProbability > 0.6 ? 3 : 4; 
    double xend = newX + Evis * TMath::Cos(phinew); 
    double yend = newY + Evis * TMath::Sin(phinew); 

    density->Fill(newX, newY); 
    energy->Fill(newX, newY, Evis); 
    if (meshE>0)
    {
      menergy->Fill(newX, newY, meshE); 
      mdensity->Fill(newX,newY);
      ratio->Fill(newX,newY,Evis/meshE); 
    }


    int binx = rose->GetXaxis()->FindBin(newX) -1; 
    int biny = rose->GetYaxis()->FindBin(newY) -1; 

    double angle_norm = DmtpcMath::normalizeAngle(phinew, TMath::Pi()) * 180 / TMath::Pi(); 
    rose_phis_all[binx][biny]->Fill(angle_norm);
    if (fitProbability > 0.8 && chisq / ndof < 10)
      rose_phis_cut[binx][biny]->Fill(angle_norm);

//    TArrow * arrow = new TArrow(newX,newY, xend,yend); 
//    arrow->SetLineColor(color); 
//    arrow->Draw(); 
  }

  densityc->cd(); 
  density->Draw("colz"); 

  mdensityc->cd(); 
  mdensity->Draw("colz"); 


  energy->Divide(density); 
  energyc->cd(); 
  energy->Draw("colz"); 

  menergy->Divide(mdensity); 
  menergyc->cd(); 
  menergy->Draw("colz"); 

  ratioc->cd(); 
  ratio->Divide(mdensity); 
  ratio->Draw("colz"); 


  rosec->cd(); 
  rose->Draw(); 
  double scale = 6; 

  for (unsigned i = 0; i < nxbins; i++)
  {
    for (unsigned j = 0; j < nybins; j++)
    {
       TH1S * phi = rose_phis_all[i][j];
       TH1S * phi_cut = rose_phis_cut[i][j]; 

       for (unsigned k = 1; k <= phi->GetNbinsX(); k++)
       {
          TCrown *crown = new TCrown(rose->GetXaxis()->GetBinCenter(i+1), 
                   rose->GetYaxis()->GetBinCenter(j+1), 
                   0, scale * phi->GetBinContent(k), 
                   phi->GetXaxis()->GetBinLowEdge(k), 
                   phi->GetXaxis()->GetBinLowEdge(k+1)); 

          crown->SetFillColor(4); 
          crown->Draw();  

          TCrown * cutcrown = new TCrown(rose->GetXaxis()->GetBinCenter(i+1), 
                   rose->GetYaxis()->GetBinCenter(j+1), 
                   0, scale * phi_cut->GetBinContent(k), 
                   phi->GetXaxis()->GetBinLowEdge(k) + 0.1 * phi->GetXaxis()->GetBinWidth(k), 
                   phi->GetXaxis()->GetBinLowEdge(k+1) - 0.1 * phi->GetXaxis()->GetBinWidth(k)); 

          cutcrown->SetFillColor(2); 
          cutcrown->Draw();  
       }

       htxy->SetBinContent(i+1,j+1, 1 - phi->Integral(1,phi->GetNbinsX()/2.) / phi->Integral()); 
    }
  }




  

  TCanvas *htxyc = new TCanvas("htxyc","htxyc",1000,1000); 
  htxy->Draw("colz"); 



  /////////////////////////////////////slides go here///////////////////////

const char * overview[] = {  
  "Ran fitting code on 60Torr AmBe playlists from intersection of Jeremy's + Shawn's rough nuclear recoil cuts (Shawn is CCD only, Jeremy is mostly charge)", 
  "5368 events total.", 
  "Analysis done on individual images, so some duplicate tracks in intersection",
  "\\alert{Caveat:} Gains only approximately known (Shawn will do analysis of alpha data)",
  "Fit summary images available at http://fritz.lns.mit.edu/fits/ for each track (organized by run number)",
  "Jeremy + Igal working on expected angular spectrum from MC. Lead brick is ~1 MFP. Chamber wall is ~0.1 MFP. "
  }; 

//const char * overview[] = {"Uses low threshold $S\\cap J$ playlist", 
//                           "Each cluster is fit using DmtpcProjection and SRIMLineFit()",
//                          "Individual fit plots are available at \\url{http://fritz.lns.mit.edu/teh_fits/}" }; 



const char * mtg4sh_013013[] = 
{
  "Big discrepancy between Evis and Emesh (Slides 2,3)", 
  "How Emesh impacts various things (Slides 4,5,6,7,8,9,10,11,12,13). Not all events have an associated Emesh because they weren't in Jeremy's tree. ", 
  "range / $\\sigma > ~3$ appears to be important. Can we optimize this parameter? (Slides 14,15,16,17). Also, longer ranges at higher $\\sigma$ even though using fitRange (so $\\sigma$ is automatically convolved in)? Not sure if this is a reconstruction artifact or something meaningful." 
}; 

//  slides.AddListSlide("Overview","overview",6, overview); 
  slides.AddListSlide("Comments","comments",3, mtg4sh_013013); 

 // slides.AddPlotSlide(" $\\chi^2/ndof$ ","chisq", chisqc); 
  //slides.AddPlotSlide("Raw $\\phi$ (blue is fit, red is getPhi2)","rawPhi", rawPhi); 
  //slides.AddPlotSlide("Combined $\\phi$ (blue is fit, red is getPhi2)","combined", combined); 
 // slides.AddPlotSlide("fit $\\phi$ components", "split", split); 
//  slides.AddPlotSlide("Fit $\\phi$ vs. implied probability", "phiprob", phiVsProb); 
  slides.AddPlotSlide("Evis vs. Emesh", "veme", visEmeshE); 
  slides.AddPlotSlide("Evis/Emesh", "ratio", ratioc); 
 // slides.AddPlotSlide("Emesh vs. fit range", "rve", rvec); 
  slides.AddPlotSlide("Evis vs. implied probability", "eprob", EVsProb); 
  slides.AddPlotSlide("Emesh vs. implied probability", "meprob", mEVsProb); 
//  slides.AddPlotSlide("Fit $\\phi$ in probability slices ($R >0.8$, $G \\in (0.6,0.8]$, $B \\leq 0.6$), ", "cutPhi", cutPhi); 
//  slides.AddPlotSlide("Combined Fit $\\phi$ in probability slices ($R >0.8$, $G \\in (0.6,0.8]$, $B \\leq 0.6$), ", "cutPhiAll", cutPhiAll); 
  slides.AddPlotSlide("Vis Energy spectrum in probability slices ($R >0.8$, $G \\in (0.6,0.8]$, $B \\leq 0.6$), ", "cutE", cutE); 
  slides.AddPlotSlide("Mesh Energy spectrum in probability slices ($R >0.8$, $G \\in (0.6,0.8]$, $B \\leq 0.6$), ", "cutmE", cutmE); 
  slides.AddPlotSlide("H-T vs Evis in probability slices ($R >0.8$, $G \\in (0.6,0.8]$, $B \\leq 0.6$), ", "htp", htpc); 
  slides.AddPlotSlide("H-T vs Emesh in probability slices ($R >0.8$, $G \\in (0.6,0.8]$, $B \\leq 0.6$), ", "htpm", htpmc); 
  slides.AddPlotSlide("Combined H-T vs Evis in probability slices ($R >0.8$, $G \\in (0.6,0.8]$, $B \\leq 0.6$), ", "htap", htpac); 
  slides.AddPlotSlide("Combined H-T vs Emesh in probability slices ($R >0.8$, $G \\in (0.6,0.8]$, $B \\leq 0.6$), ", "htmap", htpmac); 
//  slides.AddPlotSlide("H-T vs range in probability slices ($R >0.8$, $G \\in (0.6,0.8]$, $B \\leq 0.6$), ", "htpr", htrpc); 
//  slides.AddPlotSlide("Combined H-T vs range in probability slices ($R >0.8$, $G \\in (0.6,0.8]$, $B \\leq 0.6$), ", "htrap", htprac); 
//  slides.AddPlotSlide("H-T vs $\\sigma$ in probability slices ($R >0.8$, $G \\in (0.6,0.8]$, $B \\leq 0.6$), ", "htps", htspc); 
//  slides.AddPlotSlide("Combined H-T vs $\\sigma$ in probability slices ($R >0.8$, $G \\in (0.6,0.8]$, $B \\leq 0.6$), ", "htsap", htpsac); 
  slides.AddPlotSlide("Fit $\\sigma$ vs. Evis", "sige", sigec); 
  slides.AddPlotSlide("Fit $\\sigma$ vs. Emesh", "sigme", sigmec); 
  slides.AddPlotSlide("Fit $r/\\sigma$ vs. probability", "rsVsProb", rsVsProb); 
  slides.AddPlotSlide("H-T vs $r/\\sigma$ in probability slices ($R >0.8$, $G \\in (0.6,0.8]$, $B \\leq 0.6$), ", "htprs", htrspc); 
  slides.AddPlotSlide("Combined H-T vs $r/\\sigma$ in probability slices ($R >0.8$, $G \\in (0.6,0.8]$, $B \\leq 0.6$), ", "htrsap", htprsac); 
//  slides.AddPlotSlide("Fit $\\sigma$", "sig", sigc); 
  slides.AddPlotSlide("Fit $\\sigma$ vs. fit Range", "sigr", sigrc); 
  slides.AddPlotSlide("Fit $\\sigma$ in. chargE slices (R $<$ 100 keV,, G $\\in (100,200]$ keV, B $>$ 200 keV", "sigeslicec", sigeslicec); 
//  slides.AddPlotSlide("$\\phi$ vs. implied probability", "phiVsProb", phiVsProb); 
//  slides.AddPlotSlide("HT vs P in energy slices ($R >150$, $G \\in (75,150]$, $B \\leq 75), ", "hte", htec); 
//  slides.AddPlotSlide("HT vs P in energy slices ($R >150$, $G \\in (75,150]$, $B \\leq 75), ", "htea", hteac); 

  //slides.AddPlotSlide("Track density", "density", densityc); 
 // slides.AddPlotSlide("Track density with meshE", "mdensity", mdensityc); 
 // slides.AddPlotSlide("Average track meshE", "menergy", menergyc); 
//  slides.AddPlotSlide("Average track Evis", "energy", energyc); 
//  slides.AddPlotSlide("Binned $\\phi$ vs. position", "rose", rosec); 
//  slides.AddPlotSlide("H-T fraction vs. position", "htxy", htxyc); 
  slides.Close(); 
} 


