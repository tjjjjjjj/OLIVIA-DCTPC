TCanvas * phiVsX(TH2 * hist /* phi on y axis*/ , double correctPhi = 0, double range_factor = 1.5 , bool draw_fits = true) 
{
  TCanvas * c = new TCanvas(TString::Format("%s_c",hist->GetName()), TString::Format("#phi vs. %s", hist->GetXaxis()->GetTitle()), 1800 ,1000);

  c->Divide(3,2); 
  c->cd(1); 
  hist->DrawCopy("colz"); 
  c->cd(2); 


  int ncolors = gStyle->GetNumberOfColors(); 
  int dz= hist->GetMaximum() - hist->GetMinimum(); 

  double Xmax = hist->GetXaxis()->GetXmax(); 
  double Xrange =range_factor * Xmax; 
  gPad->Range(-Xrange, -Xrange, Xrange, Xrange) ;

//  ax->Draw("same"); 
  for (int i = 1; i <= hist->GetNbinsX(); i++)
  {
    for (int j = 1; j <= hist->GetNbinsY(); j++)
    {
      double val = hist->GetBinContent(i,j); 
      if (val == 0) continue; 
      int color = gStyle->GetColorPalette(((val-hist->GetMinimum())/dz)  * double(ncolors)-1); 

      double phi_min = hist->GetYaxis()->GetBinLowEdge(j) * 180 / TMath::Pi();  
      double phi_max = hist->GetYaxis()->GetBinLowEdge(j+1) * 180 / TMath::Pi() ; 

      double X_min = hist->GetXaxis()->GetBinLowEdge(i) ; 
      double X_max = hist->GetXaxis()->GetBinLowEdge(i+1) ; 

      TCrown * crown = new TCrown(0,0, X_min,X_max, phi_min,phi_max); 
      crown->SetFillColor(color); 
      crown->Draw(); 
    }
  }

  TGaxis * r_axis = new TGaxis(0,0, Xmax * cos(correctPhi), Xmax * sin(correctPhi),0,Xmax); 
  r_axis->SetNdivisions(hist->GetNbinsX()); 
  r_axis->SetLabelSize(0.02); 
  r_axis->Draw(); 

  TArrow * arrow = new TArrow(Xmax * cos(correctPhi), Xmax * sin(correctPhi), range_factor * Xmax * cos(correctPhi), range_factor * Xmax * sin(correctPhi)); 
  arrow->Draw(); 

  TText * r_label = new TText(1.04*Xmax * cos(correctPhi) , 1.02 * Xmax * sin(correctPhi) ,  hist->GetXaxis()->GetTitle()); 
  r_label->SetTextAngle(180 * correctPhi / TMath::Pi()); 
  r_label->SetTextSize(0.02); 
  r_label->SetTextAlign(11); 
  r_label->SetTextColor(2); 
  r_label->Draw(); 

  TLatex * title = new TLatex(0, Xrange, hist->GetTitle()); 
  title->SetTextSize(0.05); 
  title->SetTextAlign(23); 

  title->Draw(); 
  
  for (int i = 0; i < hist->GetNbinsY(); i++)
  {
    int denom = hist->GetNbinsY(); 
    int num = i; 
    double phi = (TMath::Pi() * 2 *num) / denom; 

    TLatex * angle_label = new TLatex (1.01 *Xmax * cos(phi) , 1.01*Xmax * sin(phi) , TString::Format("%d^{o}",num* 360 / denom)); 
    angle_label->SetTextSize(0.02); 
    angle_label->SetTextColor(22); 
    angle_label->SetTextAngle(180 *phi / TMath::Pi() - 90); 
    angle_label->Draw(); 
  }

  double xmin = cos(correctPhi) > 0  ? -Xrange: 0.95 * Xrange;
  double xmax = cos(correctPhi) > 0  ? -Xrange + 0.05 : Xrange;

//  TPaletteAxis * palette = new TPaletteAxis( xmin, -Xrange, xmax, Xrange, 0); 
//  palette->Draw(); 

  gPad->SetEditable(0); 
  c->cd(3); 

  TH1I all("allhist", "allhist", hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax()); 
  TH1I same("samehist", "samehist", hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax()); 

  for (int i = 1; i<= hist->GetNbinsX(); i++)
  {
    int total = 0; 
    int sameside = 0; 
    for (int j = 1; j<= hist->GetNbinsY(); j++)
    {
       total += hist->GetBinContent(i,j);  

       if (cos(correctPhi - hist->GetYaxis()->GetBinCenter(j)) > 0) 
       {
          sameside += hist->GetBinContent(i,j); 
       }
    }

    all.SetBinContent(i,total); 
    same.SetBinContent(i,sameside); 
  }

  TH2C axismaker("axismaker","H-T fraction", 10,hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(),10, 0.4, 1); 
  axismaker.GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle()); 
  axismaker.GetYaxis()->SetTitle("Fraction correct (same hemisphere)"); 
  axismaker.GetYaxis()->SetNdivisions(406); 
  axismaker.DrawCopy(); 
  TGraphAsymmErrors*  g = new TGraphAsymmErrors(all.GetNbinsX()); 
  g->Divide(&same,&all); 
  g->SetMarkerColor(4); 
  g->Draw("cpsame"); 
  g->SetTitle("Head-Tail Fraction"); 


  vector<double> vm_fracs(hist->GetNbinsX()); 
  vector<double> vm_fracs_errors(hist->GetNbinsX()); 
  vector<double> vm_mean(hist->GetNbinsX()); 
  vector<double> vm_mean_errors(hist->GetNbinsX()); 
  vector<double> vm_k(hist->GetNbinsX()); 
  vector<double> vm_k_errors(hist->GetNbinsX()); 
  vector<double> vm_X(hist->GetNbinsX()); 
  vector<double> vm_X_errors(hist->GetNbinsX()); 

  //do fits to von mises distribution for each slice of X 
  
  c.cd(4); 
  TF1 * vmf = new TF1("vmf",DmtpcMath::vonMisesDistHT, 0, 2*TMath::Pi(), 4);

  TCanvas * fits = draw_fits ? new TCanvas(TString::Format("fits_%s",hist->GetName()),"fits", 1800,1000): 0; 
  if (draw_fits)
  {
    fits->Divide((hist->GetNbinsX()+1)/4,4); 
    fits->cd(1); 
  }
  for (int i = 1; i<= hist->GetNbinsX(); i++)
  {
    TH1 * proj = hist->ProjectionY(TString::Format("_py%s_%d",hist->GetName(),i),i,i); 
    vmf->SetParameter(0,0.75); 
    vmf->SetParError(0,0.25); 
    vmf->SetParLimits(0,0.5,1); 
    vmf->SetParameter(1,proj->Integral()); 
    vmf->SetParError(1,10); 
    vmf->SetParLimits(1,0,10*proj->Integral()); 
    vmf->SetParameter(2,correctPhi); 
    vmf->SetParError(2,0.05); 
    vmf->SetParameter(3,8); 
    vmf->SetParError(3,2); 
    proj->Fit(vmf,"LIM"); 

    vm_fracs[i-1] = vmf->GetParameter(0); 
    vm_fracs_errors[i-1] = vmf->GetParError(0); 
    vm_mean[i-1] = DmtpcMath::normalizeAngle(vmf->GetParameter(2)); 
    vm_mean_errors[i-1] = vmf->GetParError(2); 
    vm_k[i-1] = fabs(vmf->GetParameter(3)); 
    vm_k_errors[i-1] = vmf->GetParError(3); 
    vm_X[i-1] = hist->GetBinCenter(i); 
    vm_X_errors[i-1] = hist->GetBinWidth(i)/2; 


    if (draw_fits)
    {
      fits->cd(i); 
      proj->Draw("pol"); 
//      vmf.DrawCopy("same"); 
    }


  }

  TGraphErrors * gfracs = new TGraphErrors(hist->GetNbinsX(), &(vm_X[0]), &(vm_fracs[0]), &(vm_X_errors[0]),&(vm_fracs_errors[0])); 
  gfracs->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle()); 
  gfracs->GetYaxis()->SetTitle("Fraction correct (from Von-Mises fit)"); 
  gfracs->Draw("acp"); 

  c.cd(5); 
  TGraphErrors*  gmeans = new TGraphErrors(hist->GetNbinsX(), &(vm_X[0]), &(vm_mean[0]), &(vm_X_errors[0]),&(vm_mean_errors[0])); 
  gmeans->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle()); 
  gmeans->GetYaxis()->SetTitle("mean direction (from Von-Mises fit)"); 
  gmeans->Draw("acp"); 

  c.cd(6); 
  TGraphErrors*  gk = new TGraphErrors(hist->GetNbinsX(), &(vm_X[0]), &(vm_k[0]), &(vm_X_errors[0]),&(vm_k_errors[0])); 
  gmeans->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle()); 
  gmeans->GetYaxis()->SetTitle("mean k (from Von-Mises fit)"); 
  gmeans->Draw("acp"); 


  return c; 

}
