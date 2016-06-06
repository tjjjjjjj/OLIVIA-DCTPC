void make_cutoff_plots(TString cam_id, bool ci=false, bool save = false, TString format = ".pdf")
{

   gStyle->SetOptTitle(1);
   gStyle->SetOptStat(1);
   TString eval = !ci ? TString("weights/Eval_")+cam_id+TString("_worms.root") 
                      :TString("weights/Eval_")+cam_id+TString("_worms_ci.root"); 

   TCanvas *  c = new TCanvas("cutoff_c","cutoff_c",600,1000); 
   
   TFile f(eval); 
   c->Divide(1,3); 
   c->cd(1); 
   c->GetPad(1)->SetLogy(); 
   TH1F * sig_value = (TH1F*) f.Get("sig_value"); 
   sig_value->SetTitle("Signal and Background Classifier Distributions"); 
   sig_value->GetXaxis()->SetTitle("Classifier Value"); 
   sig_value->GetYaxis()->SetTitle("N"); 
   sig_value->SetLineColor(3); 
   sig_value->DrawCopy(); 
   TH1F * bg_value = (TH1F*) f.Get("bg_value"); 
   bg_value->SetLineColor(2); 
   bg_value->DrawCopy("same"); 

   c->cd(2); 
   c->GetPad(2)->SetLogy(); 

   TGraph * sig_efficiency = (TGraph*) f.Get("signal_efficiency"); 
   sig_efficiency->SetLineColor(3); 
   sig_efficiency->Draw("alp"); 
   sig_efficiency->SetTitle("Signal and Background Efficiency"); 
   sig_efficiency->GetXaxis()->SetTitle("Classifier Value"); 
   sig_efficiency->GetYaxis()->SetTitle("Efficiency"); 
   TGraph * background_eff_cmp = (TGraph*) f.Get("background_eff_cmp"); 
   background_eff_cmp->SetLineColor(2); 
   background_eff_cmp->Draw("lpsame"); 

   c->cd(3); 
   TH2F * value_v_E = (TH2F*) f.Get("value_v_E"); 
   value_v_E->SetTitle("Signal Classifier Value vs. Energy"); 
   value_v_E->GetYaxis()->SetTitle("Signal Classifier Value"); 
   value_v_E->GetXaxis()->SetTitle("Signal Energy (adu)"); 
   value_v_E->DrawCopy("colz"); 

   if (save)
   {
    c->SaveAs(cam_id + (ci ? TString("_ci_") : TString("_")) + TString("_classifier")+format); 
   }

   f.Close(); 
}
void save_all_cutoff_plots(TString format=".pdf")
{
  make_cutoff_plots("100439",false,true,format); 
  make_cutoff_plots("100439",true,true,format); 
  make_cutoff_plots("081264",false,true,format); 
  make_cutoff_plots("081264",true,true,format); 
  make_cutoff_plots("A80334",false,true,format); 
  make_cutoff_plots("A80334",true,true,format); 
}


void save_all_parameter_plots(TString format=".pdf")
{
  make_parameter_plots("100439",false,true,format); 
  make_parameter_plots("100439",true,true,format); 
  make_parameter_plots("081264",false,true,format); 
  make_parameter_plots("081264",true,true,format); 
  make_parameter_plots("A80334",false,true,format); 
  make_parameter_plots("A80334",true,true,format); 
}

void make_parameter_plots(TString cam_id, bool ci = false, bool save = false, TString format = ".pdf" )
{

   TString tmva = !ci ? TString("weights/TMVA_")+cam_id+TString("_worms.root") 
                      :TString("weights/TMVA_")+cam_id+TString("_worms_ci.root"); 

   TCanvas * c = new TCanvas("par_c","par_c", 1600,400); 

   c->Divide(7,2); 

   TFile f(tmva); 

   TDirectoryFile * d = (TDirectoryFile*) f.Get("InputVariables_Id"); 

   c->cd(1); 
   TH1F * log_E___Signal_Id = (TH1F*) d->Get("log_E___Signal_Id"); 
   log_E___Signal_Id->SetLineColor(3); 
   log_E___Signal_Id->DrawCopy(); 
   c->cd(2); 
   TH1F * log_range___Signal_Id = (TH1F*) d->Get("log_range___Signal_Id"); 
   log_range___Signal_Id->SetLineColor(3); 
   log_range___Signal_Id->DrawCopy(); 
   c->cd(3); 
   TH1F * log_cluster_rms___Signal_Id = (TH1F*) d->Get("log_cluster_rms___Signal_Id"); 
   log_cluster_rms___Signal_Id->SetLineColor(3); 
   log_cluster_rms___Signal_Id->DrawCopy(); 
   c->cd(4); 
   TH1F * log_maxpixel___Signal_Id = (TH1F*) d->Get("log_maxpixel___Signal_Id"); 
   log_maxpixel___Signal_Id->SetLineColor(3); 
   log_maxpixel___Signal_Id->DrawCopy(); 
   c->cd(5); 
   TH1F * neighbors__Signal_Id = (TH1F*) d->Get("neighbors__Signal_Id"); 
   neighbors__Signal_Id->SetLineColor(3); 
   neighbors__Signal_Id->DrawCopy(); 
   c->cd(6); 
   TH1F * log_npixel___Signal_Id = (TH1F*) d->Get("log_npixel___Signal_Id"); 
   log_npixel___Signal_Id->SetLineColor(3); 
   log_npixel___Signal_Id->DrawCopy(); 
   c->cd(7); 
   TH1F * log_npixel_red___Signal_Id = (TH1F*) d->Get("log_npixel_red___Signal_Id"); 
   log_npixel_red___Signal_Id->SetLineColor(3); 
   log_npixel_red___Signal_Id->DrawCopy(); 

   c->cd(8); 
   TH1F * log_E___Background_Id = (TH1F*) d->Get("log_E___Background_Id"); 
   log_E___Background_Id->SetLineColor(2); 
   log_E___Background_Id->DrawCopy(); 
   c->cd(9); 
   TH1F * log_range___Background_Id = (TH1F*) d->Get("log_range___Background_Id"); 
   log_range___Background_Id->SetLineColor(2); 
   log_range___Background_Id->DrawCopy(); 
   c->cd(10); 
   TH1F * log_cluster_rms___Background_Id = (TH1F*) d->Get("log_cluster_rms___Background_Id"); 
   log_cluster_rms___Background_Id->SetLineColor(2); 
   log_cluster_rms___Background_Id->DrawCopy(); 
   c->cd(11); 
   TH1F * log_maxpixel___Background_Id = (TH1F*) d->Get("log_maxpixel___Background_Id"); 
   log_maxpixel___Background_Id->SetLineColor(2); 
   log_maxpixel___Background_Id->DrawCopy(); 
   c->cd(12); 
   TH1F * neighbors__Background_Id = (TH1F*) d->Get("neighbors__Background_Id"); 
   neighbors__Background_Id->SetLineColor(2); 
   neighbors__Background_Id->DrawCopy(); 
   c->cd(13); 
   TH1F * log_npixel___Background_Id = (TH1F*) d->Get("log_npixel___Background_Id"); 
   log_npixel___Background_Id->SetLineColor(2); 
   log_npixel___Background_Id->DrawCopy(); 
   c->cd(14); 
   TH1F * log_npixel_red___Background_Id = (TH1F*) d->Get("log_npixel_red___Background_Id"); 
   log_npixel_red___Background_Id->SetLineColor(2); 
   log_npixel_red___Background_Id->DrawCopy(); 

   if (save)
   {
      c->SaveAs(cam_id + (ci ? TString("_ci_") : TString("_")) + TString("_parameters")+format); 
   }

}
