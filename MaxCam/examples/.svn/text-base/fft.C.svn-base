#include "../DmtpcSkimDataset.hh"
#include "../DmtpcSkimEvent.hh"
#include "TH2D.h"
#include "TH1F.h"

bool show_init = false; 
TCanvas * c; 
TTree * intree; 
TH2D * mag0; 
TH2D * ph0; 
TH2D * mag1; 
TH2D * ph1; 
TFile * f = NULL;

void load_file(char * file)
{
 if (f!=NULL) 
 {
     f->Close(); 
     f->Delete();
 }
 f = new TFile(file); 
 mag0 = 0; 
 ph0=0;
 mag1 = 0; 
 ph1=0;
 intree = (TTree*) f->Get("fft"); 
 intree->SetBranchAddress("mag0",&mag0); 
 intree->SetBranchAddress("phase0",&ph0); 
 intree->SetBranchAddress("mag1",&mag1); 
 intree->SetBranchAddress("phase1",&ph1); 
}




void progression(char * prefix, int n, int * runs, int * events)
{
 
  TCanvas* pc = new TCanvas("pc","pc",800,800);   
  pc->Divide(4,n); 
  
  for (int i = 0; i < n; i++)
  {
    TString fs = prefix; 
    fs+=runs[i]; 
    fs+=".root"; 
    TString label = "Run "; 
    label+=runs[i]; 
    label+= " Event "; 
    label+= events[i]; 

    cout << label << endl; 

    load_file(fs.Data()); 
    cout << "here" << endl; 
    intree->GetEntry(events[i]); 
    cout << "her2e" << endl; 

    pc->cd(4*i+1); 
    mag0->Draw("colz"); 
    mag0->GetXaxis()->SetTitle(TString("|FFT(cam0)| ")+label); 

    pc->cd(4*i+2); 
    ph0->Draw("colz"); 
    ph0->GetXaxis()->SetTitle(TString("#phi(FFT(cam0)) ")+label); 

    pc->cd(4*i+3); 
    mag1->Draw("colz"); 
    mag1->GetXaxis()->SetTitle(TString("|FFT(cam1)| ")+label); 

    pc->cd(4*i+4); 
    ph1->Draw("colz"); 
    ph1->GetXaxis()->SetTitle(TString("#phi(FFT(cam1)) ")+label); 
  }

}

void show(int i)
{
  if (!show_init)
  {
    c = new TCanvas("ffts","ffts",800,800); 
    c->Divide(2,2); 
    show_init = true; 
  }

  if (f==NULL)
  {
    cout << "NO FILE LOADED!!!" << endl;
  }
  intree->GetEntry(i); 
  c->cd(1);
  mag0->Draw("colz"); 
  c->cd(2); 
  ph0->Draw("colz"); 
  c->cd(3);
  mag1->Draw("colz"); 
  c->cd(4); 
  ph1->Draw("colz"); 
}


void doit(int ev=10)
{
  DmtpcSkimDataset ds; 
  ds.openRootFile("../../AnalysisFramework/v3/skim/dmtpc_10L_02400skim.root");
  ds.getEvent(ev); 
  TH2F * cleaned = fft_clean(ds.event(),1,128); 

  TCanvas * c = new TCanvas("test","test",800,400); 
  c->Divide(2,1); 
  c->cd(1); 
  ds.event()->cluster(1)->getImage()->DrawCopy("colz"); 
  c->cd(2); 
  cleaned->DrawCopy("colz"); 
  c->Update(); 
}

TH2F * fft_clean(const DmtpcSkimEvent * event, int cam, int bin_width)
{

    TH2F * img = event->cluster(cam)->getImage(); 
    TH2D * mag = new TH2D("mag_temp","mag_temp",
                          img->GetNbinsX(), 
                          img->GetXaxis()->GetXmin(),
                          img->GetXaxis()->GetXmax(),
                          img->GetNbinsY(), 
                          img->GetYaxis()->GetXmin(),
                          img->GetYaxis()->GetXmax()
            ); 
    TH2D * ph = new TH2D("ph_temp","ph_temp",
                          img->GetNbinsX(), 
                          img->GetXaxis()->GetXmin(),
                          img->GetXaxis()->GetXmax(),
                          img->GetNbinsY(), 
                          img->GetYaxis()->GetXmin(),
                          img->GetYaxis()->GetXmax()
            ); 

    mag = (TH2D*) img->FFT((TH1*) mag, "MAG"); 
    TVirtualFFT * fft = TVirtualFFT::GetCurrentTransform(); 
    
    
    double * re = new double[img->GetNbinsX()*img->GetNbinsY()+2]; 
    double * im = new double[img->GetNbinsX()*img->GetNbinsY()+2]; 

    fft->GetPointsComplex(re,im);   

    cout << re[50] << " " << im[50] << endl; 

    ph = (TH2D*) img->FFT((TH1*) ph, "MAG"); 

//    mag->DrawCopy("colz"); 

//    return NULL;

    for (int c = 1; c <= img->GetNbinsX()/bin_width; c++)
    {
        TH1D * proj = (TH1D*) (mag->ProjectionY("_py",(c-1)*bin_width+1,c*bin_width))->Clone(); 
        proj->Scale(1./bin_width); 

        TSpectrum spec(20); 

        spec.Search(proj); 

        int npeaks = spec.GetNPeaks(); 

        double * params = new double[3*npeaks+1]; 
        
        TString fstr = "[0] + "; 
        params[0] = 0; 

        for (int i = 0; i < npeaks; i++)
        {
          params[3*i+1] = spec.GetPositionY()[i];  
          params[3*i+2] = spec.GetPositionX()[i];  
          params[3*i+3] = 3; 

          fstr+= "gaus(";
          fstr+=3*i+1; 
          fstr+= ")"; 
          if (i < npeaks-1) fstr+= " + "; 
        }

        TF1 * fitfunc = new TF1("fitfunc",fstr,proj->GetXaxis()->GetXmin(),
                                 proj->GetXaxis()->GetXmax()); 


        fitfunc->SetParameters(params); 

        for (int i = 0; i < npeaks; i++)
        {
          fitfunc->SetParLimits(3*i+1,0.9*params[3*i+1],1.1*params[3*i+1]); 
          fitfunc->SetParLimits(3*i+2,params[3*i+2]-4,params[3*i+2]+4); 
          fitfunc->SetParLimits(3*i+3,0,10); 
        }

        proj->Fit(fitfunc); 

        for (int x = (c-1)* bin_width+1 ; x <= (c)*bin_width; x++)
        {
          for (int y = 1; y <= proj->GetNbinsX(); y++)
          {
            double val = fitfunc->Eval(proj->GetBinCenter(y)) - fitfunc->GetParameter(0); 
            double ang = ph->GetBinContent(x,y); 

            re[x+y*(proj->GetNbinsX()+2)] -= val * TMath::Cos(ang); 
            im[x+y*(proj->GetNbinsX()+2)] -= val * TMath::Sin(ang); 
            
          }
        } 


        proj->Delete(); 
    }

    cout << re[50] << " " << im[50] << endl; 

    int dim_arr[2]; 
    dim_arr[0] = mag->GetNbinsX();
    dim_arr[1] = mag->GetNbinsY();

    TVirtualFFT * ifft = TVirtualFFT::FFT(2,dim_arr,"C2R M K"); 
    ifft->SetPointsComplex(re,im); 
    ifft->Transform();

    TH2F * ret = new TH2F("fft_out","fft_out",img->GetNbinsX(), 
                          img->GetXaxis()->GetXmin(),
                          img->GetXaxis()->GetXmax(),
                          img->GetNbinsY(), 
                          img->GetYaxis()->GetXmin(),
                          img->GetYaxis()->GetXmax());
              
    ret =  (TH2F*)TH2F::TransformHisto(ifft,(TH1*)ret,"RE"); 

    ret->Scale(1./(ret->GetNbinsY()*ret->GetNbinsY()));
    mag->Delete();
    ifft->Delete();
    ph->Delete(); 
    delete im; 
    delete re; 
    return ret; 

}


void do_fft(char * infile, char * outfile, int nevents=-1)
{
  DmtpcSkimDataset ds; 
  ds.openRootFile(infile); 
  ds.tree()->SetBranchStatus("*trigger*",0); 
  TFile ff(outfile,"RECREATE"); 
  TH2D * mag0 = new TH2D("mag0","mag0",256,0,256,256,0,256); 
  TH2D * phase0 = new TH2D("phase0","phase0",256,0,256,256,0,256); 
  TH2D * mag1 = new TH2D("mag1","mag1",256,0,256,256,0,256); 
  TH2D * phase1 = new TH2D("phase1","phase1",256,0,256,256,0,256); 
  TTree * outtree = new TTree("fft","fft"); 
  outtree->Branch("mag0",&mag0); 
  outtree->Branch("phase0",&phase0); 
  outtree->Branch("mag1",&mag1); 
  outtree->Branch("phase1",&phase1); 
  if (nevents < 0) nevents = ds.nevents(); 
  for (int i = 0; i < nevents; i++)
  {
    ds.getEvent(i); 
    mag0 = (TH2D*) ds.event()->cluster(0)->getImage()->FFT((TH1*)mag0,"MAG"); 
    phase0 = (TH2D*) ds.event()->cluster(0)->getImage()->FFT((TH1*)phase0,"PH"); 
    mag1 = (TH2D*) ds.event()->cluster(1)->getImage()->FFT((TH1*)mag1,"MAG"); 
    phase1 = (TH2D*) ds.event()->cluster(1)->getImage()->FFT((TH1*)phase1,"PH"); 
    outtree->Fill(); 
  }
  outtree->Write();
  ff.Close(); 
}


