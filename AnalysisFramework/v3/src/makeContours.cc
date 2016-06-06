#include "TApplication.h"
#include "contourPicker.hh"
#include "AnalysisResult.hh"
#include "../../../MaxCam/DmtpcSkimEvent.hh"
#include "TChain.h"
#include "recoilEnergy.hh"
#include "TCanvas.h"
#include "TFile.h"
#include "style.hh"
#include "TH2F.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include <malloc.h>
#include <iostream>

#include "TMarker.h"
#include "TPolyLine.h"


double rms_limit = 80;
double max_limit = 300; 
int nbins = 50; 
double E_limit = 300; 
double pct = 0.68; 

char * merged_dir = "/net/zwicky/dmtpc/analysis/v3/ntuples/";
char * pass_prefix = "/net/zwicky/dmtpc/analysis/v3/pass/pass";

char * get_full_name(char * base)
{
  char * ret = (char*)malloc(strlen(merged_dir) + strlen(base) + 1); 
  strcpy(ret,merged_dir);
  strcpy(ret+strlen(merged_dir),base);
  return ret; 
}

char * concat(char * a, int b) 
{

  TString str(a);
  str+=b;
  char * ret = (char *) malloc(strlen(str.Data())+1); 
  strcpy(ret,str.Data()); 
  return ret;
}


int get_min_cam(int run)
{

  if (run > 100421005 && run <= 100421010) return 1;
  if (run > 100421015 && run <= 100421020) return 1;
  if (run > 100429005 && run <= 100429010) return 1; 
  if (run > 100502000 && run <= 100502004) return 1; 
  return 0; 
}

int convex_hull(int n, double * xin, double * yin, double * xout, double * yout)
{
  double * ysort = (double*)malloc(sizeof(double) * n); 
  double * xsort = (double*)malloc(sizeof(double) * n); 
  double * temp = (double*)malloc(sizeof(double) * n); 
  int * yind = (int*)malloc(sizeof(int) * n); 
  int * tempind = (int*)malloc(sizeof(int) * n); 

  //Sort by y

  TMath::Sort(n,yin,yind,kFALSE); 
  for (int i = 0; i < n; i++)
  {
    ysort[i] = yin[yind[i]]; 
    xsort[i] = xin[yind[i]]; 
  }

  //Find leftmost min

  bool done = false; 
  int j = 1; 
  temp[0] = xsort[0]; 
  while(!done)
  {
      if(ysort[j]!=ysort[j-1]) done=true;
      else
      {
	 temp[j]=xsort[j];
	 j++;
      }
   }
   //cout << j << endl;

   TMath::Sort(j,temp,tempind,kFALSE);

   for(int i=0; i<j; i++)
   {
      xsort[i]=temp[tempind[i]];
   }

   //cout << xsort[0] << "," << ysort[0] << endl;

   int nhull=1;
   double x,y;
   x=xsort[0];
   y=ysort[0];
   xout[0]=xsort[0];
   yout[0]=ysort[0];
   double prevang=-1;

   TString testp;

   done = false;

   while(!done)
   {
      for(int i=0; i<n; i++)
      {
	 if(xsort[i]==x && ysort[i]==y) temp[i]=TMath::TwoPi();
	 else
	 {
	    temp[i]=atan2(ysort[i]-y,xsort[i]-x);
	    if(temp[i] < 0) temp[i]=temp[i]+TMath::TwoPi();
	 }
      }
      
      TMath::Sort(n,temp,tempind,kFALSE);
      
      bool done2=false;
      int k=0;
      while(!done2)
      {
	 if(temp[tempind[k]] > prevang){done2=true;}
	 else k++;
      }

      done2=false;
      int kk=1+k;
      while(!done2)
      {
	 if(kk > n-2){done2=true; break;}
	 if(temp[tempind[kk]]==temp[tempind[k]]){kk++;}
	 else done2=true;
      }

      double dist = 0;
      double ind = k;
      for(int i=k; i<kk; i++)
      {
	 double distt = sqrt(pow(xsort[tempind[i]]-x,2)+
			     pow(ysort[tempind[i]]-y,2));
	 if(distt > dist){ ind=i; dist=distt;}
      }

      k=int(ind);

//      cout << prevang << ":" << temp[tempind[k]] << ":" << k << endl;

      if(xsort[tempind[k]]==xout[0] && ysort[tempind[k]]==yout[0])
      {
	 done=true;
      }
      else
      {
	 prevang = atan2(y-ysort[tempind[k]],x-xsort[tempind[k]]);
	 xout[nhull]=xsort[tempind[k]];
	 yout[nhull]=ysort[tempind[k]];
	 x=xsort[tempind[k]];
	 y=ysort[tempind[k]];
	 nhull++;
	 
//	 cout << x << "," << y << ";" << nhull << endl;
      }
     

   }

   //close hull
   xout[nhull]=xout[0];
   yout[nhull]=yout[0];
   nhull++;

   
   return nhull;
}

int main (int nargs, char ** args)
{

  TApplication * app = new TApplication("make_contours",0,0); 
  AnalysisStyle::setStyle(); 

  int cam = 0;
  if (nargs > 1) cam = atoi(args[1]);
  if (nargs > 2) nbins = atoi(args[2]); 
  if (nargs > 3) pct = atof(args[3]); 

  TChain * mc = new TChain("skim"); 
  DmtpcSkimEvent * ev = 0; 
  mc->SetBranchAddress("event",&ev); 
  if (cam == 0)
  {
    mc->Add(get_full_name("merged100429002_100429005.root"));
    mc->Add(get_full_name("merged100502005_100502008.root"));
  }
  else
  {

    mc->Add(get_full_name("merged100429006_100429010.root"));
    mc->Add(get_full_name("merged100502001_100502004.root"));
  }

  TH2F * rmsplot = new TH2F("rmsplot","rmsplot",nbins,0,E_limit,nbins, 0,rms_limit); 
  TH2F * maxplot = new TH2F("maxplot","maxplot",nbins,0,E_limit,nbins, 0,max_limit); 

  int n = mc->GetEntries(); 
  double * Er = (double*) malloc(sizeof(double) * n); 
  double * rms = (double*) malloc(sizeof(double) * n); 
  double * max = (double*) malloc(sizeof(double) * n); 

  int j = 0; 
  for (int i = 0; i < n; i++)
  {
    mc->GetEntry(i); 
    
    if (ev->ntracks(0) == 0)
    {
      continue;
    }
    Er[j] = (double) RecoilEnergy::getRecoilEnergy(ev->E(0,0),cam); 
    rms[j] = (double) ev->cluster_rms(0,0); 
    max[j] = (double) ev->maxpixel(0,0); 
    rmsplot->Fill(Er[j],rms[j]);
    maxplot->Fill(Er[j],max[j]);
    j++; 
  }

  n = j; 



  int rms_max_bin = rmsplot->GetMaximumBin(); 
  int rms_max_binx = rms_max_bin % (rmsplot->GetNbinsX() + 2);
  int rms_max_biny = ((rms_max_bin - rms_max_binx)/(rmsplot->GetNbinsX()+2)) % (rmsplot->GetNbinsY() + 2);
  double rms_center_x = rmsplot->GetXaxis()->GetBinCenter(rms_max_binx); 
  double rms_center_y = rmsplot->GetYaxis()->GetBinCenter(rms_max_biny); 

  int max_max_bin = maxplot->GetMaximumBin(); 
  int max_max_binx = max_max_bin % (maxplot->GetNbinsX() + 2);
  int max_max_biny = ((max_max_bin - max_max_binx)/(maxplot->GetNbinsX()+2)) % (maxplot->GetNbinsY() + 2);
  double max_center_x = maxplot->GetXaxis()->GetBinCenter(max_max_binx); 
  double max_center_y = maxplot->GetYaxis()->GetBinCenter(max_max_biny); 

  double * rms_dists = (double*) malloc(sizeof(double)*n); 
  double * max_dists = (double*) malloc(sizeof(double)*n); 

  for (int i = 0; i < n; i++)
  {
    double dErms = (Er[i] - rms_center_x); 
    double dEmax = (Er[i] - max_center_x); 
    double dr = (rms[i] - rms_center_y); 
    double dm = (max[i] - max_center_y); 

    double weight_rms = 1./(rmsplot->GetBinContent(rmsplot->FindBin(Er[i],rms[i]))); 

    double weight_max = 1./(maxplot->GetBinContent(maxplot->FindBin(Er[i],max[i]))); 

    if (Er[i] > E_limit || rms[i] > rms_limit) weight_rms = 1e200;  
    if (Er[i] > E_limit || max[i] > max_limit) weight_max = 1e200;  

    rms_dists[i] = weight_rms * TMath::Sqrt(dErms*dErms + dr*dr);     
    max_dists[i] = weight_max * TMath::Sqrt(dEmax*dEmax + dm*dm);          
  } 

 
  int * rms_dist_srt = (int*) malloc(sizeof(int)*n);   
  int * max_dist_srt = (int*) malloc(sizeof(int)*n);   

  TMath::Sort(n,rms_dists,rms_dist_srt,false); 
  TMath::Sort(n,max_dists,max_dist_srt,false); 

  int nCut = (int)(pct*n); 

  double rms_cut_dist = rms_dists[rms_dist_srt[nCut]]; 
  double max_cut_dist = max_dists[max_dist_srt[nCut]]; 
  double rms_in[2][nCut]; 
  double max_in[2][nCut]; 


  int rms_i = 0; 
  int max_i = 0; 
  for (int i = 0; i < n; i++)
  {
    if (rms_dists[i] < rms_cut_dist)
    {
      rms_in[0][rms_i] = Er[i]; 
      rms_in[1][rms_i++] = rms[i]; 
    }

    if (max_dists[i] < max_cut_dist)
    {
      max_in[0][max_i] = Er[i]; 
      max_in[1][max_i++] = max[i]; 
    }
  }

  double rms_hull[2][10000]; 
  double max_hull[2][10000]; 

  cout << "nCut: " << nCut << endl; 
  int n_rms_hull = convex_hull(nCut, rms_in[0], rms_in[1], rms_hull[0], rms_hull[1]); 
  int n_max_hull = convex_hull(nCut, max_in[0], max_in[1], max_hull[0], max_hull[1]); 

  cout << "n_rms_hull: " << n_rms_hull << endl;
  cout << "n_max_hull: " << n_max_hull << endl;


  TCanvas * rmscanvas = new TCanvas("rms","rms",500,500); 
  rmscanvas->cd(); 
  rmsplot->Draw("colz"); 
  TMarker *  rms_cent = new TMarker(rms_center_x, rms_center_y,5); 
  rms_cent->Draw("same"); 
  TPolyLine * rms_line = new TPolyLine(n_rms_hull, rms_hull[0], rms_hull[1]); 
  rms_line->Draw("same"); 

  TCanvas * maxcanvas = new TCanvas("max","max",500,500); 
  maxcanvas->cd(); 
  maxplot->Draw("colz"); 
  TMarker *  max_cent = new TMarker(max_center_x, max_center_y,5); 
  max_cent->Draw("same"); 
  TPolyLine * max_line = new TPolyLine(n_max_hull, max_hull[0], max_hull[1]); 
  max_line->Draw("same"); 

  app->Run(); 


}

