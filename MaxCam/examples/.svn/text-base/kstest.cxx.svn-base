#include <iostream>
#include <math.h>
#include "../DmtpcDataset.hh"
#include "TH2F.h"
using namespace std;

TH2F* TH2Fks(DmtpcDataset* d, int cameraNum, int start, int stop, float sigmas, int iter, float change){

  //determine size of images
  d->getEvent(start);
  TH2F* pixhisto=d->event()->ccdData(cameraNum);
  int xpix=pixhisto->GetNbinsX();
  int ypix=pixhisto->GetNbinsY();

  TH2F *kshisto = new TH2F("kshisto", "kshisto", xpix, 0, xpix, ypix, 0, ypix);

  //to be used in the next loop
  int n = stop-start+1;
  double dp[n];
  sigmas = sqrt(sigmas);
	
  //loop over all pixels of all images, ks clipping

   for (int ii=0; ii<xpix; ii++){
     cout<<"xpix"<<endl;
     for (int jj=0; jj<ypix; jj++){
       cout<<"ypix"<<endl;
       float sum=0;
       float sumsq=0;
       float m;
       float changetest;
       int k = 0, r;

       for (int i=start; i<=stop; i++){
	 //cout<<"image"<<endl;
	 d->getEvent(i);
	 TH2F* dhisto= d->event()->ccdData(cameraNum);

	 double content=dhisto->GetBinContent(ii, jj);
	 sum=sum+content;
	 sumsq=sumsq+sqrt(content);
	 dp[i-start]=content;
     }
 
	float s = sigmas * (sumsq / (1.0 * n) - sqrt(sum / (1.0 * n)));
       
	if (!(n%2)){
	  m=(dp[int(n/2.0-1)]+dp[int(n/2.0)])/2;
	}
	else {
	  m=dp[int(n/2.0-0.5)];
	}
  
	do {
	  cout<<"do"<<endl;
		iter --;
		sum = 0.0;
		sumsq = 0.0;
		r = k;
		k = 0;
		changetest=m;

		for (int i=0; i<n; i++) {
			if (sqrt(dp[i] - m) < s) {
				sum += dp[i];
				sumsq += sqrt(dp[i]);
				k++;
			}
		}
		if (k == 0){
			break;
		}
		m = sum / (1.0 * (k));
	    
		s = sigmas * (sumsq / (1.0 * (k)) - sqrt(m));
	} while ((iter > 0) && (k != r) && (fabs(changetest-m)>change));

     kshisto->SetBinContent(ii, jj, m);
   }
  }
  
  return kshisto;
}

int kstest(){
  DmtpcDataset e;
  e.openRootFile("proem_a0.root");
  DmtpcDataset* d;
  d=&e;

  TH2F* khisto=TH2Fks(d, 0, 1, 9, 2, 30, 0.01);

  return 0;
}
