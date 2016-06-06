//=================================================
// class Tools
//=================================================
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <iostream>

#include "Tools.hh"

using namespace std;

//ClassImp(Tools)

Tools::Tools()
{
  // constructor: not much to do here... 
}

void 
Tools::Draw4Hist(TH1F* h1,TH1F* h2,TH1F* h3,TH1F* h4,TString myCanName="myGenCan")
{
  TCanvas* myCan = new TCanvas(myCanName,myCanName, 600,900);; 
  myCan->Clear();
  myCan->Divide(1,2);
  
  h1->SetLineColor(1);
  h2->SetLineColor(2);
  
  h3->SetLineColor(1);
  h4->SetLineColor(2);
  
  myCan->cd(1);
  h1->Draw("histo");

  h2->Draw("SAME,histo");
    
  myCan->cd(2);
  h3->Draw("histo");
  h4->Draw("SAME,histo");

}

void
Tools::DrawTwoHistosNorm(const TH1F* h1, const TH1F* h2, TString canvasName ) 
{ 

  // Draw 2 histos on top of each other normalizing to the same numb of entries 

  double nE1 = h1->GetEntries(); 
  double nE2 = h2->GetEntries();
  double norm = nE1/nE2; 

  double max1 = h1->GetMaximum();
  double max2 = h2->GetMaximum();

  max2 *= norm; 

  double themax = max1;
  if (max2>themax) themax = max2;
 
  TH1F* hsig1 = (TH1F*)h1->Clone(); hsig1->Sumw2(); 
  TH1F* hsig2 = (TH1F*)h2->Clone(); hsig2->Sumw2(); 
  hsig1->SetMaximum(themax*1.1);

  hsig1->SetLineColor(2); hsig1->SetStats(false); 
  hsig2->SetLineColor(1); hsig2->SetStats(false); 

  // Draw histos 
  // -----------  
  TCanvas* canvas = new TCanvas(canvasName,canvasName,400,400); // open a window
  canvas->Divide(1,1);
  canvas->cd(1);
  hsig1->Draw("histo");
  hsig2->Scale(norm);
  hsig2->Draw("histo,same");

  TString name = canvasName + ".eps"; 
  canvas->Print(name);

}


void
Tools::OptimizeCut2(const TH1F* h1, const TH1F* h2, float lumi1=1.0 ,float lumi2=1.0,TString myCan="testCan")
{ 

  TH1F* hsig1o = (TH1F*)h1->Clone(); hsig1o->Sumw2(); 
  TH1F* hsig2o = (TH1F*)h2->Clone(); hsig2o->Sumw2(); 

  TH1F* hsig1 = (TH1F*)h1->Clone(); hsig1->Sumw2();
  TH1F* hsig2 = (TH1F*)h2->Clone(); hsig2->Sumw2();

  double nsig, nbkg;
  nsig = 0;
  nbkg = 0;
  for (int i = 0; i <= h1->GetNbinsX()+1; i++) {
    nsig += h1->GetBinContent(i)/lumi1;
    nbkg += h2->GetBinContent(i)/lumi2;
    double signif = 0;
    if (nsig+nbkg > 0) signif = nsig/sqrt(nsig+nbkg);
    hsig1->SetBinContent(i,signif);
  }
  nsig = 0;
  nbkg = 0;
  for (int i = h1->GetNbinsX()+1; i >= 0; i--) {
    nsig += h1->GetBinContent(i)/lumi1;
    nbkg += h2->GetBinContent(i)/lumi2;
    double signif = 0;
    if (nsig+nbkg > 0) signif = nsig/sqrt(nsig+nbkg);
    hsig2->SetBinContent(i,signif);
  }

  Draw4Hist(hsig1o,hsig2o,hsig1,hsig2,myCan); 
  
  PrintSignificance(hsig1,hsig2);  

//   TString name = cutname + "_optimize.eps"; 
//   myCan->Print(name);

}


void
Tools::OptimizeAndCompareWData(const TH1F* h1, const TH1F* h2, const TH1F* h3, const TH1F* h4,TString myCanName)
{ 
  // all histos are expected to be normalized to the same luminosity
  // h1= signal MC 
  // h2= generic MC 
  // h3= offpeak DATA
  // h4 = onpeak DATA

  TH1F* hsig1o = (TH1F*)h1->Clone(); hsig1o->Sumw2(); 
  TH1F* hsig2o = (TH1F*)h2->Clone(); hsig2o->Sumw2(); 

  TH1F* hsig1 = (TH1F*)h1->Clone(); hsig1->Sumw2();
  TH1F* hsig2 = (TH1F*)h2->Clone(); hsig2->Sumw2();

  double nsig, nbkg;
  nsig = 0;
  nbkg = 0;
  for (int i = 0; i <= h1->GetNbinsX()+1; i++) {
    nsig += h1->GetBinContent(i);
    nbkg += h2->GetBinContent(i);
    double signif = 0;
    if (nsig+nbkg > 0) signif = nsig/sqrt(nsig+nbkg);
    hsig1->SetBinContent(i,signif);
  }
  nsig = 0;
  nbkg = 0;
  for (int i = h1->GetNbinsX()+1; i >= 0; i--) {
    nsig += h1->GetBinContent(i);
    nbkg += h2->GetBinContent(i);
    double signif = 0;
    if (nsig+nbkg > 0) signif = nsig/sqrt(nsig+nbkg);
    hsig2->SetBinContent(i,signif);
  }


  PrintSignificance(hsig1,hsig2);  
  //  Draw4Hist(hsig1o,hsig2o,hsig1,hsig2,myCanName); 
  

  TCanvas* myCan = new TCanvas(myCanName,myCanName,600,900);
  myCan->Clear();  myCan->Divide(2,2);
  
  hsig1o->SetLineColor(2);
  hsig2o->SetLineColor(1);
  
  hsig1->SetLineColor(2);
  hsig2->SetLineColor(1);
  

  // autoscaling y axis
  double ymax = hsig1o->GetMaximum(); if (ymax<hsig2o->GetMaximum()) ymax = hsig2o->GetMaximum();  hsig1o->SetMaximum(ymax*1.1);
  myCan->cd(1);  hsig1o->Draw("histo");  hsig2o->Draw("SAME,histo");
    
  ymax = hsig1->GetMaximum(); if (ymax<hsig2->GetMaximum()) ymax = hsig2->GetMaximum();  hsig1->SetMaximum(ymax*1.1);
  myCan->cd(3);  hsig1->Draw("histo");  hsig2->Draw("SAME,histo");


  //  TH1F* h3copy = (TH1F*)h3->Clone(); h3copy->Sumw2();  h3copy->SetLineColor(3); h3copy->SetMarkerColor(3);// off peak not used
  TH1F* h4copy = (TH1F*)h4->Clone(); h4copy->Sumw2();  h4copy->SetLineColor(4); h4copy->SetMarkerColor(4);

  // Normalize data to the integral of the background MC 
  double intBack = hsig2o->Integral(); double intData = h4copy->Integral();   h4copy->Scale(intBack/intData);

  // set max and draw 
  ymax = hsig2o->GetMaximum(); if (ymax<h4copy->GetMaximum()) ymax = h4copy->GetMaximum();  hsig2o->SetMaximum(ymax*1.1);
  myCan->cd(2);  hsig2o->Draw("histo");  hsig1o->SetFillColor(2); hsig1o->Draw("histo,same");  h4copy->Draw("SAME");

  // now take the ratio of data and MC 
  myCan->cd(4); 

  TH1F* hdcopy = (TH1F*)h4copy->Clone("hdcopy");
  hdcopy->Divide(hsig2o);
//   hdcopy->SetMaximum(1.5); 
//   hdcopy->SetMinimum(0.5); 

  hdcopy->SetStats(false); 
  hdcopy->SetTitle("DATA/MC"); 
  hdcopy->Draw();



  TString name = myCanName + "_optimize.eps"; 
  myCan->Print(name);

} 



void
Tools::OptimizeCut3(const TH1F* h1, const TH1F* h2, const TH1F* h3, float lumi1=1.0 ,float lumi2=1.0,float lumi3=1.0,
                   TString myCan="testCan", bool draw=false )
{ 

  TH1F* h1copy = (TH1F*)h1->Clone(); h1copy->Sumw2(); 
  TH1F* h2copy = (TH1F*)h2->Clone(); h2copy->Sumw2(); 
  TH1F* h3copy = (TH1F*)h3->Clone(); h3copy->Sumw2(); 

  TH1F* hsig1 = (TH1F*)h1->Clone(); hsig1->Sumw2();
  TH1F* hsig2 = (TH1F*)h2->Clone(); hsig2->Sumw2();
  TH1F* hsig3 = (TH1F*)h3->Clone(); hsig3->Sumw2();

  double nsig, nbkg;
  nsig = 0;
  nbkg = 0;
  for (int i = 0; i <= h1->GetNbinsX()+1; i++) {
    nsig += h1->GetBinContent(i)/lumi1;
    nbkg += h2->GetBinContent(i)/lumi2;
    nbkg += h3->GetBinContent(i)/lumi3;
    double signif = 0;
    if (nsig+nbkg > 0) signif = nsig/sqrt(nsig+nbkg);
    hsig1->SetBinContent(i,signif);
  }
  nsig = 0;
  nbkg = 0;
  for (int i = h1->GetNbinsX()+1; i >= 0; i--) {
    nsig += h1->GetBinContent(i)/lumi1;
    nbkg += h2->GetBinContent(i)/lumi2;
    nbkg += h3->GetBinContent(i)/lumi3;
    double signif = 0;
    if (nsig+nbkg > 0) signif = nsig/sqrt(nsig+nbkg);
    hsig2->SetBinContent(i,signif);
  }

  // renormalize histos 
   h1copy->Scale(1./lumi1); 
   h2copy->Scale(1./lumi2); 
   h3copy->Scale(1./lumi3); 
  
  // add backgrounds 
  h2copy->Add(h3copy); 

  if(draw) Draw4Hist(h1copy,h2copy,hsig1,hsig2,myCan); 

//   double bestS = hsig1->GetMaximum(); 
//   double bestS2 = hsig2->GetMaximum(); 
//   cout << " Maximum sensitivity: integrating L->R: " << bestS << "; integrating R->L: " << bestS2<< endl; 

  PrintSignificance(hsig1,hsig2);  

  TCanvas* c1 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(myCan);
  TString name = myCan + ".eps"; 

  c1->Print(name);

} 



void
Tools::RejVsEff(int nSets, const TH1F* h1, const TH1F* h2, const TH1F* h3,const TH1F* h4, float lumi1=1.0 ,float lumi2=1.0,float lumi3=1.0,float lumi4=1.0,TString myCan="testCan", bool draw=false )
{ 

  double n1 = 0;
  double n2 = 0;
  double n3 = 0;
  double n4 = 0;

  int nBins1 = h1->GetNbinsX(); 
  int nBins2 = h2->GetNbinsX(); 
  int nBins3 = 0; 
  if (h3!=NULL)  nBins3= h3->GetNbinsX(); 

  int nBins4 = 0; 
  if (h4!=NULL)  nBins4= h4->GetNbinsX(); 

  if (nBins1!= nBins2) {
    cout << " >>>> EffVsEff Problems: histos have different binning!" << endl; 
    return; 
  }

  double nEnt1 = h1->GetEntries(); 
  double nEnt2 = h2->GetEntries(); 
  double nEnt3 = 0;  if (h3!=NULL)  nEnt3 =h3->GetEntries(); 
  double nEnt4 = 0;  if (h4!=NULL)  nEnt4 =h4->GetEntries(); 

  //  cout << " nbins = " << nBins1; 

  double *eff_1, *eff_2, *eff_3, *eff_4;
  eff_1 = new double[nBins1+2];
  eff_2 = new double[nBins1+2];
  eff_3 = new double[nBins1+2];
  eff_4 = new double[nBins1+2];


  for (int i = h1->GetNbinsX()+1; i>=0; i--) {           // nbins+1 to cover overflows 
    n1 += h1->GetBinContent(i);    eff_1[i] = n1 / nEnt1; 
    n2 += h2->GetBinContent(i);    eff_2[i] = 1. - n2 / nEnt2; 
    if (h3!=NULL) {n3 += h3->GetBinContent(i); eff_3[i] = 1. - n3 / nEnt3;   }
    if (h4!=NULL) {n4 += h4->GetBinContent(i); eff_4[i] = 1. - n4 / nEnt4;   }
  }

  TCanvas* myCanvas = new TCanvas(myCan,myCan, 600,900);

  TGraph* gr12 = new TGraph(nBins1,eff_1,eff_2); 
  gr12->SetMarkerStyle(21); 
  gr12->SetMarkerColor(2); 
  gr12->GetXaxis()->SetTitle("Selection Efficiency");
  gr12->GetYaxis()->SetTitle("Backg Rejection"); 
 
  TGraph* gr13 = new TGraph(nBins1,eff_1,eff_3);  
  if (h3!=NULL) {
    gr13->SetMarkerStyle(23);   
    gr13->SetMarkerColor(4); 
  }

  TGraph* gr14 = new TGraph(nBins1,eff_1,eff_4);  
  if (h4!=NULL) {
    gr14->SetMarkerStyle(24);   
    gr14->SetMarkerColor(5); 
  }

  TMultiGraph* mg  = new TMultiGraph(); 
  mg->Add(gr12); 
  if (h3!=NULL) mg->Add(gr13); 
  if (h4!=NULL) mg->Add(gr14); 
  if(draw) mg->Draw("ALP*"); 

  mg->GetXaxis()->SetLimits(0.4,0.9);
  mg->SetMaximum(1.0);
  mg->SetMinimum(0.86);


  TLegend* legend = new TLegend(0.8,0.8,0.9,0.9);
  legend->AddEntry(gr12,"graph 1","P"); // "L" for line, "F" for fill color/style
  legend->AddEntry(gr13,"graph 2","P");
  legend->AddEntry(gr14,"graph 3","P");
  legend->Draw();

  TCanvas* c1 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(myCanvas); 
  TString name = myCan     + ".eps"; 
  c1->Print(name);

  delete[] eff_1;
  delete[] eff_2;
  delete[] eff_3;
  delete[] eff_4;

} 


void
Tools::RejVsEffComp(const TH1F* h1, const TH1F* h2, const TH1F* h3, const TH1F* h4,TString myCan="testCan", bool draw=false )
  // Comapre performance of 2 different NN 
{ 
  double n1 = 0;
  double n2 = 0;
  double n3 = 0;
  double n4 = 0;

  int nBins1 = h1->GetNbinsX(); 
  int nBins2 = h2->GetNbinsX(); 
  int nBins3 = h3->GetNbinsX(); 
  int nBins4 = h4->GetNbinsX(); 

  double nEnt1 = h1->GetEntries(); 
  double nEnt2 = h2->GetEntries(); 
  double nEnt3 = h3->GetEntries(); 
  double nEnt4 = h4->GetEntries(); 

  if (nBins1!=nBins2 ||nBins1!=nBins3 ||nBins1!=nBins4  ) {
    cout << " RejVsEffComp: Error: histos with different binnings " << endl; 
    return ; 
  }

  double *eff_1, *eff_2, *eff_3, *eff_4;
  eff_1 = new double[nBins1+2];
  eff_2 = new double[nBins1+2];
  eff_3 = new double[nBins1+2];
  eff_4 = new double[nBins1+2];


  for (int i = h1->GetNbinsX()+1; i>=0; i--) {
    n1 += h1->GetBinContent(i);
    n2 += h2->GetBinContent(i);
    eff_1[i] = n1 / nEnt1; 
    eff_2[i] = 1. - n2 / nEnt2; 

    n3 += h3->GetBinContent(i);
    n4 += h4->GetBinContent(i);
    eff_3[i] = n3 / nEnt3; 
    eff_4[i] = 1. - n4 / nEnt4; 
  }

  TCanvas* myCanvas = new TCanvas(myCan,myCan, 600,900);

  TGraph* gr12 = new TGraph(nBins1,eff_1,eff_2); 
  gr12->SetMarkerColor(2); 
  gr12->GetXaxis()->SetTitle("Selection Efficiency");
  gr12->GetYaxis()->SetTitle("Backg Rejection"); 
 
  TGraph* gr34 = new TGraph(nBins3,eff_3,eff_4);  
  gr34->SetMarkerColor(4); 

  TMultiGraph* mg  = new TMultiGraph(); 
  mg->Add(gr12); 
  mg->Add(gr34); 

  if(draw) mg->Draw("ALP*"); 

  mg->GetXaxis()->SetLimits(0.4,0.9);
  mg->SetMaximum(1.0);
  mg->SetMinimum(0.86);

  TLegend* legend = new TLegend(0.8,0.8,0.9,0.9);
  legend->AddEntry(gr12,"graph 1-2","P"); // "L" for line, "F" for fill color/style
  legend->AddEntry(gr34,"graph 3-4","P");
  legend->Draw();

  TCanvas* c1 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(myCanvas); 
  TString name = myCan     + ".eps"; 
  c1->Print(name);

  delete[] eff_1;
  delete[] eff_2;
  delete[] eff_3;
  delete[] eff_4;



} 


TGraph* 
Tools::RejVsEffUnit(const TH1F* h1, const TH1F* h2)
  // Return the graph of Eff Back vs Eff Signal 
{ 
  double n1 = 0;
  double n2 = 0;

  int nBins1 = h1->GetNbinsX(); 
  int nBins2 = h2->GetNbinsX(); 

  double nEnt1 = h1->GetEntries(); 
  double nEnt2 = h2->GetEntries(); 

  if (nBins1!=nBins2  ) {
    cout << " RejVsEffComp: Error: histos with different binnings " << endl; 
    return 0 ; 
  }

  double *eff_1, *eff_2; 
  eff_1 = new double[nBins1+2];
  eff_2 = new double[nBins1+2];


  for (int i = h1->GetNbinsX()+1; i>=0; i--) {
    n1 += h1->GetBinContent(i);
    n2 += h2->GetBinContent(i);
    eff_1[i] = n1 / nEnt1; 
    eff_2[i] = 1. - n2 / nEnt2; 
  }

  TGraph* gr12 = new TGraph(nBins1,eff_1,eff_2); 


  delete[] eff_1;
  delete[] eff_2;

  return gr12; 

} 


//DrawNHist: takes any number of histograms and colors and plots them on top of eachother
//Format for arguments: number of histograms (int), histogram 1(TH1F), histogram 2(TH1F), etc..., color for histogram 1(int), color for histogram 2(int), etc...
//Note that the number of histograms and colors must match, and they must be passed in the corresponding order. 
void 
Tools::DrawNHist(Int_t count, ...) 
{
  const Int_t count2=count;
  double max=0;
  double tempmax=0;
  Int_t colors[count2];
  va_list list2;
  va_start ( list2, count );
  //finds max. of the histograms or puts colors into colorArray
  for(Int_t q=0; q<2*count; q++)
    {
      if(q<count)
	{
	  tempmax = va_arg( list2, TH1F* )->GetMaximum();
	  max = TMath::Max(max , tempmax);
	}
      else
	{
 	  colors[q-count]=(int)va_arg(list2, int);
	}
    }
  for(Int_t i=0; i<count; i++)
    //    cout<<colors[i]<<endl;
  
  va_end ( list2 );
  //set hist colors, rescale Y axis, draw.
  va_start(list2, count);
  for(Int_t y=0; y<count; y++)
    {
      TH1F *tmpHist = va_arg( list2, TH1F* );
      tmpHist->SetLineColor(colors[y]);
      tmpHist->GetYaxis()->SetRangeUser(0, 1.1*max);
      if(y==0)
	tmpHist->Draw();
      else
	tmpHist->Draw("same");
    }
  va_end ( list2 );    
}

//DrawNHist2: Takes arrays of histograms and histogram colors and plots them. Basically same as DrawNHist, but histograms and colors must first be put into arrays. This is probably safter than DrawNHist.
//Arguments: array of histograms(TH1F), array of colors (int), size of arrays (int)
//Dimention of histogram and color arrays must match
void 
Tools::DrawNHist2(TH1F* histArray[], Int_t colorArray[], Int_t arraySize=1) 
{

  double max=0;
  double tempmax=0;
  //finds max of histograms
  for(Int_t i=0; i<arraySize; i++)
    {
      tempmax = histArray[i]->GetMaximum();
      max = TMath::Max(max , tempmax);
    }
  //set hist colors, rescale Y axis, draw.
  for(Int_t y=0; y<arraySize; y++)
    {
      TH1F *tmpHist = histArray[y];
      tmpHist->SetLineColor(colorArray[y]);
      tmpHist->GetYaxis()->SetRangeUser(0, 1.1*max);
      if(y==0)
	tmpHist->Draw();
      else
	tmpHist->Draw("same");
    }
}



//OptimizeCutN2: same as OptimizeCutN, except takes histograms and luminosities in arrays. Probably safter than OptimizeCutN.
//Arguments: array of histograms with signal as first entry (TH1F), array of luminosities with signal luminosity as first entry (double), signal color (int), background color (int), draw or no (bool), number of histograms or dimension of arrays (int).
//note that signal must be first in both arrays, and their dimensions must match. 
void 
Tools::OptimizeCutN2(TH1F* histoArray[], double lums[], Int_t sigColor=1, Int_t bkgColor=2, bool draw=false, Int_t numHists=1, TString myCanName="testCan")
{ 

  Int_t colors[2];
  colors[0]=sigColor;
  colors[1]=bkgColor;
  double refLumi = 300.; 

  //makes canvas
  TCanvas* myCan = new TCanvas(myCanName, myCanName, 600,900);; 
  myCan->Clear();
  myCan->Divide(1,2);

  TH1F* hsignalo = (TH1F*)histoArray[0];
  TH1F* hsignal = (TH1F*)hsignalo->Clone();
  hsignal->Sumw2();
  Int_t numBins=hsignalo->GetNbinsX();
  double axisLow=hsignalo->GetXaxis()->GetXmin();
  double axisHigh=hsignalo->GetXaxis()->GetXmax();
  TH1F* totalsigLeft = new TH1F("totalsigLeft", "Total Significance", numBins, axisLow, axisHigh);
  TH1F* totalsigRight = new TH1F("totalsigRight", "Total Significance", numBins, axisLow, axisHigh);
  TH1F* hbackgrounds = new TH1F("hbackgrounds", "Total Background", numBins, axisLow, axisHigh);


  for(Int_t k=1; k<numHists; k++)
    {
      TH1F *hbkgo = (TH1F*)histoArray[k]; 
      TH1F* hbkg = (TH1F*)hbkgo->Clone(); 
      hbkg->Sumw2();
      hbkg->Scale(refLumi/lums[k]);
      hbackgrounds->Add(hbkg);  

      double nsig, nbkg;
      nsig = 0;
      nbkg = 0;
      for (int i = 0; i <= hsignal->GetNbinsX()+1; i++) {
	nbkg += hbkgo->GetBinContent(i)*refLumi/lums[k];
	totalsigLeft->SetBinContent(i, totalsigLeft->GetBinContent(i)+nbkg);
	if(k==numHists-1) 
	  {
	    nsig += hsignalo->GetBinContent(i)*refLumi/lums[0];
	    double signif = 0;
	    if (nsig+totalsigLeft->GetBinContent(i) > 0) signif=nsig/sqrt(nsig+totalsigLeft->GetBinContent(i));
	    totalsigLeft->SetBinContent(i, signif);
	  }
      }
      nsig = 0;
      nbkg = 0;
      for (int i = hsignal->GetNbinsX()+1; i >= 0; i--) {
	nbkg += hbkgo->GetBinContent(i)*refLumi/lums[k];
	totalsigRight->SetBinContent(i, totalsigRight->GetBinContent(i)+nbkg);
	if(k==numHists-1) 
	  {
	    nsig += hsignalo->GetBinContent(i)*refLumi/lums[0];
	    double signif = 0;
	    if (nsig+totalsigRight->GetBinContent(i) > 0) signif=nsig/sqrt(nsig+totalsigRight->GetBinContent(i));
	    totalsigRight->SetBinContent(i, signif);
	  }
      }
    }

  hsignal->Scale(refLumi/lums[0]);
  
  TH1F* histPair[2];
  histPair[0]=hsignal;
  histPair[1]=hbackgrounds;
  
  TH1F* sigPair[2];
  sigPair[0]=totalsigLeft;
  sigPair[1]=totalsigRight;
  
  if(draw) 
    {
      myCan->cd(1);
      DrawNHist2(histPair, colors, 2); 
      myCan->cd(2);
      DrawNHist2(sigPair, colors, 2);	 
    }  
  
//   double bestS = (((double)(totalsigLeft->GetMaximumBin())/numBins))*(axisHigh-axisLow)+axisLow; 
//   double bestS2 = (((double)(totalsigRight->GetMaximumBin())/numBins))*(axisHigh-axisLow)+axisLow; 
//   cout << " Maximum sensitivity: integrating L->R: " << bestS << "; integrating R->L: " << bestS2<< endl; 

  PrintSignificance(totalsigLeft,totalsigRight );  

  TCanvas* c1 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(myCan);
  TString name = myCanName + ".eps"; 

  c1->Print(name);
}


// ------------------------ test -------------------------
//Adapted from OptimizeCutN: Takes histograms for signal and any number of backgrounds with luminosities, computes significances, and plots the distributions along with significances. Uses DrawNHist
//Arguments: plot or no (bool), total number of signal and background histograms(int), canvas name (TString), signal color (int), background color (int), signal lumiosity (double), background 1 luminosity (double), bkg 2 luminosity (double), etc..., signal histogram (TH1F), background histogram 1(TH1F), background histogram 2(TH1F), etc...
//Only takes one signal histogram. total number of histograms and luminosities must match. 
// ===============================
// NB: first histo is always data! 
// ===============================
void 
Tools::CompareData_MC(bool draw=false, const Int_t numHists=2, TString myCanName="testCan", Int_t sigColor=1, Int_t bkgColor=2, ... )
{ 
  double lums[numHists];

  //makes canvas
  TCanvas* myCan = new TCanvas(myCanName, myCanName, 600,900);; 
  myCan->Clear();
  //  myCan->Divide(1,1);

  va_list list;
  va_start(list, numHists);
  
  TH1F* hsignalo = new TH1F;
  TH1F* hsignal = new TH1F;
  
  TH1F* hbackgrounds = new TH1F;

  for(Int_t k=0; k<numHists-1; k++) {
    if (k==0) {
      // signal==DATA 
      hsignalo = (TH1F*)va_arg(list, TH1F*);
      hsignal = (TH1F*)hsignalo->Clone(); 
      hsignal->Sumw2();      
      lums[k]=(double)va_arg(list, double);
      cout << "Luminosity = " <<   lums[k] << endl; 
    }

    // if it's not signal, it must be background...  
    TH1F* hbkgo = (TH1F*)va_arg(list, TH1F*); 
    lums[k]=(double)va_arg(list, double);
    cout << "Luminosity = " <<   lums[k] << endl; 

    TH1F* hbkg = (TH1F*)hbkgo->Clone(); 
    hbkg->Sumw2();
    hbkg->Scale(lums[0]/lums[k]);
    hbackgrounds->Add(hbkg);
  }
  
  va_end(list);
  
  if(draw) 
    {
      //      myCan->cd(1);
      DrawNHist(2, hsignal, hbackgrounds, sigColor, bkgColor);       

      hsignal->Draw();  hbackgrounds->Draw("same"); 


    } 

  TCanvas* c1 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(myCan);
  TString name = myCanName + ".eps"; 
  c1->Print(name);
  
} 


// ===========================================================

//OptimizeCutN: Takes histograms for signal and any number of backgrounds with luminosities, computes significances, and plots the distributions along with significances. Uses DrawNHist
//Arguments: plot or no (bool), total number of signal and background histograms(int), canvas name (TString), signal color (int), background color (int), signal lumiosity (double), background 1 luminosity (double), bkg 2 luminosity (double), etc..., signal histogram (TH1F), background histogram 1(TH1F), background histogram 2(TH1F), etc...
//Only takes one signal histogram. total number of histograms and luminosities must match. 
void 
Tools::OptimizeCutN(bool draw=false, const Int_t numHists=2, TString myCanName="testCan", Int_t sigColor=1, Int_t bkgColor=2, ... )
{ 

  va_list list;

  double lums[numHists];
  double refLumi = 300.; 

  //makes canvas
  TCanvas* myCan = new TCanvas(myCanName, myCanName, 600,900);; 
  myCan->Clear();
  myCan->Divide(1,2);
  
  va_start(list, numHists);

  TH1F* hsignalo = new TH1F;
  TH1F* hsignal = new TH1F;
  
  TH1F* totalsigLeft=0;  //new TH1F;
  TH1F* totalsigRight=0; //new TH1F;

  TH1F* hbackgrounds = new TH1F;

  Int_t numBins=0;
  double axisLow=0;
  double axisHigh=0; 

  for(Int_t k=0; k<2*numHists; k++)
    {
      //gets luminosities
      if(k<numHists)
	{
	  lums[k]=(double)va_arg(list, double);
	}
      //gets signal histogram
      else if(k==numHists)
	{

	  hsignalo = (TH1F*)va_arg(list, TH1F*);

	  hsignal = (TH1F*)hsignalo->Clone(); 
	  hsignal->Sumw2();
	  numBins=hsignalo->GetNbinsX();
	  axisLow=hsignalo->GetXaxis()->GetXmin();
	  axisHigh=hsignalo->GetXaxis()->GetXmax();
	  totalsigLeft = new TH1F("totalsigLeft", "Total Significance", numBins, axisLow, axisHigh);
	  totalsigRight = new TH1F("totalsigRight", "Total Significance", numBins, axisLow, axisHigh);
	  hbackgrounds = new TH1F("hbackgrounds", "Total Background", numBins, axisLow, axisHigh);


	}
      //gets background histograms and computes significances, integrating from left and right
      else 
	{  
	  TH1F* hbkgo = (TH1F*)va_arg(list, TH1F*); 
	  TH1F* hbkg = (TH1F*)hbkgo->Clone(); 
	  hbkg->Sumw2();
	  hbkg->Scale(refLumi/lums[k-numHists]);
	  hbackgrounds->Add(hbkg);
	  
	  double nsig, nbkg;
	  nsig = 0;
	  nbkg = 0;
	  for (int i = 0; i <= hsignal->GetNbinsX()+1; i++) {
	    nbkg += hbkgo->GetBinContent(i)*refLumi/lums[k-numHists];
	    totalsigLeft->SetBinContent(i, totalsigLeft->GetBinContent(i)+nbkg);
	    if(k==2*numHists-1) 
	      {
		nsig += hsignalo->GetBinContent(i)*refLumi/lums[0];
		double signif = 0;
		if (nsig+totalsigLeft->GetBinContent(i) > 0) signif=nsig/sqrt(nsig+totalsigLeft->GetBinContent(i));
		totalsigLeft->SetBinContent(i, signif);
	      }
	  }
	  nsig = 0;
	  nbkg = 0;
	  for (int i = hsignal->GetNbinsX()+1; i >= 0; i--) {
	    nbkg += hbkgo->GetBinContent(i)*refLumi/lums[k-numHists];
	    totalsigRight->SetBinContent(i, totalsigRight->GetBinContent(i)+nbkg);
	    if(k==2*numHists-1) 
	      {
		nsig += hsignalo->GetBinContent(i)*refLumi/lums[0];
		double signif = 0;
		if (nsig+totalsigRight->GetBinContent(i) > 0) signif=nsig/sqrt(nsig+totalsigRight->GetBinContent(i));
		totalsigRight->SetBinContent(i, signif);
	      }
	  }
	  
	}
    }

  
  va_end(list);
  
  //draws plots
  hsignal->Scale(refLumi/lums[0]);
  
  if(draw) 
    {
      myCan->cd(1);
      DrawNHist(2, hsignal, hbackgrounds, sigColor, bkgColor); 
      myCan->cd(2);
      DrawNHist(2, totalsigLeft, totalsigRight, sigColor, bkgColor);
      
    } 
  //outputs values with maximum signifcance
  // Calculate the cut values: 

  PrintSignificance(totalsigLeft , totalsigRight);  

//   int bin1 = totalsigLeft->GetMaximumBin(); 
//   int bin2 = totalsigRight->GetMaximumBin();

//   double bestS = ((float(totalsigLeft->GetMaximumBin())/numBins))*(axisHigh-axisLow)+axisLow; 
//   double bestS2 = ((float(totalsigRight->GetMaximumBin())/numBins))*(axisHigh-axisLow)+axisLow; 

//   cout << " Optimal cut value: integrating L->R: " << bestS << "; integrating R->L: " << bestS2<< endl;      
//   cout << " Maximum sensitivity: integrating L->R: " << totalsigLeft->GetBinContent(bin1)  
//        << "; integrating R->L: " << totalsigRight->GetBinContent(bin2) << endl;      


  TCanvas* c1 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(myCan);
  TString name = myCanName + ".eps"; 

  c1->Print(name);
  
} 

void 
Tools::PrintSignificance(TH1F* hsig1 ,TH1F* hsig2 )  
{ 
  int bin1 = hsig1->GetMaximumBin(); 
  int bin2 = hsig2->GetMaximumBin();

  int numBins=hsig1->GetNbinsX();
  double axisLow=hsig1->GetXaxis()->GetXmin();
  double axisHigh=hsig1->GetXaxis()->GetXmax();
  
  double bestS = ((float(hsig1->GetMaximumBin())/numBins))*(axisHigh-axisLow)+axisLow; 
  double bestS2 = ((float(hsig2->GetMaximumBin())/numBins))*(axisHigh-axisLow)+axisLow; 

  // NB: these errors are not correct!!! One should just do error propagation based on the # of 
  // candidates (purely statistical) GS. 

//   double Sig2 = hsig2->GetBinContent(bin2); 
//   double Sig2p =0.; if(bin2>0) Sig2p= hsig2->GetBinContent(bin2-1); 
//   double Sig2m =0.; if(bin2< hsig2->GetNbinsX()) Sig2m= hsig2->GetBinContent(bin2+1); 

//   double error2 = 0.; 
//   if (abs(Sig2p-Sig2) > abs(Sig2-Sig2m)) 
//     error2= fabs(Sig2p-Sig2);
//   else
//     error2= fabs(Sig2m-Sig2);
//  error2 = error2/2.; 


//   double Sig1 = hsig1->GetBinContent(bin1); 
//   double Sig1p =0.; if(bin1>0) Sig1p= hsig1->GetBinContent(bin1-1); 
//   double Sig1m =0.; if(bin1< hsig1->GetNbinsX()) Sig1m= hsig1->GetBinContent(bin1+1); 

//   double error1 = 0.; 
//   if (abs(Sig1p-Sig1) > abs(Sig1-Sig1m)) 
//     error1= fabs(Sig1p-Sig1);
//   else
//     error1= fabs(Sig1m-Sig1);
//   error1 = error1/2.; 

  cout << " Optimal cut value: integrating L->R: " << bestS << "; integrating R->L: " << bestS2<< endl;      
  cout << " Maximum sensitivity: integrating L->R: " << hsig1->GetBinContent(bin1) // <<" +/- " << error1  
       << "; integrating R->L: " << hsig2->GetBinContent(bin2) // <<" +/- " << error2 
       << endl;      

}
