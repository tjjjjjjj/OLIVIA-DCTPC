// Mainframe macro generated from application: /afs/lns.mit.edu/public/sl54-32/bin/root.exe
// By ROOT version 5.32/00 on 2012-03-26 16:08:15



#ifndef ROOT_TGDockableFrame
#include "TGDockableFrame.h"
#endif
#ifndef ROOT_TGMdiDecorFrame
#include "TGMdiDecorFrame.h"
#endif
#ifndef ROOT_TG3DLine
#include "TG3DLine.h"
#endif
#ifndef ROOT_TGMdiFrame
#include "TGMdiFrame.h"
#endif
#ifndef ROOT_TGMdiMainFrame
#include "TGMdiMainFrame.h"
#endif
#ifndef ROOT_TGMdiMenu
#include "TGMdiMenu.h"
#endif
#ifndef ROOT_TGListBox
#include "TGListBox.h"
#endif
#ifndef ROOT_TGNumberEntry
#include "TGNumberEntry.h"
#endif
#ifndef ROOT_TGScrollBar
#include "TGScrollBar.h"
#endif
#ifndef ROOT_TGComboBox
#include "TGComboBox.h"
#endif
#ifndef ROOT_TGuiBldHintsEditor
#include "TGuiBldHintsEditor.h"
#endif
#ifndef ROOT_TGuiBldNameFrame
#include "TGuiBldNameFrame.h"
#endif
#ifndef ROOT_TGFrame
#include "TGFrame.h"
#endif
#ifndef ROOT_TGMenu
#include "TGMenu.h"
#endif
#ifndef ROOT_TGFileDialog
#include "TGFileDialog.h"
#endif
#ifndef ROOT_TGShutter
#include "TGShutter.h"
#endif
#ifndef ROOT_TGButtonGroup
#include "TGButtonGroup.h"
#endif
#ifndef ROOT_TGCanvas
#include "TGCanvas.h"
#endif
#ifndef ROOT_TGFSContainer
#include "TGFSContainer.h"
#endif
#ifndef ROOT_TGuiBldEditor
#include "TGuiBldEditor.h"
#endif
#ifndef ROOT_TGColorSelect
#include "TGColorSelect.h"
#endif
#ifndef ROOT_TGButton
#include "TGButton.h"
#endif
#ifndef ROOT_TRootContextMenu
#include "TRootContextMenu.h"
#endif
#ifndef ROOT_TGFSComboBox
#include "TGFSComboBox.h"
#endif
#ifndef ROOT_TGLabel
#include "TGLabel.h"
#endif
#ifndef ROOT_TGMsgBox
#include "TGMsgBox.h"
#endif
#ifndef ROOT_TRootGuiBuilder
#include "TRootGuiBuilder.h"
#endif
#ifndef ROOT_TGTab
#include "TGTab.h"
#endif
#ifndef ROOT_TGListView
#include "TGListView.h"
#endif
#ifndef ROOT_TGSplitter
#include "TGSplitter.h"
#endif
#ifndef ROOT_TGStatusBar
#include "TGStatusBar.h"
#endif
#ifndef ROOT_TGListTree
#include "TGListTree.h"
#endif
#ifndef ROOT_TGuiBldGeometryFrame
#include "TGuiBldGeometryFrame.h"
#endif
#ifndef ROOT_TGToolTip
#include "TGToolTip.h"
#endif
#ifndef ROOT_TGToolBar
#include "TGToolBar.h"
#endif
#ifndef ROOT_TRootEmbeddedCanvas
#include "TRootEmbeddedCanvas.h"
#endif
#ifndef ROOT_TCanvas
#include "TCanvas.h"
#endif
#ifndef ROOT_TGuiBldDragManager
#include "TGuiBldDragManager.h"
#endif

#include "Riostream.h"

#include "GuiFirFrame.hh"
#include "TApplication.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TMath.h"
#include <vector>
#include <iostream>
#include "kernels.cc"
using std::vector;
using std::cout;
using std::cerr;
using std::endl;


ClassImp(GuiFirFrame);
   static const int BLACKMAN=0x0;
   static const int BHARRIS=0x1;
   static const int BNUTALL=0x2;
   static const int BOHMANN=0x3;
   static const int COSINE=0x4;
   static const int EXPOCONV=0x5;
   static const int EXPO=0x6;
   static const int FERMICONV=0x7;
   static const int FLATTOP=0x8;
   static const int FUN=0x9;
   static const int GAUSCONV=0xA;
   static const int GAUS=0xB;
   static const int HAMMING=0xC;
   static const int HANN=0xD;
   static const int KAISER=0xE;
   static const int LANCZOS=0xF;
   static const int NUTALL=0x10;
   static const int RAISECOS=0x11;
   static const int RECT=0x12;
   static const int SG=0x13;
   static const int TRIANGLE=0x14;
   static const int TRICONV=0x15;
   static const int USER=0x16;

GuiFirFrame::GuiFirFrame(const TGWindow* win,unsigned int w, unsigned int h)

 : TGMainFrame(win,10,10,kMainFrame | kVerticalFrame)
{
  fNCoeff = 11;
  fNPar = 10;
  fNbin = 1000;
  SetGUI(w,h);
  gROOT->SetStyle("Pub");
  gStyle->SetHistLineColor(kGreen+2);
  gStyle->SetHistLineWidth(2);
  gStyle->SetTitleOffset(1,"y");
  gStyle->SetTitleOffset(1,"x");
  gStyle->SetLabelSize(0.03,"x");
  gStyle->SetLabelSize(0.03,"y");
  gStyle->SetTitleSize(0.04,"x");
  gStyle->SetTitleSize(0.04,"y");
  coeffH = new TH1D("coeffH" ,"coeffs;Coefficient;Value",fNCoeff,0,fNCoeff);  
  magH = new TH1D("magH","Magnitude;Frequency [Bin^{-1}];Magnitude",fNbin+1,0,0.5+0.5/(fNbin));
  phH = new TH1D("phH","Phase;Frequency [Bin^{-1}];Phase",fNbin+1,0,0.5+0.5/(fNbin));
  reH = new TH1D("reH","Real;Frequency [Bin^{-1}];Real",fNbin+1,0,0.5+0.5/(fNbin));
  imH = new TH1D("imH","Imaginary;Frequency [Bin^{-1}];Imaginary",fNbin+1,0,0.5+0.5/(fNbin));

  fParSetButton->Connect("Clicked()","GuiFirFrame",this,"SetParameter()");
  fCoeffSetButton->Connect("Clicked()","GuiFirFrame",this,"SetCoefficient()");
  fCoeffUpdateButton->Connect("Clicked()","GuiFirFrame",this,"UpdateCoefficients()");
  fCalcFreqButton->Connect("Clicked()","GuiFirFrame",this,"UpdateFrequency()");
//  fFilterLabel->Connect("Clicked()","GuiFirFrame",this,"SetFilter()");
  fFilterComboBox->Connect("Selected(Int_t)","GuiFirFrame",this,"SetFilter()");
  Connect("CloseWindow()","GuiFirFrame",this,"Terminate()");
  

  fFileMenu->Connect("Activated(Int_t)","GuiFirFrame",this,"HandleFileMenu(int)" );
  fViewMenu->Connect("Activated(Int_t)","GuiFirFrame",this,"HandleViewMenu(int)" );
  fOptionMenu->Connect("Activated(Int_t)","GuiFirFrame",this,"HandleOptionMenu(int)" );
  
  SetFilter();

  fFilterFun = 0;
}

void GuiFirFrame::Terminate()
{
  gApplication->Terminate(0);
}

GuiFirFrame::~GuiFirFrame()
{

  if (fFilterFun) delete fFilterFun;
  delete coeffH;
  delete magH;
  delete reH;
  delete imH;

  delete fMagECan; 
  delete fPhaseECan;
  delete fRealECan;
  delete fImagECan;
  delete fCoeffECan;
  delete fMagFrame;
  delete fPhaseFrame;
  delete fRealFrame;
  delete fImagFrame;
  delete fCoeffFrame;
  delete fHistoTab;
  delete fParEntry;
  delete fParNumberEntry;
  delete fParSetButton; 
  delete fCoeffSetButton;
  delete fCoeffNumberEntry;
  delete fCoeffEntry;
  delete fCoeffUpdateButton; 
  delete fFilterComboBox;
  delete fFilterLabel;
  delete fNumCoeffEntry;
  delete fCoeffsLabel;
  delete fFilterProps;
  delete fGuiCompFrame;

}

void GuiFirFrame::HandleFileMenu(int x)
{
  switch(x)
  {
    case 0:
      gApplication->Terminate(0);
      break;
    case 1: 
      break;
    case 2:
      break;
    case 3:
      break;
    case 4:
      break;
    default:
      break;

  }
}

void GuiFirFrame::HandleViewMenu(int x)
{
  switch(x)
  {
    case 0:
      fHistoTab->SetTab(0);
      break;
    case 1:
      fHistoTab->SetTab(1);
      break;
    case 2:
      fHistoTab->SetTab(2);
      break;
    case 3:
      fHistoTab->SetTab(3);
      break;
    case 4:
      fHistoTab->SetTab(4);
      break;
    default:
      break;
  }
}

void GuiFirFrame::HandleOptionMenu(int x)
{

  switch(x)
  {
    case 0:
      if (fOptionMenu->IsEntryChecked(x)){
        fOptionMenu->UnCheckEntry(x);
        fAutoUpdate = false;
      }else{
        fOptionMenu->CheckEntry(x);
        fAutoUpdate = true;
        UpdateCoefficients();
      }
      break;
    case 1:
      if (fOptionMenu->IsEntryChecked(x)){
        fOptionMenu->UnCheckEntry(x);
        fAutoFreq = false; 
      }else{
        fOptionMenu->CheckEntry(x);
        fAutoFreq = true;
        UpdateFrequency();
      }

      break;
    case 2:
      if (fOptionMenu->IsEntryChecked(x)){
        fOptionMenu->UnCheckEntry(x);
        fMagECan->GetCanvas()->SetLogz(false);
      }else{
        fOptionMenu->CheckEntry(x);
        fMagECan->GetCanvas()->SetLogz(true);
      }
      fMagECan->GetCanvas()->Update();
      break;
    default:
      break;
  }
}

void GuiFirFrame::SetFilter(bool update)
{
  
  fNCoeff = fNumCoeffEntry->GetIntNumber();
  if (fCoeffs.size()<(unsigned int)fNCoeff) fCoeffs.resize(fNCoeff);
  if (!update) fNCoeff = 11;
  int sel = fFilterComboBox->GetSelected();
  vector<double> vec;
  vector<double> par;
  TGString text;
  switch(sel)
  {
   case BLACKMAN:
     fNPar = 1;
     if (fPars.size()<(unsigned int)fNPar) fPars.resize(fNPar);
     if (!update) fPars[0] = 0.16;
     vec = blackmanWin((fNCoeff)/2,fPars[0]);
     text = "FILTER: \nBlackman\n";
     text.Append(TString::Format("LENGTH: %i\n",vec.size()));
     text.Append("PARAMETERS: 1\n");
     text.Append(TString::Format("0: %4.2f\n",fPars[0]));
     fFilterProps->SetText(text);
     break;
   case BHARRIS:
     fNPar = 0;
     vec = blackmanHarrisWin((fNCoeff)/2);
     text = "FILTER: \nBlackman-Harris\n";
     text.Append(TString::Format("LENGTH: %i\n",vec.size()));
     text.Append("PARAMETERS: 0\n");
     fFilterProps->SetText(text);
     break;
   case BNUTALL:
     fNPar = 0;
     vec = blackmanNutallWin((fNCoeff)/2);
     text = "FILTER: \nBlackman-Nutall\n";
     text.Append(TString::Format("LENGTH: %i\n",vec.size()));
     text.Append("PARAMETERS: 0\n");
     fFilterProps->SetText(text);
     break;
   case BOHMANN:
     fNPar = 0;
     vec = bohmanWin((fNCoeff)/2);
     text = "FILTER: \nBohman\n";
     text.Append(TString::Format("LENGTH: %i\n",vec.size()));
     text.Append("PARAMETERS: 0\n");
     fFilterProps->SetText(text);
     break;
   case COSINE:
     fNPar = 0;
     vec = cosineWin((fNCoeff)/2);
     text = "FILTER: \nCosine\n";
     text.Append(TString::Format("LENGTH: %i\n",vec.size()));
     text.Append("PARAMETERS: 0\n");
     fFilterProps->SetText(text);
     break;
   case EXPOCONV:
     fNPar = 2;
     if (fPars.size()<(unsigned int)fNPar) fPars.resize(fNPar);
     if (!update) {fPars[0] = 4; fPars[1] = 0.01;}
     vec = expoConv(fPars[0],fPars[1]);
     text = "FILTER: \nExponential Conv.\n";
     text.Append(TString::Format("LENGTH: %i\n",vec.size()));
     text.Append("PARAMETERS: 2\n");
     text.Append(TString::Format("0: Decay: %4.2f\n",fPars[0]));
     text.Append(TString::Format("1: Cutoff: %4.2f\n",fPars[1]));
     fFilterProps->SetText(text);
     break;
   case EXPO:
     fNPar = 1;
     if (fPars.size()<(unsigned int)fNPar) fPars.resize(fNPar);
     if (!update) {fPars[0] = 4;}
     vec = expoWin((fNCoeff)/2,fPars[0]);
     text = "FILTER: \nExponential\n";
     text.Append(TString::Format("LENGTH: %i\n",vec.size()));
     text.Append("PARAMETERS: 1\n");
     text.Append(TString::Format("0: Decay: %4.2f\n",fPars[0]));
     fFilterProps->SetText(text);
     break;
   case FERMICONV:
     fNPar = 3;
     if (fPars.size()<(unsigned int)fNPar) fPars.resize(fNPar);
     if (!update){fPars[0] = 3; fPars[1] = 5; fPars[2] = 0.01;}
     vec = fermiConv(fPars[0],fPars[1],fPars[2]);
     text = "FILTER: \nFermi Conv.\n";
     text.Append(TString::Format("LENGTH: %i\n",vec.size()));
     text.Append("PARAMETERS: 3\n");
     text.Append(TString::Format("0: x_{0}: %4.2f\n",fPars[0]));
     text.Append(TString::Format("1: T: %4.2f\n",fPars[1]));
     text.Append(TString::Format("2: Cutoff: %4.2f\n",fPars[2]));
     fFilterProps->SetText(text);
     break;
   case FLATTOP:
     fNPar = 0;
     vec = flatTopWin((fNCoeff)/2);
     text = "FILTER: \nFlat-Top\n";
     text.Append(TString::Format("LENGTH: %i\n",vec.size()));
     text.Append("PARAMETERS: 0\n");
     fFilterProps->SetText(text);
     break;
   case FUN:
     return;
     //Not done yet
     break;
   case GAUSCONV:
     fNPar = 2;
     if (fPars.size()<(unsigned int)fNPar) fPars.resize(fNPar);
     if (!update){fPars[0] = 2; fPars[1] = 0.01;}
     vec = gausConv(fPars[0],fPars[1]);
     text = "FILTER: \nGaussian Conv.\n";
     text.Append(TString::Format("LENGTH: %i\n",vec.size()));
     text.Append("PARAMETERS: 2\n");
     text.Append(TString::Format("0: Sigma: %4.2f\n",fPars[0]));
     text.Append(TString::Format("1: Cutoff: %4.2f\n",fPars[1]));
     fFilterProps->SetText(text);
     break;
   case GAUS:
     fNPar = 1;
     if (fPars.size()<(unsigned int)fNPar) fPars.resize(fNPar);
     if (!update){fPars[0] = 2;}
     vec = gausWin((fNCoeff)/2,fPars[0]);
     text = "FILTER: \nGaussian\n";
     text.Append(TString::Format("LENGTH: %i\n",vec.size()));
     text.Append("PARAMETERS: 1\n");
     text.Append(TString::Format("0: Sigma: %4.2f\n",fPars[0]));
     fFilterProps->SetText(text);
     break;
   case HAMMING:
     fNPar = 0;
     vec = hammingWin((fNCoeff)/2);
     text = "FILTER: \nHamming\n";
     text.Append(TString::Format("LENGTH: %i\n",vec.size()));
     text.Append("PARAMETERS: 0\n");
     fFilterProps->SetText(text);
     break;
   case HANN:
     fNPar = 0;
     vec = hannWin((fNCoeff)/2);
     text = "FILTER: \nHann\n";
     text.Append(TString::Format("LENGTH: %i\n",vec.size()));
     text.Append("PARAMETERS: 0\n");
     fFilterProps->SetText(text);
     break;
   case KAISER:
     fNPar = 1;
     if (fPars.size()<(unsigned int)fNPar) fPars.resize(fNPar);
     if (!update){fPars[0] = 3;}
     vec = kaiserWin((fNCoeff)/2,fPars[0]);
     break;
   case LANCZOS:
     fNPar = 0;
     vec = lanczosWin((fNCoeff)/2);
     text = "FILTER: \nLanczos\n";
     text.Append(TString::Format("LENGTH: %i\n",vec.size()));
     text.Append("PARAMETERS: 0\n");
     fFilterProps->SetText(text);
     break;
   case NUTALL:
     fNPar = 0;
     vec = nutallWin((fNCoeff)/2);
     text = "FILTER: \nNutall\n";
     text.Append(TString::Format("LENGTH: %i\n",vec.size()));
     text.Append("PARAMETERS: 0\n");
     fFilterProps->SetText(text);
     break;
   case RAISECOS:
     fNPar = 5;
     if (fPars.size()<(unsigned int)fNPar) fPars.resize(fNPar);
     if (!update){fPars[0]  = 0.5;
     fPars[1] = 0.5;
     fPars[2] = 0; fPars[3] = 0; fPars[4] = 0;}
     vec = raisedCosWin((fNCoeff)/2,fPars[0],fPars[1],fPars[2],fPars[3],fPars[4]);
     text = "FILTER: \nRaised Cosine\n";
     text.Append(TString::Format("LENGTH: %i\n",vec.size()));
     text.Append("PARAMETERS: 5\n");
     text.Append(TString::Format("0: a0: %4.2f\n",fPars[0]));
     text.Append(TString::Format("1: a1: %4.2f\n",fPars[1]));
     text.Append(TString::Format("2: a2: %4.2f\n",fPars[2]));
     text.Append(TString::Format("3: a3: %4.2f\n",fPars[3]));
     text.Append(TString::Format("4: a4: %4.2f\n",fPars[4]));
     fFilterProps->SetText(text);
     break;
   case RECT:
     fNPar = 0;
     vec = rectWin((fNCoeff)/2);
     text = "FILTER: \nRectangle\n";
     text.Append(TString::Format("LENGTH: %i\n",vec.size()));
     text.Append("PARAMETERS: 0\n");
     fFilterProps->SetText(text);
     break;
   case SG:
     fNPar = 2;
     if (fPars.size()<(unsigned int)fNPar) fPars.resize(fNPar);
     if (!update){fPars[0] = 4; fPars[1] = 0;   } 
     vec = savitzkyGolay((int)fPars[0],(fNCoeff)/2,(int)fPars[1]);
     text = "FILTER: \nSavitzky-Golay\n";
     text.Append(TString::Format("LENGTH: %i\n",vec.size()));
     text.Append("PARAMETERS: 2\n");
     text.Append(TString::Format("0: Order: %i\n",(int)fPars[0]));
     text.Append(TString::Format("1: Deriv.: %i\n",(int)fPars[1]));
     fFilterProps->SetText(text);
     
     break;
   case TRIANGLE:
     fNPar = 0;
     vec = triangleWin((fNCoeff)/2);
     text = "FILTER: \nTriangle\n";
     text.Append(TString::Format("LENGTH: %i\n",vec.size()));
     text.Append("PARAMETERS: 0\n");
     fFilterProps->SetText(text);
     break;
   case TRICONV:
     fNPar = 0;
     vec = triangleConv((fNCoeff)/2);
     text = "FILTER: \nTriangle Conv.\n";
     text.Append(TString::Format("LENGTH: %i\n",vec.size()));
     text.Append("PARAMETERS: 0\n");
     fFilterProps->SetText(text);
     break;
   case USER:
     fNPar = 0; 
     text = "FILTER: \nUser-Defined\n";
     text.Append(TString::Format("LENGTH: %i\n",fNCoeff));
     text.Append("PARAMETERS: 0\n");
     fFilterProps->SetText(text);
     return;
     break;
   default:
     break;
  }

  if (sel!=USER){ fNCoeff = vec.size();
    fNumCoeffEntry->SetIntNumber(fNCoeff);
    if (fCoeffs.size()<vec.size()) fCoeffs.resize(vec.size());
    for (unsigned int i = 0; i < vec.size(); i++)
      fCoeffs[i] = vec[i];
  }
 
  if (sel==GAUSCONV||sel==FERMICONV||sel==EXPOCONV)
  {
    fNumCoeffEntry->SetState(false);
  }else fNumCoeffEntry->SetState(true);

  SetAutoUpdate();
  if (fAutoUpdate){
    UpdateCoefficients(false);
  }
  if (fAutoFreq)
    UpdateFrequency();
}

void GuiFirFrame::UpdateFilter()
{
  SetFilter(true);
}

void GuiFirFrame::UpdateCoefficients()
{UpdateCoefficients(true);}

void GuiFirFrame::UpdateCoefficients(bool updateFil)
{
 // cout <<"Updating coefficients"<<endl;
  if (updateFil) UpdateFilter();
  coeffH->SetBins(fNCoeff,0,fNCoeff);
  for (int i = 0; i < fNCoeff; i++)
  {
    coeffH->SetBinContent(i+1,fCoeffs[i]);
  }
 // cout <<"Drawing "<<endl; 
  fCoeffECan->GetCanvas()->cd();
  coeffH->Draw();
  fCoeffECan->GetCanvas()->Update();
}

void GuiFirFrame::UpdateFrequency()
{
  //cout <<"Updating frequency stuff"<<endl;
  magH->SetBins(fNbin+1,0,0.5*(1+1./fNbin));
  reH->SetBins(fNbin+1,0,0.5*(1+1./fNbin));
  imH->SetBins(fNbin+1,0,0.5*(1+1./fNbin));
  phH->SetBins(fNbin+1,0,0.5*(1+1./fNbin));
  SetLogPlot();

  for (int i = 1; i<=fNbin+1; i++)
  {
    double re=0,im=0,f = magH->GetXaxis()->GetBinLowEdge(i)*TMath::TwoPi();;
    for (int j = 0; j<fNCoeff; j++)
    {
      re += fCoeffs[fNCoeff-j-1]*cos( -(fNCoeff-j-1)*f );
      im += fCoeffs[fNCoeff-j-1]*sin( -(fNCoeff-j-1)*f );
    }
    reH->SetBinContent(i,re);
    imH->SetBinContent(i,im);
    magH->SetBinContent(i,sqrt(re*re+im*im));
    phH->SetBinContent(i,atan2(im,re));
  }
  fMagECan->GetCanvas()->cd();
  magH->Draw();
  fMagECan->GetCanvas()->Update();
  fPhaseECan->GetCanvas()->cd();
  phH->Draw();
  fPhaseECan->GetCanvas()->Update();
  fRealECan->GetCanvas()->cd();
  reH->Draw();
  fRealECan->GetCanvas()->Update();
  fImagECan->GetCanvas()->cd();
  imH->Draw();
  fImagECan->GetCanvas()->Update();
  fCoeffECan->GetCanvas()->cd();

}

void GuiFirFrame::SetLogPlot()
{
  if (fOptionMenu->IsEntryChecked(2))
  {
    fMagECan->GetCanvas()->SetLogy(true);
  }else fMagECan->GetCanvas()->SetLogy(false);
}

void GuiFirFrame::SetAutoUpdate()
{
  if (fOptionMenu->IsEntryChecked(0))
  {
    fAutoUpdate = true;
  }else fAutoUpdate = false;
}

void GuiFirFrame::SetAutoFreq()
{
  if (fOptionMenu->IsEntryChecked(1))
  {
    fAutoFreq = true;
  }else fAutoFreq = false;
}


void GuiFirFrame::SetCoefficient()
{
  SetAutoUpdate();
  SetAutoFreq();
  int ncoeff = fCoeffNumberEntry->GetIntNumber();
  double coeff= fCoeffEntry->GetNumber();

  if (ncoeff>=(int)fCoeffs.size()) fCoeffs.resize(ncoeff+1);
  fCoeffs[ncoeff] = coeff;

  if (fAutoUpdate)
    UpdateCoefficients();

  if (fAutoFreq)
    UpdateFrequency();
}

void GuiFirFrame::SetParameter()
{
  SetAutoUpdate();
  int npar = fParNumberEntry->GetIntNumber();
  double param = fParEntry->GetNumber();
  
  if (npar>=(int)fPars.size()) fPars.resize(npar+1);
  fPars[npar] = param;
 
  if (fAutoUpdate)
    UpdateCoefficients();
  if (fAutoFreq)
    UpdateFrequency();

}

void GuiFirFrame::SetGUI(unsigned int w,unsigned int h)
{

   // main frame
  
   SetName("FIR Filter Analyzer");

   //
   fMenuBar = new TGMenuBar(this,400,20, kHorizontalFrame);

   fFileMenu = new TGPopupMenu(gClient->GetDefaultRoot());
/*
   fFileMenu->AddEntry("&Save Histograms...",1);
   fFileMenu->AddEntry("&Read Ascii File...",2);
   fFileMenu->AddEntry("&About...",3);
   fFileMenu->AddEntry("&Help...",4);
   fFileMenu->AddSeparator();
*/
   fFileMenu->AddEntry("&Quit",0);   
   fViewMenu = new TGPopupMenu(gClient->GetDefaultRoot());
   fViewMenu->AddEntry("&Coefficients...",0);
   fViewMenu->AddEntry("&Magnitude...",1);
   fViewMenu->AddEntry("&Phase...",2);
   fViewMenu->AddEntry("&Real Part...",3);
   fViewMenu->AddEntry("&Imaginary Part...",4);
   fOptionMenu = new TGPopupMenu(gClient->GetDefaultRoot());
   fOptionMenu->AddEntry("Auto Update &Coeffs",0);
   fOptionMenu->AddEntry("Auto Update &Frequency Plots",1);
   fOptionMenu->AddEntry("Magnitude in &Log Scale",2);

   fMenuBar->AddPopup("&File",fFileMenu,new TGLayoutHints(kLHintsTop | kLHintsLeft,0,4,0));
   fMenuBar->AddPopup("&View",fViewMenu,new TGLayoutHints(kLHintsTop | kLHintsLeft,0,4,0));
   fMenuBar->AddPopup("&Option",fOptionMenu,new TGLayoutHints(kLHintsTop | kLHintsLeft,0,4,0));
   AddFrame(fMenuBar);

   // composite frame
   fGuiCompFrame = new TGHorizontalFrame(this,546,435);//kVerticalFrame);
   fGuiCompFrame->SetName("fGuiCompFrame");
   fGuiCompFrame->SetLayoutBroken(kFALSE);

   fLeftFrame = new TGVerticalFrame(fGuiCompFrame,140,434);
   fRightFrame = new TGVerticalFrame(fGuiCompFrame,412,435);
   fTabFrame = new TGHorizontalFrame(fRightFrame,412,384);
   fBottomButtons = new TGHorizontalFrame(fRightFrame,412,54);   


   // tab widget
   fHistoTab = new TGTab(fTabFrame,432,384);

   // container of "Coefficients"
   fCoeffFrame = fHistoTab->AddTab("Coefficients");
   fCoeffFrame->SetLayoutManager(new TGVerticalLayout(fCoeffFrame));

   // embedded canvas
   fCoeffECan = new TRootEmbeddedCanvas(0,fCoeffFrame,424,355);
   Int_t wfCoeffECan = fCoeffECan->GetCanvasWindowId();
   TCanvas *ccoeff = new TCanvas("ccoeff", 9, 9, wfCoeffECan);
   fCoeffECan->AdoptCanvas(ccoeff);
   fCoeffFrame->AddFrame(fCoeffECan, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,2,2,2,2));


   // container of "Magnitude"
   fMagFrame = fHistoTab->AddTab("Magnitude");
   fMagFrame->SetLayoutManager(new TGVerticalLayout(fMagFrame));

   // embedded canvas
   fMagECan = new TRootEmbeddedCanvas(0,fMagFrame,424,355);
   Int_t wfMagECan = fMagECan->GetCanvasWindowId();
   TCanvas *cmag = new TCanvas("cmag", 9, 9, wfMagECan);
   fMagECan->AdoptCanvas(cmag);
   fMagFrame->AddFrame(fMagECan, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,2,2,2,2));


   // container of "Phase"
   fPhaseFrame = fHistoTab->AddTab("Phase");
   fPhaseFrame->SetLayoutManager(new TGVerticalLayout(fPhaseFrame));

   // embedded canvas
   fPhaseECan = new TRootEmbeddedCanvas(0,fPhaseFrame,424,355);
   Int_t wfPhaseECan = fPhaseECan->GetCanvasWindowId();
   TCanvas *cphase = new TCanvas("cphase", 9, 9, wfPhaseECan);
   fPhaseECan->AdoptCanvas(cphase);
   fPhaseFrame->AddFrame(fPhaseECan, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,2,2,2,2));


   // container of "Real"
   fRealFrame = fHistoTab->AddTab("Real");
   fRealFrame->SetLayoutManager(new TGVerticalLayout(fRealFrame));

   // embedded canvas
   fRealECan = new TRootEmbeddedCanvas(0,fRealFrame,424,355);
   Int_t wfRealECan = fRealECan->GetCanvasWindowId();
   TCanvas *creal = new TCanvas("creal", 9, 9, wfRealECan);
   fRealECan->AdoptCanvas(creal);
   fRealFrame->AddFrame(fRealECan, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,2,2,2,2));


   // container of "Imaginary"
   fImagFrame = fHistoTab->AddTab("Imaginary");
   fImagFrame->SetLayoutManager(new TGVerticalLayout(fImagFrame));

   // embedded canvas
   fImagECan = new TRootEmbeddedCanvas(0,fImagFrame,424,355);
   Int_t wfImagECan = fImagECan->GetCanvasWindowId();
   TCanvas *cimag = new TCanvas("cimag", 9, 9, wfImagECan);
   fImagECan->AdoptCanvas(cimag);
   fImagFrame->AddFrame(fImagECan, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,2,2,2,2));


   fHistoTab->SetTab(0);

//Histo bin label
   fHistoBinLabel = new TGLabel(fBottomButtons,"Histogram Bins");
   fHistoBinLabel->SetTextJustify(36);
   fHistoBinLabel->SetMargins(0,0,0,0);
   fHistoBinLabel->SetWrapLength(-1);
   fBottomButtons->AddFrame(fHistoBinLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY,2,2,2,2));
   //fHistoBinLabel->MoveResize(171,412,96,18);

//Hist bin entry
   fHistoBinEntry = new TGNumberEntry(fBottomButtons, (Double_t) fNbin,9,-1,(TGNumberFormat::EStyle) 0,(TGNumberFormat::EAttribute) 2);
   fHistoBinEntry->SetName("fHistoBinEntry");
   fBottomButtons->AddFrame(fHistoBinEntry, new TGLayoutHints(kLHintsLeft | kLHintsCenterY));
   //fHistoBinEntry->MoveResize(272,410,80,22);

//Frequency response button
   fCalcFreqButton = new TGTextButton(fBottomButtons,"Calculate Freqency Response",-1,TGTextButton::GetDefaultGC()(),TGTextButton::GetDefaultFontStruct(),kRaisedFrame);
   fCalcFreqButton->SetTextJustify(36);
   fCalcFreqButton->SetMargins(0,0,0,0);
   fCalcFreqButton->SetWrapLength(-1);
   fCalcFreqButton->Resize(184,24);
   fBottomButtons->AddFrame(fCalcFreqButton, new TGLayoutHints(kLHintsRight | kLHintsBottom,2,2,2,2));
   //fCalcFreqButton->MoveResize(360,410,184,24);
//Filter label
   fFilterLabel = new TGLabel(fLeftFrame,"Set Filter");
   fFilterLabel->SetTextJustify(36);
   fFilterLabel->SetMargins(0,0,0,0);
   fFilterLabel->SetWrapLength(-1);
   fFilterLabel->Resize(112,24);
   fLeftFrame->AddFrame(fFilterLabel, new TGLayoutHints(kLHintsExpandX | kLHintsTop,2,2,2,2));
  //fFilterLabel->MoveResize(0,53,112,24);

   // combo box
   fFilterComboBox = new TGComboBox(fLeftFrame,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
   fFilterComboBox->SetName("fFilterComboBox");
   fFilterComboBox->AddEntry("Blackman",BLACKMAN);
   fFilterComboBox->AddEntry("Blackman-Harris",BHARRIS);
   fFilterComboBox->AddEntry("Blackman-Nutall",BNUTALL);
   fFilterComboBox->AddEntry("Bohman",BOHMANN);
   fFilterComboBox->AddEntry("Cosine",COSINE);
   fFilterComboBox->AddEntry("Expo. Conv.",EXPOCONV);
   fFilterComboBox->AddEntry("Exponential",EXPO);
   fFilterComboBox->AddEntry("Fermi Conv.",FERMICONV);
   fFilterComboBox->AddEntry("Flat Top",FLATTOP);
   fFilterComboBox->AddEntry("Function",FUN);
   fFilterComboBox->AddEntry("Gauss. Conv.",GAUSCONV);
   fFilterComboBox->AddEntry("Gaussian",GAUS);
   fFilterComboBox->AddEntry("Hamming",HAMMING);
   fFilterComboBox->AddEntry("Hann",HANN);
   fFilterComboBox->AddEntry("Kaiser",KAISER);
   fFilterComboBox->AddEntry("Lanczos",LANCZOS);
   fFilterComboBox->AddEntry("Nutall",NUTALL);
   fFilterComboBox->AddEntry("Raised Cosine",RAISECOS);
   fFilterComboBox->AddEntry("Rectangle",RECT);
   fFilterComboBox->AddEntry("Savitzky-Golay",SG);
   fFilterComboBox->AddEntry("Triangle",TRIANGLE);
   fFilterComboBox->AddEntry("Triangle Conv.",TRICONV);
   fFilterComboBox->AddEntry("User",USER);
   fFilterComboBox->Resize(112,22);
   fFilterComboBox->Select(GAUS);
   fLeftFrame->AddFrame(fFilterComboBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fFilterComboBox->Resize(112,22);
   //fFilterComboBox->MoveResize(0,80,112,22);
   fCoeffsFrame = new TGHorizontalFrame(fLeftFrame,106,24);
//Coeffs label
   fCoeffsLabel = new TGLabel(fCoeffsFrame,"Coeffs");
   fCoeffsLabel->SetTextJustify(36);
   fCoeffsLabel->SetMargins(0,0,0,0);
   fCoeffsLabel->SetWrapLength(-1);
   fCoeffsFrame->AddFrame(fCoeffsLabel, new TGLayoutHints(kLHintsLeft | kLHintsExpandX,2,2,2,2));
   fCoeffsLabel->Resize(48,24);
   //fCoeffsLabel->MoveResize(0,106,48,24);
   

//Coefficients
   fNumCoeffEntry = new TGNumberEntry(fCoeffsFrame, (Double_t) fNCoeff,6,-1,(TGNumberFormat::EStyle) 0,(TGNumberFormat::EAttribute) 2);
   fNumCoeffEntry->SetName("fNumCoeffEntry");
   fNumCoeffEntry->Resize(59,22);
   fCoeffsFrame->AddFrame(fNumCoeffEntry, new TGLayoutHints(kLHintsExpandX | kLHintsRight));
   //fNumCoeffEntry->MoveResize(52,108,59,22);
   fLeftFrame->AddFrame(fCoeffsFrame,new TGLayoutHints(kLHintsExpandX | kLHintsLeft));

//Parameter set button
   fParSetButton = new TGTextButton(fLeftFrame,"Set Parameter",-1,TGTextButton::GetDefaultGC()(),TGTextButton::GetDefaultFontStruct(),kRaisedFrame);
   fParSetButton->SetTextJustify(36);
   fParSetButton->SetMargins(0,0,0,0);
   fParSetButton->SetWrapLength(-1);
   fParSetButton->Resize(112,24);
   fLeftFrame->AddFrame(fParSetButton, new TGLayoutHints(kLHintsExpandX | kLHintsTop,2,2,2,2));
  // fParSetButton->MoveResize(0,130,112,24);
   fParFrame = new TGHorizontalFrame(fLeftFrame,106,24);
//Parameter number entry
   fParNumberEntry = new TGNumberEntry(fParFrame, (Double_t) 0,3,-1,(TGNumberFormat::EStyle) 0,(TGNumberFormat::EAttribute) 1);
   fParNumberEntry->SetName("fParNumberEntry");
   fParFrame->AddFrame(fParNumberEntry, new TGLayoutHints(kLHintsExpandX | kLHintsLeft));
   fParNumberEntry->Resize(40,22);
   //fParNumberEntry->MoveResize(0,157,40,22);
//Parameter entry
   fParEntry = new TGNumberEntry(fParFrame, (Double_t) 0,8,-1,(TGNumberFormat::EStyle) 5);
   fParEntry->SetName("fParEntry");
   fParFrame->AddFrame(fParEntry, new TGLayoutHints(kLHintsRight | kLHintsTop));
   fParEntry->Resize(65,22);
   //fParEntry->MoveResize(40,157,72,22);
   fLeftFrame->AddFrame(fParFrame,new TGLayoutHints(kLHintsExpandX | kLHintsLeft,2,2,2,2));
//Coefficient set button
   fCoeffSetButton = new TGTextButton(fLeftFrame,"Set Coefficient",-1,TGTextButton::GetDefaultGC()(),TGTextButton::GetDefaultFontStruct(),kRaisedFrame);
   fCoeffSetButton->SetTextJustify(36);
   fCoeffSetButton->SetMargins(0,0,0,0);
   fCoeffSetButton->SetWrapLength(-1);
   fCoeffSetButton->Resize(112,24);
   fLeftFrame->AddFrame(fCoeffSetButton, new TGLayoutHints(kLHintsExpandX | kLHintsTop,2,2,2,2));
   //fCoeffSetButton->MoveResize(0,180,112,24);
   fCoeffsFrame2 = new TGHorizontalFrame(fLeftFrame,106,24);
//Coefficient number entry
   fCoeffNumberEntry = new TGNumberEntry(fCoeffsFrame2, (Double_t) 0,3,-1,(TGNumberFormat::EStyle) 0,(TGNumberFormat::EAttribute) 1);
   fCoeffNumberEntry->SetName("fCoeffNumberEntry");
   fCoeffsFrame2->AddFrame(fCoeffNumberEntry, new TGLayoutHints(kLHintsLeft | kLHintsExpandX,2,2,2,2));
   fCoeffNumberEntry->Resize(40,22);  
 //fCoeffNumberEntry->MoveResize(0,205,40,22);

//Coefficient entry
   fCoeffEntry = new TGNumberEntry(fCoeffsFrame2, (Double_t) 0,8,-1,(TGNumberFormat::EStyle) 5);
   fCoeffEntry->SetName("fCoeffEntry");
   fCoeffsFrame2->AddFrame(fCoeffEntry, new TGLayoutHints(kLHintsLeft | kLHintsExpandX,2,2,2,2));
   fCoeffEntry->Resize(65,22);
   //fCoeffEntry->MoveResize(40,205,72,22);
   fLeftFrame->AddFrame(fCoeffsFrame2,new TGLayoutHints(kLHintsLeft | kLHintsLeft));

//Coefficient update 
   fCoeffUpdateButton = new TGTextButton(fLeftFrame,"Update Coefficients",-1,TGTextButton::GetDefaultGC()(),TGTextButton::GetDefaultFontStruct(),kRaisedFrame);
   fCoeffUpdateButton->SetTextJustify(36);
   fCoeffUpdateButton->SetMargins(0,0,0,0);
   fCoeffUpdateButton->SetWrapLength(-1);
   fCoeffUpdateButton->Resize(112,24);
   fLeftFrame->AddFrame(fCoeffUpdateButton, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   //fCoeffUpdateButton->MoveResize(0,248,112,24);

   ULong_t ucolor;        // will reflect user color changes
   gClient->GetColorByName("#ffffff",ucolor);


   gClient->GetColorByName("#ffffff",ucolor);
   //
   //Filter properties
   fFilterProps = new TGLabel(fLeftFrame,"fFilterProps",TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),kSunkenFrame,ucolor);
   fFilterProps->SetTextJustify(36);
   fFilterProps->SetMargins(0,0,0,0);
   fFilterProps->SetWrapLength(-1);
   fFilterProps->Resize(112,134);
   fLeftFrame->AddFrame(fFilterProps, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,2,2,2,2));
   //fFilterProps->MoveResize(0,273,112,134);

   fGuiCompFrame->AddFrame(fLeftFrame,new TGLayoutHints(kLHintsExpandY,2,2,2,2));
   fTabFrame->AddFrame(fHistoTab,new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,2,2,2,2));
   fRightFrame->AddFrame(fTabFrame, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,2,2,2,2));
   fRightFrame->AddFrame(fBottomButtons,new TGLayoutHints(kLHintsExpandX |kLHintsBottom,2,2,2,2));
   fGuiCompFrame->AddFrame(fRightFrame,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,2,2,2,2));  
   AddFrame(fGuiCompFrame, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
//   fGuiCompFrame->Resize(GetDefaultSize());

   fGuiCompFrame->MoveResize(0,0,546,435); 

 //  Layout();
   
/*
   SetMWMHints(kMWMDecorAll,
                        kMWMFuncAll,
                        kMWMInputModeless);
*/
   MapSubwindows();

   Resize(GetDefaultSize());
   MapWindow();
   Resize(w,h);
}  
