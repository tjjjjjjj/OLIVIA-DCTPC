#include "../../MaxCam/DmtpcSkimDataset.hh"
#include "../../MaxCam/MaxCamPulse.hh"
#include <list>
#include <algorithm>
using std::pair;
using std::make_pair;
using std::list;
using std::copy;
#include "TCanvas.h"
#include "TArrow.h"
#include "TTree.h"

DmtpcSkimDataset *sd;
int ii;
int jj;

TCanvas* c1;
TCanvas* c2;

TTree* analysistree;
Double_t qanode, vpmt, qpmt;
Bool_t _shouldDraw;

void init(TString filename="/net/zwicky/dmtpc/production/skimout_4sh/skim/dmtpc_4sh_00147skim.root",Bool_t shouldDraw=false);
void next(void);
void print(void);
void draw(void);
void drawscope(void);
void drawpeakfinding(char* channel="B1_CHB");
void doanalysis(void);
void loop(Int_t n=0);

void init(TString filename,Bool_t shouldDraw){
  
  sd=new DmtpcSkimDataset();
  sd->openRootFile(filename);
  cout << "filename=" << filename << endl;
  ii=1;
  jj=0;
  cout << "ii=" << ii << endl;
  cout << "jj=" << jj << endl;
  sd->getEvent(ii);

  analysistree=new TTree("analysis","analysis tree for 4-shooter pmts");
  analysistree->Branch("qanode",&qanode,"qanode/D");
  analysistree->Branch("qpmt",&qpmt,"qpmt/D");
  analysistree->Branch("vpmt",&vpmt,"vpmt/D");

  _shouldDraw=shouldDraw;
  if(_shouldDraw) c1=new TCanvas("c1","c1",0,0,1000,800);
  if(_shouldDraw) c2=new TCanvas("c2","c2",0,0,1000,800);
  if(_shouldDraw) c1->Divide(2,1);

}

void next(void){
  //  cout << "sd->event()->ntriggers()=" << sd->event()->ntriggers() << endl;
  // have to do, otherwise we pick up datafiles
  // that have no trigger groups; I didn't 
  // do the algebra of this loop very well...
  do{
    if(jj+1<(sd->event()->ntriggers())) ++jj;
    else {
      jj=0;
      ++ii;
    }
    sd->getEvent(ii);
  } while (sd->event()->ntriggers()==0);
  draw();
}

void print(void){
  cout << "ii=" << ii << " jj=" << jj << endl;
  c1->Update();
}

void draw(void){
  cout << "*** Event: " << ii << " Trigger: " << jj << endl;
  
  if(_shouldDraw) drawscope();
  if(_shouldDraw) drawpeakfinding("B1_CHB");

  doanalysis();
}

void loop(Int_t n){
  for(Int_t ii=0; ii<n; ++ii) next();
}

void drawpeakfinding(char* channel){
  //  cout << "sd->event()->trigger_group(jj)=" << sd->event()->trigger_group(jj) << endl;
  MaxCamWaveformTools b1chb_wf(sd->event()->trigger_group(jj)->getWaveform("B1_CHB")->getWaveform(),false,10,1);
  //  b1chb_wf.print();
  
  //  cout << "b1chb_wf.getSizeOfPulseList()=" << b1chb_wf.getSizeOfPulseList() << endl;
  list<MaxCamPulse> L;
  b1chb_wf.getCopyOfPulseList(L);
  //  cout << "L.size()=" << L.size() << endl;
  
  if(L.size()>0){
    TH1F *wfi=(TH1F*)sd->event()->trigger_group(jj)->getWaveform(channel)->getWaveform()->Clone();
    TH1F *wfi_filled_pulse=(TH1F*)sd->event()->trigger_group(jj)->getWaveform(channel)->getWaveform()->Clone();
    wfi_filled_pulse->SetLineStyle(2);
    wfi_filled_pulse->SetLineColor(kGreen);
    wfi_filled_pulse->SetFillColor(kGreen);
    
    Int_t npeak=1;
    for( list< MaxCamPulse >::iterator iter=L.begin(); iter!=L.end(); ++iter){
      cout << "(*iter).getPulseHeight()=" 
	   << (*iter).getPulseHeight() 
	   << endl;
      
      c2->cd();
      Int_t halfbinrange=100;
      cout << "*** Event: " << ii << " Trigger: " << jj << " Peak: " << npeak << " (of " << L.size() << ")" << endl;
      wfi->GetXaxis()->SetRange((*iter).getBin()-halfbinrange,(*iter).getBin()+halfbinrange);
      wfi_filled_pulse->GetXaxis()->SetRange((*iter).getPulseStartBin(),(*iter).getPulseEndBin());
      wfi->Draw();
      wfi_filled_pulse->Draw("same");
      wfi->Draw("same");
      
      // arrow pointing to the found channel 1 pulse
      TArrow* thisPulse = new TArrow((*iter).getPulseHeightTime(),(*iter).getPulseHeight()-0.06,(*iter).getPulseHeightTime(),(*iter).getPulseHeight()-0.01);
      thisPulse->SetLineColor(kRed);
      thisPulse->SetLineWidth(3);
      thisPulse->Draw();
      
      TLine *thisThreshold=new TLine(wfi->GetBinCenter((*iter).getBin()-halfbinrange),b1chb_wf.getBaseline()-b1chb_wf.getThreshold(),wfi->GetBinCenter((*iter).getBin()+halfbinrange),b1chb_wf.getBaseline()-b1chb_wf.getThreshold());
      thisThreshold->SetLineWidth(1);
      thisThreshold->SetLineStyle(2);
      thisThreshold->SetLineColor(kBlue);
      thisThreshold->Draw();
      
      TLine *thisBaseline=new TLine(wfi->GetBinCenter((*iter).getBin()-halfbinrange),b1chb_wf.getBaseline(),wfi->GetBinCenter((*iter).getBin()+halfbinrange),b1chb_wf.getBaseline());
      thisBaseline->SetLineWidth(1);
      thisBaseline->SetLineStyle(2);
      thisBaseline->SetLineColor(kMagenta);
      thisBaseline->Draw();
      
      c2->Update();
      
      getchar();
      ++npeak;
    }
  }
  
  L.erase(L.begin(),L.end());
  c1->Update();
}

void doanalysis(void){
  
  // the pmt waveform
  MaxCamWaveformTools b1chb_wf(sd->event()->trigger_group(jj)->getWaveform("B1_CHB")->getWaveform(),false,10,1);
  list<MaxCamPulse> L;
  b1chb_wf.getCopyOfPulseList(L);
  list< MaxCamPulse >::iterator iter=L.begin();
  //  cout << "void doanalysis(void): (*iter).getPulseIntegral()=" << (*iter).getPulseIntegral() << endl;
  L.erase(L.begin(),L.end());
  
  vpmt=(*iter).getPulseHeight();
  qpmt=(*iter).getPulseIntegral()*1.e9;
  
  // the anode waveform
  MaxCamWaveformTools b1cha_wf(sd->event()->trigger_group(jj)->getWaveform("B1_CHA")->getWaveform(),false,10,1);
  //  cout << "b1cha_wf.getPeakHeight()=" << b1cha_wf.getPeakHeight() << endl;
  
  qanode=b1cha_wf.getPeakHeight();
  
  analysistree->Fill();
}

void drawscope(void){
  // had to change line 435 in DmtpcSkimEvent.hh...figure out why...
  c1->cd(1);
  sd->event()->trigger_group(jj)->getWaveform("B1_CHA")->getWaveform()->Draw();
  
  c1->cd(2);
  sd->event()->trigger_group(jj)->getWaveform("B1_CHB")->getWaveform()->SetLineColor(kRed);
  sd->event()->trigger_group(jj)->getWaveform("B1_CHB")->getWaveform()->Draw();
  
  c1->Update();
}

// drawing code --> Kelsey
// rise times with current skim files, not this code --> Kelsey
// check in code --> how to run to Kelsey
