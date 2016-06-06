#include "DmtpcDataConverter.hh"

#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include <memory>
#include <iostream>


/** Utlity casting functions */
static inline short sh(float val)
{
  val+=0.5;
  unsigned short uval = (unsigned short) val;
  return *( (short*) &(uval) ); 
}

static inline char ch(float val, unsigned char zero)
{
  unsigned char uch =  zero  + (unsigned char) (val > 0 ? val + 0.5 : val -0.5); 
  return *(char*) &uch; 
}

static inline float fl(short val)
{
  return (float)  *( (unsigned short*) &val); 

}

void DmtpcDataConverter::convertDataFile(const char * in, const char * out, bool overwrite) 
{
  TFile *in_file = new TFile(in); 
  TTree * intree = (TTree*) in_file->Get("dmtpc"); 
  DmtpcEvent * in_ev = new DmtpcEvent; 
  intree->SetBranchAddress("event",&in_ev); 
  TFile * out_file = new TFile(out, overwrite ? "RECREATE" : "CREATE" ); 


  /* Copy everything directly except for event tree */
  TIter next(in_file->GetListOfKeys());  
  while (TKey * k = (TKey*) next())
  {
    if (k->GetName() == "dmtpc") 
      continue; 

    TObject * copy = k->ReadObj()->Clone(); 
    copy->Write(k->GetName()); 
  } 

  /* Convert event tree */
  out_file->cd(); 
  TTree * outtree = new TTree("dmtpc","dmtpc"); 
  DmtpcEvent * out_ev = new DmtpcEvent; 
  outtree->Branch("event","DmtpcEvent",&out_ev,32000,0); 

  for (int i = 0; i < intree->GetEntries(); i++)
  {
    cout << "Converting event " << i << endl; 
    intree->GetEntry(i); 
    out_ev = DmtpcDataConverter::convertHistDataTypes(in_ev); 
    outtree->GetBranch("event")->SetAddress(&out_ev); 
    outtree->Fill(); 
    out_ev->Delete(); 
  }

  outtree->Write(); 
  outtree->Delete();
  
 // out_file->Close(); 
 // in_file->Close(); 

}

DmtpcEvent * DmtpcDataConverter::convertHistDataTypes(DmtpcEvent * old, DmtpcEvent * placement_ptr)
{

  DmtpcEvent * ev; 
 
  if (placement_ptr)
    ev = new(placement_ptr) DmtpcEvent(*old,false);  
  else
    ev = new DmtpcEvent(*old,false);  
  

  for (int i = 0; i < old->rawCcdData()->GetEntries(); i++)
  {
    ccdReduce((TH2F*) old->rawCcdData(i),old->rawCcdData(i)->GetName(),(TH2S*) (*(ev->rawCcdData()))[i]); 
  } 


  for (int i = 0; i < old->rawScopeData()->GetEntries(); i++)
  {
    //scopeReduce((TH1F*) old->rawScopeData(i), old->scopeDataInfo(i),old->rawScopeData(i)->GetName(),(ScopeWaveformData*)(*(ev->rawScopeData()))[i]); 
  }

 
  for (int i = 0; i < old->rawOverscan()->GetEntries(); i++)
  {
    ccdReduce((TH2F*)old->rawOverscan(i),old->rawOverscan(i)->GetName(),(TH2S*) (*(ev->rawOverscan()))[i]); 
  }
 
  return ev; 

}


TH2S * DmtpcDataConverter::ccdReduce(const TH2F * old, const char * name, TH2S * placement_ptr) 
{

  const char * new_name = name ? name : old->GetName(); 

  TH2S * reduced; 
  if (placement_ptr)
  {
    reduced =  new (placement_ptr) TH2S(new_name,new_name,old->GetNbinsX(), old->GetXaxis()->GetXmin(), old->GetXaxis()->GetXmax(),
                                           old->GetNbinsY(), old->GetYaxis()->GetXmin(), old->GetYaxis()->GetXmax()); 
  }
  else
  {
    reduced =  new TH2S(new_name,new_name,old->GetNbinsX(), old->GetXaxis()->GetXmin(), old->GetXaxis()->GetXmax(),
                                           old->GetNbinsY(), old->GetYaxis()->GetXmin(), old->GetYaxis()->GetXmax()); 

  }

  for (int i = 1; i <= old->GetNbinsX(); i++)
  {
    for (int j = 1; j <= old->GetNbinsY(); j++)
    {
      reduced->SetBinContent(i,j,sh(old->GetBinContent(i,j))); 
    }
  }

  //cout << placement_ptr<< endl; 
  //cout << reduced<< endl; 
  return reduced;
}

TH2F * DmtpcDataConverter::ccdExpand(const TH2S * reduced, const char * name, TH2F * placement_ptr) 
{
 
  const char * new_name = name ? name : reduced->GetName(); 
  TH2F * expanded; 
  if (placement_ptr)
  {
    expanded =  new (placement_ptr) TH2F(new_name,new_name,reduced->GetNbinsX(), reduced->GetXaxis()->GetXmin(), reduced->GetXaxis()->GetXmax(),
                                         reduced->GetNbinsY(), reduced->GetYaxis()->GetXmin(), reduced->GetYaxis()->GetXmax()); 
  }
  else
  {
    expanded =  new TH2F(new_name,new_name,reduced->GetNbinsX(), reduced->GetXaxis()->GetXmin(), reduced->GetXaxis()->GetXmax(),
                                         reduced->GetNbinsY(), reduced->GetYaxis()->GetXmin(), reduced->GetYaxis()->GetXmax()); 
  }

  for (int i = 1; i <= reduced->GetNbinsX(); i++)
  {
    for (int j = 1; j <= reduced->GetNbinsY(); j++)
    {
      expanded->SetBinContent(i,j,fl((short)reduced->GetBinContent(i,j))); 
    }
  }

  if (placement_ptr) assert(expanded == placement_ptr); 

  return expanded; 
}


ScopeWaveformData * DmtpcDataConverter::scopeReduce(const TH1F *old, const ScopeDataInfo * sinfo, const char * name, TH1F * placement_ptr) 
{

  const char * new_name = name ? name : old->GetName(); 
  float vstep = sinfo->getVoltageStep(); 
  float vmin = sinfo->getVoltageMin(); 
  float vmax = sinfo->getVoltageMax(); 
  float vzero = (vmin + vmax)/2; 
  unsigned char czero = vzero == 0 ? 128 : 0; 
  ScopeWaveformData * reduced;
  double timestamp = old->GetBinContent(0); 
  if (placement_ptr)
    reduced =  new (placement_ptr) ScopeWaveformData(new_name,old->GetNbinsX(), old->GetXaxis()->GetXmin(), old->GetXaxis()->GetXmax(),timestamp, vstep, czero); 
  else
    reduced =  new ScopeWaveformData(new_name,old->GetNbinsX(), old->GetXaxis()->GetXmin(), old->GetXaxis()->GetXmax(), timestamp, vstep, czero); 

  for (int i = 1; i <= old->GetNbinsX(); i++)
  {
    reduced->SetBinContent(i,(old->GetBinContent(i)-vzero)/vstep,czero); 
//    cout << old->GetBinContent(i) << " : " << 0 + (unsigned char) reduced->GetBinContent(i) << endl; 
  }

  return reduced; 
}

TH1F * DmtpcDataConverter::scopeExpand(const ScopeWaveformData * reduced, const char * name, TH1F * placement_ptr)
{

  const char * new_name = name ? name : reduced->GetName(); 
  TH1F * expanded; 
  if (placement_ptr)
    expanded =  new (placement_ptr) TH1F(new_name,new_name,reduced->GetNbinsX(), reduced->GetXaxis()->GetXmin(), reduced->GetXaxis()->GetXmax()); 
  else
    expanded =  new TH1F(new_name,new_name,reduced->GetNbinsX(), reduced->GetXaxis()->GetXmin(), reduced->GetXaxis()->GetXmax()); 

  for (int i = 1; i <= reduced->GetNbinsX(); i++)
  {
    expanded->SetBinContent(i,((reduced->GetBinContent(i)) - reduced->getZero())*reduced->getleveltoVolt()); 
  }

  expanded->SetBinContent(0,reduced->getTimeStamp()); 


  return expanded; 

}
