#include "DmtpcStitchGUIFrame.hh"
#include "../DmtpcCameraMap.hh"
#include "../DmtpcDataset.hh"
#include "../MaxCamImageTools.hh"
#include "TF1.h" 
#include "TLine.h" 
#include "TGLayout.h"
#include "TEllipse.h"
#include "TGButtonGroup.h"
#include "TThread.h" 
#include "../MaxCamConfig.hh"
#include "TCanvas.h"
#include "TGFileDialog.h"
#include <iostream>


using namespace std; 

static const double DEFAULT_LENS_CORRECTION[] =  {1,0,3.2e-8,0,0} ; 
static const char * filetypes_raw[] = { "Raw Files","*.root",0,0};
static const char * filetypes_map[] = { "Map Files","*.map",0,0};
static const char * filetypes_overlay[] = { "Overlay Files","*.overlay",0,0};
static const char * filetypes_stitch[] = { "Stitch Files","*stitch.root",0,0};
static const char * INIT_DIR = "/net/zwicky/esata01/dmtpc/data/4sh/"; 
static const char * STITCH_DIR = "/net/zwicky/esata01/dmtpc/production/stitch/"; 
static const char * MAP_DIR = "/net/zwicky/esata01/dmtpc/production/map/"; 
static const char * CAMERA_MAP = "/net/zwicky/esata01/dmtpc/production/map/default.map"; 

/*** M#!@$EW)$*(!@##*@!#@!#@!# **/ 
inline static void ThreadFix()
{
#if ROOT_VERSION_CODE == ROOT_VERSION(5,32,0)
   void ** tsdptr  = TThread::Tsd(0,ROOT::kDirectoryThreadSlot); 
   if (tsdptr) *tsdptr = gROOT; 
#endif
}

ClassImp(DmtpcStitchGUIFrame); 

DmtpcStitchGUIFrame::DmtpcStitchGUIFrame(const TGWindow * p, UInt_t w, UInt_t h) : TGMainFrame(p,w,h) 
{

  mapname = strdup(CAMERA_MAP); 
  doGUI(); 
  SetWindowName("Stitch GUI"); 
  MapSubwindows();
  Resize(GetDefaultSize()); 
  MapWindow(); 



  stitch = 0; 
  stitched = 0; 
  filename = 0; 
  for (int i = 0; i < 4; i++)
  {
    summed[i] = 0; 
  }
  SetPreLoad(); 

}

void DmtpcStitchGUIFrame::SaveStitched()
{
  if (!stitch) return; 

  TGFileInfo fis; 
  fis.fFileTypes=filetypes_raw; 
  fis.fOverwrite = 1; 
  fis.fFilename = strdup("stitched.root"); 
  new TGFileDialog(gClient->GetRoot(),this, kFDSave, &fis); 
  TFile out (fis.fFilename,"RECREATE"); 
  out.cd(); 
  TH2 *stitchedcpy = (TH2*) stitched->Clone("stitched"); 
  stitchedcpy->Write(); 
  out.Close(); 
}


void DmtpcStitchGUIFrame::LoadStitch()
{
  TGFileInfo fio; 
  fio.fFileTypes=filetypes_stitch; 
  fio.fIniDir= strdup(STITCH_DIR); 
  new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fio); 
  if (fio.fFilename == NULL || strcmp(fio.fFilename,"")==0) return; 
  cout << "Opening stitch file from" << fio.fFilename << endl; 

  TFile in(fio.fFilename); 
  if (stitch) stitch->Delete(); 
  stitch = (Dmtpc4ShooterStitcher*) in.Get("stitch"); 
  in.Close(); 

  SetTrained(); 
  DrawStitched(); 
}

void DmtpcStitchGUIFrame::SaveStitch()
{
  if (!stitch) return; 

  TGFileInfo fis; 
  fis.fFileTypes=filetypes_stitch; 
  fis.fOverwrite = 1; 
  fis.fIniDir = strdup(STITCH_DIR); 
  char buf[64]; 
  sprintf(buf,"dmtpc_4sh_%05dstitch.root",runNumber); 
  fis.fFilename = strdup(buf); 
  new TGFileDialog(gClient->GetRoot(),this, kFDSave, &fis); 
  cout <<"Saving to " << fis.fFilename << endl; 

  stitch->Dump(); 
  TFile f(fis.fFilename,"RECREATE");
  f.cd(); 
  for (int i = 0; i < 4; i++)
  {
//    stitch->median(i)->SetDirectory(&f); 
  }
  stitch->Write("stitch"); 
  f.Close(); 
}

void *  DmtpcStitchGUIFrame::doLoading(void *p)
{
  ThreadFix(); 
//  TThread::Lock(); 
  DmtpcStitchGUIFrame * f = (DmtpcStitchGUIFrame *) p; 
  DmtpcDataset d; 
  d.openRootFile(f->filename); 

  d.getEvent(0); 
  f->runNumber = d.event()->runNumber(); 
  f->ncam = d.event()->ccdData()->GetEntries(); 


  for (int c = 0; c < f->ncam; c++)
  {

    char buf[8]; 
    sprintf(buf,"sum%d",c); 
    gROOT->cd(); 

    if (f->summed[c])
    {
      delete f->summed[c]; 
    }

    f->summed[c] = NEW_HIST2D_WITH_SAME_SIZE(d.event()->ccdData(c),TH2D,buf); 
    f->cameraNames[c] = d.event()->ccdConfig(c)->serialNumber; 
    f->summed[c]->SetTitle(f->cameraNames[c]); 
    f->ccdnames[c]->SetText(f->cameraNames[c].Data()); 
  }

 // TThread::UnLock(); 
  
  for (int i = 0; i < d.tree()->GetEntries(); i++)
  {
    d.getEvent(i); 

    for (int c = 0; c < f->ncam; c++)
    {
      f->summed[c]->Add(d.event()->ccdData(c)); 
    }
  }


  int max_i = 0; 
  double max_integral = -1; 
  for (int c = 0; c < f->ncam; c++)
  {
    f->summed[c]->Add(d.getBiasFrame(c+1), -d.tree()->GetEntries()); 
    f->summed[c]->Scale(1./d.tree()->GetEntries()); 

    if (f->summed[c]->Integral() > max_integral)
    {
      max_i = c; 
      max_integral = f->summed[c]->Integral(); 
    }
  }

//  TThread::Lock(); 

  f->ccdWithLed->SetNumber(max_i); 

  f->DrawRaw(); 

  f->SetUnbusy(); 
  f->SetPreTrained(); 

//  TThread::UnLock(); 
  cout << "Done" << endl; 

  gROOT->cd(); 
  return 0; 
}

void *  DmtpcStitchGUIFrame::doAdjusting(void *p)
{
  ThreadFix(); 
//  TThread::Lock(); 
  DmtpcStitchGUIFrame *f = (DmtpcStitchGUIFrame*) p; 

  double weights[4]; 
  double scales[4]; 
  double rots[4]; 
  double xcenters[4]; 
  double ycenters[4]; 

  vector<const TH2*> images; 
  for (int i = 0; i < f->ncam; i++)
  {
    weights[i] = f->weights[i]->GetNumber(); 
    scales[i] = f->scale[i]->GetNumber(); 
    xcenters[i] = f->xcenter[i]->GetNumber(); 
    ycenters[i] = f->ycenter[i]->GetNumber(); 
    rots[i] = f->rot[i]->GetNumber(); 
    images.push_back(f->summed[i]); 
  }

  f->stitch->setWeights(weights); 
  f->stitch->setScales(scales,false); 
  f->stitch->setCenters(xcenters,ycenters,false); 
  f->stitch->setRotations(rots,true); 
//  f->stitched->SetName("old_stitched"); 
//  TThread::UnLock(); 

  TH2 * new_stitched = f->stitch->stitch(&images); 
  cout << "Done Stitching!" <<endl; 

  TThread::Lock(); 
  if (f->stitched) delete f->stitched; 
  f->stitched = new_stitched; 

  f->DrawStitched(); 
  f->SetUnbusy(); 
  TThread::UnLock(); 

  return 0; 
}

void *  DmtpcStitchGUIFrame::doTraining(void *p)
{

 //Fix for bug in 5.32
  ThreadFix(); 

  DmtpcStitchGUIFrame *f = (DmtpcStitchGUIFrame*) p; 

//  TThread::Lock(); 
  if (f->stitch) f->stitch->Delete(); 


  f->stitch = new Dmtpc4ShooterStitcher("stitch"); 

  DmtpcCameraMap map; 
  map.loadMap(f->mapname); 


  vector<const TH2*> images; 
  TString rotation_guess[4]; 

  for (int i = 0; i < f->ncam; i++)
  {
    images.push_back(f->summed[i]); 
    rotation_guess[i] = map.getRot(f->cameraNames[i]);
  }


  double lens_coeff[LENS_CORRECTION_ORDER+1]; 
  int max_order = 0; 
  for (int i = 0; i < LENS_CORRECTION_ORDER+1; i++)
  {
     lens_coeff[i] = f->lensCorrection[i]->GetNumber(); 
     if (lens_coeff[i]!=0) max_order = i; 
  }
  DmtpcLensCorrection * l = new DmtpcLensCorrection("lens",max_order, lens_coeff); 


  f->stitch->setLensCorrection(l); 
  f->stitch->setBlurLevel(f->blurLevel->GetNumber()); 
  f->stitch->setEdgeLowThreshold(f->edgeLowThreshold->GetNumber()); 
  f->stitch->setEdgeHighThreshold(f->edgeHighThreshold->GetNumber()); 
  f->stitch->setMinEdgeNeighbors(f->minEdgeNeigbhors->GetIntNumber()); 
  f->stitch->setLinearHoughRbins(f->linHoughRbins->GetIntNumber()); 
  f->stitch->setLinearHoughThetabins(f->linHoughThetabins->GetIntNumber()); 
  f->stitch->setLinearHoughMinVotes(f->linHoughMinVotes->GetIntNumber()); 
  f->stitch->setSpacerJoinThresholds(f->spacerJoinRthresh->GetNumber(), f->spacerJoinThetathresh->GetNumber()); 
  f->stitch->setCircularHoughFirstPassNbins(f->circHoughFirstPassNbinsX->GetIntNumber(), 
                                         f->circHoughFirstPassNbinsY->GetIntNumber(), 
                                         f->circHoughFirstPassNbinsZ->GetIntNumber()); 

  f->stitch->setCircularHoughSecondPassNbins(f->circHoughSecondPassNbinsX->GetIntNumber(), 
                                          f->circHoughSecondPassNbinsY->GetIntNumber(), 
                                          f->circHoughSecondPassNbinsZ->GetIntNumber()); 

  f->stitch->setCircularHoughMins( f->circHoughMinX->GetNumber(), f->circHoughMinY->GetNumber(), f->circHoughMinZ->GetNumber()); 
  f->stitch->setCircularHoughMaxs( f->circHoughMaxX->GetNumber(), f->circHoughMaxY->GetNumber(), f->circHoughMaxZ->GetNumber()); 

  f->stitch->setNWidthsSecondPass( f->nwidthsSecondPass->GetIntNumber()); 
  f->stitch->setNSpectrPeaksR( f->nSpectrPeaksR->GetIntNumber()); 

  f->stitch->setMedianNbins( f->medianNBins->GetIntNumber()); 
  f->stitch->setMedianNIter( f->medianNiter->GetIntNumber()); 
  f->stitch->setCCDWithLED( f->ccdWithLed->GetIntNumber()); 
  f->stitch->setLEDBorderWidth( f->LEDBorderWidth->GetIntNumber()); 
  f->stitch->setLEDThresh( f->LEDThresh->GetNumber()); 
//  f->stitch->setImageHighThresh( f->ImageHighThresh->GetNumber()); 

//  TThread::UnLock(); 
  f->stitch->setScaleMethod((Dmtpc4ShooterStitcher::SCALE_METHOD) f->radioChoice); 
  f->stitch->train(&images,rotation_guess, f->cameraNames,true); 

//  TThread::Lock();
//  f->stitched->SetName("old_stitched"); 
//  TThread::UnLock(); 

  cout << "Done Training!" <<endl; 
  TH2 * new_stitched = f->stitch->stitch(&images); 
  cout << "Done Stitching!" <<endl; 
 


  for (int i = 0; i < f->ncam; i++)
  {
    f->weights[i]->SetNumber(1.); 
    f->xcenter[i]->SetNumber(f->stitch->xCenter(i)); 
    f->ycenter[i]->SetNumber(f->stitch->yCenter(i)); 
    f->scale[i]->SetNumber(f->stitch->getScale(i)); 
    f->rot[i]->SetNumber(f->stitch->getRotation(i)); 
  }
  for (int i = f->ncam; i < 4; i++)
  {
    f->ccdnames[i]->SetText("NULL"); 
  }


  TThread::Lock(); 
 
  if (f->stitched) delete f->stitched; 
  f->stitched = new_stitched; 
  cout << "HERE" << endl; 
  f->DrawStitched(); 
  cout << "THERE" << endl; 
  f->SetUnbusy(); 
  TThread::UnLock(); 

  f->SetTrained(); 
  return 0;
}

void DmtpcStitchGUIFrame::Train()
{

 
 SetBusy("Training... Be patient! (See Terminal for Output)");
 TThread * t = new TThread(doTraining, this); 
 t->Run(); 
//doTraining(this); 
}

void DmtpcStitchGUIFrame::SaveOverlays() 
{

  if (!stitch) return; 

  for (int i = 0; i < ncam; i++)
  {
    TGFileInfo fis; 
    fis.fFileTypes=filetypes_overlay; 
    fis.fOverwrite = 1; 
    char buf[64]; 
    sprintf(buf,"%s.overlay", cameraNames[i].Data()); 
    fis.fFilename = strdup(buf); 
    new TGFileDialog(gClient->GetRoot(),this, kFDSave, &fis); 
    if (fis.fFilename) 
    {
      stitch->writeOverlayFile(fis.fFilename, i); 
    }
    else
    {
      break; 
    }
  }
}




void DmtpcStitchGUIFrame::DrawRaw() 
{
  if (!summed[0]) return; 
  canvas->GetCanvas()->Clear(); 
  canvas->GetCanvas()->Divide(2,2); 
  for (int i = 0; i < ncam; i++)
  {
    if (summed[i] == 0) continue; 

    canvas->GetCanvas()->cd(i+1); 
    summed[i]->DrawCopy("colz"); 
  }

  if (stitch)
  {
    if (drawRings->IsDown()) drawRingsMulti(); 
    if (drawSpacers->IsDown()) drawSpacersMulti(); 
  }
 canvas->GetCanvas()->Update(); 
}

void DmtpcStitchGUIFrame::DrawStitched() 
{
  canvas->GetCanvas()->Clear(); 
  canvas->GetCanvas()->cd(); 
  if (stitched) 
  {
     stitched->SetMaximum(LEDThresh->GetNumber()/10.); 
     stitched->DrawCopy("colz");   

     if (drawRings->IsDown())
     {

      double inner_radius = 0; 
      double outer_radius = 0; 
      for (int i = 0; i < ncam; i++)
      {
        inner_radius += stitch->innerRadius(i); 
        outer_radius += stitch->outerRadius(i); 
      }

      inner_radius/= ncam; 
      outer_radius/= ncam; 

      TEllipse * inner = new TEllipse(0,0, inner_radius, inner_radius); 
      inner->SetFillStyle(0); 
      inner->SetLineColor(0); 
      inner->SetLineWidth(2); 
      TEllipse * outer = new TEllipse(0,0, outer_radius, outer_radius); 
      outer->SetFillStyle(0); 
      outer->SetLineColor(0); 
      outer->SetLineWidth(3); 
      inner->Draw("same"); 
      outer->Draw("same"); 
     }

     if (drawSpacers->IsDown())
     {
        for (int i = 0; i < ncam; i++)
        {
          double c = stitch->getScale(i); 
          double th = stitch->getRotation(i); 
          double x0 = stitch->xCenter(i); 
          double y0 = stitch->yCenter(i); 
        
          for (unsigned s = 0; s < stitch->getNSpacers(i); s++ )
          {
            double m = stitch->getSpacerSlope(i,s); 
            double b = stitch->getSpacerIntercept(i,s); 
            TF1 f("f","1/(cos([1]) + [0]*sin([1])) * (x * ([0] * cos([1]) - sin([1])) - [2] * [3] + [0]*[3]*[4] + [3]*[5])",-c * (x0*cos(th) + y0*sin(th)), c * ((1024-x0)*cos(th) + (1024-y0)*sin(th)) ); 
            f.SetParameters(m,th,y0,c,x0,b); 
            f.SetLineColor(20+i); 
            f.SetLineWidth(2); 
            f.SetLineStyle(7); 
            f.DrawCopy("same"); 
          }
        }
     }

  }

  canvas->GetCanvas()->Update(); 
}

void DmtpcStitchGUIFrame::DrawEdges() 
{
  canvas->GetCanvas()->Clear(); 
  if (stitch)
  {
    canvas->GetCanvas()->Divide(2,2); 
    for (int i = 0; i < ncam; i++)
    {
       canvas->GetCanvas()->cd(i+1); 
       stitch->edge(i)->DrawCopy("colz"); 
    }

    if (drawRings->IsDown()) drawRingsMulti(); 
    if (drawSpacers->IsDown()) drawSpacersMulti(); 
  }
  canvas->GetCanvas()->Update(); 
}


void DmtpcStitchGUIFrame::DrawPolar()
{

  canvas->GetCanvas()->Clear(); 
  if (stitch)
  {
    canvas->GetCanvas()->Divide(2,2); 
    for (int i = 0; i < ncam; i++)
    {
      canvas->GetCanvas()->cd(i+1); 
      stitch->polar(i)->DrawCopy("colz"); 
      if (drawRings->IsDown())
      {
          TLine * inner = new TLine(stitch->innerRadius(i),0,stitch->innerRadius(i),1024); 
          TLine * outer = new TLine(stitch->outerRadius(i),0,stitch->outerRadius(i),1024); 
          inner->SetLineWidth(2); 
          inner->SetLineStyle(7); 
          inner->Draw("same"); 
          outer->SetLineWidth(2); 
          outer->SetLineStyle(7); 
          outer->Draw("same"); 
      }
    }
  }



  canvas->GetCanvas()->Update(); 

}

void DmtpcStitchGUIFrame::DrawMedians() 
{
  canvas->GetCanvas()->Clear(); 
  if (stitch)
  {
    canvas->GetCanvas()->Divide(2,2); 
    for (int i = 0; i < ncam; i++)
    {
      canvas->GetCanvas()->cd(i+1); 
      stitch->median(i)->DrawCopy("colz"); 
    }

    if (drawRings->IsDown()) drawRingsMulti(); 
    if (drawSpacers->IsDown()) drawSpacersMulti(); 
  }
  canvas->GetCanvas()->Update(); 

}


void DmtpcStitchGUIFrame::SetPreLoad() 
{
  adjustButton->SetEnabled(false); 
  trainButton->SetEnabled(false); 
  saveStitchButton->SetEnabled(false); 
  saveStitchedButton->SetEnabled(false); 
  saveOverlayButton->SetEnabled(false); 
  loadStitchButton->SetEnabled(false); 
  busy->SetText("Must load a file"); 
}

void DmtpcStitchGUIFrame::SetPreTrained() 
{
  adjustButton->SetEnabled(false); 
  saveStitchButton->SetEnabled(false); 
  saveStitchedButton->SetEnabled(false); 
  saveOverlayButton->SetEnabled(false); 
  loadStitchButton->SetEnabled(true); 
}

void DmtpcStitchGUIFrame::SetTrained() 
{
  adjustButton->SetEnabled(true); 
  saveStitchButton->SetEnabled(true); 
  saveStitchedButton->SetEnabled(true); 
  saveOverlayButton->SetEnabled(true); 
}


void DmtpcStitchGUIFrame::Map()
{
  
  TGFileInfo fio;  
  fio.fFileTypes = filetypes_map; 
  fio.fIniDir = strdup(MAP_DIR); 

  new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fio); 
  
  if (fio.fFilename == NULL || strcmp(fio.fFilename,"")==0) return; 
  free(mapname); 
  mapname = strdup(fio.fFilename); 
  mapLabel->SetText(fio.fFilename); 
  Layout(); 
  gSystem->ProcessEvents();
}


void DmtpcStitchGUIFrame::Load()
{
  TGFileInfo fio; 
  fio.fFileTypes=filetypes_raw; 
  fio.fIniDir= strdup(INIT_DIR); 
  new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fio); 
  if (fio.fFilename == NULL || strcmp(fio.fFilename,"")==0) return; 
  cout << "Opening " << fio.fFilename << endl; 
  if (filename) free(filename); 
  filename = strdup(fio.fFilename); 
 
  fileLabel->SetText(fio.fFilename); 
  Layout(); 
  gSystem->ProcessEvents();

  SetBusy("Loading and summing File!"); 

  TThread * t = new TThread(doLoading, this); 
  t->Run(); 
}


void DmtpcStitchGUIFrame::Adjust()
{

  SetBusy("Adjusting and stitching File!"); 
  TThread * t = new TThread(doAdjusting, this); 
  t->Run(); 

}



void DmtpcStitchGUIFrame::SetBusy( const char * text)
{
  
  trainButton->SetEnabled(false); 
  loadButton->SetEnabled(false); 
  mapButton->SetEnabled(false); 
  loadStitchButton->SetEnabled(false); 
  adjustButton->SetEnabled(false); 
  busy->SetText(text); 
  gSystem->ProcessEvents(); 
}

void DmtpcStitchGUIFrame::SetUnbusy()
{
  
  trainButton->SetEnabled(true); 
  loadButton->SetEnabled(true); 
  adjustButton->SetEnabled(true); 
  mapButton->SetEnabled(true); 
  loadStitchButton->SetEnabled(true); 
  busy->SetText("Idle..."); 
  gSystem->ProcessEvents(); 
}




void DmtpcStitchGUIFrame::doGUI()
{


  TGHorizontalFrame * base = new TGHorizontalFrame(this,1200,800); 


  TGVerticalFrame * left = new TGVerticalFrame(base,400,800); 


  /*******Train Frame *********/
  trainFrame = new TGGroupFrame(left,"Train"); 
  
  TGGroupFrame * loadFrame = new TGGroupFrame(trainFrame,"File Choice"); 
  TGHorizontalFrame * h0 = new TGHorizontalFrame(loadFrame,400,40); 
  loadButton = new TGTextButton(h0,"Data File..."); 
  loadButton->Connect("Clicked()","DmtpcStitchGUIFrame",this,"Load()"); 
  fileLabel = new TGLabel(h0,"no file loaded "); 
  h0->AddFrame(loadButton); 
  h0->AddFrame(fileLabel, new TGLayoutHints(kLHintsExpandX,3,3,3,3)); 
  loadFrame->AddFrame(h0); 

  TGHorizontalFrame * h01 = new TGHorizontalFrame(loadFrame,400,40); 
  mapButton = new TGTextButton(h01,"Choose Map..."); 
  mapButton->Connect("Clicked()","DmtpcStitchGUIFrame",this,"Map()"); 
  mapLabel = new TGLabel(h01,mapname); 
  h01->AddFrame(mapButton); 
  h01->AddFrame(mapLabel); 
  loadFrame->AddFrame(h01); 

  TGHorizontalFrame * h02 = new TGHorizontalFrame(loadFrame,400,40); 
  loadStitchButton = new TGTextButton(h02,"Load Stitch File..."); 
  loadStitchButton->Connect("Clicked()","DmtpcStitchGUIFrame",this,"LoadStitch()"); 
  h02->AddFrame(loadStitchButton); 
  loadFrame->AddFrame(h02); 



  trainFrame->AddFrame(loadFrame, new TGLayoutHints(kLHintsExpandX)); 

  //lens frame
  TGGroupFrame * lensFrame = new TGGroupFrame(trainFrame,"Lens Correction Coefficients"); 

  TGHorizontalFrame * h1 = new TGHorizontalFrame(lensFrame,400,40); 
  for (int i = 0; i <=4; i++)
  {
    lensCorrection[i] = new TGNumberEntry(h1,  DEFAULT_LENS_CORRECTION[i],5,-1, TGNumberFormat::kNESReal); 
    lensCorrection[i]->Resize(60,15); 
    h1->AddFrame(lensCorrection[i]); 
  }
  lensFrame->AddFrame(h1); 
  trainFrame->AddFrame(lensFrame, new TGLayoutHints(kLHintsExpandX)); 

  //edge frame
  TGGroupFrame * edgeFrame = new TGGroupFrame(trainFrame,"Edge Detection Settings"); 
  TGHorizontalFrame *h2 = new TGHorizontalFrame(edgeFrame,400,40); 
  h2->AddFrame(new TGLabel(h2, "Blur: ")); 
  blurLevel = new TGNumberEntry(h2, 2.4,4,-1, TGNumberFormat::kNESReal); 
  blurLevel->Resize(40,15); 
  h2->AddFrame(blurLevel); 
  h2->AddFrame(new TGLabel(h2, " Low Thr: ")); 
  edgeLowThreshold = new TGNumberEntry(h2,-0.5,4,-1, TGNumberFormat::kNESReal); 
  edgeLowThreshold->Resize(40,15); 
  h2->AddFrame(edgeLowThreshold); 
  h2->AddFrame(new TGLabel(h2, " High Thr: ")); 
  edgeHighThreshold = new TGNumberEntry(h2, 0.1,4,-1, TGNumberFormat::kNESReal); 
  edgeHighThreshold->Resize(40,15); 
  h2->AddFrame(edgeHighThreshold); 
  h2->AddFrame(new TGLabel(h2, " Min Ngbrs: ")); 
  minEdgeNeigbhors = new TGNumberEntry(h2,100,4,-1, TGNumberFormat::kNESReal); 
  minEdgeNeigbhors->Resize(40,15); 
  h2->AddFrame(minEdgeNeigbhors); 
  edgeFrame->AddFrame(h2); 
  trainFrame->AddFrame(edgeFrame, new TGLayoutHints(kLHintsExpandX)); 

  //join frame
  TGGroupFrame * joinFrame = new TGGroupFrame(trainFrame,"Spacer Join Settings"); 
  TGHorizontalFrame *h3 = new TGHorizontalFrame(joinFrame,400,40); 
  h3->AddFrame(new TGLabel(h3, " R Thr: ")); 
  spacerJoinRthresh = new TGNumberEntry(h3, 40,4,-1, TGNumberFormat::kNESReal); 
  spacerJoinRthresh->Resize(40,15); 
  h3->AddFrame(spacerJoinRthresh); 
  h3->AddFrame(new TGLabel(h3, " Theta Thr: ")); 
  spacerJoinThetathresh = new TGNumberEntry(h3,0.4,4,-1, TGNumberFormat::kNESReal); 
  spacerJoinThetathresh->Resize(40,15); 
  h3->AddFrame(spacerJoinThetathresh); 
  joinFrame->AddFrame(h3); 
  trainFrame->AddFrame(joinFrame, new TGLayoutHints(kLHintsExpandX)); 

  //linear hough frame 
  TGGroupFrame * linhoughframe = new TGGroupFrame(trainFrame,"Linear Hough Transform Settings"); 
  TGHorizontalFrame *h4 = new TGHorizontalFrame(linhoughframe,400,40); 
  h4->AddFrame(new TGLabel(h4,"NBinsR: ")); 
  linHoughRbins = new TGNumberEntry(h4,1024,5,-1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive); 
  linHoughRbins->Resize(60,15); 
  h4->AddFrame(linHoughRbins); 
  h4->AddFrame(new TGLabel(h4," NBinsTheta: ")); 
  linHoughThetabins = new TGNumberEntry(h4,1024,5,-1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive); 
  linHoughThetabins->Resize(60,15); 
  h4->AddFrame(linHoughThetabins); 
  h4->AddFrame(new TGLabel(h4," MinVotes: ")); 
  linHoughMinVotes = new TGNumberEntry(h4,400,4,-1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive); 
  linHoughMinVotes->Resize(60,15); 
  h4->AddFrame(linHoughMinVotes); 
  linhoughframe->AddFrame(h4); 
  trainFrame->AddFrame(linhoughframe, new TGLayoutHints(kLHintsExpandX)); 


  TGGroupFrame * circhoughframe = new TGGroupFrame(trainFrame,"Circular Hough Transform Settings"); 
  TGHorizontalFrame * h5 = new TGHorizontalFrame(circhoughframe, 400,200); 
  h5->AddFrame(new TGLabel(h5, "1st Pass Nbins(xyr): ")); 
  circHoughFirstPassNbinsX = new TGNumberEntry(h5,40,4,-1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive);
  circHoughFirstPassNbinsX->Resize(40,15); 
  h5->AddFrame(circHoughFirstPassNbinsX); 
  circHoughFirstPassNbinsY = new TGNumberEntry(h5,40,4,-1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive);
  circHoughFirstPassNbinsY->Resize(40,15); 
  h5->AddFrame(circHoughFirstPassNbinsY); 
  circHoughFirstPassNbinsZ = new TGNumberEntry(h5,40,4,-1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive);
  circHoughFirstPassNbinsZ->Resize(40,15); 
  h5->AddFrame(circHoughFirstPassNbinsZ); 
  circhoughframe->AddFrame(h5); 
  TGHorizontalFrame * h6 = new TGHorizontalFrame(circhoughframe, 400,200); 
  h6->AddFrame(new TGLabel(h6, "2nd Pass Nbins(xyr): ")); 
  circHoughSecondPassNbinsX = new TGNumberEntry(h6,40,4,-1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive);
  circHoughSecondPassNbinsX->Resize(60,15); 
  h6->AddFrame(circHoughSecondPassNbinsX); 
  circHoughSecondPassNbinsY = new TGNumberEntry(h6,40,4,-1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive);
  circHoughSecondPassNbinsY->Resize(60,15); 
  h6->AddFrame(circHoughSecondPassNbinsY); 
  circHoughSecondPassNbinsZ = new TGNumberEntry(h6,2400,5,-1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive);
  circHoughSecondPassNbinsZ->Resize(60,15); 
  h6->AddFrame(circHoughSecondPassNbinsZ); 
  circhoughframe->AddFrame(h6); 
  TGHorizontalFrame * h7 = new TGHorizontalFrame(circhoughframe, 400,200); 
  h7->AddFrame(new TGLabel(h7, "Min(xyr): ")); 
  circHoughMinX = new TGNumberEntry(h7,-200,4,-1,  TGNumberFormat::kNESReal);
  circHoughMinX->Resize(60,15); 
  h7->AddFrame(circHoughMinX); 
  circHoughMinY = new TGNumberEntry(h7,800, 4,-1, TGNumberFormat::kNESReal);
  circHoughMinY->Resize(60,15); 
  h7->AddFrame(circHoughMinY); 
  circHoughMinZ = new TGNumberEntry(h7,800, 4,-1, TGNumberFormat::kNESReal);
  circHoughMinZ->Resize(60,15); 
  h7->AddFrame(circHoughMinZ); 
  circhoughframe->AddFrame(h7); 
   TGHorizontalFrame * h77 = new TGHorizontalFrame(circhoughframe, 400,200); 
  h77->AddFrame(new TGLabel(h77, "Max(xyr): ")); 
  circHoughMaxX = new TGNumberEntry(h77,200, 4,-1,TGNumberFormat::kNESReal);
  circHoughMaxX->Resize(60,15); 
  h77->AddFrame(circHoughMaxX); 
  circHoughMaxY = new TGNumberEntry(h77,1200,4,-1, TGNumberFormat::kNESReal);
  circHoughMaxY->Resize(60,15); 
  h77->AddFrame(circHoughMaxY); 
  circHoughMaxZ = new TGNumberEntry(h77,1600,4,-1, TGNumberFormat::kNESReal);
  circHoughMaxZ->Resize(60,15); 
  h77->AddFrame(circHoughMaxZ); 
  circhoughframe->AddFrame(h77); 
  TGHorizontalFrame * h8 = new TGHorizontalFrame(circhoughframe, 400,200); 
  h8->AddFrame(new TGLabel(h8, "NWidthsSecondPass: ")); 
  nwidthsSecondPass = new TGNumberEntry(h8,1,2,-1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive); 
  nwidthsSecondPass->Resize(40,15); 
  h8->AddFrame(nwidthsSecondPass); 
  h8->AddFrame(new TGLabel(h8, " NSpectrPeaksR: ")); 
  nSpectrPeaksR = new TGNumberEntry(h8,4,2,-1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive); 
  nSpectrPeaksR->Resize(40,15); 
  h8->AddFrame(nSpectrPeaksR); 
  circhoughframe->AddFrame(h8); 
  trainFrame->AddFrame(circhoughframe, new TGLayoutHints(kLHintsExpandX)); 


  TGGroupFrame * ledframe = new TGGroupFrame(trainFrame,"LED Settings"); 

  TGHorizontalFrame * h9 = new TGHorizontalFrame(ledframe, 400,200); 
  h9->AddFrame(new TGLabel(h9, "CCD #: ")); 
  ccdWithLed = new TGNumberEntry(h9,0,1,-1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive); 
  ccdWithLed->Resize(40,15); 
  h9->AddFrame(ccdWithLed); 
  h9->AddFrame(new TGLabel(h9, " Border: ")); 
  LEDBorderWidth = new TGNumberEntry(h9,30,3,-1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative); 
  LEDBorderWidth->Resize(40,15); 
  h9->AddFrame(LEDBorderWidth); 
  h9->AddFrame(new TGLabel(h9, " Thr: ")); 
  LEDThresh = new TGNumberEntry(h9,6000,5,-1, TGNumberFormat::kNESReal); 
  LEDThresh->Resize(60,15); 
  h9->AddFrame(LEDThresh); 
  ledframe->AddFrame(h9); 
  trainFrame->AddFrame(ledframe, new TGLayoutHints(kLHintsExpandX)); 

  TGGroupFrame * medianframe = new TGGroupFrame(trainFrame, "Median Settings"); 
  TGHorizontalFrame *h10 = new TGHorizontalFrame(medianframe, 400,200); 
  h10->AddFrame(new TGLabel(h10, "Nbins: ")); 
  medianNBins = new TGNumberEntry(h10,3,2,-1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive); 
  medianNBins->Resize(40,15); 
  h10->AddFrame(medianNBins); 
  h10->AddFrame(new TGLabel(h10, " NIter: ")); 
  medianNiter = new TGNumberEntry(h10,1,2,-1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive); 
  medianNiter->Resize(40,15); 
  h10->AddFrame(medianNiter); 
//  h10->AddFrame(new TGLabel(h10, " Thr: ")); 
 // ImageHighThresh = new TGNumberEntry(h10,-1,2, TGNumberFormat::kNESReal); 
 // ImageHighThresh->Resize(40,15); 
 // h10->AddFrame(ImageHighThresh); 
  medianframe->AddFrame(h10); 
  trainFrame->AddFrame(medianframe, new TGLayoutHints(kLHintsExpandX)); 


  TGHButtonGroup * radiusframe = new TGHButtonGroup(trainFrame, "Scale Settings"); 
  radiusframe->SetBorderDrawn(true); 
  radiusframe->AddFrame(new TGRadioButton(radiusframe,"Fixed",0)); 
  radiusframe->AddFrame(new TGRadioButton(radiusframe,"Inner",1)); 
  radiusframe->AddFrame(new TGRadioButton(radiusframe,"Outer",2)); 
  radiusframe->AddFrame(new TGRadioButton(radiusframe,"Both",3)); 
  radiusframe->SetButton(0); 

  radiusframe->Connect("Clicked(Int_t)","DmtpcStitchGUIFrame",this,"HandleRadio(Int_t)"); 
  trainFrame->AddFrame(radiusframe); 

  trainButton= new TGTextButton(trainFrame,"Train"); 
  trainButton->Connect("Clicked()", "DmtpcStitchGUIFrame",this,"Train()"); 
  trainButton->SetBackgroundColor(0x0099dd); 
  trainButton->Resize(60,80); 
  trainFrame->AddFrame(trainButton, new TGLayoutHints(kLHintsExpandX)); 

  left->AddFrame(trainFrame); 


  /*** Adjust Frame ***/
  adjustframe = new TGGroupFrame(left,"Adjust"); 
  TGMatrixLayout * adjusttbl = new TGMatrixLayout(adjustframe, 7, 5,3); 
  adjustframe->SetLayoutManager(adjusttbl); 

  adjustframe->AddFrame(new TGLabel(adjustframe,"Property")); 

  for (int i = 0; i < 4; i++)
  {
    char buf[8]; 
    sprintf(buf,"CCD %d", i); 
    ccdnames[i] = new TGLabel(adjustframe,buf); 
    adjustframe->AddFrame(ccdnames[i]); 
  }

  adjustframe->AddFrame(new TGLabel(adjustframe,"Rotation:  ")); 
  for (int i = 0; i < 4; i++)
  {
     rot[i] = new TGNumberEntry(adjustframe,0,7); 
     rot[i]->Resize(80,15); 
     adjustframe->AddFrame(rot[i]); 
  }

  adjustframe->AddFrame(new TGLabel(adjustframe,"X Center: ")); 
  for (int i = 0; i < 4; i++)
  {
     xcenter[i] = new TGNumberEntry(adjustframe,0,7); 
     xcenter[i]->Resize(80,15); 
     adjustframe->AddFrame(xcenter[i]); 
  }

  adjustframe->AddFrame(new TGLabel(adjustframe,"Y Center: ")); 
  for (int i = 0; i < 4; i++)
  {
     ycenter[i] = new TGNumberEntry(adjustframe,0,7); 
     ycenter[i]->Resize(80,15); 
     adjustframe->AddFrame(ycenter[i]); 
  }

  adjustframe->AddFrame(new TGLabel(adjustframe,"Scales:    ")); 
  for (int i = 0; i < 4; i++)
  {
     scale[i] = new TGNumberEntry(adjustframe,0,7); 
     scale[i]->Resize(80,15); 
     adjustframe->AddFrame(scale[i]); 
  }

  adjustframe->AddFrame(new TGLabel(adjustframe,"Weights:   ")); 
  for (int i = 0; i < 4; i++)
  {
     weights[i] = new TGNumberEntry(adjustframe,1,7); 
     weights[i]->Resize(80,15); 
     adjustframe->AddFrame(weights[i]); 
  }
  adjustButton= new TGTextButton(adjustframe,"    Adjust   "); 
//  adjustButton->SetBackgroundColor(0x22ff22); 
  adjustButton->Connect("Clicked()", "DmtpcStitchGUIFrame",this,"Adjust()"); 
  adjustframe->AddFrame(adjustButton); 
  
  left->AddFrame(adjustframe, new TGLayoutHints(kLHintsRight | kLHintsExpandX,5,5,5,5)); 
  base->AddFrame(left, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 5,5,5,5)); 

  TGVerticalFrame * right = new TGVerticalFrame(base,800,800); 
  TGGroupFrame * canvasframe = new TGGroupFrame(right,"View"); 


  TGHorizontalFrame * buttonFrame = new TGHorizontalFrame(canvasframe); 
  
  TGButton * raw_button = new TGTextButton(buttonFrame,"Raw"); 
  TGButton * stitch_button = new TGTextButton(buttonFrame,"Stitched"); 
  TGButton * edge_button = new TGTextButton(buttonFrame,"Edge"); 
  TGButton * median_button = new TGTextButton(buttonFrame,"Cleaned"); 
  TGButton * polar_button = new TGTextButton(buttonFrame,"Polar"); 
  drawSpacers = new TGCheckButton(buttonFrame,"Draw Spacers"); 
  drawRings = new TGCheckButton(buttonFrame,"Draw Rings"); 

  buttonFrame->AddFrame(raw_button); 
  buttonFrame->AddFrame(stitch_button); 
  buttonFrame->AddFrame(edge_button); 
  buttonFrame->AddFrame(median_button); 
  buttonFrame->AddFrame(polar_button); 
  buttonFrame->AddFrame(drawSpacers); 
  buttonFrame->AddFrame(drawRings); 

  raw_button->Connect("Clicked()","DmtpcStitchGUIFrame",this,"DrawRaw()"); 
  stitch_button->Connect("Clicked()","DmtpcStitchGUIFrame",this,"DrawStitched()"); 
  edge_button->Connect("Clicked()","DmtpcStitchGUIFrame",this,"DrawEdges()"); 
  median_button->Connect("Clicked()","DmtpcStitchGUIFrame",this,"DrawMedians()"); 
  polar_button->Connect("Clicked()","DmtpcStitchGUIFrame",this,"DrawPolar()"); 

  TGHorizontalFrame * saveFrame = new TGHorizontalFrame(canvasframe); 

  saveOverlayButton = new TGTextButton(saveFrame,"Save Overlays"); 
  saveStitchedButton = new TGTextButton(saveFrame,"Save Stitched Image"); 
  saveStitchButton = new TGTextButton(saveFrame,"Save Stitch Information"); 
  saveStitchButton->SetBackgroundColor(0xff4444); 
  saveFrame->AddFrame(saveStitchedButton); 
  saveFrame->AddFrame(saveOverlayButton); 
  saveFrame->AddFrame(saveStitchButton); 
  saveStitchButton->Connect("Clicked()","DmtpcStitchGUIFrame",this,"SaveStitch()"); 
  saveStitchedButton->Connect("Clicked()","DmtpcStitchGUIFrame",this,"SaveStitched()"); 
  saveOverlayButton->Connect("Clicked()","DmtpcStitchGUIFrame",this,"SaveOverlays()"); 


  canvas = new TRootEmbeddedCanvas("ECanvas",canvasframe,800,800); 
  canvasframe->AddFrame(canvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,5,5,5,5)); 

  canvasframe->AddFrame(buttonFrame, new TGLayoutHints(kLHintsExpandX,5,5,5,5)); 
  canvasframe->AddFrame(saveFrame, new TGLayoutHints(kLHintsExpandX,5,5,5,5)); 

  right->AddFrame(canvasframe, new TGLayoutHints(kLHintsRight | kLHintsExpandX | kLHintsExpandY,5,5,5,5)); 

  base->AddFrame(right, new TGLayoutHints(kLHintsLeft | kLHintsExpandY | kLHintsExpandX, 5,5,5,5)); 
  AddFrame(base, new TGLayoutHints(kLHintsExpandY | kLHintsExpandX)); 

  busy= new TGStatusBar(this); 
  AddFrame(busy, new TGLayoutHints(kLHintsExpandX)); 
}


void DmtpcStitchGUIFrame::drawSpacersMulti()
{
  for (int i = 0; i < ncam; i++)
  {
    canvas->GetCanvas()->cd(i+1);

    for (unsigned s = 0; s < stitch->getNSpacers(i); s++)
    {
      TF1 f1("f1","pol1",0,1024); 
      f1.SetParameter(0,stitch->getSpacerIntercept(i,s)); 
      f1.SetParameter(1,stitch->getSpacerSlope(i,s)); 
      f1.DrawCopy("lsame"); 
    }
   
  }
}

void DmtpcStitchGUIFrame::drawRingsMulti()
{
  for (int i = 0; i < ncam; i++)
  {
    canvas->GetCanvas()->cd(i+1); 
    TEllipse * ei = new TEllipse(stitch->xCenter(i), stitch->yCenter(i), stitch->innerRadius(i), stitch->innerRadius(i)); 
    ei->SetFillStyle(0); 
    ei->SetLineColor(3); 
    ei->SetLineWidth(2); 
    TEllipse *eo = new TEllipse(stitch->xCenter(i), stitch->yCenter(i), stitch->outerRadius(i), stitch->outerRadius(i)); 
    eo->SetFillStyle(0); 
    eo->SetLineColor(4); 
    eo->SetLineWidth(2); 
    ei->Draw("same"); 
    eo->Draw("same");  
  }
}

void DmtpcStitchGUIFrame::HandleRadio(Int_t id)
{
     radioChoice =  id; 
}
