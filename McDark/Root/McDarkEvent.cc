// @(#)root/test:$Name:  $:$Id: McDarkEvent.cc,v 1.1 2008/08/18 09:31:25 tcsahin Exp $
// Author: Rene Brun   19/08/96

#include "TRandom.h"
#include "TDirectory.h"
#include "TProcessID.h"

#include "McDarkEvent.hh"


//ClassImp(EventHeader)
ClassImp(McDarkEvent)
//ClassImp(Track)
//ClassImp(HistogramManager)

//TClonesArray *Event::fgTracks = 0;
//TH1F *Event::fgHist = 0;

//______________________________________________________________________________
McDarkEvent::McDarkEvent() //: fIsValid(kFALSE)
{
   // Create an Event object.
   // When the constructor is invoked for the first time, the class static
   // variable fgTracks is 0 and the TClonesArray fgTracks is created.

   /*if (!fgTracks) fgTracks = new TClonesArray("Track", 1000);
   fTracks = fgTracks;
   fHighPt = new TRefArray;
   fMuons  = new TRefArray;
   fNtrack = 0;
   fH      = 0;
   Int_t i0,i1;
   for (i0 = 0; i0 < 4; i0++) {
      for (i1 = 0; i1 < 4; i1++) {
         fMatrix[i0][i1] = 0.0;
      }
   }
   for (i0 = 0; i0 <10; i0++) fMeasures[i0] = 0;
   for (i0 = 0; i0 <20; i0++) fType[i0] = 0;
   fClosestDistance = 0;
   fEventName = 0;
   fWebHistogram.SetAction(this);*/

	//Set r, p, and E;
/*	x_i = 0; y_i = 0; z_i = 0;
	x_f = 0; y_f = 0; z_f = 0;
	p_x = 0; p_y = 0; p_z = 0;
	E = -1;
	evt = -1;
	id = -1;
	track = -1;
	ccdImages = new TObjArray; */

	Reset();
}

//______________________________________________________________________________
McDarkEvent::~McDarkEvent()
{
   delete ccdImages;
/*   Clear();
   if (fH == fgHist) fgHist = 0;
   delete fH; fH = 0;
   delete fHighPt; fHighPt = 0;
   delete fMuons;  fMuons = 0;
   delete [] fClosestDistance;
   if (fEventName) delete [] fEventName;
*/
}

//______________________________________________________________________________

void McDarkEvent::Build(Int_t EventNo, 
			Double_t newXi, Double_t newYi, Double_t newZi,
			Double_t newXf, Double_t newYf, Double_t newZf,
			Double_t newPx, Double_t newPy, Double_t newPz,
			Double_t newE,
			Int_t newID, Int_t newTrackID, TObjArray* newCCD) {
	SetEvt(EventNo);
	SetX_i(newXi);
	SetY_i(newYi);
	SetZ_i(newZi);
	SetX_f(newXf);
	SetY_f(newYf);
	SetZ_f(newZf);
	SetP_x(newPx);
	SetP_y(newPy);
	SetP_z(newPz);
	SetE(newE);
	SetID(newID);
	SetTrackID(newTrackID);
	SetCCDImages(newCCD);
}

/*
void Event::Build(Int_t ev, Int_t arg5, Float_t ptmin) {
  fIsValid = kTRUE;
  char etype[20];
  Float_t sigmat, sigmas;
  gRandom->Rannor(sigmat,sigmas);
  Int_t ntrack   = Int_t(arg5 +arg5*sigmat/120.);
  Float_t random = gRandom->Rndm(1);

  //Save current Object count
  Int_t ObjectNumber = TProcessID::GetObjectCount();
  Clear();
  fHighPt->Delete();
  fMuons->Delete();
  
  Int_t nch = 15;
  if (ev >= 100)   nch += 3;
  if (ev >= 10000) nch += 3;
  if (fEventName) delete [] fEventName;
  fEventName = new char[nch];
  sprintf(fEventName,"Event%d_Run%d",ev,200);
  sprintf(etype,"type%d",ev%5);
  SetType(etype);
  SetHeader(ev, 200, 960312, random);
  SetNseg(Int_t(10*ntrack+20*sigmas));
  SetNvertex(Int_t(1+20*gRandom->Rndm()));
  SetFlag(UInt_t(random+0.5));
  SetTemperature(random+20.);

  for(UChar_t m = 0; m < 10; m++) {
     SetMeasure(m, Int_t(gRandom->Gaus(m,m+1)));
  }
  for(UChar_t i0 = 0; i0 < 4; i0++) {
    for(UChar_t i1 = 0; i1 < 4; i1++) {
       SetMatrix(i0,i1,gRandom->Gaus(i0*i1,1));
    }
  }

  fTriggerBits.SetBitNumber((UInt_t)(64*gRandom->Rndm(1)));
  fTriggerBits.SetBitNumber((UInt_t)(64*gRandom->Rndm(1)));
  fTriggerBits.SetBitNumber((UInt_t)(64*gRandom->Rndm(1)));

  //  Create and Fill the Track objects
  for (Int_t t = 0; t < ntrack; t++) AddTrack(random,ptmin);
  
  //Restore Object count 
  //To save space in the table keeping track of all referenced objects
  //we assume that our events do not address each other. We reset the 
  //object count to what it was at the beginning of the event.
  TProcessID::SetObjectCount(ObjectNumber);
} */ 

//______________________________________________________________________________
/*Track *Event::AddTrack(Float_t random, Float_t ptmin)
{
   // Add a new track to the list of tracks for this event.
   // To avoid calling the very time consuming operator new for each track,
   // the standard but not well know C++ operator "new with placement"
   // is called. If tracks[i] is 0, a new Track object will be created
   // otherwise the previous Track[i] will be overwritten.

   TClonesArray &tracks = *fTracks;
   Track *track = new(tracks[fNtrack++]) Track(random);
   //Save reference to last Track in the collection of Tracks
   fLastTrack = track;
   //Save reference in fHighPt if track is a high Pt track
   if (track->GetPt() > ptmin)   fHighPt->Add(track);
   //Save reference in fMuons if track is a muon candidate
   if (track->GetMass2() < 0.11) fMuons->Add(track);
   return track;
}*/

//______________________________________________________________________________

//void Event::Clear(Option_t * /*option*/)
//{

//   ccdImages->Clear("C"); //will also call Track::Clear

   // fTracks->Clear("C") 
   // fHighPt->Delete();
   // fMuons->Delete();
   // fTriggerBits.Clear();
//}

//______________________________________________________________________________
void McDarkEvent::Reset()
{
	x_i = 0; y_i = 0; z_i = 0;
	x_f = 0; y_f = 0; z_f = 0;
	p_x = 0; p_y = 0; p_z = 0;
	E = -1;
	evt = -1;
	id = -1;
	track = -1;
	ccdImages = new TObjArray;
}

//______________________________________________________________________________
/*void Event::SetHeader(Int_t i, Int_t run, Int_t date, Float_t random)
{
   fNtrack = 0;
   fEvtHdr.Set(i, run, date);
   if (!fgHist) fgHist = new TH1F("hstat","Event Histogram",100,0,1);
   fH = fgHist;
   fH->Fill(random);
}*/

//______________________________________________________________________________

/*
void Event::SetMeasure(UChar_t which, Int_t what) {
   if (which<10) fMeasures[which] = what;
}
*/

//______________________________________________________________________________

/*void Event::SetRandomVertex() {
   // This delete is to test the relocation of variable length array
   if (fClosestDistance) delete [] fClosestDistance;
   if (!fNvertex) {
      fClosestDistance = 0;
      return;
   }
   fClosestDistance = new Double32_t[fNvertex];
   for (Int_t k = 0; k < fNvertex; k++ ) {
      fClosestDistance[k] = gRandom->Gaus(1,1);
   }
}*/

//______________________________________________________________________________

/*Track::Track(const Track &orig) : TObject(orig)
{
   // Copy a track object

   fPx = orig.fPx;
   fPy = orig.fPy;
   fPz = orig.fPx; 
   fRandom = orig.fRandom;
   fMass2 = orig.fMass2;
   fBx = orig.fBx;
   fBy = orig.fBy;
   fMeanCharge = orig.fMeanCharge;
   fXfirst = orig.fXfirst;
   fXlast  = orig.fXlast;
   fYfirst = orig.fYfirst;
   fYlast  = orig.fYlast;
   fZfirst = orig.fZfirst;
   fZlast  = orig.fZlast;
   fCharge = orig.fCharge;

   fVertex[0] = orig.fVertex[0];
   fVertex[1] = orig.fVertex[1];
   fVertex[2] = orig.fVertex[2];
   fNpoint = orig.fNpoint;
   fNsp = orig.fNsp;
   if (fNsp) {
      fPointValue = new Double32_t[fNsp];
      for(int i=0; i<fNsp; i++) {
         fPointValue[i] = orig.fPointValue[i];
      }
   } else {
      fPointValue = 0;
   }
   fValid  = orig.fValid;

   fTriggerBits = orig.fTriggerBits;

}

//______________________________________________________________________________
Track::Track(Float_t random) : TObject(),fTriggerBits(64)
{
   // Create a track object.
   // Note that in this example, data members do not have any physical meaning.

   Float_t a,b,px,py;
   gRandom->Rannor(px,py);
   fPx = px;
   fPy = py;
   fPz = TMath::Sqrt(px*px+py*py);
   fRandom = 1000*random;
   if (fRandom < 10) fMass2 = 0.106;
   else if (fRandom < 100) fMass2 = 0.8;
   else if (fRandom < 500) fMass2 = 4.5;
   else if (fRandom < 900) fMass2 = 8.9;
   else  fMass2 = 9.8;
   gRandom->Rannor(a,b);
   fBx = 0.1*a;
   fBy = 0.1*b;
   fMeanCharge = 0.01*gRandom->Rndm(1);
   gRandom->Rannor(a,b);
   fXfirst = a*10;
   fXlast  = b*10;
   gRandom->Rannor(a,b);
   fYfirst = a*12;
   fYlast  = b*16;
   gRandom->Rannor(a,b);
   fZfirst = 50 + 5*a;
   fZlast  = 200 + 10*b;
   fCharge = Double32_t(Int_t(3*gRandom->Rndm(1)) - 1);

   fTriggerBits.SetBitNumber((UInt_t)(64*gRandom->Rndm(1)));
   fTriggerBits.SetBitNumber((UInt_t)(64*gRandom->Rndm(1)));
   fTriggerBits.SetBitNumber((UInt_t)(64*gRandom->Rndm(1)));

   fVertex[0] = gRandom->Gaus(0,0.1);
   fVertex[1] = gRandom->Gaus(0,0.2);
   fVertex[2] = gRandom->Gaus(0,10);
   fNpoint = Int_t(60+10*gRandom->Rndm(1));
   fNsp = Int_t(3*gRandom->Rndm(1));
   if (fNsp) {
      fPointValue = new Double32_t[fNsp];
      for(int i=0; i<fNsp; i++) {
         fPointValue[i] = i+1;
      }
   } else {
      fPointValue = 0;
   }
   fValid  = Int_t(0.6+gRandom->Rndm(1));
}
*/
//______________________________________________________________________________
//void Track::Clear(Option_t * /*option*/)
//{
  // fTriggerBits.Clear(); 
  // delete [] fPointValue; 
  // fPointValue=0; 
//}

//______________________________________________________________________________

/*
HistogramManager::HistogramManager(TDirectory *dir)
{
   // Create histogram manager object. Histograms will be created
   // in the "dir" directory.

   // Save current directory and cd to "dir".
   TDirectory *saved = gDirectory;
   dir->cd();

   fNtrack      = new TH1F("hNtrack",    "Ntrack",100,575,625);
   fNseg        = new TH1F("hNseg",      "Nseg",100,5800,6200);
   fTemperature = new TH1F("hTemperature","Temperature",100,19.5,20.5);
   fPx          = new TH1F("hPx",        "Px",100,-4,4);
   fPy          = new TH1F("hPy",        "Py",100,-4,4);
   fPz          = new TH1F("hPz",        "Pz",100,0,5);
   fRandom      = new TH1F("hRandom",    "Random",100,0,1000);
   fMass2       = new TH1F("hMass2",     "Mass2",100,0,12);
   fBx          = new TH1F("hBx",        "Bx",100,-0.5,0.5);
   fBy          = new TH1F("hBy",        "By",100,-0.5,0.5);
   fMeanCharge  = new TH1F("hMeanCharge","MeanCharge",100,0,0.01);
   fXfirst      = new TH1F("hXfirst",    "Xfirst",100,-40,40);
   fXlast       = new TH1F("hXlast",     "Xlast",100,-40,40);
   fYfirst      = new TH1F("hYfirst",    "Yfirst",100,-40,40);
   fYlast       = new TH1F("hYlast",     "Ylast",100,-40,40);
   fZfirst      = new TH1F("hZfirst",    "Zfirst",100,0,80);
   fZlast       = new TH1F("hZlast",     "Zlast",100,0,250);
   fCharge      = new TH1F("hCharge",    "Charge",100,-1.5,1.5);
   fNpoint      = new TH1F("hNpoint",    "Npoint",100,50,80);
   fValid       = new TH1F("hValid",     "Valid",100,0,1.2);

   // cd back to original directory
   saved->cd();
}

//______________________________________________________________________________
HistogramManager::~HistogramManager()
{
   // Clean up all histograms.

   // Nothing to do. Histograms will be deleted when the directory
   // in which tey are stored is closed.
}

//______________________________________________________________________________
void HistogramManager::Hfill(Event *event)
{
   // Fill histograms.

   fNtrack->Fill(event->GetNtrack());
   fNseg->Fill(event->GetNseg());
   fTemperature->Fill(event->GetTemperature());

   for (Int_t itrack = 0; itrack < event->GetNtrack(); itrack++) {
      Track *track = (Track*)event->GetTracks()->UncheckedAt(itrack);
      fPx->Fill(track->GetPx());
      fPy->Fill(track->GetPy());
      fPz->Fill(track->GetPz());
      fRandom->Fill(track->GetRandom());
      fMass2->Fill(track->GetMass2());
      fBx->Fill(track->GetBx());
      fBy->Fill(track->GetBy());
      fMeanCharge->Fill(track->GetMeanCharge());
      fXfirst->Fill(track->GetXfirst());
      fXlast->Fill(track->GetXlast());
      fYfirst->Fill(track->GetYfirst());
      fYlast->Fill(track->GetYlast());
      fZfirst->Fill(track->GetZfirst());
      fZlast->Fill(track->GetZlast());
      fCharge->Fill(track->GetCharge());
      fNpoint->Fill(track->GetNpoint());
      fValid->Fill(track->GetValid());
   }
}
*/
