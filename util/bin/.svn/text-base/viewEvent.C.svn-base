// run in batch mode!! 
// root -l -b viewEvent.C 

// this macro creates a pdf containing the CCD image and associated waveforms of events specified in events.dat 
{
  string DETECTOR="BIGDCTPC";//BIGDCPTC or LITTLEDCTPC
  string HIGHRES="YES"; //Create high res CCD image (which may be slower on older GPUs)
  
  double lightCalib = 0.48; 
  double anodeCalib= 23.2;
  double meshCalib = 37.8;
      
  gStyle->SetPalette(1,0); 
  gSystem->Load("libMaxCam");
  gSystem->Load("libWaveformTools.so");
  gStyle->SetOptStat(0);
  double stops[5] = {0,0.34,0.61,0.84,1.0};
  double red[5] = {0.0,0.0,0.87,1.0,0.51};
  double green[5] = {0.0,0.81,1.00,0.2,0.0};
  double blue[5] = {0.51,1.0,0.12,0.0,0.0};
  TColor::CreateGradientColorTable(5,stops,red,green,blue,255);
  gStyle->SetNumberContours(255);

  std::stringstream sstm;
  string skimend = "skim.root";
  string origend = ".root";
  string origfilename;
  string skimfilename;
  vector<double> noiseVector;
  DmtpcSkimDataset d;

TRandom3* r = new TRandom3();

int run,event;
int lineno=0;
int totallines=0;
string line;
ifstream infile; 

int index=-1;
double bestDiff=1e9;
double diff;

TLatex *tex1 = new TLatex();
tex1->SetTextSize(0.03);
tex1->SetTextColor(1);
tex1->SetNDC();

TLatex *tex2 = new TLatex();
tex2->SetTextSize(0.05);
tex2->SetTextColor(1);
tex2->SetNDC();

TLatex *tex3 = new TLatex();
tex3->SetTextSize(0.07);
tex3->SetTextColor(1);
tex3->SetNDC();

TLine l;
l.SetLineStyle(1);
double phi,TrackX, TrackY,range_ccd,TrackXStart,TrackYStart,TrackXEnd,TrackYEnd,Trackrms;
//count number of lines first
infile.open("events.dat");
    while ( getline (infile,line) )
    {    	
    	istringstream ss(line);
        totallines++;
    }   
infile.close();   

infile.open("events.dat");
  if (infile.is_open())
  {   
    while ( getline (infile,line) )
    {    	
    	istringstream ss(line);
        ss >> run >> event;
   lineno++;

   if(run==0)
   continue;

   cout <<"Run: "<< run <<" Event: "<<event<<endl;
   sstm.str("");
   string origfile = "/net/nudsk0001/d00/scratch/dctpc_tmp/bigdctpc_data/BigDCTPC_run_";
   string skimfile = "/net/nudsk0001/d00/scratch/dctpc_tmp/bigdctpc_skim/BigDCTPC_run_";
   
   if(DETECTOR=="LITTLEDCTPC")
   {
   origfile = "/net/nudsk0001/d00/scratch/dctpc_tmp/data/dmtpc_DC_";
   skimfile = "/net/nudsk0001/d00/scratch/dctpc_tmp/skim/dmtpc_DC_";
   lightCalib = .15;
   anodeCalib = 6.5;
   meshCalib = 32.;
   }
  
   
   if (run<10000){ origfile+="0"; skimfile+="0"; }
   if (run<1000){ origfile+="0"; skimfile+="0"; }
   if (run<100){ origfile+="0"; skimfile+="0"; }
   if (run<10){ origfile+="0"; skimfile+="0"; }
   sstm << origfile << run << origend;
   origfilename = sstm.str();
   sstm.str("");
   sstm << skimfile << run << skimend;
   skimfilename = sstm.str();
//    cout << skimfilename << endl;
//    cout << origfilename << endl;
   ifstream ifile(origfilename.c_str());
   if(!ifile)
   {
   cout << origfilename <<" does not exist; skipping"<<endl;
	continue;
}
  // Make DmtpcSkimDataset class instance (data holder)
  DmtpcSkimDataset d;
  // Load skim data file ... this contains reconstructed information & ccd image
  d.openRootFile(skimfilename.c_str());
  // Load original data file ... this is used as an input to reconstruction, and called "original" ... contains charge read-out waveforms
  d.loadDmtpcEvent(true,origfilename.c_str());
  DmtpcSkimEvent* ev=d.event(); 
  // Let's see how many events we have in this file.
  Long_t nEvents=d.tree()->GetEntriesFast();

    d.getEvent(event);
    TCanvas *c=new TCanvas("c","",600,500);
    c->SetGrid();
    TH2F* hImage = (TH2F*)(d.event()->cluster(0)->getImage());// access ccd
    //rebin to reduce the size of the output file!!
    
    if(HIGHRES=="NO")
    {
    hImage->Rebin2D(2,2);
    hImage->SetMinimum(-100);
    hImage->SetMaximum(500); 
    }
    else
    {
    hImage->SetMinimum(-50);
    hImage->SetMaximum(250);
    }  
    
    hImage->SetTitle(Form("Run %d, Event %d;",run,event));
        
    hImage->Draw("COLZ");
    tex1->DrawLatex(.72,.97,Form("Raw E: %4.2f keV",ev->E(0,0)*lightCalib*.56));    
    tex1->DrawLatex(.72,.94,Form("Range: %4.2f pix",ev->range(0,0)));
    tex1->DrawLatex(.72,.91,Form("1#sigma width: %4.2f pix",sqrt(ev->transverse_moment(0,2,0)))); 
    
phi = ev->phi(0,0);
phi = atan2(sin(phi),cos(phi)) * 180./ 3.1416; 
TrackX = ev->x(0,0);
TrackY = ev->y(0,0);
range_ccd=ev->range(0,0); 
// TrackXStart = TrackX + (0.5*range_ccd*cos(phi*TMath::Pi()/180.));
// TrackYStart = TrackY + (0.5*range_ccd*sin(phi*TMath::Pi()/180.));
// TrackXEnd = TrackX - (0.5*range_ccd*cos(phi*TMath::Pi()/180.));
// TrackYEnd = TrackY - (0.5*range_ccd*sin(phi*TMath::Pi()/180.));

TrackXStart = ev->xbegin(0,0);
TrackYStart = ev->ybegin(0,0);
TrackXEnd = ev->xend(0,0);
TrackYEnd = ev->yend(0,0);

//Trackrms=sqrt(ev->transverse_moment(0,2,0));
cout<<ev->npixel(0,0)<<" "<<ev->range(0,0)<<endl;
Trackrms=4.*ev->npixel(0,0)/ev->range(0,0);

TLine ll;
ll.SetLineStyle(1);
ll.SetLineWidth(3);
TLine ll2;
ll2.SetLineColor(2);
ll2.SetLineStyle(1);
ll2.SetLineWidth(3);

double deltax=2.*Trackrms*cos(1.5708+phi*TMath::Pi()/180.);
double deltay=2.*Trackrms*sin(1.5708+phi*TMath::Pi()/180.);
double xstart1=TrackXStart+deltax+(.5*deltax*cos(1.5708+phi*TMath::Pi()/180.));
double ystart1=TrackYStart+deltay+(.5*deltay*sin(1.5708+phi*TMath::Pi()/180.));
double xstart2=TrackXStart-deltax-(.5*deltax*cos(1.5708+phi*TMath::Pi()/180.));
double ystart2=TrackYStart-deltay-(.5*deltay*sin(1.5708+phi*TMath::Pi()/180.));
double xend1=TrackXEnd+deltax+(.5*deltax*cos(1.5708+phi*TMath::Pi()/180.));
double yend1=TrackYEnd+deltay+(.5*deltay*sin(1.5708+phi*TMath::Pi()/180.));
double xend2=TrackXEnd-deltax-(.5*deltax*cos(1.5708+phi*TMath::Pi()/180.));
double yend2=TrackYEnd-deltay-(.5*deltay*sin(1.5708+phi*TMath::Pi()/180.));

// ll.DrawLine(TrackXStart,TrackYStart,TrackXEnd,TrackYEnd);
ll.DrawLine(xstart1,ystart1,xstart2,ystart2);
ll2.DrawLine(xend1,yend1,xend2,yend2);

    c->Print("CCDandWF.pdf(");

    // Access # of waveforms stored
     Int_t nWFs = (Int_t)(d.orig_event()->scopeData()->GetEntries());

     TCanvas *c1=new TCanvas("c","",600,500);
  
  
  
bestDiff=1e9;
     
if(DETECTOR=="BIGDCTPC")    
{ 
     c1->Divide(1,2);
       
       if(d.orig_event()->scopeDataInfo(0)==NULL){continue;}//no waveforms!!
       
    for(int j=0; j<nWFs/3; j++){
     

      TObjArray* arr = ev->waveform_vectors();
      CspWfVector* anode = (CspWfVector*) arr->At(0);
	  FastWfVector* mesh = (FastWfVector*) arr->At(1);
	  CspWfVector* veto = (CspWfVector*) arr->At(2); 
    
    
    for(int k=0;k<ev->ntracks(0);k++)
	  {	 
	  if (ev->E(0,k)<1){continue;} //Cut if E is below 1 (measured in ADU)
	  if (ev->range(0,k)<=0){continue;} //Cut if it has a zero range
	  if (ev->maxpixel(0,k)/ev->E(0,k)>0.25){continue;}
	  
    
      double E1 = ev->E(0,k)*lightCalib;
      double E2 = (anode->at(j,0).getPeak()*1000.)*anodeCalib;
      diff=TMath::Abs( (E1 - E2)/(E1 + E2) );

    if(fabs(diff)<fabs(bestDiff)) 
	{  
	bestDiff=diff;
    index=j;
    }
  
      }
  
	}

	  TH1F* hMesh  = (TH1F*)d.orig_event()->scopeData(1,index);
      hMesh->SetTitle(Form("Run %d, Event %d Mesh Waveform; Time; Amplitude",run,event));
      TH1F* hAnode = (TH1F*)d.orig_event()->scopeData(0,index);
      hAnode->SetTitle(Form("Run %d, Event %d Anode Waveform; Time; Amplitude",run,event));
      
      c1->cd(1);
      c1->cd(1)->SetGrid();
      //rebins reduce the size of the output file!!
      hAnode->Rebin(50);      
      hAnode->Draw();

      tex2->DrawLatex(.55,.83,Form("Raw E (amplitude): %4.2f keV",anode->at(index,0).getPeak()*1000.*anodeCalib*0.56));
	  tex2->DrawLatex(.55,.78,Form("Rise time: %4.2f samp",anode->at(index,0).getRise0())); 	

      l.DrawLine(anode->at(index,0).getStartBin(),0,anode->at(index,0).getStartBin(),3000.);
      l.DrawLine(anode->at(index,0).getPeakBin(),0,anode->at(index,0).getPeakBin(),3000.);
             
      c1->cd(2); 
      c1->cd(2)->SetGrid();
      hMesh->Rebin(50); 
      hMesh->Draw();
      
      tex2->DrawLatex(.55,.83,Form("Raw E (area): %4.2f keV",mesh->at(index,0).getIntegral()*meshCalib*0.56));
      tex2->DrawLatex(.55,.78,Form("Rise time: %4.2f samp",mesh->at(index,0).getRise0()));
      tex2->DrawLatex(.55,.73,Form("Total time: %4.2f samp",mesh->at(index,0).getPeakBin()+mesh->at(index,0).getFall10()-mesh->at(index,0).getStartTime()));

      l.DrawLine(mesh->at(index,0).getStartBin(),0,mesh->at(index,0).getStartBin(),3000.);
      l.DrawLine(mesh->at(index,0).getPeakBin(),0,mesh->at(index,0).getPeakBin(),3000.);      	 l.DrawLine(mesh->at(index,0).getPeakBin()+mesh->at(index,0).getFall10(),0,mesh->at(index,0).getPeakBin()+mesh->at(index,0).getFall10(),3000.);
      
      
  //     
       //cout<<mesh->at(index,0).getEndBin()<<" "<<mesh->at(index,0).getPeakBin()+mesh->at(index,0).getFall10()<<endl;
      
       
}


if(DETECTOR=="LITTLEDCTPC")    
{ 
     c1->Divide(1,3);

     if(d.orig_event()->scopeDataInfo(0)==NULL){continue;}//no waveforms!!
       
    for(int j=0; j<nWFs/4; j++){
 

      TObjArray* arr = ev->waveform_vectors();
      FastWfVector* mesh = (FastWfVector*) arr->At(0);
	  CspWfVector* anode = (CspWfVector*) arr->At(1);
	  CspWfVector* veto = (CspWfVector*) arr->At(2);
       
    for(int k=0;k<ev->ntracks(0);k++)
	  {	 
	  
 	  if (ev->E(0,k)<1){continue;} //Cut if E is below 1 (measured in ADU)
	  if (ev->range(0,k)<=0){continue;} //Cut if it has a zero range
	  if (ev->maxpixel(0,k)/ev->E(0,k)>0.25){continue;}
	  
      double E1 = ev->E(0,k)*lightCalib;
      double E2 = (anode->at(j,0).getPeak()*1000.)*anodeCalib;
      diff=TMath::Abs( (E1 - E2)/(E1 + E2) );

    if (fabs(diff)<fabs(bestDiff)) 
	{  	
	bestDiff=diff;
    index=j;
    }

  }
  
	}

	  TH1F* hMesh  = (TH1F*)d.orig_event()->scopeData(0,index);
      hMesh->SetTitle(Form("Run %d, Event %d Mesh Waveform; Time; Amplitude",run,event));
      TH1F* hAnode = (TH1F*)d.orig_event()->scopeData(1,index);
      hAnode->SetTitle(Form("Run %d, Event %d Anode Waveform; Time; Amplitude",run,event));
      TH1F* hVeto = (TH1F*)d.orig_event()->scopeData(2,index);
      hVeto->SetTitle(Form("Run %d, Event %d Veto Waveform; Time; Amplitude",run,event));

      c1->cd(1);
      c1->cd(1)->SetGrid();
      //rebins reduce the size of the output file!!
      hAnode->Rebin(32);      
      hAnode->Draw();
      tex3->DrawLatex(.55,.8,Form("Raw E (amplitude): %4.2f keV",anode->at(index,0).getPeak()*1000.*anodeCalib));
      tex3->DrawLatex(.55,.73,Form("Rise time: %4.2f microsec",anode->at(index,0).getRise0()*1000000.)); 
 
l.DrawLine(hAnode->GetBinCenter(anode->at(index,0).getStartBin()/32.),-3000.,hAnode->GetBinCenter(anode->at(index,0).getStartBin()/32.),3000.); 
 l.DrawLine(hAnode->GetBinCenter(anode->at(index,0).getPeakBin()/32.),-3000.,hAnode->GetBinCenter(anode->at(index,0).getPeakBin()/32.),3000.); 
 




        
      c1->cd(2); 
      c1->cd(2)->SetGrid();
    
      hMesh->Rebin(32); 
      hMesh->Draw();
      tex3->DrawLatex(.55,.8,Form("Raw E (area): %4.2f keV",mesh->at(index,0).getIntegral()*meshCalib));
      tex3->DrawLatex(.55,.73,Form("Rise time: %4.2f microsec",mesh->at(index,0).getRise0()*1000000.));
      tex3->DrawLatex(.55,.66,Form("Total time: %4.2f microsec",((hMesh->GetBinCenter(mesh->at(index,0).getPeakBin()/32.))+(mesh->at(index,0).getFall10())-hMesh->GetBinCenter(mesh->at(index,0).getStartBin()/32.))*1000000.));
      
l.DrawLine(hMesh->GetBinCenter(mesh->at(index,0).getStartBin()/32.),-3000.,hMesh->GetBinCenter(mesh->at(index,0).getStartBin()/32.),3000.); 
l.DrawLine(hMesh->GetBinCenter(mesh->at(index,0).getPeakBin()/32.),-3000.,hMesh->GetBinCenter(mesh->at(index,0).getPeakBin()/32.),3000.);  
//l.DrawLine(hMesh->GetBinCenter(mesh->at(index,0).getEndBin()/32.),-3000.,hMesh->GetBinCenter(mesh->at(index,0).getEndBin()/32.),3000.); 
l.DrawLine((hMesh->GetBinCenter(mesh->at(index,0).getPeakBin()/32.))+(mesh->at(index,0).getFall10()),-3000.,(hMesh->GetBinCenter(mesh->at(index,0).getPeakBin()/32.))+(mesh->at(index,0).getFall10()),3000.); 




// (hMesh->GetBinCenter(mesh->at(index,0).getPeakBin())+mesh->at(index,0).getFall10())/32.


// (mesh->at(index,0).getPeakBin()+mesh->at(index,0).getFall10())


      c1->cd(3); 
      c1->cd(3)->SetGrid();
      hVeto->Rebin(32); 
      hVeto->Draw();
}
   
  if(lineno!=totallines)
  c1->Print("CCDandWF.pdf");
  
  if(lineno==totallines)
  c1->Print("CCDandWF.pdf)");

}
	}
  }  
  
infile.close();

}
