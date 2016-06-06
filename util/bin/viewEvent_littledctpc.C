

{
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
  
//   c->SetRightMargin(0.08);
//   c->SetLeftMargin(0.13);
//   c->SetBottomMargin(0.15);
//   std::string orig_file="/net/nudsk0001/d00/scratch/dctpc_tmp/data/dmtpc_DC_01646.root";
//   std::string skim_file="/net/nudsk0001/d00/scratch/dctpc_tmp/skim/dmtpc_DC_01646skim.root";


 std::stringstream sstm;
  string origfile = "/net/nudsk0001/d00/scratch/dctpc_tmp/bigdctpc_data/BigDCTPC_run_0";
  string skimfile = "/net/nudsk0001/d00/scratch/dctpc_tmp/bigdctpc_skim/BigDCTPC_run_0";
  string skimend = "skim.root";
  string origend = ".root";
  string origfilename;
  string skimfilename;
  vector<double> noiseVector;

  // Make DmtpcSkimDataset class instance (data holder)
  DmtpcSkimDataset d;

int run,event;
int lineno=0;
int totallines=0;
string line;
ifstream infile; 

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

   cout <<"Run: "<< run <<" Event: "<<event<<endl;
   sstm.str("");
   sstm << origfile << run << origend;
   origfilename = sstm.str();
   sstm.str("");
   sstm << skimfile << run << skimend;
   skimfilename = sstm.str();
   cout << skimfilename << endl;
   cout << origfilename << endl;


  // Make DmtpcSkimDataset class instance (data holder)
  DmtpcSkimDataset d;
  // Load skim data file ... this contains reconstructed information & ccd image
  d.openRootFile(skimfilename.c_str());
  // Load original data file ... this is used as an input to reconstruction, and called "original" ... contains charge read-out waveforms
  d.loadDmtpcEvent(true,origfilename.c_str());

  // Let's see how many events we have in this file.
  Long_t nEvents=d.tree()->GetEntriesFast();
  std::cout
    << Form("Found %d events loaded...", nEvents)
    << std::endl
    << std::endl;

  // Let's loop over first 20 events.
  std::cout
    << "Start an event loop!"
    << std::endl;



    d.getEvent(event);
    TCanvas *c=new TCanvas("c","",600,500);
    // Access ccd image (TH2F) & store image!
    TH2F* hImage = (TH2F*)(d.event()->cluster(0)->getImage());// access ccd
    hImage->SetTitle(Form("Run %d, Event %d;",run,event));

    //new TCanvas;
    hImage->SetMinimum(-50);
    hImage->SetMaximum(250);
    hImage->Draw("COLZ");
    c->Print("CCDandWF.pdf(");

    // Access # of waveforms stored
     Int_t nWFs = (Int_t)(d.orig_event()->scopeData()->GetEntries());
     if(nWFs)
     std::cout
       << Form(" found %d waveforms (=%d times charge read-out triggered)", nWFs, nWFs/3)       << std::endl;


TCanvas *c1=new TCanvas("c1","",600,500);
c1->Divide(nWFs/3,3);

    // Access waveforms (TH1F) & store image!
    for(int j=0; j<nWFs/3; j++){

      TH1F* hMesh  = (TH1F*)d.orig_event()->scopeData(1,j);
      hMesh->SetTitle(Form("Run %d, Event %d Mesh Waveform %d/%d; Time; Amplitude",run,event,j,nWFs/3));
      TH1F* hAnode = (TH1F*)d.orig_event()->scopeData(0,j);
      hAnode->SetTitle(Form("Run %d, Event %d Anode Waveform %d/%d; Time; Amplitude",run,event,j,nWFs/3));
      TH1F* hVeto = (TH1F*)d.orig_event()->scopeData(2,j);
      hVeto->SetTitle(Form("Run %d, Event %d Veto Waveform %d/%d; Time; Amplitude",run,event,j,nWFs/3));
      
      c1->cd(j+1);          
      hAnode->Draw();
      c1->cd(j+1+(nWFs/3));
      hMesh->Draw();
      c1->cd(j+1+(2*(nWFs/3)));
      hVeto->Draw();

// 
//       new TCanvas;
//       hMesh->Draw();
//       c->SaveAs("mesh_waveform.gif+50");
//       new TCanvas;
//       hAnode->Draw();
//       c->SaveAs("anode_waveform.gif+50");

  }
  if(lineno!=totallines)
  c1->Print("CCDandWF.pdf");
  else
  c1->Print("CCDandWF.pdf)");
 
 
 


	}
  }
  
  
infile.close();


}
