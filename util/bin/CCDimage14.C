void CCDimage14(int firstrun, int lastrun)
{
  gStyle->SetPalette(1,0); 
  gSystem->Load("libMaxCam");
  gSystem->Load("libWaveformTools.so");

  std::stringstream sstm;

  TH2F * sum = 0; 
  TH2I * count = 0;
  TH2F * p=0;
         sum = new TH2F("alpha_sum","alpha_sum",256,0,1024,256,0,1024); 
        count = new TH2I("alpha_count","alpha_count",256,0,1024,256,0,1024);
  const int c = 0;
 
  gROOT->cd(); 

  int runnum;
  int setnum;
  
    if		(firstrun>=889 && lastrun<=954)		{setnum=0;}
    else if	(firstrun>=955 && lastrun<=3868)	{setnum=1;}
    else if	(firstrun>=3876 && lastrun<=5789)	{setnum=2;}
	else if	(firstrun>=5792)					{setnum=3;}
	else
	{
	gApplication->Terminate();
	exit;}
  
  
for (int x = firstrun; x <= lastrun; x++){
	runnum=x;
	
       string origfile = "/net/nudsk0001/d00/scratch/dctpc_tmp/bigdctpc_data/BigDCTPC_run_";
       string skimfile = "/net/nudsk0001/d00/scratch/dctpc_tmp/bigdctpc_skim/BigDCTPC_run_";  
      
//            string origfile = "/net/hisrv0001/home/spitzj/myOutput_";
//      string skimfile = "/net/hisrv0001/home/spitzj//myOutput_";
      
      string skimend = "skim.root";
      string origend = ".root";
      string origfilename;
      string skimfilename;
      sstm.str("");
      if (x<10000){ origfile+="0"; skimfile+="0"; }
      if (x<1000){ origfile+="0"; skimfile+="0"; }
      if (x<100){ origfile+="0"; skimfile+="0"; }
      if (x<10){ origfile+="0"; skimfile+="0"; }
      sstm << origfile << x << origend;
      origfilename = sstm.str();
      sstm.str("");
      sstm << skimfile << x << skimend;
      skimfilename = sstm.str();
      cout << origfilename << endl;
      ifstream ifile(origfilename.c_str());
      ifstream ifile2(skimfilename.c_str());
      if(!ifile)
	  continue;
      if(!ifile2)
	continue;
  	DmtpcSkimDataset d;
  	d.openRootFile(skimfilename.c_str());    
  	d.loadDmtpcEvent(true,origfilename.c_str());
  	for (int i = 0; i < d.nevents(); i ++){
        if(i%100==0){
			cout<<"Event: "<<i<<endl;
		}

	//cout<<"here"<<endl
    d.getEvent(i);    
    //cout<<"here2"<<endl;
    for (int t =0 ; t < d.event()->ntracks(c); t++)
    {
      if(d.event()->maxpixel(0,t)>150){continue;}
      if(d.event()->spark(0)){continue;} 
      p = (TH2F*)d.event()->cluster(c)->getImage(); 
      vector<int> clust = d.event()->cluster(c)->getCluster(t); 
     
    	for (vector<int>::iterator it = clust.begin(); it!=clust.end(); it++)
		{
	  		sum->SetBinContent(*it, sum->GetArray()[*it] + p->GetArray()[*it]);
	  		count->SetBinContent(*it, count->GetArray()[*it] + 1); 
		}    
    }
  }
 }

  TH2F * normalized = (TH2F*) sum->Clone("normalized"); 
  normalized->SetName("normalized");
  normalized->Divide(count); 
  normalized->SetStats(0);
  normalized->Draw("COLZ");
  c1->Print(Form("CCDimage14_set_%d_runs_%d_%d.pdf",setnum,firstrun,lastrun));

  TFile *f=new TFile(Form("CCDimage14_set_%d_runs_%d_%d.root",setnum,firstrun,lastrun),"RECREATE");
  normalized->Write();
  f->Close();

  gApplication->Terminate();
}
