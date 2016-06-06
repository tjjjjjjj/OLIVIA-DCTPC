// Note that the 'skimmed' runs need to exist before this macro will run

//The official list of physics runs will soon be made available. For now, it is ok to assume that all Big DCTPC runs >=955 are physics runs. The Big DCTPC logbook can be viewed at ($BIGDCTPC_RUNLOG): /net/hisrv0001/home/spitzj/runlog

// To run over runs 2100-2105 (for example), inclusive: 
// root
// .L analysis_125CF4_bigdctpc_alpha.C
// analysis_125CF4_bigdctpc_alpha(2100,2105) 
  
// Alternatively, this program can be run with the run_analysis_alpha.C script:
// root .x run_analysis_alpha.C

// This program's output is a file called outtree_firstrun_lastrun.root which contains a TTree of all the events that pass the cuts. 
// All of the cuts that are utilized are explained in arXiv:1108.4894

void analysis_125CF4_bigdctpc_alpha(int firstrun, int lastrun){//inclusive 
  
  //It is necessary to employ calibration constants in creating the reduced file because we only associate one waveform to one CCD track per event. The matching is based off of energy reconstruction agreement. 	
  double lightCalib = .48;
  double anodeCalib = 19.37;
  double meshCalib = 19.37;
  
  gSystem->Load("libMaxCam");
  gSystem->Load("libWaveformTools.so");

  DmtpcSkimDataset d;
  std::stringstream sstm;
  int passTrack = 0;
  int passTrig = 0;
  int passAll = 0;
  int ntotalTrack = 0;
  int ntrig=0;
  int ntotalTrig = 0;
  int lastSpark = 0; 
  int cutnum=0;
  int spark, edge, burnin, last_spark, runnum, runnum2, evnum, neighbors, Tracknpixel;
  int seqnum=-1;
  int seqnum2=-1;
  int seqcounter=0;
  int setnum=-1;
  int sernum=-1;
  int timenow,time_start;
  double xx, yy, diff, E1, E2, E3, phi;
  double range_ccd=-1.;
  double range_ccd_diffused=-1.;
  double Etrack=-1.;
  double Etrig=-1.;
  double Emesh=-1.;
  double TrackX=-1.;
  double TrackY=-1.;
  double Trackskewness=-1.;
  double Trackmaxpixel=-1.;
  double bestDiff=-1.;
  double phi=-10.;
  double rr=-1.;
  double anodeRMS=0.;
  double meshRMS=0.;
  double cluster_mean=0.;
  double cluster_rms=0.;
  int exposure=0;
  double Trackrms=0.;
  double mesh_peak=0.;
  double mesh_fastpeak=0.;
  double mesh_R0=0.;
  double mesh_R10=0.;
  double mesh_R25=0.;
  double mesh_R50=0.;
  double mesh_R75=0.;
  double mesh_R90=0.;
  double mesh_peaktime=0.;

  TFile *histfile=new TFile(Form("outtree_%d_%d.root",firstrun,lastrun), "RECREATE");
  tree = new TTree("dctpc_eventinfo", "Event info");
  tree->Branch("RunNum", &runnum, "runnum/I");
  tree->Branch("SetNum", &setnum, "setnum/I");
  tree->Branch("SequenceNum", &seqnum, "seqnum/I");
  tree->Branch("SeriesNum", &sernum, "sernum/I");
  tree->Branch("EventNum", &evnum, "evnum/I");
  tree->Branch("Edge", &edge, "edge/I");
  tree->Branch("BurnIn", &burnin, "burnin/I");
  tree->Branch("LastSpark", &last_spark, "last_spark/I");
//   tree->Branch("PassTrack", &passTrack, "passTrack/I");
//   tree->Branch("PassTrig", &passTrig, "passTrig/I");
//   tree->Branch("PassAll", &passAll, "passAll/I");
  tree->Branch("Ntrig", &ntrig, "ntrig/I");
  tree->Branch("Etrack_kev", &Etrack, "Etrack/D");
  tree->Branch("Etrig_kev", &Etrig, "Etrig/D");
  tree->Branch("Emesh_kev", &Emesh, "Emesh/D"); 
  tree->Branch("Cluster_mean_adu", &cluster_mean, "cluster_mean/D");
  tree->Branch("Cluster_rms_adu", &cluster_rms, "cluster_rms/D");
  tree->Branch("Neighbors", &neighbors, "neighbors/I");
  tree->Branch("Phi_deg", &phi, "phi/D"); 
  tree->Branch("TrackX_pix", &TrackX, "TrackX/D"); 
  tree->Branch("TrackY_pix", &TrackY, "TrackY/D");
  tree->Branch("TrackRMS_pix", &Trackrms, "Trackrms/D");
  tree->Branch("Track_maxpixel_adu", &Trackmaxpixel, "Trackmaxpixel/D");
  tree->Branch("Track_npixel", &Tracknpixel, "Tracknpixel/I");
  tree->Branch("Track_skewness", &Trackskewness, "Trackmskewness/D");
  tree->Branch("Rangeccd_pix", &range_ccd, "range_ccd/D");
  tree->Branch("Rangeccd_diffused_pix", &range_ccd_diffused, "range_ccd_diffused/D");
//   tree->Branch("CutNum", &cutnum, "cutnum/I");
  tree->Branch("AnodeRMS_adu", &anodeRMS, "anodeRMS/D");
  tree->Branch("MeshRMS_adu", &meshRMS, "meshRMS/D");
  tree->Branch("Mesh_peak_adu", &mesh_peak, "mesh_peak/D");
  tree->Branch("Mesh_fastpeak_adu", &mesh_fastpeak, "mesh_fastpeak/D");
  tree->Branch("Mesh_R0time_samp", &mesh_R0, "mesh_R0/D");
  tree->Branch("Mesh_R10time_samp", &mesh_R10, "mesh_R10/D");
  tree->Branch("Mesh_R25time_samp", &mesh_R25, "mesh_R25/D");
  tree->Branch("Mesh_R50time_samp", &mesh_R50, "mesh_R50/D");
  tree->Branch("Mesh_R75time_samp", &mesh_R75, "mesh_R75/D");
  tree->Branch("Mesh_R90time_samp", &mesh_R90, "mesh_R90/D");
  tree->Branch("Mesh_peaktime_samp", &mesh_peaktime, "mesh_peaktime/D");
  tree->Branch("Timenow_sec",&timenow, "timenow/I");
  tree2 = new TTree("dctpc_runinfo", "Run info");
  tree2->Branch("RunNum", &runnum, "runnum2/I");
  tree2->Branch("Exposure_sec", &exposure, "exposure/I");
  tree2->Branch("Time_startofrun_sec", &time_start, "time_start/I");
  
//start run loop
for (int x = firstrun; x <= lastrun; x++)
{
      runnum=x;
           
      if(runnum>=889&&runnum<=891)
      seqnum=0;
      
      if(runnum>=893&&runnum<=896)
      seqnum=1;
      
      if(runnum>=897&&runnum<=904)
      seqnum=2;
      
      if(runnum>=910&&runnum<=917)
      seqnum=3;
      
      if(runnum>=919&&runnum<=925)
      seqnum=4;
      
      if(runnum>=926&&runnum<=938)
      seqnum=5;
      
      if(runnum>=939&&runnum<=954)
      seqnum=6;    
      
      if(runnum>=955&&runnum<=1645)
      seqnum=7;
      
      if(runnum>=1646&&runnum<=2205)
      seqnum=8;
      
      if(runnum>=2206&&runnum<=2648)
      seqnum=9; 
      
      if(runnum>=2650&&runnum<=3868)
      seqnum=10; 
        
      if(runnum>=3876)
      seqnum=11; 
      
              
      if(runnum>=889&&runnum<=954)
      setnum=0;
      
      if(runnum>=955&&runnum<=3868)
      setnum=1;
      
      if(runnum>=3876)
      setnum=2;
      
      
      if(seqnum!=seqnum2)
      seqcounter=0;
      
      seqnum2=seqnum;
    
      
      string origfile = "/net/nudsk0001/d00/scratch/dctpc_tmp/bigdctpc_data/BigDCTPC_run_";
      string skimfile = "/net/nudsk0001/d00/scratch/dctpc_tmp/bigdctpc_skim/BigDCTPC_run_";
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
      if(!ifile)
	  continue;
    
      DmtpcSkimDataset d;
      d.openRootFile(skimfilename.c_str());     
      d.loadDmtpcEvent(true,origfilename.c_str());
      Long_t nEvents=d.tree()->GetEntriesFast();
	  runnum2=x; 
	  exposure=0;	
	  lastSpark=0;

if(seqcounter%10==0)
sernum++;

seqcounter++;

//start event loop 
for(int i=0; i<nEvents; i++){    

if(i%100==0)
cout<<"Event: "<<i<<" / "<<nEvents<<endl;
    
evnum=i;	
DmtpcSkimEvent* ev=d.event(); //creating an event object		
d.getEvent(i);

if(i==0)
time_start=ev->timenow(0);

	if(ev->spark(0)){lastSpark = ev->eventNumber(); continue;} 
	
	for(int k=0;k<ev->ntracks(0);k++)
	  {	
	  if (ev->maxpixel(0,k)>600){lastSpark = ev->eventNumber();continue;}//maxpixel has to be less than 600ADU. also, the event might have been a spark!
	  }
	if(evnum-lastSpark<7){continue;} // make sure that at least 6 events have passed after a spark

    exposure++;//increase the integrated exposure as long as event is not a spark or close to a spark in time

	if(ev->ntracks(0)==0){continue;}		

	Etrack=0.;
	Etrig=0.;
	Emesh=0.;
	passTrack=0;
	passTrig=0;
	passAll=0;
	timenow = ev->timenow(0);
	bestDiff=1e9;

	//loop over tracks in each event	
	  for(int k=0;k<ev->ntracks(0);k++)
	  {	    
	    if (ev->E(0,k)<1){continue;} //Cut if E is below 1 (measured in ADU)
	    //if (ev->nburnin(0,k)){cutnum = 41;  continue;} //Cut if there is a burnin
	    // if (ev->edge(0,k)){continue;} //Cut if track reaches the edge
	    if (ev->range(0,k)<=0){continue;} //Cut if it has a zero range
// 	    if (ev->x(0,k)<48) continue; //this and the following->
// 	    if (ev->y(0,k)<48) continue; //instructions make sure that->
// 	    if (ev->x(0,k)>976) continue; //the track does not cross->
// 	    if (ev->y(0,k)>976)continue; //the edge
	    if (ev->maxpixel(0,k)/ev->E(0,k)>0.25){continue;} //To eliminate noise artifacts the max pixel has to be less than 25% of the total signal      
           
	    passTrig=0;
	    passTrack++;

	    if(passTrack==1)
	    {
		E1 = ev->E(0,k)*lightCalib;
		Etrack = E1;
		phi = ev->phi(0,k);
		phi = atan2(sin(phi),cos(phi)) * 180./ 3.1416;   
		range_ccd=ev->range(0,k);
		range_ccd_diffused=ev->diffusedRange(0,k);		         	      	    
	    TrackX = ev->x(0,k)-512;
	    TrackY = ev->y(0,k)-512; 
        Trackskewness=ev->skewness(0,k);
        edge=ev->edge(0,k);
        burnin=ev->nburnin(0,k);
        last_spark=evnum-lastSpark;
        cluster_mean=ev->cluster_mean(0,k);
        cluster_rms=ev->cluster_rms(0,k);
        neighbors=ev->neighbors(0,k);
        Trackmaxpixel=ev->maxpixel(0,k);
        Tracknpixel=ev->npixel(0,k);
        Trackrms=sqrt(ev->transverse_moment(0,2,k));      
	    }
	    else continue;

if(d.orig_event()->scopeDataInfo(0)==NULL){continue;}//no waveforms!!
    TObjArray* arr = ev->waveform_vectors();
    CspWfVector* anode = (CspWfVector*) arr->At(0);
	FastWfVector* mesh = (FastWfVector*) arr->At(1);
	CspWfVector* veto = (CspWfVector*) arr->At(2); 

	    //loop over triggers
	    ntrig=mesh->size();
	    for (int nt = 0; nt < mesh->size(); nt++)
	    {
	    diff=1e9;
	    meshRMS=mesh->at(nt).getRMS();
	    anodeRMS=anode->at(nt).getRMS();

	    passTrig++;//increment number of tracks that pass the trigger cuts
	    E2 = (anode->at(nt,0).getPeak()*1000.)*anodeCalib;
	    diff = E1 - E2;
	    
	    	if (fabs(diff)<fabs(bestDiff)) 
	    	{	    	
			bestDiff = diff;
			Etrig = E2;			
			mesh_peak=mesh->at(nt,0).getPeak();
			mesh_fastpeak=mesh->at(nt,0).getFastPeak();
			mesh_R0=mesh->at(nt,0).getFastRise0();
			mesh_R10=mesh->at(nt,0).getFastRise10();
			mesh_R25=mesh->at(nt,0).getFastRise25();
			mesh_R50=mesh->at(nt,0).getFastRise50();
			mesh_R75=mesh->at(nt,0).getFastRise75();
			mesh_R90=mesh->at(nt,0).getFastRise90();
			mesh_peaktime=mesh->at(nt,0).getFastPeakTime();
			Emesh=mesh->at(nt,0).getIntegral()*meshCalib;
	    	}
	    	
	    	
	    } // end trigger loop 
	    ntotalTrig+=passTrig;
	  }//end track loop

	ntotalTrack+=passTrack;
	
	if (passTrack!=1){continue;} // events are unlikely to have more than 1 track
	if (passTrig==0){continue;} 

	passAll=1;
		
	tree->Fill();	
	
}//end event loop

tree2->Fill();
cout <<"Total seconds of exposure: "<<exposure<<endl;
}//end file loop



histfile->Write();
histfile->Close();
gApplication->Terminate(); 
}
