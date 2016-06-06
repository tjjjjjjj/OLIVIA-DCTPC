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
#include <iomanip>;
#include <sys/stat.h>;
void analysis_125CF4_bigdctpc_alpha(int firstrun, int lastrun){//inclusive 
  
  //It is necessary to employ calibration constants in creating the reduced file because we only associate one waveform to one CCD track per event. The matching is based off of energy reconstruction agreement. 	

  double lightCalib = 0.48;
  double anodeCalib= 23.2;
  double meshCalib = 37.8;
  const int maxentries=100;
  
  int BATCH=0;//0 if running interactively, 1 if running batch
  
  gSystem->Load("libMaxCam");
  gSystem->Load("libWaveformTools.so");

  DmtpcSkimDataset d;
  DmtpcSkimDataset d2;
  std::stringstream sstm;
  std::stringstream sstmxxx;
  int passTrack = 0;
  int passTrig = 0;
  int passAll = 0;
  int ntotalTrack = 0;
  int ntrack=0;
  int ntrig=0;
  int ntotalTrig = 0;
  int lastSpark = -1; 
  int cutnum=0;
  int spark, edge, burnin, last_spark, runnum, runnum2, evnum, neighbors, Tracknpixel;
  double image_mean=0.;
  double image_rms=0.;
  int pixels_killed=0;
  int seqnum=-1;
  int seqnum2=-1;
  int seqcounter=0;
  int setnum=-1;
  int timenow,time_start;
  int next_spark=10000;
  
  double pressure=0.;
  double voltage_amp=0.;
  double voltage_drift=6000.;
  double xx, yy, diff, E1, E2, E3, phi;
  double range_ccd=-1.;
  double Etrack=-1.;
  double Etrig=-1.;
  double Emesh=-1.;
  double TrackX=-10000.;
  double TrackY=-10000.;
  double TrackXStart=-10000.;
  double TrackYStart=-10000.;
  double TrackXEnd=-10000.;
  double TrackYEnd=-10000.;
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
  int expose=0;
  double Trackrms=0.;
  double Trackwidth=0.;
  double mesh_peak=0.;
  double mesh_start=0.;	
  double anode_start=0.;
  double mesh_base=0.;
  double anode_base=0.; 
  double mesh_max;
  double anode_max;
  double veto_peak=0.;
  double mesh_width=0.;
  double anode_R0=0.;
  double mesh_R0=0.;
  double mesh_R10=0.;
  double mesh_R25=0.;
  double mesh_R50=0.;
  double mesh_R75=0.;
  double mesh_R90=0.;
  double mesh_F0=0.;
  double mesh_F10=0.;
  double mesh_F25=0.;
  double mesh_F50=0.;
  double mesh_F75=0.;
  double mesh_F90=0.;
  double mesh_peaktime=0.;
  double triggertimestamp[maxentries];
  int triggerindex;
  int badrun; 
  int badrunflag;
  int totaltrack=0;
  int totaltrig=0;
  double totalmesh_allwf=0.;
  double totalanode_allwf=0.;


long int orig_filesize=-1;
long int new_filesize=0;
long int orig_filesize2=-1;
long int new_filesize2=0;

      string origfile;
      string skimfile;

      string outfile1;
      std::stringstream out1;
      out1.str("");
      if (firstrun<10000){ outfile1+="0";}
      if (firstrun<1000){ outfile1+="0"; }
      if (firstrun<100){ outfile1+="0";}
      if (firstrun<10){ outfile1+="0"; }
      out1 << outfile1 << firstrun ;
  
      string outfile2;
      std::stringstream out2;
      out2.str("");
      if (lastrun<10000){ outfile2+="0";}
      if (lastrun<1000){ outfile2+="0"; }
      if (lastrun<100){ outfile2+="0";}
      if (lastrun<10){ outfile2+="0"; }
      out2 << outfile2 << lastrun ;
      
      string out1filename=out1.str();
      string out2filename=out2.str();
  
  TFile *histfile=new TFile(Form("outtree_%s_%s.root",out1filename.c_str(),out2filename.c_str()), "RECREATE");
  tree = new TTree("dctpc_eventinfo", "Event info");
  tree->Branch("RunNum", &runnum, "runnum/I");
  tree->Branch("SetNum", &setnum, "setnum/I");
  tree->Branch("SequenceNum", &seqnum, "seqnum/I");
  tree->Branch("ExposureInRun_sec", &expose, "expose/I");
  tree->Branch("EventNum", &evnum, "evnum/I");
  tree->Branch("Image_mean_ccdadu", &image_mean, "image_mean/D");
  tree->Branch("Image_rms_ccdadu", &image_rms, "image_rms/D");
  tree->Branch("Edge", &edge, "edge/I");
  tree->Branch("BurnIn", &burnin, "burnin/I");
  tree->Branch("Pixels_killed", &pixels_killed, "pixels_killed/I");
  tree->Branch("LastSpark", &last_spark, "last_spark/I");
  tree->Branch("NextSpark", &next_spark, "next_spark/I");
  tree->Branch("Ntrack", &ntrack, "ntrack/I");
  tree->Branch("Ntrig", &ntrig, "ntrig/I");
  tree->Branch("Etrack_kev", &Etrack, "Etrack/D");
  tree->Branch("Etrig_kev", &Etrig, "Etrig/D");
  tree->Branch("Emesh_kev", &Emesh, "Emesh/D"); 
  tree->Branch("Track_mean_ccdadu", &cluster_mean, "cluster_mean/D");
  tree->Branch("Track_rms_ccdadu", &cluster_rms, "cluster_rms/D"); 
  tree->Branch("Track_x_pix", &TrackX, "TrackX/D"); 
  tree->Branch("Track_y_pix", &TrackY, "TrackY/D");
  tree->Branch("Track_x_start_pix", &TrackXStart, "TrackXStart/D"); 
  tree->Branch("Track_y_start_pix", &TrackYStart, "TrackYStart/D");
  tree->Branch("Track_x_end_pix", &TrackXEnd, "TrackXEnd/D"); 
  tree->Branch("Track_y_end_pix", &TrackYEnd, "TrackYEnd/D");
  tree->Branch("Track_range_pix", &range_ccd, "range_ccd/D");
  tree->Branch("Track_fitwidth_pix", &Trackrms, "Trackrms/D");
  tree->Branch("Track_width_pix", &Trackwidth, "Trackwidth/D");
  tree->Branch("Track_maxpixel_ccdadu", &Trackmaxpixel, "Trackmaxpixel/D");
  tree->Branch("Track_neighbors", &neighbors, "neighbors/I");
  tree->Branch("Track_pixels", &Tracknpixel, "Tracknpixel/I");
  tree->Branch("Track_phi_deg", &phi, "phi/D");
  tree->Branch("Track_skewness", &Trackskewness, "Trackmskewness/D");
  tree->Branch("Anode_rms_v", &anodeRMS, "anodeRMS/D");
  tree->Branch("Mesh_rms_v", &meshRMS, "meshRMS/D");
  tree->Branch("Mesh_peak_v", &mesh_peak, "mesh_peak/D");
  tree->Branch("Mesh_base_v", &mesh_base, "mesh_base/D");
  tree->Branch("Anode_base_v", &anode_base, "anode_base/D");
  tree->Branch("Mesh_max_v", &mesh_max, "mesh_max/D");
  tree->Branch("Anode_max_v", &anode_max, "anode_max/D");
  tree->Branch("Veto_peak_v", &veto_peak, "veto_peak/D");
  tree->Branch("Mesh_starttime_samp", &mesh_start, "mesh_start/D");
  tree->Branch("Anode_starttime_samp", &anode_start, "anode_start/D");
  tree->Branch("Mesh_peaktime_samp", &mesh_peaktime, "mesh_peaktime/D");
  tree->Branch("Mesh_totaltime_samp", &mesh_width, "mesh_width/D");
  tree->Branch("Anode_R0time_samp", &anode_R0, "anode_R0/D");
  tree->Branch("Mesh_R0time_samp", &mesh_R0, "mesh_R0/D");
  tree->Branch("Mesh_R10time_samp", &mesh_R10, "mesh_R10/D");
  tree->Branch("Mesh_R25time_samp", &mesh_R25, "mesh_R25/D");
  tree->Branch("Mesh_R50time_samp", &mesh_R50, "mesh_R50/D");
  tree->Branch("Mesh_R75time_samp", &mesh_R75, "mesh_R75/D");
  tree->Branch("Mesh_R90time_samp", &mesh_R90, "mesh_R90/D");
  tree->Branch("Mesh_F0time_samp", &mesh_F0, "mesh_F0/D");
  tree->Branch("Mesh_F10time_samp", &mesh_F10, "mesh_F10/D");
  tree->Branch("Mesh_F25time_samp", &mesh_F25, "mesh_F25/D");
  tree->Branch("Mesh_F50time_samp", &mesh_F50, "mesh_F50/D");
  tree->Branch("Mesh_F75time_samp", &mesh_F75, "mesh_F75/D");
  tree->Branch("Mesh_F90time_samp", &mesh_F90, "mesh_F90/D");  
  tree->Branch("Timenow_sec",&timenow, "timenow/I");
  tree->Branch("Triggerindex",&triggerindex, "triggerindex/I");
  tree->Branch("Triggertimestamp_samp",&triggertimestamp, "triggertimestamp[ntrig]/D");
  tree->Branch("Mesh_total_allwf_kev", &totalmesh_allwf, "totalmesh_allwf/D");
  tree->Branch("Anode_total_allwf_kev", &totalanode_allwf, "totalanode_allwf/D");
    
    
    
  tree2 = new TTree("dctpc_runinfo", "Run info");
  tree2->Branch("RunNum", &runnum, "runnum2/I");
  tree2->Branch("SetNum", &setnum, "setnum/I");
  tree2->Branch("SequenceNum", &seqnum, "seqnum/I");
  tree2->Branch("Exposure_sec", &exposure, "exposure/I");
  tree2->Branch("Totaltracks", &totaltrack, "totaltrack/I");
  tree2->Branch("Totaltrigs", &totaltrig, "totaltrig/I");
  tree2->Branch("Time_startofrun_sec", &time_start, "time_start/I");
  tree2->Branch("Pressure_torr", &pressure, "pressure/D");
  tree2->Branch("Voltage_amp_volts", &voltage_amp, "voltage_amp/D");
  tree2->Branch("Voltage_drift_volts", &voltage_drift, "voltage_drift/D");
//start run loop

for (int x = firstrun; x <= lastrun; x++)
{


badrunflag=0;
totaltrack=0;
totaltrig=0;
next_spark=10000;
string line;
ifstream infile; 
infile.open("./BadRuns_list_BigDCTPC.dat");
  if (infile.is_open())
  {   
    while ( getline (infile,line) )
    {    	
    
    	istringstream ss(line);
        ss >> badrun;
         
		if(x==badrun)
		{
		badrunflag=1;
		}
	}
  }
infile.close();

if(badrunflag==1)
continue;


      runnum=x;
      
      voltage_drift=6000.;
           
      if(runnum>=889&&runnum<=891)
      {
      seqnum=0;
      lightCalib = .48 * 0.7663;
      anodeCalib = 23.2 * 0.771;
      meshCalib = 37.8 * 0.769;
      }
      
      if(runnum>=893&&runnum<=896)
      {
      seqnum=1;
      lightCalib = .48 * 0.89;
      anodeCalib = 23.2 * 0.8861;
      meshCalib = 37.8 * 0.8853;
      }
      
      if(runnum>=897&&runnum<=904)
      {
      seqnum=2;
      lightCalib = .48 * 0.8887;
      anodeCalib = 23.2 * 0.8849;
      meshCalib = 37.8 * 0.881;
      }
      
      if(runnum>=910&&runnum<=917)
      {
      seqnum=3;
      lightCalib = .48 * 0.8644;
      anodeCalib = 23.2 * 0.857;
      meshCalib = 37.8 * 0.856;
      }
      
      if(runnum>=919&&runnum<=925)
      {
      seqnum=4;
      lightCalib = .48 * 0.7078;
      anodeCalib = 23.2 * 0.715;
      meshCalib = 37.8 * 0.7134;
      }
      
      if(runnum>=926&&runnum<=938)
      {
      seqnum=5;
      lightCalib = .48 * 0.968;
      anodeCalib = 23.2 * 0.8401;
      meshCalib = 37.8 * 0.8392;
      }
      
      if(runnum>=939&&runnum<=954)
      {
      seqnum=6;   
      lightCalib = .48 * 0.9719;
      anodeCalib = 23.2 * 0.8291;
      meshCalib = 37.8 * 0.8284;
      }
      
      if(runnum>=955&&runnum<=1645)
      seqnum=7;
      
      if(runnum>=1646&&runnum<=2205)
      seqnum=8;
      
      if(runnum>=2206&&runnum<=2648)
      seqnum=9; 
      
      if(runnum>=2650&&runnum<=3868)
      seqnum=10; 
        
      if(runnum>=3876&&runnum<=5789)
      seqnum=11; 
      
      if(runnum>=5792&&runnum<=7520)
      seqnum=12;
      
      if(runnum>=7530&&runnum<=9143)
      seqnum=13;

      if(runnum>=9147&&runnum<=12833)
      seqnum=14;
      
      if(runnum>=12837&&runnum<=13311)
      seqnum=15;
      
      if(runnum>=13315&&runnum<=15161)
      seqnum=16;
      
      if(runnum>=15162&&runnum<=16073)
      seqnum=17;

      if(runnum>=16074&&runnum<=16175)
      seqnum=18;

      if(runnum>=16150&&runnum<=17275)
      seqnum=19;	
      
      if(runnum>=17276&&runnum<=17334)
      seqnum=20;
      
      if(runnum>=17336&&runnum<=17349)
      seqnum=21;
      
      if(runnum>=17351&&runnum<=17366)
      seqnum=22;
      
      if(runnum>=17369&&runnum<=17370)
      seqnum=23;

      if(runnum>=17528 && runnum<=17628)
      seqnum=24;

      if(runnum>=17630 && runnum<=17723)
      seqnum=25;

      if(runnum>=17725 && runnum<=18829)
      seqnum=26;
      
      if(runnum>=18836 && runnum<=19675)
      seqnum=27;

      if(runnum>=19678 && runnum<=20232)
      seqnum=28;

      if(runnum>=20234 && runnum<=21625)
      seqnum=29;	

      if(runnum>=21733 && runnum<=25399)
      seqnum=30;
      
      if(runnum>=25432 && runnum<=26854)
      seqnum=31;
      
      if(runnum>=26855)
      seqnum=32;

      if(runnum>=889&&runnum<=954)
      setnum=0;
      
      if(runnum>=955&&runnum<=3868)
      setnum=1;
      
      if(runnum>=3876&&runnum<=5789)
      setnum=2;

      if(runnum>=5792&&runnum<=17334)
      {
      setnum=3;
      voltage_drift=7500.;
      }
      
      if(runnum>=5792&&runnum<=17334)
      {
      setnum=3;
      voltage_drift=7500.;
      }
      
      if(runnum>=17336&&runnum<=17366)
      {
      setnum=4;
      voltage_drift=6200.;
      }
      
      if(runnum>=17369&&runnum<=17370)
      {
      setnum=5;
      voltage_drift=4480.;
      }

      if(runnum>=17528&&runnum<=26854)
      {
      setnum=6;
      voltage_drift=5500.;
      }
      
      if(runnum>=26855)
      {
      setnum=7;
      voltage_drift=5500.;
      }
      
      if(setnum==-1||seqnum==-1)
      continue;

      
      if(seqnum!=seqnum2)
      seqcounter=0;
      
      seqnum2=seqnum;
    
      if(BATCH==0)
      {
      origfile = "/net/nudsk0001/d00/scratch/dctpc_tmp/bigdctpc_data/BigDCTPC_run_";
      skimfile = "/net/nudsk0001/d00/scratch/dctpc_tmp/bigdctpc_skim/BigDCTPC_run_";
      }
      else
      {
      origfile = "./BigDCTPC_run_";
      skimfile = "./BigDCTPC_run_";   
      }
//       origfile = "/mnt/hadoop/mitgroups/dchooz/dctpc_tmp/bigdctpc_data/BigDCTPC_run_";
//       skimfile = "/mnt/hadoop/mitgroups/dchooz/dctpc_tmp/bigdctpc_skim/BigDCTPC_run_";

//uncomment all this for batch running:

      
      string skimend = "skim.root";
      string origend = ".root";
      string origfilename;
      string skimfilename;
      sstm.str("");
      sstmxxx.str("");
      string xxx;
      if (x<10000){ origfile+="0"; skimfile+="0"; xxx+="0";}
      if (x<1000){ origfile+="0"; skimfile+="0"; xxx+="0";}
      if (x<100){ origfile+="0"; skimfile+="0"; xxx+="0";}
      if (x<10){ origfile+="0"; skimfile+="0"; xxx+="0";}
      sstm << origfile << x << origend;
      sstmxxx << xxx << x;
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

if(BATCH==1)
{
   stringstream t;
   t<<gSystem->GetFromPipe(Form("du -sb /mnt/hadoop/mitgroups/dchooz/dctpc_tmp/bigdctpc_skim/%s |awk '{print $1}'",skimfilename.c_str()));
   t >> orig_filesize;
   
   stringstream t;
   t<<gSystem->GetFromPipe(Form("du -sb %s |awk '{print $1}'",skimfilename.c_str()));
   t >> new_filesize;

   stringstream t;
   t<<gSystem->GetFromPipe(Form("du -sb /mnt/hadoop/mitgroups/dchooz/dctpc_tmp/bigdctpc_data/%s |awk '{print $1}'",origfilename.c_str()));
   t >> orig_filesize2;
   
   stringstream t;
   t<<gSystem->GetFromPipe(Form("du -sb %s |awk '{print $1}'",origfilename.c_str()));
   t >> new_filesize2;

if(orig_filesize!=new_filesize)
continue;

if(orig_filesize2!=new_filesize2)
continue;

}


      DmtpcSkimDataset d;
      try
      {
      d.openRootFile(skimfilename.c_str());     
      d.loadDmtpcEvent(true,origfilename.c_str());
      }
      catch(const std::exception& e)
      {
	  continue;	
      }
      DmtpcSkimDataset d2;
      try
      {
      d2.openRootFile(skimfilename.c_str());     
      d2.loadDmtpcEvent(true,origfilename.c_str());
      }
      catch(const std::exception& e)
      {
	  continue;	
      }
      
      
      Long_t nEvents=d.tree()->GetEntriesFast();
	  runnum2=x; 
	  exposure=0;	
	  expose=0;
	  lastSpark=-1;
	  
std::string s = sstmxxx.str();

stringstream ff;
ff<<gSystem->GetFromPipe(Form("grep run_%s /net/hisrv0001/home/spitzj/runlog|awk '{print $46}'",s.c_str()));
ff >> pressure;

stringstream fff;
fff<<gSystem->GetFromPipe(Form("grep run_%s /net/hisrv0001/home/spitzj/runlog|awk '{print $30}'",s.c_str()));
fff >> voltage_amp;

voltage_amp*=1000.;
// cout<<pressure<<" "<<voltage_drift<<" "<<voltage_amp<<endl;

seqcounter++;

std::vector<int> sparkvector;
for(int i=0; i<nEvents; i++){ 

DmtpcSkimEvent* ev2=d2.event(); //creating an event object		
d2.getEvent(i);

if(ev2->spark(0)){sparkvector.push_back(ev2->eventNumber());lastSpark = ev2->eventNumber();continue;} 
	
ntrack=ev2->ntracks(0);
for(int k=0;k<ev2->ntracks(0);k++)
{	
	  if (ev2->maxpixel(0,k)>1000||ev2->npixel(0,k)>1500.){
	  sparkvector.push_back(ev2->eventNumber());
	  lastSpark = ev2->eventNumber();
	  continue;}
}

if(i-lastSpark>=7)
expose++;

}

lastSpark = -1;

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
	
	ntrack=ev->ntracks(0);
	for(int k=0;k<ev->ntracks(0);k++)
	  {	
	  if (ev->maxpixel(0,k)>1000||ev->npixel(0,k)>1500.){lastSpark = ev->eventNumber();continue;}//maxpixel has to be less than 1000ADU. also, the event might have been a spark! also, sparks tend to have a lot of pixels.
	  }
	if(evnum-lastSpark<7){continue;} // make sure that at least 6 events have passed after a spark

	for(int a=0;a<sparkvector.size();a++)
	{
	if(i<sparkvector[a])
	{
	next_spark=sparkvector[a]-evnum;
	break;
	}
	next_spark=10000;
	}
	

    exposure++;//increase the integrated exposure as long as the event is not a spark or close to a spark in time

	//if(ev->ntracks(0)==0){continue;}		

	Etrack=-1.;
	Etrig=-1.;
	Emesh=-1.;
	passTrack=0;
	passTrig=0;
	passAll=0;
	timenow = ev->timenow(0);
	bestDiff=1e9;
    Etrack = -1.;
    phi = -1.;   
    range_ccd=-1.;	         	      	    
    TrackX = -1.;
    TrackY = -1.; 
    Trackskewness=-1.;
    edge=-1;
    burnin=-1;
    cluster_mean=-1.;
    cluster_rms=-1.;
    neighbors=-1;
    Trackmaxpixel=-1.;
    Tracknpixel=-1;
    Trackrms=-1.;  
    Trackwidth=-1.;
    TrackXStart = -1.;
    TrackYStart = -1.;
    TrackXEnd = -1.;
    TrackYEnd = -1.;

image_mean=ev->image_mean(0);
image_rms=ev->image_rms(0);
last_spark=evnum-lastSpark;
//cout<<image_rms<<endl;
pixels_killed=ev->pixels_killed(0);

	diff=1e9;	
	
	if(d.orig_event()->scopeDataInfo(0)==NULL){continue;}//no waveforms!!

    TObjArray* arr = ev->waveform_vectors();
    CspWfVector* anode = (CspWfVector*) arr->At(0);
	FastWfVector* mesh = (FastWfVector*) arr->At(1);
	CspWfVector* veto = (CspWfVector*) arr->At(2); 

totalmesh_allwf=0.;
totalanode_allwf=0.;

totaltrig+=mesh->size();
totaltrack+=ev->ntracks(0);

	    //loop over triggers
	    ntrig=mesh->size();
	for (int nt = 0; nt < mesh->size(); nt++)
	{	
	ScopeWaveformData *hwf=d.orig_event()->rawScopeData(nt);	
	triggertimestamp[nt]=hwf->getTimeStamp();

	totalmesh_allwf+=mesh->at(nt,0).getIntegral()*meshCalib;
	totalanode_allwf+=(anode->at(nt,0).getPeak()*1000.)*anodeCalib;
	
	//loop over tracks in each event
	  for(int k=0;k<ev->ntracks(0);k++)
	  {	 
 
 
	    if (ev->E(0,k)<1){continue;} //Cut if E is below 1 (measured in ADU)
	    if (ev->range(0,k)<=0){continue;} //Cut if it has a zero range
	    if (ev->maxpixel(0,k)/ev->E(0,k)>0.25){continue;} //To eliminate noise artifacts the max pixel has to be less than 25% of the total signal      
           
	    passTrig=0;


	    passTrig++;//increment number of tracks that pass the trigger cuts
	    E1 = ev->E(0,k)*lightCalib;
	    E2 = (anode->at(nt,0).getPeak()*1000.)*anodeCalib;
	    diff=TMath::Abs( (E1 - E2)/(E1 + E2) );
	    
	    
	    //cout<<ev->minorAxis(0,k)<<endl;
	    	if (fabs(diff)<fabs(bestDiff)) 
	    	{	
			Etrack = E1;
			phi = ev->phi(0,k);
			phi = atan2(sin(phi),cos(phi)) * 180./ 3.1416;   
			range_ccd=ev->range(0,k);	         	      	    
	    	TrackX = ev->x(0,k)-512;
	    	TrackY = ev->y(0,k)-512; 
        	Trackskewness=ev->skewness(0,k);
        	edge=ev->edge(0,k);
        	burnin=ev->nburnin(0,k);
        	cluster_mean=ev->cluster_mean(0,k);
        	cluster_rms=ev->cluster_rms(0,k);
        	neighbors=ev->neighbors(0,k);
        	Trackmaxpixel=ev->maxpixel(0,k);
        	Tracknpixel=ev->npixel(0,k);
        	Trackrms=sqrt(ev->transverse_moment(0,2,k)); 
        	Trackwidth=4*ev->npixel(0,k)/ev->range(0,k);
        	//cout<<ev->transverse_moment(0,1,k)<<" "<<ev->transverse_moment(0,2,k)<<endl; 
        	 
        	 
        	TrackXStart = (double)ev->xbegin(0,k)-512;
        	TrackYStart = (double)ev->ybegin(0,k)-512;
        	TrackXEnd = (double)ev->xend(0,k)-512;
        	TrackYEnd = (double)ev->yend(0,k)-512;
	    	
			Etrig = E2;			
			mesh_peak=mesh->at(nt,0).getPeak();
			veto_peak=veto->at(nt,0).getPeak();
			anode_R0=anode->at(nt,0).getRise0();
			mesh_R0=mesh->at(nt,0).getRise0();
			mesh_R10=mesh->at(nt,0).getRise10();
			mesh_R25=mesh->at(nt,0).getRise25();
			mesh_R50=mesh->at(nt,0).getRise50();
			mesh_R75=mesh->at(nt,0).getRise75();
			mesh_R90=mesh->at(nt,0).getRise90();
			mesh_F0=mesh->at(nt,0).getFall0();
			mesh_F10=mesh->at(nt,0).getFall10();
			mesh_F25=mesh->at(nt,0).getFall25();
			mesh_F50=mesh->at(nt,0).getFall50();
			mesh_F75=mesh->at(nt,0).getFall75();
			mesh_F90=mesh->at(nt,0).getFall90();
			mesh_width=mesh->at(nt,0).getRise0()+mesh->at(nt,0).getFall0();
			mesh_peaktime=mesh->at(nt,0).getSlowPeakTime();
			Emesh=mesh->at(nt,0).getIntegral()*meshCalib;
			meshRMS=mesh->at(nt).getRMS();
	        anodeRMS=anode->at(nt).getRMS();
	        mesh_start=mesh->at(nt,0).getStartTime();	
	        anode_start=anode->at(nt,0).getStartTime();	
	        mesh_base=mesh->at(nt).getBase();
	        anode_base=anode->at(nt).getBase();
	        mesh_max=mesh->at(nt).getWfMax();
	        anode_max=anode->at(nt).getWfMax();
	        triggerindex=nt;
	        bestDiff = diff;	        
	    	}
	    		    	
	    } // end track loop 
	   
	 //    cout<<evnum<<" "<<diff<<" "<<bestDiff<<" "<<ev->ntracks(0)<<endl;
	    
	   if(ev->ntracks(0)==0) 
	    {
	    E1 = mesh->at(nt,0).getIntegral()*meshCalib;
	    E2 = (anode->at(nt,0).getPeak()*1000.)*anodeCalib;
	    diff=TMath::Abs( (E1 - E2)/(E1 + E2) );
	    
	   
	    
	    	if (fabs(diff)<fabs(bestDiff)) 
	    	{	
	    	// cout<<evnum<<endl;
	    	
	    	Etrack = -1.;
			phi = -1.;   
			range_ccd=-1.;	         	      	    
	    	TrackX = -1.;
	    	TrackY = -1.; 
        	Trackskewness=-1.;
        	edge=-1;
        	burnin=-1;
        	cluster_mean=-1.;
        	cluster_rms=-1.;
        	neighbors=-1;
        	Trackmaxpixel=-1.;
        	Tracknpixel=-1;
        	Trackrms=-1.;  
        	Trackwidth=-1.;
        	TrackXStart = -1.;
        	TrackYStart = -1.;
        	TrackXEnd = -1.;
        	TrackYEnd = -1.;

			Etrig = E2;			
			mesh_peak=mesh->at(nt,0).getPeak();
			veto_peak=veto->at(nt,0).getPeak();
			anode_R0=anode->at(nt,0).getRise0();
			mesh_R0=mesh->at(nt,0).getRise0();
			mesh_R10=mesh->at(nt,0).getRise10();
			mesh_R25=mesh->at(nt,0).getRise25();
			mesh_R50=mesh->at(nt,0).getRise50();
			mesh_R75=mesh->at(nt,0).getRise75();
			mesh_R90=mesh->at(nt,0).getRise90();
			mesh_F0=mesh->at(nt,0).getFall0();
			mesh_F10=mesh->at(nt,0).getFall10();
			mesh_F25=mesh->at(nt,0).getFall25();
			mesh_F50=mesh->at(nt,0).getFall50();
			mesh_F75=mesh->at(nt,0).getFall75();
			mesh_F90=mesh->at(nt,0).getFall90();
			mesh_width=mesh->at(nt,0).getRise0()+mesh->at(nt,0).getFall0();
			mesh_peaktime=mesh->at(nt,0).getSlowPeakTime();
			Emesh=mesh->at(nt,0).getIntegral()*meshCalib;
			meshRMS=mesh->at(nt).getRMS();
	        anodeRMS=anode->at(nt).getRMS();
	        mesh_start=mesh->at(nt,0).getStartTime();	
	        anode_start=anode->at(nt,0).getStartTime();
	        mesh_base=mesh->at(nt).getBase();
	        anode_base=anode->at(nt).getBase();
	        mesh_max=mesh->at(nt).getWfMax();
	        anode_max=anode->at(nt).getWfMax();
	        triggerindex=nt;	
	        bestDiff = diff;		
	    	} 
	    
	    }
    
	    ntotalTrig+=passTrig;
	  }//end trigger loop

	ntotalTrack+=passTrack;
	
	//if (passTrack!=1){continue;} // events are unlikely to have more than 1 track
// 	if (passTrig==0){continue;} 

	passAll=1;
		//cout<<Etrack<<" "<<Etrig<<" "<<Emesh<<endl;
	tree->Fill();	
	
}//end event loop

tree2->Fill();
cout <<"Total seconds of exposure: "<<exposure<<" "<<expose<<endl;
}//end file loop


histfile->Write();
histfile->Close();
if(BATCH==1)
gApplication->Terminate(); 
}


