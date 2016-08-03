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

void analysis_125CF4_dctpc_LI8_sim(int firstrun, int lastrun){//inclusive 
 
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

  double leftint=0;
  double leftint2=0;
  double rightint=0;
  double rightint2=0;
  int jmin=0;
  int jbragg=0;
  int jterm=0;
  int jterm2=0;
  int jterm3=0;
  int origin=0;
  double rec_SD=0;
  int wfd_delta=0;
  int peak1=0;
  int peak2=0;
  double peak1val=0;
  double peak2val=0;
  int half1=0;
  int half2=0;
  double termdist=0;
  int maxdevloc=0;
  double rms_left=0;
  double rms_right=0;
  double rms_outer=0;
  double rms_full=0;
  
  int estimate=0;
  int estimate_error=0;
  double delta_E=0;
  double delta_E_true=0;
  double delta_E_err=0;
  double delta_E_cheat=0;

  double truth_phi_deg=-10.;
  double truth_x_start_mm=0.;
  double truth_y_start_mm=0.;
  double truth_x_start_pix=0.;
  double truth_y_start_pix=0.;
  double truth_z_start_mm=0.;
  int truth_sequence=-1;
  double truth_x_end_mm=0.;
  double truth_y_end_mm=0.;
  double truth_x_end_pix=0.;
  double truth_y_end_pix=0.;
  double truth_z_end_mm=0.;
  double truth_theta_deg=-10.;
  double truth_particlee_kev=0.;
  double truth_particlee_kev2=0.;
  double truth_depositede_ccdadu=-100.;
  double truth_range_pix=-1.;
  double truth_range_mm=-1.;
  double truth_zrange_mm=-1.;
  double truth_lengthcalibration=0.;
  double truth_pressure_torr=0.;
  double truth_gain=0.;
  double truth_noise=0.;
  double truth_bias=0.;

  double truth_phi_deg2=-10.;
  double truth_x_start_mm2=0.;
  double truth_y_start_mm2=0.;
  double truth_x_start_pix2=0.;
  double truth_y_start_pix2=0.;
  double truth_z_start_mm2=0.;
  int truth_sequence2=-1;
  double truth_x_end_mm2=0.;
  double truth_y_end_mm2=0.;
  double truth_x_end_pix2=0.;
  double truth_y_end_pix2=0.;
  double truth_z_end_mm2=0.;
  double truth_theta_deg2=-10.;
  double truth_range_pix2=-1.;
  double truth_range_mm2=-1.;
  double truth_zrange_mm2=-1.;


long int orig_filesize=-1;
long int new_filesize=0;
long int orig_filesize2=-1;
long int new_filesize2=0;

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

  tree->Branch("wfr_LeftInt", &leftint, "leftint/D");
  tree->Branch("wfr_LeftInt2", &leftint2, "leftint2/D");
  tree->Branch("wfr_RightInt", &rightint, "rightint/D");
  tree->Branch("wfr_RightInt2", &rightint2, "rightint2/D");
  tree->Branch("wfr_jMin", &jmin, "jmin/I");
  tree->Branch("wfr_jBragg", &jbragg, "jbragg/I");
  tree->Branch("wfr_jTerm", &jterm, "jterm/I");
  tree->Branch("wfr_jTerm2", &jterm2, "jterm2/I");
  tree->Branch("wfr_jTerm3", &jterm3, "jterm3/I");
  tree->Branch("wfr_Origin", &origin, "origin/I");
  tree->Branch("wfr_RecSD", &rec_SD, "rec_SD/D");
  tree->Branch("wfr_WfdDelta", &wfd_delta, "wfd_delta/I");
  tree->Branch("wfr_Peak1", &peak1, "peak1/I");
  tree->Branch("wfr_Peak2", &peak2, "peak2/I");
  tree->Branch("wfr_Peak1Val", &peak1val, "peak1val/D");
  tree->Branch("wfr_Peak2Val", &peak2val, "peak2val/D");
  tree->Branch("wfr_Half1", &half1, "half1/I");
  tree->Branch("wfr_Half2", &half2, "half2/I");
  tree->Branch("wfr_TermDist", &termdist, "termdist/D");
  tree->Branch("wfr_MaxDevLoc", &maxdevloc, "maxdevloc/I");
  tree->Branch("wfr_RMSLeft", &rms_left, "rms_left/D");
  tree->Branch("wfr_RMSRight", &rms_right, "rms_right/D");
  tree->Branch("wfr_RMSFull", &rms_full, "rms_full/D");
  tree->Branch("wfr_RMSOuter", &rms_outer, "rms_outer/D");
  tree->Branch("wfr_Origin_rec", &estimate, "estimate/I");
  tree->Branch("wfr_Origin_err", &estimate_error, "estimate_error/I");
  tree->Branch("wfr_DE",&delta_E,"delta_E/D");
  tree->Branch("TRUTH_DE",&delta_E_true,"delta_E_true/D");
  tree->Branch("wfr_DE_err",&delta_E_err,"delta_E_err/D");
  tree->Branch("wfr_DE_cheat",&delta_E_cheat,"delta_E_cheat/D");

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



  
  tree->Branch("TRUTH_sequence", &truth_sequence, "truth_sequence/I");
  tree->Branch("TRUTH_phi_deg", &truth_phi_deg, "truth_phi_deg/D");
  tree->Branch("TRUTH_x_start_mm", &truth_x_start_mm, "truth_x_start_mm/D");
  tree->Branch("TRUTH_y_start_mm", &truth_y_start_mm, "truth_y_start_mm/D");
  tree->Branch("TRUTH_x_start_pix", &truth_x_start_pix, "truth_x_start_pix/D");
  tree->Branch("TRUTH_y_start_pix", &truth_y_start_pix, "truth_y_start_pix/D");
  tree->Branch("TRUTH_z_start_mm", &truth_z_start_mm, "truth_z_start_mm/D");
  tree->Branch("TRUTH_x_end_mm", &truth_x_end_mm, "truth_x_end_mm/D");
  tree->Branch("TRUTH_y_end_mm", &truth_y_end_mm, "truth_y_end_mm/D");
  tree->Branch("TRUTH_x_end_pix", &truth_x_end_pix, "truth_x_end_pix/D");
  tree->Branch("TRUTH_y_end_pix", &truth_y_end_pix, "truth_y_end_pix/D");
  tree->Branch("TRUTH_z_end_mm", &truth_z_end_mm, "truth_z_end_mm/D");
  tree->Branch("TRUTH_theta_deg", &truth_theta_deg, "truth_theta_deg/D");
  tree->Branch("TRUTH_particleE_kev", &truth_particlee_kev, "truth_particlee_kev/D");
  tree->Branch("TRUTH_range_pix", &truth_range_pix, "truth_range_pix/D");
  tree->Branch("TRUTH_range_mm", &truth_range_mm, "truth_range_mm/D");
  tree->Branch("TRUTH_zrange_mm", &truth_zrange_mm, "truth_zrange_mm/D");


  //tree->Branch("TRUTH_sequence2", &truth_sequence2, "truth_sequence2/I");
  tree->Branch("TRUTH_phi_deg2", &truth_phi_deg2, "truth_phi_deg2/D");
  //tree->Branch("TRUTH_x_start_mm2", &truth_x_start_mm2, "truth_x_start_mm2/D");
  //tree->Branch("TRUTH_y_start_mm2", &truth_y_start_mm2, "truth_y_start_mm2/D");
  //tree->Branch("TRUTH_x_start_pix2", &truth_x_start_pix2, "truth_x_start_pix2/D");
  //tree->Branch("TRUTH_y_start_pix2", &truth_y_start_pix2, "truth_y_start_pix2/D");
  //tree->Branch("TRUTH_z_start_mm2", &truth_z_start_mm2, "truth_z_start_mm2/D");
  tree->Branch("TRUTH_x_end_mm2", &truth_x_end_mm2, "truth_x_end_mm2/D");
  tree->Branch("TRUTH_y_end_mm2", &truth_y_end_mm2, "truth_y_end_mm2/D");
  tree->Branch("TRUTH_x_end_pix2", &truth_x_end_pix2, "truth_x_end_pix2/D");
  tree->Branch("TRUTH_y_end_pix2", &truth_y_end_pix2, "truth_y_end_pix2/D");
  tree->Branch("TRUTH_z_end_mm2", &truth_z_end_mm2, "truth_z_end_mm2/D");
  tree->Branch("TRUTH_theta_deg2", &truth_theta_deg2, "truth_theta_deg2/D");
  tree->Branch("TRUTH_particleE_kev2", &truth_particlee_kev2, "truth_particlee_kev2/D");
  tree->Branch("TRUTH_range_pix2", &truth_range_pix2, "truth_range_pix2/D");
  tree->Branch("TRUTH_range_mm2", &truth_range_mm2, "truth_range_mm2/D");
  tree->Branch("TRUTH_zrange_mm2", &truth_zrange_mm2, "truth_zrange_mm2/D");




  tree->Branch("TRUTH_depositedE_ccdadu", &truth_depositede_ccdadu, "truth_depositede_ccdadu/D");

  tree->Branch("TRUTH_lengthcalibration", &truth_lengthcalibration, "truth_lengthcalibration/D");
  tree->Branch("TRUTH_pressure_torr", &truth_pressure_torr, "truth_pressure_torr/D");
  tree->Branch("TRUTH_gain", &truth_gain, "truth_gain/D");
  tree->Branch("TRUTH_noise_ccdadu", &truth_noise, "truth_noise/D");
  tree->Branch("TRUTH_bias", &truth_bias, "truth_bias/D");
    
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

totaltrack=0;
totaltrig=0;


      if(BATCH==0)
      {
      string origfile = "/net/hisrv0001/home/spitzj/tj/DCTPC_soft/MaxCam/Simulations/v1/dmtpc_mc_";
      string skimfile = "/net/hisrv0001/home/spitzj/tj/DCTPC_soft/MaxCam/Simulations/v1/skim/dmtpc_mc_";    
	  }
// 	  if(BATCH==1)
// 	  {
// 	  string origfile = "./myOutput_";
//       string skimfile = "./myOutput_"; 
// 	  }

 
 //cout << origfilename << endl;
 
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
      //cout << origfilename << endl;
      ifstream ifile(origfilename.c_str());
      ifstream ifile2(skimfilename.c_str());
      if(!ifile)
	continue;
      if(!ifile2)
	continue;
	
      DmtpcSkimDataset d;
      d.openRootFile(skimfilename.c_str());     
      d.loadDmtpcEvent(true,origfilename.c_str());
      
      DmtpcMCDataset mc;
      mc.loadFile(origfilename.c_str());

      Long_t nEvents=d.tree()->GetEntriesFast();
	  runnum2=x;  
	  exposure=0;	
	  expose=0;
	  lastSpark=-1;

//start event loop 
for(int i=0; i<nEvents; i++){    

// if(i%100==0)
cout<<"Event: "<<i<<" / "<<nEvents<<endl;
    
evnum=i;	
DmtpcSkimEvent* ev=d.event(); //creating an event object		
d.getEvent(i);

 mc.getEvent(i);
 truth_phi_deg=mc.getPhi()*180./3.14159;
 truth_x_start_mm=mc.getX(true);
 truth_y_start_mm=mc.getY(true);
 truth_x_start_pix=mc.getX()-512.;
 truth_y_start_pix=mc.getY()-512.;
 truth_z_start_mm=mc.getZ();
 truth_sequence=mc.getSequence();
 truth_theta_deg=mc.getTheta()*180./3.14159;
 truth_particlee_kev=mc.getParticleE();
 truth_range_pix=mc.getRange();
 truth_range_mm=mc.getRange()*mc.getLengthcal();
 truth_zrange_mm=mc.getRangeZ();
 truth_x_end_mm=truth_x_start_mm+(truth_range_mm*TMath::Sin(mc.getTheta())*TMath::Cos(mc.getPhi()));
 truth_y_end_mm=truth_y_start_mm+(truth_range_mm*TMath::Sin(mc.getTheta())*TMath::Sin(mc.getPhi()));
 truth_z_end_mm=truth_z_start_mm+(truth_range_mm*TMath::Cos(mc.getTheta()));
 truth_x_end_pix=truth_x_start_pix+(truth_range_pix*TMath::Sin(mc.getTheta())*TMath::Cos(mc.getPhi()));
 truth_y_end_pix=truth_y_start_pix+(truth_range_pix*TMath::Sin(mc.getTheta())*TMath::Sin(mc.getPhi()));

 truth_phi_deg2=mc.getPhi2()*180./3.14159;
 truth_x_start_mm2=truth_x_start_mm;//mc.getX2(true);
 truth_y_start_mm2=truth_y_start_mm;//mc.getY2(true);
 truth_x_start_pix2=truth_x_start_pix;//mc.getX2()-512.;
 truth_y_start_pix2=truth_y_start_pix;//mc.getY2()-512.;
 truth_z_start_mm2=truth_z_start_mm;//mc.getZ2();
 //truth_sequence2=mc.getSequence2();
 truth_theta_deg2=mc.getTheta2()*180./3.14159;
 truth_particlee_kev2=mc.getParticleE2();
 truth_range_pix2=mc.getRange2();
 truth_range_mm2=mc.getRange2()*mc.getLengthcal();
 truth_zrange_mm2=mc.getRangeZ2();
 truth_x_end_mm2=truth_x_start_mm2+(truth_range_mm2*TMath::Sin(mc.getTheta2())*TMath::Cos(mc.getPhi2()));
 truth_y_end_mm2=truth_y_start_mm2+(truth_range_mm2*TMath::Sin(mc.getTheta2())*TMath::Sin(mc.getPhi2()));
 truth_z_end_mm2=truth_z_start_mm2+(truth_range_mm2*TMath::Cos(mc.getTheta2()));
 truth_x_end_pix2=truth_x_start_pix2+(truth_range_pix2*TMath::Sin(mc.getTheta2())*TMath::Cos(mc.getPhi2()));
 truth_y_end_pix2=truth_y_start_pix2+(truth_range_pix2*TMath::Sin(mc.getTheta2())*TMath::Sin(mc.getPhi2()));


truth_depositede_ccdadu=mc.getE();

truth_lengthcalibration=mc.getLengthcal();
truth_pressure_torr=mc.getPressure();
truth_gain=mc.getGain();
truth_noise=mc.getNoise();
truth_bias=mc.getBias();


// cout<<i<<endl;
// cout<<truth_x_start_pix<<" "<<truth_y_start_pix<<endl;
//cout<<truth_x_end_pix<<" "<<truth_y_end_pix<<endl;


if(i==0)
time_start=ev->timenow(0);
	

    exposure++;//increase the integrated exposure as long as the event is not a spark or close to a spark in time

	//if(ev->ntracks(0)==0){continue;}		
    E1=-1.;
	Etrack=-1.;
	phi=-10.;
    range_ccd=-1.; 
    TrackX=-10000.;
    TrackY=-10000.;
    Trackskewness=-1.;
    edge=-10;
    burnin=-10.;
    cluster_mean=0.;
    cluster_rms=0.;
    neighbors=-10;
    Trackmaxpixel=-1.;
    Tracknpixel=-1;
    Trackrms=0.;
    TrackXStart=-10000.;
    TrackYStart=-10000.;
    TrackXEnd=-10000.;
    TrackYEnd=-10000.;
     
	passTrack=0;
	passTrig=0;
	passAll=0;
	timenow = ev->timenow(0);

image_mean=ev->image_mean(0);
image_rms=ev->image_rms(0);
pixels_killed=ev->pixels_killed(0);

diff=1e9;	

totaltrack+=ev->ntracks(0);
ntrack=ev->ntracks(0);

//waveform filling code:

TObjArray* arr = ev->waveform_vectors();
FastWfVector* mesh = (FastWfVector*) arr->At(0);
ntrig=mesh->size();

	for (int nt = 0; nt < mesh->size(); nt++)
	{
	  mesh_peak=mesh->at(nt,0).getPeak();
	 
	  leftint=mesh->at(nt).getLeftInt();
	  leftint2=mesh->at(nt).getLeftInt2();
	  rightint=mesh->at(nt).getRightInt();
	  rightint2=mesh->at(nt).getRightInt2();
	  jmin=mesh->at(nt).getjMin();
	  jbragg=mesh->at(nt).getjBragg();
	  jterm=mesh->at(nt).getjTerm();
	  jterm2=mesh->at(nt).getjTerm2();
	  jterm3=mesh->at(nt).getjTerm3();
	  origin=mesh->at(nt).getOrigin();
	  rec_SD=mesh->at(nt).getRecSD();
	  wfd_delta=mesh->at(nt).getWfdDelta();
	  peak1=mesh->at(nt).getPeak1();
	  peak2=mesh->at(nt).getPeak2();
	  peak1val=mesh->at(nt).getPeak1Val();
	  peak2val=mesh->at(nt).getPeak2Val();
	  half1=mesh->at(nt).getHalf1();
	  half2=mesh->at(nt).getHalf2();
	  termdist=mesh->at(nt).getTermDist();
	  maxdevloc=mesh->at(nt).getMaxDevLoc();
	  rms_left=mesh->at(nt).getRMSLeft();
	  rms_right=mesh->at(nt).getRMSRight();
	  rms_outer=mesh->at(nt).getRMSOuter();
	  rms_full=mesh->at(nt).getRMSFull();
	  estimate = jterm - 2;
	  estimate_error=estimate-origin;
	  delta_E=fabs(leftint-rightint);
	  delta_E_cheat=fabs(leftint2-rightint2);
          delta_E_true=fabs(truth_particlee_kev-truth_particlee_kev2);
	  delta_E_err= delta_E - delta_E_true;


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
	  Emesh=mesh->at(nt,0).getIntegral()*meshCalib;
	  mesh_width=mesh->at(nt,0).getRise0()+mesh->at(nt,0).getFall0();
	  mesh_peaktime=mesh->at(nt,0).getSlowPeakTime();
	  meshRMS=mesh->at(nt).getRMS();
	  mesh_start=mesh->at(nt,0).getStartTime();
	  mesh_base=mesh->at(nt).getBase();
	  mesh_max=mesh->at(nt).getWfMax();
	}
	
//end of waveform filling code
	
	//loop over tracks in each event
	  for(int k=0;k<ev->ntracks(0);k++)
	  {	 

	        E1 = ev->E(0,k)*lightCalib;
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
        	TrackXStart = (double)ev->xbegin(0,k)-512;
        	TrackYStart = (double)ev->ybegin(0,k)-512;
        	TrackXEnd = (double)ev->xend(0,k)-512;
        	TrackYEnd = (double)ev->yend(0,k)-512;
	    	
	    		    	
	    } // end track loop 


	ntotalTrack+=passTrack;
	
	//if (passTrack!=1){continue;} // events are unlikely to have more than 1 track
// 	if (passTrig==0){continue;} 

	passAll=1;
		//cout<<Etrack<<" "<<Etrig<<" "<<Emesh<<endl;
	tree->Fill();	
	
}//end event loop

tree2->Fill();
}//end file loop


histfile->Write();
histfile->Close();
gApplication->Terminate(); 
}


