void analysis_125CF4_bigdctpc_neutron(int firstrun, int lastrun){//inclusive 

  // CF252 (neutron) runs: 1040-1158
  // Note that the 'skimmed' runs need to exist before this macro will run
  // To run over all CF252 runs (1040-1158, inclusive): 
  // root
  // .L analysis.C
  // analysis(1040, 1158) 

  //For running at dcnode, the input file location needs to be changed

  //The program creates some relevant histograms and a TTree. The TTree is only relevant if this program is run in parallel for processing a large number of runs. 
  // All of the cuts that are utilized are explained in arXiv:1108.4894

  double lightCalib = .449;
  double anodeCalib = 17.54;
  double meshCalib = 6.55;
  double mmperpixel = 166.6/1024;
  bool useRadial=true;
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1,0); 
  gSystem->Load("libMaxCam");
  gSystem->Load("libWaveformTools.so");
  std::stringstream sstm;
  vector<double> noiseVector;
  vector<double> eventnums;
  vector<int> runnums;
  DmtpcSkimDataset d;
  int totalevents = 0;
  int passTrack = 0;
  int passTrig = 0;
  int passAll = 0;
  int ntotalTrack = 0;
  int ntrig=0;
  int ntotalTrig = 0;
  int lastSpark = -100; 
  int cutnum=0;
  int spark, runnum, evnum;
  int timenow;
  double xx, yy, diff, E1, E2, E3, phi;
  double range_ccd=-1.;
  double Etrack=-1.;
  double Etrig=-1.;
  double Emesh=-1.;
  double TrackX=-1.;
  double TrackY=-1.;
  double bestDiff=-1.;
  double phi=-10.;
  double r2=-1.;
  double exposure=0.;
  double anodeRMS=0.;
  double meshRMS=0.;
  TH1D *hist_diff=new TH1D("Charge-light diff","Charge-light diff",50,-300,300);
  
  TH1D *hist_energy_ccd=new TH1D("CCD energy","CCD energy",50,10,1000);
  TH1D *hist_energy_anode=new TH1D("Anode energy","Anode energy",50,10,1000);
  TH1D *hist_energy_mesh=new TH1D("Mesh energy","Mesh energy",50,10,1000);
  TH1D *hist_phi=new TH1D("phi","phi",50,0.,180.);
  TH1D *hist_cosphi=new TH1D("cosphi","cosphi",50,-1.,1.);
  TH2D *hist_pos=new TH2D("pos","pos",128,-512,512,128,-512,512);
  TH2D *hist_energy=new TH2D("Energy","Energy",50,10,1000,50,10,1000);
  TH2D *hist_energy2=new TH2D("Anode v Mesh","Anode v Mesh",50,10,1000,50,10,1000);
  TH2D *hist_energy3=new TH2D("CCD v Mesh","CCD v Mesh",50,10,1000,50,10,1000);
  TH2D *hist_energy_range_ccd=new TH2D("Energy Range CCD","Energy Range CCD",70,10,1000,40,0,15);

  TH2D *  hist_energy_range_anode=new TH2D("Energy Range Anode","Energy Range Anode",70,10,1000,40,0,15);
    TH2D *  hist_energy_range_mesh=new TH2D("Energy Range Mesh","Energy Range Mesh",70,10,1000,40,0,15);
  TH1F* ehist = new TH1F("ehist",";E_{anode} [keVee];N [arb. units]",100,0,1000);
  TH1F* eh1 = new TH1F("eh1","",50,0,1000);
  TH1F* effHist = new TH1F("effHist",";E_{anode} [keVee];Efficiency",50,0,1000);
  TH1F* effHist2 = new TH1F("effHist2",";E_{anode} [keVee];Efficiency",50,0,1000);
  TH1D *hist_cutnum = new TH1D("Cut number","Cut number", 50, 0, 50);
  TFile *histfile=new TFile(Form("outtree_%d_%d.root",firstrun,lastrun), "RECREATE");
  tree = new TTree("dctpc", "Event info");
  tree->Branch("RunNum", &runnum, "runnum/I");
  tree->Branch("EventNum", &evnum, "evnum/I");
  tree->Branch("Spark", &spark, "spark/I");
  tree->Branch("PassTrack", &passTrack, "passTrack/I");
  tree->Branch("PassTrig", &passTrig, "passTrig/I");
  tree->Branch("PassAll", &passAll, "passAll/I");
  tree->Branch("Ntrig", &ntrig, "ntrig/I");
  tree->Branch("Etrack_kev", &Etrack, "Etrack/D");
  tree->Branch("Etrig_kev", &Etrig, "Etrig/D"); 
  tree->Branch("Emesh_kev", &Emesh, "Emesh/D");
  tree->Branch("deltaE_kev", &bestDiff, "bestDiff/D");
  tree->Branch("Phi_deg", &phi, "phi/D"); 
  tree->Branch("TrackX_pix", &TrackX, "TrackX/D"); 
  tree->Branch("TrackY_pix", &TrackY, "TrackY/D");
  tree->Branch("Rangeccd_mm", &range_ccd, "range_ccd/D");
  tree->Branch("r2_mm2", &r2, "r2/D");
  tree->Branch("exposure_s", &exposure, "exposure/D");
  tree->Branch("cutnum", &cutnum, "cutnum/I");
  tree->Branch("anodeRMS", &anodeRMS, "anodeRMS/D");
  tree->Branch("meshRMS", &meshRMS, "meshRMS/D");
  tree->Branch("timenow",&timenow, "timenow/I");
  //This bit attaches/connects certain variables to "branches" of the generalized tree


  ////////////////////////////////////////File Name Creation///////////////////////////////////
  for (int x = firstrun; x <= lastrun; x++)
    {
      runnum=x;
      if (x == 10002 || x==10003) {continue;}
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
      /////////////////////////////////////////////////////////////////////////////////////////////


      //////////////////////////////Creation of objects from Skim/Data Files///////////////////////       
      DmtpcSkimDataset d;  
      d.openRootFile(skimfilename.c_str());     
      d.loadDmtpcEvent(true,origfilename.c_str());
      Long_t nEvents=d.tree()->GetEntriesFast();
      //above prepares the "objects" so they are ready for analysis 

      //start event loop 
      for(int i=0; i<nEvents; i++){    
	if(i%100==0)
	  cout<<"Event: "<<i<<" / "<<nEvents<<endl;
    
	d.getEvent(i);
	DmtpcSkimEvent* ev=d.event(); //creating an event object
	evnum=i;
	Etrack=0.;
	Etrig=0.;
	Emesh=0.;
	passTrack=0;
	passTrig=0;
	//cutnum=0;
	passAll=0;
	spark=0;
	exposure=0;
	ntrig=0;
	r2=0.;
	TrackX=0.;
	TrackY=0.;
	range_ccd=0.;
	anodeRMS=0.;
	meshRMS=0.;
	timenow = ev->timenow(0);
	//cout << timenow << endl;
	//Above instantiates all needed variabes

	if(ev->spark(0)){ lastSpark = ev->eventNumber();spark=1;cutnum = 7; hist_cutnum->Fill(cutnum); tree->Fill();continue;} 
	//is it a spark? if so update the tree and move on to the next
	if(ev-lastSpark<6){spark=1;tree->Fill(); continue;} // make sure that at least 6 events have passed after a spark before considering the data
	totalevents++;
	exposure=1;

	if(ev->ntracks(0)==0){tree->Fill(); continue;}


	range_ccd=0.; 
	bestDiff=1e9;

	//loop over tracks in each event
	for(int k=0;k<ev->ntracks(0);k++)
	  {	    

	    if (ev->E(0,k)<1){cutnum=40; hist_cutnum->Fill(cutnum); continue;} //Cut if E is below 1 (measured in ADU)
	    if (ev->nburnin(0,k)){cutnum = 41; hist_cutnum->Fill(cutnum); continue;} //Cut if there is a burnin
	    if (ev->edge(0,k)){cutnum = 50; hist_cutnum->Fill(cutnum); continue;} //Cut if track reaches the edge
	    if (ev->range(0,k)<=0){cutnum=42; hist_cutnum->Fill(cutnum); continue;} //Cut if it has a zero range
	    if (ev->range(0,k)>160){cutnum=43; hist_cutnum->Fill(cutnum); continue;} //also cut if it is too big! 
	    if (ev->maxpixel(0,k)>600){cutnum=44; hist_cutnum->Fill(cutnum); continue;} //maxpixel has to be less than 250ADU
	    if (ev->x(0,k)<48) continue; //this and the following->
	    if (ev->y(0,k)<48) continue; //instructions make sure that->
	    if (ev->x(0,k)>976) continue; //the track does not cross->
	    if (ev->y(0,k)>976){cutnum = 48; hist_cutnum->Fill(cutnum); continue;} //the edge
	    if (ev->maxpixel(0,k)/ev->E(0,k)>0.25){cutnum=49; hist_cutnum->Fill(cutnum); continue;} //To eliminate noise artifacts the max pixel has to be less than 25% of the total signal 
	      
    	      	    
	    xx = ev->x(0,k)-512;
	    yy = ev->y(0,k)-512;      
	    //if (yy >760-512&&yy<780-512) continue;//what is this?
	    r2 = (xx*xx+yy*yy)* 16.66/1024*16.66/1024; 
        
	    passTrig=0;
	    passTrack++;
         
	    if(passTrack==1)
	      {

		E1 = ev->E(0,k)*lightCalib;
		//if (useRadial) E1 *= pow(1+ r2/(16.123*16.123),3.5857/2 ); 
		Etrack = E1;
		phi = ev->phi(0,k);
		phi = atan2(sin(phi),cos(phi)) * 180./ 3.1416;   
	  	TrackX = xx;
		TrackY = yy;
		range_ccd=ev->range(0,k)*mmperpixel;
		eh1->Fill(Etrack);
	      }
	    else continue;
             
	    TObjArray* arr = ev->waveform_vectors();
	    FastWfVector* anode = (FastWfVector*) arr->At(0);
	    CspWfVector* mesh = (CspWfVector*) arr->At(1);
	    CspWfVector* veto = (CspWfVector*) arr->At(2);
	     	  
	    //loop over triggers
	    ntrig=mesh->size();
	    
	    for (int nt = 0; nt < mesh->size(); nt++){
  //cout<<cutnum<<endl;
	      diff=1e9;
	      anodeRMS=mesh->at(nt).getRMS();
	      meshRMS=anode->at(nt).getRMS();
	      // if (mesh->at(nt).getRMS() >0.5){cutnum=1; continue;}//for neutron calib run: >0.002
// 	      if (anode->at(nt).getRMS()>40.){cutnum=2; continue;}//for neutron calib run: >0.003
// 	      if (fabs(mesh->at(nt,0).getFastRise0())>0.5){cutnum=3; continue;}
// 	      if (fabs(mesh->at(nt,0).getFastRise10())>0.5){cutnum=4; continue;}
// 	      if (fabs(mesh->at(nt,0).getFastRise25())>0.5){cutnum=5; continue;}
// 	      if (fabs(mesh->at(nt,0).getFastRise50())>0.5){cutnum=6; continue;}
// 	      if (fabs(mesh->at(nt,0).getFastRise75())>0.5){cutnum=7; continue;}  
// 	      if (fabs(mesh->at(nt,0).getFastRise90())>0.5){cutnum=8; continue;}
// 	      if (fabs(mesh->at(nt,0).getRise0())>0.5){cutnum=9; continue;}
// 	      if (fabs(mesh->at(nt,0).getRise10())>0.5){cutnum=10; continue;}
// 	      if (fabs(mesh->at(nt,0).getRise25())>0.5){cutnum=11; continue;}
// 	      if (fabs(mesh->at(nt,0).getRise50())>0.5){cutnum=12; continue;}
// 	      if (fabs(mesh->at(nt,0).getRise75())>0.5){cutnum=13; continue;}
// 	      if (fabs(mesh->at(nt,0).getRise90())>0.5){cutnum=14; continue;}
// 	      if (fabs(mesh->at(nt,0).getSlowFall0())>0.5){cutnum=15; continue;}
// 	      if (fabs(mesh->at(nt,0).getSlowFall10())>0.5){cutnum=16; continue;}
// 	      if (fabs(mesh->at(nt,0).getSlowFall25())>0.5){cutnum=17; continue;}
// 	      if (fabs(mesh->at(nt,0).getSlowFall50())>0.5){cutnum=18; continue;}
// 	      if (fabs(mesh->at(nt,0).getSlowFall75())>0.5){cutnum=19; continue;}
// 	      if (fabs(mesh->at(nt,0).getSlowFall90())>0.5){cutnum=20; continue;}
// 	      if (fabs(mesh->at(nt,0).getFall0())>0.5){cutnum=21; continue;}
// 	      if (fabs(mesh->at(nt,0).getFall10())>0.5){cutnum=22; continue;}
// 	      if (fabs(mesh->at(nt,0).getFall25())>0.5){cutnum=23; continue;}
// 	      if (fabs(mesh->at(nt,0).getFall50())>0.5){cutnum=24; continue;}
// 	      if (fabs(mesh->at(nt,0).getFall75())>0.5){cutnum=25; continue;}
// 	      if (fabs(mesh->at(nt,0).getFall90())>0.5){cutnum=26; continue;}  
// 	      if (mesh->at(nt).getWfMax() >0.39 || mesh->at(nt).getWfMin() < -0.39){cutnum=27; continue;}
// 	      if (anode->at(nt).getWfMax() >0.195 || anode->at(nt).getWfMin() < -0.195){cutnum=28; continue;}
// 	      if (veto->at(nt).getWfMax() >0.195 || veto->at(nt).getWfMin() < -0.195){cutnum=29; continue;}
// 	      if (mesh->at(nt,0).getIntegral()<0.5){cutnum=30; continue;}
// 	      if (anode->at(nt,0).getPeak()<0.025){cutnum=31; continue;}
// 	      if (anode->at(nt,0).getPeak() < 3 * veto->at(nt,0).getPeak()){cutnum=32; hist_cutnum->Fill(cutnum); continue;}//r_veto>
// 	      if (veto->at(nt,0).getRise10() - veto->at(nt,0).getRise90()>800e-9){cutnum=33; hist_cutnum->Fill(cutnum); continue;}//t_veto
// 	      //if (mesh->at(nt,0).getPeak() > 3 * anode->at(nt,0).getPeak()){cutnum=34; hist_cutnum->Fill(cutnum); continue;} //r_mesh < 3
// 	      //if (mesh->at(nt,0).getPeak() < 1.3 * anode->at(nt,0).getPeak()){cutnum=35; hist_cutnum->Fill(cutnum); continue;} //r_mesh > 1.3
// 	      if (anode->at(nt,0).getRise10() - anode->at(nt,0).getRise90()>1.25e-6){cutnum=36; hist_cutnum->Fill(cutnum); continue;}//t_anode<1250ns
        
	      passTrig++;//increment number of tracks that pass the trigger cuts
    	    
	      E2 = (anode->at(nt,0).getPeak()*1000)*anodeCalib;//*chargeCalib //creation of E_anode
	      E3 = (mesh->at(nt,0).getPeak()*1000)*meshCalib;//*meshCalib
	      //diff = E1 - fitPoly12Percent(E2);
	      diff = E1 - E2;
	      //cout <<"E1: "<< E1<< endl << "E2: " << E2 << endl << "E3: " << E3 << endl << "diff: " << diff << endl;
	      if (fabs(diff)<fabs(bestDiff)) {
		bestDiff = diff;
		Etrig = E2;
		Emesh = E3;
	      }
	      //cout <<"Event: "<<i<< " Etrack: " << Etrack << endl << "Etrig: " << Etrig << endl << "Emesh: "<<Emesh<< endl <<"bestdiff: " << bestDiff<< endl;
	    } // end trigger loop 
	    ntotalTrig+=passTrig;

	  }//end track loop

	ntotalTrack+=passTrack;
	
	if (passTrack!=1){tree->Fill(); continue;} //neutron events are unlikely to have more than 1 track
	//	cout << "here0" <<endl;
	if (passTrig==0){tree->Fill(); continue;} 
	//	cout << "here1" <<endl;
	//if (Etrig < 250 && bestDiff > 75){cutnum=37;hist_cutnum->Fill(cutnum);tree->Fill(); continue;}
	//	cout << "here2" <<endl;	
	//if (Etrig>=250&&Etrig<500&& bestDiff>125){cutnum=38; hist_cutnum->Fill(cutnum); tree->Fill(); continue;}
	//	cout << "here3" <<endl;	
	//if (Etrig>=500&& bestDiff>150){cutnum=39; hist_cutnum->Fill(cutnum); tree->Fill(); continue;}
	//	cout << "here4" <<endl;        
	if(Etrack<10 && Etrack<10){tree->Fill(); continue;}
	passAll=1;
	//	cout << "here5" <<endl;
	effHist->Fill(Etrack); //fill the effHist if all cuts are passed.
        
	  
	//cout<<i<<" Energies: "<<Etrack<<" "<<Etrig<<" "<<phi<<endl;
	hist_energy->Fill(Etrack,Etrig);
	hist_energy2->Fill(Etrig,Emesh);
	hist_energy3->Fill(Etrack,Emesh);
	hist_energy_anode->Fill(Etrig);
	hist_energy_mesh->Fill(Emesh);
	hist_energy_ccd->Fill(Etrack);
	hist_diff->Fill(bestDiff);      
	hist_energy_range_ccd->Fill(Etrack,range_ccd);  
	hist_energy_range_anode->Fill(Etrig,range_ccd);  
	hist_energy_range_mesh->Fill(Emesh,range_ccd);
	hist_phi->Fill(abs(phi)); 
	hist_cosphi->Fill(cos(phi));
	hist_pos->Fill(TrackX,TrackY); 
	//	cout << "Anode: " << Etrig << endl << "CCD: " << Etrack << endl;
	eventnums.push_back(i);
	runnums.push_back(x);
	  
		
	tree->Fill();	
	
      }//end event loop

    }//end file loop

  for (int i = 1; i<=effHist->GetNbinsX(); i++)
    {
      double n1 = eh1->GetBinContent(i);
      double n2 = effHist->GetBinContent(i);
    
      if (!n1) continue;
      effHist2->SetBinContent(i,n2/n1);
      effHist2->SetBinError(i,sqrt(n2/n1 * (1-n2/n1)/n1));
    }

  cout <<"Total seconds of exposure: "<<totalevents<<endl;
 
 TLine line;
line.SetLineWidth(.1);
line.SetLineStyle(2);
 
  new TCanvas; 
  hist_energy->SetTitle("CCD Energy vs. Anode Energy");
  hist_energy->SetXTitle("E_{CCD} ");
  hist_energy->SetYTitle("E_{anode} ");
  hist_energy->Draw("COLZ");
  line.DrawLine(0,0,1000,1000);
  
  new TCanvas; 
  hist_energy2->SetTitle("Anode Energy vs. Mesh Energy");
  hist_energy2->SetXTitle("E_{anode} ");
  hist_energy2->SetYTitle("E_{mesh} ");
  hist_energy2->Draw("COLZ");
  line.DrawLine(0,0,1000,1000);
  
  new TCanvas; 
  hist_energy3->SetTitle("CCD Energy vs. Mesh Energy");
  hist_energy3->SetXTitle("E_{CCD} ");
  hist_energy3->SetYTitle("E_{mesh} ");
  hist_energy3->Draw("COLZ");
  line.DrawLine(0,0,1000,1000);


  new TCanvas;
  hist_cutnum->SetTitle("Cut Number");
  hist_cutnum->SetXTitle("Cut Number");
  hist_cutnum->Draw();
  for (i=1; i<51; i++)
    {cout << "Cut" << i -1 << ": "<< hist_cutnum->GetBinContent(i) << "\n";}
  new TCanvas; 
  hist_energy_range_ccd->SetTitle("Range vs. CCD Energy");
  hist_energy_range_ccd->SetXTitle("E_{CCD} ");
  hist_energy_range_ccd->SetYTitle("2D range");
  hist_energy_range_ccd->Draw("COLZ");
  
  new TCanvas; 
  hist_energy_range_anode->SetTitle("Range vs. Anode Energy");
  hist_energy_range_anode->SetXTitle("E_{anode} ");
  hist_energy_range_anode->SetYTitle("2D range ");
  hist_energy_range_anode->Draw("COLZ");
  
    new TCanvas; 
  hist_energy_range_mesh->SetTitle("Range vs. Mesh Energy");
  hist_energy_range_mesh->SetXTitle("E_{mesh}");
  hist_energy_range_mesh->SetYTitle("2D range");
  hist_energy_range_mesh->Draw("COLZ");

  new TCanvas;
  hist_energy_anode->Draw();
  hist_energy_anode->SetTitle("Anode Energy");
  hist_energy_anode->SetXTitle("E_{anode}");
 
  new TCanvas;
  hist_energy_mesh->Draw();
  hist_energy_mesh->SetTitle("Mesh Energy");
  hist_energy_mesh->SetXTitle("E_{mesh} ");	
 
  new TCanvas;
  hist_energy_ccd->Draw(); 
  hist_energy_ccd->SetTitle("CCD Energy");
  hist_energy_ccd->SetXTitle("E_{CCD} ");

  new TCanvas;
  hist_diff->Draw();
  hist_diff->SetTitle("CCD and Anode #Delta_{energy}");
  hist_diff->SetXTitle("E_{CCD} - E_{anode} ");
  hist_diff->Fit("gaus");


  new TCanvas;
  hist_phi->Draw(); 
  hist_phi->SetTitle("Track angle");
  hist_phi->SetXTitle("|#phi| (degrees)");

  new TCanvas;
  hist_cosphi->Draw(); 
  hist_cosphi->SetTitle("Track angle");
  hist_cosphi->SetXTitle("cos(#phi)");	

  new TCanvas;
  hist_pos->Draw("colz"); 
  hist_pos->SetTitle("Track position");
  hist_pos->SetXTitle("X (pixels)");
  hist_pos->SetYTitle("Y (pixels)");

 //  new TCanvas;
//   effHist2->Draw("e1");

  for(int k = 0; k < eventnums.size(); k++){
    cout <<runnums.at(k)<<": "<< eventnums.at(k) <<endl; }

  // hist_energy->Write();
  // hist_energy_range_ccd->Write();
  // hist_energy_anode->Write();
  // hist_energy_ccd->Write();
  // hist_diff->Write();
  // hist_phi->Write();
  // hist_pos->Write();
  // effHist->Write();
  // eh1->Write();
  histfile = tree->GetCurrentFile();
  histfile->Write();
  histfile->Close();
  
  gApplication->Terminate(); 
}

double fitPoly12Percent(double E){
  return sqrt(0.9741*0.9741*E*E+113.835*113.835)-62.675;
  //return -11.50+39.02*E+0.03459*E*E;
  //return 469 + 10*E + 0.548 *E*E -3.39e-3*E*E*E+7.5e-6*E*E*E*E;
}

