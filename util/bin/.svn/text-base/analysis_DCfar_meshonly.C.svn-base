void analysis(int firstrun, int lastrun){//inclusive 

// CF252 (neutron) runs: 1041-1158
// Note that the 'skimmed' runs need to exist before this macro will run
//CF252 runs w/ 100% CF4 (878-1040, inclusive)
// To run over all CF252 runs w/ 12.5% He (1041-1158, inclusive): 
// root
// .L analysis_DCfar_meshonly.C
// analysis(1040, 1158) 

//For running at dcnode, the input file location needs to be changed

//The program creates some relevant histograms and a TTree. The TTree is only relevant if this program is run in parallel for processing a large number of runs. 
// All of the cuts that are utilized are explained in arXiv:1108.4894

//double lightCalib = 9.19;//light calibration factor  (12.5%, default)
//double chargeCalib = 823./4440;//anode calibration factor  (12.5%, default)

double lightCalib = 72235./4440;//light calibration factor (CF4)
double chargeCalib = 1635./4440;//mesh calibration factor (CF4), just a guess right now
double mmperpixel = 166.6/1024;
bool useRadial=true; 
gROOT->SetStyle("Plain");
gStyle->SetPalette(1,0); 
gSystem->Load("libMaxCam");
gSystem->Load("libWaveformTools.so");
std::stringstream sstm;
vector<double> noiseVector;
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
double xx, yy, diff, E1, E2, phi;
double range_ccd=-1.;
double Etrack=-1.;
double Etrig=-1.;
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
TH1D *hist_energy_anode=new TH1D("Mesh energy","Mesh energy",50,10,1000);
TH1D *hist_phi=new TH1D("phi","phi",50,-180.,180.);
TH2D *hist_pos=new TH2D("pos","pos",128,-512,512,128,-512,512);
TH2D *hist_energy=new TH2D("Energy","Energy",50,10,1000,50,0,1000);
TH2D *hist_energy_range_ccd=new TH2D("Energy Range CCD","Energy Range CCD",50,10,1000,50,0,20);
TH1F* ehist = new TH1F("ehist",";E_{anode} [keVee];N [arb. units]",100,0,1000);
TH1F* eh1 = new TH1F("eh1","",50,0,1000);
TH1F* effHist = new TH1F("effHist",";E_{anode} [keVee];Efficiency",50,0,1000);
TH1F* effHist2 = new TH1F("effHist2",";E_{anode} [keVee];Efficiency",50,0,1000);
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

for (int x = firstrun; x <= lastrun; x++)
{
runnum=x;
string origfile = "/net/nudsk0001/d00/scratch/dctpc_tmp/data/dmtpc_DC_";//change me if running at DCNODE
string skimfile = "/net/nudsk0001/d00/scratch/dctpc_tmp/skim/dmtpc_DC_";//change me if running at DCNODE
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

	//start event loop 
    for(int i=0; i<nEvents; i++){    
    // if(i%100==0)
//     cout<<"Event: "<<i<<" / "<<nEvents<<endl;
    
	d.getEvent(i);
	DmtpcSkimEvent* ev=d.event();
	evnum=i;
	Etrack=0.;
	Etrig=0.;
	passTrack=0;
	passTrig=0;
	cutnum=0;
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
	
	if(ev->spark(0)){ lastSpark = ev->eventNumber();spark=1;tree->Fill();continue;}
    if(ev-lastSpark<6){spark=1;tree->Fill(); continue;}
    totalevents++;
    exposure=1;

	if(ev->ntracks(0)==0){tree->Fill(); continue;}

	range_ccd=0.; 
	bestDiff=1e9;

		//loop over tracks in each event
		for(int k=0;k<ev->ntracks(0);k++)
	    {	    
		if (ev->E(0,k)<1) continue;
		if (ev->nburnin(0,k)) continue;
		if (ev->edge(0,k)) continue;
		if (ev->range(0,k)<=0) continue;
		if (ev->range(0,k)>160) continue;
		if (ev->maxpixel(0,k)>250) continue;
		if (ev->x(0,k)<50) continue;//normally < 24
		if (ev->y(0,k)<50) continue;//normally < 24
		if (ev->x(0,k)>974) continue;//normally > 1000
		if (ev->y(0,k)>974) continue;//normally > 1000
		if (ev->maxpixel(0,k)/ev->E(0,k)>0.25) continue;
		if(ev->pixels_killed(0)==0) continue;//added by spitz. gets rid of sparks.
    	if(ev->image_mean(0)>1005.) continue;//added by spitz. gets rid of sparks
    	      	    
		xx = ev->x(0,k)-512;
        yy = ev->y(0,k)-512;      
        //if (yy >760-512&&yy<780-512) continue;//what is this?
        r2 = (xx*xx+yy*yy)* 16.66/1024*16.66/1024; 
        
        passTrig=0;
        passTrack++;
         
		if(passTrack==1)
		{
		E1 = ev->E(0,k)/lightCalib;
// 		cout<<x<<" "<<i<<endl;
		if (useRadial) E1 *= pow(1+ r2/(16.123*16.123),3.5857/2 ); 
		Etrack = E1;
		phi = ev->phi(0,k);
        phi = atan2(sin(phi),cos(phi)) * 180./ 3.1416;   
	  	TrackX = xx;
		TrackY = yy;
		range_ccd=ev->range(0,k)*mmperpixel;
		eh1->Fill(Etrack);
		}
		else
		{
		cutnum=-1;
		continue;
		}
             
		TObjArray* arr = ev->waveform_vectors();
		FastWfVector* mesh = (FastWfVector*) arr->At(0);
		CspWfVector* anode = (CspWfVector*) arr->At(1);
		CspWfVector* veto = (CspWfVector*) arr->At(2);
	     	  
	    //loop over triggers
	    ntrig=mesh->size();
	    
		for (int nt = 0; nt < mesh->size(); nt++){	  
		diff=1e9;
		anodeRMS=mesh->at(nt).getRMS();
		meshRMS=anode->at(nt).getRMS();
		//cout<<mesh->at(nt).getRMS()<<endl;
	    if (mesh->at(nt).getRMS() >0.009){cutnum=1; continue;}//for neutron calib run: >0.002
	    //if (anode->at(nt).getRMS()>0.0006){cutnum=2; continue;}//for neutron calib run: >0.003
	    if (fabs(mesh->at(nt,0).getFastRise0())>0.5){cutnum=3; continue;}
	    if (fabs(mesh->at(nt,0).getFastRise10())>0.5){cutnum=4; continue;}
	    if (fabs(mesh->at(nt,0).getFastRise25())>0.5){cutnum=5; continue;}
	    if (fabs(mesh->at(nt,0).getFastRise50())>0.5){cutnum=6; continue;}
	    if (fabs(mesh->at(nt,0).getFastRise75())>0.5){cutnum=7; continue;}  
	    if (fabs(mesh->at(nt,0).getFastRise90())>0.5){cutnum=8; continue;}
	    if (fabs(mesh->at(nt,0).getRise0())>0.5){cutnum=9; continue;}
	    if (fabs(mesh->at(nt,0).getRise10())>0.5){cutnum=10; continue;}
	    if (fabs(mesh->at(nt,0).getRise25())>0.5){cutnum=11; continue;}
	    if (fabs(mesh->at(nt,0).getRise50())>0.5){cutnum=12; continue;}
	    if (fabs(mesh->at(nt,0).getRise75())>0.5){cutnum=13; continue;}
	    if (fabs(mesh->at(nt,0).getRise90())>0.5){cutnum=14; continue;}
	    if (fabs(mesh->at(nt,0).getSlowFall0())>0.5){cutnum=15; continue;}
	    if (fabs(mesh->at(nt,0).getSlowFall10())>0.5){cutnum=16; continue;}
	    if (fabs(mesh->at(nt,0).getSlowFall25())>0.5){cutnum=17; continue;}
	    if (fabs(mesh->at(nt,0).getSlowFall50())>0.5){cutnum=18; continue;}
	    if (fabs(mesh->at(nt,0).getSlowFall75())>0.5){cutnum=19; continue;}
	    if (fabs(mesh->at(nt,0).getSlowFall90())>0.5){cutnum=20; continue;}
	    if (fabs(mesh->at(nt,0).getFall0())>0.5){cutnum=21; continue;}
	    if (fabs(mesh->at(nt,0).getFall10())>0.5){cutnum=22; continue;}
	    if (fabs(mesh->at(nt,0).getFall25())>0.5){cutnum=23; continue;}
	    if (fabs(mesh->at(nt,0).getFall50())>0.5){cutnum=24; continue;}
	    if (fabs(mesh->at(nt,0).getFall75())>0.5){cutnum=25; continue;}
	    if (fabs(mesh->at(nt,0).getFall90())>0.5){cutnum=26; continue;}  
	    //spitz comments this out:
	    //if (mesh->at(nt).getWfMax() >0.39 || mesh->at(nt).getWfMin() < -0.39){cutnum=27; continue;}
	    
	    //if (anode->at(nt).getWfMax() >0.195 || anode->at(nt).getWfMin() < -0.195){cutnum=28; continue;}
	    //if (veto->at(nt).getWfMax() >0.195 || veto->at(nt).getWfMin() < -0.195){cutnum=29; continue;}
	    //spitz default is 0.5:
	    if (mesh->at(nt,0).getIntegral()<0.5){cutnum=30; continue;}
	    
	    //if (anode->at(nt,0).getPeak()<0.01){cutnum=31; continue;}
	    //if (anode->at(nt,0).getPeak() < 3 * veto->at(nt,0).getPeak()){cutnum=32; continue;}
	    //if (veto->at(nt,0).getRise10() - veto->at(nt,0).getRise90()>800e-9){cutnum=33; continue;}
	   // if (mesh->at(nt,0).getPeak() > 3 * anode->at(nt,0).getPeak()){cutnum=34;continue;}
	    //if (mesh->at(nt,0).getPeak() < 1.3 * anode->at(nt,0).getPeak()){cutnum=35; continue;}
	    //if (anode->at(nt,0).getRise10() - anode->at(nt,0).getRise90()>1.25e-6){cutnum=36; continue;}

        passTrig++;    	    
		//E2 = (anode->at(nt,0).getPeak()*1000)/chargeCalib;//default, based on anode
		E2 = (mesh->at(nt,0).getPeak()*1000)/chargeCalib;
		//diff = E1 - fitPoly12Percent(E2); //default
		diff = E1 - E2;
		//cout<<"Pass trig!"<<endl;
			if (fabs(diff)<fabs(bestDiff)) {
			bestDiff = diff;
			Etrig = E2;
			}
	  
		} // end trigger loop 
		//cout<<"Cut num: "<<cutnum<<" "<<x<<" "<<i<<endl;
		}//end track loop

	ntotalTrig+=passTrig;
	ntotalTrack+=passTrack;

	if (passTrack!=1){tree->Fill(); continue;}
    if (passTrig==0){tree->Fill(); continue;}  
    // cout<<x<<" "<<i<<endl;
//     cout<<"Etrack: "<<Etrack<<" Etrig: "<<Etrig<<endl;
      
 //    if (Etrig < 250 && abs(bestDiff) > 75){cutnum=37;tree->Fill(); continue;}
//     if (Etrig>=250&&Etrig<500&& abs(bestDiff)>125){cutnum=38;tree->Fill(); continue;}
//     if (Etrig>=500&& abs(bestDiff)>150){cutnum=39;tree->Fill(); continue;}
    passAll=1;
    effHist->Fill(Etrack);
        
		if(Etrack>0 && Etrig>0)
		{
		
		cout<<x<<" "<<i<<" Energies: "<<Etrack<<" "<<Etrig<<" "<<phi<<endl;
		hist_energy->Fill(Etrack,Etrig);   
		hist_energy_anode->Fill(Etrig);
		hist_energy_ccd->Fill(Etrack);
		hist_diff->Fill(bestDiff);      
		hist_energy_range_ccd->Fill(Etrack,range_ccd);  
		hist_phi->Fill(phi); 
		hist_pos->Fill(TrackX,TrackY); 
		}
		
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
 
new TCanvas; 
hist_energy->SetTitle("CCD Energy vs. Mesh Energy");
hist_energy->SetXTitle("E_{CCD} (keV)");
hist_energy->SetYTitle("E_{mesh} (keV)");
hist_energy->Draw("COLZ");

new TCanvas; 
hist_energy_range_ccd->SetTitle("Range vs. CCD Energy");
hist_energy_range_ccd->SetXTitle("E_{CCD} (keV)");
hist_energy_range_ccd->SetYTitle("2D range (mm)");
hist_energy_range_ccd->Draw("COLZ");

new TCanvas;
hist_energy_anode->Draw();
hist_energy_anode->SetTitle("Mesh Energy");
hist_energy_anode->SetXTitle("E_{mesh} (keV)");
 
new TCanvas;
hist_energy_ccd->Draw(); 
hist_energy_ccd->SetTitle("CCD Energy");
hist_energy_ccd->SetXTitle("E_{CCD} (keV)");

new TCanvas;
hist_diff->Draw();
hist_diff->SetTitle("CCD and Mesh #Delta_{energy}");
hist_diff->SetXTitle("E_{CCD} - E_{mesh} (keV)");
hist_diff->Fit("gaus");

new TCanvas;
hist_phi->Draw(); 
hist_phi->SetTitle("Track angle");
hist_phi->SetXTitle("#phi (degrees)");

new TCanvas;
hist_pos->Draw("colz"); 
hist_pos->SetTitle("Track position");
hist_pos->SetXTitle("X (pixels)");
hist_pos->SetYTitle("Y (pixels)");

new TCanvas;
effHist2->Draw("e1");

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
  
//gApplication->Terminate(); 
}

double fitPoly12Percent(double E){
  return sqrt(0.9741*0.9741*E*E+113.835*113.835)-62.675;
  //return -11.50+39.02*E+0.03459*E*E;
  //return 469 + 10*E + 0.548 *E*E -3.39e-3*E*E*E+7.5e-6*E*E*E*E;
}

