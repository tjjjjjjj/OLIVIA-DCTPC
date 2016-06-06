///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
////Use this file to simulate alpha and neutron tracks               //
////To change projectile distributions:                              //
////Look at generateAlphas and generateNeutrons                      //
////                                                                 //
////Questions/suggestions? jplopez@mit.edu                           //
////01/20/09                                                         //
///////////////////////////////////////////////////////////////////////

#include "../MaxCamMC.hh"
#include "../MaxCamImageTools.hh"
#include "../DmtpcEvent.hh"
#include "../DmtpcDataset.hh"
#include "../MaxCamChannel.hh"
#include "../MaxCamConfig.hh"

#include "TH2.h"
#include "TSystem.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TVector3.h"
#include "TObjString.h"
#include "TROOT.h"
#include "TDatime.h"
#include "TTimeStamp.h"

#include "makeSimTools.cc"
#include "SimGenerator.hh"

#include <vector>
#include <iostream>
#include <fstream>
#include <assert.h>
using namespace std;

MaxCamMC *mc = 0;
TTree *sim;
DmtpcDataset *dataset;
float simx,simy,simz,simlength,simEnergy,simEscint;
float simmass, simphi,simtheta,simprojE,simIntegral;
float simprojtheta,simprojphi,simprojmass;
TDatime *simtime;
int simeventnum;

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
////Init()///////////////////////////////////////////////////////////////////
////Initialize most MaxCamMC parameters                                    //
////Create trees for holding images, simulation info                       //
////                                                                       //
////Default Parameters:                                                    //
////Gain            :   17.0                                               //
////Noise           :   12.0                                               //
////Width           :   160                                                //
////Pressure        :   75                                                 //
////Spacer Diameter :   0.5                                                //
////Spacer Spacing  :   20.0                                               //
////Diff. Cont.     :   0.2704                                             //
////Diff. Dz        :   0.0076                                             //
////Things to add: anodeHV, meshHV: not really needed                      //
/////////////////////////////////////////////////////////////////////////////

void Init(float gain, float noise,TString recoil = "F",TString gas = "CF4_100Torr",
	  int npix = 1024,int nbins = 256, float width = 0.143*1024, 
	  float pressure = 75,float spacerdiam = 0.5, 
	  float spacerspacing = 20.0,float diffconst = 0.2704,
	  float diffdz = 0.0076){

  /////////////////////////////////////////////////////////////////////////////
  ////Parameters://////////////////////////////////////////////////////////////
  ////Gain          :   ADC / per keV                                        //
  ////Noise         :   ADC (counts)                                         //
  ////filename      :   .root file                                           // 
  ////npix          :   # of pixels in one dimension (assume square CCD)     //
  ////nbins         :   # of bins in one dimension                           //
  ////width         :   roi size in mm                                       //
  ////pressure      :   in torr                                              //
  ////spacerdiam    :   diameter of spacers, in mm                           //
  ////spacerspacing :   separation between spacers, in mm                    //
  ////diffconst     :   Diffusion constant term                              //
  ////diffdz        :   Diffusion dz term                                    //
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  gSystem->Setenv("MCTABLES","../tables");
  TString fname = "mcrun.root";
  mc = new MaxCamMC(fname);
  gRandom = new TRandom3(0);
  //Set image histograms
  TString filename = "SRIM_" + recoil + "_in_" + gas;
  ifstream instr("../tables/" + filename);
  cout <<"Filling SRIM table "<< filename<<endl;
  if (!instr.is_open()){
    cout << "Recoil type not found!"<<endl;
    assert(0);
  }

  mc->fillSrimTable(filename);
 
  mc->setTrackImage(nbins,0.,width,nbins,0.,width);//Track Image
  mc->setWireImage(nbins,0.,width,nbins,0.,width);//Anode Image
  mc->setCCDImage(nbins,0,npix,nbins,0,npix);      //CCD Image

  //Set parameters
  mc->setPressure(pressure);
  mc->setPhotonsPerkeV(gain);
  mc->setPixelsPermm(nbins/width);
  mc->setNuclScintillation(0.3);
  mc->setNoiseADC(noise);
  mc->setDiffusionConstTerm(diffconst);
  mc->setDiffusionDzTerm(diffdz);

  //Add spacers
  mc->clearWireList();
  mc->clearSpacerList();
  mc->setSpacerDiameter(spacerdiam);
  int n_spacers =(int)(width/spacerspacing + 1);
  for (int i = 0; i < n_spacers; i++) mc->addSpacer(i*spacerspacing);

}//end Init()

void dataSetup(float gain = 3.47, float noise = 7,TString filename = "example.root",
	  int npix = 1024,int nbins = 256, float width = 0.143*1024, 
	  float pressure = 75,float spacerdiam = 0.5, 
	  float spacerspacing = 20.0,float diffconst = 0.2704,
	  float diffdz = 0.0076){

  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  ////Set up data file                                                         //
  ///////////////////////////////////////////////////////////////////////////////

  dataset = new DmtpcDataset();
  dataset->createRootFile((const char*) filename,"recreate");
  dataset->setComment("Radial gain effects added 4/21/09");
  dataset->setLocation("26-456");
  dataset->setKeyword("Simulations");
  dataset->comment()->Write();
  dataset->keyword()->Write();
  dataset->location()->Write();
  //dataset->tree()->SetCircular(10);

  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  ////Set up experimentConfig values. (copied from DmtpcRun.cc::L44-74)        //
  ///////////////////////////////////////////////////////////////////////////////

  int nPars = 0;
  MaxCamChannel* driftHV   = new MaxCamChannel("driftHV",    "Drift Voltage", 0,  0);
  driftHV->currentValue = 3;
  driftHV->setValue = 3;
  new( (*dataset->event()->experimentConfig())[nPars++] ) MaxCamChannel(*driftHV);
  delete driftHV;
  
  MaxCamChannel* anodeHV   = new MaxCamChannel("anodeHV",    "Anode Voltage", 1,  1);
  anodeHV->currentValue = 0.68;
  anodeHV->setValue = 0.68;
  new( (*dataset->event()->experimentConfig())[nPars++] ) MaxCamChannel(*anodeHV);
  delete anodeHV;
  
  MaxCamChannel* driftI    = new MaxCamChannel("driftI",     "Drift current", 2, -1);
  driftI->currentValue = 0;
  driftI->setValue = 0;
  new( (*dataset->event()->experimentConfig())[nPars++] ) MaxCamChannel(*driftI);
  delete driftI;

  MaxCamChannel* temp0   = new MaxCamChannel("temp0",    "Temperature point 0", 4,  -1);
  temp0->currentValue = 0;
  temp0->setValue = 0;
  new( (*dataset->event()->experimentConfig())[nPars++] ) MaxCamChannel(*temp0);
  delete temp0;
  
  MaxCamChannel* anodeI    = new MaxCamChannel("anodeI",     "Anode current", 5, -1);
  anodeI->currentValue = 0;
  anodeI->setValue = 0;
  new( (*dataset->event()->experimentConfig())[nPars++] ) MaxCamChannel(*anodeI);
  delete anodeI;

  MaxCamChannel* pressure0 = new MaxCamChannel("pressure",   "Gas pressure",  -1, -1);    
  pressure0->currentValue = pressure;
  pressure0->setValue = pressure;
  new( (*dataset->event()->experimentConfig())[nPars++] ) MaxCamChannel(*pressure0);
  delete pressure0;

  //////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////
  ////Set up MaxCamConfig file for the dataset.                                   //
  //////////////////////////////////////////////////////////////////////////////////

  MaxCamConfig *config0 = new MaxCamConfig("ccdConfig","CCD Camera configuration");
  config0->cameraID = 1;
  config0->row_width = npix;
  config0->img_rows = npix;
  config0->hbin = npix / nbins;
  config0->vbin = npix / nbins;
  config0->ul_x = 0;
  config0->ul_y = 0;
  config0->lr_x = npix;
  config0->lr_y = npix;
  config0->CCDTemp = -20;
  config0->CCDTempSet = -20;
  config0->exposureTime = 1000;
  config0->frameType = 0;
  config0->nFlushes = -1;
  config0->bitDepth = 65535;
  config0->daqTime = -1;
  new( (*dataset->event()->ccdConfig())[0]) MaxCamConfig(*config0);
  delete config0;

  //////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////
  ////Create a 1 entry tree holding all of the input parameters.                  //
  //////////////////////////////////////////////////////////////////////////////////

  int par_1[2];
  float par_2[8];
  TTree *params = new TTree("params","params");
  params->Branch("par_1",&par_1,"Npix/I:Nbins",1000);
  params->Branch("par_2",&par_2,
		 "gain/F:ADCnoise:ImageWidth:Pressure:DiffConst:DiffDz",2000);

  par_1[0] = npix;
  par_1[1] = nbins;

  par_2[0] = gain;
  par_2[1] = noise;
  par_2[2] = width;
  par_2[3] = pressure;
  par_2[4] = diffconst;
  par_2[5] = diffdz;
  par_2[6] = spacerspacing;
  par_2[8] = spacerdiam;
  params->Fill();

  //////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////
  ////Create a tree to hold the simulated data for comparison w/ reconstruction   //
  //////////////////////////////////////////////////////////////////////////////////

  sim = new TTree("Simulation","Simulation");
  sim->Branch("mass",&simmass,"mass/F");
  sim->Branch("E",&simEnergy,"E/F");
  sim->Branch("Escint",&simEscint,"Escint/F");
  sim->Branch("phi",&simphi,"phi/F");
  sim->Branch("theta",&simtheta,"theta/F");
  sim->Branch("x",&simx,"x/F");
  sim->Branch("y",&simy,"y/F");
  sim->Branch("z",&simz,"z/F");
  sim->Branch("length",&simlength,"length/F");
  sim->Branch("time","TDatime",&simtime);
  sim->Branch("eventnum",&simeventnum,"eventnum/I");
  sim->Branch("projE",&simprojE,"projE/F");
  sim->Branch("Integral",&simIntegral,"Integral/F");
  sim->Branch("projPhi",&simprojphi,"projPhi/F");
  sim->Branch("projTheta",&simprojtheta,"projTheta/F");
  sim->Branch("projMass",&simprojmass,"projMass/F");
  // sim->SetCircular(100);
  //sim->Branch("simdata",&simdata,"mass/F:E:Escint:phi:theta:x:y:z:length:time",32000);

  //////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////
  ////Make bias frame for camera 1 with correct noise level                       //
  //////////////////////////////////////////////////////////////////////////////////

  TH2F *bias = new TH2F("biasFrame1","biasFrame1",nbins,0,npix,nbins,0,npix);
  for (int i = 1; i <= nbins; i++){
    for (int j = 1; j<= nbins;j++){
      
      bias->SetBinContent(i,j,gRandom->Gaus(0,noise));

    }
  }
  bias->Write();

}//End dataSetup()

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////Generate a blank image                                                        //
////////////////////////////////////////////////////////////////////////////////////
void generateEmptyEvent(){

  double noise = mc->getNoiseADC();
  int xbins = mc->getCCDImage()->GetNbinsX();
  int ybins = mc->getCCDImage()->GetNbinsY();
  //Gaussian noise only
  for (int i = 1; i<= xbins;i++){
    for (int j = 1; j<= ybins;j++){
      mc->getCCDImage()->SetBinContent(i,j,gRandom->Gaus(0,noise));
    }
  }
  simx = 0;
  simy = 0;
  simz = 0;
  simlength = 0;
  simEnergy = 0;
  simEscint = 0;
  simmass = 0;
  simphi = 0;
  simtheta = 0;
  simtime = new TDatime();
  simprojE = 0;
  simIntegral = 0;
  simprojtheta = 0;
  simprojphi = 0;
  simprojmass = 0;
  sim->Fill();
}//End generateBlankEvent()

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
////Generate events of various kinds                                               //
/////////////////////////////////////////////////////////////////////////////////////
bool Event(TString type, TString type2, double min, double max, double zmin = 0, double zmax =250,
	   double wimpMass = 100e6, double emin = 0){

  SimGenerator *simgen = new SimGenerator();
  simgen->setRecoilMass(19e6);
  simgen->setProjectileMass(wimpMass);//changed automatically if not wimp event
  simgen->setPositionMaxValues(max,max,zmax);
  simgen->setPositionMinValues(min,min,zmin);
  simgen->setMinEnergy(emin);
  //ar: random alpha
  //af: fixed position/direction alpha
  //cf: Cf-252 neutron
  //nr: random neutron
  //w:  WIMP
  if (type == "ar"){
    simgen->setMaxEnergy(6000);
    simgen->setMinEnergy(25);

  }
  else if (type == "nr"){
    simgen->setMaxEnergy(5000);
    simgen->setMinEnergy(25);
  }
  else if (type == "af"){
    simgen->setThetaMax(5);//Collimation: in degrees
    simgen->setRecoilEnergy(5350);
    simgen->setProjectileEnergy(5350);
    simgen->setSourcePosition(15,15,50);
    simgen->setProjectileVector(1,1,0);//y-direction
  }
  else if (type == "cf"){
    simgen->setSourcePosition(-1000,75,125);
    simgen->setMaxEnergy(10000);
    simgen->setMinEnergy(0);
  }
  else if (type == "nf"){
    simgen->setSourcePosition(-1000,75,125);
    simgen->setProjectileEnergy(800);
    simgen->setMaxEnergy(10000);
    simgen->setMinEnergy(0);
  }

  bool isGenerated =  simgen->generateRecoil(type,type2);
  //bool isGenerated = 1;
  if (isGenerated){
    TVector3 projVect = simgen->projectileVector().Unit();
    TVector3 recVect = simgen->recoilVector().Unit();
    TVector3 posVect = simgen->recoilPosition();
    //Fix position of random energy alpha events
    double projMass = simgen->getProjectileMass();
    double recoilMass = simgen->getRecoilMass();
    double time = simgen->getTime();
    double projEn = simgen->getProjectileEnergy();
    double recoilEn = simgen->getRecoilEnergy();

    //Run event
    mc->setProjectile(projEn*projVect.X(),projEn*projVect.Y(),projEn*projVect.Z(),projMass);
    mc->setRecoil(recoilEn*recVect.X(),recoilEn*recVect.Y(),recoilEn*recVect.Z(),recoilMass);
    mc->setRecoilCoord(posVect.X(),posVect.Y(),posVect.Z());
    mc->event();

    //Save data in simdata tree
    //time
    simtime = new TDatime(*makeTDatime(time));
    //Mass of recoil particle in keV
    simmass = recoilMass;
    //Recoil energy in keV
    simEnergy = recoilEn;
    //Energy of track in counts (Escint - spacers)
    simIntegral = mc->getAnodeImage()->Integral();
    //Scintillated Energy (Erecoil * Quenching)
    simEscint = mc->getNumPhotons();
    //x position (mm from lower left corner, horizontal on TH2F)
    simx = posVect.X();
    //y position (vertical axis when TH2F is drawn)
    simy = posVect.Y();
    //z position/drift length: perpendicular to imaging plane
    simz = posVect.Z();
    //rough estimate of projected track length
    simlength = getTrackLength(mc->getTrackImage(),recVect.Phi());
    //recoil azimuthal angle: tan(phi) = vy/vx
    simphi = recVect.Phi();
    //recoil angle from z: cos(theta) = vz/|v|
    simtheta = recVect.Theta();
    //projectile direction, kinetic energy and mass
    simprojphi = projVect.Phi();
    simprojtheta = projVect.Theta();
    simprojE = projEn;
    simprojmass = projMass;
    //Fill
    sim->Fill();
  }
 
  return isGenerated;
}//End Event()

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
////Save CCD image and fill dataset                         //
//////////////////////////////////////////////////////////////

void saveEvent(){


  mc->getCCDImage()->SetName("ccd_0");
  TH2F* imagehist = (TH2F*)mc->getCCDImage()->Clone();
  dataset->event()->timeStamp()->Set(simtime->GetYear(),simtime->GetMonth(),simtime->GetDay(),simtime->GetHour(),simtime->GetMinute(),simtime->GetSecond());
  dataset->event()->UTCtimeStamp()->Set(simtime->GetYear(),simtime->GetMonth(),simtime->GetDay(),simtime->GetHour(),simtime->GetMinute(),simtime->GetSecond(),0,1,0);
  //Events start from 1.
  dataset->event()->setRunNumber(1);
  dataset->event()->increaseEventNumber();
  //syntax to save TH2F in a TClonesArray
  new( (*dataset->event()->ccdData())[0]) TH2F(*imagehist);
  dataset->fill();
  delete imagehist;
}//End saveEvent()

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
////Shell to run all simulation functions                   //
////Initialize, generate and save events, write file        //
//////////////////////////////////////////////////////////////

void RunAll(int NEvents,TString opt = "w",TString opt2 = "F",TString opt3 = "CF4",
	    TString filename ="example.root",double gain = 2.5, double noise = 7,
	    double wimpMass = 100e6,
	    int npix = 1024,int nbins = 256, float width = 0.143*1024, 
	    float pressure = 75,float spacerdiam = 1., 
	    float spacerspacing = 20.0,float diffconst = 0.2704,
	    float diffdz = 0.0076){
  

  cout << "Initializing..."<<endl;
  dataSetup(gain, noise,filename,npix,nbins,width,pressure,
       spacerdiam,spacerspacing,diffconst,diffdz);
  cout << "Loading tables, Setting up parameters"<<endl;
  Init(gain, noise,opt2,opt3,npix,nbins,width,pressure,
        spacerdiam,spacerspacing,diffconst,diffdz);

  gROOT->cd();
  cout << "Running Events" << endl;

  if (opt.Contains("em") || opt.Contains("bl")){
    for (int i = 0; i <NEvents; i++){
      generateEmptyEvent();
      saveEvent();
    }
  }
  else  {
  for (int i = 0; i < NEvents; i++){
    simeventnum = i+1;
    //Run event. If successful, save.
    if (i % 20 == 0){
      generateEmptyEvent();
      saveEvent();
    }
    else if(Event(opt,opt2,0,width,0,250)) saveEvent();
    else cout << "WARNING: Failed to generate event!!"<<endl;
    if (i % 10 == 0) cout << "Event " << i+1 << " of "<<NEvents<<endl;
  }
}
  cout<<"Done.  Saving."<<endl;

  dataset->write();
  dataset->file()->Close();
  cout <<"Finished." << endl;

}//End RunAll()

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////makeSimulations()///////////////////////////////////////////////////////
////Use this to run events from the terminal                              //
////Outputs .root file and .txt file                                      //
////Parameters same as in Init().                                         //
////NEvents   :   Number of events each for alpha and neutron             //
////Date      :   Date in YYMMDD format                                   //
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

void makeSimulations(int NEvents = 100, char* date = "000000",char* particle = "w",
		     char* recoil = "F",char* gas = "CF4",
		     float gain = 3.47, float noise = 7,float wimpMass = 100e6,
		     int npix = 1024,int nbins = 256,float width = 0.143*1024,
		     float pressure = 75,float spacerdiam = 1.,
		     float spacerspacing = 20.0,float diffconst = 0.2704,
		     float diffdz = 0.0076){

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////To get the default diffusion values:                                  //
  //// Use diffMeasurement.cxx on dmtpc_run00009.root and                   //
  ////dmtpc_run00010.root.  These were alpha measurements at 3 kV, 75 torr  //
  ////then let diffconst = sigma0^2, and diffdz = diffusion^2               //
  ////////////////////////////////////////////////////////////////////////////

  TString particle0 = particle;
  particle0.ToLower();

    bool done = false;
    TString file = "dmtpc_run";
    file += date;
    file += "_";
    TString textfile;
    int count = 0;
    while(!done){
      count++;
      textfile = file;
      textfile += count;
      textfile += "_sim.txt";

      ifstream ifstr = ifstream(textfile);
      if (!ifstr.is_open()) done = true;

    }
    TString rootfile = file;
    rootfile += count;
    rootfile += "_sim.root";
    makeTextFile(NEvents,particle0,recoil,gas,textfile,gain,noise,wimpMass,npix,nbins,width,
	         pressure,spacerdiam,spacerspacing,diffconst,diffdz);
    RunAll(NEvents,particle0,recoil,gas,rootfile,gain,noise,wimpMass,npix,nbins,width,pressure,
	   spacerdiam,spacerspacing,diffconst,diffdz);

  cout <<"Closing ROOT." << endl;

}//End makeSimulations()
