#include "SimScope.hh"
#include "ScopeDataInfo.hh"
#include "RunGenerator.hh"

#include "../../../MaxCamChannel.hh"
#include "../../../DmtpcGainMap.hh"
#include "../../../MaxCamConfig.hh"
#include "../../../DmtpcStringTools.hh"

#include "TDatime.h"
#include "TROOT.h"

//new
//#include "../../../ScopeHandler.hh"


#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <string.h>
using std::istringstream;
using std::endl;
using std::cin;
using std::cout;
using std::flush;
using std::cerr;
using std::vector;

using namespace std;

ClassImp(RunGenerator);

RunGenerator::RunGenerator()
{
  fNEvents = 1;
  fDatapath = "output/data/";
  fTextpath = "output/summary/";
  fFilename = "example.root";
  fTextfile = "example.txt";
  fFiletag = "dmtpc_mc_";
  fRnd = new TRandom3(0);
  fTime = new TTimeStamp();
  fSimTime = new TDatime();
  fRunNumber = 12345;
  fEventGen = new EventGenerator();
  fSim = NULL;
  fData = NULL;
  fSimTrueClusterArray=0;
  fSaveTrueClusters=false;
  fVerbose = 0;
}

RunGenerator::~RunGenerator()
{

  delete fEventGen;
  delete fRnd;
  delete fTime;
  delete fSimTime;
//  if (fData) delete fData;
  if (fSim) delete fSim;
}

void
RunGenerator::dataSetup(){

  cout <<"Creating dataset"<<endl;
  fData = new DmtpcDataset();
  fData->createRootFile((const Char_t*) fFilename,"recreate");
  fData->setLocation("26-456");
  fData->setKeyword("Simulations");
  fData->setComment("Put any comments here");
  fData->comment()->Write();
  fData->keyword()->Write();
  fData->location()->Write();

  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  ////Set up experimentConfig values. (copied from DmtpcRun.cc::L44-74)        //
  ///////////////////////////////////////////////////////////////////////////////
  cout <<"Saving experimentConfig channels"<<endl;
  Int_t nPars = 0;
  MaxCamChannel* driftHV   = new MaxCamChannel("driftHV",    "Drift Voltage", 0,  0);
  driftHV->currentValue = fEventGen->chamber()->getDriftVoltage()/1000;
  driftHV->setValue = fEventGen->chamber()->getDriftVoltage()/1000;
  new( (*fData->event()->experimentConfig())[nPars++] ) MaxCamChannel(*driftHV);
  delete driftHV;
  
  MaxCamChannel* anodeHV   = new MaxCamChannel("anodeHV",    "Anode Voltage", 1,  1);
  anodeHV->currentValue = fEventGen->chamber()->getAnodeVoltage()/1000;
  anodeHV->setValue = fEventGen->chamber()->getAnodeVoltage()/1000;
  new( (*fData->event()->experimentConfig())[nPars++] ) MaxCamChannel(*anodeHV);
  delete anodeHV;
  
  MaxCamChannel* driftI    = new MaxCamChannel("driftI",     "Drift current", 2, -1);
  driftI->currentValue = 0;
  driftI->setValue = 0;
  new( (*fData->event()->experimentConfig())[nPars++] ) MaxCamChannel(*driftI);
  delete driftI;

  MaxCamChannel* temp0   = new MaxCamChannel("temp0",    "Temperature point 0", 4,  -1);
  temp0->currentValue = fEventGen->chamber()->getTemperature()-273.15;
  temp0->setValue = fEventGen->chamber()->getTemperature()-273.15;
  new( (*fData->event()->experimentConfig())[nPars++] ) MaxCamChannel(*temp0);
  delete temp0;
  
  MaxCamChannel* anodeI    = new MaxCamChannel("anodeI",     "Anode current", 5, -1);
  anodeI->currentValue = 0;
  anodeI->setValue = 0;
  new( (*fData->event()->experimentConfig())[nPars++] ) MaxCamChannel(*anodeI);
  delete anodeI;

  MaxCamChannel* pressure0 = new MaxCamChannel("pressure",   "Gas pressure",  -1, -1);
  pressure0->currentValue = fEventGen->chamber()->getPressure();
  pressure0->setValue = fEventGen->chamber()->getPressure();
  new( (*fData->event()->experimentConfig())[nPars++] ) MaxCamChannel(*pressure0);
  delete pressure0;

  /////////////////////////////////////////////
  // Set up ScopeConfig for each PMT channel //
  /////////////////////////////////////////////
  //Int_t npmt = fEventGen->pmt()->GetEntries();
  cout << "Setting up ScopeConfig" << endl;
  cout << "not yet implemented" << endl;
  Int_t ncam = fEventGen->camera()->GetEntries();  
  /* this has been moved to save event because of the stupid way our code stores
     this information
  //////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////
  ////Set up MaxCamConfig file for the dataset.                                   //
  //////////////////////////////////////////////////////////////////////////////////
  cout <<"Setting MaxCamConfig"<<endl;
  for (Int_t i = 0; i < ncam; i++){
    MaxCamConfig *config0 = new MaxCamConfig("ccdConfig","CCD Camera configuration");
    config0->serialNumber = fEventGen->camera(i)->getSerialNumber();
    config0->cameraID = fEventGen->camera(i)->getCameraNumber();
    config0->row_width = fEventGen->camera(i)->getPixelsX();
    config0->img_rows = fEventGen->camera(i)->getPixelsY();
    config0->hbin = fEventGen->camera(i)->getPixelPerBin();
    config0->vbin = fEventGen->camera(i)->getPixelPerBin();
    config0->ul_x = 0;
    config0->ul_y = 0;
    config0->lr_x = fEventGen->camera(i)->getPixelsX();
    config0->lr_y = fEventGen->camera(i)->getPixelsY();
    config0->CCDTemp = fEventGen->chamber()->getTemperature() - 273.15;
    config0->CCDTempSet = fEventGen->chamber()->getTemperature()-273.15;
    config0->exposureTime = 1000;//in ms
    config0->frameType = 0;
    config0->nFlushes = -1;
    config0->bitDepth = 65535;
    config0->daqTime = -1;
    new( (*fData->event()->ccdConfig())[i]) MaxCamConfig(*config0);
    delete config0;
  }
  */
  fSimIntegral.resize(ncam);
  fSimEscint.resize(ncam);
  fTotalGain.resize(ncam);
  Float_t SimInt[ncam];
  Float_t SimEsc[ncam];
  Float_t TotGain[ncam];
  //////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////
  ////Create a tree to hold the simulated data for comparison w/ reconstruction   //
  //////////////////////////////////////////////////////////////////////////////////
  cout <<"Creating simulated data tree"<<endl;
  fSim = new TTree("Simulation","Simulation");
  fSim->Branch("ncamera",&ncam,"ncamera/I");
  fSim->Branch("mass",&fSimMass,"mass/F");
  fSim->Branch("Sequence",&fSeq,"sequence/I");
  fSim->Branch("E",&fSimEnergy,"E/F");
  fSim->Branch("Escint",&SimEsc,"Escint[ncamera]/F");
  fSim->Branch("phi",&fSimPhi,"phi/F");
  fSim->Branch("theta",&fSimTheta,"theta/F");
  fSim->Branch("x",&fSimX,"x/F");
  fSim->Branch("y",&fSimY,"y/F");
  fSim->Branch("z",&fSimZ,"z/F");
  fSim->Branch("length",&fSimLength,"length/F");
  fSim->Branch("zlength",&fSimZLength,"zlength/F");
  fSim->Branch("deltaZ",&fSimDeltaZ,"deltaZ/F");
  fSim->Branch("time","TDatime",&fSimTime);
  fSim->Branch("eventnum",&fSimEventNum,"eventnum/I");
  fSim->Branch("projE",&fSimProjE,"projE/F");
  fSim->Branch("Integral",&SimInt,"Integral[ncamera]/F");
  fSim->Branch("projPhi",&fSimProjPhi,"projPhi/F");
  fSim->Branch("projTheta",&fSimProjTheta,"projTheta/F");
  fSim->Branch("projMass",&fSimProjMass,"projMass/F");
  fSim->Branch("totalGain",TotGain,"totalGain[ncamera]/F");
  fSim->Branch("trackPhi",&fTrackPhi,"trackPhi/F");
  fSim->Branch("trackTheta",&fTrackTheta,"trackTheta/F");
  fSim->Branch("cosCygnus",&fCosCygnus,"cosCygnus/F");
  cout << "fSaveTrueClusters=" << fSaveTrueClusters << endl;
  if(fSaveTrueClusters){
    fSimTrueClusterArray=new TObjArray();
    fSim->Branch("trueClusterArray","TObjArray", &fSimTrueClusterArray, 32000, 0);
  }

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////Create a tree holding all of the camera parameters
  ///////////////////////////////////////////////////////////////////////
  cout <<"Creating camera tree and bias frames"<<endl;

  Int_t camnum;
  Int_t npix[2], nbin[2];
  Float_t pos[2], width[2];
  Float_t emgain, noisefac,dark, bias0, read, totgain, gain_landau_1, gain_landau_2;
  TString ser;
  TTree* _camtree = new TTree("camera","camera");
  _camtree->Branch("cameraNumber",&camnum,"cameraNumber/I");
  _camtree->Branch("pixels",&npix,"NpixX/I:NpixY");
  _camtree->Branch("bins",&nbin,"NbinX/I:NbinY");
  _camtree->Branch("width",&width,"widthX/F:widthY");
  _camtree->Branch("cameraPosition",&pos,"camX/F:camY");
  _camtree->Branch("noiseFactor",&noisefac,"noiseFactor/F");
  _camtree->Branch("darkCurrent",&dark,"darkCurrent/F");
  _camtree->Branch("bias",&bias0,"bias/F");
  _camtree->Branch("readNoise",&read,"readNoise/F");
  _camtree->Branch("gain",&totgain,"gain/F");
  _camtree->Branch("gain_landau_1",&gain_landau_1,"gain_landau_1/F");
  _camtree->Branch("gain_landau_2",&gain_landau_2,"gain_landau_2/F");
  _camtree->Branch("emGain",&emgain,"emGain/F");
  
  for (Int_t j = 0; j < ncam; j++){
    camnum = fEventGen->camera(j)->getCameraNumber();
    npix[0] = fEventGen->camera(j)->getPixelsX();
    npix[1] = fEventGen->camera(j)->getPixelsY();
    nbin[0] = fEventGen->camera(j)->getBinsX();
    nbin[1] = fEventGen->camera(j)->getBinsY();
    pos[0] = fEventGen->camera(j)->getPositionX();
    pos[1] = fEventGen->camera(j)->getPositionY();
    width[0] = fEventGen->camera(j)->getWidthX();
    width[1] = fEventGen->camera(j)->getWidthY();
    emgain = fEventGen->camera(j)->getEMGain();
    noisefac = fEventGen->camera(j)->getNoiseFactor();
    dark = fEventGen->camera(j)->getDarkCurrent();
    bias0 = fEventGen->camera(j)->getBias();
    read = fEventGen->camera(j)->getReadNoise();
    totgain = fEventGen->camera(j)->getGain();
    gain_landau_1 = fEventGen->camera(j)->getGainLandau_1();
    gain_landau_2 = fEventGen->camera(j)->getGainLandau_2();
    _camtree->Fill();
  /////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////
  //////Create a bias frame w/ correct noise level
  /////////////////////////////////////////////////////////////////////
    cout <<"Creating bias frame"<<endl;
    TString name = "biasFrame";
    Int_t num = fEventGen->camera(j)->getCameraNumber()+1;
    name += num;
    TH2F *biasfr = (TH2F*) fEventGen->camera(j)->getCCDImage()->Clone(name);

    //for(int i=1;i<=256;i++)
    //for(int j=1;j<=256;j++)
    //cout<<biasfr->GetBinContent(i,j)<<endl; 

    biasfr->Reset();
    Int_t xbin = biasfr->GetNbinsX();
    Int_t ybin = biasfr->GetNbinsY();
    //Double_t M = fEventGen->camera(j)->getEMGain();
    //Double_t F = fEventGen->camera(j)->getNoiseFactor();
    //Double_t Nd = fEventGen->camera(j)->getDarkCurrent();
    //Double_t Nr = fEventGen->camera(j)->getReadNoise();
    //Double_t sigma = sqrt(Nd*M*M*F*F + Nr*Nr);
    Double_t bias = fEventGen->camera(j)->getBias();
    for (Int_t i = 1; i <= xbin; i++){
      for (Int_t k = 1; k <= ybin; k++){
        biasfr->SetBinContent(i,k,bias);
	//std::cout<<bias<<std::endl;
      }
    }
    biasfr->Write();
    biasfr->Delete();
  }
  _camtree->Write();
  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  ////Create a tree holding all of the chamber parameters
  ////////////////////////////////////////////////
  cout <<"Creating chamber parameter tree"<<endl;
  Float_t height = fEventGen->chamber()->getHeight();
  Float_t drift = fEventGen->chamber()->getDriftLength();
  Float_t pres = fEventGen->chamber()->getPressure();
  Float_t elperkeV = fEventGen->chamber()->getElectronPerkeV();
  Float_t driftV = fEventGen->chamber()->getDriftVoltage();
  Float_t anodeV = fEventGen->chamber()->getAnodeVoltage();
  Float_t temp = fEventGen->chamber()->getTemperature();
  Float_t domu = fEventGen->chamber()->getDOverMu();
  Float_t diffc = fEventGen->chamber()->getDiffusionConstantTerm();
  Float_t diffdz = fEventGen->chamber()->getDiffusionDzTerm();
  Float_t elifetime = fEventGen->chamber()->getElectronLifetime();
  Float_t elecScint = fEventGen->chamber()->getElectricScintillation();
  Float_t nuclScint = fEventGen->chamber()->getNuclearScintillation();
  Float_t driftVel = fEventGen->chamber()->getDriftVelocity();
  Float_t spacw = fEventGen->chamber()->getSpacerWidth();
  Float_t atten = fEventGen->chamber()->getAttenuation();
  TObjString* spacAxis = new TObjString(fEventGen->chamber()->getSpacerAxis());
  Float_t spacersp = fSpacerSpacing;
  
  TTree* _chamtree = new TTree("chamber","chamber");
  _chamtree->Branch("height",&height,"height/F");
  _chamtree->Branch("driftLength",&drift,"driftLength/F");
  _chamtree->Branch("driftVelocity",&driftVel,"driftVelocity/F");
  _chamtree->Branch("attenuation",&atten,"attenuation/F");
  _chamtree->Branch("electronPerkeV",&elperkeV,"electronPerkeV/F");
  _chamtree->Branch("driftVoltage",&driftV,"driftVoltage/F");
  _chamtree->Branch("anodeVoltage",&anodeV,"anodeVoltage/F");
  _chamtree->Branch("temperature",&temp,"temperature/F");
  _chamtree->Branch("pressure",&pres,"pressure/F");
  _chamtree->Branch("DoverMu",&domu,"DoverMu/F");
  _chamtree->Branch("diffConst",&diffc,"diffConst/F");
  _chamtree->Branch("diffDz",&diffdz,"diffDz/F");
  _chamtree->Branch("elifetime",&elifetime,"elifetime/F");
  _chamtree->Branch("elecScint",&elecScint,"elecScint/F");
  _chamtree->Branch("nuclScint",&nuclScint,"nuclScint/F");
  _chamtree->Branch("spacerWidth",&spacw,"spacerWidth/F");
  _chamtree->Branch("spacerAxis",&spacAxis);
  _chamtree->Branch("spacerSpacing",&spacersp,"spacerSpacing/F");
  _chamtree->Fill();
  _chamtree->Write();
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  ////////Create a file containing other run info.
  ////////////////////////////////////////////////////////
  cout <<"Creating run parameter tree"<<endl;
  TObjString* recoilPart = new TObjString( fEventGen->particleGen()->getRecoilParticle());
  TObjString* projPart = new TObjString( fEventGen->particleGen()->getProjectileParticle());
  TObjString* rectype = new TObjString( fEventGen->particleGen()->getRecoilType());
  TObjString* enopt = new TObjString( fEventGen->particleGen()->getEnergyOption());
  TObjString* posopt = new TObjString(fEventGen->particleGen()->getPositionOption());
  TObjString* timopt = new TObjString( fEventGen->particleGen()->getTimeOption());
  TObjString* diropt = new TObjString( fEventGen->particleGen()->getDirectionOption());
  TObjString* specopt = new TObjString( fEventGen->particleGen()->getSpecialOption());
  TObjString* endfspec = new TObjString( fEventGen->particleGen()->getSpectrum());
  TObjString* totalcs = new TObjString( fEventGen->particleGen()->getCrossSection());
  TObjString* elastscat = new TObjString( fEventGen->particleGen()->getScattering());
  Int_t nevents = fNEvents;
  Float_t thetamax = fEventGen->particleGen()->getThetaMax();
  Float_t timres = fEventGen->getTimeResolution();
  TObjString* srim = new TObjString( fEventGen->getSrimName());

  TTree* _runtree = new TTree("runInfo","runInfo");
  _runtree->Branch("recoilType",&rectype);
  _runtree->Branch("recoilPart",&recoilPart);
  _runtree->Branch("projPart",&projPart);
  _runtree->Branch("srim",&srim);
  _runtree->Branch("endfSpectrum",&endfspec);
  _runtree->Branch("endfElastScat",&elastscat);
  _runtree->Branch("endfTotalCS",&totalcs);
  _runtree->Branch("NEvents",&nevents,"NEvents/I");
  _runtree->Branch("thetaMax",&thetamax,"thetaMax/F");
   _runtree->Branch("energyOption",&enopt);
  _runtree->Branch("positionOption",&posopt);
  _runtree->Branch("directionOption",&diropt);
  _runtree->Branch("timeOption",&timopt);
  _runtree->Branch("specialOption",&specopt);
  _runtree->Branch("timeResolution",&timres,"timeResolution/F");
  _runtree->Fill();
  _runtree->Write();
  gROOT->cd();
  cout <<"All data trees created"<<endl;

  delete _camtree;
  delete _chamtree;
  delete _runtree;
}

void
RunGenerator::emptyEvent()
{

  Int_t ncam = fEventGen->camera()->GetEntries();
  for (Int_t i = 0; i < ncam; i++) {

    fEventGen->camera(i)->emptyImage();

    fSimEscint[i] = 0;
    fSimIntegral[i] = 0;
    fTotalGain[i] = 0;
  }
  //fEventGen->getElectronicSignal()->Reset();
  fEventGen->particleGen()->generateTime();
  fSimTime = fEventGen->particleGen()->getTDatime();
  fSimMass = 0;
  fSeq=0;
  fSimEnergy = 0;
  fSimPhi = 0;
  fSimTheta = 0;
  fSimX = 0;
  fSimY = 0;
  fSimZ = 0;
  fSimLength = 0;
  fSimZLength = 0;
  fSimProjE = 0;
  fSimProjPhi = 0;
  fSimProjTheta = 0;
  fSimProjMass = 0;
  fTrackPhi = 0;
  fTrackTheta = 0;
  fCosCygnus = 0;
}

void
RunGenerator::recoilEvent()
{
  if (fVerbose > 2) cout << "recoilEvent() begin " << endl;
//Run Event with EventGenerator
  fEventGen->runEvent(true);
  if (fVerbose > 2) cout << "recoilEvent() begin2 " << endl;
//Load tree parameters into the proper places
  fSimMass = fEventGen->particleGen()->getRecoilMass();
  
  //correct for index corresponding to actual sequence number
  if(fEventGen->particleGen()->getRunType()=="BigDCTPC_near")
  fSeq = 12+fEventGen->particleGen()->getSeqNum(); 
  if(fEventGen->particleGen()->getRunType()=="BigDCTPC_far")
  fSeq = 24+fEventGen->particleGen()->getSeqNum();
  if(fEventGen->particleGen()->getRunType()=="LittleDCTPC_far")
  fSeq = 1+fEventGen->particleGen()->getSeqNum();
  
  fSimEnergy = fEventGen->particleGen()->getRecoilEnergy();
  for (Int_t i = 0; i < fEventGen->camera()->GetEntries(); i++){
    fSimEscint[i] = fEventGen->camera(i)->getCountBeforeSpacers();
    fSimIntegral[i] = fEventGen->camera(i)->getFinalCount();
    fTotalGain[i] = fEventGen->camera(i)->getGain()*fEventGen->camera(i)->getEMGain();
  }
  fSimPhi = fEventGen->particleGen()->recoilVector().Phi();
  fSimTheta = fEventGen->particleGen()->recoilVector().Theta();
  fSimX = fEventGen->particleGen()->recoilPosition().X();
  fSimY = fEventGen->particleGen()->recoilPosition().Y();
  fSimZ = fEventGen->particleGen()->recoilPosition().Z();
  fSimProjE = fEventGen->particleGen()->getProjectileEnergy();

  fSimProjPhi = fEventGen->particleGen()->projectileVector().Phi();
  fSimProjTheta = fEventGen->particleGen()->projectileVector().Theta();
  fSimProjMass = fEventGen->particleGen()->getProjectileMass();
  fSimTime = fEventGen->particleGen()->getTDatime();
  fSimLength = fEventGen->getLength();
  fSimZLength = fEventGen->getZLength();
  fTrackPhi = fEventGen->getTrackPhi();
  //  cout << "fEventGen->getTrueClusterArray()=" << fEventGen->getTrueClusterArray() << endl;
  //  cout << "fEventGen->getTrueClusterArray()->GetEntries()=" << fEventGen->getTrueClusterArray()->GetEntries() << endl;
  fSimTrueClusterArray = fEventGen->getTrueClusterArray();
  fTrackTheta = fEventGen->getTrackTheta();
  fCosCygnus = fEventGen->particleGen()->getCosCygnus();
}



void
RunGenerator::findFileNames()
{

  Int_t runNum = 1;
  fFilename = fDatapath+fFiletag+"0000";
  fFilename += runNum;
  fFilename += ".root";
  while (!gSystem->AccessPathName(fFilename)){
    runNum++;
    fFilename = fDatapath+fFiletag;
    if (runNum < 10) fFilename+="0000";
    else if (runNum<100) fFilename+="000";
    else if (runNum<1000) fFilename+="00";
    else if (runNum<10000) fFilename+="0";
    fFilename+= runNum;
    fFilename+= ".root";
  }
  cout <<"Unused filename found: "<<fFilename<<endl;
  fTextfile = fFilename;
  fTextfile.ReplaceAll(".root",".sum");
  fTextfile.ReplaceAll(fDatapath,fTextpath);
  fRunNumber = runNum;
  
/*Old date based file name code 
  Int_t date = fTime->GetDate();
  date = date % 20000000;
  TString name;
  if (date < 100000) name = "dmtpc_run0";
  else name = "dmtpc_run";
  date *= 1000;
  date++;
  fFilename = fDatapath;
  fFilename += name;
  fFilename += date;
  fFilename += ".root";
  while (!gSystem->AccessPathName(fFilename)){
    date++;
    fFilename = fDatapath;
    fFilename += name;
    fFilename += date;
    fFilename += ".root";
  }
  fTextfile = fFilename;
  cout <<"Unused filename found: "<<fTextfile<<endl;
  fTextfile.ReplaceAll(".root",".sum");
  fTextfile.ReplaceAll(fDatapath,fTextpath);
  fRunNumber = date;
*/

}

void 
RunGenerator::saveEvent()
{
  //cout << "saveEvent() " << endl;
  fData->event()->timeStamp()->Set(fSimTime->GetYear(),fSimTime->GetMonth(),
				   fSimTime->GetDay(),fSimTime->GetHour(),
				   fSimTime->GetMinute(),fSimTime->GetSecond());

  fData->event()->UTCtimeStamp()->Set(fSimTime->GetYear(),fSimTime->GetMonth(),
				     fSimTime->GetDay(),fSimTime->GetHour(),
				     fSimTime->GetMinute(),fSimTime->GetSecond(),0,1,0);

  fData->event()->setRunNumber(fRunNumber);
  fData->event()->increaseEventNumber();

  ////////////////////////////////////////////////////////////////////////
  //  Write the PMT data to file
  /* FIXME: hardcoded number of triggers 
     -- can there be more than 1 trigger per event? */
  Int_t ntrig   = 1;
  for (Int_t ipmt=0; ipmt<fEventGen->pmt()->GetEntries(); ipmt++) {
    for (Int_t itrig=0; itrig<ntrig; itrig++) {
      //cout << "found a pmt:  " << endl;
      new(  (*fData->event()->rawScopeData())[itrig] ) ScopeWaveformData(*fEventGen->pmt(ipmt)->wf(itrig));
    } // loop over triggers
  } // loop over PMTs
  
  
  //cout<<"NENTRIES SCOPE&&&&&&&&&&&&&&&&&&&&&: "<<fEventGen->scope()->GetEntries()<<endl;
  
  //Store the waveform data TJ
  new(  (*fData->event()->rawScopeData())[0] ) ScopeWaveformData(*fEventGen->scope(0)->wf(0));

// new( (*fData->event()->scopeDataInfo())[0] ) ScopeDataInfo(dataInfo(0));

cout<<"size of data info: "<<dataInfo(0)<<endl;

//addScopeDataInfo(0);
getScopeDataInfo()->setNTriggers(1);
cout<<"Get N triggers: "<<getScopeDataInfo()->getNTriggers()<<endl;

// new( (*fData->event()->scopeDataInfo())[0] ) ScopeDataInfo(*dataInfo(0));

new( (*fData->event()->scopeDataInfo())[0] ) ScopeDataInfo(*getScopeDataInfo());


////////////////////////////////////////////////////////////////////
  // Store the CCD data
  for (Int_t i = 0; i < fEventGen->camera()->GetEntries(); i++){
    TH2S* image = fEventGen->camera(i)->getRawCCDImage(); 

    new( (*fData->event()->ccdData())[i]) TH2S(*image);
    delete image;
  }

  //////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////
  ////Set up MaxCamConfig file for the dataset.                                   //
  //////////////////////////////////////////////////////////////////////////////////
  Int_t ncam = fEventGen->camera()->GetEntries();  
  //  cout <<"Setting MaxCamConfig"<<endl;
  for (Int_t i = 0; i < ncam; i++){
    MaxCamConfig *config0 = new MaxCamConfig("ccdConfig","CCD Camera configuration");
    config0->serialNumber = fEventGen->camera(i)->getSerialNumber();
    config0->cameraID = fEventGen->camera(i)->getCameraNumber();
    config0->row_width = fEventGen->camera(i)->getPixelsX();
    config0->img_rows = fEventGen->camera(i)->getPixelsY();
    config0->hbin = fEventGen->camera(i)->getPixelPerBin();
    config0->vbin = fEventGen->camera(i)->getPixelPerBin();
    config0->ul_x = 0;
    config0->ul_y = 0;
    config0->lr_x = fEventGen->camera(i)->getPixelsX();
    config0->lr_y = fEventGen->camera(i)->getPixelsY();
    config0->CCDTemp = fEventGen->chamber()->getTemperature() - 273.15;
    config0->CCDTempSet = fEventGen->chamber()->getTemperature()-273.15;
    config0->exposureTime = 1000;//in ms
    config0->frameType = 0;
    config0->nFlushes = -1;
    config0->bitDepth = 65535;
    config0->daqTime = -1;
    new( (*fData->event()->ccdConfig())[i]) MaxCamConfig(*config0);
    delete config0;
  }

  fData->fill();

  Float_t simint[ncam];
  Float_t simesc[ncam];
  Float_t totgain[ncam];
  fSim->SetBranchAddress("ncamera",&ncam);
  fSim->SetBranchAddress("Integral",&simint);
  fSim->SetBranchAddress("Escint",&simesc);
  fSim->SetBranchAddress("totalGain",&totgain);
  for (Int_t i = 0; i<ncam;  i++){
    simint[i] = fSimIntegral[i];
    simesc[i] = fSimEscint[i];
    totgain[i] = fTotalGain[i];
  }
  
  //  fSimTrueCluster=new MaxCamClusterImage();
  //  fSimTrueCluster=new EventGenerator();
  //  cout << "fSimTrueCluster=" << fSimTrueCluster << endl;
  //  fSim->GetBranch("trueClusterArray")->SetAddress(&fSimTrueClusterArray);
  //  fSimTrueClusterArray->AddAtAndExpand(new MaxCamClusterImage(new TH2F(),new TDatime()),0);
  //  cout << "fSimTrueClusterArray->GetEntries()=" << fSimTrueClusterArray->GetEntries() << endl;
  fSim->Fill();
  // this pointer gets recreated everytime EventGenerator iterates, so delete it
  // so we don't eat up to much memory...
  if(fSimTrueClusterArray) fSimTrueClusterArray->Delete();
   
  fData->event()->isCleared = false;
  fData->clearEventMemory();
}

//Int_t
//RunGenerator::getNlevels(Int_t scope) {
//  return (Int_t)(TMath::Power(2,fEventGen->scope(scope)->getNbitsPerSamp()));
//}
//
//Float_t
//RunGenerator::getVoltageStep(Int_t scope, Int_t chan)
//{
//  Int_t nlevels = getNlevels(scope);
//  Float_t vmin  = fEventGen->scope(scope)->chan(chan)->getVoltageMin();
//  Float_t vmax  = fEventGen->scope(scope)->chan(chan)->getVoltageMax();
//  return (vmax-vmin)/nlevels;
//}

//Float_t
//RunGenerator::levelOfVolts(Int_t scope, Int_t chan, Float_t level_volts) 
//{
//  Float_t vmin  = fEventGen->scope(scope)->chan(chan)->getVoltageMin();
//  int level     = (level_volts-vmin)/getVoltageStep(scope, chan);
//
//  Int_t nlevels = getNlevels(scope);
//  if (level<0) level=0;
//  if (level > nlevels) level=nlevels-1;
//    
//  return level;
//}

void 
RunGenerator::endRun()
{
  //  cout <<"Writing simulation values"<<endl;
  fData->file()->cd();
  fSim->Write();
  delete fSim;
  fSim = NULL;
  cout <<"Writing event data"<<endl;
  fData->write();
  gROOT->cd();
}

void
RunGenerator::readRunParameters(TString filename)
{
  cout <<"Reading run parameters."<<endl;
  ifstream ifstr(filename);
  if (!ifstr.is_open()){
    cerr <<"File: "<<filename<<" not found!"<<endl;
    assert(0);
  }
  string line;
  TString name;
  TString a;
  Double_t m;
  while (!ifstr.eof()){
    getline(ifstr,line);
    if (line[0]=='#' || line=="") continue;
    istringstream linestr(line);
    linestr >> name;
    //cout << name <<endl;
    if(name == "NumberOfEvents"){
      Int_t num;
      linestr >> num;
      fNEvents = num;
    }
    else if (name == "FileTag"){
      linestr>>a;
      fFiletag = a;
    }
    else if (name == "DataFile"){
      linestr >> a;
      fFilename = a;
    }
    else if (name == "TextFile"){
       linestr >> a;
       fTextfile = a;
    }
    else if(name=="DataPath"){
      linestr >> a;
      fDatapath = a;
    }
    else if (name=="TextPath"){
      linestr >> a;
      fTextpath = a;
    }
    else if(name =="EventType"){
      linestr >> a;
      fEventGen->particleGen()->setRecoilType(a);
    }
    else if(name =="ProjectileParticle"){
      linestr >> a;
      fEventGen->particleGen()->setProjectileParticle(a);
    }
    else if(name =="ProjectileMass"){
      linestr >> m;
      fEventGen->particleGen()->setProjectileMass(m);
    }
    else if(name =="RecoilParticle"){
      linestr >> a;
      fEventGen->particleGen()->setRecoilParticle(a);
    }
    else if(name =="RecoilMass"){
      linestr >> m;
      fEventGen->particleGen()->setRecoilMass(m);
    }
    else if(name =="SRIMfile"){
      linestr>>a;
      fEventGen->setSrimName(a);
    }
    else if (name =="RunType"){
      linestr>>a;
      fEventGen->particleGen()->setRunType(a);
    }
    else if (name =="Config"){
      linestr>>a;
      fEventGen->particleGen()->setConfigType(a);
    }
    else if (name =="Location"){
      linestr>>a;
      fEventGen->particleGen()->setLocation(a);
    }
    else if(name =="EnergyOption"){
      linestr>>a;
      fEventGen->particleGen()->setEnergyOption(a);
    }
    else if (name == "DecayChain")
    {
      linestr >> a ; 
      TString key; 
      linestr >> key; 
      fEventGen->particleGen()->setDecayChain(a,key); 
    }
    else if(name =="MinEnergy"){
      linestr>>m;
      fEventGen->particleGen()->setMinEnergy(m);
    }
    else if(name =="MaxEnergy"){
      linestr>>m;
      fEventGen->particleGen()->setMaxEnergy(m);
    }
    else if(name =="FixEnergy"){
      linestr>>m;
      fEventGen->particleGen()->setProjectileEnergy(m);
      // in case user wants a gaussian
      fEventGen->particleGen()->setProjectileGaussianEnergy(m);
    }
    else if(name =="EnergySpread"){
      linestr>>m;
      fEventGen->particleGen()->setProjectileGaussianEnergySpread(m);
    }
   else if(name =="PositionOption"){
      linestr>>a;
      fEventGen->particleGen()->setPositionOption(a);
    }
    else if(name =="MinPositionLimits"){
      Double_t x,y,z;
      linestr>>x>>y>>z;
      fEventGen->particleGen()->setPositionMinValues(x,y,z);
    }
    else if(name =="MaxPositionLimits"){
      Double_t x,y,z;
      linestr>>x>>y>>z;
      fEventGen->particleGen()->setPositionMaxValues(x,y,z);
    }
    else if(name =="RingPositionValues")
    {
      Int_t n; 
      Double_t rin,rout, x0,y0, zmin,zmax,dz, ztop, fractop, zbottom, fracbottom; 
      linestr>>n>>x0>>y0>>rin>>rout>>zmin>>zmax>>dz >> ztop >> fractop >> zbottom >> fracbottom;
      fEventGen->particleGen()->setRingValues(n,x0,y0,rin,rout,zmin,zmax,dz, ztop, fractop, zbottom, fracbottom); 
    }
    else if(name =="FixPosition"){
      Double_t x,y,z;
      linestr>>x>>y>>z;
      fEventGen->particleGen()->setRecoilPosition(x,y,z);
    }
    else if(name =="TimeOption"){
      linestr >> a;
      fEventGen->particleGen()->setTimeOption(a);
    }
    else if(name =="MinTime"){
      Int_t year, month, day, hour, min, sec;
      linestr >> year >> month>>day>>hour>>min>>sec;
      fEventGen->particleGen()->setBeginTime(year,month,day,hour,min,sec);
    }
    else if(name =="MaxTime"){
      Int_t year, month, day, hour, min, sec;
      linestr >> year >> month>>day>>hour>>min>>sec;
      fEventGen->particleGen()->setEndTime(year,month,day,hour,min,sec);
    }
    else if(name =="FixTime"){
      Int_t year, month, day, hour, min, sec;
      linestr >> year >> month>>day>>hour>>min>>sec;
      TTimeStamp stamp(year,month,day,hour,min,sec);
      fEventGen->particleGen()->setTime(&stamp);
    }
    else if(name == "TimeStep"){
      linestr>>m;
      fEventGen->particleGen()->setTimeStep(m);
    }
    else if(name =="DirectionOption"){
      linestr>>a;
      fEventGen->particleGen()->setDirectionOption(a);
    }
    else if(name =="SourcePosition"){
      Double_t x,y,z;
      linestr>>x>>y>>z;
      fEventGen->particleGen()->setSourcePosition(x,y,z);
    }
    else if(name =="SourceDirection"){
      Double_t x,y,z;
      linestr>>x>>y>>z;
      fEventGen->particleGen()->setSourceDir(x,y,z);
    }
    else if(name =="ThetaMax"){
      linestr>>m;
      fEventGen->particleGen()->setThetaMax(m);
    }
    else if(name =="SpecialOption"){
      linestr>>a;
      fEventGen->particleGen()->setSpecialOption(a);
    }
    else if(name =="FissionSpectrumENDF"){
      linestr>>a;
      fEventGen->particleGen()->setSpectrum(a);
    }
    else if(name =="TotalCS_ENDF"){
      linestr>>a;
      fEventGen->particleGen()->setCrossSection(a);
    }
    else if(name =="ElastScatterDCS_ENDF"){
      linestr>>a;
      fEventGen->particleGen()->setScattering(a);
    }
    else if(name =="TimeResolution"){
      linestr >> m;
      fEventGen->setTimeResolution(m);
    }
    else if(name =="WimpEnergyOption"){
       linestr>>a;
       fEventGen->particleGen()->setTheoryEnergyOption(a);
    }
    else if(name =="SaveTrueClusters"){
       linestr>>a;
       setSaveTrueClusters(a);
    }
    else cout << "Keyword " << name <<" in " << filename << " not found."<<endl;
  }
  if (fEventGen->particleGen()->getTimeOption().Contains("ser")){
    fEventGen->particleGen()->setTime(fEventGen->particleGen()->getBeginTime());
  }
}

void 
RunGenerator::readScopeProperties(TString filename)
{
  cout <<"Reading Scope properties." << endl;
  ifstream ifstr(filename);
  cout << filename << endl;
  string line;
  TString name;
  Int_t nscope;
  while(!ifstr.eof()) {
    getline(ifstr, line);
    DmtpcStringTools::ltrim(line);
    if (line[0]=='#' || line=="") continue;

    istringstream linestr(line);
    linestr >> name >> nscope;//Keyword  ScopeNumber  Value(s)

    Int_t iscope=0;
    for (;iscope<fEventGen->scope()->GetEntries(); iscope++) {
      if (nscope==fEventGen->scope(iscope)->getScopeNumber()) break;
    }
    if (iscope==fEventGen->scope()->GetEntries()) {
      cout << "RunGenerator: found new Scope, no=" << nscope << endl;
      // scope number not in list of existing scopes:  add and assign a new number
      fEventGen->scope()->Add(new SimScope());
      fEventGen->scope(nscope)->setScopeNumber(nscope);

    }

    if (name == "BoardType"){
      TString btype;
      linestr >> btype;
      fEventGen->scope(nscope)->setBoardType(btype);
      
      
    } else if (name == "SerialNumber"){
      TString ser;
      linestr >> ser;
      fEventGen->scope(nscope)->setSerialNumber(ser);
    } else if (name == "ClockRate"){
      Float_t cr;
      linestr >> cr;
      fEventGen->scope(nscope)->setClockRate(cr);
    } else if (name == "RecordPreSize"){
      Int_t sz;
      linestr >> sz;
      fEventGen->scope(nscope)->setRecordPreSize(sz);
    } else if (name == "RecordLength"){
      Int_t sz;
      linestr >> sz;
      fEventGen->scope(nscope)->setRecordLength(sz);
    } else if (name == "ChannelName"){
      cout << "found ChannelName" << endl;
      Int_t chan;
      TString cname;
      linestr >> chan >> cname;
      prepScopeChannel(nscope, chan);
      fEventGen->scope(nscope)->chan(chan)->setName(cname);
    } else if (name == "ConnectedDevice"){
      cout << "found ConnectedDevice" << endl;
      Int_t chan;
      TString cdname;
      linestr >> chan >> cdname;
      prepScopeChannel(nscope, chan);  

 cout << "setConnectedDevice: " << cdname << endl;
      fEventGen->scope(nscope)->chan(chan)->setConnectedDevice(cdname);
      cout << "post setConnectedDevice: " << cdname << endl;      
    } else if ((name == "VoltageMin")||(name=="VoltageMax")) {
      cout << "found VoltageMin or VoltageMax: " << name << endl;
      Int_t chan;
      Float_t val;
      linestr >> chan >> val;
      prepScopeChannel(nscope, chan);
      //TJ
      
      //These lines do not have any effect on the output after compiling?

      

      if (name == "VoltageMin") 
	fEventGen->scope(nscope)->chan(chan)->setVoltageMin(val);
      else
	fEventGen->scope(nscope)->chan(chan)->setVoltageMax(val);
    } else{
      if (name!="ScopeNumber") 
        cout <<"Keyword " << name <<" in "<< filename<<" not found."<<endl;
    }
  }
}

void
RunGenerator::prepScopeChannel(Int_t scope, Int_t chan) 
{
  // if this channel doesn't exist, create it
  /* FIXME:  assumes here that the scope already exists... */
  Int_t ichan=0;
  for (;ichan<fEventGen->scope(scope)->chans()->GetEntries(); ichan++) {
    if (chan==fEventGen->scope(scope)->chan(ichan)->getChanNumber()) break;
  }
  if (ichan==fEventGen->scope(scope)->chans()->GetEntries()) {
    cout << "RunGenerator: found new Scope Channel, no=" << chan << endl;
    // scope chan number does not yet exist:  add and assign a new number
    scopeinfo=new ScopeDataInfo();
    fEventGen->scope(scope)->chans()->Add(scopeinfo);
    
    addScopeDataInfo(0);
    fEventGen->scope(scope)->chan(chan)->setNumber(chan);

      
  }
}


void 
RunGenerator::readPMTProperties(TString filename)
{
  cout <<"Reading PMT properties." << endl;
  ifstream ifstr(filename);
  cout << filename << endl;
  string line;
  TString name;
  Int_t npmt;
  while(!ifstr.eof()) {
    getline(ifstr, line);
    DmtpcStringTools::ltrim(line);
    if (line[0]=='#' || line=="") continue;

    istringstream linestr(line);
    linestr >> name >> npmt;//Keyword  CameraNumber  Value(s)

    Int_t ipmt=0;
    for (;ipmt<fEventGen->pmt()->GetEntries(); ipmt++) {
      if (npmt==fEventGen->pmt(ipmt)->getPMTNumber()) break;
    }
    if (ipmt==fEventGen->pmt()->GetEntries()) {
      cout << "RunGenerator: found new PMT, no=" << npmt << endl;
      // pmt number not in list of existing pmts:  add and assign a new number
      fEventGen->pmt()->Add(new SimPMT());
      fEventGen->pmt(npmt)->setPMTNumber(npmt);
    }

    if (name == "PMTPosition"){
      Double_t x,y;
      linestr >> x >> y;
      fEventGen->pmt(npmt)->setPosition(x,y);
    } else if (name == "SerialNumber"){
      TString ser;
      linestr >> ser;
      fEventGen->pmt(npmt)->setSerialNumber(ser);
    } else if (name == "PMTDiameter"){
      Double_t diam;
      linestr >> diam;
      fEventGen->pmt(npmt)->setDiameter(diam);
    } else if (name == "PMTDistance"){
      Double_t dist;
      linestr >> dist;
      fEventGen->pmt(npmt)->setDistance(dist);
    } else if (name == "PMTPosition"){
      Double_t x,y;
      linestr >> x >> y;
      fEventGen->pmt(npmt)->setPosition(x,y);
    } else{
      if (name!="PMTNumber") 
        cout <<"Keyword " << name <<" in "<< filename<<" not found."<<endl;
    }
  }
}

void
RunGenerator::readCameraProperties(TString filename)
{
  cout <<"Reading camera properties from: "<<filename<<endl;
  ifstream ifstr(filename);
  cout << filename << endl;
  string line;
  TString name;
  TString a;
  TString b;
  Double_t m;
  Int_t i;
  Int_t ncam;
  while (!ifstr.eof()){
    getline(ifstr,line);
    if (line[0]=='#' || line=="") continue;
    istringstream linestr(line);
    linestr >> name >> ncam;//Keyword  CameraNumber  Value(s)
    //cout << fEventGen->camera()->GetEntries() << endl;

    int icam=0; 
    for (; icam<fEventGen->camera()->GetEntries(); icam++) {
        if (ncam==fEventGen->camera(icam)->getCameraNumber()) break;
    }
    if (icam==fEventGen->camera()->GetEntries()) {
        cout << "RunGenerator: found new camera, no=" << ncam << endl;
        // camera number not in list of existing cameras: add and assign new number
        fEventGen->camera()->Add(new SimCamera());
        fEventGen->camera(ncam)->setCameraNumber(ncam);
    }
    
    //cout << name << endl;
    if (name == "CameraPosition"){
      Double_t x,y;
      linestr >> x >> y;
      fEventGen->camera(ncam)->setPosition(x,y);
    }
    else if (name == "SerialNumber"){
      TString ser;
      linestr >> ser;
      fEventGen->camera(ncam)->setSerialNumber(ser);
    }
    else if (name == "CameraBins"){
      Int_t x,y;
      linestr >> x >> y;
      fEventGen->camera(ncam)->setBins(x,y);
    }
    else if (name =="ImageWidths"){
      Double_t x,y;
      linestr>>x>>y;
      fEventGen->camera(ncam)->setWidths(x,y);
    }
    else if (name =="PixelsPerBin"){
      linestr >> i;
      fEventGen->camera(ncam)->setPixelPerBin(i);
    }
    else if (name == "ReadNoise"){
      linestr >> m;
      fEventGen->camera(ncam)->setReadNoise(m);
    }
    else if (name == "Bias"){
      linestr >> m;
      fEventGen->camera(ncam)->setBias(m);
    }
    else if (name =="Gain"){        
      linestr >> m;
      fEventGen->camera(ncam)->setGain(m);
    }
    else if (name =="Gain_Landau"){
      double landau2;
      linestr >> m >> landau2;
      fEventGen->camera(ncam)->setGain_Landau(m,landau2);
    }
    else if (name == "RadialParams"){
      double radial1,radial2,radial3;
      linestr >> radial1 >> radial2 >> radial3;
      fEventGen->camera(ncam)->setRadialParams(radial1,radial2,radial3);
    }
    else if (name == "PhotonsADU"){
      double photonsadu;
      linestr >> photonsadu;
      fEventGen->camera(ncam)->setPhotonsADU(photonsadu);
    }
    else if (name == "GainMap")
    {
      linestr >> a >> b;  
      TFile f(a); 
      gROOT->cd(); 
      TH2F * img = (TH2F*) ((DmtpcGainMap*)f.Get(b))->getGainMap();
      cout << img->Integral() << endl; 
      fEventGen->camera(ncam)->setGainMap(img); 
      fEventGen->camera(ncam)->normalizeGainMap(); 
      f.Close(); 
    }
    else if (name == "ActiveRegionUnits")
    {
      TString a;
      linestr >> a;
      fEventGen->camera(ncam)->setActiveRegionUnits(a);
    }
    else if (name == "ActiveRegion")
    {
      Int_t npoints; 
      Double_t ax, ay; 
      linestr >> npoints; 

      TCutG cut("active", npoints); 

      for (int i = 0; i < npoints; i++)
      {
        linestr >> ax >> ay; 
        cut.SetPoint(i, ax,ay); 
      }
      fEventGen->camera(ncam)->setActiveRegion(&cut); 

    }
    else if (name =="NoiseFactor"){
      linestr >>m;
      fEventGen->camera(ncam)->setNoiseFactor(m);
    }
    else if (name =="DarkCurrent"){
      linestr >> m;
      fEventGen->camera(ncam)->setDarkCurrent(m);
    }
    else if (name =="EMGain"){
      linestr >>m;
      fEventGen->camera(ncam)->setEMGain(m);
    }
    else{
      if (name!="CameraNumber") 
        cout <<"Keyword " << name <<" in "<< filename<<" not found."<<endl;
    }

  }
}

void
RunGenerator::readDetectorProperties(TString filename)
{
  cout <<"Reading chamber properties."<<endl;
  ifstream ifstr(filename);
  string line;
  TString name;
  TString a;
  Double_t m;
  Double_t diffdz = 0;
  Double_t elifetime = 0;
  Int_t usediffdz = 0;
  while (!ifstr.eof()){
    getline(ifstr,line);
    istringstream linestr(line);
    if (line[0]=='#' || line =="") continue;
    linestr >> name;
    //cout << name << endl;
    if (name == "SpacerDiameter"){
      linestr >> m;
      fEventGen->chamber()->setSpacerWidth(m);
    }
    else if (name == "SpacerSpacing"){
      linestr >>m;
      fSpacerSpacing = m;
    }
    else if (name == "Temperature"){
      linestr >> m;
      fEventGen->chamber()->setTemperature(m);
    }
    else if (name == "SpacerAxis"){
      linestr >> a;
      fEventGen->chamber()->setSpacerAxis(a);
    } 
    else if (name == "Height"){
      linestr >> m;
      fEventGen->chamber()->setHeight(m);
    }
    else if (name == "DriftLength"){
      linestr >> m;
      fEventGen->chamber()->setDriftLength(m);
    }
    else if (name == "Pressure"){
      linestr >> m;
      fEventGen->chamber()->setPressure(m);
    }
    else if (name == "DriftVoltage"){
      linestr >> m;
      fEventGen->chamber()->setDriftVoltage(m);
    }
    else if (name == "AnodeVoltage"){
      linestr >> m;
      fEventGen->chamber()->setAnodeVoltage(m);
    }
    else if (name == "DiffusionConstantTerm"){
      linestr >> m;
      fEventGen->chamber()->setDiffusionConstantTerm(m);
    }
    else if (name == "DiffusionDzTerm"){
      linestr >> diffdz >> usediffdz;
    }
    else if (name == "ElectronLifetime"){
      linestr >> elifetime;
      fEventGen->chamber()->setElectronLifetime(elifetime);   
    }
    else if (name == "ElectricScintillation"){
      linestr >> m;
      fEventGen->chamber()->setElectricScintillation(m);
    }
    else if (name == "NuclearScintillation"){
      linestr >> m;
      fEventGen->chamber()->setNuclearScintillation(m);
    }
    else if (name == "Attenuation"){
      linestr >> m;
      fEventGen->chamber()->setAttenuation(m);
    }
    else if (name == "ElectronPerkeV"){
      linestr >> m;
      fEventGen->chamber()->setElectronPerkeV(m);
    }
    else if (name == "DriftVelocity"){
      linestr >> m;
      fEventGen->chamber()->setDriftVelocity(m);
    }
    else if (name == "LongDiffusionConstTerm"){
      linestr >> m;
      fEventGen->chamber()->setLongDiffusionConstantTerm(m);
    }
    else if (name == "LongDiffusionDzTerm"){
      linestr >> m;
      fEventGen->chamber()->setLongDiffusionDzTerm(m);
    }
    else if (name == "OrientationAngle"){
      linestr >> m;
      fEventGen->chamber()->setOrientationAngle(m);
      fEventGen->particleGen()->setRotationAngle(m);
    }
    else if (name == "TopOrBottomTPC"){
      Bool_t tb;
      linestr >> tb;
      fEventGen->particleGen()->setTopOrBottomTPC(tb);
    }
    else cout << "Keyword " << name <<" in "<< filename <<" not found."<<endl;
  }
  if (usediffdz != 0) fEventGen->chamber()->setDiffusionDzTerm(diffdz);
  else fEventGen->chamber()->setDOverMu();
}

void
RunGenerator::readParametersFromFile(TString filename){
  cout <<"Reading from "<<filename<<endl;
  ifstream ifstr(filename);
  string name;
  istringstream iss((string)filename);
  while(!iss.eof()){
    getline(iss,name,'/');
  }
  TString path = filename;
  path.ReplaceAll(name,"");
  TString file;
  ifstr >> file;
  readRunParameters(path+file);
  ifstr>>file;
  readDetectorProperties(path+file);
  while(!ifstr.eof()){
    ifstr>>file;
    file = file.Strip(TString::kLeading);
    cout << "----- " << file << endl;
    if ((file.BeginsWith("#")) || (file.Length() <= 0)) continue;
    
    cout<<"FILE!!!! " <<path+file<<endl;
    
    if (file.Contains("PMT")) readPMTProperties(path+file);
    //else if (file.Contains("Scope")) 
    cout<<"here1"<<endl;
    if (file.Contains("scope"))
    readScopeProperties(path+file);
    cout<<"here2"<<endl;
    if (file.Contains("camera"))
    readCameraProperties(path+file);
    //if (file.Length() > 0) readCameraProperties(path+file);
  }
}

void
RunGenerator::outputTextFile(Bool_t fromFile, TString filename){
  ofstream ofstr(fTextfile);
  //cout << "fTextfile = " << fTextfile << endl;
  //cout << "fromFile = " << fromFile << endl;
  //cout << "filename = " << filename << endl;
  ofstr <<"#Run Summary for file: "<<endl << fFilename << endl;
  ofstr <<"#==================================="<<endl<<endl;
  if (fromFile){
    ifstream ifstr(filename);
    string name;
    istringstream iss((string)filename);
    while(!iss.eof()){
      getline(iss,name,'/');
    }
    TString path = filename;
    path.ReplaceAll(name,"");
    cout << "path = " << path << endl;
    TString file;

    // Run Properties
    ifstr >> file;
    cout <<"opening file: "<<path+file<<endl;
    ifstream filestr(path+file);
    string line;
    while (!filestr.eof()){
      getline(filestr,line);
      //cout << line << endl;
      ofstr << line <<endl;
      if (line=="") break;
    }
    ofstr <<endl <<"#Chamber Properties: "<<endl;
    ofstr <<"#==================================="<<endl;
    filestr.clear();
    filestr.close();

    // Chamber Properties
    ifstr >> file;
    filestr.open(path+file);
    while (!filestr.eof()){
      getline(filestr,line);
      ofstr << line <<endl;
      if (line=="") break;
    }

    // Camera, PMT and Scope properties
    while (!ifstr.eof()){
      ifstr >> file;
      file = file.Strip(TString::kLeading);
      if ((file.BeginsWith("#")) || (file=="")) continue;
      if (file.Contains("Scope")) {
	ofstr <<endl<<"#Scope Properties: "<<endl;
	ofstr <<"#==================================="<<endl;
      } else if (file.Contains("PMT")) {
	ofstr <<endl<<"#PMT Properties: "<<endl;
	ofstr <<"#==================================="<<endl;
      } else {  // camera (to maintain backwards compatibility)
	ofstr <<endl<<"#Camera Properties: "<<endl;
	ofstr <<"#==================================="<<endl;
      }
      filestr.clear();
      filestr.close();
      filestr.open(path+file);
      while (!filestr.eof()){
        getline(filestr,line);
        if (line=="")break;
        ofstr << line <<endl;        
      }
    }
    ifstr.clear();
    ifstr.close();
  }
  else{
    cout <<"Outputting text file from MC objects not yet completed!"<<endl;
    ofstr <<"Run Parameters"<<endl;
    ofstr <<"==================================="<<endl;
    ofstr <<"NumberOfEvents:\t\t"<< fNEvents<<endl;
    ofstr <<"EventType:\t\t"<< fEventGen->particleGen()->getRecoilType()<<endl;
    ofstr <<"ProjectileParticle:\t"<<endl;
    ofstr <<"ProjectileMass:\t\t"<<endl;
    ofstr <<"RecoilParticle:\t\t"<<endl;
    ofstr <<"RecoilMass:\t\t"<<endl;
    ofstr <<"SRIMfile:\t\t"<<endl;
    ofstr <<"EnergyOption:\t\t"<<endl;
    ofstr <<"MinEnergy:\t\t"<<endl;
    ofstr <<"MaxEnergy:\t\t"<<endl;
    ofstr <<"FixEnergy:\t\t"<<endl;
    ofstr <<"PositionOption:\t\t"<<endl;
    ofstr <<"MinPositionLimits:\t"<<endl;
    ofstr <<"MaxPositionLimits:\t"<<endl;
    ofstr <<"FixPosition:\t\t"<<endl;
    ofstr <<"TimeOption:\t\t"<<endl;
    ofstr <<"MinTime:\t\t\t"<<endl;
    ofstr <<"MaxTime:\t\t\t"<<endl;
    ofstr <<"FixTime:\t\t\t"<<endl;
    ofstr <<"DirectionOption:\t\t"<<endl;
    ofstr <<"SourcePosition:\t\t"<<endl;
    ofstr <<"SourceDirection:\t\t"<<endl;
    ofstr <<"ThetaMax:\t\t"<<endl;
    ofstr <<"SpecialOption:\t\t"<<endl;
    ofstr <<"FissionSpectrumENDF:\t"<<endl;
    ofstr <<"TotalCS_ENDF:\t\t" <<endl;
    ofstr <<"ElastScatterDCS_ENDF:\t"<<endl;
  }
  ofstr.close();
}

void
RunGenerator::setSpacers()
{

  TString axis = fEventGen->chamber()->getSpacerAxis();
  axis.ToLower();
  Double_t max = 0,min = 0,mid = 0;
  if (axis =="x" || axis =="y"){
    if (axis == "x"){
      Double_t width = fEventGen->camera(0)->getWidthY(), y = fEventGen->camera(0)->getPositionY();
      min = y - width/2.; 
      max = y + width/2.;
      for (Int_t i = 1; i < fEventGen->camera()->GetEntries(); i++){
        width = fEventGen->camera(i)->getWidthY();
        y = fEventGen->camera(i)->getPositionY();
        min = TMath::Min(min,y-width/2.);
        max = TMath::Max(max,y+width/2.);       
      }
      mid = 0.5*(min+max);
    }
    else if (axis == "y"){
      Double_t width = fEventGen->camera(0)->getWidthX(), x = fEventGen->camera(0)->getPositionX();
      min = x - width/2.; 
      max = x + width/2.;
      for (Int_t i = 1; i < fEventGen->camera()->GetEntries(); i++){
        width = fEventGen->camera(i)->getWidthX();
        x = fEventGen->camera(i)->getPositionX();
        min = TMath::Min(min,x-width/2.);
        max = TMath::Max(max,x+width/2.);
      }
      mid = 0.5*(min+max);
    }
    cout <<"Max: "<<max<<" Min: "<<min<<endl;
    Double_t position = 0;
    fEventGen->chamber()->addSpacer(mid+position);
    while (mid+position+fSpacerSpacing <= max){
      position += fSpacerSpacing;
      fEventGen->chamber()->addSpacer(mid+position);
      fEventGen->chamber()->addSpacer(mid-position);
      cout <<"Setting spacer at "<<mid+position<<endl;
    }

  }
}

void
RunGenerator::initialize()
{
  cout << "=============================" << endl;
  cout << "======  Initializing  =======" <<endl;
  cout << "=============================" << endl;
  for (Int_t i = 0; i < fEventGen->camera()->GetEntries(); i++) fEventGen->camera(i)->setImage();
  fEventGen->setSrim();
//  fEventGen->setElectronicSignal();

  // Setup PMT waveforms...
  for (Int_t ipmt=0; ipmt<fEventGen->pmt()->GetEntries(); ipmt++) {
    TString devName = fEventGen->pmt(ipmt)->getSerialNumber();
    // Determine which scopenum and channelnum this pmt is connected to
    // loop over all scopes and channels and look for the PMT device name
    bool matchFound = false;
    Int_t iscope=0;
    Int_t ichan =0;
    for (; iscope<fEventGen->scope()->GetEntries(); iscope++) {
      cout << "iscope = " << iscope << endl;
      for (; ichan<fEventGen->scope(iscope)->chans()->GetEntries(); ichan++) {
	cout << "ichan = " << ichan << endl;
	if (devName==fEventGen->scope(iscope)->chan(ichan)->getConnectedDevice()) {
	  cout << "found ConnectedDevice: " << devName << " --> scope, chan="<<iscope<<", " << ichan << endl;
	  matchFound=true;
	  break;
	}
      }
      if (matchFound) break;  // needed to break out of the nested for loop.
    }
    cout << "got here" << endl;
    TString msg = "No correspondence between PMT and Scope found for PMT name = ";
    msg+= devName;
    if (!matchFound) assert(!msg);

    // determine which scope the PMT is paired to
    // pass the scope, the scope number and the channel number to the PMT
    cout << "assign scope" << endl;
    cout << "iscope, ichan = " << iscope << ", " << ichan << endl;
    fEventGen->pmt(ipmt)->assignScope(fEventGen->scope(iscope), ichan);
    cout << "set wfs" << endl;
    fEventGen->pmt(ipmt)->setWfs();
  }

    fEventGen->scope(0)->assignScope(fEventGen->scope(0), 0);
    //TJ
    fEventGen->scope(0)->setWfs();
    //fEventGen->scope(0)->setWfs(fEventGen->scope(0), 0);


  fEventGen->setPressure(fEventGen->chamber()->getPressure());
  if (fEventGen->particleGen()->getSpecialOption() == "endf"){
    fEventGen->particleGen()->setENDFfiles();
  }
  setSpacers();
  dataSetup();


  cout << "=============================" << endl;
  cout << "====== end of init    =======" <<endl;
  cout << "=============================" << endl;
}

void
RunGenerator::run(Bool_t empty)
{

  cout <<"Beginning event generation."<<endl;
  if (!empty){
    cout << "Requested " << fNEvents << " events." << endl;
    cout << "Progress: " << flush;
    for (Int_t i = 1; i <= fNEvents; i++){
      fSimEventNum = i;
      if (i % 1 == 0) 
        //cout <<"Event "<< i<<" of "<<fNEvents<<endl;
	cout << i << " " << flush;
      recoilEvent();
      saveEvent();
    }
    
    cout << endl;
  }
  else {
    for (Int_t i = 1; i <= fNEvents; i++){
      fSimEventNum = i;
      if (i % 20 == 0) cout << "Event "<< i <<" of " <<fNEvents << endl;
      
      emptyEvent();
      saveEvent();
    }
  }
  endRun();
}




void 
RunGenerator::normalizeGainMaps()
{

  for (int i = 0; i<fEventGen->camera()->GetEntries(); i++)
  {
    fEventGen->camera(i)->normalizeGainMap();
  }
}



void
RunGenerator::addScopeDataInfo(int ii) { _scopeDataInfo.push_back( new ScopeDataInfo() );
 }
 
ScopeDataInfo*
RunGenerator::dataInfo(int ii) { return _scopeDataInfo.at(ii); }