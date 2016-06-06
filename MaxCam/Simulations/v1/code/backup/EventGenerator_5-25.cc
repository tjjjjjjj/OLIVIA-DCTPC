#include "TCanvas.h" // debugging..

#include "../../../MaxCamSRIM.hh"
#include "../../../MaxCamConfig.hh"
#include "../../../DmtpcDataset.hh"
#include "EventGenerator.hh"

#include "TMath.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TString.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TVector3.h"
#include "TObject.h"
#include "TAxis.h"
#include "TH1.h"
#include "TF1.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <assert.h>

using namespace std;

ClassImp(EventGenerator);



double noise(TRandom3* rand, double SD) //Gaussian random value
{
  return rand->Gaus(0,SD);
}
   
double timenoise(double a1,double a2,double sz) //uncertainty of electron arrival time
{
  if (sz < 0) sz = 0;
  return a1 + a2 *pow(sz,0.5); 
}

double j(double t, double sz, double vd, double dt) //expected electron arrival time
{
  return t/dt + sz/(vd * dt);
}

double v(double En)
{
  if (En <= 0) En = 1.e-20;
  return pow(2*En/(3.728e6 / 9.e22),0.5); //magic number is alpha particle mass over c^2.
}
  
int step(int Nsamples, int& jmin, double& En, double dx, double& t, double dt, double a1, double a2, int i, double vd, int jmax, int jmax_init, vector<double>& Q, vector<double>& fRecoilZ, vector<double>& fRecoilEn, TRandom3* rand) //Computes when the electrons hit the mesh.
{

  double sz; //z-position
  double jD; //expected arrival time
  int jI; //rounded to int
  double tSD; //uncertainty of arrival time
  double El; //energy loss
  int need_space = 0; //vectors are almost out of allocated space

  sz = fRecoilZ.at(i);
  El = fRecoilEn.at(i);
  //cout << En << " En, " << El << " El." << '\n';
 
  t += dx/v(En);

  tSD = timenoise(a1,a2,sz);
  jD = j(t,sz,vd,dt);

  for (int samples = 0; samples < Nsamples; samples++)//randomly drops electrons onto the mesh at times around jI
    {
      jI = int(jD + 0.5 + noise(rand,tSD));
      if (jI < 0) jI = 1;
      jmin = (jI < jmin) ? jI : jmin;
      if (jI > (jmax - jmax_init)) need_space = 1;
  
      Q.at(jI)+=El/Nsamples;
    }  
  En -= fRecoilEn.at(i);
return need_space;
}

double vecsum(vector<double>& B, int risesteps)
{
  double sum = 0;
  for(int r=0;r<risesteps;r++)
    sum+=B.at(r);
  return sum;
}

void decay(double Tr, double Td, double dt, vector<double>& Q, vector<double>& R, int jmax) //simulates rise, decay times
{
  int risesteps = int(Tr/dt) + 1;
  double DecayConst = std::exp(-dt/Td);
  vector<double> B(risesteps,0.); //used as a "buffer" for the rise time

    for(int i=1;i<jmax;i++)
      {
	B.at(i%risesteps) = Q.at(i);
	R.at(i) = DecayConst*R.at(i-1) + vecsum(B,risesteps)/risesteps;
      }
}

void noisify(vector<double>& R, vector<double>& S, TRandom3* rand,double SD, int jmax) //adds final layer of random noise to data
{
  for(int i=0;i<jmax;i++)
{
    S.at(i) = R.at(i) + noise(rand, SD);
}
}

double round(double x, double res)//Rounds x to the nearest multiple of res
{
  if (x==0)
    return 0;
  else
    return double(int(x/res + ( x>0 ? .5 : -.5 ) ))*res;
}

void
EventGenerator::clusterfunction(double En, double dx, vector<double>& fRecoilZ,vector<double>& fRecoilEn)
{ 
  TRandom3 *rand = new TRandom3;
  /*
  const double Tr = 5.e-9; //rise time (linear, seconds)
  const double Td = 500.e-9; //decay time (exponential, seconds)
  const double Scale = .2; //scales final waveform
  const double vd = 1.646e7; //drift velocity (mm/s)
  const double scoperes = .016*256; //output voltage resolution
  const double SD = (.75*scoperes)/Scale; //final noise strength
  const double rate = 250.e6; //sampling rate (Hz)
  const double dt = 1./rate; //time resolution (s)
  const double a1 = 2.; //time noise strength (position independent)
  const double a2 = 1.; //(position dependent)
  const int Nsamples = 10; //raising this number increases the accuracy of the Gaussian time-uncertainty. 1 might be fine for small step size.
  const int zero_offset = 50; //fine-tunes the centering of the waveform on the final plot
  */
  
  double Td = scope(0)->getDecayTime();
  double Tr = scope(0)->getRiseTime();
  double vd = scope(0)->getDriftVel();
  double Scale = scope(0)->getScale();
  double scoperes = scope(0)->getScopeRes();
  double SD = scope(0)->getVoltageNoise() * scoperes/Scale;
  double rate = scope(0)->getClockRate();
  double dt = 1./rate;
  double a1 = scope(0)->getTimeNoise1();
  double a2 = scope(0)->getTimeNoise2();
  int zero_offset = scope(0)->getZeroOffset();
  int Nsamples = scope(0)->getNSamples();

  int need_space;
 
  const int jmax_init = 20000; //size of memory allocated when vectors need more space. Should be much larger than the time uncertainty caused by random noise.
  const int imin = 0; //imin corresponds to t=0
  const int imax = imin + fRecoilZ.size();
  const int j_init = int(j(0,fRecoilZ.at(0),vd,dt)); // first electron detected (approximately)
  int jmax = jmax_init + j_init; //initial length of output vectors.
  int jmin = j_init; //first electron detected (actual)
  vector<double> Q(jmax,0); //electrons hit mesh
  vector<double> R(jmax,0); //signal recorded
  vector<double> S(jmax,0); //random noise added
  double t = 0.; //time elapsed
 
   for(int i=imin;En>fMinEnergy;i++) //main simulation; generates Q (electron arrival at mesh)
 {
   need_space = step(Nsamples,jmin,En,dx,t,dt,a1,a2,i,vd,jmax,jmax_init,Q,fRecoilZ,fRecoilEn,rand); //step(...) is ugly but I didn't know classes when I wrote it
      if(need_space == 1) //allocates extra space to vectors if necessary
	{
	  jmax+=jmax_init;
	  Q.reserve(jmax);
	  R.reserve(jmax);
	  S.reserve(jmax);
	  for(int k = 0;k<jmax_init;k++)
	    {
	      Q.push_back(0.);
	      R.push_back(0.);
	      S.push_back(0.);
	    }
	}
    }
    
  decay(Tr,Td,dt,Q,R,jmax); //generates R (output after rise/decay)
  noisify(R,S,rand,SD,jmax); //generates S (output after random noise)
  
  int plotmax = scope(0)->getRecordLength()+scope(0)->getRecordPreSize();
  int presize = scope(0)->getRecordPreSize();
  int plotmin = presize + zero_offset; 
 
 //Discretizes the waveform to the scope's resolution, and centers the waveform spike near t=0
 for (int j = plotmin; j < plotmax; j++) 
   scope(0)->wf(0)->SetBinContent(j,round(Scale*S.at(j - (plotmin - jmin)),scoperes));
 for (int j = 1; j < plotmin; j++)
   scope(0)->wf(0)->SetBinContent(j,round(Scale*rand->Gaus(0,SD),scoperes));

}




EventGenerator::EventGenerator()
{

  
  hist=new TH1D("","",2000,0,2000);
  //Constructor
  fPartGen         = new ParticleGenerator();
  fChamber         = new SimChamber();
  fCamera          = new TObjArray();
  fPMT             = new TObjArray();
  fScope           = new TObjArray();
  fStepSize        = 0.0005;//mm 0.005 is default. this is an important number!! 
  fMinEnergy       = 0.2;
  //fSrimName        = "SRIM_He_in_75TorrCF4_525TorrHe_Aug2014";
  fSrimName        = "SRIM_Ne_in75TorrCF4_525TorrNe_compoundcorr_975_June2015";   
//  fSignal	   = NULL;
  //fPMTSignal       = NULL;
  fSrim 	   = NULL;
  fRnd             = new TRandom3(0);
  fUseSpacers      = false;
  fUseRadialEffect = true;
  fUseTPCImageEffect = true;
  fUseLandauEffect = false;
  fUseDoubleAlpha = true;
  //fLongStraggle    = NULL;// new TF1("longstraggle", "1/sqrt(2*TMath::Pi()*[0]^2)*exp(-(x-[1])^2/(2*[0]^2))");
  //fLongStraggle    = new TF1("longstraggle", "1/sqrt(2*TMath::Pi()*[0]^2)*exp(-(x-[1])^2/(2*[0]^2))"); // NULL;
  //fLongStraggle->SetParNames("sigma","t0");
  //fLongStraggle->SetNpx(1000);
  fTrackPhi = 0;
  fTrackTheta = 0;
  fTrueClusterArray = new TObjArray();
  fTrueClusterArray->SetOwner(kTRUE);
  fRecoilZ_min = 1e10;
  fRecoilZ_max = -1e10;
 // straggle=0.;
  fStraggle=0.;
  
  //hist=new TH1D("","",100,-10,10);


}

EventGenerator::~EventGenerator()
{
  //Destructor
  //hist->Write();
  //f->Close();
  delete fChamber;
  for (int i = 0; i < fCamera->GetEntries(); i++){
    delete fCamera->At(i);
  }
  delete fCamera;
  for (int i = 0; i < fPMT->GetEntries(); i++) delete fPMT->At(i);
  delete fPMT;
  for (int i = 0; i < fScope->GetEntries(); i++) delete fScope->At(i);
  delete fScope;
  if(fSrim) delete fSrim;
  delete fRnd;
  
  TFile *test=new TFile("test.root", "RECREATE"); 
hist->Write();
test->Close();

  
}

SimCamera*
EventGenerator::camera(int i)
{
  int ncam = fCamera->GetEntries(); 
  SimCamera* cam;
  if (i < ncam) cam = (SimCamera*) fCamera->At(i);
  else{
    cout << "Invalid camera number. There are "<<ncam <<" cameras. Returning null pointer." <<endl;
    cam = NULL;
  }
  return cam;
}

SimPMT*
EventGenerator::pmt(int i) 
{
  int npmt = fPMT->GetEntries();
  SimPMT* pmt;
  if (i < npmt) pmt = (SimPMT*) fPMT->At(i);
  else {
    cout << "Invalid PMT number.  There are "<<npmt<<" pmts.  Returning null pointer." << endl;
    pmt = NULL;
  }
  return pmt;
}

SimScope*
EventGenerator::scope(int i) 
{
  int nscope = fScope->GetEntries();
  SimScope* scope;
  if (i < nscope) scope = (SimScope*) fScope->At(i);
  else {
    cout << "Invalid Scope number: " << i << ".  There are "<<nscope<<" scopes.  Returning null pointer." << endl;
    scope = NULL;
  }
  return scope;
}

void
EventGenerator::reset()
{

  //cout<<"Hi TJ"<<endl;
  //cout << "reset()" << endl;
  int ncam = fCamera->GetEntries();
  for (int i = 0; i < ncam; i++) 
    camera(i)->resetImage();

  //cout << "reset PMTs" << endl;
  for (int ipmt=0; ipmt<fPMT->GetEntries(); ipmt++)
    pmt(ipmt)->resetWfs();

// for (int iscope=0; iscope<fScope->GetEntries(); iscope++)
  //scope(0)->resetWfs();

  //fSignal->Reset();
}

void
EventGenerator::setSrim()
{
  //creates the srim object once a name and the environmental variable
  //"MCTABLES" is set
  cout <<"Setting SRIM table " << fSrimName<<endl;
  fSrim = new MaxCamSRIM(fSrimName);
  fSrim->setStopping(fChamber->getPressure());
  //cout<<fChamber->getPressure()<<endl;

}

void
EventGenerator::setPressure(double pres)
{
  //sets the pressure needed for the srim table
  fChamber->setPressure(pres);
  fSrim->setStopping(pres);
}

double
EventGenerator::applyDiffusion(double x,double z)
{
  //Gaussian diffusion based on our diffusion constants
  double dz = z;
  if (dz<0) dz = 0;
  double sigma = sqrt(fChamber->getDiffusionConstantTerm()+dz*fChamber->getDiffusionDzTerm());
  //std::cout<<z<<" "<<dz<<" "<<fChamber->getDiffusionDzTerm()<<" "<<fChamber->getDiffusionConstantTerm()<<" "<<x<<" "<<sigma<<" "<<fRnd->Gaus(x,sigma)<<std::endl;
  return fRnd->Gaus(x,sigma);
}

double
EventGenerator::applyElectronLifetime(double x,double z)
{
  //Gaussian diffusion based on our diffusion constants
  double dz = z;
  // if (dz<0) dz = 0;
  if(dz<0. || dz>528.) x=0.;
  double signalleftover = x*exp(-dz/fChamber->getElectronLifetime());
  //std::cout<<x<<" "<<signalleftover<<" "<<exp(-z/fChamber->getElectronLifetime())<<" "<<z<<std::endl;
  return signalleftover;
}


Int_t
EventGenerator::getSequence()
{

//Big DCTPC at near hall exposure
//Sequences 12-20 (inclusive)
double exposure_near[9]={1286022.,299446.,654818.,380117.,1394001.,695633.,55359.,855258.,44421.};

//Big DCTPC at far hall exposure 
//Sequences 24-30 (inclusive)
double exposure_far[7]={88271.,81060.,940424.,637316.,465074.,1161811.,2175558.};

//Little DCTPC at far hall exposure 
//Sequences 2-18 (inclusive)
double exposure_littledctpc_far[17]={345130., 87757., 248539., 367407., 364018., 89584., 251442., 324605., 358127., 359194., 357272., 365276., 354349., 370056., 125739., 103122., 1111687.};

double totalexp=0.;
double exp=0.;

	if(particleGen()->getRunType()=="BigDCTPC_near")  
	{
	for(int i=0;i<=9;i++)
	totalexp+=exposure_near[i];
	double num=fRnd->Uniform(0,totalexp);

	for(int i=0;i<9;i++)
	{
	exp+=exposure_near[i];
	if(num<=exp)
	{
	return i;
	}
	}
	}

	if(particleGen()->getRunType()=="BigDCTPC_far")  
	{
	for(int i=0;i<=7;i++)
	totalexp+=exposure_far[i];
	double num=fRnd->Uniform(0,totalexp);
	for(int i=0;i<7;i++)
	{
	exp+=exposure_far[i];
	if(num<=exp)
	{
	return i;	
	}
	}
	}
	
	if(particleGen()->getRunType()=="LittleDCTPC_far")  
	{
	for(int i=0;i<=17;i++)
	totalexp+=exposure_littledctpc_far[i];
	double num=fRnd->Uniform(0,totalexp);
	for(int i=0;i<17;i++)
	{
	exp+=exposure_littledctpc_far[i];
	if(num<=exp)
	{
	return i;	
	}
	}
	}

return -1;	
}


void 
EventGenerator::applyGainLandauEffect(int ncam)
{
  int nbinsx = camera(ncam)->getBinsX();
  int nbinsy = camera(ncam)->getBinsY();
  double factor = fRnd->Landau(camera(ncam)->getGainLandau_1(),camera(ncam)->getGainLandau_2());
  for (int i = 1; i<= nbinsx; i++){
    for (int j = 1; j<= nbinsy; j++){
      double signal = camera(ncam)->getCCDImage()->GetBinContent(i,j);
      camera(ncam)->getCCDImage()->SetBinContent(i,j,signal*factor);
    }
  }
}

void
EventGenerator::applyRadialEffect(int ncam)
{
  //Once image is created and spacer effects are done, use this to apply the radial effect due to geometric acceptance from the mesh onto the ccd
  // double D = fChamber->getHeight();
  int nbinsx = camera(ncam)->getBinsX();
  int nbinsy = camera(ncam)->getBinsY();
  double x0 = camera(ncam)->getPositionX();
  double y0 = camera(ncam)->getPositionY();
  for (int i = 1; i<= nbinsx; i++){
    //double x = camera(ncam)->getBinPositionX(i);
    double x = (i*4)-512;
    for (int j = 1; j<= nbinsy; j++){
      //double y = camera(ncam)->getBinPositionY(j);
      double y = (j*4)-512;  
      double r = sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0));    
      double rr=pow(r,2)*pow((camera(ncam)->getWidthX()/camera(ncam)->getPixelsX())/10.,2);//convert to mm^2
	double factor= camera(ncam)->getRadialParams_1()/pow(camera(ncam)->getRadialParams_2()+rr,camera(ncam)->getRadialParams_3());
      double signal = camera(ncam)->getCCDImage()->GetBinContent(i,j);

      camera(ncam)->getCCDImage()->SetBinContent(i,j,signal*factor);
    }
  }
}


void
EventGenerator::applyTPCImageEffect(int ncam)
{

int nbinsx = camera(ncam)->getBinsX();
  int nbinsy = camera(ncam)->getBinsY();
  double x0 = camera(ncam)->getPositionX();
  double y0 = camera(ncam)->getPositionY();
  for (int i = 1; i<= nbinsx; i++){
    //double x = camera(ncam)->getBinPositionX(i);
    double x = (i*4)-512;
    for (int j = 1; j<= nbinsy; j++){
      //double y = camera(ncam)->getBinPositionY(j);
      double y = (j*4)-512;  
      double r = sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0));

		if(r>600 && particleGen()->getRunType()=="BigDCTPC_near")
        camera(ncam)->getCCDImage()->SetBinContent(i,j,0.);
	}}

}

void
EventGenerator::applySpacers(int ncam){
  //Once an image is created and spacer positions are saved, use this to add spacer effects
  vector<double> _spacers = fChamber->getSpacers();
  TString _spacerAxis = fChamber->getSpacerAxis();
  double _spacerDiameter = fChamber->getSpacerWidth();
  _spacerAxis.ToLower();
  int size = _spacers.size();
  TVectorT<double> spacerlow;
  TVectorT<double> spacerhigh;
  spacerlow.ResizeTo(size);
  spacerhigh.ResizeTo(size);
  int nbins = 0;
  int nopp =0 ;
  
  
  
  if (size > 0 && (_spacerAxis.Contains("x") || _spacerAxis.Contains("y"))){
    if (_spacerAxis.Contains("x")){
      nbins = camera(ncam)->getCCDImage()->GetNbinsY();
      nopp = camera(ncam)->getCCDImage()->GetNbinsX();
    }
    else if (_spacerAxis.Contains("y")){
      nbins = camera(ncam)->getCCDImage()->GetNbinsX();
      nopp = camera(ncam)->getCCDImage()->GetNbinsY();
    }
    for (int i = 0; i < size; i++){
		
      spacerlow[i] = _spacers[i] - _spacerDiameter / 2;
      spacerhigh[i] = _spacers[i] + _spacerDiameter/2;
      for (int j = 1; j <= nbins; j++){
	double binlow = 0, binup = 0;
	if (_spacerAxis.Contains("x")){
	  binlow = camera(ncam)->getBinPositionY(j,"lo");
	  binup = camera(ncam)->getBinPositionY(j,"hi");
	}
	else if (_spacerAxis.Contains("y")){
	  binlow = camera(ncam)->getBinPositionX(j,"lo");
	  binup = camera(ncam)->getBinPositionX(j,"hi");
	}
	double redfac=1;
	//fitted spacers to 1- 0.85*exp(-|x|^1.5 / 18), x in mm
	double spacermid = _spacers[i];
	double sFactor = 0.3784;//empirical scaling factor
	double spacerdiam = sFactor*_spacerDiameter;
	if (spacerlow[i] > binlow && spacerhigh[i] < binup){
	  redfac = 1-0.85*0.3018*spacerdiam / (binup - binlow);

	}
	
	
	
	
// 	if ((spacerlow[i]<binlow && spacerhigh[i] > binup) ||
// 		 (spacerlow[i] > binlow && spacerlow[i] < binup && spacerhigh > binup)||
// 		 (spacerhigh[i] >binlow && spacerhigh[i] < binup && spacerlow[i] < binlow)){
	
	  
	  double tempsum = 0;
	  int counter = 0;
	  double stepsize = (binup-binlow)/10;
	  for(double mm=binlow; mm<=binup;mm+=stepsize){
	    tempsum += exp(-pow(fabs((mm-spacermid)/spacerdiam),1.5)/18);
	    counter++;
	  }
	  double Sum = tempsum/counter;
	  redfac = 1-0.75*Sum;//spitz. nominally 0.55. but, this is totally arbitrary

	//}
	
	// cout<<redfac<<" "<<1-0.55*Sum<<endl;
	
	if (redfac<1){
	  for (int k = 1; k<= nopp; k++){
	    if (_spacerAxis.Contains("x")){
	  
	      double contents = camera(ncam)->getCCDImage()->GetBinContent(k,j)*redfac;
	      camera(ncam)->getCCDImage()->SetBinContent(k,j,contents);
	      //std::cout<<redfac<<std::endl;
	    }
	    else if (_spacerAxis.Contains("y")){
	      double contents = camera(ncam)->getCCDImage()->GetBinContent(j,k)*redfac;
	      camera(ncam)->getCCDImage()->SetBinContent(j,k,contents);
	      //std::cout<<contents<<std::endl;
	    }
	  }//for (int k = 0; k<= ...: loop over second axis
	}//if (redfac<1)
      }//for (int j = 1; j<= ...:loop over first axis

    }//for (int i = 0; i<size...: loop over spacers
  }//if (axis is x or y)
}//applySpacers()

bool
EventGenerator::propagateRecoil()
{

  //Propagate a recoil product by a step size
  double dR = getStepSize();
  assert(fTempRecoil.Mag());
  assert(fSrim->getStoppingVsEnergy());
  double E = fTempEnergy;
  //get straggling lengths

  //fSrim->getTable()->Scan("","");

  //fSrim->getRangeVsEnergy()->Print(); 
  // fSrim->getStoppingVsEnergy(false)->Print();

  double srim_range = fSrim->ROOTv24TGraphEval(fSrim->getRangeVsEnergy(),E); //  fSrim->getRangeVsEnergy()->Eval(E);
  //double lstraggle = fRnd->Gaus(1,fSrim->getLStraggleVsEnergy()->Eval(E)/srim_range);
  
//   cout<<fSrim->getTStraggleVsEnergy()->Eval(E)<<" "<<E<<" "<<srim_range<<endl;
//   double lstraggle = fRnd->Gaus(1,fSrim->getLStraggleVsEnergy()->Eval(E)/srim_range);
//   double tstraggle = fRnd->Gaus(0,fSrim->getTStraggleVsEnergy()->Eval(E)/srim_range)/sqrt(2);
//   double nstraggle = fRnd->Gaus(0,fSrim->getTStraggleVsEnergy()->Eval(E)/srim_range)/sqrt(2);

  
  double srim_totalrange= fSrim->ROOTv24TGraphEval(fSrim->getRangeVsEnergy(),fPartGen->getRecoilEnergy());
 
 

double lstraggle = 1.+(fRnd->Gaus(0,fSrim->getLStraggleVsEnergy()->Eval(E))/sqrt(srim_range/dR)); 
double tstraggle = fRnd->Gaus(0,fSrim->getTStraggleVsEnergy()->Eval(E))/sqrt(srim_range/dR);
double nstraggle = tstraggle;
   
  //cout<<E<<" "<<fSrim->getTStraggleVsEnergy()->Eval(E)<<" "<<srim_totalrange<<" "<<dR<<endl;
  // cout<<straggle<<endl;
  
  fStraggle+=tstraggle;
  
//   cout <<E<< " fSrim->getRangeVsEnergy()->Eval(E)=" << fSrim->getStoppingVsEnergy()->Eval(E) << " "<< fSrim->ROOTv24TGraphEval(fSrim->getStoppingVsEnergy(),E)<<endl;
  //cout << "fSrim->ROOTv24TGraphEval(fSrim->getRangeVsEnergy(),E)=" << fSrim->ROOTv24TGraphEval(fSrim->getRangeVsEnergy(),E) << endl;
  //getchar();
  

  //make a vector using the straggling lengths and renormalize
  TVector3 longVec = fTempRecoil;
  TVector3 tranVec = longVec.Orthogonal();
  TVector3 normVec = longVec.Cross(tranVec);
  TVector3 stepVec = longVec*lstraggle + tranVec * tstraggle + normVec * nstraggle;
  stepVec = stepVec.Unit();
  //get step stopping energy
  double dE, dE_elec, dE_scint, omegaB;
  //double skew=((getSkew()/fPartGen->getRecoilEnergy())*E)+(1-(getSkew()/2.));
  
  
  //fPartGen->getRecoilEnergy()
//use this for the full hist. much slower!
//   dE=fSrim->getdedx(E)*dR;
//   dE_elec=dE;
//   dE_scint=dE;
  
//He correction  
double x_He=1.;
double corr_He=1.;
if( (E/4.) < (75.*2.) )
{
x_He=(E/4.)/(25.*2.);
corr_He=0.5*((1.36*sqrt(x_He))-(0.16*pow(x_He,1.5)));
}

//C correction  
double x_C=1.;
double corr_C=1.;
if( (E/4.) < (75.*6.) )
{
x_C=(E/4.)/(25.*6.);
corr_C=0.5*((1.36*sqrt(x_C))-(0.16*pow(x_C,1.5)));
}

//F correction  
double x_F=1.;
double corr_F=1.;
if( (E/4.) < (75.*9.) )
{
x_F=(E/4.)/(25.*9.);
corr_F=0.5*((1.36*sqrt(x_F))-(0.16*pow(x_F,1.5)));
}

//cout<<E<<" "<<corr_F<<endl;


// double omegaB2_He_alt=corr_He*0.157*(0.2415*0.000465)*(fStepSize*1000.*100.)*pow(2.,2.)*2/4.; 
// double omegaB2_C_alt=corr_C*0.157*(0.1035*0.000465)*(fStepSize*1000.*100.)*pow(2.,2.)*6/12.;
// double omegaB2_F_alt=corr_F*0.157*(0.3333*0.000465)*(fStepSize*1000.*100.)*pow(2.,2.)*9/19.;
// double omegaB_alt=sqrt(omegaB2_He_alt+omegaB2_C_alt+omegaB2_F_alt);

//http://www.slideshare.net/ysina/energy-loss-and-energy-straggling-a-presentation-by-younes-sina
double omegaB2_He=corr_He*0.26*(0.5834*28.963)*(fStepSize/10.)*pow(2.,2.)*2; 
double omegaB2_C=corr_C*0.26*(0.0833*28.963)*(fStepSize/10.)*pow(2.,2.)*6;
double omegaB2_F=corr_F*0.26*(0.3333*28.963)*(fStepSize/10.)*pow(2.,2.)*9; 
 
 
// double corr_He= 
 
 omegaB=sqrt(omegaB2_He+omegaB2_C+omegaB2_F);



 //cout<<skew<<" "<<getSkew()<<" "<<E<<endl;
//use this for poor dE/dx estimate. much faster!  
// dE = (fSrim->getStoppingVsEnergy()->Eval(E)*dR)+fRnd->Gaus(0.0,omegaB);
dE = (fSrim->getStoppingVsEnergy()->Eval(E)*dR);
// if(dE<0.)
// {
// dE=0.;
// }
dE_elec = fSrim->getStoppingVsEnergy(false)->Eval(E)*dR;
dE_scint = (dE-dE_elec)*fChamber->getNuclearScintillation() + dE_elec*fChamber->getElectricScintillation();
  

// std::cout<<E<<" "<<dE/dR<<" "<<fSrim->getStoppingVsEnergy()->Eval(E)*dR<<" "<<dE<<" "<<omegaB<<" "<<std::endl;
  
  //Save (x,y,n) in TVectorTs to be later entered into images
  fRecoilX.push_back(fTempPosition.X());
  fRecoilY.push_back(fTempPosition.Y());
  fRecoilZ.push_back(fTempPosition.Z());
  if (fTempPosition.Z() > fRecoilZ_max) {
    fRecoilZ_max = fTempPosition.Z();
  }
  if (fTempPosition.Z() < fRecoilZ_min) {
    fRecoilZ_min = fTempPosition.Z();
  }

  fRecoilEn.push_back(dE_scint);//default!!
  //std::cout<<std::endl;
  //fTempRecoil.Print("");
  fTempRecoil = fTempRecoil*fTempEnergy-stepVec*dE;//doesn't seem correct to me
  
  //fTempRecoil = (fTempRecoil-stepVec)*(fTempEnergy-dE);
  //stepVec.Print("");
  //fTempRecoil.Print("");
  
  TVector3 testrecoil;
  testrecoil = (fTempRecoil-stepVec); 
  //testrecoil.Print("");
  
  //fTempEnergy = fTempRecoil.Mag();
  //fTempRecoil = fTempRecoil.Unit();
  fTempEnergy = E-dE;
  fTempRecoil = testrecoil.Unit();
  fTempPosition += stepVec*dR*fLengthScaleFactor;
  return E-dE > fMinEnergy && E-dE > dE;
}

void
EventGenerator::driftElectrons()
{
  //Drift electrons held in each element of vector objects
  //This is done simply by applying the diffusion function for each element
  double xdiff, ydiff;
  int count = fRecoilX.size();
  for (int i = 0; i < count; i++){
    xdiff = applyDiffusion(fRecoilX[i],fRecoilZ[i]);
    ydiff = applyDiffusion(fRecoilY[i],fRecoilZ[i]);
    fRecoilX[i] = xdiff;
    fRecoilY[i] = ydiff;
  }
}

void
EventGenerator::makeImage(int ncam)
{
  //Puts the vectors (x,y,n) onto the ccd image TH2F
  int count = fRecoilX.size();
  double countBfSpacers = 0;
  double totalen = 0;
  int pixperbin = camera(ncam)->getPixelPerBin();
  //cout << fRecoilX[count-1] << "," << fRecoilY[count-1] << endl;
  
  // for storing the cluster in the image in a MaxCamClusterImage object that
  // if the user wants, gets stored in the Simulation tree...
  vector<int> cluster;
  double gain_fluc=1.;
 
 // double gain_fluc=applyGainFluctuation();
//12   13   14   15   16  17  18  19  20  
// 0.65,0.56,0.56,0.74,.65,.83,.90,.65,.65
 
 double gain_factor_near[9]={1./0.65,1./0.56,1./0.56,1./0.74,1./0.65,1./0.83,1./0.90,1./0.65,1./0.65};
 
 double gain_factor_far[7]={1./.52,1./.50,1./.46,1./.52,1./.52,1./.52,1./.52};
 
 double gain_factor_littledctpc_far[17]={1./0.955114,1./1.06241,1./1.07105,1./1.07669,1./1.02827,1./1.15574,1./1.18828,1./1.08779,1./0.998239,1./1.05118,1./1.03128,1./1.02639,1./0.992746,1./0.925486,1./1.1664,1./1.05791,1./0.989169};
 
  if(particleGen()->getRunType()=="BigDCTPC_near")
  gain_fluc=gain_factor_near[fPartGen->getSeqNum()];
  
  if(particleGen()->getRunType()=="BigDCTPC_far")
  gain_fluc=gain_factor_far[fPartGen->getSeqNum()];
 
  if(particleGen()->getRunType()=="LittleDCTPC_far")
  gain_fluc=gain_factor_littledctpc_far[fPartGen->getSeqNum()];
 
  for (int i = 0; i < count; i++){
    if (camera(ncam)->isInImage(fRecoilX[i],fRecoilY[i])){
      if (fRecoilZ[i] < 0 || fRecoilZ[i] > fChamber->getDriftLength()) continue; //We exit the detector. 
      double xpixel = camera(ncam)->getPixelX(fRecoilX[i]);
      double ypixel = camera(ncam)->getPixelY(fRecoilY[i]);
      int xbin = static_cast<int> (xpixel/pixperbin);
      int ybin = static_cast<int> (ypixel/pixperbin);
      if (xpixel / pixperbin - xbin != 0) xbin++;
      if (ypixel / pixperbin - ybin != 0) ybin++;
      //double npho = fRnd->Poisson(camera(ncam)->getGain(xbin,ybin)*fRecoilEn[i]);

//cout<<"seq num "<<fPartGen->getSeqNum()<<" "<<1./gain_fluc<<endl;

      double npho = camera(ncam)->getGain(xbin,ybin)*fRecoilEn[i];
      
      // cout<<gain_fluc<<endl;
      
      //cout<<"1 "<<npho<<" "<<fRecoilZ[i]<<endl;
      //cout<<"1.5 "<<npho/applyElectronLifetime(npho,fRecoilZ[i])<<endl;
      npho=applyElectronLifetime(npho,fRecoilZ[i]);
      //cout<<"2 "<<npho<<" "<<fRecoilZ[i]<<endl;
      
      //gain_fluc=1./.65;
      npho*=gain_fluc;
      
      //cout<<"3 "<<npho<<" "<<fRecoilZ[i]<<endl;
//       if (camera(ncam)->getEMGain() > 1.0) 
//       npho = fRnd->Poisson(npho*camera(ncam)->getEMGain());
      totalen += fRecoilEn[i];
      countBfSpacers += npho;
      //cout<<camera(ncam)->getCCDImage()->GetBinContent(xbin,ybin)<<endl;
      camera(ncam)->getCCDImage()->SetBinContent(xbin,ybin,camera(ncam)->getCCDImage()->GetBinContent(xbin,ybin)+npho);
//       std::cout<<xbin<<" "<<ybin<<" "<<camera(ncam)->getCCDImage()->GetBinContent(xbin,ybin)+npho<<std::endl;
      // add this bin to the cluster
      cluster.push_back(camera(ncam)->getCCDImage()->GetBin(xbin,ybin));
    }
  }



for (int i=1;i<=camera(ncam)->getCCDImage()->GetNbinsX();i++)
{
for (int j=1;j<=camera(ncam)->getCCDImage()->GetNbinsY();j++)
{

if(camera(ncam)->getCCDImage()->GetBinContent(i,j)<=0)
continue;
 

 
double binc=fRnd->Poisson(camera(ncam)->getPhotonsADU()*camera(ncam)->getCCDImage()->GetBinContent(i,j))/camera(ncam)->getPhotonsADU();//~1.3 photons/ADU (poisson stats should go according to #photons collected i think)



camera(ncam)->getCCDImage()->SetBinContent(i,j,binc);

}
}

  
  // done making the pre-noise and bias image
  
  // push this image onto the stack of pre-noise and bias images
  //  cout << "camera(ncam)->getCCDImage()=" << camera(ncam)->getCCDImage() << endl;
  //  cout << "camera(ncam)->GetMean()=" << camera(ncam)->getCCDImage()->GetMean() << endl;
  // maybe not good to make a copy, but...also, the TDatime is just empty right now...
  MaxCamClusterImage *this_cluster=new MaxCamClusterImage(new TH2F(*(camera(ncam)->getCCDImage())), new TDatime());
  // add the cluster to it
  // there are duplicates in our bin vector; kill them
  
  vector<int>::iterator it;
  /*
  cout << "cluster contains:";
  Int_t counter=0;
  for (it=cluster.begin(); it!=cluster.end(); ++it){
    cout << " " << *it;
    ++counter;
    if(counter>100) break;
  }
  cout << endl;
  */

  sort(cluster.begin(),cluster.end());

  /*
  cout << "post-sort cluster contains:";
  counter=0;
  for (it=cluster.begin(); it!=cluster.end(); ++it){
    cout << " " << *it;
    ++counter;
    if(counter>100) break;
  }
  cout << endl;
  */

  // using default comparison:
  it = unique (cluster.begin(), cluster.end()); // 10 20 30 20 10 30 20 20 10
                                                  //                ^
  
  cluster.resize( it - cluster.begin() );       // 10 20 30 20 10
  
  // using predicate comparison:
  unique (cluster.begin(), cluster.end());   // (no changes)

  /*
  cout << "post-uniqueing cluster contains:";
  counter=0;
  for (it=cluster.begin(); it!=cluster.end(); ++it){
    cout << " " << *it;
    ++counter;
    if(counter>100) break;
  }
  cout << endl;
  */

  this_cluster->addCluster(cluster);
  fTrueClusterArray->AddAtAndExpand(this_cluster,ncam);
  //  cout << "fTrueClusterArray->At(ncam)->getImage()->GetMean()=" 
  //       << ((MaxCamClusterImage*)fTrueClusterArray->At(ncam))->getImage()->GetMean() 
  //       << endl;
  //  cout << "fTrueClusterArray->GetEntries()=" << fTrueClusterArray->GetEntries() << endl;

  camera(ncam)->setCountBeforeSpacers(countBfSpacers);
}

//void
//EventGenerator::makeElectronicSignal()
//{
//  int count = fRecoilX.size();
//  double elperkeV = fChamber->getElectronPerkeV();
//  double atten = fChamber->getAttenuation();
//  double vDrift = fChamber->getDriftVelocity();
//  double diffconst = fChamber->getLongDiffusionConstantTerm();
//  double diffdz = fChamber->getLongDiffusionDzTerm();
//  for (int i = 0; i<count; i++){
//    double z = fRecoilZ[i];
//    double E = fRecoilEn[i];
//    //Poisson distribution of electrons at mesh
//    double n = fRnd->PoissonD(E*elperkeV*exp(-1*atten*z));
//    double t0 = z/vDrift;
//    double sigma = sqrt(diffconst+z*diffdz);
//    fLongStraggle->SetParameters(sigma,t0);
//    //
//    double xmin = t0 - 7*sigma;
//    double xmax = t0 + 7*sigma;
//    fLongStraggle->SetRange(xmin,xmax);
//    fSignal->FillRandom("longstraggle",(int)n);
//  }
//  fSignal->Scale(-1);
//
//}

//void
//EventGenerator::setElectronicSignal()
//{
//  double min = 0;
//  double max = 1.5*fChamber->getDriftLength()/fChamber->getDriftVelocity();
//  int nbins = (int) (max-min) / ((int) fTimeResolution) + 1;
// max = nbins*fTimeResolution;
//  fSignal = new TH1F("signal","Electronic Readout",nbins,min,max); 
//  fSignal->GetYaxis()->SetTitle("Q/e");
//  fSignal->GetXaxis()->SetTitle("Time [ns]");
//}


//void
//EventGenerator::applyPMTQE() {
//  fPhotonSpectrum->Multiply(fPMT->getPMTQE());
//}

Double_t
EventGenerator::solidAngleSphericalCap(Double_t D, Double_t h) {
  // Solid angle of cone with vertex angle 2*theta
  //    = 2*pi*(1-cos(theta))
  //
  // For a disk of diameter D that sits a distance h above a point
  // tan(theta) = 0.5*D/d
  //
  return TMath::TwoPi()*(1.0-h/TMath::Sqrt(h*h+0.25*D*D));
}

Double_t
EventGenerator::acceptance(Double_t D, Double_t h) {
  return 0.25*TMath::InvPi()*solidAngleSphericalCap(D,h);
}

void
EventGenerator::makePMTSignal(Int_t ipmt, Int_t itrig)
{
  //cout << "makePMTSignal()" << endl;
    
  Int_t count = fRecoilX.size();
  Double_t elperkeV  = fChamber->getElectronPerkeV();
  Double_t atten     = fChamber->getAttenuation(); 
  Double_t vDrift    = fChamber->getDriftVelocity();
  Double_t diffconst = fChamber->getLongDiffusionConstantTerm();
  Double_t diffdz    = fChamber->getLongDiffusionDzTerm();

  /* FIXME:  photons per electron is pressure dependent... add to config file */
  Double_t photonsPerElectron = 0.3;

  Double_t photonPerKeV  = elperkeV*photonsPerElectron;
  Double_t pmtDiameter   = pmt(ipmt)->getDiameter(); 
  Double_t distToPMT     = pmt(ipmt)->getDistance();
  Double_t pmtAcceptance = acceptance(pmtDiameter, distToPMT);
  //0.0625*pmtDiameter*pmtDiameter/(distToPMT*distToPMT); // simple approx.
  //cout << "pmtAcceptance = " << pmtAcceptance << endl;

  /* FIXME: should account for 1/r^2 and projected size of PMT face for points that are not sub-PMT */

  /* FIXME:  spectral information ignored */
  /* FIXME:  PMT gain ignored */
  /* FIXME:  PMT QE ignored */
  /* FIXME:  Chamber window transmission (wl dependent) ignored */
  /* FIXME:  Amplification noise in PMT ignored */
  /* FIXME:  Digitization noise ignored (amplifier noise, scope noise etc.) */

  Double_t photonsPerKeVAtPMT = photonPerKeV*pmtAcceptance;

  // debug
  cout << "count = " << count << endl; 
  //TGraph *gz = new TGraph();  gz->SetTitle("energy deposition vs. z;z [mm]; dE/dz [keV/mm]");
  //TGraph *gx = new TGraph(); gx->SetTitle("energy deposition vs. x;x [mm]; dE/dx [keV/mm]");
  //TGraph *gy = new TGraph(); gy->SetTitle("energy deposition vs. y;y [mm]; dE/dy [keV/mm]");
  for (int i = 0; i<count; i++){
    double z = fRecoilZ[i];
    double E = fRecoilEn[i];
    //Poisson distribution of electrons at mesh (attachment included)
    //What about noise in the amplification process?? Landau?  Polya?
    double n  = fRnd->PoissonD(E*photonsPerKeVAtPMT*exp(-atten*z));
    double t0 = (z-fRecoilZ_min)/vDrift;
    double sigma = sqrt(diffconst+z*diffdz);
    /* FIXME: add check that itrig is not larger than length of SimPMT->wfs() */
    /* FIXME: longitudinal diffusion should be inverse gaussian, not gaussian */
    /* FIXME: data should be stored as TH1C (8-bit)... */
    for (int ipt=0; ipt<n; ipt++) {
      pmt(ipmt)->wf(itrig)->Fill(fRnd->Gaus(t0, sigma));
    }
    //gz->SetPoint(i, z, E);
  }
  
  //TCanvas *cc = new TCanvas("cc","cc", 800, 800); // debug
  //cc->Divide(1,2);
  //cc->cd(1);
  //gz->Draw("al");
  //cc->cd(2);
  //pmt(ipmt)->wf(itrig)->Draw();
  //cc->Update();
  //getchar();
}

void
EventGenerator::runEvent(bool resetImage)
{
  //Runs an event: Use only once all parameters and options are set
  //i.   If desired, reset image
  //ii.  Make recoil
  //iii. Propagate Recoil
  //iv.  Drift Electrons
  //v.  Create Image
  //vi.   Apply Spacers
  //vii.  Apply Radial Effect
  //viii. Apply Noise
  if (resetImage) reset();
  fPartGen->generateRecoil();
  
  //cout<<"run event"<<endl;
  
  fPartGen->setSeqNum(getSequence());
  
  fTempRecoil = fPartGen->recoilVector();
  fTempPosition = fPartGen->recoilPosition();
  fTempEnergy = fPartGen->getRecoilEnergy();
  //Generate straggling
  double E = fTempEnergy;
  double Initial_Energy = E;
  double R = fSrim->getRangeVsEnergy()->Eval(E);
//   double lstraggle = fRnd->Gaus(R,fSrim->getLStraggleVsEnergy()->Eval(E));
//   double tstraggle = fRnd->Gaus(0,fSrim->getTStraggleVsEnergy()->Eval(E));
//   double nstraggle = fRnd->Gaus(0,fSrim->getTStraggleVsEnergy()->Eval(E));
  
  double lstraggle= R;
  double tstraggle=0.;
  double nstraggle=0.;
  
  TVector3 longVec = fTempRecoil;
  TVector3 tranVec = longVec.Orthogonal();
  TVector3 normVec = longVec.Cross(tranVec);
  TVector3 stepVec = longVec*lstraggle + tranVec*tstraggle+normVec*nstraggle;
  fLengthScaleFactor = stepVec.Mag() / R;//Rescale all lengths by this: track loses same amount of energy but over a different total length
  fTempRecoil = stepVec.Unit();
  fTrackPhi = stepVec.Phi();
  fTrackTheta = stepVec.Theta();
//   fTempEnergy = (E*lstraggle/R);//definitely not really correct but
//   //included in MaxCamMC

   //setSkew(-TMath::Abs(fRnd->Gaus(0.0,.2)));
  //cout<<fStraggle<<endl;
  //hist->Fill(fStraggle);
//   new TCanvas;
//   hist->Draw();
  
  
  
  fStraggle=0.;
  
  
  while(propagateRecoil()) {}
  
 //////////////////////////
 
if (fUseDoubleAlpha)
 {
  

  fTempRecoil = fPartGen->recoilVector2();
  fTempPosition = fPartGen->recoilPosition();
  fTempEnergy = fPartGen->getRecoilEnergy2();
  //Generate straggling
  E = fTempEnergy;
  
  R = fSrim->getRangeVsEnergy()->Eval(E);
//   double lstraggle = fRnd->Gaus(R,fSrim->getLStraggleVsEnergy()->Eval(E));
//   double tstraggle = fRnd->Gaus(0,fSrim->getTStraggleVsEnergy()->Eval(E));
//   double nstraggle = fRnd->Gaus(0,fSrim->getTStraggleVsEnergy()->Eval(E));
  
  lstraggle= R;
  tstraggle=0.;
  nstraggle=0.;
  
  longVec = fTempRecoil;
  tranVec = longVec.Orthogonal();
  normVec = longVec.Cross(tranVec);
  stepVec = longVec*lstraggle + tranVec*tstraggle+normVec*nstraggle;
  fLengthScaleFactor = stepVec.Mag() / R;//Rescale all lengths by this: track loses same amount of energy but over a different total length
  fTempRecoil = stepVec.Unit();
  fTrackPhi = stepVec.Phi();
  fTrackTheta = stepVec.Theta();
//   fTempEnergy = (E*lstraggle/R);//definitely not really correct but
//   //included in MaxCamMC

   //setSkew(-TMath::Abs(fRnd->Gaus(0.0,.2)));
  //cout<<fStraggle<<endl;
  //hist->Fill(fStraggle);
//   new TCanvas;
//   hist->Draw();
  
  
  
  fStraggle=0.;
  
  
  while(propagateRecoil()) {}
}  
  
  ////////////////////
  
  
  
 
  
  
  Int_t npoints = fRecoilX.size();
//get length of track
  fLength = sqrt(pow(fRecoilX[0] - fRecoilX[npoints-1],2)+pow(fRecoilY[0]-fRecoilY[npoints-1],2) + pow(fRecoilZ[0] - fRecoilZ[npoints-1],2));
  fZLength = sqrt(pow(fRecoilZ[0] -fRecoilZ[npoints-1],2));
  
  
  // std::cout<<"Z Length: "<<fZLength<<std::endl;
//  makeElectronicSignal();

  /////////////////////////////////////////////////////////
  // Generate PMT signals
  for (int ipmt=0; ipmt<fPMT->GetEntries(); ipmt++) {
    Int_t itrig=0;
    makePMTSignal(ipmt,itrig);
  }
  /////////////////////////////////////////////////////////

  driftElectrons();
  // make the cluster array the correct size to hold the images
  for (int i = 0; i < fCamera->GetEntries(); i++){
    makeImage(i);
    if (fUseSpacers) applySpacers(i);
    camera(i)->setCountBeforeRadialEffect(camera(i)->getCCDImage()->Integral());
    if (fUseLandauEffect) applyGainLandauEffect(i);
    if (fUseRadialEffect) applyRadialEffect(i);
    
    if (fUseTPCImageEffect) applyTPCImageEffect(i);
    
    camera(i)->setFinalCount(camera(i)->getCCDImage()->Integral());
    // save MaxCamClusterImage here, before noise is applied 
    camera(i)->applyNoise();
  }
  
  clusterfunction(Initial_Energy, fStepSize, fRecoilZ, fRecoilEn);

  fRecoilX.clear();
  fRecoilY.clear();
  fRecoilZ.clear();
  fRecoilEn.clear();
}





