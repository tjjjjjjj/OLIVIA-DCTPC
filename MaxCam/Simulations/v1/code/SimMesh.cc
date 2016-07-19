#include "SimMesh.hh"
#include "TObject.h"
#include "SimScope.hh"
#include "EventGenerator.hh"
#include "TVector.h"

ClassImp(SimMesh)

  SimMesh::SimMesh() {/*empty constructor*/}

SimMesh::SimMesh(int plot_offset_EG, int imin_EG, double Initial_Energy_EG, double StepSize_EG, SimScope* scope_EG, double fMinEnergy_EG, vector<double>& fRecoilZ_EG, vector<double>& fRecoilEn_EG, SimChamber* fChamber_EG, int New_Plot_EG)
{//Constructor
  plotname = "null";
  new_plot = New_Plot_EG;
  mass = 3.727e6 / (9e22); //MAGIC NUMBER! alpha mass in keV/c^2, c in mm/s
  fRecoilEn = fRecoilEn_EG; //Energy loss vector
  fRecoilZ = fRecoilZ_EG; //Position vector
  fMinEnergy = fMinEnergy_EG; //Energy threshold for end of simulation
  scope = scope_EG; //pointer to scope information and waveform histogram
  fChamber = fChamber_EG; //contains variables relevant to the TPC
  Td = scope->getDecayTime(); 
  Tr = scope->getRiseTime();
  vd = scope->getDriftVel();
  Scale = scope->getScale(); //final scale factor for the y-axis
  scoperes = scope->getScopeRes(); //resolution of scope
  SD = scope->getVoltageNoise() * scoperes/Scale; //amplitude of y-axis noise
  atten = scope->getAttenuation();
  rate = scope->getClockRate();
  dt = 1./rate;
  a1 = scope->getTimeNoise1(); //smearing of electron signal due to uncertainty in electronics
  a2 = scope->getTimeNoise2(); //smearing of electron signal due to diffusion
  zero_offset = scope->getZeroOffset(); //centers plot on x-axis
  zero_offset_var = scope->getZeroOffsetVar(); //adjusts zero_offset based on distance from mesh
  Nsamples = scope->getNSamples(); //precision of electron smearing
  vertical_offset = scope->getVerticalOffset() * 1/Scale; //centers plot on y-axis
  jmax_init = 20000; //MAGIC NUMBER! value isn't very important. Could theoretically cause out-of-range errors if it is too low.
  dx = StepSize_EG; //step size of simulation
  Initial_Energy = Initial_Energy_EG; //inital alpha energy
  En = Initial_Energy; //current alpha energy
  rand = new TRandom3(0); //RNG
  imin = imin_EG; //start point of position vector
  imax = fRecoilZ.size(); //end of position vector
  j_init = int(j(fRecoilZ.at(imin))); //expected arrival time of first emitted electron
  jmin = j_init; //actual arrival time of first arriving electron
  zero_offset += int(pow(fRecoilZ.at(imin),0.5)*zero_offset_var*a2);
  jmax = jmax_init + j_init; //space allocated to Q,R
  Q = vector<double>(jmax,0); //vector that stores electron arrival times
  R = vector<double>(jmax,0); //vector that stores electron signal after rise/decay
  t = 0; //current time
  plot_offset = plot_offset_EG; //used for centering histogram

  double sz = fRecoilZ.at(imin);
  double tSD = timenoise(sz);
  MakeGaussian(ker,tSD,Nsamples);

  int runnumberfound = 0;
  runnum = 0;

  std::stringstream nametmp;
  nametmp.str("");

  for (int n = 0; runnumberfound == 0; n++)
    {
      std::stringstream nametmp2;
      nametmp2.str("");
      nametmp2 << "waveforms/wfs" << n << ".root";

      TString nametmp3 = nametmp2.str();

      TFile* tmpfile = new TFile(nametmp3,"read");
      //tmpfile->Open();

      if(!tmpfile->IsOpen())
      {
	cout << "Ignore this error message. Idk how to make it go away but it's not a problem.";
	runnumberfound = 1;
        runnum = n-1;
	if (new_plot == 1) runnum += 1;
      }	
      tmpfile->Close();
    }
  nametmp << "waveforms/wfs" << runnum << ".root";
  rootfilename = nametmp.str();
}

SimMesh::~SimMesh(){/*destructor*/}


void SimMesh::simulate()// simulates waveforms on an initialized instance of SimMesh
{
  TH1D* firsthist = new TH1D("first_"+plotname,"",16384,0,16384);
  for(int i = imin; En>fMinEnergy && i < imax-1 && fRecoilZ.at(i)>0 && fRecoilZ.at(i)<fChamber->getDriftLength(); i++)
    //while the particle still has energy, the index i is in range, and the particle is in the chamber
    {
      final_position = i+1;
      SimMesh::step(i,firsthist);
      if(need_space == 1) SimMesh::allocate();
    }

  TString fileoption = (new_plot == 1) ? "recreate" : "update";
  TFile* waveforms = new TFile(rootfilename,fileoption);
  firsthist->Write();
  delete firsthist;
  waveforms->Close();

  SimMesh::decay();
  SimMesh::plot();
}

double SimMesh::normaldist(double peak, double x, double SD)
{return exp((x-peak)*(peak-x)/(2*SD*SD)) / (SD * pow(3.141592,0.5));}

void SimMesh::MakeGaussian(vector<double>& kernel, double SD, double acc)
{
  kernel.clear();
  int size = int(acc*SD+1);
  for (int i = -size; i < size + 1; i++)
    kernel.push_back((4*normaldist(0.,double(i),SD)+normaldist(0.,double(i)+.5,SD)+normaldist(0.,double(i)-.5,SD))/6.);
  normalize(ker);
}

void SimMesh::normalize(vector<double>& vec)
{
  double sum = vecsum(vec,vec.size());
  if (sum == 0) sum = 1;
  for (int i = 0; i < vec.size(); i++)
    vec.at(i) /= sum;
}

void SimMesh::step(int i, TH1D* firsthist) //runs one step of electron ionization and drift
{
  double sz; //z-position                                    
  double jD; //expected arrival time of electrons   
  int jI; //jD rounded to int                                    
  double tSD; //uncertainty of arrival time             
  double El; //energy loss at current step 
  need_space = 0; //records whether vectors are almost out of allocated space  
  sz = fRecoilZ.at(i);
  //if (i%1000 == 0) cout << sz << " | " << i << "|" << vel() << endl;
  El = fRecoilEn.at(i);
  t += dx/vel();
  tSD = timenoise(sz);
  jD = j(sz);
  jI = int(jD+.5);
  if (i%20 == 0) MakeGaussian(ker,tSD,Nsamples);
  int vecsize = (ker.size()-1)/2;

  jmin = (jI - vecsize < jmin) ? jI-vecsize : jmin;
  if (jI + vecsize + 1000 >= jmax) need_space = 1;
  firsthist->SetBinContent(jI,firsthist->GetBinContent(jI)+El);
  
  for (int i = -vecsize; i < vecsize+1; i++)
    {
      if (jI + i > 0)
	Q.at(jI+i)+= exp(-1. * sz * atten) * El * ker.at(i+vecsize); //record electrons' signal
    }
  En -= /*pow(10,1.5)*/El;
}

double SimMesh::round(double x, double res) //rounds x to nearest multiple of res
{double offset;
  if (x <0.)
    offset = -.5; //needed to actually round x, instead of just truncating it
  if (x==0.)
    offset =  0.;
  if (x >0.)
    offset =  .5;
  return res*(double(int(x/res + offset)));}

double SimMesh::vecsum(vector<double>& B, int risesteps) //computes sum of vector's elements
{
  double sum = 0;
  for(int r=0;r<risesteps;r++)
    sum+=B.at(r);
  return sum;
}

void SimMesh::allocate() //allocates space to Q,R
{
  jmax+=jmax_init;
  Q.reserve(jmax);
  R.reserve(jmax);
  for(int k = 0;k<jmax_init;k++)
    {
      Q.push_back(0.);
      R.push_back(0.);
    }
}

void SimMesh::decay() //simulates rise, decay times
{
  TString fileoption = (new_plot == 1) ? "update" : "update";
  TFile* waveforms = new TFile(rootfilename,fileoption);
  TH1D* bragghist = new TH1D("bragg_"+plotname,"",16384,0,16384);
  for(int j=0; j<16384; j++)
    bragghist->SetBinContent(j,Q.at(j));
  bragghist->Write();
  delete bragghist;
  waveforms->Close();


  int risesteps = int(Tr/dt) + 1;
  double DecayConst = std::exp(-dt/Td);
  vector<double> B(risesteps,0.); //used as a "buffer" for the rise time  
  for(int i=1;i<jmax;i++)
    {
      B.at(i%risesteps) = Q.at(i); //stores the risesteps most recent elements of Q
      R.at(i) = DecayConst*R.at(i-1) + vecsum(B,risesteps)/risesteps;
    }
}

void SimMesh::noisify() //adds random noise to waveform. [Currently not used. Has been replaced by plotfinish(). ]
{
  for(int i=1;i<jmax;i++)
    R.at(i) = R.at(i) + noise(SD);
}

void SimMesh::plot() //plots results
{
  int presize = scope->getRecordPreSize();
  int plotmax = presize + scope->getRecordLength();
  int plotmin = presize + zero_offset;

  if (new_plot == 1) //Re-centers and resets waveform for the first particle in a batch
    {
      for (int j = 0; j < plotmax; j++)
        scope->wf(0)->SetBinContent(j,0.);

      plot_offset = jmin - plotmin;
    }
  TH1D* hist = new TH1D(plotname,"",16384,0,16384);

  for (int j = 1; j < plotmax; j++)
    {  
      if (j + plot_offset > 0 && j + plot_offset < R.size()) 
	{
	  scope->wf(0)->SetBinContent(j,scope->wf(0)->GetBinContent(j)+Scale*R.at(j + plot_offset));
	  hist->SetBinContent(j,Scale*(vertical_offset+R.at(j+plot_offset)));
	}
    }

  TH1I* bizarre_workaround = new TH1I("plot_offset","",1,0,1);
  TH1D* bizarre_workaround2 = new TH1D("vertical_offset","",1,0,1);
  bizarre_workaround->SetBinContent(0,plot_offset);
  bizarre_workaround2->SetBinContent(0,Scale*vertical_offset);

  TFile *waveforms = new TFile(rootfilename, "update");
  bizarre_workaround->Write();
  bizarre_workaround2->Write();
  hist->Write();
  delete bizarre_workaround;
  delete bizarre_workaround2;
  delete hist;
  waveforms->Close();
}

void SimMesh::plotfinish() //Shifts, adds noise to, and discretizes the waveform.
{
  int plotmax = scope->getRecordPreSize()+scope->getRecordLength();
  
  TH1D* histprefinal = new TH1D("hist0","",16384,0,16384);
  TH1D* histfinal = new TH1D("hist1","",16384,0,16384);

  for (int j = 1; j < plotmax; j++)
    {
      histprefinal->SetBinContent(j,scope->wf(0)->GetBinContent(j) + Scale*vertical_offset);
      scope->wf(0)->SetBinContent(j,scope->wf(0)->GetBinContent(j)+Scale*(vertical_offset+noise(SD)));
      scope->wf(0)->SetBinContent(j,round(scope->wf(0)->GetBinContent(j),scoperes));
      histfinal->SetBinContent(j,scope->wf(0)->GetBinContent(j));
    }
  TFile* waveforms = new TFile(rootfilename, "update");
  histfinal->Write();
  histprefinal->Write();
  delete histfinal;
  delete histprefinal;
  waveforms->Close();
}
