#include "SimMesh.hh"
#include "TObject.h"
#include "SimScope.hh"
#include "EventGenerator.hh"


ClassImp(SimMesh)

  SimMesh::SimMesh() {/*empty constructor*/};

SimMesh::SimMesh(int plot_offset_EG, int imin_EG, double Initial_Energy_EG, double StepSize_EG, SimScope* scope_EG, double fMinEnergy_EG, vector<double>& fRecoilZ_EG, vector<double>& fRecoilEn_EG, SimChamber* fChamber_EG)
{//Constructor
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

}

SimMesh::~SimMesh(){/*destructor*/}


void SimMesh::simulate()// simulates waveforms on an initialized instance of SimMesh
{
  for(int i = imin; En>fMinEnergy && i < imax-1 && fRecoilZ.at(i)>0 && fRecoilZ.at(i)<fChamber->getDriftLength(); i++)
    //while the particle still has energy, the index i is in range, and the particle is in the chamber
    {
      final_position = i+1;
      SimMesh::step(i);
      if(need_space == 1) SimMesh::allocate();
    }
  SimMesh::decay();
  //SimMesh::noisify();
  SimMesh::plot();
}

void SimMesh::step(int i) //runs one step of electron ionization and drift
{
  double sz; //z-position                                    
  double jD; //expected arrival time of electrons   
  int jI; //jD rounded to int                                    
  double tSD; //uncertainty of arrival time             
  double El; //energy loss at current stop                                
  need_space = 0; //records whether vectors are almost out of allocated space  
  sz = fRecoilZ.at(i);
  El = fRecoilEn.at(i);
  t += dx/vel();
  tSD = timenoise(sz);
  jD = j(sz);
  
  for (int samples = 0; samples < Nsamples; samples++){
    //randomly drops electrons onto the mesh at times around jI
    
    jI = int(round(jD + noise(tSD),1));
    if (jI < 1) jI = 1; //prevents underflow
    jmin = ((jI < jmin) ? jI : jmin); //keeps track of first signal for plotting purposes
    if (jI > (jmax - jmax_init)) need_space = 1; //Q,R are almost out of space
    
    Q.at(jI)+=El/Nsamples; //record electrons' signal
  }
  En -= fRecoilEn.at(i);
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

  if (plot_offset == -1) //Resets waveform to 0 for first particle in a batch
    {
      for (int j = 0; j < plotmax; j++)
        scope->wf(0)->SetBinContent(j,0.);
    }

  if (plot_offset == -1) //-1 is used as an input when calling the function to tell the program that this is the first particle in an event. plot_offset is calculated automatically in this case. For additional particles, call the SimMesh constructor with the plot_offset value from the first particle.
    {
      plot_offset = jmin - plotmin;
      if (plot_offset == -1) plot_offset = 0;
    }
  for (int j = 1; j < plotmax; j++)
    {  
      if (j + plot_offset > 0 && j + plot_offset < R.size()) 
	scope->wf(0)->SetBinContent(j,scope->wf(0)->GetBinContent(j)+Scale*R.at(j + plot_offset));
    }
}

void SimMesh::plotfinish() //Shifts, adds noise to, and discretizes the voltage.
{
  int plotmax = scope->getRecordPreSize()+scope->getRecordLength();

  for (int j = 1; j < plotmax; j++)
    {
      scope->wf(0)->SetBinContent(j,scope->wf(0)->GetBinContent(j)+Scale*(vertical_offset+noise(SD)));
      scope->wf(0)->SetBinContent(j,round(scope->wf(0)->GetBinContent(j),scoperes));
    }
}
