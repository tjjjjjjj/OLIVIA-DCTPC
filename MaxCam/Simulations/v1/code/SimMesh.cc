#include "SimMesh.hh"
#include "TObject.h"
#include "SimScope.hh"
#include "EventGenerator.hh"


ClassImp(SimMesh)

  SimMesh::SimMesh() {/*empty constructor*/};
SimMesh::SimMesh(double Initial_Energy_EG, double StepSize_EG, SimScope* scope_EG, int nscope_EG, double fMinEnergy_EG, vector<double>& fRecoilZ_EG, vector<double>& fRecoilEn_EG, SimChamber* fChamber_EG)
{
  //Constructor
  mass = 3.727e6 / (9e22); // alpha in keV/c^2, c in mm/s
  fRecoilEn = fRecoilEn_EG;
  fRecoilZ = fRecoilZ_EG;
  fMinEnergy = fMinEnergy_EG;
  nscope = 0; //I think this is never used
  scope = scope_EG;
  fChamber = fChamber_EG;
  Td = scope->getDecayTime();
  Tr = scope->getRiseTime();
  vd = scope->getDriftVel();
  Scale = scope->getScale();
  scoperes = scope->getScopeRes();
  SD = scope->getVoltageNoise() * scoperes/Scale;
  rate = scope->getClockRate();
  dt = 1./rate;
  a1 = scope->getTimeNoise1();
  a2 = scope->getTimeNoise2();
  zero_offset = scope->getZeroOffset();
  zero_offset_var = scope->getZeroOffsetVar();
  zero_offset += int(pow(fRecoilZ.at(0),0.5)*zero_offset_var);
  Nsamples = scope->getNSamples();
  vertical_offset = scope->getVerticalOffset() * 1/Scale;
  jmax_init = 20000; //Magic number since value isn't very important. Could theoretically cause out-of-range errors if it is too low.
  dx = StepSize_EG;
  Initial_Energy = Initial_Energy_EG;
  En = Initial_Energy;
  rand = new TRandom3(0);
  imin = 0;
  imax = imin + fRecoilZ.size();
  j_init = int(j(fRecoilZ.at(0)));
  jmax = jmax_init + j_init;
  Q = vector<double>(jmax,0);
  R = vector<double>(jmax,0);
    t = 0;
}
    SimMesh::~SimMesh(){}

    void SimMesh::simulate()
    {
      for(int i = imin; En>fMinEnergy && i < imax-1 && fRecoilZ.at(i)>0 && fRecoilZ.at(i)<fChamber->getDriftLength(); i++)
	{
	  SimMesh::step(i);
	  if(need_space == 1) SimMesh::allocate();
	}
      SimMesh::decay();
      SimMesh::noisify();
      SimMesh::plot();
    }
  void SimMesh::step(int i)
  {
    double sz; //z-position                                    
    double jD; //expected arrival time   
    int jI; //rounded to int                                    
    double tSD; //uncertainty of arrival time             
    double El; //energy loss                                
    need_space = 0; //vectors are almost out of allocated space  
    sz = fRecoilZ.at(i);
    El = fRecoilEn.at(i);
    En -= fRecoilEn.at(i);
    t += dx/vel();
    tSD = timenoise(sz);
    jD = j(sz);
    
    for (int samples = 0; samples < Nsamples; samples++)//randomly drops electrons onto the mesh at times around jI                                     
      {
	jI = int(round(jD + noise(tSD),1));
	if (jI < 1) jI = 1;
	jmin = (jI < jmin) ? jI : jmin;
	if (jI > (jmax - jmax_init)) need_space = 1;

	Q.at(jI)+=El/Nsamples;
      }
    En -= fRecoilEn.at(i);
  }

double SimMesh::round(double x, double res)
{double offset;
  if (x <0.)
    offset = -.5;
  if (x==0.)
    offset =  0.;
  if (x >0.)
    offset =  .5;
  return res*(double(int(x/res + offset)));}



double SimMesh::vecsum(vector<double>& B, int risesteps)
{
  double sum = 0;
  for(int r=0;r<risesteps;r++)
    sum+=B.at(r);
  return sum;
}


  void SimMesh::allocate()
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

  void SimMesh::decay()
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

  void SimMesh::noisify()
  {
    for(int i=0;i<jmax;i++)
	R.at(i) = R.at(i) + noise(SD);
  }

  void SimMesh::plot()
  {
    int presize = scope->getRecordPreSize();
    int plotmax = presize + scope->getRecordLength();
    int plotmin = presize + zero_offset;

    //Discretizes the waveform to the scope's resolution, and centers the waveform spike near t=0                                                     
    for (int j = plotmin; j < plotmax; j++)
      { scope->wf(0)->SetBinContent(j,round(Scale*(vertical_offset+R.at(j - (plotmin - jmin))),scoperes));
	//cout << endl << round(Scale*(vertical_offset + R.at(j+jmin-plotmin)),scoperes);
	//cout << endl << '#' << Scale* (vertical_offset+R.at(j+jmin-plotmin));
	/*if (scope->wf(0)->GetBinContent(j)!=0)
	  cout << endl << scope->wf(0)->GetBinContent(j);
	if (j%500==0)
	cout << "500 more...";*/
      }
      for (int j = 1; j < plotmin; j++)
	scope->wf(0)->SetBinContent(j,round(Scale*(vertical_offset+noise(SD)),scoperes));

 
  }
