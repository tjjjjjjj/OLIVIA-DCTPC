#include "SimMesh.hh"
#include "TObject.h"
#include "SimScope.hh"
#include "EventGenerator.hh"


ClassImp(SimMesh)

  SimMesh::SimMesh(double Initial_Energy_EG, double StepSize_EG, SimScope* scope_EG, int nscope_EG, EventGenerator* fEventGen_EG)
{
  //Constructor
  mass = 3.727e6 / (9 * 10^22); // alpha in keV/c^2, c in mm/s
  fEventGen = fEventGen_EG;
  nscope = 0; //magic number for now, since we only have one
  scope = scope_EG;
  Td = scope(0)->getDecayTime();
  Tr = scope(0)->getRiseTime();
  vd = scope(0)->getDriftVel();
  Scale = scope(0)->getScale();
  scoperes = scope(0)->getScopeRes();
  SD = scope(0)->getVoltageNoise() * scoperes/Scale;
  rate = scope(0)->getClockRate();
  dt = 1./rate;
  a1 = scope(0)->getTimeNoise1();
  a2 = scope(0)->getTimeNoise2();
  zero_offset = scope(0)->getZeroOffset();
  Nsamples = scope(0)->getNSamples();
  vertical_offset = scope(0)->getVerticalOffset() * 1/Scale;
  jmax_init = 20000; //Magic number, not in .tmp files since value isn't very important. Could theoretically cause out-of-range errors if it is too low.
  Q = new vector<double>(jmax_init,0);
  R = new vector<double>(jmax_init,0);
  dx = StepSize_EG;
  Initial_Energy = Initial_Energy_EG;
  rand = new TRandom3(0);
  imin = 0;
  imax = imin + fEventGen->fRecoilZ.size();
  j_init = int(j(0,fEventGen->fRecoilZ.at(0));
  t = 0;
}
    SimMesh::~SimMesh()
    {
    }


    void SimMesh::simulate()
    {
      for(int i = imin; fEventGen->fRecoilEn.at(i)>fEventGen->fMinEnergy && i < imax-1 && fEventGen->fRecoilZ.at(i)>0 && fEventGen->fRecoilZ.at(i)<fChamber->getDriftLength(); i++)
	{
	  SimMesh::step();
	  if(need_space == 1) SimMesh::allocate();
	}
      SimMesh::decay();
      SimMish::noisify();
      SimMesh::plot();
    }
  void SimMesh::step()
  {
    double sz; //z-position                                    
    double jD; //expected arrival time   
    int jI; //rounded to int                                    
    double tSD; //uncertainty of arrival time             
    double El; //energy loss                                
    need_space = 0; //vectors are almost out of allocated space  
    sz = fEventGen->fRecoilZ.at(i);
    El = fEventGen->fRecoilEn.at(i);
    t += dx/v(En);
    tSD = timenoise(sz);
    jD = j(sz);

    for (int samples = 0; samples < Nsamples; samples++)//randomly drops electrons onto the mesh at times around jI                                     
      {
	jI = int(jD + 0.5 + noise(tSD));
	if (jI < 1) jI = 1;
	jmin = (jI < jmin) ? jI : jmin;
	if (jI > (jmax - jmax_init)) need_space = 1;

	Q.at(jI)+=El/Nsamples;
      }
    En -= fRecoilEn.at(i);
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
	R.at(i) = R.at(i) + noise(rand, SD);
  }

  void SimMesh::plot()
  {
    int plotmax = scope(0)->getRecordLength()+scope(0)->getRecordPreSize();
    int presize = scope(0)->getRecordPreSize();
    int plotmin = presize + zero_offset;

    //Discretizes the waveform to the scope's resolution, and centers the waveform spike near t=0                                                     
  for (int j = plotmin; j < plotmax; j++)
    scope(0)->wf(0)->SetBinContent(j,round(Scale*(vertical_offset+S.at(j - (plotmin - jmin))),scoperes));
 for (int j = 1; j < plotmin; j++)
   scope(0)->wf(0)->SetBinContent(j,round(Scale*(vertical_offset+rand->Gaus(0,SD)),scoperes));

 
  }
