#include "recoilEnergy.hh"
#include "TH2F.h"
#include "TFile.h"
#include "TRint.h"
#include "TMath.h"

static bool gainMapInit = false; 
static TH2F * gainMap[2]; 
static std::string gainMapDir = "." ;
static TFile * f[2]; 

static void loadGainMaps()
{
  f[0] = TFile::Open((gainMapDir+"/gainMap0.root").c_str());  
  f[1] = TFile::Open((gainMapDir+"/gainMap1.root").c_str());  
  gainMap[0] = (TH2F*) f[0]->Get("image"); 
  gainMap[1] = (TH2F*) f[1]->Get("image"); 
  gainMapInit = true; 
}

double RecoilEnergy::visibleToRecoil(const double E)
{
  if (E>200) return E/0.76879; 
  return 4.46388 + 1.65231*E - 0.00286719*E*E + E*E*E*(5.8344e-6) - E*E*E*E*(4.22845e-9);
}

void RecoilEnergy::setGainMapDirectory(std::string dir)
{
  gainMapDir = dir;
  if (gainMapInit)
  {
    gainMapInit==false; 
    f[0]->Close(); 
    f[1]->Close();   
    delete f[0];
    delete f[1];
  }
}


double RecoilEnergy::getRecoilEnergy(DmtpcSkimEvent * ev, 
                                     int c, int t, bool gain)
{
  if (gain)
    return getRecoilEnergy(ev->E(c,t),c,ev->x(c,t), ev->y(c,t)); 
  else 
    return getRecoilEnergy(ev->E(c,t),c); 
}

double RecoilEnergy::getVisibleEnergy(DmtpcSkimEvent * ev, 
                                      int c, int t, bool gain)
{

  double E_cal = c == 0 ? 9.5 : 12.9; 

  if (gain) 
  {
    if (!gainMapInit) loadGainMaps(); 
    double cal = gainMap[c]->GetBinContent(TMath::Nint(ev->x(c,t)/8), TMath::Nint(ev->y(c,t)/8)); 
    if (cal!=0) E_cal = cal; 
  } 

  return ev->E(c,t)/E_cal; 

}

double RecoilEnergy::getRecoilEnergy(double E, int cam, double x, double y)
{
  if (!gainMapInit) loadGainMaps(); 

  double gain = gainMap[cam]->GetBinContent(TMath::Nint(x/8),
                                            TMath::Nint(y/8)); 
  if (gain ==0) return getRecoilEnergy(E,cam); 
  double vis_kev = E/gain; 

  return visibleToRecoil(vis_kev); 
}

double RecoilEnergy::getRecoilEnergy(double E, int cam)
{
  double E_cal = cam == 0 ? 9.5 : 12.9; 

  double vis_kev = E/E_cal; 
  return visibleToRecoil(vis_kev); 
}
