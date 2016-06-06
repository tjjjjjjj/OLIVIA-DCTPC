#include "DmtpcProjection.hh"
#include "TH2.h" 
#include "TH1.h"
#include "TF1.h"
#include "MaxCamClusterImage.hh"
#include "Math/GSLMCIntegrator.h"
#include "TPaveText.h"
#include "TCanvas.h" 
#include "DmtpcGainMap.hh"
#include "TGraphErrors.h"
#include "DmtpcSkimDataset.hh"
#include "Math/Minimizer.h" 
#include "Math/Functor.h" 
#include "Math/Factory.h" 
#include "DmtpcSkimEvent.hh"
#include "MaxCamSRIM.hh"
#include "TMinuit.h"
#include <algorithm>


using namespace std; 

ClassImp(DmtpcProjection); 

struct srimLineFitIntegralParams
{
  double E; 
  double sigma; 
  const DmtpcProjection * proj; 
  enum which {ALL, POS, NEG} which; 
}; 

double SRIMLineFitFnIntegral(double *x, size_t ndim, void * p)
{
  srimLineFitIntegralParams * params = (srimLineFitIntegralParams*) p; 
  assert(ndim=3); 

  double x0 = x[0]; 
  double y0 = x[1]; 
  double y1 = x[2]; 
  if (params->which == srimLineFitIntegralParams::NEG && y1 > y0) return 0; 
  else if (params->which == srimLineFitIntegralParams::POS && y1 < y0) return 0; 
  double srimLineFitParams[5] = { params->E,x0,y0,y1,params->sigma}; 
  return  exp(-params->proj->SRIMLineFitFn(srimLineFitParams) / 2); 
}


DmtpcProjection::~DmtpcProjection()
{
  if (_long_profile) _long_profile->Delete(); 
  if (_long_profile_interpolated) _long_profile_interpolated->Delete(); 
  if (_tranv_profile) _tranv_profile->Delete(); 
  if (_EProfile) _EProfile->Delete(); 
  if (_EFitProfile) _EFitProfile->Delete(); 
}

DmtpcProjection::DmtpcProjection()
{
  _long_profile = 0; 
  _long_profile_interpolated = 0; 
  _tranv_profile = 0; 
  _EProfile = 0; 
  _EFitProfile = 0; 
}

void EfitFunc(int & npar, double * /*grad*/, double &ans, double * pars, int /*flag*/)
{
  double likely = 0; 
  DmtpcProjection * obj = (DmtpcProjection*) gMinuit->GetObjectFit(); 

  const TH1 * p = obj->getLongitudinalProfile(); 
  const TH1 * ip = obj->getInterpolatedLongitudinalProfile(); 

  int njumps = obj->getNJumps(); 
  
  TGraphErrors * g; 

  g = obj->getEFitProfile();

  double offset = pars[npar-njumps-1]; 

  double * jumpE = &(pars[npar - njumps]); 
  double runningE = offset; 
  double runningE_err = fabs(offset); 

  vector<double> E; 
  vector<double> dE; 
  vector<double> Eerr; 
  vector<double> dEerr; 
    
  int ndof = 0; 

  int ibin = 1; 
  int ijump = 0; 

//    std::cout << std::endl <<"Start jumping! " << std::endl;
    while (ijump <= njumps)
    {

//      std::cout <<" ijump: " << ijump <<  " of " << njumps << std::endl; 
      if (ijump!=njumps)
      {
//          std::cout << "Expected minbin: " << obj->getJumpLowBin(ijump) << std::endl; 
      }
     
      while (ijump < njumps ? (ibin < obj->getJumpLowBin(ijump))
                            : ibin <= ip->GetNbinsX())
      {
      
        if (ibin < obj->getNStartIgnore()) 
        {
           ibin++; 
//           std::cout << "start ignore" << std::endl; 
           continue; 
        }

        if (ibin > p->GetNbinsX()-obj->getNEndIgnore())
        {
//          std::cout << "end ignore" << std::endl; 
          break;
        }

//        std::cout << " ibin: " << ibin << std::endl; 
        double val = p->GetBinContent(ibin); 
        double err = p->GetBinError(ibin); 
        dE.push_back(val);
        dEerr.push_back(err); 

        runningE += val; 
        runningE_err += err*err; 
        E.push_back(runningE);
        Eerr.push_back(runningE_err);
        ibin++; 
      }

      if (ijump == njumps) break; 

      runningE += jumpE[ijump]; 
      runningE_err += jumpE[ijump]; 
//      std::cout << "Jumping! With value " << jumpE[ijump] << " to "; 
      ibin = obj->getJumpHighBin(ijump++) + 1; 
 //     std:: cout << ibin << std::endl; 
//
    }
    

    g->Clear(); 
    g->Set(E.size()); 

    for (int i = 0; i < int(E.size()); i++)
    {
      g->SetPoint(i,E[i],dE[i]);
      g->SetPointError(i,sqrt(Eerr[i]),dEerr[i]);
    }
    while(g->GetN() > int( E.size()))
    {
      g->RemovePoint(E.size()); 
    }


  TF1 * func = obj->getUserFunc(); 
  if (g->GetN() > 1) func->SetRange(0,g->GetX()[g->GetN()-1]); 
  func->SetParameters(pars); 

  for (int i = 0; i < g->GetN(); i++)
  {
     double y = g->GetY()[i]; 
     double dy = g->GetEY()[i]; 
     double x = g->GetX()[i]; 
     double dx = g->GetEX()[i]; 
     double Y = func->Eval(x); 
     double dY = func->Derivative(x); 
     double error2 = dy*dy + dx*dx*dY*dY;
     likely += pow(Y-y,2)/error2; 
     ndof++; 
     // likely += log(TMath::Landau(Y,y,sqrt(error2),true); 
  }

  ans =  likely;
  func->SetNDF(ndof - npar); 
}

void fitfunc(int & npar, double * /*g*/, double &f, double* pars, int /*flag*/)
{
  double likely = 0; 
  DmtpcProjection * obj = (DmtpcProjection*) gMinuit->GetObjectFit(); 

  int ndof = 0; 


  int jumpnum = 0; 
  int max = obj->getLongitudinalProfile()->GetNbinsX() - obj->getNEndIgnore(); 
//  std::cout << max << std::endl; 
  for (int i = 1; i <= max; i++)
  {
     
      double x = obj->getLongitudinalProfile()->GetBinCenter(i); 

      if (obj->getNJumps() > jumpnum && i >= obj->getJumpLowBin(jumpnum))
      {
        i = obj->getJumpHighBin(jumpnum++); 
        continue; 
      }

      TF1 * userfunc = obj->getUserFunc();

      userfunc->SetParameters(pars); 

      double lambda = userfunc->Eval(x);
      double error = obj->getLongitudinalProfile()->GetBinError(i); 
      double k = obj->getLongitudinalProfile()->GetBinContent(i); 

//      std::cout << error << std::endl; 
      //chi^2
      likely += pow(lambda - k,2)/(error*error); 
      ndof++; 
  }

  obj->getUserFunc()->SetNDF(ndof-npar); 

  f = likely;
}

TF1 * DmtpcProjection::EFit(TF1 * function, TF1 * rangefn ,
                            double jump_limits, double jump_error,  vector<string> *  minuit_cmds ) 
{



  _userfunc = function; 
  TMinuit * m = new TMinuit; 
  //m->SetPrintLevel(-1); 
  gMinuit->SetObjectFit(this); 
  m->SetFCN(&EfitFunc); 
  int err; 
  double l_lim,u_lim; 
  int parcount = 0;

  for (int i = 0; i < _userfunc->GetNpar(); i++)
  {
    _userfunc->GetParLimits(i,l_lim,u_lim); 
    m->mnparm(parcount++,_userfunc->GetParName(i), _userfunc->GetParameter(i), _userfunc->GetParError(i), l_lim, u_lim, err); 
  }

  if (rangefn == 0)
  {
    m->mnparm(parcount++,"offset",0, fabs(_long_profile->Integral()/10.),0,0,err); 
  }
  else
  {
    m->mnparm(parcount++,"offset",rangefn->Integral(-1e99,_long_profile->GetBinLowEdge(1),(const double *) NULL,0.001), fabs(_long_profile->Integral()/10.),0,0,err); 
  }

  for (int i = 0; i < _njumps; i++)
  {
    int minbin = getJumpLowBin(i);
    int maxbin = getJumpHighBin(i);
    double jump;
    if (rangefn== 0)
    {
      jump = fabs(getInterpolatedLongitudinalProfile()->Integral(minbin,maxbin)); 
    }
    else
    {
      double minval = _long_profile->GetBinLowEdge(minbin); 
      double maxval = _long_profile->GetBinLowEdge(maxbin+1); 
      jump = rangefn->Integral(minval,maxval,(const double *) NULL,0.001); 
    }
      
    if (jump > _long_profile->Integral()) jump = _long_profile->Integral()/10; 
    char buf[24]; 
    sprintf(buf,"jump%d",i); 
    m->mnparm(parcount++,buf, fabs(jump), fabs(jump*jump_error), fabs((1-jump_limits)*jump), fabs(jump*(1+jump_limits)), err); 
  }

  if (!minuit_cmds)
  {
//    m->mncomd("SEEK 100",err); 
    m->mncomd("MINIMIZE 500",err); 
  }
  else
  {
    for (unsigned i = 0; i < minuit_cmds->size(); i++)
    {
      m->mncomd(minuit_cmds->at(i).c_str(),err); 
    }
  }

  for (int i = 0; i < _userfunc->GetNpar(); i++)
  {
    double val, err; 
    m->GetParameter(i,val,err); 
    _userfunc->SetParameter(i, val); 
    _userfunc->SetParError(i, err); 
  }

  double val,edm,erddef;
  int npars,maxpar,status; 

  m->mnstat(val,edm,erddef,npars,maxpar,status); 

  _userfunc->SetChisquare(val); 
  m->Delete(); 
  return _userfunc; 
}

TF1 * DmtpcProjection::longitudinalFit(TF1 * function, vector<string> * minuit_cmds)       
{


  _userfunc = function; 
  TMinuit * m = new TMinuit; 
  //m->SetPrintLevel(1); 
  gMinuit->SetObjectFit(this); 
  m->SetFCN(&fitfunc); 
  int err; 
  double l_lim, u_lim; 

  for (int i = 0; i < _userfunc->GetNpar(); i++)
  {
    _userfunc->GetParLimits(i,l_lim,u_lim); 
    m->mnparm(i,_userfunc->GetParName(i), _userfunc->GetParameter(i), _userfunc->GetParError(i), l_lim, u_lim, err); 
    if (_userfunc->GetParError(i) == 0) 
    {
      m->FixParameter(i); 
    }
  }

  if (!minuit_cmds)
  {
    m->mncomd("MINIMIZE 500",err); 
  }
  else
  {
    for (unsigned i = 0; i < minuit_cmds->size(); i++)
    {
      m->mncomd(minuit_cmds->at(i).c_str(),err); 
    }
  }

  for (int i = 0; i < _userfunc->GetNpar(); i++)
  {
    double val, err; 
    m->GetParameter(i,val,err); 
    _userfunc->SetParameter(i, val); 
    _userfunc->SetParError(i, err); 
  }

  double val,edm,erddef;
  int npars,maxpar,status; 

  m->mnstat(val,edm,erddef,npars,maxpar,status); 

  _userfunc->SetChisquare(val); 
  m->Delete(); 
  return _userfunc; 
}


int DmtpcProjection::getJumpLowBin(int i) const 
{ 
  int sign = 1;//_flip_spacers ? -1 : 1; 
  return _long_profile->FindBin(_long_spacers[i] -  sign*getLongitudinalSpacerWidth(i));
}

int DmtpcProjection::getJumpHighBin(int i) const 
{ 
  int sign =1;// _flip_spacers ? -1 : 1; 
  return _long_profile->FindBin(_long_spacers[i]  + sign * getLongitudinalSpacerWidth(i));
}

DmtpcProjection::DmtpcProjection(DmtpcSkimEvent * e, int c, int track, const DmtpcGainMap * gm, bool flip, double sigma, int ninterp)
{
  init(e->cluster(c), track, e->phi(c,track), e->x(c,track), e->y(c,track), gm, flip, sigma, ninterp); 
}

DmtpcProjection::DmtpcProjection(const MaxCamClusterImage * clust, int track, double phi, double x, double y,  const DmtpcGainMap * gm, bool flip, double sigma, int ninterp)
{ 
  init(clust, track, phi, x, y,  gm, flip, sigma, ninterp); 
}

//Constructor
void DmtpcProjection::init(const MaxCamClusterImage * clust, int track, double phi, double x, double y,  const DmtpcGainMap * gm, bool flip, double sigma, int /*ninterp*/)
{

    _start_cutoff = false; 
    _end_cutoff = false; 

   _phi = phi ;
   double phi_t = _phi + M_PI/2;
   _theta_param = -1; 

   if (flip)
   {
     _phi += M_PI; 
     phi_t += M_PI; 
   }

   _spacer_sigma = sigma; 

  _long_profile = clust->projectClusterInterpolate(track, _phi,"bicubic",gm, 0.2, false); 
  _tranv_profile = clust->projectClusterInterpolate(track, phi_t,"bicubic",gm, 0.2, false); 
  
  _width = 4*_tranv_profile->GetRMS(); 
  //double bin_width = _tranv_profile->GetBinWidth(1); 

  //Figure out where we intersect the spacers
 

  double m = tan(_phi); 
  double b = -m * x + y; 

  double long_min, long_max, tran_min, tran_max; 

  clust->getMinMaxPosition(track,_phi, long_min,long_max); 
  clust->getMinMaxPosition(track,phi_t , tran_min,tran_max); 
  

  
  _flip_spacers = cos(_phi) < 0; 
  for (int ii = 0; gm && ii < gm->getNSpacers(); ii++)
  {

    int i = _flip_spacers ? gm->getNSpacers() - ii - 1: ii; 

    double mm = gm->getSpacerSlope(i); 
    double bb = gm->getSpacerIntercept(i); 
    
    if (m == mm) continue; 
    double xint = (bb - b) / (m-mm); 
    double yint = m * xint + b; 

    double plong = xint * cos(_phi) + yint* sin(_phi); 
    if (plong < long_max && plong > long_min)
    {
      double cos_ang =fabs(1+m*mm) / sqrt((1+m*m)*(1+mm*mm)); 
      double sin_ang = fabs(m - mm) / sqrt((1+m*m)*(1+mm*mm)); 
      double tan_ang = sin_ang / cos_ang; 
      _long_spacers.push_back(plong); 
      _long_spacer_widths.push_back(sigma * gm->getSpacerWidth(i) / sin_ang + _width / tan_ang); 
    }

    double ptran = xint * cos(phi_t) + yint* sin(phi_t); 
    if (ptran < tran_max && ptran > tran_min)
    {
//      double sin_ang = fabs(m - mm) / sqrt((1+m*m)*(1+mm*mm)); 
      double cos_ang =(1+m*mm) / sqrt((1+m*m)*(1+mm*mm)); 
      _tranv_spacers.push_back(ptran); 
      _tranv_spacer_widths.push_back(gm->getSpacerWidth(i) / cos_ang); 
    }
  }

  _njumps = _long_spacers.size(); 
  int max_interpolated = _long_profile->GetNbinsX(); 

  if (_njumps > 0 && getJumpHighBin(_njumps-1) > _long_profile->GetNbinsX())
  {
    max_interpolated = getJumpLowBin(--_njumps); 
  }


  //_long_profile_interpolated = (TH1*) _long_profile->Clone("long_prof_interpolated"); 
  char buf[128]; 
  sprintf(buf,"interpolated_%s",_long_profile->GetName()); 
  _long_profile_interpolated = new TH1F(buf,buf, max_interpolated,_long_profile->GetBinLowEdge(1), _long_profile->GetBinLowEdge(max_interpolated+1)); 

  for (int i = 1; i <= max_interpolated; i++)
  {
    _long_profile_interpolated->SetBinContent(i, _long_profile->GetBinContent(i));  
    _long_profile_interpolated->SetBinError(i, _long_profile->GetBinError(i));  
  }


/*** 

  //Interpolate with pol2 over spacers
  for (int i = 0; i < _njumps; i++)
  {

     if (max_interpolated < 2 * ninterp)
     {
       break; 
     }
     vector<double> x; 
     vector<double> dE; 
     vector<double> dE_error; 
    
     int minbin = getJumpLowBin(i); 

     for (int bin = minbin-(ninterp+1); bin < minbin; bin++)
     {
      if (bin > 0)
      {
        x.push_back(_long_profile->GetBinCenter(bin)); 
        dE.push_back(_long_profile->GetBinContent(bin)); 
        dE_error.push_back(_long_profile->GetBinError(bin)); 
      }
     }

     int maxbin = getJumpHighBin(i); 
     for (int bin = maxbin+1; bin < maxbin+(ninterp+1); bin++)
     {
      if (bin <= _long_profile->GetNbinsX())
      {
        x.push_back(_long_profile->GetBinCenter(bin)); 
        dE.push_back(_long_profile->GetBinContent(bin)); 
        dE_error.push_back(_long_profile->GetBinError(bin)); 
      }
     }
   
     TF1 f("f","pol2",x[0],x[x.size()-1]);     
     TGraphErrors g(x.size(), &x[0], &dE[0],0,&dE_error[0]); 
     g.Fit(&f,""); 


     for (int bin = minbin; bin <= maxbin; bin++)
     {
       double newval = f.Eval(_long_profile->GetBinCenter(bin));
//       cout << newval << endl; 
       _long_profile_interpolated->SetBinContent(bin, newval); 
       _long_profile_interpolated->SetBinError(bin, _long_profile->GetBinError(bin) * 2); 
     }

  }

  ***/
  vector<double> E, dE, Eerr, dEerr; 

  double runningE = 0; 
  double runningdE = 0; 

  for (int i = 1; i <= _long_profile_interpolated->GetNbinsX(); i++)
  {
    double val = _long_profile_interpolated->GetBinContent(i);
    dE.push_back(val); 
    double err = _long_profile_interpolated->GetBinError(i);
    
    dEerr.push_back(err); 
    runningE += val; 
    runningdE = sqrt(runningdE * runningdE + err*err);
    E.push_back(runningE); 
    Eerr.push_back(runningdE); 
  }

  _EProfile = new TGraphErrors(E.size(), &(E[0]), &(dE[0]), &(Eerr[0]), &(dEerr[0])); 
  _EFitProfile = (TGraphErrors*)_EProfile->Clone(); 
}

static TCanvas * line_fit_canvas = 0; 
static TCanvas * srim_fit_canvas = 0; 

DmtpcProjection::SRIMFitParams * DmtpcProjection::doSRIMFit(DmtpcSRIMProjection::PARTICLE type, double lengthcal, double gain, double pressure, const char * opt) const
{

  TH1* tran = (TH1*)  _tranv_profile->Clone("tran"); 
  TF1 * gaus = new TF1("fgaus","gaus(0)", tran->GetXaxis()->GetXmin(), tran->GetXaxis()->GetXmax()); 
  gaus->SetParameters(tran->GetMaximum(), tran->GetBinCenter(tran->GetMaximumBin()), tran->GetRMS()); 
  tran->Fit(gaus,"IG"); 

  double sigma = gaus->GetParameter(2); 
  double E  = gaus->GetParameter(0) * sqrt(2*M_PI) * sigma; 
  MaxCamSRIM srim(DmtpcSRIMProjection::getFileName(type)); 
  srim.setPressure(pressure); 
  _srim = &srim; 

  double binning = _long_profile->GetBinWidth(1); 
  double xmin = _long_profile->GetXaxis()->GetXmin(); 
  double xmax = _long_profile->GetXaxis()->GetXmax(); 
  double offset =  (xmax - xmin) > 1.5 * sigma ? 0.5 * sigma : binning/2; 
  double range = xmax - xmin  - 2 * offset;

  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad"); 
  min->SetMaxFunctionCalls(10000); 
  min->SetPrintLevel(100); 
  

  ROOT::Math::Functor f(this, &DmtpcProjection::SRIMfitFn, 7); 
  min->SetFunction(f); 

  min->SetFixedVariable(0,"gain", gain); 
  min->SetFixedVariable(1, "lengthcal",lengthcal); 
  min->SetFixedVariable(2,"E",E); 
  min->SetLimitedVariable(3,"x0",xmin+offset, 0.5 * sigma, xmin - sigma, xmax + sigma); 
  min->SetLimitedVariable(4,"range",range, 0.5 * sigma, 0, xmax - xmin); 
  min->SetFixedVariable(5,"sigma",sigma); 
  min->SetFixedVariable(6,"pressure",pressure); 
  
  SRIMFitParams * ret = new SRIMFitParams; 


  ret->success = false; 
  if (xmax!= xmin && sigma !=0)
  {
    ret->success = min->Minimize(); 
  }


  const double * answer = min->X(); 
  ret->chisq = min->MinValue(); 
  ret->gain = answer[0]; 
  ret->sin_theta = answer[1]; 
  ret->lengthcal = answer[2]; 
  ret->E = answer[3]; 
  ret->x0 = answer[4]; 
  ret->range = answer[5]; 
  ret->sigma = answer[6]; 
  ret->type = type; 
  ret->pressure = answer[7]; 
  ret->ndof = _long_profile->GetNbinsX() - min->NFree() + 2; 

  ret->phi = _phi; 

  if (strchr(opt,'D')) 
  {
    if (!srim_fit_canvas) srim_fit_canvas = new  TCanvas("srim_c","srim_c",1600,800);  
    srim_fit_canvas->Clear(); 
    srim_fit_canvas->Divide(2,1); 
    srim_fit_canvas->cd(1); 

    DmtpcSRIMProjection p(_srim, ret->sigma, ret->pressure, 0, ret->range  * ret->lengthcal / ret->sin_theta); 
    p.setGain(ret->gain); 
    p.setLengthCal(ret->lengthcal); 
    p.setSinTheta(ret->sin_theta); 
    p.setOffset(ret->x0); 

    double * pdfx = new double[_long_profile->GetNbinsX()]; 
    double * pdfy = new double[_long_profile->GetNbinsX()]; 

    for (int ii = 1; ii <= _long_profile->GetNbinsX(); ii++)
    {
      std::cout << ii << std::endl; 
      pdfx[ii-1] = _long_profile->GetBinCenter(ii);  
      pdfy[ii-1] = p.E(_long_profile->GetBinLowEdge(ii), _long_profile->GetBinLowEdge(ii+1)); 
    }
      
    _long_profile->DrawCopy(); 
    TGraph  * g = new TGraph(_long_profile->GetNbinsX(), pdfx, pdfy); 
    g->SetMarkerColor(2); 
    g->SetMarkerStyle(4); 
    g->Draw("psame"); 

    //delete pdfx; 
    //delete pdfy; 

    srim_fit_canvas->cd(2); 
    tran->SetTitle("Transverse Fit"); 
    tran->DrawCopy(); 
    gaus->DrawCopy("lpsame"); 
    srim_fit_canvas->Update(); 
  }


  delete tran; 
  delete gaus; 
  delete min; 

  return ret; 
}


DmtpcProjection::SRIMLineFitParams * 
DmtpcProjection::doSRIMLineFit(MaxCamSRIM * srim, double lengthcal, double gain, const char * opt) const
{
  TH1* tran = (TH1*)  _tranv_profile->Clone("tran"); 
  TF1 * gaus = new TF1("fgaus","gaus(0)", tran->GetXaxis()->GetXmin(), tran->GetXaxis()->GetXmax()); 
  gaus->SetParameters(tran->GetMaximum(), tran->GetBinCenter(tran->GetMaximumBin()), tran->GetRMS()); 
  tran->Fit(gaus,TString("I") + (strchr(opt,'V') ? TString("V") : TString("Q")) ); 

  double sigma = gaus->GetParameter(2); 
  double E  = gaus->GetParameter(0) * sqrt(2*M_PI) * sigma; 
  double binning = _long_profile->GetBinWidth(1); 
  double Ephys = E/gain / binning; 

  double xmin = _long_profile->GetXaxis()->GetXmin(); 
  double xmax = _long_profile->GetXaxis()->GetXmax(); 
  double offset =  (xmax - xmin) > 1.5 * sigma ?  sigma : binning/2; 
  double range_3d = srim->getRangeVsEnergy(false)->Eval(Ephys) / lengthcal + binning;
  double range_guess = (xmax - xmin - 2 * offset); 
  if (range_guess < binning) range_guess = binning;

  double srimheight = srim->getEnergyVsRange(false)->Eval(lengthcal * binning) * gain; 
  double maxheight = E/binning*2; 
  double startheight = E/range_guess; 
  double minheight = srim->getEnergyVsRange(false)->Eval(lengthcal * binning/2) * gain; 
  
  if (startheight < minheight)
  {
    startheight = srimheight;  
  }

  if (strchr(opt,'V'))
  {
    std::cout << "range_3d: " << range_3d << std::endl; 
    std::cout << "range_guess: " << range_guess << std::endl; 
    std::cout << "minheight: " << minheight << std::endl; 
    std::cout << "srimheight: " << srimheight << std::endl; 
    std::cout << "maxheight: " << maxheight << std::endl; 
    std::cout << "startheight: " << startheight << std::endl; 
    std::cout << "offset: " << offset << std::endl; 
    std::cout << "Ephys: " << Ephys << std::endl; 
    std::cout << "E: " << E << std::endl; 
    std::cout << "xmin: " << xmin << std::endl; 
    std::cout << "xmax: " << xmax << std::endl; 
    std::cout << "sigma: " << sigma << std::endl; 
  }


  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad"); 
  min->SetMaxFunctionCalls(10000); 
  

  ROOT::Math::Functor f(this, &DmtpcProjection::SRIMLineFitFn, 5); 
  min->SetFunction(f); 

  min->SetFixedVariable(0,"E", E); 
  min->SetLimitedVariable(1,"x0",xmin+offset, 0.5 * sigma, xmin - sigma, xmax + sigma); 
  min->SetLimitedVariable(2,"y0",startheight, minheight/4, minheight, maxheight); 
  min->SetLimitedVariable(3,"y1",startheight, minheight/4, minheight, maxheight); 
  min->SetFixedVariable(4,"sigma",sigma); 
  
  SRIMLineFitParams * ret = new SRIMLineFitParams; 

  ret->success = min->Minimize(); 

  if ((strchr(opt,'V')))
  {
    min->PrintResults(); 
  }

  const double * answer = min->X(); 

  double x0 = answer[1]; 
  double y0 = answer[2]; 
  double y1 = answer[3]; 

  min->GetMinosError(1, ret->x0ErrorLow,ret->x0ErrorHigh); 
  min->GetMinosError(2, ret->y0ErrorLow,ret->y0ErrorHigh); 
  min->GetMinosError(3, ret->y1ErrorLow,ret->y1ErrorHigh); 
  ret->x0 = x0; 
  ret->x0Error = min->Errors()[1]; 
  ret->y0 = y0;
  ret->y0Error = min->Errors()[2]; 
  ret->y1 = y1;
  ret->y1Error = min->Errors()[3]; 

  ret->chisq = min->MinValue(); 
  ret->gain = gain; 
  ret->sinThetaMin = srimheight/ TMath::Min(y0,y1);

  ret->lengthcal = lengthcal; 
  ret->E = E/binning; 
  ret->EPhys = Ephys; 
  ret->EError = gaus->GetParError(0) / binning; 
  ret->EPhysError = gaus->GetParError(0) / binning / gain; 
  ret->range = 2 * E / (y0+y1);  
  ret->rangeError = ret->range *  2 * sqrt( pow(ret->EError,2) / pow(ret->E,2) + (pow(ret->y0Error,2) + pow(ret->y1Error,2)) / pow(y0 + y1,2)); 
  ret->sinThetaSlope = ret->range / range_3d; 
  ret->sigma = sigma; 
  ret->sigmaError = gaus->GetParError(2); 
  ret->ndof = _long_profile->GetNbinsX() - min->NFree() + 2; 
  ret->transverseChisq = gaus->GetChisquare(); 
  ret->transverseNdof = gaus->GetNDF(); 

  double slope = (y1*y1 - y0*y0)/(2*E); 
  double x1 = x0 + ret->range; 
  double b =  y0 - x0 * slope; 
  ret->phi = _phi;   
  if (y0<y1)  ret->phi = DmtpcMath::normalizeAngle(ret->phi + M_PI) ; 

  if (strchr(opt,'D')) 
  {
    if (!line_fit_canvas) line_fit_canvas = new  TCanvas("line_c","line_c",1800,1000);  
    line_fit_canvas->Clear(); 
    line_fit_canvas->Divide(3,2); 
    line_fit_canvas->cd(1); 

    double max_y = TMath::Max(y0, y1); 
    max_y = TMath::Max(max_y, _long_profile->GetMaximum()); 
    max_y+=10; 
    max_y *= 1.1; 

    double min_y = TMath::Min(y0,y1); 
    min_y = TMath::Min(min_y, _long_profile->GetMinimum()); 
    min_y-=10; 
    min_y *= 1.1; 


    TH2C axis_maker("axis_maker","Longitudinal Fit", _long_profile->GetNbinsX(),xmin,xmax, 10, min_y, max_y); 
    axis_maker.SetStats(0); 
    axis_maker.DrawCopy(); 
    _long_profile->DrawCopy("lpsame");    
    TF1 func("func", &DmtpcMath::lineSegmentConvolvedWithGaussian, xmin, xmax, 5); 
    double line_params[]  = {slope,b,x0,x1,sigma}; 
    func.SetParameters(line_params); 
    func.DrawCopy("lpsame"); 
    TF1 func2("func2","pol1", x0,x1); 
    func2.SetParameters(b,slope); 
    func2.SetLineColor(3); 
    func2.DrawCopy("lpsame"); 

    double * pdfx = new double[_long_profile->GetNbinsX()]; 
    double * pdfy = new double[_long_profile->GetNbinsX()]; 

    for (int ii = 1; ii <= _long_profile->GetNbinsX(); ii++)
    {
      pdfx[ii-1] = _long_profile->GetBinCenter(ii);  
      pdfy[ii-1] = DmtpcMath::integralOfLineSegmentConvolvedWithGaussian(_long_profile->GetXaxis()->GetBinLowEdge(ii), _long_profile->GetXaxis()->GetBinLowEdge(ii+1), line_params)/binning; 

      //double x = pdfx[ii-1]; 
      //cout << DmtpcMath::lineSegmentConvolvedWithGaussian(&x, line_params) << endl;; 
    }
      
    TGraph  * g = new TGraph(_long_profile->GetNbinsX(), pdfx, pdfy); 
    g->SetMarkerColor(4); 
    g->Draw("psame"); 

    //delete pdfx; 
    //delete pdfy; 

    line_fit_canvas->cd(2); 
    
    TPaveText * pt = new TPaveText(0.05,0.2,0.95,0.95); 
    pt->AddText( TString::Format("#chi^2/ndof : %f/%d  = %f", ret->chisq, ret->ndof, ret->chisq/ret->ndof)); 
    pt->AddText( TString::Format("x0: %f +/-%f (%f,%f)", ret->x0, ret->x0Error, ret->x0ErrorLow, ret->x0ErrorHigh)); 
    pt->AddText( TString::Format("y0: %f +/-%f (%f,%f)", ret->y0, ret->y0Error, ret->y0ErrorLow, ret->y0ErrorHigh)); 
    pt->AddText( TString::Format("y1: %f +/-%f (%f,%f)", ret->y1, ret->y1Error, ret->y1ErrorLow, ret->y1ErrorHigh)); 
    pt->AddText( TString::Format("E: %f +/- %f adu (%f +/- %f keVee)", ret->E, ret->EError, ret->EPhys, ret->EPhysError)); 
    pt->AddText( TString::Format("Range: %f +/- %f pixels (%f +/- %f mm)", ret->range, ret->rangeError, ret->range * lengthcal, ret->rangeError * lengthcal)); 
    pt->AddText( TString::Format("sin #theta (from range): %f  (from min): %f", ret->sinThetaSlope, ret->sinThetaMin)); 

    pt->Draw(); 
    line_fit_canvas->cd(3); 
    tran->SetTitle("Transverse Fit"); 
    tran->DrawCopy(); 
    gaus->DrawCopy("lpsame"); 


    double sigmas[] = {0.5,1,1.5,2,2.5,3}; 
    int nsigmas = sizeof(sigmas) / sizeof(double); 
    unsigned npts = 50; 

    line_fit_canvas->cd(4); 
    double error_def = min->ErrorDef(); 

    TH2C cont0axis("cont0axis","cont0axis", 10, minheight, ret->y0 + 5 * ret->y0Error , 10, minheight, ret->y1 + 5 * ret->y1Error); 
    cont0axis.SetTitle("y0 vs y1"); 
    cont0axis.SetStats(0); 
    cont0axis.GetXaxis()->SetTitle("y0"); 
    cont0axis.GetYaxis()->SetTitle("y1"); 
    cont0axis.DrawCopy(); 

    for (int sig = 0; sig < nsigmas; sig++)
    {
      double cx[npts+1]; 
      double cy[npts+1]; 
      min->SetErrorDef(sigmas[nsigmas - sig - 1]); 
      min->Contour(2,3,npts,cx,cy); 
      cx[npts] = cx[0]; 
      cy[npts] = cy[0]; 
      TGraph * cf = new TGraph(npts+1,cx,cy); 
      cf->SetLineColor(nsigmas - sig); 
      cf->Draw("csame"); 
 
    }

    TF1 div("div","pol1", ret->y0 - 10 * ret->y0Error, ret->y0 + 10 * ret->y0Error); 
    div.SetParameters(0,1); 
    div.SetLineStyle(2); 
    div.SetLineColor(1); 
    div.DrawCopy("lsame");

    line_fit_canvas->cd(5); 

    TH2C cont1axis("cont1axis","cont1axis", 10, ret->x0 - 5*ret->x0Error, ret->x0 + 5 * ret->x0Error, 10, minheight, ret->y1 + 5 * ret->y1Error); 
    cont1axis.SetTitle("x0 vs y1"); 
    cont1axis.SetStats(0); 
    cont1axis.GetXaxis()->SetTitle("x0"); 
    cont1axis.GetYaxis()->SetTitle("y1"); 
    cont1axis.DrawCopy(); 

    for (int sig = 0; sig < nsigmas; sig++)
    {
      double cx[npts+1]; 
      double cy[npts+1]; 
      min->SetErrorDef(sigmas[nsigmas  - sig- 1]); 
      min->Contour(1,3,npts,cx,cy); 
      cx[npts] = cx[0]; 
      cy[npts] = cy[0]; 
      TGraph * cf = new TGraph(npts+1,cx,cy); 
      cf->SetLineColor(nsigmas - sig); 
      cf->Draw("csame"); 
      
    }

    line_fit_canvas->cd(6); 

    TH2C cont2axis("cont2axis","cont2axis", 10, ret->x0 - 5*ret->x0Error, ret->x0 + 5 * ret->x0Error, 10, minheight, ret->y0 + 5 * ret->y0Error); 
    cont2axis.SetTitle("x0 vs y0"); 
    cont2axis.SetStats(0); 
    cont2axis.GetXaxis()->SetTitle("x0"); 
    cont2axis.GetYaxis()->SetTitle("y0"); 
    cont2axis.DrawCopy(); 

    for (int sig = 0; sig < nsigmas; sig++)
    {
      double cx[npts+1]; 
      double cy[npts+1]; 
      min->SetErrorDef(sigmas[nsigmas-sig - 1]); 
      min->Contour(1,2,npts,cx,cy); 
      cx[npts] = cx[0]; 
      cy[npts] = cy[0]; 
      TGraph * cf = new TGraph(npts+1,cx,cy); 
      cf->SetLineColor(nsigmas-sig); 
      cf->Draw("csame"); 
    }

    min->SetErrorDef(error_def); 
  }

  if (strchr(opt,'V'))
  {
    std::cout <<  "Chisq/ndof: " << ret->chisq / ret->ndof << std::endl <<
                  "x0: " << ret->x0 << " +/-" << ret->x0Error << std::endl <<
                  "y0: " << ret->y0 << " +/-" << ret->y0Error << std::endl <<
                  "y1: " << ret->y1 << " +/-" << ret->y1Error << std::endl <<
                  "m: " << slope << std::endl <<
                  "range: " << ret->range << std::endl <<
                  "E: " << ret->E << std::endl <<
                  "phi " << ret->phi << std::endl <<
                  "sin_theta_min: " << ret->sinThetaMin << std::endl <<
                  "sin_theta_slope: " << ret->sinThetaSlope << std::endl; 
  }


  //do integrals! 
  
  if (strchr(opt,'I'))
  {
    ROOT::Math::GSLMCIntegrator gator(ROOT::Math::MCIntegration::kVEGAS,0,0,5000);  

    double gator_maxes[] = { 
          TMath::Min(xmax+sigma, ret->x0 + 4 * (ret->x0ErrorHigh > 0 ? ret->x0ErrorHigh : ret->x0Error)) , 
          TMath::Min(maxheight, ret->y0 + 4 * ( ret->y0ErrorHigh > 0 ? ret->y0ErrorHigh : ret->y0Error)),
          TMath::Min(maxheight, ret->y1 + 4 * ( ret->y1ErrorHigh > 0 ? ret->y1ErrorHigh : ret->y1Error))
          } ; 
    double gator_mins[] = { 
          TMath::Max(xmin-sigma, ret->x0 + 4 * (ret->x0ErrorLow < 0 ? ret->x0ErrorLow : -ret->x0Error)), 
          TMath::Max(minheight, ret->y0 + 4 * ( ret->y0ErrorLow < 0 ? ret->y0ErrorLow : -ret->y0Error)),
          TMath::Max(minheight, ret->y1 + 4 * ( ret->y1ErrorLow < 0 ? ret->y1ErrorLow : -ret->y1Error))
          } ; 

    srimLineFitIntegralParams gator_pars; 
    gator_pars.E = E; 
    gator_pars.sigma = ret->sigma; 
    gator_pars.proj = this; 
    gator_pars.which = srimLineFitIntegralParams::ALL; 

    //doing total likely integral 
    ret->likelyIntegralAll = gator.Integral(&SRIMLineFitFnIntegral,3,gator_mins, gator_maxes, (void *) &gator_pars); 

    if (strchr(opt,'V')) std::cout << " total integral: " << ret->likelyIntegralAll << std::endl; 
    gator_pars.which = ret->y0 >= ret->y1 ? srimLineFitIntegralParams::NEG : srimLineFitIntegralParams::POS; 
    ret->likelyIntegralSame = gator.Integral(&SRIMLineFitFnIntegral,3,gator_mins, gator_maxes, (void *) &gator_pars); 

    if (strchr(opt,'V')) std::cout << " same side integral: " << ret->likelyIntegralSame << std::endl; 

    ret->prob = ret->likelyIntegralSame / ret->likelyIntegralAll; 

    if (strchr(opt,'V')) std::cout << " implied probability: " << ret->prob << std::endl; 
  }

  else 
  {
    ret->likelyIntegralAll = 0; 
    ret->likelyIntegralSame = 0; 
    ret->prob = 0; 
  }

  delete tran; 
  delete gaus; 
  delete min; 

  return ret; 
}
double DmtpcProjection::SRIMfitFn(const double *xx) const
{
    
  double gain = xx[0]; 
  double lengthcal = xx[1];
  double E = xx[2]; 
  double offset = xx[3];  
  double range = xx[4]; 
  double sigma = xx[5]; 
  double pressure = xx[6]; 

  double sin_theta = range / _srim->getRangeVsEnergy(false)->Eval(E); 

  DmtpcSRIMProjection p(_srim, sigma, pressure, 0, range  * lengthcal / sin_theta); 
  double binning = _long_profile->GetBinWidth(1); 
  p.setGain(gain); 
  p.setLengthCal(lengthcal); 
  p.setSinTheta(sin_theta); 
  p.setOffset(offset); 

  double chisq = 0; 
  for (int i = 1; i <= _long_profile->GetNbinsX(); i++)
  {
    double val_hist = _long_profile->GetBinContent(i); 
    double err_hist = _long_profile->GetBinError(i); 
    double val_start = _long_profile->GetBinLowEdge(i); 
    double val_end = _long_profile->GetBinLowEdge(i+1); 
    double val_fn = p.E(val_start, val_end);
    chisq += pow((val_hist-val_fn) / err_hist,2); 
  }


  double xmin = _long_profile->GetXaxis()->GetXmin(); 
  double xmax = _long_profile->GetXaxis()->GetXmax(); 

  //we want the function to be close to 0 outside the interval... 
  for (int i = 1; i <= 10; i++)
  {
    double eval_low  = (xmin - binning*i - offset) * lengthcal * sin_theta;  
    double eval_high  = (xmax + binning*i - offset) * lengthcal * sin_theta;  
    double val_low = p.dE(eval_low) * gain * binning  * lengthcal / sin_theta;  
    double val_high = p.dE(eval_high) * gain * binning  * lengthcal / sin_theta;  

    chisq += pow(val_low/_long_profile->GetBinError(1),2); 
    chisq += pow(val_high/_long_profile->GetBinError(_long_profile->GetNbinsX()),2); 
  }

  return chisq; 

}



DmtpcProjection::LineFitParams * DmtpcProjection::doLineFit(const char * opt) const
{

  TH1* tran = (TH1*)  _tranv_profile->Clone("tran"); 
  TF1 * gaus = new TF1("fgaus","gaus(0)", tran->GetXaxis()->GetXmin(), tran->GetXaxis()->GetXmax()); 
  gaus->SetParameters(tran->GetMaximum(), tran->GetBinCenter(tran->GetMaximumBin()), tran->GetRMS()); 
  tran->Fit(gaus,"IG"); 

  LineFitParams * ret = new LineFitParams(); 

  double sigma = gaus->GetParameter(2); 
  double E  = gaus->GetParameter(0) * sqrt(2*M_PI) * sigma; 

  ret->sigma = sigma; 
  ret->sigma_error = gaus->GetParErrors()[2]; 
  ret->transverse_chisq = gaus->GetChisquare(); 
  ret->transverse_ndof = gaus->GetNDF(); 

  std::cout << "Sigma: " << sigma << std::endl; 

  //fit left half and right half of same sigma gaussian onto longitudinal projection
  //
  double xmin = _long_profile->GetXaxis()->GetXmin(); 
  double xmax = _long_profile->GetXaxis()->GetXmax(); 

  double binning = _long_profile->GetBinWidth(1); 

  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad"); 

  min->SetMaxFunctionCalls(100000); 
  min->SetPrintLevel(100); 

  ROOT::Math::Functor f(this, &DmtpcProjection::lineFitFn, 5); 
  min->SetFunction(f); 

  double offset =  (xmax - xmin) > 1.5 * sigma ? 0.5 * sigma : binning/2; 
  min->SetLimitedVariable(0, "phi", 0, 0.25, -M_PI/2., M_PI/2.); 
//  min->SetLimitedVariable(1, "E",  E,0.01*E ,0.75*E, 1.5*E); 
  min->SetFixedVariable(1, "E",  E); 
  min->SetLimitedVariable(2, "x0", xmin + offset, 0.5 * sigma, xmin - sigma, xmax + sigma); 
  min->SetLimitedVariable(3, "range", xmax - 2*offset - xmin, 0.5 * sigma, 0, xmax-xmin); 
  min->SetFixedVariable(4, "sigma", sigma); 

  ret->success = (xmin!= xmax && sigma > 0 ) ? min->Minimize(): 0; 

  const double * answer = min->X(); 
  ret->alpha_error = min->Errors()[0]; 
  ret->alpha=answer[0]; 

  ret->chisq = min->MinValue(); 
  ret->slope = tan(answer[0]); 
  ret->offset= answer[2]; 
  ret->offset_error = min->Errors()[2]; 
  ret->ndof = _long_profile->GetNbinsX() - min->NFree() + 2; 

//  std::cout << f

  ret->phi = _phi; 
  if (ret->slope > 0)  ret->phi = DmtpcMath::normalizeAngle(ret->phi + M_PI) ; 
  ret->range = answer[3]; 
  ret->range_error = min->Errors()[3]; 

  double x0 = answer[2]; 
  double x1 = answer[3] + answer[2]; 
  ret->E = answer[1]; 
  double b = (ret->E - 0.5 * ret->slope *(pow(x1,2) - pow(x0,2))) / (x1-x0); 
  ret->b = b; 
  if (strchr(opt,'D')) 
  {
    if (!line_fit_canvas) line_fit_canvas = new  TCanvas("line_c","line_c",1600,800);  
    line_fit_canvas->Clear(); 
    line_fit_canvas->Divide(2,1); 
    line_fit_canvas->cd(1); 

    double max_y = TMath::Max(x0*ret->slope+b, x1*ret->slope+b); 
    max_y = TMath::Max(max_y, _long_profile->GetMaximum()); 
    max_y+=10; 
    max_y *= 1.1; 

    double min_y = TMath::Min(x0*ret->slope+b, x1*ret->slope+b); 
    min_y = TMath::Min(min_y, _long_profile->GetMinimum()); 
    min_y-=10; 
    min_y *= 1.1; 


    TH2C axis_maker("axis_maker","Longitudinal Fit", _long_profile->GetNbinsX(),xmin,xmax, 10, min_y, max_y); 
    axis_maker.DrawCopy(); 
    _long_profile->DrawCopy("plsame");    
    TF1 func("func", &DmtpcMath::lineSegmentConvolvedWithGaussian, xmin, xmax, 5); 
    double line_params[]  = {ret->slope,b,x0,x1,sigma}; 
    func.SetParameters(line_params); 
    func.DrawCopy("lpsame"); 
    TF1 func2("func2","pol1", x0,x1); 
    func2.SetParameters(b,ret->slope); 
    func2.SetLineColor(3); 
    func2.DrawCopy("lpsame"); 

    double * pdfx = new double[_long_profile->GetNbinsX()]; 
    double * pdfy = new double[_long_profile->GetNbinsX()]; 

    for (int ii = 1; ii <= _long_profile->GetNbinsX(); ii++)
    {
      pdfx[ii-1] = _long_profile->GetBinCenter(ii);  
      pdfy[ii-1] = DmtpcMath::integralOfLineSegmentConvolvedWithGaussian(_long_profile->GetXaxis()->GetBinLowEdge(ii), _long_profile->GetXaxis()->GetBinLowEdge(ii+1), line_params)/binning; 

      //double x = pdfx[ii-1]; 
      //cout << DmtpcMath::lineSegmentConvolvedWithGaussian(&x, line_params) << endl;; 
    }
      
    TGraph  * g = new TGraph(_long_profile->GetNbinsX(), pdfx, pdfy); 
    g->SetMarkerColor(4); 
    g->Draw("psame"); 

    //delete pdfx; 
    //delete pdfy; 

    line_fit_canvas->cd(2); 
    tran->SetTitle("Transverse Fit"); 
    tran->DrawCopy(); 
    gaus->DrawCopy("lpsame"); 
  }


  delete tran; 
  delete gaus; 
  delete min; 

  return ret; 
}

double DmtpcProjection::SRIMLineFitFn(const double * p) const 
{

  double E = p[0]; 
  double x0 = p[1]; 
  double y0 = p[2]; 
  double y1 = p[3]; 
  double sigma = p[4]; 
  double m = (y1*y1-y0*y0)/(2*E); 

  double range =2*E/(y0+y1); 
  double x1 = x0 + range; 
  double binning = _long_profile->GetBinWidth(1); 
  double b =  y0 - m * x0; 

  double pars[] = { m, b, x0,x1,sigma}; 

  double chisq = 0; 
  for (int i = 1; i <= _long_profile->GetNbinsX(); i++)
  {
    double val_hist = _long_profile->GetBinContent(i); 
    double lim_low = _long_profile->GetBinLowEdge(i); 
    double lim_high = lim_low + binning; 
    double err_hist = _long_profile->GetBinError(i); 
    double val_fn = DmtpcMath::integralOfLineSegmentConvolvedWithGaussian(lim_low,lim_high, pars)  / binning; 
    chisq += pow((val_hist-val_fn) / err_hist,2); 
  }

  //we want the function to be close to 0 outside the interval... 
  
  double xmin = _long_profile->GetXaxis()->GetXmin(); 
  double xmax = _long_profile->GetXaxis()->GetXmax(); 
  chisq += pow(DmtpcMath::integralOfLineSegmentConvolvedWithGaussian(xmin - 10 * sigma, xmin,pars)/ (_long_profile->GetBinError(1)  * (10*sigma/binning)),2); 
  chisq += pow(DmtpcMath::integralOfLineSegmentConvolvedWithGaussian(xmax, xmax +  10*sigma,pars) / (_long_profile->GetBinError(_long_profile->GetNbinsX()) * (10*sigma/binning)),2); 

  return chisq;  
}




double DmtpcProjection::lineFitFn(const double * p) const 
{

  double binning = _long_profile->GetBinWidth(1); 
  double m = tan(p[0]); 
  double E = p[1]; 
  double start = p[2]; 
  double range = p[3]; 
  double end = start + range; 
  double sigma = p[4]; 
  double b =  (E- 0.5 * m * (end*end - start*start)) / (end - start); 

  double pars[] = { m, b, start,end,sigma}; 

  double chisq = 0; 
  for (int i = 1; i <= _long_profile->GetNbinsX(); i++)
  {
    double val_hist = _long_profile->GetBinContent(i); 
    double lim_low = _long_profile->GetBinLowEdge(i); 
    double lim_high = lim_low + binning; 
    double err_hist = _long_profile->GetBinError(i); 
//    double val_center = _long_profile->GetBinCenter(i); 
//    double val_fn = DmtpcMath::lineSegmentConvolvedWithGaussian(&val_center,pars) ; 
    double val_fn = DmtpcMath::integralOfLineSegmentConvolvedWithGaussian(lim_low,lim_high, pars)  / binning; 
    chisq += pow((val_hist-val_fn) / err_hist,2); 
  }


 //We don't want the line to cross the x axis! 
  double xint = -b/m;   
  if (xint > start && xint < end) 
  {
    if (m > 0)
    {
      chisq -= 10 * (xint-start) * (m * start + b); 
    }
    else
    {
      chisq -= 10 * (end-xint) * (m * end + b); 
    }
  }

  //we want the function to be close to 0 outside the interval... 
  
  double xmin = _long_profile->GetXaxis()->GetXmin(); 
  double xmax = _long_profile->GetXaxis()->GetXmax(); 
  chisq += pow(DmtpcMath::integralOfLineSegmentConvolvedWithGaussian(xmin - 3* sigma, xmin,pars),2); /// _long_profile->GetBinError(1),2); 
  chisq += pow(DmtpcMath::integralOfLineSegmentConvolvedWithGaussian(xmax, xmax + 3*sigma,pars),2);// / _long_profile->GetBinError(_long_profile->GetNbinsX()),2); 

  return chisq;  
}


