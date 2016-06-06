#include "DmtpcTrackFitter.hh" 

#include "DmtpcMath.hh"
#include "DmtpcRootTools.hh" 
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include "MaxCamImageTools.hh"
#include "TPaveText.h"
#include "Math/Factory.h"
#include "Math/GSLMCIntegrator.h"
#include "TRandom3.h"
#include "limits.h"

static const size_t NPAR = 7; 
static const int PAR_E = 0; 
static const int PAR_X0 = 1; 
static const int PAR_Y0 = 2; 
static const int PAR_PHI = 3; 
static const int PAR_Z0 = 4; 
static const int PAR_DELTA_Z = 5; 
static const int PAR_SIGMA = 6; 

ClassImp(DmtpcTrackFitter); 
ClassImp(DmtpcTrackFitter::Result); 
ClassImp(DmtpcTrackFitter::Param); 

struct srimLineFitIntegralParams
{
   const DmtpcTrackFitter * fitter; 
   enum which {ALL, SAME,OPPOSITE} which; 
   double phi; 
   double x0; 
   double y0; 
   double sigma; 
   double E; 
   bool flip; 
}; 

static double SRIMLineFitFnIntegral(double *x, size_t ndim, void * p)
{
  srimLineFitIntegralParams * params = (srimLineFitIntegralParams*) p; 
  double phi = params->phi; 
  double z1 = x[2]; 
  double z0 = x[1]; 
  double dz = z1-z0; 
  double offset = x[0]; 


  double fpar[NPAR]; 
  fpar[PAR_E] = params->E;
  fpar[PAR_X0] = params->x0 + offset * cos(phi); 
  fpar[PAR_Y0] = params->y0 + offset * sin(phi); 
  fpar[PAR_PHI] = params->phi; 
  fpar[PAR_Z0] = z0; 
  fpar[PAR_DELTA_Z] = dz; 
  fpar[PAR_SIGMA] = params->sigma; 


  bool same_side = (params->flip == (dz > 0));
  if (params->which == srimLineFitIntegralParams::SAME && !same_side) return 0; 
  else if (params->which == srimLineFitIntegralParams::OPPOSITE && same_side) return 0; 
  return exp(-params->fitter->SRIMLineFit2DFn(fpar)/2); 
}



void DmtpcTrackFitter::fillResultParams()
{

  const double * vals = min->X(); 
  const double * errors = min->Errors(); 

  //set values

  result->fitE.val = vals[PAR_E]; 
  result->E.val = 0; 
  result->E.err = 0; 

  result->x0.val = vals[PAR_X0]; 
  result->y0.val = vals[PAR_Y0]; 
  result->phi.val = vals[PAR_PHI]; 
  result->z0.val = vals[PAR_Z0]; 
  result->delta_z.val = vals[PAR_DELTA_Z]; 
  result->z1.val= vals[PAR_Z0] + vals[PAR_DELTA_Z]; 
  result->sigma.val = vals[PAR_SIGMA]; 

  result->fitE.err = errors[PAR_E]; 
  result->x0.err = errors[PAR_X0]; 
  result->y0.err = errors[PAR_Y0]; 
  result->phi.err = errors[PAR_PHI]; 
  result->z0.err = errors[PAR_Z0]; 
  result->delta_z.err = errors[PAR_DELTA_Z]; 
  result->z1.err= sqrt(pow(errors[PAR_Z0],2) + pow(errors[PAR_DELTA_Z],2)); 
  result->sigma.err = errors[PAR_SIGMA]; 

  result->fitE.zeroErr(); 
  result->x0.zeroErr(); 
  result->y0.zeroErr(); 
  result->phi.zeroErr(); 
  result->z0.zeroErr(); 
  result->z1.zeroErr(); 
  result->z1.zeroErr(); 
  result->sigma.zeroErr(); 

  //do minos if always minos or z0, z1 close enough to limits 
  result->minos = always_minos 
    || fabs(result->z0.val-minheight)/result->z0.err < minos_thresh
    || fabs(result->z1.val-minheight)/result->z1.err < minos_thresh
    || fabs(result->z0.val-maxheight)/result->z0.err < minos_thresh
    || fabs(result->z1.val-maxheight)/result->z1.err < minos_thresh;

  if (result->minos)
  {
    min->GetMinosError(PAR_E,result->fitE.errDn, result->fitE.errUp); 
  //  min->GetMinosError(PAR_X0,result->x0.errDn, result->x0.errUp); 
  //  min->GetMinosError(PAR_Y0,result->y0.errDn, result->y0.errUp); 
    min->GetMinosError(PAR_Z0,result->z0.errDn, result->z0.errUp); 
    min->GetMinosError(PAR_DELTA_Z,result->delta_z.errDn, result->delta_z.errUp); 
    result->z1.errDn = -DmtpcMath::sumError(result->delta_z.errDn, result->z0.errDn); 
    result->z1.errUp = -DmtpcMath::sumError(result->delta_z.errUp, result->z0.errUp); 
 //   min->GetMinosError(PAR_PHI,result->phi.errDn, result->phi.errUp); 
 //   min->GetMinosError(PAR_SIGMA,result->sigma.errDn, result->sigma.errUp); 
  }


  result->fitE.fixErr(); 
  result->x0.fixErr(); 
  result->y0.fixErr(); 
  result->phi.fixErr(); 
  result->z0.fixErr(); 
  result->z1.fixErr(); 
  result->delta_z.fixErr(); 
  result->sigma.fixErr(); 

  result->range.val = 2 * binwidth * result->fitE.val / (2*result->z0.val + result->delta_z.val); 
  //incorporate covariance of z0 and z1... but not with E since E error is probably relatively small 
  result->range.err = result->range.val *  2 * sqrt( pow(result->fitE.err,2) / pow(result->fitE.val,2) + (2*pow(result->z0.err,2) + pow(result->delta_z.err,2) + 2*min->CovMatrix(PAR_Z0,PAR_DELTA_Z)) / pow(2*result->z0.val + result->delta_z.val,2)); 

  //TODO: figure out a better way to set these
  result->range.errDn = result->range.err; 
  result->range.errUp = result->range.err; 



  result->chisq = min->MinValue(); 
  result->ndof = fithist->GetNbinsX()*fithist->GetNbinsY() - min->NFree(); 
  result->gain = gain; 


  result->flip = result->delta_z.val > 0; 

  result->htphi.err = result->phi.err; 
  result->htphi.errUp = result->phi.errUp; 
  result->htphi.errDn = result->phi.errDn; 
  result->htphi.val = result->flip ? DmtpcMath::normalizeAngle(result->phi.val + M_PI) : result->phi.val; 



}


DmtpcTrackFitter::~DmtpcTrackFitter()
{
  if (fitfn) delete fitfn; 
  if (wrapped_fitfn) delete wrapped_fitfn; 
  if (fit_integrator) delete fit_integrator; 
  if (result) delete result; 
  if (arrow) delete arrow; 
  if (canvas) delete canvas; 
  if (fithist) delete fithist; 
  if (local_gainmap) delete local_gainmap; 
  delete min; 
}



DmtpcTrackFitter::DmtpcTrackFitter(MaxCamSRIM * sr, double rangecal, double g, const TH2* gainmap, double min_g, double rn)
  : minfn(this, &DmtpcTrackFitter::SRIMLineFit2DFn,NPAR) 
{
  gain = g; 
  gainMap = gainmap;
  lengthcal = rangecal; 
  srim = sr; 
  result = 0; 
  min_gain = min_g; 
  srim = sr; 
  read_noise = rn; 
  min = ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad"); 
  fitfn = new TF2("trackfitfn",&DmtpcMath::lineSegmentConvolvedWithGaussian2D,0,1,0,1,NPAR); 
  wrapped_fitfn = new ROOT::Math::WrappedMultiFunction<TF2&> (*fitfn,2); 
  fit_integrator = new ROOT::Math::AdaptiveIntegratorMultiDim(*wrapped_fitfn, 1e-5,1e-4,100); 
  min->SetFunction(minfn); 
  min->SetPrintLevel(10); 
  arrow = new TArrow(); 
  zero_out = false; 
  padding = 3; 
  fithist = 0; 
  local_gainmap = 0; 
  integral_nsigma=3; 
  integral_reltol=0.05; 
  integral_abstol = 0; 
  integral_ncalls = 10000; 
  canvas = 0; 
  do_integral = false; 
  always_minos = false; 
  minos_thresh = 1; 
  quickeval = false; 
}




double DmtpcTrackFitter::SRIMLineFit2DFn(const double * p) const
{

  double E = p[PAR_E] * binwidth * binwidth; 
  double x0 = p[PAR_X0]; 
  double y0 = p[PAR_Y0]; 
  double phi = p[PAR_PHI];
  double raw_z0 = p[PAR_Z0]; 
  double raw_z1 = p[PAR_Z0] + p[PAR_DELTA_Z]; 
  double z0 = p[PAR_Z0] * binwidth; 
  double z1 = p[PAR_DELTA_Z] * binwidth + z0; 


  double chisq = 0; 

  //penalize for minheight, maxheight
  
  if (raw_z0 < minheight)
  {
    chisq += pow(raw_z0-minheight,2) / minheight; 
    if (raw_z0 < 0) 
    {
      chisq-= 10 * raw_z0;   
    }
  }

  if (raw_z1 < minheight)
  {
    chisq += pow(raw_z1-minheight,2) / minheight; 
    if (raw_z1 < 0) 
    {
      chisq -= 10 * raw_z1;   
    }
  }

  if (raw_z1 > maxheight)
  {
    chisq += pow(raw_z1-maxheight,2) / maxheight; 
  }

  if (raw_z0 > maxheight)
  {
    chisq += pow(raw_z1-maxheight,2) / maxheight; 
  }


  double sigma = p[PAR_SIGMA]; 

  double range =  2*E/(z0+z1);  
  double m =  (z1 - z0)/range; 
  double x1 = x0 + range * cos(phi); 
  double y1 = y0 + range * sin(phi); 

  double pars[] = {m,z0,x0,y0,x1,y1,sigma}; 

  if (DmtpcRootTools::checkNaN(NPAR,p)) 
  {
    std::cout << "p is NaN: ";  
    DmtpcRootTools::printVec(NPAR,p);
    return 1e99; 
  }

  if (DmtpcRootTools::checkNaN(NPAR,pars)) 
  {
    std::cout << "pars is NaN: ";  
    DmtpcRootTools::printVec(NPAR,pars);
    return 1e99 ; 
  }

  if (E<=0) 
  {
    std::cout << "E is Non-Positive: ";  
    DmtpcRootTools::printVec(NPAR,pars);
    return 1e99 ; 
  }




  fitfn->SetRange(fithist->GetXaxis()->GetXmin(), fithist->GetXaxis()->GetXmax(), fithist->GetYaxis()->GetXmin(), fithist->GetYaxis()->GetXmax()) ; 
  fitfn->SetParameters(pars); 

  double bin2 = binwidth*binwidth; 
  for (int i = 1; i <= fithist->GetNbinsX(); i++)
  {
    for (int j = 1; j <= fithist->GetNbinsY(); j++)
    {

      double valhist = fithist->GetBinContent(i,j); 
      
      double xlow = fithist->GetXaxis()->GetBinLowEdge(i); 
      double ylow = fithist->GetYaxis()->GetBinLowEdge(j); 
      double xhigh = xlow + binwidth;
      double yhigh = ylow + binwidth; 
      double err2 = read_noise * read_noise + valhist/binwidth; 
 //     double valfn = fitfn->Integral(xlow,xhigh,ylow,yhigh); 
      //double valfn = fitfn->Eval((xlow+xhigh)/2, (ylow+yhigh)/2); 
      double vlow[2] = {xlow,ylow}; 
      double vhigh[2] = {xhigh,yhigh}; 
     // double valfn = fit_integrator->Integral(vlow,vhigh)/bin2 * local_gainmap->GetBinContent(i,j); 
      double valfn = (quickeval ? fitfn->Eval((xlow+xhigh)/2,(ylow+yhigh)/2) : fit_integrator->Integral(vlow,vhigh)/bin2) * local_gainmap->GetBinContent(i,j); 
      //cout <<valfn << endl; 
      chisq += pow(valfn-valhist,2)/err2; 
    }
  }

  return chisq; 
}


int DmtpcTrackFitter::fit(DmtpcSkimEvent * ev, int c, int t) 
{
  vector<int> pix = ev->cluster(c)->getCluster(t); 
  return fit(ev->cluster(c)->getImage(), &pix, ev->phi(c,t), ev->range(c,t), sqrt(ev->transverse_moment(c,2,t)));
}

int DmtpcTrackFitter::fit(const TH2 * img, const vector<int> * pix, double phi_guess, double range_guess, double sigma_guess)
{
  binwidth = img->GetXaxis()->GetBinWidth(1); 

  //assign not ridiculous values if nothing specified 
  if (phi_guess == DBL_MAX)
  {
    static TRandom3 rand; 
    phi_guess =rand.Rndm() * 2 * TMath::Pi(); 
  }
  if (range_guess == DBL_MAX)
  {
    range_guess = binwidth * sqrt(pix->size())+1; 
  }
  if (sigma_guess == DBL_MAX)
  {
    sigma_guess = 2 * binwidth; 
  }


  double x0,y0; 

  double max = -DBL_MAX; 
  double E = 0;

  double xmin,xmax,ymin,ymax; 
  int xminbin = INT_MAX, yminbin = INT_MAX, xmaxbin = -INT_MAX,ymaxbin = -INT_MAX; 

  for (unsigned i = 0; i < pix->size(); i++)
  {
    int b = pix->at(i); 
    double val = img->GetBinContent(b);
    int bx,by,bz; 
    img->GetBinXYZ(b,bx,by,bz); 

    if (bx < xminbin)
    {
      xminbin = bx; 
      xmin = img->GetXaxis()->GetBinLowEdge(bx); 
    }
    if (bx > xmaxbin)
    {
      xmaxbin = bx; 
      xmax = img->GetXaxis()->GetBinLowEdge(bx)+binwidth; 
    }
    if (by < yminbin)
    {
      yminbin = by; 
      ymin = img->GetYaxis()->GetBinLowEdge(by); 
    }
    if (by > ymaxbin)
    {
      ymaxbin = by; 
      ymax = img->GetYaxis()->GetBinLowEdge(by)+binwidth; 
    }


    if (val > max)
    {
      max = val; 
      img->GetBinXYZ(b,bx,by,bz); 
      x0 = img->GetXaxis()->GetBinCenter(bx); 
      y0 = img->GetXaxis()->GetBinCenter(by); 
    }

    /*
    E +=  !gainMap ? val 
          : gainMap->GetBinContent(b) > min_gain ? val / gainMap->GetBinContent(b) 
          : 0; 
          */
    E+= val; 
  }

  if (xminbin - padding >= 1) 
  {
    xminbin -= padding; 
    xmin -= padding * binwidth; 
  }
  else 
  {
    xmin -= (xminbin-1)*binwidth;
    xminbin = 1; 
  }

  if (yminbin - padding >= 1) 
  {
    yminbin -= padding; 
    ymin -= padding * binwidth; 
  }
  else 
  {
    ymin -= (yminbin-1)*binwidth;
    yminbin = 1; 
  }

  if (xminbin + padding <= img->GetNbinsX()) 
  {
    xmaxbin += padding; 
    xmax += padding * binwidth; 
  }
  else 
  {
    xmax += (img->GetNbinsX() - xmaxbin)*binwidth;
    xmaxbin = img->GetNbinsX(); 
  }

  if (yminbin + padding <= img->GetNbinsY()) 
  {
    ymaxbin += padding; 
    ymax += padding * binwidth; 
  }
  else 
  {
    ymax += (img->GetNbinsY() - ymaxbin)*binwidth;
    ymaxbin = img->GetNbinsY(); 
  }

  if (fithist) delete fithist; 
  if (local_gainmap) delete local_gainmap; 

  fithist = new TH2F("fithist","fithist", xmaxbin-xminbin + 1, xmin,xmax, ymaxbin-yminbin +1,ymin,ymax); 
  local_gainmap = new TH2F("gm_local","gm_local", xmaxbin-xminbin + 1, xmin,xmax, ymaxbin-yminbin +1,ymin,ymax); 


  if (zero_out)
  {
    for (unsigned i = 0; i < pix->size(); i++)
    {
        int bin = pix->at(i); 
        int ix,iy,iz; 
        img->GetBinXYZ(bin,ix,iy,iz); 

        fithist->SetBinContent(ix-xminbin+1, iy-yminbin+1, img->GetBinContent(ix,iy)); 
        local_gainmap->SetBinContent(ix-xminbin+1, iy-yminbin+1, !gainMap ? 1 : gainMap->GetBinContent(ix,iy)); 

    }

  }
  else
  {
    for (int ix = xminbin; ix <=xmaxbin; ix++)
    {
      for (int iy = yminbin; iy <=ymaxbin; iy++)
      {
        fithist->SetBinContent(ix-xminbin+1, iy-yminbin+1, img->GetBinContent(ix,iy)); 
        local_gainmap->SetBinContent(ix-xminbin+1, iy-yminbin+1, !gainMap ? 1 : gainMap->GetBinContent(ix,iy)); 
      }
    }
  }




//  double x1 = x0 + range_guess * cos(phi_guess); 
//  double y1 = x0 + range_guess * cos(phi_guess); 

  double srimheight = srim->getEnergyVsRange(false)->Eval(lengthcal * binwidth) * gain; 
  maxheight = E/2; 
  double startheight = srim->getStoppingVsEnergy(false)->Eval(E/gain)*gain*lengthcal*binwidth; 
  double start_delta = srim->getEnergyVsRange(false)->Eval(lengthcal * binwidth) * gain - startheight; 
  if (start_delta > 0) start_delta = -startheight/2; 

  minheight = srim->getEnergyVsRange(false)->Eval(lengthcal * binwidth/2) * gain * 0.75 ; 

  if (startheight < minheight) startheight = srimheight; 


//  min->SetVariable(PAR_E,"E",E,sqrt(E)); 
  cout <<"start E: " << E << endl; 
  min->SetLimitedVariable(PAR_E,"E", E, sqrt(E),0,2*E); 
  min->SetLimitedVariable(PAR_X0,"x0",x0, 0.5 * sigma_guess, xmin - sigma_guess, xmax + sigma_guess); 
  min->SetLimitedVariable(PAR_Y0,"y0",y0, 0.5 * sigma_guess, ymin - sigma_guess, ymax + sigma_guess); 
  min->SetVariable(PAR_PHI,"phi",phi_guess,TMath::Pi()/16); 
  min->SetVariable(PAR_Z0,"z0",startheight,minheight/4); 
  min->SetVariable(PAR_DELTA_Z,"delta_z",start_delta,minheight/4); 
  min->SetLimitedVariable(PAR_SIGMA,"sigma",sigma_guess, sqrt(sigma_guess),0, binwidth * sqrt(pow(fithist->GetNbinsX(),2.) + pow(fithist->GetNbinsY(),2.)) ); 

  if (!result) result = new Result; 
  result->success = min->Minimize(); 
  min->Hesse(); 

  if (verbose)
  {
    min->PrintResults(); 
  }


  fillResultParams(); 

  double m = (result->delta_z.val*binwidth) / (result->range.val); 
  double x1 = result->x0.val + result->range.val * cos(result->phi.val); 
  double y1 = result->y0.val + result->range.val * sin(result->phi.val); 
  double pars[] = {m,result->z0.val*binwidth,result->x0.val,result->y0.val,x1,y1,result->sigma.val}; 
  fitfn->SetParameters(pars); 

  double bin2 = binwidth * binwidth; 
  for (int i = 1; i <=fithist->GetNbinsX(); i++)
  {
    for (int j = 1; j <=fithist->GetNbinsY(); j++)
    {
      double xlow = fithist->GetXaxis()->GetBinLowEdge(i); 
      double ylow = fithist->GetYaxis()->GetBinLowEdge(j); 
      double xhigh = xlow + binwidth;
      double yhigh = ylow + binwidth; 
      double vlow[2] = {xlow,ylow}; 
      double vhigh[2] = {xhigh,yhigh}; 
      double valfn = (quickeval ? fitfn->Eval((xlow+xhigh)/2,(ylow+yhigh)/2) : fit_integrator->Integral(vlow,vhigh)/bin2) * local_gainmap->GetBinContent(i,j); 
      result->E.val+= valfn; 
    }
  }

  result->Ephys.val = result->E.val / gain; 
  result->E.err = result->fitE.err * result->E.val / result->fitE.val; 
  result->Ephys.err = result->E.err/gain; 
  result->E.errDn =result->fitE.errDn * result->E.val / result->fitE.val; 
  result->E.errUp =result->fitE.errUp * result->E.val / result->fitE.val; 
  result->Ephys.errDn =result->E.errDn / gain;
  result->Ephys.errUp =result->E.errUp / gain; 



  


  //compute quick probability estimate based on overlap of z0, z1 multiplied by probability that phi is on same side
  double zscore = fabs(result->delta_z.val); 
  zscore /=  (result->delta_z.val > 0 ? -result->delta_z.errDn : result->delta_z.errUp); 
  cout <<"z: " << zscore << endl; 
  result->probFast = 0.5  *( 1 + TMath::Erf(zscore/sqrt(2.))*TMath::Erf(M_PI/result->phi.err/sqrt(2.))); ; 

  result->fullProb = do_integral; 

  if (do_integral)
  {

#ifdef HAVE_MATHMORE
    ROOT::Math::GSLMCIntegrator gator(ROOT::Math::MCIntegration::kVEGAS,integral_abstol,integral_reltol,integral_ncalls);  

    // parameters for integral: 

    // xoffset 
    // z0 
    // z1
     
    double gator_maxes[3] = 
    {
      result->sigma.val, 
      maxheight, 
      maxheight

    }; 

    double gator_mins[3] = 
    {
      -result->sigma.val, 
      minheight,
      minheight
    }; 



    srimLineFitIntegralParams gator_pars; 
    gator_pars.fitter = this; 
    gator_pars.phi = result->phi.val;
    gator_pars.flip = result->flip; 
    gator_pars.x0 = result->x0.val; 
    gator_pars.y0 = result->y0.val; 
    gator_pars.sigma = result->sigma.val; 
    gator_pars.E = result->E.val; 


    if (verbose)
    {
      for (size_t par = 0; par < 3; par++)
      {
        std::cout << "par " << par <<" integration limits:  (" << gator_mins[par] << " , " << gator_maxes[par] << ")"  <<std::endl; 
      }
    }

    //do  total likely integral 
    bool q = quickeval; 
    quickeval = true; 
    gator_pars.which = srimLineFitIntegralParams::SAME;  
    result->likelyIntegralSame = gator.Integral(&SRIMLineFitFnIntegral,3,gator_mins, gator_maxes, (void *) &gator_pars); 
    result->likelyIntegralSameError = gator.Sigma(); 
    gator_pars.which = srimLineFitIntegralParams::ALL; 
    result->likelyIntegralTotal = gator.Integral(&SRIMLineFitFnIntegral,3,gator_mins, gator_maxes, (void *) &gator_pars); 
    result->likelyIntegralTotalError = gator.Sigma(); 
    result->prob = result->likelyIntegralSame / result->likelyIntegralTotal; 
    result->probError = DmtpcMath::quotientError(result->likelyIntegralSame,result->likelyIntegralTotal,result->likelyIntegralSameError, result->likelyIntegralTotalError); 
    quickeval = q; 

#else
    std::cerr <<"Integration not supported without HAVE_MATHMORE" <<std::endl; 
    assert(0); 
#endif
  }
  else
  {
    result->likelyIntegralTotal = 0; 
    result->likelyIntegralTotalError = 0; 
    result->likelyIntegralSame = 0; 
    result->likelyIntegralSameError = 0; 
    result->prob = -1; 
    result->probError = 0; 
  }


  if (draw)
  {

    fitfn->SetParameters(pars); 

    gStyle->SetOptStat(0); 
    if (!canvas) canvas = new TCanvas("trackfitter_c","Track Fitter Canvas", 1800,1000); 
    canvas->Clear(); 
    canvas->Divide(3,2); 

    canvas->cd(1); 
    fithist->SetTitle("Data Hist"); 
    fithist->Draw("colz"); 

    canvas->cd(2); 
    TH2 * fnhist = (TH2*) fithist->Clone(); 
    fnhist->SetTitle("Fit Function"); 



    TH2 * residualhist = (TH2*) fithist->Clone(); 

    double bin2 = binwidth * binwidth; 
    for (int i = 1; i <=fithist->GetNbinsX(); i++)
    {
      for (int j = 1; j <=fithist->GetNbinsY(); j++)
      {
        double xlow = fithist->GetXaxis()->GetBinLowEdge(i); 
        double ylow = fithist->GetYaxis()->GetBinLowEdge(j); 
        double xhigh = xlow + binwidth;
        double yhigh = ylow + binwidth; 
//        fnhist->SetBinContent(i,j,fitfn->Integral(xlow,xhigh,ylow,yhigh)); 
//        fnhist->SetBinContent(i,j,fitfn->Eval((xlow+xhigh)/2,(ylow+yhigh)/2)); 
        double vlow[2] = {xlow,ylow}; 
        double vhigh[2] = {xhigh,yhigh}; 
        double valfn = (quickeval ? fitfn->Eval((xlow+xhigh)/2, (ylow+yhigh)/2) : fit_integrator->Integral(vlow,vhigh)/bin2) * local_gainmap->GetBinContent(i,j); 
        fnhist->SetBinContent(i,j,valfn); 
        double err2 = valfn + read_noise*read_noise; 
        residualhist->SetBinContent(i,j,pow(valfn-fithist->GetBinContent(i,j),2)/err2); 
      }
    }

    new (arrow) TArrow( result->x0.val, result->y0.val, x1, y1, 0.1, result->flip ? "<" : ">"); 

    fnhist->SetName("Fitfn"); 
    fnhist->SetMaximum(fithist->GetMaximum()); 
    fnhist->DrawCopy("colz"); 
    arrow->Draw(); 
    canvas->cd(3); 

    residualhist->SetTitle("#chi^{2} value"); 
    residualhist->SetName("chisq"); 
    residualhist->DrawCopy("colz"); 

    delete fnhist; 
    delete residualhist; 

    canvas->cd(4); 

    
    TH2 * rotated = MaxCamImageTools::rotateInterpolate(fithist, result->phi.val, result->x0.val, result->y0.val,0, "bicubic"); 
    TH2 * rotated_gm = MaxCamImageTools::rotateInterpolate(local_gainmap, result->phi.val, result->x0.val, result->y0.val,0, "bicubic"); 
    TH2 * norml =  (TH2*) local_gainmap->Clone("norml"); 

    for (int i = 1; i <= norml->GetNbinsX(); i++)
    {
      for (int j = 1; j <= norml->GetNbinsY(); j++)
      {
        norml->SetBinContent(i,j,1); 
      }
    }
 
    TH2 * rotated_norml = MaxCamImageTools::rotateInterpolate(norml,result->phi.val, result->x0.val, result->y0.val,0,"bicubic"); 

   
    TH1 * longi = rotated->ProjectionX("long"); 
    TH1 * longi_gm = rotated_gm->ProjectionX("long_gm"); 
    TH1 * longi_norml = rotated_norml->ProjectionX("long_norml"); 
    longi_gm->Divide(longi_norml); 
    double average_gain = longi_gm->Integral()/longi_gm->GetNbinsX(); 
    std::cout << "average gain: " << average_gain << std::endl; 

    TF1 func1d("func1d", &DmtpcMath::lineSegmentConvolvedWithGaussian, longi->GetXaxis()->GetXmin(), longi->GetXaxis()->GetXmax(), 5); 

    double b = result->z0.val- m*result->x0.val/binwidth;  
    double line_params[]  = {m/binwidth*average_gain,b*average_gain,result->x0.val,result->x0.val+result->range.val,result->sigma.val}; 
    func1d.SetParameters(line_params); 

    longi->SetMaximum(TMath::Max(result->z0.val,result->z1.val) * 1.25*average_gain); 
    longi->DrawCopy(); 
    func1d.DrawCopy("lpsame"); 
    TF1 funcline("funcline","pol1",result->x0.val, result->x0.val+result->range.val); 
    funcline.SetParameters(b*average_gain,m/binwidth*average_gain); 
    funcline.SetLineColor(3); 
    funcline.DrawCopy("lpsame"); 
    double * pdfx = new double[longi->GetNbinsX()]; 
    double * pdfy = new double[longi->GetNbinsX()]; 

    double line_params_unaverage[] = {m/binwidth,b,result->x0.val,result->x0.val+result->range.val,result->sigma.val}; 
    for (int ii = 1; ii <= longi->GetNbinsX(); ii++)
    {
      pdfx[ii-1] = longi->GetBinCenter(ii);  
      pdfy[ii-1] = DmtpcMath::integralOfLineSegmentConvolvedWithGaussian(longi->GetXaxis()->GetBinLowEdge(ii), longi->GetXaxis()->GetBinLowEdge(ii+1), line_params_unaverage)/binwidth*longi_gm->GetBinContent(ii); 
    }
      
    TGraph  * g = new TGraph(longi->GetNbinsX(), pdfx, pdfy); 
    g->SetMarkerColor(4); 
    g->Draw("*psame"); 


    double pts_x[2]={result->x0.val, result->x0.val + result->range.val}; 
    double pts_y[2]={result->z0.val*average_gain, result->z1.val*average_gain}; 
    double pts_x_err[2]={result->x0.err, sqrt(pow(result->x0.err,2) + pow(result->range.err,2))}; 
    double pts_y_err_dn[2]={fabs(result->z0.errDn)*average_gain, fabs(result->z1.errDn)*average_gain}; 
    double pts_y_err_up[2]={result->z0.errUp*average_gain, result->z1.errUp*average_gain}; 

    TGraphAsymmErrors * pts = new TGraphAsymmErrors(2,pts_x,pts_y,pts_x_err,pts_x_err,pts_y_err_dn,pts_y_err_up); 
    pts->SetLineColor(8); 
    pts->SetMarkerColor(3); 
    pts->Draw("psame"); 

    canvas->cd(5); 


    TH1 * tranv = rotated->ProjectionY("tranv"); 
    TH1 * tranv_gm = rotated_gm->ProjectionY("tranv_gm"); 
    TH1 * tranv_norml = rotated_norml->ProjectionY("tranv_normal"); 
    tranv_gm->Divide(tranv_norml); 
    double * pdf_t_x = new double [tranv->GetNbinsX()]; 
    double * pdf_t_y = new double [tranv->GetNbinsX()]; 

    TF1 gaus("gaus","gaus", tranv->GetXaxis()->GetXmin(), tranv->GetXaxis()->GetXmax()); 
    gaus.SetParameters(average_gain*min->X()[PAR_E] *binwidth / sqrt(2*TMath::Pi()) / result->sigma.val,  result->y0.val, result->sigma.val); 
    for (int ii = 1; ii <= tranv->GetNbinsX(); ii++)
    {
      pdf_t_x[ii-1] = tranv->GetBinCenter(ii); 
      pdf_t_y[ii-1] = (TMath::Erf((tranv->GetBinLowEdge(ii+1) - result->y0.val)/(result->sigma.val * sqrt(2))) - TMath::Erf((tranv->GetBinLowEdge(ii) - result->y0.val)/(result->sigma.val * sqrt(2))) ) * min->X()[PAR_E] / (sqrt(2*TMath::Pi())) * tranv_gm->GetBinContent(ii); 
    }

    TGraph * gt = new TGraph(tranv->GetNbinsX(), pdf_t_x, pdf_t_y); 
    gt->SetMarkerColor(4); 
    tranv->DrawCopy(); 
    gt->Draw("psame"); 

    gaus.DrawCopy("lpsame"); 

    canvas->cd(6); 

    TPaveText * pt = new TPaveText(0.05,0.2,0.95,0.95); 
    pt->AddText( TString::Format("#chi^{2}/ndof : %f/%d  = %f", result->chisq, result->ndof, result->chisq/result->ndof)); 
    //pt->AddText( TString::Format("x0: %f +/-%f (%f,%f)", result->x0.val, result->x0.err, result->x0.errDn, result->x0.errUp)); 
    pt->AddText( TString::Format("x0: %f +/-%f", result->x0.val, result->x0.err)); 
    //pt->AddText( TString::Format("y0: %f +/-%f (%f,%f)", result->y0.val, result->y0.err, result->y0.errDn, result->y0.errUp)); 
    pt->AddText( TString::Format("y0: %f +/-%f", result->y0.val, result->y0.err)); 
    pt->AddText( TString::Format("z0: %f +/-%f (%f,%f)", result->z0.val, result->z0.err, result->z0.errDn, result->z0.errUp)); 
    pt->AddText( TString::Format("#Delta z: %f +/-%f (%f,%f)", result->delta_z.val, result->delta_z.err, result->delta_z.errDn, result->z1.errUp)); 
    pt->AddText( TString::Format("z1: %f +/-%f (%f,%f)", result->z1.val, result->z1.err, result->z1.errDn, result->z1.errUp)); 
    //pt->AddText( TString::Format("phi: %f +/-%f (%f,%f)", result->htphi.val, result->phi.err, result->phi.errDn, result->phi.errUp)); 
    pt->AddText( TString::Format("phi: %f +/-%f", result->htphi.val, result->phi.err)); 
   // pt->AddText( TString::Format("#sigma: %f +/-%f (%f,%f)", result->sigma.val, result->sigma.err, result->sigma.errDn, result->sigma.errUp)); 
    pt->AddText( TString::Format("#sigma: %f +/-%f pixels", result->sigma.val, result->sigma.err)); 
    pt->AddText( TString::Format("E: %f +/- %f adu (%f +/- %f keVee)", result->E.val, result->E.err, result->Ephys.val, result->Ephys.err)); 
    pt->AddText( TString::Format("Range: %f +/- %f pixels (%f +/- %f mm)", result->range.val, result->range.err, result->range.val * lengthcal, result->range.err * lengthcal)); 
    pt->AddText( TString::Format("Minos?: %s  LikelihoodIntegral?: %s", result->minos ? "true" : "false", do_integral ? "true" : "false"));
     pt->AddText(TString::Format("Z-test Probability: %f", result->probFast));  
    if (do_integral)
    {
      pt->AddText(TString::Format("Implied Probability: %f +/- %f", result->prob, result->probError));  
    }

    pt->Draw(); 

    delete pdfx;
    delete pdfy;
    delete pdf_t_x;
    delete pdf_t_y;
    delete tranv;
    delete longi; 
    delete longi_gm; 
    delete longi_norml; 
    delete tranv_gm; 
    delete tranv_norml; 
    delete rotated; 
    delete rotated_gm; 
    //delete g; 
  }


  if (verbose)
  {
    std::cout <<  "Chisq/ndof: " << result->chisq / result->ndof << std::endl <<
                  "x0: " << result->x0.val << " +/-" << result->x0.err << std::endl <<
                  "y0: " << result->y0.val << " +/-" << result->y0.err << std::endl <<
                  "z0: " << result->z0.val << " +/-" << result->z0.err << std::endl <<
                  "deltaz: " << result->delta_z.val << " +/-" << result->delta_z.err << std::endl <<
                  "z1: " << result->z1.val << " +/-" << result->z1.err << std::endl <<
                  "range: " << result->range.val << " +/-" << result->range.err << std::endl <<
                  "E: " << result->E.val << std::endl <<
                  "phi " << result->htphi.val << std::endl <<
                  "fast probability " << result->probFast << endl; 

    if (do_integral)
    {

     std::cout << " total integral: " << result->likelyIntegralTotal << " +/- " << result->likelyIntegralTotalError << 
                  " same side integral: " << result->likelyIntegralSame << " +/- " << result->likelyIntegralSameError << 
                  " implied probability: " << result->prob <<  " +/- " << result->probError << std::endl; 
    }

 }



  return 0; 

}
