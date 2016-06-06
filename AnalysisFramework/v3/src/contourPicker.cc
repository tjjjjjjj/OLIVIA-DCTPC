#include  "contourPicker.hh"
#include <iostream>

using std::cout;
using std::endl;



ContourPicker::ContourPicker(int n, 
                             double * levels, double tol)
{
  setLevels(n,levels); 
  setTol(tol); 
}

ContourPicker::~ContourPicker()
{}


static double getPortionAbove(const TH2 * hist, const double z)
{
  double ret = 0; 
  for (int x = 1; x <= hist->GetNbinsX(); x++)
  {
    for (int y = 1; y <= hist->GetNbinsY(); y++)
    {
      double val = hist->GetBinContent(x,y); 
      if (val >= z)
      { 
        ret+= val -z;
      }
    }
  }

  return ret; 
}


TH2F * ContourPicker::makeContourHist (const TH2 * hist, char  * name)
{
  double z = hist->GetMaximum();  
  double fraction = 0.; 
  double total = hist->Integral(); 

  TH2F * ret = new TH2F(name,name,
                        hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(),
                        hist->GetNbinsY(), hist->GetYaxis()->GetXmin(), hist->GetYaxis()->GetXmax());  
  

  while (fraction < 1) 
  {
    z-= _tol; 
    double portion = getPortionAbove(hist,z); 
    fraction = portion / total;  
    for (int x = 1; x <= hist->GetNbinsX(); x++)
    {
      for (int y = 1; y <= hist->GetNbinsY(); y++)
      {
        if (ret->GetBinContent(x,y) == 0 && hist->GetBinContent(x,y) > z) ret->SetBinContent(x,y,fraction); 
      }
    }
  }

  return ret; 
}

double * ContourPicker::getContours(const TH2 * hist)
{
  double total = hist->Integral(); 
  if (CNTRDBG) cout << "Integral : " << total << endl; 
  double fraction = 0.;  
  double z = hist->GetMaximum(); 
  double * ret = new double[_n]; 

  for (int i = 0; i < _n ; i++)
  {
    double level = _levels[i];      

    while (fraction < level)
    {
      z-=_tol;  
      double portion = getPortionAbove(hist,z); 
  //    if (CNTRDBG) cout << " Portion above " << z << " is " << portion << endl; 
      fraction = portion/total; 
    }
    if (CNTRDBG) cout << "Level " << level << " value is " << z << endl;    
    ret[_n - 1 - i] = z  ;
  }

  return ret; 
}

void ContourPicker::drawHistogramContours(TH2 * hist, const bool same)
{
  double * contours = getContours(hist); 
  hist->SetContour(_n, contours);  
  TString options = "cont z list";  
  if (same) options = "same " + options; 
  hist->Draw(options); 
}
