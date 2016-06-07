#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "GausConvSpline.hh" 
#include "DmtpcMath.hh"
#include "TVector3.h"
#include "TTimeStamp.h"
#include "TGraph.h"
#include "DmtpcProjection.hh"
#include "TCanvas.h"
#include "TCutG.h"
#include "TH1.h"
#include <vector>
#include <set> 
#include <cmath>
#include <iostream>
#include <limits.h>
#include "MaxCamClusterImage.hh"
#include "McWimp.hh"
#include "MaxCamImageTools.hh"
#include "Math/Minimizer.h" 
#include "Math/Functor.h" 
#include "Math/Factory.h" 

using namespace std;

MaxCamClusterImage::MaxCamClusterImage()
{
   _nbinsx=0;
   _nbinsy=0;
   _nbins=0;
   image_rms = -1;
}

MaxCamClusterImage::MaxCamClusterImage(TH2* image, TDatime* time)
{
  setTime(time);  
  _image =  image;  

   //   cout << "There are " << _pixels.size() << " entries in the cluster \n";

   _nbinsx = _image->GetNbinsX()+2;
   _nbinsy = _image->GetNbinsY()+2;
   _nbins = _nbinsx*_nbinsy;

   image_rms = -1;
}

MaxCamClusterImage::MaxCamClusterImage(TH2* image, TTimeStamp* time)
{
   setTime(time); 
   _image =  image;  

   //   cout << "There are " << _pixels.size() << " entries in the cluster \n";

   _nbinsx = _image->GetNbinsX()+2;
   _nbinsy = _image->GetNbinsY()+2;
   _nbins = _nbinsx*_nbinsy;

   image_rms = -1;
}

MaxCamClusterImage::~MaxCamClusterImage()
{
   delete _image;
   //_pixels.clear();
   //_pixelsred.clear();
}


void MaxCamClusterImage::addCluster(vector < int > newpx, vector<int> redpx,  vector<vector<char> > * frac)
{
  if (redpx.size() > 0)
  {
    _pixelsred.push_back(redpx);
  }
  else
  {
    _pixelsred.push_back(newpx);
  }

  _pixels.push_back(newpx);  
}

double MaxCamClusterImage::getIntegral(int i) const
{
   vector <int> px = getCluster(i);
   
   double sum=0;
   for(int j=0; j<int(px.size());j++)
   {
      sum += _image->GetBinContent(px[j]);
   }

   return sum;
}

double MaxCamClusterImage::getIntegralWithGainMap(int i, const DmtpcGainMap* gainMap, double mingain) const
{
   if(!gainMap || gainMap->getGainMap()->GetNbinsX() != _image->GetNbinsX() ||
      gainMap->getGainMap()->GetNbinsY() != _image->GetNbinsY())
   {
      cout << "Gain Map does not exist or does match image dimensions; proceeding with non-weighted integral." << endl;
      return getIntegral(i);
   }
   else
   {
      vector <int> px = getCluster(i);
      
      double sum=0;
      for(int j=0; j<int(px.size());j++)
      {
        double g = gainMap->getGainMap()->GetBinContent(px[j]); 
        if (g > mingain)
          sum += _image->GetBinContent(px[j])/(g);
      }
      
      return sum;
   }
}
   
double MaxCamClusterImage::getLength(int i, double &x1, double &y1, double &x2, double &y2) const
{
   vector <int> px = getClusterRed(i);
   
   double length = 0;
   for(int k=0; k<int(px.size()); k++)
   {
      for(int j=i+1; j<int(px.size());j++)
      {
	 int ibinx = px[k]%_nbinsx;
	 int ibiny = ((px[k]-ibinx)/_nbinsx)%_nbinsy;
	 int jbinx = px[j]%_nbinsx;
	 int jbiny = ((px[j]-jbinx)/_nbinsx)%_nbinsy;

	 double ix = _image->GetXaxis()->GetBinCenter(ibinx);
	 double iy = _image->GetYaxis()->GetBinCenter(ibiny);
	 double jx = _image->GetXaxis()->GetBinCenter(jbinx);
	 double jy = _image->GetYaxis()->GetBinCenter(jbiny);

	 double testlen = sqrt((ix-jx)*(ix-jx)+(iy-jy)*(iy-jy));

	 if(testlen > length)
	 {
	    length=testlen;
	    x1=ix;
	    x2=jx;
	    y1=iy;
	    y2=jy;
	 }

      }
   }

   return length;
}

double MaxCamClusterImage::getDiffusedLength(int i, double &x1, double &y1, double &x2, double &y2) const
{
   vector <int> px = getCluster(i);
   
   double length = 0;
   for(int k=0; k<int(px.size()); k++)
   {
      for(int j=i+1; j<int(px.size());j++)
      {
	 int ibinx = px[k]%_nbinsx;
	 int ibiny = ((px[k]-ibinx)/_nbinsx)%_nbinsy;
	 int jbinx = px[j]%_nbinsx;
	 int jbiny = ((px[j]-jbinx)/_nbinsx)%_nbinsy;

	 double ix = _image->GetXaxis()->GetBinCenter(ibinx);
	 double iy = _image->GetYaxis()->GetBinCenter(ibiny);
	 double jx = _image->GetXaxis()->GetBinCenter(jbinx);
	 double jy = _image->GetYaxis()->GetBinCenter(jbiny);

	 double testlen = sqrt((ix-jx)*(ix-jx)+(iy-jy)*(iy-jy));

	 if(testlen > length)
	 {
	    length=testlen;
	    x1=ix;
	    x2=jx;
	    y1=iy;
	    y2=jy;
	 }

      }
   }

   return length;
}

bool MaxCamClusterImage::isInCluster(int i, int bin, bool red) const
{
   vector <int> px = red ? getClusterRed(i) : getCluster(i);
   for(int j=0; j<int(px.size()); j++)
   {  
      if(px[j]==bin) return true;
   }
   return false;
}

bool MaxCamClusterImage::isContainedInRegion(int itrk, TCutG* region) const 
{
  vector<int> pix = getCluster(itrk);
  for (unsigned int ii=0; ii<pix.size(); ii++) {
    int xbin, ybin, junk;
    getImage()->GetBinXYZ(pix.at(ii), xbin, ybin, junk);
    double xpix = getImage()->GetXaxis()->GetBinCenter(xbin);
    double ypix = getImage()->GetXaxis()->GetBinCenter(ybin);
    if (!region->IsInside(xpix, ypix)) return false;
  }
  return true;
}

bool MaxCamClusterImage::intersectsRegion(int itrk, TCutG* region) const
{
  vector<int> pix = getCluster(itrk);
  for (unsigned int ii=0; ii<pix.size(); ii++) {
    int xbin, ybin, junk;
    getImage()->GetBinXYZ(pix.at(ii), xbin, ybin, junk);
    double xpix = getImage()->GetXaxis()->GetBinCenter(xbin);
    double ypix = getImage()->GetXaxis()->GetBinCenter(ybin);
    if (region->IsInside(xpix, ypix)) return true;
  }
  return false;
}

void MaxCamClusterImage::getXY(int i, double &x, double &y) const
{

   x = 0;
   y = 0; 

   vector <int> px = getCluster(i);
   int xx,yy;
   double xw,yw;

   int nx=_image->GetNbinsX()+2;
   int ny=_image->GetNbinsY()+2;

   for(int j=0; j<int(px.size()); j++)
   {
      xx = px[j]%nx;
      yy = ((px[j]-xx)/nx)%ny;
      xw=_image->GetXaxis()->GetBinCenter(xx);
      yw=_image->GetYaxis()->GetBinCenter(yy);
      double content = _image->GetBinContent(px[j]);
      x+=content*xw;
      y+=content*yw;

   }
   double integral = getIntegral(i); 
   x/=integral;
   y/=integral;

}

double MaxCamClusterImage::getPhi(int i) const
{
   double x,y,x1,x2,y1,y2,length;
   length=getLength(i,x1,y1,x2,y2);
   getXY(i,x,y);

   bool isOne;
   if(pow(x-x1,2)+pow(y-y1,2) < pow(x-x2,2)+pow(y-y2,2)){isOne=true;}
   else {isOne=false;}

   double phi;

   if(isOne)
   {
      double deltax = x2-x1;
      double deltay = y2-y1;
      phi=atan2(deltay,deltax);
   }
   else
   {
      double deltax = x1-x2;
      double deltay = y1-y2;
      phi=atan2(deltay,deltax);
   }

   return phi;

}

double MaxCamClusterImage::getPhi2(int i, int nerode,int ndilate ) const {
  //returns angle in xy plane, no directionality yet
 
  vector <int> px = getCluster(i);
  if (nerode || ndilate)
  {
    std::set<int> moar_pix; 
    moar_pix.insert(px.begin(), px.end()); 

    for (int i = 0; i < nerode; i++)
    {
       MaxCamImageTools::erode(_image, &moar_pix);
    }

    for (int i = 0; i < ndilate; i++)
    {
       MaxCamImageTools::dilate(_image, &moar_pix);
    }


    px.clear();
    for (std::set<int>::iterator it = moar_pix.begin(); it!= moar_pix.end(); it++)
    {
      px.push_back(*it); 
    }
  }

  int size = (int) px.size();
  double meanx = 0, meany = 0, xsq = 0, ysq = 0, xy = 0, sum = 0;

  for (int j = 0; j < size; j++){

    int xbin = px[j]%_nbinsx;
    int ybin = px[j]/_nbinsx;
    double content = _image->GetBinContent(xbin,ybin);
    double xpos = _image->GetXaxis()->GetBinCenter(xbin);
    double ypos = _image->GetYaxis()->GetBinCenter(ybin);

    meanx += content*xpos;
    meany += content*ypos;
    xsq += content*xpos*xpos;
    ysq   += content*ypos*ypos;
    xy    += content*xpos*ypos;
    sum   += content;

  }
  meanx /= sum; meany /= sum; xsq /= sum; ysq /= sum; xy /= sum;

  double varx = xsq - meanx*meanx;//x-variance
  double vary = ysq - meany*meany;//y-variance
  double cov  = xy - meanx*meany;//covariance
  
  //unnormalized eigenvectors of the matrix I = [vary, - cov; - cov, varx] are in the form (A,B+-)
  double A      = 2 * cov;
  double Bplus  = vary - varx - sqrt((varx - vary)*(varx - vary) + 4.*cov*cov);
  double Bminus = vary - varx + sqrt((varx - vary)*(varx - vary) + 4.*cov*cov);
  double phi1 = atan2(Bplus,A);
  double phi2 = atan2(Bminus,A);

  //cout << "Phi 1: " << phi1*180/TMath::Pi() << "\t Phi 2: " << phi2*180/TMath::Pi() << endl;

  double var1 = varx * cos(phi1)*cos(phi1) + vary * sin(phi1)*sin(phi1) + cov * sin(2.*phi1);
  double var2 = varx * cos(phi2)*cos(phi2) + vary * sin(phi2)*sin(phi2) + cov * sin(2.*phi2);

  double phi;
  if (var1 < var2) phi = phi2;
  else phi = phi1;

  if(phi > TMath::Pi()/2) phi = phi-TMath::Pi();

  //cout << "Phi : " << phi*180/TMath::Pi() << endl;

  //double skewness = getSkewness(i,phi);
  double skewness = getAsymmetry(i,phi);

  //cout << "Skew: " << skewness << endl;
  if(skewness > 0)
  {
     if(phi <= 0)				
	phi=TMath::Pi()+phi;
     else if(phi > 0)
	phi=phi-TMath::Pi();
  }

//  cout << "New Phi: " << phi*180/TMath::Pi() << endl;
  
  return phi;

}




void MaxCamClusterImage::getEllipseAxes(int i, double& a, double& b)
{
   vector <int> px = getCluster(i);
  int size = (int) px.size();
  double meanx = 0, meany = 0, xsq = 0, ysq = 0, xy = 0, sum = 0;

  for (int j = 0; j < size; j++){

    int xbin = px[j]%_nbinsx;
    int ybin = px[j]/_nbinsx;
    double content = _image->GetBinContent(xbin,ybin);
    double xpos = _image->GetXaxis()->GetBinCenter(xbin);
    double ypos = _image->GetYaxis()->GetBinCenter(ybin);

    meanx += content*xpos;
    meany += content*ypos;
    xsq += content*xpos*xpos;
    ysq   += content*ypos*ypos;
    xy    += content*xpos*ypos;
    sum   += content;

   

  }
  meanx /= sum; meany /= sum; xsq /= sum; ysq /= sum; xy /= sum;

  double varx = xsq - meanx*meanx;//x-variance
  double vary = ysq - meany*meany;//y-variance
  double cov  = xy - meanx*meany;//covariance
  cout<<"covariance "<<cov<<endl;
  a = (varx+vary)/2 + sqrt((varx-vary)*(varx-vary)/4-cov*cov);
  //a = (varx+vary) + sqrt((varx-vary)*(varx-vary)-(4*cov*cov)); 
  a = sqrt(a);
  
  b = (varx+vary)/2 - sqrt((varx-vary)*(varx-vary)/4-cov*cov);
  //b = (varx+vary) - sqrt((varx-vary)*(varx-vary)-(4*cov*cov));  
  b = sqrt(b);
  cout<<a<<" "<<b<<endl;

}


//Calculates phi by finding angle with maximum sigma.  Similar precision to getPhi2 but slightly less accurate.
double MaxCamClusterImage::getPhi3(int i) const{
   vector <int> px = getCluster(i);
   int size     = (int) px.size();
   double meanx = 0, meany = 0, xsq = 0, ysq = 0, xy = 0, sum = 0;
   
   for (int i = 0; i < size; i++){
      
      int xbin       = px[i]%(_nbinsx);
      int ybin       = px[i]/(_nbinsx);
      double content = _image->GetBinContent(xbin,ybin);
      double xpos    = _image->GetXaxis()->GetBinCenter(xbin);
      double ypos    = _image->GetYaxis()->GetBinCenter(ybin);
      
      meanx += content*xpos;
      meany += content*ypos;
      xsq   += content*xpos*xpos;
      ysq   += content*ypos*ypos;
      xy    += content*xpos*ypos;
      sum   += content;
   }
   
   meanx /= sum; meany /= sum; xsq /= sum; ysq /= sum; xy /= sum;
   
   double varx = xsq - meanx*meanx;//x-variance
   double vary = ysq - meany*meany;//y-variance
   double cov  = xy - meanx*meany;//covariance
   double phi1 = 0.5*atan2(2.*cov,varx-vary);//see getPhi() for explanation of this
   double phi2 = phi1 + 0.5*TMath::Pi();
   double var1 = varx * cos(phi1)*cos(phi1) + vary * sin(phi1)*sin(phi1) + cov * sin(2.*phi1);
   double var2 = varx * cos(phi2)*cos(phi2) + vary * sin(phi2)*sin(phi2) + cov * sin(2.*phi2);
   
   double phi;
   if (var1 < var2){
      phi       = phi2;
      //    rmswidth  = sqrt(var1);
      //   rmslength = sqrt(var2);
   }
   else {
      phi       = phi1;
      //    rmswidth  = sqrt(var2);
      //    rmslength = sqrt(var1);
   }
   
   return phi;
}

double MaxCamClusterImage::getLength2(int i, double theta, double pixpermm){
  //calculates length of cluster along on axis at angle theta above x.
  vector <int> px = getCluster(i);
  int size = (int) px.size();
  double minx = _image->GetXaxis()->GetBinLowEdge(1);
  double miny = _image->GetYaxis()->GetBinLowEdge(1);
  
  double mindist = 0;
  double maxdist = 0;

  for (int i = 0; i < size; i++){
    //save only the two extreme values
    int xbin = px[i] % _nbinsx;
    int ybin = px[i] / _nbinsx;
    double dist = (_image->GetXaxis()->GetBinCenter(xbin) - minx)*cos(theta)/pixpermm + 
      (_image->GetYaxis()->GetBinCenter(ybin) - miny)*sin(theta)/pixpermm;

    if (dist < mindist || i == 0) mindist = dist;
    if (dist > maxdist || i == 0) maxdist = dist;
  }
  double length = maxdist - mindist;
  return length;
}

double MaxCamClusterImage::getSkewness(int i,double theta){

   vector <int> px = getCluster(i);
   int size = (int) px.size();
   //x*cos(theta) + y*sin(theta)
   double mean = 0;
   double total = 0;

   double offset = 0;
  
   for (int j = 0; j < size; j++){
      
      int xbin = px[j] % _nbinsx;
      int ybin = ((px[j]-xbin)/_nbinsx)%_nbinsy;

      float distance=xbin*cos(theta)+ybin*sin(theta);
      
      float content = _image->GetBinContent(px[j])+offset;

      if(content >0)
      {
	 mean += content*distance;
	 total += content;
      }
   }

   //cout << " w*x: " << mean << endl;
   //cout << " w: " << total << endl;
   mean /=total;
   
   //cout << "Mean: " << mean << endl;
   
   double moment3 = 0;
   double moment2 = 0;
   for (int j = 0; j < size; j++){
      
      int xbin = px[j] % _nbinsx;
      int ybin = ((px[j]-xbin)/_nbinsx)% _nbinsy;
      
      double dist = xbin*cos(theta) + ybin*sin(theta);

      double content = _image->GetBinContent(px[j])+offset;
      if(content > 0)
      {	 
	 moment3 += content*(dist-mean)*(dist-mean)*(dist-mean);
	 //moment3 += content*(dist)*(dist)*(dist);
	 moment2 += content*(dist-mean)*(dist-mean); 
      }
   }

   moment3 /= total;
   moment2 /= total;
   //cout << "second moment around mean: " << moment2 << endl;
   //cout << "sigma: " << sqrt(moment2) << endl;

   double skewness = moment3 / pow(moment2,1.5);

   //cout << "Moment 2: " << moment2 << "\t Moment 3: " 
	//<< moment3 << "\t Skewness: " << skewness << endl;
      
   return skewness;
   
}

double MaxCamClusterImage::getAsymmetry(int i, double phi) const
{
   vector <int> px = getCluster(i);
   int size = (int) px.size();
   double min = 10000;
   double max = -10000;
   double xmin = 10000;
   double xmax = -10000;
   double ymin = 10000;
   double ymax = -10000;

   for(int j=0; j<size; j++)
   {

      int xbin = px[j] % _nbinsx;
      int ybin = ((px[j]-xbin)/_nbinsx)%_nbinsy;

      double distance=xbin*cos(phi)+ybin*sin(phi);
      if(distance < min)
	 min = distance;
      if(distance > max)
	 max = distance;

      if(xbin < xmin)
	 xmin = xbin;
      if(xbin > xmax)
	 xmax = xbin;

      if(ybin < ymin)
	 ymin = ybin;
      if(ybin > ymax)
	 ymax = ybin;
      
   }

   double midpoint = (min+max)/2;
   double xmidpoint = (xmin+xmax)/2;
   double ymidpoint = (ymin+ymax)/2;

  //cout << "Min: " << min << "\t Max: " << max 
	//<< "\t Mid: " << midpoint << endl;
  // cout << "XMin: " << xmin << "\t XMax: " << xmax 
	//<< "\t XMid: " << xmidpoint << endl;
  // cout << "YMin: " << ymin << "\t YMax: " << ymax 
	//<< "\t YMid: " << ymidpoint << endl;

   double above=0;
   double below=0;
   
   double xabove=0;
   double xbelow=0;
   
   double yabove=0;
   double ybelow=0;
   
   //cout << "E: " << getIntegral(i) << endl;
   for(int j=0; j<size; j++)
   {

      int xbin = px[j] % _nbinsx;
      int ybin = ((px[j]-xbin)/_nbinsx)%_nbinsy;

      double distance=xbin*cos(phi)+ybin*sin(phi);
      double content=_image->GetBinContent(px[j]);
      
      if(distance < midpoint)
	 below+=content;
      if(distance > midpoint)
	 above+=content;

      if(xbin < xmidpoint)
	 xbelow+=content;
      if(xbin > xmidpoint)
	 xabove+=content;

      if(ybin < ymidpoint)
	 ybelow+=content;
      if(ybin > ymidpoint)
	 yabove+=content;
      
   }
   

   //cout << "Above: " << above << "\t Below: " << below << endl;
   //cout << "XAbove: " << xabove << "\t XBelow: " << xbelow << endl;
   //cout << "YAbove: " << yabove << "\t YBelow: " << ybelow << endl;
   
   if(below > above)
      return -1;
   else if (above > below)
      return 1;
   else
      return 0;
   
}

bool MaxCamClusterImage::hitsInactive(int i, double rout)  const
{

  vector<int> px = getCluster(i); 
  for(vector<int>::iterator it = px.begin(); it!=px.end(); it++)
  {
      int bin = *it;  
      int xbin = bin % _nbinsx; 
      int ybin = bin / _nbinsx; 

      double x = _image->GetXaxis()->GetBinCenter(xbin); 
      double y = _image->GetYaxis()->GetBinCenter(ybin); 

      double rr = sqrt(x*x + y*y); 

      if (rr > rout) return true; 

  }

  return false; 
}


bool MaxCamClusterImage::hitsVeto(int i, double rin, double rout)  const
{

  vector<int> px = getCluster(i); 
  for(vector<int>::iterator it = px.begin(); it!=px.end(); it++)
  {
      int bin = *it;  
      int xbin = bin % _nbinsx; 
      int ybin = bin / _nbinsx; 

      double x = _image->GetXaxis()->GetBinCenter(xbin); 
      double y = _image->GetYaxis()->GetBinCenter(ybin); 

      double rr = sqrt(x*x + y*y); 

      if (rr  > rin && rr < rout) return true; 

  }

  return false; 
}

bool MaxCamClusterImage::hitsEdge(int i){

   vector <int> px = getCluster(i);
   bool hit = false;
   
   int count = 0;
   while (count < (int)px.size()){
      int bin = px[count];
      int xbin = bin % _nbinsx;
      int ybin = bin / _nbinsx; 
      
      if (xbin == 1 || xbin == (_nbinsx-2)) {
	 hit = true;
	 break;
      }
      else if (ybin == 1 || ybin == (_nbinsy-2)) {
	 hit = true;
	 break;
      }
      count++;
   }
   
   return hit;
}



void MaxCamClusterImage::changeImage(TH2* newimage)
{
   image_rms = -1;
   int xbinso = _image->GetNbinsX();
   int xbinsn = newimage->GetNbinsX();
   int ybinso = _image->GetNbinsY();
   int ybinsn = newimage->GetNbinsY();

   vector<int> newpx;
   vector < vector <int> > newbigpx;

   if (xbinso > xbinsn || ybinso > ybinsn)
   {cout << "Unable to change image, new image has fewer bins than old image" << endl; return;}

   if(xbinso == xbinsn && ybinso == ybinsn)
   {
      cout << "bins equal" << endl;
      _image =  newimage;
      return;
   }

//    if(xbinsn%xbinso !=0 || ybinsn%ybinso != 0)
//    {cout << "Unable to change image, number of bins is not evenly divisible" << endl; return;}
   else
   {
      int x,y;
      int xfactor = (xbinsn-(xbinsn%xbinso))/xbinso;
      int yfactor = (ybinsn-(ybinsn%ybinso))/ybinso;
      for(int p=0; p<int(_pixels.size()); p++)
      {
	 vector<int> px= _pixels[p];
	 newpx.clear();

	 for(int q=0; q<int(px.size()); q++)
	 {
	    x = px[q]%(xbinso+2);
	    y = ((px[q]-x)/(xbinso+2))%(ybinso+2);
	    for(int j=(x-1)*xfactor+1; j<=x*xfactor; j++)
	    {
	       for(int k=(y-1)*yfactor+1; k<=y*yfactor; k++)
	       {
		  int bin = newimage->GetBin(j,k);
		  newpx.push_back(bin);
	       }
	    }

	 }

	 newbigpx.push_back(newpx);

      }
      _pixels = newbigpx;
      _pixelsred = newbigpx;
      _image =  newimage;
      _nbinsx = _image->GetNbinsX()+2;
      _nbinsy = _image->GetNbinsY()+2;
      _nbins = _nbinsx*_nbinsy;
      return;

   }

   
}

void MaxCamClusterImage::applyRedThreshold(double threshold)
{

  vector<vector<int> > newbigpxred(_pixels.size()); 

  for (unsigned int i = 0; i < _pixels.size(); i++)
  {
    vector<int> newpxred; 

    for (unsigned int j = 0; j < _pixels[i].size(); j++)
    {

      if (_image->GetBinContent(_pixels[i][j]) > threshold)
      {
        newpxred.push_back(_pixels[i][j]); 
      }
    }

    newbigpxred[i]=newpxred; 
  }

  _pixelsred = newbigpxred; 

}


void MaxCamClusterImage::changeImageWithThreshold(TH2* newimage, double threshold)
{
  image_rms = -1;
  int xbinso = _image->GetNbinsX();
  int xbinsn = newimage->GetNbinsX();
  int ybinso = _image->GetNbinsY();
  int ybinsn = newimage->GetNbinsY();
  
  vector<int> newpx;
  vector<int> newpxred;
  vector < vector <int> > newbigpx;
  vector < vector <int> > newbigpxred;
  
  if (xbinso > xbinsn || ybinso > ybinsn)
    {cout << "Unable to change image, new image has fewer bins than old image" << endl; return;}
  
  else
    {
      int x,y;
      int xfactor = (xbinsn-(xbinsn%xbinso))/xbinso;
      int yfactor = (ybinsn-(ybinsn%ybinso))/ybinso;
      for(int p=0; p<int(_pixels.size()); p++)
	{
	  vector<int> px= _pixels[p];
	  newpx.clear();
	  newpxred.clear();

	  for(int q=0; q<int(px.size()); q++)
	    {
	      x = px[q]%(xbinso+2);
	      y = ((px[q]-x)/(xbinso+2))%(ybinso+2);
	      for(int j=(x-1)*xfactor+1; j<=x*xfactor; j++)
		{
		  for(int k=(y-1)*yfactor+1; k<=y*yfactor; k++)
		    {
		      int bin = newimage->GetBin(j,k);
		      double bincontent=newimage->GetBinContent(j,k);
		      newpx.push_back(bin);
		      if (bincontent>threshold)
			newpxred.push_back(bin);
		    }
		}
	      
	    }
	  newbigpxred.push_back(newpxred);
	  newbigpx.push_back(newpx);
	  
	}
      _pixelsred = newbigpxred;
      _pixels = newbigpx;
      _image =  newimage;
      _nbinsx = _image->GetNbinsX()+2;
      _nbinsy = _image->GetNbinsY()+2;
      _nbins = _nbinsx*_nbinsy;
      return;
      
    }
}


double MaxCamClusterImage::getCygnusAngle(int i,double nang, 
                           CAMERA_ORIENTATION orientation, 
                           double lat, double lng, double phi,
                           double theta)
{

   //Convert nang to radians... it should be entered in degrees
   nang*= TMath::Pi()/180.;

   //Get phi if not specified by function call
   if (phi==DBL_MAX) phi = getPhi2(i);

   //Get theta if not specified by function call
   if (theta==DBL_MAX) theta= getTheta(i); 
  
   //Correct for orientation by rotating around Y axis
   if (orientation == BOTTOM)
   {
     theta = TMath::Pi() - theta; 
     if (phi >0) phi = TMath::Pi()-phi;
     if (phi <0) phi -= TMath::Pi();  
   }

  //calculate the angle of the recoilw ith respect to East
   double eang = phi + nang; 

   //Use McWimp to find cygnus direction
   McWimp* mcwimp  = new McWimp();
  
   //Set Coordinates 
   mcwimp->setLat(lat);
   mcwimp->setLon(lng);
   
   //figure out time
   mcwimp->setTime(_time);
 
   //vertically oriented detector
   mcwimp->setOpt("lv");
   TVector3 dettocyg = mcwimp->cosmicGun().Unit();
   dettocyg*=-1; 
  
   TVector3 
   angindet(sin(theta)*cos(eang),sin(theta)*sin(eang),cos(theta));

   double cygang = dettocyg.Dot(angindet)/(angindet.Mag()*dettocyg.Mag());
   delete mcwimp; 
   return cygang;

}

void MaxCamClusterImage::setTime(TTimeStamp* time)
{
   double Julian = (time->GetSec()+1E-9*time->GetNanoSec())/86400.0 + 2440587.5;
   _time = Julian-2451544.0;
}

void MaxCamClusterImage::setTime(TDatime* time)
{
   TTimeStamp* ts = new TTimeStamp(time->GetYear(),
				   time->GetMonth(),
				   time->GetDay(),
				   time->GetHour(),
				   time->GetMinute(),
				   time->GetSecond(),
				   0,
				   kTRUE,
				   4*60*60);

   double Julian = (ts->GetSec()+1E-9*ts->GetNanoSec())/86400.0 + 2440587.5;
   _time = Julian-2451544.0;
}

double MaxCamClusterImage::getPhi4(int i) const 
{
   const int npx = int(_pixels[i].size());
   double xx[npx], yy[npx];
   for(int j=0; j<npx; j++)
   {
      int xj = _pixels[i][j]%_nbinsx;
      int yj = ((_pixels[i][j]-xj)/_nbinsx)%_nbinsy;
      xx[j]=_image->GetXaxis()->GetBinCenter(xj);
      yy[j]=_image->GetYaxis()->GetBinCenter(yj);
   }

   TGraph* g = new TGraph(npx,xx,yy);
   g->Fit("pol1");
   TF1* line = g->GetFunction("pol1");
   double m = line->GetParameter(1);
   double phi = atan(m);



   return phi;

}

/*
typedef struct
{
  double mean; 
  double sigma; 
  double amplitude; 
  double chisq; 
  int row; 
} phi5_fit_result_t; 


bool phi5_cmp_sigma(phi5_fit_result_t i, phi5_fit_result_t j)
{
  return i.sigma < j.sigma; 
}

bool phi5_cmp_chisq(phi5_fit_result_t i, phi5_fit_result_t j)
{
  return i.chisq < j.chisq; 
}


*/ 




//Use Radon Transform 
/*double MaxCamClusterImage::getPhi5(int i, int extrapix, const DmtpcGainMap * map, int nbins, double vector_pct) const
{
  TH2 * clusthist = getClusterHist(i,extrapix,false,map,0,0,false); 
  TH2 * clusthist_padded =  MaxCamImageTools::zeroPadSquare(clusthist); 
  TH2 * transform = MaxCamImageTools::radonTransform(clusthist_padded, nbins); 

  vector< phi5_fit_result_t> fits( transform->GetNbinsY() ); 


  TF1 f("f","gaus", transform->GetXaxis()->GetXmin(), transform->GetXaxis()->GetXmax()); 
  
  TH1 * proj = new TH1F("proj","proj",transform->GetNbinsX(), transform->GetXaxis()->GetXmin(), transform->GetXaxis()->GetXmax()); 
  for (int i = 1; i <= transform->GetNbinsY(); i++)
  {
    for (int j = 1; j <= transform->GetNbinsX(); j++)
    {
//      cout << transform->GetBinContent(j,i) << " "; 
      proj->SetBinContent(j, transform->GetBinContent(j,i)); 
    }
//    cout << endl; 

    proj->Fit(&f,"NQ"); 

    fits[i-1].amplitude = f.GetParameter(0); 
    fits[i-1].mean = f.GetParameter(1); 
    fits[i-1].sigma = f.GetParameter(2); 
    fits[i-1].chisq = f.GetChisquare() / f.GetNDF(); 
    fits[i-1].row = i; 
  }


  //figure out chisq_cutoff 

  std::sort(fits.begin(), fits.end(), phi5_cmp_chisq); 
  double chisq_cutoff = fits[fits.size()/4].chisq; 
  std::sort(fits.begin(), fits.end(), phi5_cmp_sigma); 

  int which_fit = 0; 
  for (unsigned i = 0; i < fits.size(); i++)
  {
    if (fits[i].chisq < chisq_cutoff)
    {
      which_fit = i; 
      break; 
    }
  }
  //we now have either phi + pi/2 or phi - pi/2
  
  int test_row = ((fits[which_fit].row -1+ nbins/2) % nbins) + 1; 

  double phi = transform->GetYaxis()->GetBinLowEdge(test_row); 

  for (int j = 1; j <= transform->GetNbinsX(); j++)
  {
      double val = transform->GetBinContent(j,test_row); 
      proj->SetBinContent(j, val); 
  }

  
  double total = proj->Integral(); 
  double contained = proj->GetMaximum(); 
  int minbin = proj->GetMaximumBin(); 
  int maxbin = proj->GetMaximumBin(); 

  if (vector_pct < 0 || vector_pct > 1) 
  {
    std::cerr << "Invalid vector_pct, setting to 0.9" <<std::endl; 
    vector_pct = 0.9; 
  }

  //insist on even number of bins! 
  int nbins_used = 1; 
  while (!( contained > vector_pct*total  &&  (nbins_used %2) == 0))
  {

     int leftbin = minbin - 1; 
     int rightbin = maxbin + 1; 

     //choose right 
     if (leftbin == 0 || proj->GetBinContent(rightbin) > proj->GetBinContent(leftbin))
     {
       contained += proj->GetBinContent(++maxbin); 
     }
     else 
     {
       contained += proj->GetBinContent(--minbin); 
     }

     nbins_used++; 
  }


  bool flip = proj->Integral(minbin, minbin+nbins_used/2) < contained/2; 


  if (flip)
  {
    phi += M_PI; 

    if (phi > 2*M_PI) phi-=2*M_PI; 
  }
  delete proj; 
  delete clusthist; 
  delete clusthist_padded; 
  delete transform; 

  return phi; 
}

*/

double MaxCamClusterImage::getRMS(int i, double mean)
{
  vector<int> pxs = getCluster(i); 
  vector<int>::iterator px; 
  double sum=0.; 
  for (px = pxs.begin(); px!=pxs.end(); px++)
  {
       sum += pow(_image->GetBinContent(*px) - mean,2);
  }
  return sqrt(sum/((double)pxs.size()));
}

double MaxCamClusterImage::getMax(int i, int * maxBin)
{
  vector<int> pxs = getCluster(i);
  double max = 0.;
  vector<int>::iterator px; 
  for (px=pxs.begin(); px!=pxs.end(); px++)
  {
     double bin = _image->GetBinContent(*px);
     if (bin > max)
     {        
        max = bin; 
        if (maxBin!=NULL) *maxBin = *px;  
     }
  }
  return max; 
}

void MaxCamClusterImage::getXYFromBinNo(int bin, int * x, int * y, bool undobinning) const
{
  
  int width = _image->GetNbinsX() + 2; //+2 Due to overflow/underflow bins
  int height =_image->GetNbinsY() + 2;

  //Figure out the coordinates of the bin
  
  *x = bin%width;
  *y = ((bin - *x)/width)%height; 
  
  if (undobinning)
  {
    *x = (int) _image->GetXaxis()->GetBinCenter(*x);
    *y = (int) _image->GetYaxis()->GetBinCenter(*y);
  }
}


int MaxCamClusterImage::getNumNeighbors(int i, double thresholdInRMSUnits, int maxBin)
{
  //calculate and store image rms if hasn't already been done
  if (image_rms < 0)
  {
     double mean; //dummy
     MaxCamImageTools::meanRMSNoOutliers(_image,mean,image_rms);
  }  
  double thresh = thresholdInRMSUnits * image_rms; 

  //if maxBin not specified, calculate using getMax
  if (maxBin < 0)
  {
     getMax(i,&maxBin);
  }
 //Figure out the coordinates of the maxbin
  int bin_x, bin_y; 
  
  getXYFromBinNo(maxBin, &bin_x, &bin_y); 
 

  int count = 0;

  for (int i = bin_x -1; i <= bin_x + 1; i++)
  {
      for (int j = bin_y - 1; j <= bin_y + 1; j++)
      {
          if (i==bin_x && j==bin_y) continue; 
          count += (_image->GetBinContent(i,j) > thresh ? 1 : 0); 
      }
  }
  return count; 
}


void MaxCamClusterImage::clusterBounds(int i, int * xmin, int * xmax, int * ymin, 
                                       int * ymax, int border_px) const 
{
  int minx, maxx, miny, maxy; 
  int x,y;
  maxx = 0;
  maxy = 0; 
  minx = 1<<30;
  miny = 1<<30;  
  
  for (unsigned int j = 0; j < getCluster(i).size(); j++)
  {
    int bin = getCluster(i)[j]; 
    getXYFromBinNo(bin,&x,&y,true);
     
    minx = x < minx ? x : minx;
    miny = y < miny ? y : miny;
    maxx = x > maxx ? x : maxx;
    maxy = y > maxy ? y : maxy;
    
    //cout << x << " " << y << endl;
  }
  
  minx-=border_px; 
  miny-=border_px;
  maxx+=border_px;
  maxy+=border_px; 
  
  *xmin = minx; 
  *xmax = maxx; 
  *ymin = miny; 
  *ymax = maxy; 
}


void MaxCamClusterImage::drawRegion(int minx, int maxx, int miny, int maxy)
{

  TH2 * copy = (TH2 *) _image->Clone();
  copy->SetName("drawClusterClone");
  copy->SetAxisRange(minx,maxx,"X");
  copy->SetAxisRange(miny,maxy,"Y"); 
  copy->DrawCopy("colz"); 
  gROOT->Delete("drawClusterClone"); 

}
void MaxCamClusterImage::drawCluster(int i, int border_px)
{
  
  int minx,maxx,miny,maxy; 
  clusterBounds(i,&minx,&maxx,&miny,&maxy,border_px); 
  drawRegion(minx,maxx,miny,maxy);   
}

void MaxCamClusterImage::getMinMaxPosition(int i, 
	double phi, double &minVal, double &maxVal) const
{
  vector<int>::const_iterator iter =_pixels.at(i).begin();
  for (; iter < _pixels.at(i).end(); iter++){
    int xbin = (*iter) % _nbinsx;
    int ybin = (*iter) / _nbinsx;
    //getXYFromBinNo(*iter,&xbin,&ybin,false);
    double x = _image->GetXaxis()->GetBinCenter(xbin);
    double y = _image->GetYaxis()->GetBinCenter(ybin);
    double dist = x*cos(phi)+y*sin(phi);
    if (iter == _pixels.at(i).begin()){
      minVal = dist;
      maxVal = dist;
    }
    else{
      minVal = TMath::Min(minVal,dist);
      maxVal = TMath::Max(maxVal,dist);
    }
  }
}

double MaxCamClusterImage::getMoment(int i, int n, double phi, 
		int binning, TString binType)
{
  vector<int> pix = _pixels.at(i);
  vector<int>::iterator iter = pix.begin();
  double ave =0, moment = 0, total = 0;
  if (binning == 0){
    vector<double> dist;
    vector<double> counts;
    for (; iter < pix.end(); iter++){
      int xbin = (*iter) % _nbinsx;
      int ybin = (*iter) / _nbinsx;
      double x = _image->GetXaxis()->GetBinCenter(xbin);
      double y = _image->GetYaxis()->GetBinCenter(ybin);
      double tempCount = _image->GetBinContent(*iter);
      total += tempCount;
      counts.push_back(tempCount);
      dist.push_back(x*cos(phi)+y*sin(phi));
      ave += counts[counts.size()-1]*dist[dist.size()-1];
    }
    ave /= total;
    if (n != 1){
      for (unsigned int ii = 0; ii < dist.size(); ii++){
        moment += counts[ii]*pow(dist[ii]-ave,n);
      }
      moment /= total;
    }else moment = ave;
  }else{//project onto an axis
/*
    double min = 0, max = 0;
    getMinMaxPosition(i,phi,min,max);
    TH1F hist("projHist","projection onto phi", binning,min,max);
    for (iter = pix.begin(); iter < pix.end(); iter++){
      int xbin = (*iter) % _nbinsx;
      int ybin = (*iter) / _nbinsx;
      double binCount = _image->GetBinContent(*iter);
      double x = _image->GetXaxis()->GetBinCenter(xbin);
      double y = _image->GetYaxis()->GetBinCenter(ybin);
      double dist = x*cos(phi)+y*sin(phi);
      if (dist > max || dist < min) cout <<"Warning: Pixel out of range"<<endl;
      hist.Fill(dist,binCount);


    }
    cout << min <<" "<< max<<endl;
*/
    TH1* hist = projectCluster(i,phi,binning,binType);
    ave = 0;
    for (int ii=1; ii<=hist->GetNbinsX(); ii++){
      ave+= hist->GetBinContent(ii)*hist->GetXaxis()->GetBinCenter(ii);
    }
    ave/=hist->Integral();

    if (n!=1){
      total = hist->Integral();
      for (int i = 1; i <= hist->GetNbinsX(); i++){
        moment += hist->GetBinContent(i)*pow(hist->GetXaxis()->GetBinCenter(i)-ave,n);
      }
      moment /= total;
    }else moment = ave;    
    delete hist;
  }

  return moment;
}

double
MaxCamClusterImage::getMoment2(int i, int n, double phi, int binning, char* binType)
{
  if (binning <=0) {
    cout <<"Warning: no bins used, returning 0"<<endl;
    return 0;
  }
  //Similar to getMoment() but uses arrays instead of histograms.  Should be better for memory management
  vector<int> pix = _pixels.at(i);
  vector<int>::iterator iter = pix.begin();
  double min = 0, max = 0;
  getMinMaxPosition(i,phi,min,max);
  int tempSize;
  if ((TString)binType == "pixelPerBin"){
    max += 0.5*binning;
    min -= 0.5*binning;
    tempSize = static_cast<int>((max-min)/binning);
    //tempSize = static_cast<int>((max - min) / binning);
    //max = min+binning*tempSize;
    //max += 0.5 * binning;
    //min -= 0.5 * binning;
  }else {
    tempSize = binning;
    double binwidth = (max - min ) / binning;
    max += 0.5*binwidth*binning / (binning - 1);
    min -= 0.5*binwidth*binning / (binning - 1);
  }
  const int size = tempSize;
  double bins[size];
  double binsize = (max-min)/size;

  for (int i = 0; i < size; i++) bins[i] = 0;
  double total = 0;
  //cout <<"Limits: "<<min<<" "<<max<<endl;
  for (; iter < pix.end(); iter++){
    int xbin = (*iter) % _nbinsx;
    int ybin = (*iter) / _nbinsx;
    double count = _image->GetBinContent(*iter);
    double x = _image->GetXaxis()->GetBinCenter(xbin);
    double y = _image->GetYaxis()->GetBinCenter(ybin);
    double dist = x*cos(phi)+y*sin(phi);
    if (dist < min || dist >= max)
      cout <<"Warning: bin not in histogram "<<dist<<endl;
    int bin = static_cast<int>((dist-min)/binsize);
    bins[bin]+= count;
    total += count;
  }
  double mean = 0, moment = 0;
  for (int i = 0; i < size; i++){
    mean += bins[i]/total*(min+(i+0.5)*binsize);
  }
  for (int i =0; i < size; i++){
    moment += bins[i]/total*pow(min+(i+0.5)*binsize - mean,n);
  }
  return moment;
}

void
MaxCamClusterImage::getMoments(int i,double phi,int n, int* moments, double* values,int binning, char* opt)
{
  vector<int> pix;
  TString option(opt);
  option.ToLower();
  bool useClustRed = option.Contains("r");
  bool pixPerBin = !option.Contains("t");
  if (useClustRed) pix = _pixelsred.at(i);
  else pix = _pixels.at(i);

  if (pix.size()==0){
    for (int ii = 0; ii < n; ii++) values[ii] = 0;
    return;
  }
  vector<int>::iterator iter = pix.begin();
  double min = 0, max = 0;
  getMinMaxPosition(i,phi,min,max);  int tempSize;
  if (pixPerBin){
    max += 0.5*binning;
    min -= 0.5*binning;
    tempSize = static_cast<int>((max-min)/binning);
    //tempSize = static_cast<int>((max - min) / binning);
    //max = min+binning*tempSize;
    //max += 0.5 * binning;
    //min -= 0.5 * binning;
  }else {
    tempSize = binning;
    double binwidth = (max - min ) / binning;
    max += 0.5*binwidth*binning / (binning - 1);
    min -= 0.5*binwidth*binning / (binning - 1);
  }
  const int size = tempSize;
  double bins[size];
  double binsize = (max-min)/size;
  for (int ii = 0; ii < size; ii++) bins[ii] = 0;
  double total = 0;
  //cout <<"Limits: "<<min<<" "<<max<<endl;
  for (; iter < pix.end(); iter++){
    int xbin = (*iter) % _nbinsx;    
    int ybin = (*iter) / _nbinsx;
    double count = _image->GetBinContent(*iter);
    double x = _image->GetXaxis()->GetBinCenter(xbin);
    double y = _image->GetYaxis()->GetBinCenter(ybin);
    double dist = x*cos(phi)+y*sin(phi);
    if (dist < min || dist >= max)
      cout <<"Warning: bin not in histogram "<<dist<<endl;
    int bin = static_cast<int>((dist-min)/binsize);
    bins[bin]+= count;
    total += count;
  }
  double mean = 0;
  for (int ii = 0; ii < size; ii++){
    mean += bins[ii]/total*(min+(ii+0.5)*binsize);
  }
  for (int ii = 0; ii < n; ii++){
    values[ii] = 0;
    if (moments[ii] == 1) {
      values[ii] = mean;
      continue;
    }
    for (int j = 0; j < size; j++){
      values[ii] += bins[j] / total * pow(min+(j+0.5)*binsize - mean,moments[ii]);
    }

  }
}

void
MaxCamClusterImage::getMomentsUnbinned(int i,double phi, int n, int* moments, double* values,bool useClustRed)
{
  vector<int> pix;
  if (useClustRed) pix = _pixelsred.at(i);
  else pix = _pixels.at(i);
  if (pix.size() == 0){
    for (int ii = 0; ii < n; ii++) values[ii] = 0;
    return;
  }
  vector<int>::iterator iter = pix.begin();
  double min = 0, max = 0;
  getMinMaxPosition(i,phi,min,max);  
  double total = 0;
  //cout <<"Limits: "<<min<<" "<<max<<endl;
  double mean = 0;
  for (; iter < pix.end(); iter++){
    int xbin = (*iter) % _nbinsx;
    int ybin = (*iter) / _nbinsx;
    double count = _image->GetBinContent(*iter);
    double x = _image->GetXaxis()->GetBinCenter(xbin);
    double y = _image->GetYaxis()->GetBinCenter(ybin);
    double dist = x*cos(phi)+y*sin(phi); 
    mean+= count*dist;
    total += count;
  }
  mean/= total;
  for (int ii = 0; ii < n; ii++){
    values[ii] = 0;
    if (moments[ii] == 1) {
      values[ii] = mean;
      continue;
    }
    for ( iter = pix.begin(); iter < pix.end(); iter++){
      double count = _image->GetBinContent(*iter);
      int xbin = (*iter) % _nbinsx;
      int ybin = (*iter) / _nbinsx;
      double x = _image->GetXaxis()->GetBinCenter(xbin);
      double y = _image->GetYaxis()->GetBinCenter(ybin);
      double dist = x*cos(phi)+y*sin(phi);
      values[ii] += count / total * pow(dist - mean,moments[ii]);
    }
  }
}

void MaxCamClusterImage::getRayleigh(int i, double& x, double& y, char* opt)
{
//Calculate a vector Rayleigh statistic, defined as an averaged direction of pixels in the track with respect to some center point
  TString option(opt);
  option.ToLower();
  bool verbose = option.Contains("v");
  bool unweighted = option.Contains("u");
  bool reduced = option.Contains("r");
  char centerOpt = 'a';//Unweighted average bin
  if (option.Contains("b")) centerOpt='b';//Average of min/max positions
  else if (option.Contains("c")) centerOpt='c';//Weighted Average bin

  //Get set of pixels
  vector<int> pix;
  if (reduced) {
    pix = _pixelsred.at(i);
    if (verbose) cout <<"Using reduced cluster"<<endl;
  }
  else{
    pix = _pixels.at(i);
    if (verbose) cout <<"Using full cluster"<<endl;
  }
  if (verbose) cout <<"Number of pixels: "<<pix.size()<<endl;
  //***
  //Find center position
  //***
  int size = pix.size();
  double centerX=0, centerY=0;
  if (centerOpt == 'b'){
  //Take the center to be the average of the cluster's min and max values along the two axes;
    if (verbose) cout <<"Center option: Mean of Min/Max x/y values" <<endl;
    double xmin =0, ymin = 0, xmax = 0, ymax = 0;
    getMinMaxPosition(i,0,xmin,xmax);
    getMinMaxPosition(i,TMath::Pi()/2, ymin, ymax);
    centerX = 0.5*(xmin+xmax);
    centerY = 0.5*(ymin+ymax);
  }else if (centerOpt =='c'){
    if (verbose) cout <<"Center option: average point in cluster weighted by bin content"<<endl;
    centerX = 0, centerY = 0;
    double count = 0;
  //Get a weighted average position;
    for (vector<int>::iterator it=pix.begin(); it < pix.end(); it++){
      double weight = _image->GetBinContent(*it);
      int bin_x = (*it) % _nbinsx;
      int bin_y = (*it) / _nbinsx;
      count += weight;
      centerX += weight*_image->GetXaxis()->GetBinCenter(bin_x);
      centerY += weight*_image->GetYaxis()->GetBinCenter(bin_y);
    }//Loop over pixels
    centerX /= count;
    centerY /= count;
  }else {
    if (verbose) cout <<"Center option: unweighted average cluster position"<<endl;
  //Get unweighted average position;
    for (vector<int>::iterator it=pix.begin(); it < pix.end(); it++){
      int bin_x = (*it) % _nbinsx;
      int bin_y = (*it) / _nbinsx;     
      centerX += _image->GetXaxis()->GetBinCenter(bin_x);
      centerY += _image->GetYaxis()->GetBinCenter(bin_y);
    }//Loop over pixels
    centerX /= size;
    centerY /= size;
  }//Center options
  if (verbose) cout <<"Center: (x,y) = (" <<centerX<<", "<<centerY<<")"<<endl;
  //***
  //Loop over pixels and find a type of vector Rayleigh statistic
  //***
  x = 0, y = 0;
  if (unweighted){//Each bin in the cluster gets the same weight
  //This really looks for anisotropies in the cluster shape
    if (verbose) cout <<"Using unweighted bins"<<endl;
    double count = 0;
    for (vector<int>::iterator it = pix.begin(); it < pix.end(); it++){
      int bin_x = (*it) % _nbinsx;
      int bin_y = (*it) / _nbinsx;
      double dx = _image->GetXaxis()->GetBinCenter(bin_x) - centerX;
      double dy = _image->GetYaxis()->GetBinCenter(bin_y) - centerY;
      double dist = sqrt(dx*dx+dy*dy);
      if (dist==0) continue;
      count++;
      dx /= dist;
      dy /= dist;
      x += dx;
      y += dy;
    }
    if (count >0){
      x /= count;
      y /= count;
    }
  }else{//Weight each bin by the number of counts
    double count = 0;
    if (verbose) cout <<"Using weighted bins"<<endl;
    for (vector<int>::iterator it = pix.begin(); it < pix.end(); it++){
      double weight = _image->GetBinContent(*it);
      if (weight<= 0) continue;
      weight=pow(weight,2);
      int bin_x = (*it) % _nbinsx;

      int bin_y = (*it) / _nbinsx;
      double dx = _image->GetXaxis()->GetBinCenter(bin_x) - centerX;
      double dy = _image->GetYaxis()->GetBinCenter(bin_y) - centerY;
      double dist = sqrt(dx*dx+dy*dy);
      if (dist==0) continue;
      count+= weight;
      dx /= dist;
      dy /= dist;
      x += weight*dx;
      y += weight*dy;
    }
    if (count > 0){
      x /= count;
      y /= count;
    }
  }
  if (verbose) {
    cout <<"X: "<<x<<" Y: "<<y<<endl;
    cout <<"Angle: "<<180/TMath::Pi()*atan2(y,x)<<endl;
  }

}


TH1* MaxCamClusterImage::projectCluster(int i, double phi,
	int binning, TString binType) const
{
  double minVal =0, maxVal = 0;
  getMinMaxPosition(i,phi,minVal,maxVal);
  int nbins;
//Set bins by either setting total bins or pixels per bin in projection
//Endpoints of cluster are set to the center of the end bins.
  if (binType == "pixelPerBin") {
    //nbins = static_cast<int>((maxVal - minVal) / binning);
    //maxVal = minVal+binning*tempSize;
    //maxVal += 0.5 * binning;
    //minVal -= 0.5 * binning;
    maxVal += 0.5*binning;
    minVal -= 0.5*binning;
    nbins = static_cast<int>((maxVal-minVal)/binning);
  }else {
    nbins = binning;
    double binsize = (maxVal - minVal) / binning;
    minVal -= 0.5*binsize * binning / (binning - 1.);
    maxVal += 0.5*binsize * binning / (binning - 1.);
  }

  TH1F *proj = new TH1F("projHist","projection",nbins,minVal,maxVal);
  for(vector<int>::const_iterator iter = _pixels.at(i).begin(); iter < _pixels.at(i).end(); iter++){
    int xbin = (*iter) % _nbinsx;
    int ybin = (*iter) / _nbinsy;

    double x = _image->GetXaxis()->GetBinCenter(xbin);
    double y = _image->GetYaxis()->GetBinCenter(ybin);
    proj->Fill(x*cos(phi)+y*sin(phi),_image->GetBinContent(*iter));  
  }
  return proj;
} 

void MaxCamClusterImage::mergeTracks(unsigned int i,unsigned int j)
{
  if (i == j || i >= _pixels.size() || j >= _pixels.size()){
    cout <<"Invalid track numbers"<<endl;
    return;
  }
  _pixels[i].insert(_pixels[i].end(),_pixels[j].begin(),_pixels[j].end());
  _pixelsred[i].insert(_pixelsred[i].end(),_pixelsred[j].begin(),_pixelsred[j].end());
  _pixels.erase(_pixels.begin()+j);
  _pixelsred.erase(_pixelsred.begin()+j);
}

double MaxCamClusterImage::getIntegral2(int i, int nbin)
{
  vector<int> pixel = _pixels.at(i);
  vector<int>::iterator iter = pixel.begin();
  int xmax = (*iter) % _nbinsx;
  int xmin = xmax;
  int ymax = (*iter) / _nbinsx;
  int ymin = ymax;
  for (iter = pixel.begin()+1; iter<pixel.end();iter++){
    xmax = max(xmax,(*iter) % _nbinsx);
    xmin = min(xmin,(*iter) % _nbinsx);
    ymax = max(ymax,(*iter) / _nbinsx);
    ymin = min(ymin,(*iter) / _nbinsx);
  }

  xmax = min(xmax+nbin,_nbinsx-2);
  xmin = max(xmin-nbin,1);
  ymax = min(ymax+nbin,_nbins-2);
  ymin = max(ymin-nbin,1);

  double energy = 0;
  for (int i = xmin; i <= xmax; i++){
    for (int j = ymin; j <= ymax; j++){
     for (iter = pixel.begin(); iter<pixel.end(); iter++){
       int xx = (*iter) % _nbinsx;
       int yy = (*iter) / _nbinsx;
       if (sqrt((i-xx)*(i-xx)+(j-yy)*(j-yy))<nbin){
         energy += _image->GetBinContent(i,j);
         break;
       }

      }
    }
  }
  return energy;
}

void MaxCamClusterImage::getRADec(double phi, double theta, TDatime * time, double lat, double lon, 
                                  double nang, CAMERA_ORIENTATION orient, 
                                  double & ra, double & dec, double & l, double  & b)
{

    //nang MUST BE ENTERED IN DEGREES to be consistent with simulation
    //nang is the angle of the y-axis relative to North with positive
    //angle in the counterclockwise direction

    double deg2rad = TMath::Pi()/180.0;

    nang = nang*deg2rad;

     //fix phi for bottom tpc
    if(orient==BOTTOM)
    {
       if(phi>0)
         phi=TMath::Pi()-phi;
       if(phi<0)
         phi=-1*TMath::Pi()+phi;
       theta=TMath::Pi()-theta;
    }
 
    //calculate the angle of the recoil with respect to North
    //in the eastward direction
    double azimuth = -1*TMath::Pi()/2-(phi+nang);
    double altitude = -1*(TMath::Pi()/2-theta);

    dec = asin(sin(altitude)*sin(lat*deg2rad)
               +cos(altitude)*cos(azimuth)*cos(lat*deg2rad));

    double hourangle = atan2(-sin(azimuth),
                            -cos(azimuth)*sin(lat*deg2rad)+
                            tan(altitude)*cos(lat*deg2rad));

    while(hourangle>TMath::Pi() || hourangle<-1*TMath::Pi())
    {
       if(hourangle>TMath::Pi())
         hourangle-=2*TMath::Pi();
       if(hourangle<-1*TMath::Pi())
         hourangle+=2*TMath::Pi();
    }

    double J2000 = _time;
    double T=J2000/36525.0;
    double T0=6.697374558+(2400.051337*T)+(0.000025862*T*T);
    while(T0<0 || T0>=24)
    {
       if(T0<0)
         T0=T0+24;
       if(T0>=24)
         T0=T0-24;
    }
    double sidtime = 1.002737909*(time->GetHour()+(time->GetMinute()*1/60.)
                                 +(time->GetSecond()/3600.))+T0;

    while(sidtime<0 || sidtime>=24)
    {
       if(sidtime<0)
         sidtime=sidtime+24;
       if(sidtime>=24)
         sidtime=sidtime-24;
    }

    double localdifftime = lon/15.0;
    double localsidtime = localdifftime+sidtime;

    while(localsidtime<0 || localsidtime>=24)
    {
      if(localsidtime<0)
        localsidtime=localsidtime+24;
      if(localsidtime>=24)
        localsidtime=localsidtime-24;
    }

    ra=15*localsidtime-hourangle/deg2rad;
    while(ra<=-180 || ra>180)
    {
       if(ra<=-180)
         ra+=360;
       if(ra>180)
         ra-=360;
    }

    b=asin(cos(dec)*cos(27.4*deg2rad)*cos((ra-192.25)*deg2rad) +sin(dec)*sin(27.4*deg2rad));
    l=atan2(sin(dec)-sin(b)*sin(27.4*deg2rad),cos(dec)*sin((ra-192.25)*deg2rad)*cos(27.4*deg2rad))+33*deg2rad;
    b=b/deg2rad;
    l=l/deg2rad;
    dec=dec/deg2rad;

}


TH2* MaxCamClusterImage::getClusterHist(int i, int padding,  bool red, const DmtpcGainMap * map, double mingain, const char * name, bool setzero) const
{

  int minbinx =  INT_MAX; 
  int maxbinx =  -INT_MAX; 
  int minbiny =  INT_MAX; 
  int maxbiny =  -INT_MAX; 
  double xmax = -DBL_MAX, ymax = -DBL_MAX, xmin = DBL_MAX,ymin = DBL_MAX; 

  int binx,biny; 
  double x,y; 
  vector<int> pix = red ? _pixelsred[i] : _pixels[i]; 
  for (vector<int>::const_iterator it = pix.begin(); it !=pix.end(); it++)
  {
      getXYFromBinNo(*it,&binx,&biny);     
      x = _image->GetXaxis()->GetBinCenter(binx); 
      y = _image->GetYaxis()->GetBinCenter(biny); 

      if (binx < minbinx) minbinx = binx; 
      if (binx > maxbinx) maxbinx = binx; 
      if (biny < minbiny) minbiny = biny; 
      if (biny > maxbiny) maxbiny = biny; 
      if (x < xmin) xmin = x; 
      if (x > xmax) xmax = x; 
      if (y < ymin) ymin = y; 
      if (y > ymax) ymax = y; 
  }
 
  TString n (name ? name : TString::Format("%s_clust%i",_image->GetName(),i).Data()); 

  int nbinsx = maxbinx - minbinx + 1 + 2 * padding; 
  int nbinsy = maxbiny - minbiny + 1 + 2 * padding; 

  xmin -= padding * _image->GetXaxis()->GetBinWidth(1); 
  xmax += padding * _image->GetXaxis()->GetBinWidth(1); 
  ymin -= padding * _image->GetYaxis()->GetBinWidth(1); 
  ymax += padding * _image->GetYaxis()->GetBinWidth(1); 

  TH2 * ret = DmtpcRootTools::newTH2StealType(_image, n,n,nbinsx, xmin,xmax, nbinsy, ymin,ymax); 

  for (int ix = 0; ix < nbinsx; ix++)
  {
    for (int iy = 0; iy < nbinsy; iy++)
    {
      if (!setzero || isInCluster(i,_image->GetBin(ix+minbinx,iy+minbiny), red))
      {
        double g = map ?  map->getGainMap()->GetBinContent(ix + minbinx, iy + minbiny) : 1; 
        if (ix+minbinx == 0 || ix+minbinx > _image->GetNbinsX()) continue; 
        if (iy+minbiny == 0 || iy+minbiny > _image->GetNbinsY()) continue; 
        if (g > mingain) 
          ret->SetBinContent(padding +ix+1,padding +iy+1, _image->GetBinContent(ix+minbinx, iy+minbiny)/g); 
      }
    }
  }

  return ret; 
}


TH1* MaxCamClusterImage::getLeftProj(int i, int border, bool ignoreClusters, double outlier_pct) const
{

  int xmin,xmax,ymin,ymax; 
  clusterBounds(i, &xmin,&xmax,&ymin,&ymax, border); 
  int minbiny = getImage()->GetYaxis()->FindBin(ymin); 
  int maxbiny = getImage()->GetYaxis()->FindBin(ymax); 


  TH1 * leftproj = new TH1F(TString::Format("leftproj%d",i), TString::Format("Left Projection of Cluster %d",i), 
                            maxbiny - minbiny + 1, ymin, ymax); 

  std::set<int> pix; 
  MaxCamImageTools::selectPixelsInsideBox(getImage(), &pix, 0,ymin,xmin,ymax); 


  if (ignoreClusters)
  {
    std::set<int> remove_us; 
    for (std::set<int>::iterator it = pix.begin(); it!= pix.end(); it++)
    {
      for (int c = 0; c < getNCluster(); c++)
      {
        if (isInCluster(c,*it))
        {
          remove_us.insert(*it); 
        }
      }
    }

    for (std::set<int>::iterator it = remove_us.begin(); it!= remove_us.end(); it++)
    {
      pix.erase(*it); 
    }
  }

  std::vector<double> vecpix(pix.size()); 
  int j = 0; 
  for (std::set<int>::iterator it = pix.begin(); it!= pix.end(); it++)
  {
    vecpix[j++] = getImage()->GetBinContent(*it); 
  }

  int lower_bound = int(vecpix.size() * outlier_pct); 
  int upper_bound = vecpix.size() - lower_bound; 
  std::nth_element(vecpix.begin(), vecpix.begin() + lower_bound, vecpix.end()); 
  double lower = vecpix[lower_bound];
  std::nth_element(vecpix.begin(), vecpix.begin() + upper_bound, vecpix.end()); 
  double upper = vecpix[upper_bound]; 

  std::cout << lower << " , " << upper << std::endl; 
  for (std::set<int>::iterator it = pix.begin(); it != pix.end(); it++)
  {
    double val = getImage()->GetBinContent(*it); 
    if (val < lower) continue; 
    if (val > upper) continue; 
    int bx, by, bz; 
    getImage()->GetBinXYZ(*it, bx,by,bz); 
    leftproj->SetBinContent(by - minbiny + 1, leftproj->GetBinContent(by-minbiny + 1) + val); 
  }


  return leftproj; 
}

TH1 * MaxCamClusterImage::projectClusterInterpolate(int i, double phi, const char * interpolation , const DmtpcGainMap * gainmap, double mingain, bool reset_xaxis, const char * name , double * startpos ) const
{

  //pad with 2 should be safe for rotation, I think
  TH2 * temp = getClusterHist(i,2,false,gainmap, mingain); 
  TH2 * err_temp = NEW_HIST2D_WITH_SAME_SIZE(temp,TH2F,"err_temp");

  double mean, readnoise;
  meanRmsNoClustersNoHotSingle(&mean,&readnoise); 
//  std::cout << readnoise << std::endl;


  for (int ix = 1; ix <= err_temp->GetNbinsX(); ix++)
  {
    for (int iy = 1; iy <= err_temp->GetNbinsY(); iy++)
    {
      
      //Start with readnoise
      double err2 = readnoise*readnoise; 
      //Add estimate  of poisson noise (WARNING THIS IS AN OVERESTIMATE)
      err2 += temp->GetBinContent(ix,iy) > mean ? temp->GetBinContent(ix,iy) - mean : 0 ; 
      err_temp->SetBinContent(ix,iy,err2); 
    }
  }
    




  TH2 * rotated = MaxCamImageTools::rotateInterpolate(temp, phi,-1,-1,0,interpolation); 
  TH2 * rotated_err = MaxCamImageTools::rotateInterpolate(err_temp, phi,-1,-1,0,interpolation); 

  int nbinsx = rotated->GetNbinsX(); 
  int nbinsy = rotated->GetNbinsY(); 

  double xmin=0,xmax=0; 
  int minbin = 1; 
  int maxbin = nbinsx; 

  //find minimum bins... 
  bool done = false; 
  for (int ix = 1; ix <= nbinsx; ix++)
  {
    for (int iy = 1; iy <= nbinsy; iy++)
    {
      if (rotated->GetBinContent(ix,iy) > 0)
      {
        minbin = ix;  
        xmin = rotated->GetXaxis()->GetBinLowEdge(ix); 
        done = true; 
        break; 
      }
    }
    if (done) break; 
  }

  done = false; 

  for (int ix = nbinsx; ix >= 1; ix--)
  {
    for (int iy = 1; iy <= nbinsy; iy++)
    {
      if (rotated->GetBinContent(ix,iy) > 0)
      {
        maxbin = ix;  
        xmax = rotated->GetXaxis()->GetBinLowEdge(ix) + rotated->GetXaxis()->GetBinWidth(ix); 
        done = true; 
        break; 
      }
    }
    if (done) break; 
  }
 



  TString n(name ? name : TString::Format("%s_proj%d_%f",_image->GetName(),i,phi).Data()); 

  TH1F * projection = new TH1F(n,n, maxbin-minbin+1, reset_xaxis ? 0 : xmin, reset_xaxis ? xmax-xmin : xmax); 

//  std::cout << projection->GetNbinsX() << std::endl; 

  for (int ix = minbin; ix <= maxbin; ix++)
  {
    for (int iy = 1; iy <= nbinsy; iy++)
    {
        projection->SetBinContent(ix - minbin + 1,projection->GetBinContent(ix - minbin + 1) + rotated->GetBinContent(ix,iy)); 
        projection->SetBinError(ix - minbin + 1,projection->GetBinError(ix - minbin + 1)+ rotated_err->GetBinContent(ix,iy)); //Fill error with number of things projected
    }
  }

  for (int ix = 1; ix <= projection->GetNbinsX(); ix++)
  {
    projection->SetBinError(ix, sqrt(projection->GetBinError(ix))); 
  }


  temp->Delete(); 
  err_temp->Delete(); 
  rotated->Delete(); 
  rotated_err->Delete(); 

  return projection; 
}


std::list<TGraph*> MaxCamClusterImage::getClusterBoundary(int track, int color, int linewidth, bool draw, bool red) const
{
  
  
  vector<int> pix = red ? _pixelsred[track] : _pixels[track]; 
  std::list<TGraph*> edges; 
  for (unsigned int i = 0; i < pix.size(); i++)
  {
    int bin = pix[i]; 
    //Loop over directions
    for (int dir = 0; dir < 4; dir++)
    {
      int bindir = 0; 
      if (dir == 0)
      {
       bindir = bin-1; 
      }
      else if (dir ==1)
     {
       bindir = bin - _nbinsx; 
      }
      else if (dir == 2)
      {
       bindir = bin+1; 
      }
      else if (dir ==3)
      {
       bindir = bin + _nbinsx; 
      }
 
      if (!isInCluster(track,bindir, red))
      {
        TGraph * edge = new TGraph(2); 
        edge->SetEditable(0); 
        edge->SetBit(kCannotPick); 
        int bx,by; 
        bx = bin % _nbinsx; 
        by = bin / _nbinsy; 
        double x0 = _image->GetXaxis()->GetBinLowEdge(bx); 
        double x1 = _image->GetXaxis()->GetBinLowEdge(bx+1); 
        double y0 = _image->GetYaxis()->GetBinLowEdge(by); 
        double y1 = _image->GetYaxis()->GetBinLowEdge(by+1); 
//        double x0 = xmin + (bx-1) *binsize;  
//        double x1 = xmin + (bx) *binsize;
//        double y0 = xmin + (by-1) *binsize;  
//        double y1 = xmin + (by) *binsize; 
        if (dir == 0)
        {
          edge->SetPoint(0, x0,y0); 
          edge->SetPoint(1, x0,y1); 
        }
        else if (dir == 1)
        {
          edge->SetPoint(0, x0,y0); 
          edge->SetPoint(1, x1,y0); 
        }
        else if (dir == 2)
        {
          edge->SetPoint(0, x1,y0); 
          edge->SetPoint(1, x1,y1); 
        }
        else if (dir == 3)
        {
          edge->SetPoint(0, x0,y1); 
          edge->SetPoint(1, x1,y1); 
        }
 
        edge->SetLineColor(color); 
        edge->SetLineWidth(linewidth); 
 
        if (draw) edge->Draw("l same"); 
 
        edges.push_back(edge); 
      }
    }
  }
    
  return edges; 

}

int MaxCamClusterImage::getClusterPerimeter(int track,  bool red) const
{
  vector<int> pix = red ? _pixelsred[track] : _pixels[track]; 
  int edges=0;
  for (unsigned int i = 0; i < pix.size(); i++)
  {
     int bin = pix[i]; 
    //Loop over directions
    for (int dir = 0; dir < 4; dir++)
    {
      int bindir = 0; 
      if (dir == 0)
      {
       bindir = bin-1; 
      }
      else if (dir ==1)
     {
       bindir = bin - _nbinsx; 
      }
      else if (dir == 2)
      {
       bindir = bin+1; 
      }
      else if (dir ==3)
      {
       bindir = bin + _nbinsx; 
      }
 
      if (!isInCluster(track,bindir, red))
      {
        edges++;
      }
    }
  }
    
  return edges; 

}

void MaxCamClusterImage::meanRmsNoClustersNoHotSingle(double * mean, double * rms, double factor) const
{     

  if (!_image) 
  {
    cerr << "Error in meanRmsNoClustersNoHotSingle(): No image. " <<endl; 
    return; 
  }

  double sum = 0; 
  double sum2 = 0; 
  int n = 0; 

  for (int i =1; i <= _image->GetNbinsX(); i++)
  {
    for (int j = 1; j <= _image->GetNbinsY(); j++)
    {
       
      bool breakflag = false;
      for (int c = 0; c < getNCluster(); c++)
      {
        if (isInCluster(c,_image->GetBin(i,j))) 
        {
          breakflag = true; 
          break;
        }
      }

       
     if (breakflag) continue; 

      //Check to make sure not outlier!
      
      double val = _image->GetBinContent(i,j); 
      int numneigh = 0; 
      double neighsum=0; 
      for (int ii = i-1; ii <=i+1; ii++) 
      {
         if (ii == 0 || ii == _image->GetNbinsX() + 1) continue; 
         for (int jj = j-1; jj <=j+1; jj++)
         {
           if (jj == ii) continue; 
           if (jj == 0 || jj == _image->GetNbinsY() + 1) continue; 

           numneigh++; 
           neighsum += _image->GetBinContent(ii,jj); 
         }
      }

      if (val > factor * neighsum/numneigh) continue; 

      sum += val; 
      sum2 += val*val; 
//      if (n % 1000) std::cout << n << std::endl; 
      n++; 
    }
  }

  *mean = sum/n; 
  *rms = sqrt(sum*sum-sum2)/n; 

}

void MaxCamClusterImage::changeHistType(char type)
{

  TH2 * newimage = DmtpcRootTools::newTH2StealSize(_image,type, "clustimg_short", "clustimg_short"); 
  for (int i = 1; i <= _image->GetNbinsX(); i++)
  {
    for (int j = 1; j <= _image->GetNbinsY(); j++)
    {
       newimage->SetBinContent(i,j,_image->GetBinContent(i,j)); 
    }
  }

  delete _image; 
  _image = (TH2F*) newimage; 
}

bool MaxCamClusterImage::crossesCameras (int itrk, const DmtpcStitchInfo * stitch) const
{
  vector<int> clust = getCluster(itrk); 
  set<int> used_cams; 
  for (vector<int>::iterator binit = clust.begin(); binit != clust.end(); binit++) 
  {
    for (vector<TH2F*>::const_iterator fracit = stitch->frac.begin(); fracit!= stitch->frac.end(); fracit++)
    {
      if ((*fracit)->GetBinContent(*binit) > 0)
      {
        used_cams.insert(fracit - stitch->frac.begin()); 
      }

      if (used_cams.size() > 1)
      {
        return true; 
      }
    }
  }

  return false; 
}


void MaxCamClusterImage::roundValues(bool outside_clusters, double roundTo)
{

  if (roundTo <= 0) return; 

  for (int x = 1; x <= getImage()->GetNbinsX(); x++)
  {
    for (int y = 1; y <= getImage()->GetNbinsY(); y++)
    {
      int bin = getImage()->GetBin(x,y); 
      
      bool round = true; 
      if (outside_clusters)
      {
        for (int c = 0; c < getNCluster(); c++)
        {
          if (isInCluster(c,bin))
          {
            round = false; 
            break; 
          }
        }
      }
      if (!round) continue; 

      _image->SetBinContent(bin, DmtpcMath::round(getImage()->GetBinContent(bin),roundTo)); 

    }
  }

}


void MaxCamClusterImage::morphologicalOperation(int i, int nerode, int ndilate, bool red)
{
  vector <int> px = red ? getClusterRed(i) : getCluster(i);
  if (nerode || ndilate)
  {
    std::set<int> moar_pix; 
    moar_pix.insert(px.begin(), px.end()); 

    for (int i = 0; i < nerode; i++)
    {
       MaxCamImageTools::erode(_image, &moar_pix);
    }

    for (int i = 0; i < ndilate; i++)
    {
       MaxCamImageTools::dilate(_image, &moar_pix);
    }


    px.clear();
    for (std::set<int>::iterator it = moar_pix.begin(); it!= moar_pix.end(); it++)
    {
      px.push_back(*it); 
    }
  }

  (red ? _pixelsred : _pixels)[i] = px; 

}