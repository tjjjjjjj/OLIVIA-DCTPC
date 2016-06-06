
#include "TH2F.h"
#include "TMath.h"
#include "TVector3.h"
#include "TTimeStamp.h"
#include "TGraph.h"
#include <vector>
#include <cmath>
#include <iostream>
#include "MaxCamCluster.hh"
#include "McWimp.hh"

using namespace std;

MaxCamCluster::MaxCamCluster()
{
   _nbinsx=0;
   _nbinsy=0;
   _nbins=0;
}

MaxCamCluster::MaxCamCluster(vector<int> pixels, TH2* image)
{
   _pixels = pixels;
   _image = image;

   //   cout << "There are " << _pixels.size() << " entries in the cluster \n";

   _nbinsx = _image->GetNbinsX()+2;
   _nbinsy = _image->GetNbinsY()+2;
   _nbins = _nbinsx*_nbinsy;

}

MaxCamCluster::~MaxCamCluster()
{
   delete _image;
   _pixels.clear();

}

double MaxCamCluster::getIntegral()
{
   double sum=0;
   for(int i=0; i<int(_pixels.size());i++)
   {
      sum += _image->GetBinContent(_pixels[i]);
   }

   return sum;
}

double MaxCamCluster::getLength(double &x1, double &y1, double &x2, double &y2)
{
   double length = 0;
   for(int i=0; i<int(_pixels.size()); i++)
   {
      for(int j=i+1; j<int(_pixels.size());j++)
      {
	 int ibinx = _pixels[i]%_nbinsx;
	 int ibiny = ((_pixels[i]-ibinx)/_nbinsx)%_nbinsy;
	 int jbinx = _pixels[j]%_nbinsx;
	 int jbiny = ((_pixels[j]-jbinx)/_nbinsx)%_nbinsy;

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

bool MaxCamCluster::isInCluster(int bin)
{
   for(int i=0; i<int(_pixels.size()); i++)
   {  
      if(_pixels[i]==bin) return true;
   }
   return false;
}

void MaxCamCluster::getXY(double &x, double &y)
{
   int xx,yy;
   double xw,yw;

   int nx=_image->GetNbinsX()+2;
   int ny=_image->GetNbinsY()+2;

   for(int i=0; i<int(_pixels.size()); i++)
   {
      xx = _pixels[i]%nx;
      yy = ((_pixels[i]-xx)/nx)%ny;
      xw=_image->GetXaxis()->GetBinCenter(xx);
      yw=_image->GetYaxis()->GetBinCenter(yy);
      x+=_image->GetBinContent(_pixels[i])*xw;
      y+=_image->GetBinContent(_pixels[i])*yw;
   }

   x=x/getIntegral();
   y=y/getIntegral();

}

double MaxCamCluster::getPhi()
{
   double x,y,x1,x2,y1,y2,length;
   length=getLength(x1,y1,x2,y2);
   getXY(x,y);

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

double MaxCamCluster::getPhi2(){
  //returns angle in xy plane, no directionality yet

  int size = (int) _pixels.size();
  double meanx = 0, meany = 0, xsq = 0, ysq = 0, xy = 0, sum = 0;

  for (int i = 0; i < size; i++){

    int xbin = _pixels[i]%_nbinsx;
    int ybin = _pixels[i]/_nbinsx;
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

  double var1 = varx * cos(phi1)*cos(phi1) + vary * sin(phi1)*sin(phi1) + cov * sin(2.*phi1);
  double var2 = varx * cos(phi2)*cos(phi2) + vary * sin(phi2)*sin(phi2) + cov * sin(2.*phi2);

  double phi;
  if (var1 < var2) phi = phi2;
  else phi = phi1;

  return phi;

}

//Calculates phi by finding angle with maximum sigma.  Similar precision to getPhi2 but slightly less accurate.
double MaxCamCluster::getPhi3(){

  int size     = (int) _pixels.size();
  double meanx = 0, meany = 0, xsq = 0, ysq = 0, xy = 0, sum = 0;

  for (int i = 0; i < size; i++){

    int xbin       = _pixels[i]%(_nbinsx);
    int ybin       = _pixels[i]/(_nbinsx);
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

double MaxCamCluster::getLength2(double theta, double pixpermm){
  //calculates length of cluster along on axis at angle theta above x.
  int size = (int) _pixels.size();
  double minx = _image->GetXaxis()->GetBinLowEdge(1);
  double miny = _image->GetYaxis()->GetBinLowEdge(1);

  double mindist = 0;
  double maxdist = 0;

  for (int i = 0; i < size; i++){
    //save only the two extreme values
    int xbin = _pixels[i] % _nbinsx;
    int ybin = _pixels[i] / _nbinsx;
    double dist = (_image->GetXaxis()->GetBinCenter(xbin) - minx)*cos(theta)/pixpermm + 
      (_image->GetYaxis()->GetBinCenter(ybin) - miny)*sin(theta)/pixpermm;

    if (dist < mindist || i == 0) mindist = dist;
    if (dist > maxdist || i == 0) maxdist = dist;
  }
  double length = maxdist - mindist;
  return length;
}

double MaxCamCluster::getSkewness(double theta){

  int size = (int) _pixels.size();
  //x*cos(theta) + y*sin(theta)
  double mean = 0;
  double total = 0;
  for (int i = 0; i < size; i++){

    int xbin = _pixels[i] % _nbinsx;
    int ybin = _pixels[i] / _nbinsy;
    mean += _image->GetBinContent(_pixels[i])*(xbin*cos(theta)+ ybin*sin(theta));
    total += _image->GetBinContent(_pixels[i]);

  }
  mean /=total;
  double moment3 = 0;
  double moment2 = 0;

  for (int i = 0; i < size; i++){

    int xbin = _pixels[i] % _nbinsx;
    int ybin = _pixels[i] / _nbinsy;

    double dist = xbin*cos(theta) + ybin*sin(theta);
    double content = _image->GetBinContent(_pixels[i]);
    moment3 = content*(dist-mean)*(dist-mean)*(dist-mean);
    moment2 = content*(dist-mean)*(dist-mean);

  }

  moment3 /= total;
  moment2 /= total;
  double skewness = moment3 / pow(moment2,1.5);
  return skewness;

}

bool MaxCamCluster::hitsEdge(){

  bool hit = false;

  int count = 0;
  while (count < (int)_pixels.size()){
    int bin = _pixels[count];
    int xbin = bin % _nbinsx;
    int ybin = bin / _nbinsx; 

    if ((xbin<10) || (xbin>(_nbinsx-10))) {
      hit = true;
      break;
    }
    else if ((ybin<10) || (ybin>(_nbinsy-10))) {
      hit = true;
      break;
    }
    count++;
  }

return hit;
}

void MaxCamCluster::changeImage(TH2* newimage)
{
   int xbinso = _image->GetNbinsX();
   int xbinsn = newimage->GetNbinsX();
   int ybinso = _image->GetNbinsY();
   int ybinsn = newimage->GetNbinsY();

   vector<int> newpx;

   if (xbinso > xbinsn || ybinso > ybinsn)
   {cout << "Unable to change image, new image has fewer bins than old image" << endl; return;}

   if(xbinso == xbinsn && ybinso == ybinsn)
   {
      _image = newimage;
      return;
   }
   else
   {
      int x,y;
      int xfactor = (xbinsn-(xbinsn%xbinso))/xbinso;
      int yfactor = (ybinsn-(ybinsn%ybinso))/ybinso;
      for(int i=0; i<int(_pixels.size()); i++)
      {
	 x = _pixels[i]%(xbinso+2);
	 y = ((_pixels[i]-x)/(xbinso+2))%(ybinso+2);
	 for(int j=(x-1)*xfactor+1; j<=x*xfactor; j++)
	 {
	    for(int k=(y-1)*yfactor+1; k<=y*yfactor; k++)
	    {
	       int bin = newimage->GetBin(j,k);
	       newpx.push_back(bin);
	    }
	 }
      }
      _pixels = newpx;
      _image = newimage;
      _nbinsx = _image->GetNbinsX()+2;
      _nbinsy = _image->GetNbinsY()+2;
      _nbins = _nbinsx*_nbinsy;
      return;

   }

   
}

double MaxCamCluster::getCygnusAngle(double yangN, double theta)
{
   McWimp* mcwimp  = new McWimp();
   //use cambridge coordinates
   mcwimp->setLat(42.373611);
   mcwimp->setLon(-71.110556);
   //figure out time
   mcwimp->setTime(_time);
   //vertically oriented detector
   mcwimp->setOpt("lv");
   TVector3 dettocyg = mcwimp->cosmicGun();
   dettocyg.RotateZ(yangN);

   double phi = getPhi2();
//   double theta = getTheta();
   TVector3 angindet(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
   double cygang = dettocyg.Dot(angindet)/(angindet.Mag()*dettocyg.Mag());


   return cygang;
}

void MaxCamCluster::setTime(TDatime* time)
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
   _time = Julian-2451545.0;
}

double MaxCamCluster::getPhi4()
{
   const int npx = int(_pixels.size());
   double xx[npx], yy[npx];
   for(int j=0; j<npx; j++)
   {
      int xj = _pixels[j]%_nbinsx;
      int yj = ((_pixels[j]-xj)/_nbinsx)%_nbinsy;
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
