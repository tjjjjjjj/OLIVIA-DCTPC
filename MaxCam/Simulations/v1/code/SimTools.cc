/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
////Includes various methods to be called by makeSimulations.cc        //
/////////////////////////////////////////////////////////////////////////
#include "TH2.h"
#include "TSystem.h"
#include "TString.h"
#include "TMath.h"
#include "TVector3.h"
#include "TTimeStamp.h"
#include "TDatime.h"

#include "SimTools.hh"

#include <vector>
#include <iostream>
#include <fstream>
using namespace std;


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
////Estimate the track length from the nonzero bins in the track image      //
//////////////////////////////////////////////////////////////////////////////
double 
SimTools::getTrackLength(TH2F* trackImage,double phi)
{
  double max = 0, min = 0;
  double x, y;
  double dist;
  int count = 0;
  for (int i = 1; i<=trackImage->GetNbinsX();i++){
    for (int j = 1; j<=trackImage->GetNbinsY();j++){

      if (trackImage->GetBinContent(i,j)>1.e-12){
	x = trackImage->GetXaxis()->GetBinCenter(i);
	y = trackImage->GetYaxis()->GetBinCenter(j);
	dist = x*cos(phi)+y*sin(phi);
	count++;
	if (count == 1){max = dist; min = dist;}
	else{
	  max = (dist > max) ? dist : max;
	  min = (dist < min) ? dist : min;
	}

      }
    }
  }

  double length = max - min;
  return length;
}//end getTrackLength()

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
////Make a text file giving the information about the run            //
///////////////////////////////////////////////////////////////////////

// void makeTextFile(int NEvents = 100,TString opt = "cf",TString opt2 = "F",
// 		  TString opt3 = "CF4",
// 		  TString filename = "dmtpc_run000000_0_sim.txt",
// 		  float gain = 3.47, float noise = 7, float wimpMass = 100e6,
// 		  int npix = 1024,int nbins = 256,float width = 0.143*1024,
// 		  float pressure = 75,float spacerdiam = 0.5,
// 		  float spacerspacing = 20.0,float diffconst = 0.2704,
// 		  float diffdz = 0.0076){

//   ofstream ofstr = ofstream(filename);

//   ofstr << "Today's date (YYMMDD) and run number included in filename."<<endl;
//   ofstr << "Summary of Detector Parameters Given in MaxCamMC"<<endl;
//   ofstr << "w: WIMP, cf: Cf-252 n, nr: random n, af: fixed alpha gun, ar: random alpha" << endl;
//   ofstr << "Number of Events: "         << NEvents       << endl;
//   ofstr << "Event Type: "               << opt           << endl;
//   ofstr << "Recoil Particle: "          << opt2          << endl;
//   ofstr << "Gas Type: "                 << opt3          << endl;
//   ofstr << "Wimp Mass: "                << wimpMass      << endl;
//   ofstr << "Gain (ADC / keV): "         << gain          << endl;
//   ofstr << "Noise (ADC counts): "       << noise         << endl;
//   ofstr << "Number of ccd pixels: "     << npix  << " by " << npix  << endl;
//   ofstr << "Binning: "                  << nbins << " by " << nbins << endl;
//   ofstr << "CCD roi size (mm ^ 2): "    << width << " by " << width << endl;
//   ofstr << "Chamber pressure (Torr): "  << pressure      <<endl;
//   ofstr << "Spacer Diameter (mm): "     << spacerdiam    <<endl;
//   ofstr << "Spacer Spacing (mm): "      << spacerspacing <<endl;
//   ofstr << "Diffusion constant term: "  << diffconst     <<endl;
//   ofstr << "Diffusion dz term: "        << diffdz        <<endl;
//   ofstr.close();

// }//End makeTextFile()

//From Lewin and Smith, numbers are valid for time n in days from
//1200 GMT/UT 31 Dec 1999
TTimeStamp* 
SimTools::makeTimeStamp(double t)
{
  int days = (int) t / 1;
  double seconds = (t - days) * 24*3600;
  int sec = (int) seconds / 1;
  int nanosec = (int) (sec - seconds);

  TTimeStamp *stamp = new TTimeStamp(1999,12,31+days,12,0,sec,nanosec);
  return stamp;
}

TDatime* 
SimTools::makeTDatime(double t)
{
  TTimeStamp *stamp = makeTimeStamp(t);
  TDatime *datime = convertStamptoDatime(stamp);
  return datime;
}
TDatime *
SimTools::convertStamptoDatime(TTimeStamp *stamp)
{
  unsigned int date = stamp->GetDate();
  int year = date / 10000;
  int month = (date - year*10000) / 100;
  int day = date - 10000*year - 100*month;
  unsigned int time = stamp->GetTime();
  int hour = time / 10000;
  int min = (time - hour * 10000)/100;
  int sec = time - 10000*hour - 100*min;
  TDatime *datime = new TDatime(year,month,day,hour,min,sec);
  return datime;
}
double 
SimTools::getTimeDouble(TTimeStamp *stamp)
{
  TDatime *datime = convertStamptoDatime(stamp);
  double time = getTimeDouble(datime);
  delete datime;
  return time;
}

double 
SimTools::getTimeDouble(TDatime *datime)
{
  //Should work for the years 1901 - 2099 (excluding any leap seconds)
  //Find time from New Year's, 2000
  //Works for time before, after 2000
  double time;
  int month[13];
  month[0] = 0;
  month[1] = 31;//Jan
  month[2] = 28;//Feb
  month[3] = 31;//Mar
  month[4] = 30;//Apr
  month[5] = 31;//May
  month[6] = 30;//June
  month[7] = 31;//July
  month[8] = 31;//Aug
  month[9] = 30;//Sept
  month[10] = 31;//Oct
  month[11] = 30;//Nov
  month[12] = 31;//Dec

  int Nyears = datime->GetYear() - 2000;
  int leapdays;
  if (Nyears >= 0){
    leapdays = (Nyears+3) / 4 ;//2000 is a leap year
   }
  else{
    leapdays = Nyears / 4;
  }
  time = Nyears*365 + leapdays;
  if (Nyears % 4 == 0) month[2] = 29;
  int months = datime->GetMonth();
  for (int i = 0; i< months;i++){
    time += month[i];
  }
  time += datime->GetDay()-1;
  time += datime->GetHour() / 24.;
  time += datime->GetMinute() / (24*60.);
  time += datime->GetSecond() / (24*60*60.);
  time += 0.5;//we measure time from noon, 31 Dec 99

  return time;
}

double
SimTools::findParticleMass(TString name)
{
  name.ToLower();
  double mass;
  if(name == "alpha" || name =="he") mass = 4001506;
  else if (name == "neutron") mass = 939566;
  else if (name == "h" || name == "proton") mass = 938272;
  else if (name == "c" || name == "carbon") mass = 12e6;
  else if (name == "f" || name == "fluorine") mass = 19e6;
  else if (name == "ar" || name == "argon") mass = 40e6;
  else if (name == "xe" || name == "xenon") mass = 132e6;
  else if (name == "po" || name == "polonium") mass = 218e6;
  else if (name == "rn" || name == "radon") mass = 222e6;
  else if (name == "u" || name == "uranium") mass = 238e6;
  else {
    cerr << "WARNING: Particle Not Found!!"<<endl;
    mass = -1;
  }
  return mass;
}

double
SimTools::findParticleA(TString name)
{
   name.ToLower();
   double A;
   if(name == "alpha" || name =="he") A = 4;
   else if (name == "neutron") A = 1;
   else if (name == "h" || name == "proton") A = 1;
   else if (name == "c" || name == "carbon") A = 12.01;
   else if (name == "f" || name == "fluorine") A = 19;
   else if (name == "ar" || name == "argon") A = 39.95;
   else if (name == "xe" || name == "xenon") A = 131.30;
   else if (name == "po" || name == "polonium") A = 210;
   else if (name == "rn" || name == "radon") A = 222;
   else if (name == "u" || name == "uranium") A = 238.03;
   else {
      cerr << "WARNING: Particle Not Found!!"<<endl;
      A = -1;
   }
   return A;
}
