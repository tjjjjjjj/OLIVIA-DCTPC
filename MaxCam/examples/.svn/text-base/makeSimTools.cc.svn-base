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

#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
////Estimate the track length from the nonzero bins in the track image      //
//////////////////////////////////////////////////////////////////////////////
double getTrackLength(TH2F* trackImage,double phi){

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

void makeTextFile(int NEvents = 100,TString opt = "cf",TString opt2 = "F",
		  TString opt3 = "CF4",
		  TString filename = "dmtpc_run000000_0_sim.txt",
		  float gain = 3.47, float noise = 7, float wimpMass = 100e6,
		  int npix = 1024,int nbins = 256,float width = 0.143*1024,
		  float pressure = 75,float spacerdiam = 0.5,
		  float spacerspacing = 20.0,float diffconst = 0.2704,
		  float diffdz = 0.0076){

  ofstream ofstr = ofstream(filename);

  ofstr << "Today's date (YYMMDD) and run number included in filename."<<endl;
  ofstr << "Summary of Detector Parameters Given in MaxCamMC"<<endl;
  ofstr << "w: WIMP, cf: Cf-252 n, nr: random n, af: fixed alpha gun, ar: random alpha" << endl;
  ofstr << "Number of Events: "         << NEvents       << endl;
  ofstr << "Event Type: "               << opt           << endl;
  ofstr << "Recoil Particle: "          << opt2          << endl;
  ofstr << "Gas Type: "                 << opt3          << endl;
  ofstr << "Wimp Mass: "                << wimpMass      << endl;
  ofstr << "Gain (ADC / keV): "         << gain          << endl;
  ofstr << "Noise (ADC counts): "       << noise         << endl;
  ofstr << "Number of ccd pixels: "     << npix  << " by " << npix  << endl;
  ofstr << "Binning: "                  << nbins << " by " << nbins << endl;
  ofstr << "CCD roi size (mm ^ 2): "    << width << " by " << width << endl;
  ofstr << "Chamber pressure (Torr): "  << pressure      <<endl;
  ofstr << "Spacer Diameter (mm): "     << spacerdiam    <<endl;
  ofstr << "Spacer Spacing (mm): "      << spacerspacing <<endl;
  ofstr << "Diffusion constant term: "  << diffconst     <<endl;
  ofstr << "Diffusion dz term: "        << diffdz        <<endl;
  ofstr.close();

}//End makeTextFile()

TTimeStamp* makeTimestamp(double t){

  int days = (int) t / 1;
  double seconds = (t - days) * 24*3600;
  int sec = (int) seconds / 1;
  int nanosec = (int) sec - seconds;

  TTimeStamp *stamp = new TTimeStamp(1999,12,31+days,12,0,sec,nanosec);
  return stamp;

}

TDatime* makeTDatime(double t){

  TTimeStamp *stamp = makeTimestamp(t);
  unsigned int date = stamp->GetDate();
  int year = date / 10000;
  int month = (date - year*10000) / 100;
  int day = date - 10000*year - 100*month;
  unsigned int time = stamp->GetTime();
  int hour = time / 10000;
  int min = (time - hour * 10000)/ 100;
  int sec = time - hour * 10000 - min * 100;
  TDatime *datime = new TDatime(year,month,day,hour,min,sec);
  return datime;

}
//t must be greater than zero
// vector<int> convertTime(double t){

//   if (t < 0.5)
//    double time = t - 0.5;//time in days since noon, 31/12/1999
//    int year, month, day, hour, min, sec;
//    if (time > 0){
//      int Y = (t-0.5) / 365;
//      int leapdays = Y / 4 + 1;
//      int _daysInFeb = 28;
//      int isLeapyear = 0;
//      if (Y % 4 == 0 ) {
//        isLeapyear = 1;
//        _daysInFeb = 29;
//      }
//      year = 2000+Y;
//      //Day of year = time in days - 365*full years past - # of leap days
//      int dayofyear = (int) (time - leapdays+1) %365;
// double days = time - 

//    }
//    else{
//      year = 1999;
//      month = 12;
//      day = 31;
//      double hours = t*24;
//      hour = (int) 12 + hours;
//      double mins = 60*(hours - hour + 12);
//      min = (int) mins;
//      double seconds = 60*(mins - min);
//      sec = (int) seconds;
//      nanosec = (int) 1e9 * (seconds - sec);
//    }

//    vector<int> timevect;
//    timevect.push_back(year);
//    timevect.push_back(month);
//    timevect.push_back(day);
//    timevect.push_back(hour);
//    timevect.push_back(min);
//    timevect.push_back(sec);
//    timevect.push_back(nanosec);
//    return timevect;

//  }
