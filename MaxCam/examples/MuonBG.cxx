#include "TRandom.h"
#include "TCanvas.h"

#include "TH1F.h"
#include "TF1.h"
#include "TAxis.h"

#include<iostream>


//// May be useful to access these outside scope of totalRecoils or makeHisto

TH1F recRate; // recoil rate
TH1F recIntg; //        integral
TH1F recRateX; // recoil rate after cross section
TH1F recIntgX; //        integral
TH1F neuRate; // neutron rate
TH1F neuIntg; //         integral

TCanvas c1,c2;

TF1 MuonNeutron, Xsection;

int mpb_DONT_MODIFY = 1; //  MeV per bin : change with histo titles




////////////////////////////////////////////////////////////////
////
//// totalMuonRecoils( eLow, eUp)
////
//// Make sure to run makeMuonHisto(numberEvents, siteDepth(in km.w.e.), CF4 pressure in torr ) FIRST! 
//// Example:  makeMuonHisto(5000, 3.1, 100 )  
//// 
////
//// Returns total number of recoils between two energy bounds
//// i.e. the integral of the recoil rate bewteen two energies
////
//// Input: (int) lower energy bound in MeV
////        (int) upper energy bound in MeV
////
//// Output: (float) total number of recoils per cm^2 per second
////
////////////////////////////////////////////////////////////////

float totalMuonRecoils(int eLow, int eUp) {

  int mpb = mpb_DONT_MODIFY;
  int min = (int)(eLow/mpb);
  int max = (int)(eUp/mpb);

//   double total = 0;
//   for(int i=min+1; i<=max; i++) {
//       total += recRate.GetBinContent(i)*mpb;
//   }

  double tLow, tUp;
  tLow = recIntgX.GetBinContent(min+1);
  tUp = recIntgX.GetBinContent(max+1);
  double total = tLow-tUp;

  cout << "Total recoils above " << eLow << "MeV = " << tLow << endl;
  cout << "Total recoils above " << eUp << "MeV = " << tUp << endl;
  cout << "Total recoils between " << eLow << " & " << eUp << "keV = " << total << endl;
  return total;
 
}




// double Xsection1(double x) {
//   return (pow(10,((-3-log10(x))/2*1.14+.56))-3.6);
// }

double Xsection2(double x) {
  double r = .2;
  return (3.6 + 136*exp(-pow((x-2.7e4),2)/(r*4e6)) + 47*exp(-pow((x-4.9e4),2)/(r*16e6)) + 25*exp(-pow((x-1e5),2)/(r*5.29e8)) + 8.7*exp(-pow((x-3.3e5),2)/(r*5.6e10)));
}

///////////////////////////////////////////////////
//// Makes Histograms
////
//// Input: (int) number of events to simulate
////        (double) depth of site, km.w.e.
////        (double) pressure of gas, torr
///////////////////////////////////////////////////

void makeMuonHisto(int eventN, double depth=3.1, double pressure=180) {

  double Emin   = 0;        // max and min for recoil energy range (in what unit??) 
  double Emax   = 1e2;
  double NeuMin = 10;       // max and min for neutron energy range (in what unit??) 
  double NeuMax = 3.5e3;
  double neuM   = 1e3;      // mass in MeV
  double recoilM = 19e3;


  //// Initialize all histos

  int mpb = mpb_DONT_MODIFY; // MeV per bin for recoils : change with histo titles
  int min = ((int)(Emin/mpb)) * mpb;
  int max = ((int)(Emax/mpb)) * mpb;
  int binN = (max-min)/mpb;
  int mpbN = 50; // MeV per bin for muon induced neutrons : change with histo titles
  int Nmin = ((int)(NeuMin/mpbN)) * mpbN;
  int Nmax = ((int)(NeuMax/mpbN)) * mpbN;
  int nBinN = (Nmax-Nmin)/mpbN;

  recRate = TH1F("","",binN,min,max); // former title: Event Rate per eRecoil (MeV)
  recRate.SetMinimum(0);
  recRate.GetYaxis()->SetTitle("Events (cm^-2 s^-1 (MeV)^-1)"); // change when mpb is changed
  recRate.GetXaxis()->SetTitle("Recoil Energy (MeV) (MeV bins)");      // ditto

  recIntg = TH1F("","",binN,min,max); // former title: Total Recoil Rate for thresholds (MeV)
  recIntg.SetMinimum(0);
  recIntg.GetYaxis()->SetTitle("Events over threshold (cm^-2 s^-1)");
  recIntg.GetXaxis()->SetTitle("Recoil threshold (MeV)");

  recRateX = TH1F("","",binN,min,max); // former title: Event Rate per eRecoil (MeV)
  recRateX.SetMinimum(0);
  recRateX.GetYaxis()->SetTitle("Events in m^3 (yr^-1 (MeV)^-1)"); // change when mpb is changed
  recRateX.GetXaxis()->SetTitle("Recoil Energy (MeV) (MeV bins)");      // ditto

  recIntgX = TH1F("","",binN,min,max); // former title: Total Recoil Rate for thresholds (MeV)
  recIntgX.SetMinimum(0);
  recIntgX.GetYaxis()->SetTitle("Events in m^3 over threshold (yr^-1)");
  recIntgX.GetXaxis()->SetTitle("Recoil threshold (MeV)");

  neuRate = TH1F("","",nBinN,Nmin,Nmax); // former title: Neutron Rate per energy (MeV)
  neuRate.GetYaxis()->SetTitle("Events (cm^-2 s^-1 (50MeV)^-1)"); // change when mpbN is changed
  neuRate.GetXaxis()->SetTitle("Muon Neutron Energy (MeV) (50MeV bins)");      // ditto

  neuIntg = TH1F("","",nBinN,Nmin,Nmax); // former title: Total Neutron flux for thresholds (MeV)
  neuIntg.GetYaxis()->SetTitle("Events over threshold (cm^-2 s^-1)");
  neuIntg.GetXaxis()->SetTitle("Muon Neutron threshold (MeV)");    

  // Fitting function for neutron energy distribution
  //TF1 MuonNeutron;

  // GS: please rewrite this equation so it is easy to switch between various labs 
  //   double Emean = 91.;   double a0=7.828;   // double a1=
  MuonNeutron = TF1("", "1*((exp(-7.828*x))/x + (0.324-0.641*exp(-0.014*270))*exp(-2.23*x)) + -7.505e-15 * exp(-2.831)",.01,3.5); // eq 14  Mei-Hime @GranSasso
  MuonNeutron.GetYaxis()->SetTitle("Events (cm^-2 s^-1 (50MeV)^-1)");
  MuonNeutron.GetXaxis()->SetTitle("Muon Neutron energy (keV)");  

  //TF1 Xsection;
  // Please rewrite this so one can switch betwen various elements 
  Xsection = TF1("","(x<1e-3)*(pow(10,((-3-log10(x))/2*1.14+.56))-3.6) + (x<1e6)*Xsection2(x) + ((x>=1e6) && (x<6.7e6))*(pow(10,(.56-(log10(x)-6)/.82*.32))) + (x>=6.7e6)*1.7",1e-5,3.5e9);
  Xsection.GetYaxis()->SetTitle("F n cross section (barn)");
  Xsection.GetXaxis()->SetTitle("Neutron energy (eV)");  

//   Xsection1 = TF1("","pow(10,((-3-log10(x))/2*1.14+.56))-3.6",1e-5,1e-3);
//   Xsection2 = TF1("","0 + 136*exp(-pow((x-2.7e4),2)/(2*4e6)) + 0*47*exp(-pow((x-4.9e4),2)/(2*16e6)) + 0*25*exp(-pow((x-1e5),2)/(5.29e8)) + 0*8.7*exp(-pow((x-3.3e5),2)/(2*6.2e10))",1e-5,1e6);
//   Xsection3 = TF1("","3.6-pow(10,((log10(x)-6)/.82*.32))",1e6,6.7e6);
//   Xsection4 = TF1("","1.7",6.7e6,3.5e9);


  double recoilE, Xrate, neuE;
  double rhoF = 2.4e25*1.0e-6/180.0*pressure;
  double SecYr = 365.25*24*3600;
  double vol = 1e6; // volume of m^3 detector, in cm^3
  double barn2cm = 1e-28;

  cout << "Events Simulated:" << endl;

  for (double n=1; n<=eventN; n++) {

    int nn = n;
    if (nn%5000==0) {cout << nn << "  " << flush;}
    


    neuE = MuonNeutron.GetRandom()*1e3;  // E neutron in MeV sampled from input distribution 



    //// Make Recoils

    double cosRecoil;

    while (1) {

      recoilE = gRandom->Rndm() * binN*mpb + min;

      //// Calculate cosine of recoil relative to WIMP direction
      double mRatio = recoilM/neuM; // mass ratio: target over projectile
      cosRecoil = sqrt(recoilE/(mRatio*neuE))*(1+mRatio)*0.5;     
      
      
      //// check if kinematics is valid   
      if (cosRecoil>-1 && cosRecoil<1) break;
      
    }

    Xrate = Xsection.Eval(neuE*1e6) * barn2cm * rhoF * SecYr * vol;
    // above finds rate given xsection for neuE
    // multiplies neuE by 1e3 before eval. b/c neuE is in keV and fxn is in eV

    //// Fill histos

    double shift = pow(10,((3.1-depth)/1.8));

    recRate.Fill(recoilE,(double) (1.05e-12)*shift/eventN*(max-min));
    recRateX.Fill(recoilE,(double) (1.05e-12)*Xrate*shift/eventN*(max-min));
    neuRate.Fill(neuE,(double) (3e-13)*shift/eventN*(Nmax-Nmin));

  }

  //// Integrate histos

  double Bin;
  for(int i=1; i<=(max/mpb); i++) {
    Bin = 0;
     for(int j=i; j<=(max/mpb); j++) {
       Bin += recRate.GetBinContent(j);  // sum up bins
    }
    recIntg.SetBinContent(i,Bin);
  }

  double BinX;
  for(int i=1; i<=(max/mpb); i++) {
    BinX = 0;
     for(int j=i; j<=(max/mpb); j++) {
       BinX += recRateX.GetBinContent(j);  // sum up bins
    }
    recIntgX.SetBinContent(i,BinX);
  }

  double Dis;
  for(int i=1; i<=(Nmax/mpbN); i++) {
    Dis = 0;
    for(int j=i; j<=(Nmax/mpbN); j++) {
      Dis += neuRate.GetBinContent(j);  // sum up bins
    }
    neuIntg.SetBinContent(i,Dis);
  }
    
  
  //// Draw histos 

  c1.Clear();

  c1.Divide(2,3);
  c1.cd(1);
  neuRate.Draw();
  c1.cd(2);
  neuIntg.Draw();
  c1.cd(3);
  recRate.Draw();
  c1.cd(4);
  recIntg.Draw();
  c1.cd(5);
  recRateX.Draw();
  c1.cd(6);
  recIntgX.Draw();

  cout << endl;

  c2.Clear();
  c2.cd();
  c2.SetLogx(1);
  c2.SetLogy(1);
  Xsection.Draw();
  

}
