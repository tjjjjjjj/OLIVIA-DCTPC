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

TCanvas c1, c2;

TF1 RockNeutron, Xsection;

int kpb_DONT_MODIFY = 1e1; //  keV per bin : change with histo titles




////////////////////////////////////////////////////////////////
////
//// totalRockRecoils( eLow, eUp)
////
//// Make sure to run makeRockHisto( numberEvents) FIRST!
////
//// Returns total number of recoils between two energy bounds
//// i.e. the integral of the recoil rate bewteen two energies
////
//// Input: (int) lower energy bound in keV
////        (int) upper energy bound in keV
////
//// Output: (float) total number of recoils per cm^2 per year
////
////////////////////////////////////////////////////////////////

float totalRockRecoils(int eLow, int eUp) { // eLow and eUp in keV 

  int kpb = kpb_DONT_MODIFY;
  int min = (int)(eLow/kpb);
  int max = (int)(eUp/kpb);

  double tLow, tUp;
  tLow = recIntgX.GetBinContent(min+1);
  tUp = recIntgX.GetBinContent(max+1);
  double total = tLow-tUp;

  cout << "Total recoils above " << eLow << "keV = " << tLow << endl;
  cout << "Total recoils above " << eUp << "keV = " << tUp << endl;
  cout << "Total recoils between " << eLow << " & " << eUp << "keV = " << total << endl;
  return total;
 
}




double Xsection2(double x) {
  double r = .2;
  return (3.6 + 136*exp(-pow((x-2.7e4),2)/(r*4e6)) + 47*exp(-pow((x-4.9e4),2)/(r*16e6)) + 25*exp(-pow((x-1e5),2)/(r*5.29e8)) + 8.7*exp(-pow((x-3.3e5),2)/(r*5.6e10)));
}

///////////////////////////////////////////////////
//// Makes Histograms
////
//// Input: (int) number of events to simulate
////        (double) pressure of gas, torr
////        (int) if 1 switches from CF4 to Ge
////        (int) if 1 switches from 10L to cubic m detector
///////////////////////////////////////////////////

void makeRockHisto(int eventN, double pressure=180, int Ge=0, int m3=1) {

  double Emin = 0;    // max and min for recoil energy range
  double Emax = 1e3;
  double NeuMin = 1e2;    // max and min for neutron energy range
  double NeuMax = 4.5e3;
  double neuM = 1e6;      // mass in keV
  double recoilM = 19e6;


  //// Initialize all histos

  int kpb = kpb_DONT_MODIFY; // keV per bin for recoils : change with histo titles
  int min = ((int)(Emin/kpb)) * kpb;
  int max = ((int)(Emax/kpb)) * kpb;
  int binN = (max-min)/kpb;
  int kpbN = 1e2; // keV per bin for rock neutrons : change with histo titles
  int Nmin = ((int)(NeuMin/kpbN)) * kpbN;
  int Nmax = ((int)(NeuMax/kpbN)) * kpbN;
  int nBinN = (Nmax-Nmin)/kpbN;

  recRate = TH1F("","",binN,min,max); // former title: Event Rate per eRecoil (keV)
  recRate.SetMinimum(0);
  recRate.GetYaxis()->SetTitle("Events (cm^-2 yr^-1 (10keV)^-1)"); // change when kpb is changed
  recRate.GetXaxis()->SetTitle("Recoil Energy (keV) (10keV bins)");      // ditto

  recIntg = TH1F("","",binN,min,max); // former title: Total Recoil Rate for thresholds (keV)
  recIntg.SetMinimum(0);
  recIntg.GetYaxis()->SetTitle("Events over threshold (cm^-2 yr^-1)");
  recIntg.GetXaxis()->SetTitle("Recoil threshold (keV)");

  recRateX = TH1F("","",binN,min,max); // former title: Event Rate per eRecoil (keV)
  recRateX.SetMinimum(0);
  recRateX.GetYaxis()->SetTitle("Events in m^3 (yr^-1 (10keV)^-1)"); // change when kpb is changed
  recRateX.GetXaxis()->SetTitle("Recoil Energy (keV) (10keV bins)");      // ditto

  recIntgX = TH1F("","",binN,min,max); // former title: Total Recoil Rate for thresholds (keV)
  recIntgX.SetMinimum(0);
  recIntgX.GetYaxis()->SetTitle("Events in m^3 over threshold (yr^-1)");
  recIntgX.GetXaxis()->SetTitle("Recoil threshold (keV)");

  neuRate = TH1F("","",nBinN,Nmin,Nmax); // former title: Neutron Rate per energy (keV)
  neuRate.GetYaxis()->SetTitle("Events (cm^-2 yr^-1 (50MeV)^-1)"); // change when kpbN is changed
  neuRate.GetXaxis()->SetTitle("Rock Neutron Energy (keV) (50MeV bins)");      // ditto

  neuIntg = TH1F("","",nBinN,Nmin,Nmax); // former title: Total Neutron flux for thresholds (keV)
  neuIntg.GetYaxis()->SetTitle("Events over threshold (cm^-2 yr^-1)");
  neuIntg.GetXaxis()->SetTitle("Rock Neutron threshold (keV)");  

  RockNeutron = TF1("","(x<450)*2e-2 + ((x>=450) && (x<700))*(2e-2 + (2e-3-2e-2)/(700-450)*(x-450)) + ((x>=700) && (x<1.5e3))*(2e-3 + (5e-4-2e-3)/(1.5e3-700)*(x-700)) + ((x>=1.5e3) && (x<2800))*(5e-4 + (2e-3-5e-4)/(2800-1500)*(x-1500)) + ((x>=2800) && (x<=4500))*(2e-3 + (1e-6-2e-3)/(4500-2800)*(x-2800))",100,4500);
  RockNeutron.GetYaxis()->SetTitle("Events (cm^-2 yr^-1 (50MeV)^-1)");
  RockNeutron.GetXaxis()->SetTitle("Rock Neutron energy (keV)");  

  Xsection = TF1("","(x<1e-3)*(pow(10,((-3-log10(x))/2*1.14+.56))-3.6) + (x<1e6)*Xsection2(x) + ((x>=1e6) && (x<6.7e6))*(pow(10,(.56-(log10(x)-6)/.82*.32))) + (x>=6.7e6)*1.7",1e-5,3.5e9);
  Xsection.GetYaxis()->SetTitle("F n cross section (barn)");
  Xsection.GetXaxis()->SetTitle("Neutron energy (eV)");  


  double recoilE, Xrate, neuE, rho;
  double rhoF = 2.4e25*1.0e-6/180.0*pressure;
  // double SecYr = 365.25*24*3600;
  double vol = 1e6; // volume of m^3 detector, in cm^3
  if (m3==0) {vol = 10e3;}
  double barn2cm = 1e-28;


  cout << "Events Simulated:" << endl;

  for (double n=1; n<=eventN; n++) {

    int nn = n;
    if (nn%50000==0) {cout << nn << "  " << flush;}
    


 //    //// Generate neutron energy
//     while (1) {
//       neuE = gRandom->Rndm()*(4500-100) + 100;
//       if (neuE<450) {freq = 2e-2;}
//       if (neuE>=450 && neuE<700) {freq = 2e-2 + (2e-3-2e-2)/(700-450)*(neuE-450);}
//       if (neuE>=700 && neuE<1000) {freq = 2e-3 + (5e-4-2e-3)/(1000-700)*(neuE-700);}
//       if (neuE>=1000 && neuE<1800) {freq = 5e-4 + (2e-3-5e-4)/(1800-1000)*(neuE-1000);}
//       if (neuE>=1800 && neuE<=4500) {freq = 2e-3 + (1e-6-2e-3)/(4500-1800)*(neuE-1800);}
//       //    else {freq=0;}
//       double rndm = gRandom->Rndm()*2e-2;
//       if (freq > rndm) break;
//     }



    neuE = RockNeutron.GetRandom(); 



    //// Make Recoils

    double cosRecoil;

    while (1) {

      recoilE = gRandom->Rndm() * binN*kpb + min;

      //// Calculate cosine of recoil relative to WIMP direction
      double mRatio = recoilM/neuM; // mass ratio: target over projectile
      cosRecoil = sqrt(recoilE/(mRatio*neuE))*(1+mRatio)*0.5;     
      
      
      //// check if kinematics is valid   
      if (cosRecoil>-1 && cosRecoil<1) break;
      
    }

    Xrate = Xsection.Eval(neuE*1e3) * barn2cm * rhoF * vol;

    // above finds rate given xsection for neuE
    // multiplies neuE by 1e3 before eval. b/c neuE is in keV and fxn is in eV

    //// Fill histos

    recRate.Fill(recoilE,(double) (1.5e-2)/eventN*(max-min));
    recRateX.Fill(recoilE,(double) (1.5e-2)*Xrate/eventN*(max-min));
    neuRate.Fill(neuE,(double) (.33e-2)/eventN*(Nmax-Nmin));

  }

  //// Integrate histos

  double Bin;
  for(int i=1; i<=(max/kpb); i++) {
    Bin = 0;
     for(int j=i; j<=(max/kpb); j++) {
      Bin += recRate.GetBinContent(j);
    }
    recIntg.SetBinContent(i,Bin);
  }

  double BinX;
  for(int i=1; i<=(max/kpb); i++) {
    BinX = 0;
     for(int j=i; j<=(max/kpb); j++) {
       BinX += recRateX.GetBinContent(j);  // sum up bins
    }
    recIntgX.SetBinContent(i,BinX);
  }

  double Dis;
  for(int i=1; i<=(Nmax/kpbN); i++) {
    Dis = 0;
    for(int j=i; j<=(Nmax/kpbN); j++) {
      Dis += neuRate.GetBinContent(j);
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
