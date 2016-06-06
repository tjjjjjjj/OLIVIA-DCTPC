#include "TRandom.h"
#include "TCanvas.h"

#include "TH1F.h"
#include "TF1.h"
#include "TAxis.h"

#include<iostream>

// --------------------------------------
// Useful info 
// --------------------------------------
// This code is designed to estimate the number of neutrons expected 
// in your DMTPC detector in 365 days of run time given the following 
// parameters: 
//
// N of events to simulate, 50,000 is a good number 
// Depth (in Km w e, default: 3.1 kmwe == Gram Sasso )
// mass (in grams -- default: 500, 1 m3 at 100 torr)
// Thickness of poly shield (in cm, default 10 cm) 

// How to run
// In root: 
//  .L testN.cxx+   
//  testN(50000,3.1,500,40.)


// Known problems to be fixed
// ==========================
//1) <E_neutrons> is not changin with depth as it should WIPP=86 and GSand Boulby =82 GeV: 
// does not make sense 
//2) RF of recoils is not right -- fix will give factor 2 
// 

//// May be useful to access these outside scope of totalRecoils or makeHisto

// common stuff 
const double SecYr = 365.25*24*3600;
const double barn2cm = 1e-24;
const double N_A = 6.02e23;  // Avogadro's number 

double neuM    = 1e6;      // mass in eV

// Histos for cosmic muons 
// =======================
TH1F recRate; // recoil rate
TH1F recIntg; //        integral
TH1F recRateX; // recoil rate after cross section
TH1F recIntgX; //        integral
TH1F neuRate; // neutron rate
TH1F neuIntg; //         integral
TH1F recRateGS; // recoil rates in KeV 
TCanvas c1,c2,c3;
TF1 MuonNeutron, Xsection;
int mpb_DONT_MODIFY = 1; //  MeV per bin : change with histo titles

// Histos for Rock neutrons  

TH1F recRateRock; // recoil rate
TH1F recIntRockg; //        integral
TH1F recRateRockX; // recoil rate after cross section
TH1F recIntRockgX; //        integral
TH1F neuRateRock; // neutron rate
TH1F neuIntgRock; //         integral

TCanvas c1rock; 
TF1 RockNeutron;


////////////////////////////////////////////////////////////////
//// There are 2 parts to this code
//// 
//// Neutrons from Cosmics 
//// =====================
//// Returns total number of recoils between two energy bounds
//// i.e. the integral of the recoil rate bewteen two energies
//// Output: (float) total number of recoils per YEAR (TBC) 
//// 
//// Neutrons from Rocks 
//// ===================
//// Returns total number of recoils between two energy bounds in 1 year of datataking 
////
////////////////////////////////////////////////////////////////


double Xsection2(double x) {
  double r = .2;
  return (3.6 + 136*exp(-pow((x-2.7e4),2)/(r*4e6)) + 47*exp(-pow((x-4.9e4),2)/(r*16e6)) + 25*exp(-pow((x-1e5),2)/(r*5.29e8)) + 8.7*exp(-pow((x-3.3e5),2)/(r*5.6e10)));
}



double getTotalNflux (double h0=3.1 ) {   //   h0= vertical depth in km w e -- flat overburden 
  double p0= 4.e-7; // in cm-2 s-1 
  double p1=0.86; // in km w e 
  double flux = p0 * p1 / h0 * exp(-h0/p1); 
  cout << " nutron flux for depth = " << h0 << " km.w.e. is "<< flux << " cm-2 s-1" << endl; 
  return flux; 
}

double getEmu (double h0=3.1 ) {   //  h0=3.1;// vertical depth in km w e -- flat overburden 
  double b= 0.383 ;  double epsilon = 618. ; double gamma=3.7 ; // lipari et al 
  //  double b= 0.4 ;  double epsilon = 693. ; double gamma=3.77 ; // Groom et al 
  double Emu = epsilon * (1.-exp(- b * h0))/(gamma-2.) ; // Mei Hime eq (9) 
  cout << " Average E muon  for depth = " << h0 << " km.w.e. is "<< Emu << " GeV" << endl; 
  return Emu; 
}
//

void makeMuonHisto(int eventN, double depth=3.1, double mass=500) {
  //  depth of site, km.w.e., mass in grams 
  double Emin   = 0;        // max and min for recoil energy range (in what unit??) 
  double Emax   = 1e2;
  double NeuMin = 10;       // max and min for neutron energy range (in what unit??) 
  double NeuMax = 3.5e3;

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

  recRateGS = TH1F("","",1000,0,10000); // former title: Neutron Rate per energy (keV)
  recRateGS.GetYaxis()->SetTitle("Events (cm^-2 s^-1) "); // change when mpbN is changed
  recRateGS.GetXaxis()->SetTitle("Muon Neutron Energy (keV) (10keV bins)");      // ditto

  neuIntg = TH1F("","",nBinN,Nmin,Nmax); // former title: Total Neutron flux for thresholds (MeV)
  neuIntg.GetYaxis()->SetTitle("Events over threshold (cm^-2 s^-1)");
  neuIntg.GetXaxis()->SetTitle("Muon Neutron threshold (MeV)");    

  // Fitting function for neutron energy distribution

  // common constants 
  double normalization = 1; 
  double A_mu = 1.; // 0.8e-9/4.42;

  // Site specifict constants -- use Gran Sasso as default 
  double a0=7.828;   double a1=2.23; double a2=-7.505e-15; double a3=2.831; double E_mu=270. ; // E_mu is in GeV 
  double totFlux = 2.7e-9; // cm-2 s-1 

  if( fabs(depth-3.1)<0.2 ) {
    cout << "SITE = Gran Sasso -- use default constants" << endl; 
  } else if  ( fabs(depth-2.8)<0.15 ) {
    cout << "SITE = Boulby mine " << endl; 
    a0=7.882;    a1=2.212;  a2=-2.342e-14;  a3=2.613;  E_mu=250. ; // E_mu is in GeV 
    totFlux = 4.86e-9; // cm-2 s-1 
  } else if  ( fabs(depth-1.585)<0.15 ) {
    cout << "SITE = WIPP " << endl; 
    a0=6.86;    a1=2.1;  a2=2.971e-13;  a3=2.456;  E_mu=180. ; // E_mu is in GeV 
    totFlux = 34.1e-9; // cm-2 s-1 
  } else if  ( fabs(depth-0.3)<0.1 ) {
    cout << "________ SITE = 300 mwe NB: guess many parameters !!! __________ " << endl; 
    a0=6.86;    a1=2.1;  a2=2.971e-13;  a3=2.456;  
    E_mu=getEmu(depth) ; // E_mu is in GeV and from Mei Hime formula  
    // totFlux = 3.0e-5; // cm-2 s-1 from COUPP slides 
    totFlux = getTotalNflux(depth) ; // cm-2 s-1 according to Mei Hime eq (13) 
  } else if  ( depth<=0.020 ) {
    cout << "________ Shallow site at "<< depth*1000.<< " mwe NB: guess many parameters !!! __________ " << endl; 
    a0=6.86;    a1=2.1;  a2=2.971e-13;  a3=2.456;  
    E_mu=getEmu(depth) ; // E_mu is in GeV 
    totFlux = getTotalNflux(depth) ; // cm-2 s-1 from Mei Hime eq (13)
    cout << "________ Shallow site at "<< depth*1000.<< " mwe NB: guess many parameters !!! __________ " << endl; 

  } else if  ( fabs(depth-4.3)<0.2 ) {
    cout << "________ SITE = Homestake 4.3 kmwe __________ " << endl; 
    a0=7.882;    a1=2.212;  a2=-2.342e-14;  a3=2.613; a3=2.456;  E_mu=300. ; // E_mu is in GeV change -- this is just a guess 
    totFlux = 1.0e-9; // cm-2 s-1 this is just a guess 
  } else if  ( fabs(depth-6.011)<0.1 ) {
    cout << "________ SITE = SNOLab 6.1 kmwe __________ " << endl; 
    a0=7.774;    a1=2.134;  a2=-2.939e-16;  a3=2.859;  E_mu=340. ; // E_mu is in GeV 
    totFlux = 0.054e-9; // cm-2 s-1 from Mei Hime??? 
  } else {
    cout << "************** Unknown DEPTH ***************"<< depth << endl; 
    exit(1); 
  }

  // Spectrum of neutrons from Mei-Hime eq 14
  double B_mu=0.324-0.641*exp(-0.014*E_mu); 
  char cexpress[1000];
  sprintf( cexpress, "(%f)*( exp(-(%f)*x)/x + (%f)*exp(-(%f)*x) ) + (%f)*pow(x,-(%f))", A_mu, a0, B_mu, a1, a2, a3 );
  MuonNeutron = TF1("", cexpress,.01,3.5); // eq 14  Mei-Hime @GranSasso
  MuonNeutron.GetYaxis()->SetTitle("Events (cm^-2 s^-1 (50MeV)^-1)");
  MuonNeutron.GetXaxis()->SetTitle("Muon Neutron energy (GeV)");  
  normalization = totFlux/MuonNeutron.Integral(0.001,3.5)/eventN; // 1 MeV to 3.5 GeV
  cout << " normalization = " << normalization << endl;

  //TF1 Xsection for Fluorine 
  Xsection = TF1("","(x<1e-3)*(pow(10,((-3-log10(x))/2*1.14+.56))-3.6) + (x<1e6)*Xsection2(x) + ((x>=1e6) && (x<6.7e6))*(pow(10,(.56-(log10(x)-6)/.82*.32))) + (x>=6.7e6)*1.7",1e-5,3.5e9);

  // For CF4 1 m3 
  double  A=19.0; // atomic number of the element (F=19., S=32.)  

  Xsection.GetYaxis()->SetTitle("F n cross section (barn)");
  Xsection.GetXaxis()->SetTitle("Neutron energy (eV)");  

  double recoilE,  neuE;
  double recoilM = 19e6;     // mass of F in eV
 
  cout << "Events Simulated:" << endl;

  for (int n=1; n<=eventN; n++) {
    int nn = n;
    if (nn%5000==0) {cout << nn << "  " << flush;}    
    neuE = MuonNeutron.GetRandom()*1e3;  // E neutron in MeV sampled from input distribution 
    //// Make Recoils
    double cosRecoil;
    while (1) {
      recoilE = gRandom->Rndm() * binN*mpb + min;
      //// Calculate cosine of recoil relative to WIMP direction
      double mRatio = recoilM/neuM; // mass ratio: target over projectile : this works only with F!!! 
      cosRecoil = sqrt(recoilE/(mRatio*neuE))*(1+mRatio)*0.5;                 
      //// check if kinematics is valid   
      if (cosRecoil>-1 && cosRecoil<1) break;
    }
    
    //// Fill histos

    double shift = pow(10,((3.1-depth)/1.8));
    recRate.Fill(recoilE,(double) (1.05e-12)*shift/eventN*(max-min));
    //    recRateX.Fill(recoilE,(double) (1.05e-12)*Xrate*shift/eventN*(max-min));

    double normalization2 = normalization * Xsection.Eval(neuE*1e6)*barn2cm * SecYr *(mass/A*N_A) ;
    // cout << "Xsection * barn2cm "  <<  Xsection.Eval(neuE*1e6)* barn2cm  << endl;
    recRateX.Fill(recoilE,normalization2);
 
    recRateGS.Fill(recoilE*1000.,normalization2); 

    //neuRate.Fill(neuE,(double) (3e-13)*shift/eventN*(Nmax-Nmin)); // old normalization by Albert: I don't understand it 
    neuRate.Fill(neuE,normalization);

  }

//   // back of an envelope calculations 
//   cout << endl; 
//   cout << " Simulated neutron flux  between 1 MeV and 3.5 GeV =" << neuRate.Integral()  
//        << " (cm-2 s-1 expected 8.1 10-10 cm-2 s-1 at GS from Mei-Hime table V)" << endl; 
//   double flux2=neuRate.Integral()  ; 
//   double NparticlesInTarget = mass / A *N_A ; // number of target particles 
//   cout << "NparticlesInTarget "<< NparticlesInTarget << endl; 
//   double nInteractionsPerSec2= flux2 * NparticlesInTarget * xsection ; 
//   double nInteractionsPerYear2 = nInteractionsPerSec2*SecYr;
//   cout << " GS ----- expected number of neutrons / m3 of cf4: " << nInteractionsPerYear2 << "-------" << endl; 

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

  
  c1.Clear();    c1.Divide(2,2);
  //  gPad->SetLogx();
  c1.cd(1);  gPad->SetLogy();neuRate.Draw();
  c1.cd(2);  gPad->SetLogy();  neuIntg.Draw();
  //  c1.cd(3);    recRate.Draw();  c1.cd(4);  recIntg.Draw();
  c1.cd(3);    recRateX.Draw();  c1.cd(4);  recIntgX.Draw();
  cout << endl;

  c2.Clear();  c2.cd();  c2.SetLogx(1);  c2.SetLogy(1); Xsection.Draw(); 
  c3.Clear();   c3.cd();  c3.SetLogy(1);  recRateGS.Draw();
  
}

void makeRockHisto(int eventN, double mass, double A=19 , double polythickness=10. ) {

  // ****************************************
  // All bins are 1 keV in this functions!!!!
  // All energies are given in KeVs 
  // Use spectrum from Mei Hime fig 27 which assumes 10 cm of poly 
  // then scale according to how much poly thickness in cm 
  // ****************************************

  double Emin = 0;    // max and min for recoil energy range in keV 
  double Emax = 1e3;
  double NeuMin = 1e2;    // max and min for neutron energy range in keV 
  double NeuMax = 4.5e3;

  double recoilM = A*1.e6;          // mass of F in eV (only here use eV) 

  //// Initialize all histos
  //   
  int kpbRock =  1;                 // keV per bin for recoils -- do not change or it will make mess 
  int kpbRockN = 1;                 // keV per bin for rock neutrons  -- do not change or it will make mess 

  int min = ((int)(Emin/kpbRock)) * kpbRock;
  int max = ((int)(Emax/kpbRock)) * kpbRock;
  int binN = (max-min)/kpbRock;

  int Nmin = ((int)(NeuMin/kpbRockN)) * kpbRockN;
  int Nmax = ((int)(NeuMax/kpbRockN)) * kpbRockN;
  int nBinN = (Nmax-Nmin)/kpbRockN;

  recRateRock = TH1F("","",binN,min,max); // former title: Event Rate per eRecoil (keV)
  recRateRock.SetMinimum(0);
  recRateRock.GetYaxis()->SetTitle("Events (cm^-2 yr^-1 keV^-1)"); // change when kpbRock is changed
  recRateRock.GetXaxis()->SetTitle("Recoil Energy (keV) (keV bins)");      // ditto

  recIntRockg = TH1F("","",binN,min,max); // former title: Total Recoil Rate for thresholds (keV)
  recIntRockg.SetMinimum(0);
  recIntRockg.GetYaxis()->SetTitle("Events over threshold (cm^-2 yr^-1)");
  recIntRockg.GetXaxis()->SetTitle("Recoil threshold (keV)");

  recRateRockX = TH1F("","",binN,min,max); // former title: Event Rate per eRecoil (keV)
  recRateRockX.SetMinimum(0);
  recRateRockX.GetYaxis()->SetTitle("Events in m^3 (yr^-1 (keV)^-1)"); // change when kpbRock is changed
  recRateRockX.GetXaxis()->SetTitle("Recoil Energy (keV) (keV bins)");      // ditto

  recIntRockgX = TH1F("","",binN,min,max); // former title: Total Recoil Rate for thresholds (keV)
  recIntRockgX.SetMinimum(0);
  recIntRockgX.GetYaxis()->SetTitle("Events in m^3 over threshold (yr^-1)");
  recIntRockgX.GetXaxis()->SetTitle("Recoil threshold (keV)");

  neuRateRock = TH1F("","",nBinN,Nmin,Nmax); // former title: Neutron Rate per energy (keV)
  neuRateRock.GetYaxis()->SetTitle("Events (cm^-2 yr^-1 keV-1"); // change when kpbRockN is changed
  neuRateRock.GetXaxis()->SetTitle("Rock Neutron Energy (keV) (keV bins)");      // ditto

  neuIntgRock = TH1F("","",nBinN,Nmin,Nmax); // former title: Total Neutron flux for thresholds (keV)
  neuIntgRock.GetYaxis()->SetTitle("Events over threshold (cm^-2 yr^-1)");
  neuIntgRock.GetXaxis()->SetTitle("Rock Neutron threshold (keV)");  


  RockNeutron = TF1("","(x<450)*2e-2 + ((x>=450) && (x<700))*(2e-2 + (2e-3-2e-2)/(700-450)*(x-450)) + ((x>=700) && (x<1.5e3))*(2e-3 + (5e-4-2e-3)/(1.5e3-700)*(x-700)) + ((x>=1.5e3) && (x<2800))*(5e-4 + (2e-3-5e-4)/(2800-1500)*(x-1500)) + ((x>=2800) && (x<=4500))*(2e-3 + (1e-6-2e-3)/(4500-2800)*(x-2800))",1,4500); // Hime fig 27 (in cm-2 y-1 MeV-1)
  RockNeutron.GetYaxis()->SetTitle("Events (cm^-2 yr^-1 MeV^-1)");
  RockNeutron.GetXaxis()->SetTitle("Rock Neutron energy (keV)");  
  

  // double totFlux = 3.78e-6 ;// total flux of n from rock in #/cm2/s measured at GS 
  double normalization = 1.e-3*RockNeutron.Integral(0.001,4500.)/eventN; 


  // Effect of the poly shielding 
  // parametrize effect seen in Mei-Hime Fig 28 top 
  double poly_scale = pow(10,(10-polythickness)/14.); 
  cout << " Scale factor due to poly shielding (thickness = "<< polythickness << ") = "<< poly_scale << endl; 
  normalization *= poly_scale ; // test


  if (A==19) {
    cout << " ***********Using CF4!!! ************" << endl; 
    Xsection = TF1("","(x<1e-3)*(pow(10,((-3-log10(x))/2*1.14+.56))-3.6) + (x<1e6)*Xsection2(x) + ((x>=1e6) && (x<6.7e6))*(pow(10,(.56-(log10(x)-6)/.82*.32))) + (x>=6.7e6)*1.7",1e-5,3.5e9);
    Xsection.GetYaxis()->SetTitle("F n cross section (barn)");
  }else   if (A==76) {
    cout << " ***********Using Germanium!!! ************" << endl; 
    Xsection = TF1("","(x<0.1)*(pow(10,((-1-log10(x))/4*1.82+1.07))) + ((x>=0.1) && (x<40))*(pow(10,((1.6-log10(x))/2.6*0.5+.57))) + ((x>=40) && (x<1e4))*Xsection2(x) + (x>=1e4)*(pow(10,((7.4-log10(x))/3.4*.69+0.4)))",1e-5,3.5e9);
    Xsection.GetYaxis()->SetTitle("Ge n cross section (barn)");
  } else  {
    cout << "Unknown element with A=" << A << endl; 
    exit(1); 
  } 
 Xsection.GetXaxis()->SetTitle("Neutron energy (eV)");  
  
  double recoilE,neuE;//  rho,  Xrate, ;

  cout << "Events Simulated:" << endl;

  for (int n=1; n<=eventN; n++) {
    int nn = n;
    if (nn%5000==0) {cout << nn << "  " << flush;}    
    neuE = RockNeutron.GetRandom(); // in keV 
    //// Make Recoils
    double cosRecoil;
    while (1) {
      recoilE = gRandom->Rndm() * binN*kpbRock + min;
      //// Calculate cosine of recoil relative to WIMP direction
      double mRatio = recoilM/neuM; // mass ratio: target over projectile
      cosRecoil = sqrt(recoilE/(mRatio*neuE))*(1+mRatio)*0.5;     
      //// check if kinematics is valid   
      if (cosRecoil>-1 && cosRecoil<1) break;
    }

    double normalization2 = normalization * Xsection.Eval(neuE*1e3)*barn2cm *(mass*N_A/A) ;// no need for *SecYr: this is already in y-1

    //// Fill histos

    recRateRock.Fill(recoilE,(double) (1.5e-2)/eventN*(max-min));
    recRateRockX.Fill(recoilE,normalization2);
    neuRateRock.Fill(neuE,normalization);

  }


    // back of an envelope calculations 
    cout << endl; 
    cout << " Simulated neutron flux  between 100 MeV and 4.5 GeV =" << neuRateRock.Integral(100,4500)  
	 << " (cm-2 y-1 -- expected ~1.5 10-2 cm-2 y-1 from Mei-Hime fig 27)" << endl; 
    double flux2=neuRateRock.Integral()  ; 
    double NparticlesInTarget = mass / A *N_A ; // number of target particles 
    cout << "NparticlesInTarget "<< NparticlesInTarget << endl; 
    double xsection = 5.e-24; // in cm2
    double nInteractionsPerYear= flux2 * NparticlesInTarget * xsection ; 
    cout << " GS ----- expected number of neutrons / m3 of cf4: " << nInteractionsPerYear << "-------" << endl; 


  //// Integrate histos

  double Bin;
  for(int i=1; i<=(max/kpbRock); i++) {
    Bin = 0;
     for(int j=i; j<=(max/kpbRock); j++) {
      Bin += recRateRock.GetBinContent(j);
    }
    recIntRockg.SetBinContent(i,Bin);
  }

  double BinX;
  for(int i=1; i<=(max/kpbRock); i++) {
    BinX = 0;
     for(int j=i; j<=(max/kpbRock); j++) {
       BinX += recRateRockX.GetBinContent(j);  // sum up bins
    }
    recIntRockgX.SetBinContent(i,BinX);
  }

  double Dis;
  for(int i=1; i<=(Nmax/kpbRockN); i++) {
    Dis = 0;
    for(int j=i; j<=(Nmax/kpbRockN); j++) {
      Dis += neuRateRock.GetBinContent(j);
    }
    neuIntgRock.SetBinContent(i,Dis);
  }
    
  //// Draw histos  

  c1rock.Clear(); c1rock.Divide(2,2); c1rock.SetLogy();  
  c1rock.cd(1); c1rock.SetLogy();  neuRateRock.Draw();  c1rock.cd(2);  c1rock.SetLogy(); neuIntgRock.Draw();  
  // c1rock.cd(3);  recRateRock.Draw();  c1rock.cd(4);  recIntRockg.Draw();
  c1rock.cd(3);  recRateRockX.Draw();  c1rock.cd(4);  recIntRockgX.Draw();

  //   cout << endl;
  // Plot xsection 
  //  c2rock.Clear();  c2rock.cd();  c2rock.SetLogx(1);  c2rock.SetLogy(1);  Xsection.Draw();
  

}




void  testN(int nEvents=50000,  double depth=3.1 , double mass=500, double polythickness=10. ) {

  //    gStyle->SetOptFit();
  // Inputs: 
  // =======
  // nevents= 1,000 - 10,000 ie more than enough  

  // depth in km.w.e.  3.1 = Gran Sasso, 2.8=Boulby, 1.58 = WIPP, (0.3= NUMI) 

  // mass in g: 500 g of CF4 / 1 m3 at 100 torr 

  // Generate histograms for muon histos 
  makeMuonHisto(nEvents, depth, mass );  
  
  // count number of recoils between certain energies 
  cout <<" Muon recoils between 0 and 10 MeV " << recRateGS.Integral(0,10000) << endl;
  cout <<" Muon recoils between 50 and 150 keV "<< recRateGS.Integral(5,15) << endl;
  cout <<" Muon recoils between 30 and 50 keV "<< recRateGS.Integral(3,5) << endl;
  cout <<" Muon number of recoils between 10 and 50 keV "<< recRateGS.Integral(1,5) << endl;
  cout <<" ----------------->> Muon recoils between 80 and 200 keV "<< recRateGS.Integral(8,20) << endl;

  double A=19;  // atomic mass of element considered (F=19., S=32., GE=76 ... )  
  makeRockHisto(nEvents, mass, A, polythickness) ; 


  cout << endl; 
  // count nuclear recoils  between Emin and Emax keV
  cout << "N recoils from Rock between 0 and 99999 keV = " << recRateRockX.Integral(0,10000) << endl; 
  cout << "N recoils from Rock between 50 and 150 keV = " << recRateRockX.Integral(50,150) << endl; 
  cout << "N recoils from Rock between 30 and 150 keV = " << recRateRockX.Integral(30,150) << endl; 
  cout << "N recoils from Rock between 30 and 100 keV = " << recRateRockX.Integral(30,100) << endl; 
  cout <<" ----------------->> N recoils from Rock between 80 and 200 keV "<< recRateRockX.Integral(80,200)  << endl;

}

