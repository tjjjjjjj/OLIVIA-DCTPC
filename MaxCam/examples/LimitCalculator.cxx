#include "TF1.h"
#include "TSystem.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TProfile.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLatex.h"
#include "../MaxCamWIMP.hh"
#include "../Tools.hh"

#include "TFeldmanCousins.h"
#include "math.h"
#include <iostream>
using std::cout;
using std::endl;

// 
// LimitCalculator() calculates limits (SD and SI) for various experiments
// selectExperiment: 1=CDMS, 2=Xenon, 3=DMTPC, 4=NewAge 
//
// NB: this code is not intended for publication (too many shortcuts)
// but it is useful to set QUICK limits 
// Some debugging still needed... 
// 
// Use: in root 
// .L LimitCalculator.cxx+
// LimitCalculator(3,100)
// 
// Author: G. Sciolla  
// 


double cosmoBackground ( double E1, double E2, double upTime, double activeMass, bool hasVeto, double depth ) {
  // Returns # of events from cosmic ray originated neutrons 
  // -------------------------------------------------------
  // depth in mwe
  // E1 and E2 are thresholds in keV 
  // upTime in days 
  // activeMass in kg 
  
  double recoils_keV_kg_y = 0.; 
  double x=depth/1000.; 
  recoils_keV_kg_y = 2.0*exp(-2.9*(x-0.4))+0.00025*exp(-1.3*(x-5.0)); 
  cout << " recoils_keV_kg_y = "<< recoils_keV_kg_y << endl ; 
  double background = recoils_keV_kg_y * (E2-E1) * upTime / 365. * activeMass; 
  if (hasVeto)  background =  background/100.; 
  return    background; 

}

double reducedMass (double m1, double m2){ 
  return m1*m2/(m1+m2); 
    }

double FC90(double observedEvents, double estimatedBackg ) { 

  TFeldmanCousins* fc = new TFeldmanCousins(); 
  double observedEventsD = observedEvents; 
  if ( observedEvents <= 28. ) {
    return fc->CalculateUpperLimit(observedEventsD,estimatedBackg);   
    // Call Double_t CalculateUpperLimit(Double_t Nobserved, Double_t Nbackground)
    // in http://root.cern.ch/root/html/TFeldmanCousins.html#TFeldmanCousins:TFeldmanCousins
  }  else {
    cout <<" Warning: Feldman-Cousins code is not designed for this... " << endl;  
    double signal = observedEvents-estimatedBackg; 
    double error = sqrt(observedEvents); 
    double limitRate = signal + 1.7 * error;  // NB: 1.7 is not exactly 90% CL: check; this works only if signal>0 
    return limitRate; 
  } 
}

double CL90(int observedEvents) { 

  const double vup[11]={2.3,3.89,5.32,6.68,7.99,9.27,10.53,11.77,12.99,14.21,15.41}; 
	  
  if(observedEvents>10) return observedEvents; 
  if(observedEvents<0) { 
    cout << "************* Error in 90percentRate: observedEvents is negative" << 
      observedEvents << endl; 
    return 0.; 
  }

  return vup[observedEvents]; 
}

void LimitCalculator(int selectExperiment=3, const int npoints=100, const double minDMmass=20.) {

  // Plotting options 
  bool plot100Kgy(true); 

  // Init 
  double  best_MW_SD=999.; double  best_MW_SI=999.;   double  best_MW_SD_N=999. ; 
  double sigma_SI_best=999.; double sigma_SD_best=999.; double sigma_SD_N_best=999.; 
  int nUnpairedProtons  = 0; 
  int nUnpairedNeutrons = 0; 

  // Init physics constants
  double rho_D = 0.3; // DM halo density in GeV / cm3
  double v_0 = 220;   // Average DM velocity in km/sec

 // default: CDMS Feb 2008 
 // Physics 
  TString ExpName = "CDMS 2008";  
  double Atarget =72. ; // mass target recoil in GeV c^2 (F=19, Xe=131, Ge=72)
  double spinFactor =0.065 ; // lambda^2*J(J+1) 
  double fractionUsableSD = 0.0773; // 73Ge only 7.7% of all Ge 
  // Detector   
  double activeMass = 3.75; // mass in Kg 
  double upTime = 105.;     // up time in days 
  double detEfficiency = 0.35 ; // detector efficiency 
  double E1=10.; double E2=100.0; // Lower and uper thresholds in keV
  // Data 
  double R0_measured=0.;           // Measured rate 
  double Exp_background=0.;        // Expected background 
 
 if( selectExperiment==1) { 
   nUnpairedProtons  = 0;  nUnpairedNeutrons = 1; 
 }else if( selectExperiment==2) { 
   ExpName = "Xenon 10";  
   Atarget =131. ; // mass target recoil in GeV c^2 (F=19, Xe=131, Ge=72)
   spinFactor =0.124 ; // lambda^2*J(J+1) 0.124 for Xe 129 (Xe 131 0.055)
   fractionUsableSD = 0.26; // Or 0.21 or something in between??? FIX THIS 
   activeMass = 5.4; // mass in Kg 
   upTime = 58.60;     // up time in days 
   detEfficiency = 0.19 ; // detector efficiency 
   E1=4.5;  E2=26.90; // Lower and uper thresholds in keV
   R0_measured=10.;           // Measured rate  
   Exp_background=10.;        // Expected background 
   nUnpairedProtons  = 0;  nUnpairedNeutrons = 1; 
 } else if( selectExperiment==3) { 
   ExpName = "DMTPC 0.1 kg-y";   
   Atarget =19. ; // mass target recoil in GeV c^2 (F=19, Xe=131, Ge=72)
   spinFactor =0.647 ; // lambda^2*J(J+1) 
   fractionUsableSD = 0.864; // 12(C)/(12+4*19(F)) 
   activeMass = 0.1; // mass in Kg -- 1m3@50 torr = 0.250 kg 
   upTime = 365.0;     // up time in days 
   detEfficiency = 0.99 ; // detector efficiency 
   E1=50.;  E2=150.0; // Lower and uper thresholds in keV
   R0_measured=0.;        // Measured rate  
   Exp_background = 0.01 * (E2-E1) * upTime / 365. * activeMass; 
   nUnpairedProtons  = 1;  nUnpairedNeutrons = 0; 

   // Neutrons from Rock 
   double recoilsFromNeutronsFromRock = 2.49482956860447302e-03; // events / cm^2 / year for E=50<->150; 
   double surface = 100.*100.*6.; // surface of active area in cm^2
   double recoilsFromNeutronsFromRockTotal = recoilsFromNeutronsFromRock * upTime/365. * surface ; 
   // Neutrons from cosmic muons 
   bool hasVeto= true; 
   //   double cosmoNeutronEvts = cosmoBackground ( E1, E2, upTime, activeMass, hasVeto, 2000. ); // SUSEL 2000 mwe 
   double cosmoNeutronEvts = cosmoBackground ( E1, E2, upTime, activeMass, hasVeto, 4800. ); // DUSEL 4800 mwe 
   // TOTAL Neutron background 
   cout << " ======== cosmo vs recoilsFromNeutronsFromRockTotal vs Exp_background  " << cosmoNeutronEvts 
	<< " " << recoilsFromNeutronsFromRockTotal <<" " <<Exp_background << endl;        
   // 


 } else if( selectExperiment==4) {  
   ExpName = "  NewAge surface run ";  
   Atarget =19. ; // mass target recoil in GeV c^2 (F=19, Xe=131, Ge=72)
   spinFactor =0.647 ; // lambda^2*J(J+1) 
   fractionUsableSD = 0.864; // 12(C)/(12+4*19(F)) 
   activeMass = 9.e-3;    // mass in Kg 
   upTime = 16.7;         // up time in days 
   detEfficiency = 0.40 ; // detector efficiency 
   // Test
//    E1=100.;  E2=110.0;    // Lower and uper thresholds in keV
//    R0_measured=250.;        // Measured rate  
//    Exp_background = 250.;   // just guessing...  
   // Original 
   E1=100.;  E2=400.0;    // Lower and uper thresholds in keV
   R0_measured=1686.;        // Measured rate  
   Exp_background = 1686.;   // just guessing...  
   nUnpairedProtons  = 1;  nUnpairedNeutrons = 0; 
 } else if( selectExperiment==5) {  
   ExpName = "  NewAge Tanimori paper - 3 m^3*year at 30 torrfor 1 year (>35 keV)"; 
   Atarget =19. ; // mass target recoil in GeV c^2 (F=19, Xe=131, Ge=72)
   spinFactor =0.647 ; // lambda^2*J(J+1) 
   fractionUsableSD = 0.864; // 12(C)/(12+4*19(F)) 
   activeMass = 0.45;    // mass in Kg = 3 m3 of CF4 at 30 torr
   upTime = 365.;         // up time in days 
   detEfficiency = 0.99 ; // detector efficiency 
   E1=35.;  E2=200.0;    // Lower and uper thresholds in keV
   R0_measured=0.;        // Measured rate  
   Exp_background = 0.;   // just guessing...  
   nUnpairedProtons  = 1;  nUnpairedNeutrons = 0; 
 } else if( selectExperiment==6) {  
   ExpName = " DMTPC 100 kg-y "; 
   activeMass = 100.; // mass in Kg -- 1m3@50 torr = 0.250 kg 
   upTime = 365.0;     // up time in days 
   Atarget =19. ; // mass target recoil in GeV c^2 (F=19, Xe=131, Ge=72)
   spinFactor =0.647 ; // lambda^2*J(J+1) 
   fractionUsableSD = 0.864; // 12(C)/(12+4*19(F)) 
   detEfficiency = 0.99 ; // detector efficiency 
   E1=50.;  E2=150.0; // Lower and uper thresholds in keV

   //   Exp_background = 0.01 * (E2-E1) * upTime / 365. * activeMass;  // this seems to be ok for deep site and poly=10 cm --  check 
   Exp_background = 3.6; // testN.cxx: 3 from muons at 3.1 kmwe , 0.6 from rock -- poly = 40 cm 

   R0_measured=Exp_background ;        // Measured rate  
   nUnpairedProtons  = 1;  nUnpairedNeutrons = 0; 
 } else if( selectExperiment==7 || selectExperiment==8 ) { 
     Atarget =19. ; // mass target recoil in GeV c^2 (F=19, Xe=131, Ge=72)
     spinFactor =0.647 ; // lambda^2*J(J+1) 
     fractionUsableSD = 0.864; // 12(C)/(12+4*19(F))  
     detEfficiency = 0.99 ; // detector efficiency 
     E1=50.;  E2=150.0; // Lower and uper thresholds in keV
     if( selectExperiment==7) { 
       //       ExpName = "DMTPC at SUSEL (2 kg-d)"; 
       ExpName = "DMTPC underground (2 kg-d)"; 
       activeMass = 2*(.15*.15*.25)*0.5; // mass in Kg -- 1m3@100 torr = 0.500 kg  
       upTime = 365.0;     // up time in days 
       R0_measured=0.;        // Measured rate  
       
       // NB: this calculation does not work yet!!! 
       // +++++++++++++++++++++++++++++++++++++++++ 
      // Neutrons from Rock 
       double recoilsFromNeutronsFromRock = (0.1-0.01)*1.e-3 ;        //       = recNeuFromRock(E1,E1); // in cm^-2 year-1
       double surface = 15.*15.*5.*2.; // surface of active area in cm^2
       double recoilsFromNeutronsFromRockTotal = recoilsFromNeutronsFromRock * upTime/365. * surface ; 
       // Neutrons from cosmic muons 
       bool hasVeto= true; 
       double cosmoNeutronEvts = cosmoBackground ( E1, E2, upTime, activeMass, hasVeto, 300. ); //// SUSEL 300 mwe no veto 
       // TOTAL Neutron background 
       Exp_background = cosmoNeutronEvts + recoilsFromNeutronsFromRockTotal; 
       // Exp_background = 10; // take number from Denis (no shielding) and divide it by 2 (will have direction)
       // NB: shielding reduces neutrons from rocks by ~2+ orders of magnitude
       cout << " ======== cosmo vs recoilsFromNeutronsFromRockTotal " << cosmoNeutronEvts << " " << recoilsFromNeutronsFromRockTotal << endl;        
     } else if (selectExperiment==8 ) { 
     ExpName = "DMTPC DUSEL Large Detector -- 1 tonn x 3 years @ 72 mwe "; 
     activeMass = 1000.; // mass in Kg -- 1 tonn 
     upTime = 3*365.0;     // up time in days -- 3 years 
     R0_measured=0.;        // Measured rate  
     //   Exp_background = 0.01 * (E2-E1) * upTime / 365. * activeMass; // Soudan depth, no veto  
     Exp_background = 2.5e-4 * (E2-E1) * upTime / 365. * activeMass; // 2000 mwe depth, mu veto  
     //   Exp_background = 3.e-2 * (E2-E1) * upTime / 365. * activeMass; // 2000 mwe depth, NO mu veto  
     }
     cout << "ActiveMass = "<< activeMass << " Kg "<< endl; 
     if(  Exp_background >= 1.)  R0_measured= Exp_background ; 
     nUnpairedProtons  = 1;  nUnpairedNeutrons = 0; 
 } else if( selectExperiment==99) { 
   ExpName = "TEST CDMS with 100% eff and large mass ";  
   activeMass = 0.8; // mass in Kg 
   upTime = 365;     // up time in days 
   detEfficiency = 0.99 ; // detector efficiency 
 } else { 
   cout <<" !!!!! There is no experiment # "<< selectExperiment << " will assume CDMS !!!!!"<< endl; 
   ExpName = "********************** None ******************" ; 
 } 
 cout << "Experiment selected = " << ExpName << endl; 
 double MT = Atarget*0.932;    
 
 //double limitRate = CL90(R0_measured); // set 90% CL limit no background 
  double limitRate = FC90(R0_measured, Exp_background); // set 90% CL limit 
  cout << "Observed evts=" << R0_measured << " Exp background="<< Exp_background << " --> FC 90% =" << limitRate << endl; 
  double R0_eff_corr = limitRate/detEfficiency;  // correct total rate for det efficiency 
  double exposure = activeMass*upTime; 
  cout << "Exposure = "<< exposure << " kg*days"<< endl; 

  double m_DM[npoints] ; 
  double x_SI_DM[npoints]; 
  double x_SD_DM[npoints]; 
  double x_SD_DM_N[npoints]; 

  //  float massVal[]={4 ,5,6 ,7 ,8 ,10 ,20 ,30 ,40 ,50 ,60 ,70 ,80	,90 ,100 ,200 ,300 ,500	,1000,5000,10000 }; 
  double MD; 
  for (int i=0;i<npoints;i++) { 
    // double MD =  massVal[i] ; 
    //    double MD =  minDMmass+i*1. ; 
    if(i==0) { MD=minDMmass;}
    else { MD=m_DM[i-1]*1.2; }
    m_DM[i]= MD; 


    // From partial rate to total rates 
    // cout << "Partial rate="<< R0_eff_corr << endl; 
    double c1=0.751; double c2=0.561; // with these values we will only integrate up to 0.75 or so... ???? 
    double r =4.*MD*MT/pow((MD+MT),2) ; // kinematic factor 
    double E0=MD*(1.e6)/2. * pow((230./300000.),2); // 10^6 because E0 is in keV
    
    //  cout << "E0 = " << E0<< " E1 ="<< E1<< endl; 
    //  cout << c1 << " " << c2 << " " <<  exp(-c2*E1/(E0*r)) <<" "<< exp(-c2*E2/(E0*r)) << endl; 
    
    double R0 = R0_eff_corr * c2 /c1 / (exp(-c2*E1/(E0*r)) - exp(-c2*E2/(E0*r)) ) ; 
    //    cout << "R0="<< R0<< endl; 
    
    //    // Debug: Denis's version 
    //    double Rdenis = MaxCamWIMP::dRdE_ERange(E1, E2, MD, MT); // for sigma=1pb, mass=1kg, days=1
    //    Rdenis =  Rdenis * exposure; 
    //    double sigmaDenis = R0_eff_corr/Rdenis; 
    double mup= reducedMass(0.932,MD); 
    double mun= reducedMass(0.9396,MD); 
    double mut= reducedMass(MT,MD); 
    
    // original calculation using (3.7) 
    double sigma_0 = R0 / exposure * MD * MT / 503. * 0.4/rho_D * 230/v_0 ; 
    // new calculation using (6.5) -- same results 
    //    double sigma_0 = R0 / exposure * pow(mut,2) /r / 126. ; 

    //    cout << "My sigma="<< sigma_0 <<" Denis' sigma = "<< sigmaDenis << endl; 

// correct for FF 
    double rn = pow(Atarget,0.33333); 
    double ER = (E2+2*E1)/3.; 
    double qrn = rn * sqrt(2.*MT*ER) * 1.e-3 / 0.197; 
    double F2 = 0.047; 
    if (qrn<=2.55 || qrn>=4.5) F2 = pow(sin(qrn)/qrn,2); 
    //    cout << "qrn = "<< qrn << " F2="<< F2 << endl; 
    double sigma_0_FF =  sigma_0 /F2; 
    
    // normalize to N nucleons 
    double I=pow(Atarget,2); // interaction function for spin independent x-section 
    double sigma_over_I = sigma_0_FF/I ; 

    // Spin independent x-section 
    double sigma_SI = sigma_over_I *pow(mup,2)/pow(mut,2) ; 

    // Spin dependent x-section 
    double CWp2= 0.46; // for Higssino-p 
    double CWn2= 0.34; // for Higssino-n 

    double CWN2= 1. ; 
    if(nUnpairedProtons ==1 )  CWN2= CWp2 ; 
    if(nUnpairedNeutrons ==1 ) CWN2= CWn2 ; 

    double spinFactorProton = 0.75; 
    double spinFactorNeutron = 0.75; 
    double correctionFactor = spinFactorProton/spinFactor * pow(mup/mut,2)* CWp2/CWN2;
    double correctionFactorNeutron = spinFactorNeutron/spinFactor * pow(mun/mut,2)* CWn2/CWN2;

    //    cout << correctionFactor << endl; 

    double sigma_SD =   sigma_0 * correctionFactor / fractionUsableSD; 
    double sigma_SD_N = sigma_0 * correctionFactorNeutron / fractionUsableSD; 

    //    cout <<  sigma_SD << ", "; 

//     cout << "CorrectionFactorNeutron="<< correctionFactorNeutron<<endl; 
//     cout << "CorrectionFactor for proton="<< correctionFactor<<endl; 
    // Translate pb--> cm^2: 1 pb = 10^-40 m^2 = 10^-36 cm^2 = 1000 fb 

    x_SI_DM[i]=sigma_SI/1.e36; 
    x_SD_DM[i]=sigma_SD/1.e36; 
    x_SD_DM_N[i]=sigma_SD_N/1.e36; 

    if( x_SI_DM[i]<=sigma_SI_best){ 
      best_MW_SI= m_DM[i];
      sigma_SI_best=x_SI_DM[i]; 
    }
    if( x_SD_DM[i]<= sigma_SD_best){ 
      best_MW_SD= m_DM[i];
      sigma_SD_best=x_SD_DM[i]; 
    }
    if( x_SD_DM_N[i]<= sigma_SD_N_best){ 
      best_MW_SD_N= m_DM[i];
      sigma_SD_N_best=x_SD_DM_N[i]; 
    }
    
    // Print out input for Filippini's files 
    //    cout << m_DM[i] << " " << x_SD_DM[i] << endl; 


  }

  cout << endl; 
  // Print out info of best limit 
  cout << "Best limit for sigma_SI = " << sigma_SI_best << " at mW~" << best_MW_SI << "GeV"<< endl; 
  cout << "Best limit for sigma_SD_proton = " << sigma_SD_best << " at mW~" << best_MW_SD << "GeV"<< endl; 
  cout << "Best limit for sigma_SD_neutron = " << sigma_SD_N_best << " at mW~" << best_MW_SD_N << "GeV"<< endl; 

  // plot results now 


  TGraph* graph = new TGraph(npoints, m_DM , x_SI_DM ); 
  TGraph* graphSD = new TGraph(npoints, m_DM , x_SD_DM ); 
  //  TGraph* graph = new TGraph(npoints, m_DM , m_DM ); 
  graph->SetMarkerStyle(21);
  graph->SetMarkerColor(2);
  graph->GetXaxis()->SetTitle("mass wimp (GeV)");
  graph->GetYaxis()->SetTitle("normalized SI x-section(cm-2)");

  graphSD->SetMarkerStyle(21);
  graphSD->SetMarkerColor(1);
  graphSD->GetXaxis()->SetTitle("mass wimp (GeV)");
  graphSD->GetYaxis()->SetTitle("normalized SD x-section(cm-2)");
  
  
  //    TMultiGraph* mg  = new TMultiGraph();
  //    mg->Add(graph);
  TLatex l; l.SetTextSize(0.035);  
  int verticalAxis=500; 
  if(plot100Kgy)  verticalAxis=600; 
  TCanvas* myCanvas = new TCanvas("myCan","myCan", 500, verticalAxis);
  myCanvas->SetLogy(); myCanvas->SetLogx(); 
  graph->SetMaximum(4.e-36);graph->SetMaximum(1.e-40); // GS why 2 setmax? 
  graph->Draw("ALP*");

  TCanvas* myCanvas2 = new TCanvas("myCan2","myCan2", 500,verticalAxis);
  myCanvas2->SetLogy(); myCanvas2->SetLogx(); 
  // graphSD->SetMaximum(1.);graphSD->SetMaximum(0.);
  // graphSD->Draw("ALP*");


    graphSD->SetLineWidth(3);
    graphSD->GetHistogram()->SetXTitle("WIMP mass (GeV)");
    graphSD->GetHistogram()->SetYTitle("#sigma_{p} (cm^{2})");
    if( selectExperiment==7) graphSD->SetMinimum(1e-37); // use for SUSEL 
    if( selectExperiment==3) graphSD->SetMinimum(5e-39); // use for DMTPC 0.1 kg y 
    if(plot100Kgy)  graphSD->SetMinimum(3e-45); // official when 100 kg-day is there
    graphSD->SetMaximum(4e-36);//000);
    graphSD->Draw("AL");
    double angle1=50.; 
    if(plot100Kgy) angle1=28.; 
    l.SetTextSize(0.030); // -4 
    //  l.SetTextAngle(angle1-4);     l.DrawLatex(190,1.30e-38,ExpName); // official 
    l.SetTextAngle(angle1);    l.DrawLatex(190,1.30e-37,ExpName); // SUSEL 

    // Add other experiments... 
    double x_position = 300.; 

    // COUPP 2008 - 15 keV threshold 
    //    float x_coupp[]={  4,   5,   6,   10,   20,   40,   60,   100,  200,  400,  700,  900, 1000};
    //    float y_coupp[]={ 50,   6, 2.3,  0.6, 0.36, 0.28, 0.28,  0.36, 0.59,    1,   1.7,   2,   2.3};
 float x_coupp[]={  4	,5	,6	,7	,8	,10	,20	,30	,40	,50	,60	,70	,80	,90	,100	,200	,300	,500	,1000	,5000	,10000	}; 

    float y_coupp[]={ 53.29131189495615,5.472543275525056,2.21090078408438,1.302150885904067,0.922430172157557
,0.6090229609299633,0.3282000325903187,0.282498995487303,0.2685579414388432,0.2695133159394477,0.2778315846200281,0.2903239727625864,0.3053801984483322
,0.3220604018815363,0.3398467813142201,0.5409040415082203,0.7537899968738539,1.186223883492666,2.275000882311175,1.467687092012786,21.92672598281424}; 

    for (int i=0;i<20;i++) y_coupp[i]*=1.e-36;
    TGraph* g_coupp= new TGraph(20, x_coupp, y_coupp);
    g_coupp->SetLineWidth(2);  g_coupp->SetLineStyle(2); 
    g_coupp->Draw("C");
    l.SetTextAngle(angle1); l.SetTextSize(0.025);
    l.DrawLatex(x_position,0.89e-36,"COUPP (2008)");


    // KIMS 2007 - 3409 kg-days CsI SD-proton'Cross-section data for axial interaction normalised to p 
    // data_reference='arXiv:0704.0423' // Data from Jeff Filippini 
    //      float x_kims[]={ 20,   30,   40,   50,   60,    100,  200,  300,   400,  800,   1000}; // Denis scan 
    //  float y_kims[]={ 13.7, 0.36, 0.20, 0.17, 0.16,  0.17, 0.22, 0.28,  0.36,  0.7,   0.8}; // Denis scan 

       float x_kims[]={11.2366, 12.6054, 14.1434, 16.0786, 19.0022, 21.0653, 23.9639, 27.2668, 31.8404, 36.2381, 41.2431, 46.3409, 52.7493, 60.0454, 
68.3525, 
77.8087, 
88.5755, 
100.8346,
114.7905,
130.6812,
148.7754,
169.3707,
192.8267,
219.5312,
249.934, 
284.5545,
323.9705,
368.8557,
419.9596,
478.1438,
544.3892,
619.8128,
705.686, 
803.4568,
914.7965, 1041.5653, 1185.8711, 1350.2043,
1537.31, 
1750.344,
1992.8493,
2269.0101,
2583.4402,
2941.4427,
3349.0557,
3813.1541,
4341.5653,
4943.2015,
5628.3522,
6408.3061,
7296.3428,
8307.6498,
9458.8897,  9963.435};

       float y_kims[]={ 9.6473e-35,  2.8973e-35,  9.8662e-36,  3.742e-36,  1.3692e-36,  8.8985e-37,  5.4801e-37,  3.8962e-37,  2.8203e-37,  2.3995e-37,  2.0415e-37
,  1.9e-37,  1.8004e-37,  1.7369e-37,  1.706e-37,  1.6756e-37,  1.6756e-37, 1.706e-37, 1.7369e-37
, 1.8004e-37
, 1.9e-37
, 1.9695e-37
, 2.1161e-37
, 2.2737e-37
,  2.4429e-37
, 2.6724e-37
, 2.9234e-37
, 3.2559e-37
, 3.6262e-37
, 4.0387e-37
, 4.498e-37
, 5.0096e-37
, 5.5794e-37
, 6.214e-37
, 7.0462e-37
,7.9898e-37
,8.8985e-37
,1.009e-36
, 1.1441e-36
,1.2974e-36
,1.4449e-36
,1.6384e-36
,1.8578e-36
,2.1066e-36
,2.3887e-36
,2.7086e-36
,3.0714e-36
,3.4827e-36
,4.0206e-36
,4.559e-36
,5.1695e-36
,5.968e-36
,6.7672e-36
,7.4028e-36};

    TGraph* g_kims= new TGraph(40, x_kims, y_kims);
    g_kims->SetLineWidth(2);  g_kims->SetLineStyle(2);  
    g_kims->Draw("C"); l.DrawLatex(x_position,0.32e-36,"KIMS (2007)");
    
    // 'NAIAD 2005 Final SD-proton Spin-dependent proton cross section'
    //data_reference='G. J. Alner et al. [UKDMC], Phys. Lett. B 616 (2005), 17.'
    //private_comment='Final SDp limit from NAIAD NaI experiment, obtained using DataThief from preprint by J. Filippini'
    float x_naiad[]={10.    ,14.36  ,20.62  ,30.599 ,43.94  ,63.096 ,90.603 ,134.45 ,227.58 ,425.18 ,936.33 ,2275.8 ,10000    };
    float y_naiad[]={   2.1614 , 0.85271, 0.46491, 0.33641, 0.27483, 0.26393, 0.27483, 0.32308, 0.42879, 0.69661, 1.4425, 3.3723, 13.887};
    for (int i=0;i<14;i++) y_naiad[i]*=1.e-36;
    TGraph* g_naiad= new TGraph(14, x_naiad, y_naiad);
    g_naiad->SetLineWidth(2);  g_naiad->SetLineStyle(2); 
    //    g_naiad->Draw("C");    l.DrawLatex(x_position,0.55e-36,"NAIAD (2005)");


    // Latest Picasso data // all data is wrong now: just cut and paste from naiad!!! 
    // *******************
    int nPicasso=6;
    float x_picasso[]={ 5.,   5.5,   9.,    25.,    100.,   1000 }; 
    float y_picasso[]={ 3.,   2. ,  5.e-1,   1.5e-1, 3.5e-1,  3.  };
    for (int i=0;i<nPicasso;i++) y_picasso[i]*=1.e-36;
    TGraph* g_picasso= new TGraph(nPicasso, x_picasso, y_picasso);
    g_picasso->SetLineWidth(2);  g_picasso->SetLineStyle(2); 
    g_picasso->Draw("C");    l.DrawLatex(x_position,0.55e-36,"PICASSO (2009)");


    // Theory 
    // ======
    // Ellis, Ferstl, Olive  PRD63, Fig5c tanbeta=10 - MSSM 
    int nPoints = 19; 
    float x_th[]={52,      53,  60,  70,  79,   81, 100,  130,  135, 141,  156,  161,  171, 186,  194,  200, 300, 400, 500};
    float y_th[]={0.001,  0.1, 0.6, 2.5, 3.9, 4.02, 0.6, 0.36, 0.29, 0.3, 0.25, 0.25, 0.22, 0.12, 0.12, 0.1, 0.06,0.004,  0.002};

    // Fig 7 SD 
    // Original numbers I used the paper
    // int nPoints = 16; 
    // float x_th[]={50.,   60.,  70.,  80.,   90, 100, 110, 120, 130, 140., 150, 160,   200,  300,  400,    500. };
    // float y_th[]={0.015, 0.55, 3.,   4.5,  1.5, .65, .5, .4,   .35, .3,   0.2, 0.15, 0.12, 0.007, 0.004, 0.0025    };

    // plot it on the graph 
    for (int i=0;i<nPoints;i++) y_th[i]*=1.e-36; 
    TGraph* g_th= new TGraph(nPoints, x_th, y_th);
    for (int i=0;i<g_th->GetN();i++) g_th->GetY()[i]*=1e-3; 
    g_th->SetFillColor(7);
    g_th->Draw("LF");
    l.SetTextAngle(0); l.SetTextSize(0.025); l.DrawLatex(70,1.5e-40,"MSSM");

    // More recent Ellis theory prediction CMSSM hep-ph/0001005
    int nPoints2 = 7; 
    float x_th2[]={50.,          60.,   600.,    700.,   80.,     60,   50 };
    float y_th2[]={1.e-5,        2.e-1, 5.e-5,  5.e-6, 1.e-3,  5.e-4,  1.e-5};

    for (int i=0;i<nPoints2;i++) y_th2[i]*=1.e-36; 
    TGraph* g_th2= new TGraph(nPoints2, x_th2, y_th2);
    for (int i=0;i<g_th2->GetN();i++) g_th2->GetY()[i]*=1e-3; 
    g_th2->SetFillColor(8);
    g_th2->Draw("LF");
    l.SetTextAngle(0); l.SetTextSize(0.025); l.DrawLatex(70,0.5e-41,"CMSSM");


    // DMTPC 100 Kg years at soudan depth 
    float x_dmtpc100[]={4,5,6 ,7 ,8 ,10 ,20 ,30 ,40 ,50 ,60 ,70 ,80	,90 ,100 ,200 ,300 ,500	,1000,5000,10000 }; 
    float y_dmtpc100[]={8.54258e+11, 5.5229e+06, 6163.13, 85.0696, 4.69108, 0.12806, 0.000500451, 0.000138271, 8.4499e-05, 6.7502e-05, 6.06133e-05, 5.77306e-05, 5.67993e-05, 5.69531e-05, 5.77684e-05, 7.74674e-05, 0.000102109, 0.00015369, 0.000285014 }; 
   for (int i=0;i<19;i++) y_dmtpc100[i]*=1.e-36; 
   TGraph* g_dmtpc100= new TGraph(19, x_dmtpc100, y_dmtpc100);    g_dmtpc100->SetLineWidth(3);
   l.SetTextAngle(angle1-4); l.SetTextSize(0.030);     
   if(plot100Kgy) {
     g_dmtpc100->Draw("C");    l.DrawLatex(190,1.1e-40,"DMTPC 100 kg-y");
   }
    // to here 

  cout << "CFR published results: " << endl; 
  if( selectExperiment==1) { 
    cout << "CDMS 2008 SI limit= 6.6 10-44 at 60 GeV - poisson " << endl; 
    cout << "CDMS 2008 limit on SD_n 2.7 10-38 at 60 GeV " << endl; 
  } else if ( selectExperiment==2) { 
    cout << "XENON 10 2007 SI limit= 4.5 10-44 at 30 GeV " << endl; 
    cout << "XENON 10 2007 limit on SD_n  6*10-39 at 30 GeV " << endl; 
    cout << "XENON 10 2007 limit on SD_p  5?*10-37 at 30 GeV " << endl; 
  } else if ( selectExperiment==3) { 
    cout << "CDMS 2008 SI limit= 6.6 10-44 at 60 GeV " << endl; 
    cout << "CDMS 2008 limit on SD_n 2.7 10-38 at 60 GeV " << endl; 
  } else if ( selectExperiment==4) { 
    cout << "NewAge surface run SI limit= 6.6XX 10-44 at 120? GeV " << endl; 
    cout << "NewAge surface run SD_p 2 10-32 at 120? GeV " << endl; 
  }

  cout << " Done!" << endl; 
  
// Exp results to use as reference 

// 1) CDMS 2008 
// Using 230 km/s,  Rate=0; 
// Published SI limit= 6.6 10-44 at 60 GeV (now I get 8.2 10-44 at 70GeV)
// Published: SD_n 2.7 10-38 at 90% CL Me:sigma_SD_neutron = 2.7 e-38 at mW~70 GeV -- Perfect 
// 
// 2) Xenon 10 
// SI published =    Me = 20.1485e-44 at mW~46GeV
// DSp published =   Me = 1.41167e-38 at mW~46GeV
// DSn published =   Me = 1.06015e-38 at mW~46GeV
// 
// 3) DMTPC 0.1 kg year 220 km/s 0.3 GeV/cm3
// Denis: 
// Me: sigma_SD_proton 6.7 10-39 at mW~83GeV
//
// 4) NewAge 2007 
// SD = 2x10+4 pb at 100 GeV (2 x 10-32 cm2) ---> (I find 260x-34 cm2 off by 2 orders of magnitude wrt publication)

// Giuliani does not include FF 
// CDMS   

}
