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
#include "TLegend.h"
#include "TRolke.h"
#include "TFeldmanCousins.h"
#include "math.h"
#include <iostream>
using std::cout;
using std::endl;

// 
// SurfaceRunSensitivityCF() calculates limits (SD and SI) for various experiments
// selectExperiment: 1=CDMS, 2=Xenon, 3=DMTPC, 4=NewAge 
//
// NB: this code is not intended for publication (too many shortcuts)
// but it is useful to set QUICK limits 
// Some debugging still needed... 
// 
// Use: in root 
// .L SurfaceRunSensitivityCF.cxx+
// SurfaceRunSensitivityCF(100)
// 
// Author: G. Sciolla  
//
// Modified for WR2: T. Caldwell, J. Monroe 



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
  if (observedEvents>20) fc->SetMuMax(200.);
  fc->SetMuStep(0.005);
  return fc->CalculateUpperLimit(observedEventsD,estimatedBackg); 
  
// old way: (jrm)
//
//  if ( observedEvents <= 28. ) {
//    return fc->CalculateUpperLimit(observedEventsD,estimatedBackg);   
//    // Call Double_t CalculateUpperLimit(Double_t Nobserved, Double_t Nbackground)
//    // in http://root.cern.ch/root/html/TFeldmanCousins.html#TFeldmanCousins:TFeldmanCousins
//  }  else {
//
//    cout <<" Warning: Feldman-Cousins code is not designed for this... " << endl;  
//    double signal = observedEvents-estimatedBackg; 
//    double error = sqrt(observedEvents); 
//    double limitRate = signal + 1.7 * error;  // NB: 1.7 is not exactly 90% CL: check; this works only if signal>0 
//    return limitRate; 
//  } 
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




void SurfaceRunSensitivityCF(const int npoints=100, const double minDMmass=20.) {

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
  
  double m_DM[npoints] ; 
  double x_SI_DM[7][npoints]; 
  double x_SD_DM[7][npoints]; 
  double x_SD_DM_N[7][npoints]; 
  double limit[npoints];

  for (int selectExperiment=0; selectExperiment<7; selectExperiment++) {
    
    // Init 
    double  best_MW_SD=999.; 
    double  best_MW_SI=999.;   
    double  best_MW_SD_N=999. ; 
    double sigma_SI_best=999.; 
    double sigma_SD_best=999.; 
    double sigma_SD_N_best=999.; 
    int nUnpairedProtons  = 0; 
    int nUnpairedNeutrons = 0; 

    if (selectExperiment==0) { 
      
      ExpName = "Xenon 10";  
      Atarget =131. ; // mass target recoil in GeV c^2 
      spinFactor =0.124 ; // lambda^2*J(J+1) 0.124 for Xe 129 (Xe 131 0.055)
      fractionUsableSD = 0.26; // Or 0.21 (see paper)
      activeMass = 5.4; // mass in Kg 
      upTime = 58.60;     // up time in days 
      detEfficiency = 0.19 ; // detector efficiency 
      E1=4.5;  E2=26.90; // Lower and uper thresholds in keV
      R0_measured=10.;           // Measured rate  
      Exp_background=10.;        // Expected background 
      nUnpairedProtons  = 0;  nUnpairedNeutrons = 1; 
      
    } else if (selectExperiment==1) { 
    
      ExpName = "  NewAge surface run ";  
      Atarget =19. ; // mass target recoil in GeV c^2 
      spinFactor =0.647 ; // lambda^2*J(J+1) 
      fractionUsableSD = 0.864; // 12(C)/(12+4*19(F)) 
      activeMass = 9.e-3;    // mass in Kg 
      upTime = 16.7;         // up time in days 
      detEfficiency = 0.40 ; // detector efficiency 
      E1=100.;  E2=400.0;    // Lower and uper thresholds in keV
      R0_measured=1686.;        // Measured rate  
      Exp_background = 5.5*upTime;   // just guessing...  
      Exp_background = 20.;          // just guessing...
      nUnpairedProtons  = 1;  nUnpairedNeutrons = 0; 
      
    } else if (selectExperiment==2) {  
      
      Atarget =19. ; // mass target recoil in GeV c^2 (F=19, Xe=131, Ge=72)
      spinFactor =0.647 ; // lambda^2*J(J+1) 
      fractionUsableSD = 0.864; // 12(C)/(12+4*19(F))  
      detEfficiency = 0.99 ; // detector efficiency 
      E1=50.;  E2=150.0; // Lower and uper thresholds in keV
      ExpName = "DMTPC SUSEL"; 
      activeMass = 2*(.15*.15*.25)*0.5; // mass in Kg -- 1m3@100 torr = 0.500 kg  
      upTime = 365.0;     // up time in days 
      R0_measured=0.;        // Measured rate  
      
      // Neutrons from Rock 
      double recoilsFromNeutronsFromRock = (0.1-0.01)*1.e-3 ;        //       = recNeuFromRock(E1,E1); // in cm^-2 year-1
      double surface = 15.*15.*5.*2.; // surface of active area in cm^2
      double recoilsFromNeutronsFromRockTotal = recoilsFromNeutronsFromRock * upTime/365. * surface ; 
      // Neutrons from cosmic muons 
      bool hasVeto= true; 
      double cosmoNeutronEvts = cosmoBackground ( E1, E2, upTime, activeMass, hasVeto, 300. ); //// SUSEL 300 mwe no veto 
      // TOTAL Neutron background 
      Exp_background = cosmoNeutronEvts + recoilsFromNeutronsFromRockTotal; 
      cout << " ======== cosmo vs recoilsFromNeutronsFromRockTotal " << cosmoNeutronEvts << " " << recoilsFromNeutronsFromRockTotal << endl;        
      
    } else if (selectExperiment==3) {  // DMTPC surface run (36.9 gm-day)
      
      ExpName = "DMTPC 10L Surface Run";   
      Atarget =19. ;                        // mass target recoil in GeV c^2 
      spinFactor =0.647 ;                   // lambda^2*J(J+1) 
      fractionUsableSD = 0.864;             // 12(C)/(12+4*19(F)) 
      activeMass = 0.0033;                  // total mass in Kg
      upTime = (2135795*pow(24*3600,-1))/2.;// up time in days (per cam.) 
      detEfficiency = 0.70 ;                // detector efficiency (from 100 GeV WIMP MC)
      E1=80.0;  E2=200.0;                   // Lower and upper thresholds in keV
      R0_measured=105;                      // Measured rate 
      Exp_background = 74;                  // Predicted surface n bgnd 
                                            // (from n spectrum measured by Nakamura et al., Journal of Nuclear Science and Technology 42 No. 10, 843 (2005). 
      nUnpairedProtons  = 1;  nUnpairedNeutrons = 0; 

    } else if (selectExperiment==4) {  // DMTPC surface run (36.9 gm-day), 0 exp. bgnd
      
      ExpName = "DMTPC 10L Surface Run, 0 Expected Bgnd";   
      Atarget =19. ;                        // mass target recoil in GeV c^2
      spinFactor =0.647 ;                   // lambda^2*J(J+1) 
      fractionUsableSD = 0.864;             // 12(C)/(12+4*19(F)) 
      activeMass = 0.0033;                  // total mass in Kg
      upTime = (2135795*pow(24*3600,-1))/2.;// up time in days (per cam.)
      detEfficiency = 0.70 ;                // detector efficiency (from 100 GeV WIMP MC)
      E1=80.0;  E2=200.0;                   // Lower and upper thresholds in keV
      R0_measured=105;                      // Measured rate  
      //      R0_measured=110;              // Measured rate with 100% gain map variation
      Exp_background = 0;                   // Expected bgnd (conservative case)
      nUnpairedProtons  = 1;  nUnpairedNeutrons = 0; 

    } else if (selectExperiment==5) {       // DMTPC 4-shooter at WIPP      
      ExpName = "DMTPC 20L WIPP 1-year";   
      Atarget =19. ;                        // mass target recoil in GeV c^2 
      spinFactor =0.647 ;                   // lambda^2*J(J+1) 
      fractionUsableSD = 0.864;             // 12(C)/(12+4*19(F)) 
      activeMass = 0.0081;                  // total mass in Kg 
      upTime = 365.;                        // up time in days 
      detEfficiency = 0.7 ;                 // detector efficiency 
      E1=80.;  E2=200.0;                    // Lower and upper thresholds in keV
      R0_measured=0.;                       // Measured rate  
      Exp_background = 0.;                  // Expected bgnd
      nUnpairedProtons  = 1;  nUnpairedNeutrons = 0; 

    } else if (selectExperiment==6) {  // DMTPC 1m3 at WIPP
      
      ExpName = "DMTPC 1m3 WIPP 1-year";   
      Atarget =19. ;                        // mass target recoil in GeV c^2 
      spinFactor =0.647 ;                   // lambda^2*J(J+1) 
      fractionUsableSD = 0.864;             // 12(C)/(12+4*19(F)) 
      activeMass = 0.250;                   // mass in Kg -- 1m3@50 torr = 0.250 kg 
      upTime = 365.;                        // up time in days 
      //      upTime = 365.*0.4;            // up time in days, old value 
      detEfficiency = 0.5 ;                 // detector efficiency 
      E1=50.;  E2=200.0;                    // Lower and upper thresholds in keV
      R0_measured=0.;                       // Measured rate  
      Exp_background = 0.;                  // Expected bgnd
      nUnpairedProtons  = 1;  nUnpairedNeutrons = 0; 

    } 


    double MT = Atarget*0.932;    
    double limitRate=1000;
    limitRate = FC90(R0_measured, Exp_background); // set 90% CL limit 
    cout << "Observed evts=" << R0_measured << " Exp background="<< Exp_background << " --> 90% CL =" << limitRate << endl; 
    double R0_eff_corr = limitRate/detEfficiency;  // correct total rate for det efficiency 
    double exposure = activeMass*upTime; 
    cout << "Exposure = "<< exposure << " kg*days"<< endl; 
       
    for (int i=0;i<npoints;i++) { 
      double MD =  minDMmass+i*1. ; 
      m_DM[i]= MD; 
    
      // From partial rate to total rates 
      // cout << "Partial rate="<< R0_eff_corr << endl; 
      double c1=0.751; double c2=0.561; 
      double r =4.*MD*MT/pow((MD+MT),2) ; // kinematic factor 
      double E0=MD*(1.e6)/2. * pow((230./300000.),2); // 10^6 because E0 is in keV
      double R0 = R0_eff_corr * c2 /c1 / (exp(-c2*E1/(E0*r)) - exp(-c2*E2/(E0*r)) ) ; 
      double mup= reducedMass(0.932,MD); 
      double mun= reducedMass(0.9396,MD); 
      double mut= reducedMass(MT,MD);       
      // original calculation using (3.7) 
      double sigma_0 = R0 / exposure * MD * MT / 503. * 0.4/rho_D * 230/v_0 ; 
      // new calculation using (6.5) -- same results 
      //    double sigma_0 = R0 / exposure * pow(mut,2) /r / 126. ; 
      
      // correct for SD FF 
      double rn = pow(Atarget,0.33333); 
      // find the average value of the ff over the recoil spectrum
      double ff_sum = 0;
      double pdf_sum = 0;
      double e_sum = 0;
      double F2 = 0.047; 
      for (int estep=0; estep<100; estep++) {
	double ER = E1 + estep*(E2-E1)/100.;
	double qrn = rn * sqrt(2.*MT*ER) * 1.e-3 / 0.197; 
	if (qrn<=2.55 || qrn>=4.5) F2 = pow(sin(qrn)/qrn,2); 
	double pdf=(c1*R0/(E0*r))*exp(-c2*ER/(E0*r));
	ff_sum += pdf*F2;
	pdf_sum += pdf;
	e_sum += ER*pdf;
      }
      double ff_avg=ff_sum/pdf_sum;
      double e_avg = e_sum/pdf_sum;
      if (i==0) cout << "pdf sum: " << pdf_sum << "\t" << ff_sum << "\t" << ff_avg << endl;
      if (i==0) cout << "e sum: " << e_sum << "\t" << e_avg << "\t" << E1 << "\t" << E2 << endl;
      double sigma_0_FF =  sigma_0/ff_avg; 
      // old way
      //      double ER = (E2+2*E1)/3.; // what is this? (jrm)
      //      double qrn = rn * sqrt(2.*MT*ER) * 1.e-3 / 0.197; 
      //      double F2 = 0.047; 
      //      if (qrn<=2.55 || qrn>=4.5) F2 = pow(sin(qrn)/qrn,2); 
      //      double sigma_0_FF =  sigma_0/F2; 
      
      // Spin independent x-section
      double I=pow(Atarget,2); // normalize to N nucleons 
      double sigma_over_I = sigma_0_FF/I ; // this is not the correct FF for SI (jrm)
                                           // since **SD** is calculated above 
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
      
      double sigma_SD =   sigma_0_FF * correctionFactor / fractionUsableSD; 
      double sigma_SD_N = sigma_0_FF * correctionFactorNeutron / fractionUsableSD;
      // old way
      //       double sigma_SD =   sigma_0 * correctionFactor / fractionUsableSD; 
      //       double sigma_SD_N = sigma_0 * correctionFactorNeutron / fractionUsableSD;
      x_SI_DM[selectExperiment][i]=sigma_SI/1.e36; 
      x_SD_DM[selectExperiment][i]=sigma_SD/1.e36; 
      x_SD_DM_N[selectExperiment][i]=sigma_SD_N/1.e36; 
      if( x_SI_DM[selectExperiment][i]<=sigma_SI_best){ 
	best_MW_SI= m_DM[i];
	sigma_SI_best=x_SI_DM[selectExperiment][i]; 
      }
      if( x_SD_DM[selectExperiment][i]<= sigma_SD_best){ 
	best_MW_SD= m_DM[i];
	sigma_SD_best=x_SD_DM[selectExperiment][i]; 
      }
      if( x_SD_DM_N[selectExperiment][i]<= sigma_SD_N_best){ 
	best_MW_SD_N= m_DM[i];
	sigma_SD_N_best=x_SD_DM_N[selectExperiment][i]; 
      }
      
    } // end loop over npoints
    
    // Print out info of best limit 
    cout << "experiment " << selectExperiment << endl;
    cout << "Best limit for sigma_SI = " << sigma_SI_best << " at mW~" << best_MW_SI << "GeV"<< endl; 
    cout << "Best limit for sigma_SD_proton = " << sigma_SD_best << " at mW~" << best_MW_SD << "GeV"<< endl; 
    cout << "Best limit for sigma_SD_neutron = " << sigma_SD_N_best << " at mW~" << best_MW_SD_N << "GeV"<< endl; 
    cout << endl;  
    
  } // end loop over experiments





  // plot results
  for (int i=0; i<npoints; i++) {
    limit[i]=x_SD_DM[0][i];
  }
  TGraph *graphSD_xenon = new TGraph(npoints, m_DM, limit );

  for (int i=0; i<npoints; i++) {
    limit[i]=x_SD_DM[1][i];
  }
  TGraph *graphSD_newage_surface = new TGraph(npoints, m_DM, limit );

  for (int i=0; i<npoints; i++) {
    limit[i]=x_SD_DM[3][i];
  }
  TGraph *graphSD_dmtpc_surface_bgndsub = new TGraph(npoints, m_DM, limit );  
  
  // shaded region between tgraphs
  TGraph *graphSD_dmtpc_shade=new TGraph(npoints*2);
  double lim_min[npoints], lim_max[npoints];
  for(int i=0;i<npoints;i++){
    lim_min[i]=x_SD_DM[4][i];
    lim_max[i]=x_SD_DM[3][i];
  }
  for (int j=0;j<npoints;j++) {
    graphSD_dmtpc_shade->SetPoint(j,m_DM[j],lim_max[j]);
    graphSD_dmtpc_shade->SetPoint(npoints+j,m_DM[npoints-j-1],lim_min[npoints-j-1]);
  }
  graphSD_dmtpc_shade->SetFillStyle(1001);
  graphSD_dmtpc_shade->SetFillColor(17);
  graphSD_dmtpc_shade->SetLineColor(17);

  for (int i=0; i<npoints; i++) {
    limit[i]=x_SD_DM[4][i];
  }
  TGraph *graphSD_dmtpc_surface = new TGraph(npoints, m_DM, limit );

  for (int i=0; i<npoints; i++) {
    limit[i]=x_SD_DM[5][i];
  }
  TGraph *graphSD_dmtpc_wipp = new TGraph(npoints, m_DM, limit );  


  for (int i=0; i<npoints; i++) {
    limit[i]=x_SD_DM[6][i];
  }
  TGraph *graphSD_dmtpc_1m3wipp = new TGraph(npoints, m_DM, limit );  

 
  graphSD_xenon->GetXaxis()->SetTitle("WIMP Mass (GeV)");
  graphSD_xenon->GetYaxis()->SetTitle("normalized SD x-section(cm-2)");
  graphSD_xenon->SetMarkerStyle(21);
  graphSD_xenon->SetLineColor(0);
  graphSD_newage_surface->SetMarkerStyle(21);
  graphSD_newage_surface->SetLineColor(2);
  graphSD_newage_surface->SetLineStyle(1);
  graphSD_dmtpc_surface_bgndsub->SetLineStyle(1);
  graphSD_dmtpc_surface_bgndsub->SetLineColor(1);
  graphSD_dmtpc_surface->SetLineStyle(1);
  graphSD_dmtpc_surface->SetLineColor(1);
  graphSD_dmtpc_wipp->SetLineStyle(2);
  graphSD_dmtpc_wipp->SetLineColor(1);
  graphSD_dmtpc_1m3wipp->SetLineStyle(3);
  graphSD_dmtpc_1m3wipp->SetLineColor(1);

  graphSD_xenon->SetLineWidth(4);  
  graphSD_newage_surface->SetLineWidth(4);  
  graphSD_dmtpc_surface_bgndsub->SetLineWidth(4);  
  graphSD_dmtpc_surface->SetLineWidth(4);  
  graphSD_dmtpc_wipp->SetLineWidth(4);  
  graphSD_dmtpc_1m3wipp->SetLineWidth(4);  
  
  TCanvas* myCanvas2 = new TCanvas("myCan2","myCan2", 900,700);
  myCanvas2->SetLogy(); myCanvas2->SetLogx(); 
  graphSD_xenon->GetHistogram()->SetXTitle("WIMP mass (GeV)");
  graphSD_xenon->GetHistogram()->SetYTitle("#sigma_{p} (cm^{2})");
  graphSD_xenon->SetMinimum(1e-40);
  graphSD_xenon->SetMaximum(1e-31);
  graphSD_xenon->Draw("AL");

  //  graphSD_dmtpc_surface_bgndsub->Draw("L");

  graphSD_dmtpc_1m3wipp->Draw("L");
  graphSD_dmtpc_surface->Draw("L");
  graphSD_dmtpc_wipp->Draw("L");
   



 
  // Add other experiments... 
 
  // PICASSO SD-proton (2009)
  // data_comment='Spin-dependent proton cross section'
  // data_reference='S.Archambault et al., arXiv:0907.0307.'
  // private_comment='SD-proton limit from the PICASSO experiment, provided by B. Beltran'
  // data_units = {'GeV', 'cm^2'}

  double x_picasso[]={5  ,   
		      6  ,   
		      8  ,   
		      10 ,   
		      11 ,   
		      13 ,   
		      15 ,   
		      17 ,   
		      19 ,   
		      20 ,   
		      22 ,   
		      24 ,   
		      26 ,   
		      28 ,   
		      30 ,   
		      32 ,   
		      34 ,   
		      36 ,   
		      38 ,   
		      40 ,   
		      42 ,   
		      44 ,   
		      46 ,   
		      48 ,   
		      50 ,   
		      52 ,   
		      54 ,   
		      56 ,   
		      58 ,   
		      60 ,   
		      65 ,   
		      70 ,   
		      75 ,   
		      80 ,   
		      85 ,   
		      90 ,   
		      95 ,   
		      100,   
		      200,   
		      300,   
		      400,   
		      500,   
		      600,   
		      1000};
  


  double y_picasso[] = {
    3.54466     ,  
    1.55663	    ,
    0.700941    ,
    0.376159    ,
    0.306388    ,
    0.235629    ,
    0.194705    ,
    0.177283    ,
    0.168765    ,
    0.163855    ,
    0.159189    ,
    0.158893    ,
    0.158946    ,
    0.160897    ,
    0.161201    ,
    0.166342    ,
    0.169565    ,
    0.175807    ,
    0.178214    ,
    0.183983    ,
    0.186721    ,
    0.190652    ,
    0.194429    ,
    0.200735    ,
    0.205307    ,
    0.208719    ,
    0.215641    ,
    0.219226    ,
    0.225464    ,
    0.229827    ,
    0.242663    ,
    0.25461	    ,
    0.270585    ,
    0.286068    ,
    0.298864    ,
    0.313552    ,
    0.32256	    ,
    0.338196    ,
    0.616641    ,
    0.901799    ,
    1.1817	    ,
    1.45694	    ,
    1.73862     , 
    2.899    };

  for (int i=0; i<43; i++) y_picasso[i]*=1.E-36;
  
  TGraph* g_picasso= new TGraph(43, x_picasso, y_picasso);
  g_picasso->SetLineWidth(4);  g_picasso->SetLineStyle(4);  
  g_picasso->SetLineColor(4); 
  g_picasso->Draw("C"); 



  // COUPP 2010, arXiv:1008.3518v1 Fig. 4 "COUPP this result"
  //  private communication from E. Dahl, J. Hall
  float x_coupp[]={
            8.3176377e+00,
            9.1201084e+00,
            1.0000000e+01,
            1.0964782e+01,
            1.2022644e+01,
            1.3182567e+01,
            1.4454398e+01,
            1.5848932e+01,
            1.7378008e+01,
            1.9054607e+01,
            2.0892961e+01,
            2.2908677e+01,
            2.5118864e+01,
            3.1622777e+01,
            3.9810717e+01,
            5.0118723e+01,
            6.3095734e+01,
            7.9432823e+01,
            1.0000000e+02,
            1.2589254e+02,
            1.5848932e+02,
            1.9952623e+02,
            2.5118864e+02,
            3.1622777e+02,
            1.0000000e+03,
            3.1622777e+03,
            1.0000000e+04};

  
  float y_coupp[]={ 
            1.2700098e+01,
            4.7420422,
            2.1772187,
            1.1649219,
            6.8322304e-1,
            4.3792080e-1,
            3.0238619e-1,
            2.2045415e-1,
            1.6910525e-1,
            1.3489248e-1,
            1.1152213e-1,
            9.4791725e-2,
            8.2547754e-2,
            6.4020593e-2,
            5.5102379e-2,
            5.1630249e-2,
            5.1789034e-2,
            5.4807936e-2,
            6.0431327e-2,
            6.8733429e-2,
            8.0038427e-2,
            9.4863741e-2,
            1.1396366e-1,
            1.3832620e-1,
            3.9839102e-1,
            1.2242396,
            3.8368031}; 
  
  for (int i=0;i<27;i++) y_coupp[i]*=1.e-36;
  TGraph* g_coupp= new TGraph(27, x_coupp, y_coupp);
  g_coupp->SetLineWidth(4);  g_coupp->SetLineStyle(5);
  g_coupp->SetLineColor(6);
  g_coupp->Draw("C");
  


  // KIMS 2007 - 3409 kg-days CsI SD-proton'Cross-section data for axial interaction normalised to p 
  // data_reference='arXiv:0704.0423' 
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
  g_kims->SetLineWidth(4);  g_kims->SetLineStyle(4);  
  g_kims->SetLineColor(3); 
  //  g_kims->Draw("C"); 

  
  // 'NAIAD 2005 Final SD-proton Spin-dependent proton cross section'
  //data_reference='G. J. Alner et al. [UKDMC], Phys. Lett. B 616 (2005), 17.'
  //private_comment='Final SDp limit from NAIAD NaI experiment, obtained using DataThief from preprint by J. Filippini'
  float x_naiad[]={10.    ,14.36  ,20.62  ,30.599 ,43.94  ,63.096 ,90.603 ,134.45 ,227.58 ,425.18 ,936.33 ,2275.8     };
  float y_naiad[]={   2.1614 , 0.85271, 0.46491, 0.33641, 0.27483, 0.26393, 0.27483, 0.32308, 0.42879, 0.69661, 1.4425, 3.3723};
  for (int i=0;i<12;i++) y_naiad[i]*=1.e-36;
  TGraph* g_naiad= new TGraph(12, x_naiad, y_naiad);
  g_naiad->SetLineWidth(4);  g_naiad->SetLineStyle(3);
  g_naiad->SetLineColor(4);
  //  g_naiad->Draw("C");    



  // 'Newage first underground limit'
  //data_reference="First underground results with NEWAGE-0.3a direction-sensitive dark matter detector," http://arxiv.org/abs/1002.1794
  //private_comment='data scan by James Battat 02.11.10'
  // //Abstract: A direction-sensitive dark matter search experiment at
  //Kamioka underground laboratory with the NEWAGE-0.3a detector was
  //performed. The NEWAGE- 0.3a detector is a gaseous
  //micro-time-projection chamber filled with CF4 gas at 152 Torr. The
  //fiducial volume and target mass are 20*25*31 cm3 and 0.0115 kg,
  //respectively. With an exposure of 0.524 kgdays, improved
  //spin-dependent weakly interacting massive particle (WIMP)-proton
  //cross section limits by a direction-sensitive method were achieved
  //including a new record of 5400 pb for 150 GeV/c2 WIMPs. We studied
  //the remaining background and found that ambient gamma-rays
  //contributed about one-fifth of the remaining background and
  //radioactive contaminants inside the gas chamber contributed the
  //rest.
  float x_newage_ug[] = {
    30.4789,	
    32.2107,	
    34.0408,	
    35.9749,	
    38.3707,	
    40.9261,	
    43.6516,	
    47.4242,	
    49.6592,	
    51.9996,	
    54.4503,	
    57.0164,	
    59.7035,	
    62.5173,	
    65.4636,	
    70.4693,	
    74.4732,	
    77.983 ,	
    83.1764,	
    91.2011,	
    100    ,	
    107.647,	
    118.032,	
    130.617,	
    143.219,	
    158.489,	
    173.78 ,	
    190.546,	
    231.206,	
    275.423,	
    337.287,	
    387.258,	
    432.514,	
    515.229,	
    625.173,	
    744.732,	
    847.227,	
    946.237,	
    981.748};	

  float y_newage_ug[] = {
    35527.2,
    29767.2,
    25687.4,
    21522.7,
    18033.3,
    15561.7,
    13428.8,
    11588.3,
    10607.4,
    9709.49,
    9153.53,
    8629.41,
    8135.31,
    7669.49,
    7230.34,
    7020.29,
    6618.32,
    6239.36,
    6058.1 ,
    5882.11,
    5545.3 ,
    5384.21,
    5384.21,
    5227.79,
    5227.79,
    5384.21,
    5545.3 ,
    5711.22,
    6058.1 ,
    6426.05,
    7446.68,
    8135.31,
    8887.61,
    10000  ,
    11935  ,
    13830.6,
    15109.6,
    17000.7,
    17509.4};

  for (int i=0;i<39;i++) y_newage_ug[i]*=1.e-36;
  TGraph* g_newage_ug= new TGraph(39, x_newage_ug, y_newage_ug);
  g_newage_ug->SetLineWidth(4);  g_newage_ug->SetLineStyle(1);
  g_newage_ug->SetLineColor(2);
  g_newage_ug->Draw("C");    



  
  
  // MSMM theory
  // Ellis, Ferstl, Olive, PRD63, Fig5c tanbeta=10
  // int nPoints = 19; 
  //    float x_th[]={ 53,  60,  70,  79,   81, 100,  130,  135, 141,  156,  161,  171, 186,  194,  200};
  //     float y_th[]={0.1, 0.6, 2.5, 3.9, 4.02, 0.6, 0.36, 0.29, 0.3, 0.25, 0.25, 0.22, 0.12, 0.12, 0.1};
  // Fig 7 SD 
  int nPoints = 16; 
  float x_th[]={50.,   60.,  70.,  80.,   90, 100, 110, 120, 130, 140., 150, 160,   200,  300,  400, 500 };
  float y_th[]={0.015, 0.55, 3.,   4.5,  1.5, .65, .5, .4,   .35, .3,   0.2, 0.15, 0.12, 0.007, 0.004, 0.0025    };
  
  for (int i=0;i<nPoints;i++) y_th[i]*=1.e-36; 
  TGraph* g_th= new TGraph(nPoints, x_th, y_th);
  for (int i=0;i<g_th->GetN();i++) g_th->GetY()[i]*=1e-3; 
  g_th->SetFillColor(7);
  g_th->Draw("LF");
  //  l.SetTextAngle(0); l.SetTextSize(0.025); l.DrawLatex(70,1.5e-40,"MSSM");
  



  
  TLegend *legend=new TLegend(0.7,0.7,0.9,0.9);
  legend->SetFillColor(10);
  g_newage_ug->SetName("g_newage_ug");
  legend->AddEntry("g_newage_ug","NEWAGE Underground Run 2010","l");
  graphSD_dmtpc_surface->SetName("graphSD_dmtpc_surface");
  legend->AddEntry("graphSD_dmtpc_surface","DMTPC Surface Run 2009","l");
  graphSD_dmtpc_wipp->SetName("graphSD_dmtpc_wipp");
  legend->AddEntry("graphSD_dmtpc_wipp","DMTPC in WIPP, Projected","l");
  //  g_kims->SetName("g_kims");
  //  legend->AddEntry("g_kims","KIMS 2007","l");
  g_picasso->SetName("g_picasso");
  legend->AddEntry("g_picasso","PICASSO 2009","l");
  g_coupp->SetName("g_coupp");
  legend->AddEntry("g_coupp","COUPP 2010","l");
  graphSD_dmtpc_1m3wipp->SetName("graphSD_dmtpc_1m3wipp");
  legend->AddEntry("graphSD_dmtpc_1m3wipp","DMTPC 1 m^3 in WIPP, Projected","l");
  g_th->SetName("g_th");
  legend->AddEntry("g_th","MSSM","f");
  //  legend->Draw();

  // other experiments (worse sensitivity)
  //  graphSD_dmtpc_surface_bgndsub->SetName("graphSD_dmtpc_surface_bgndsub");
  //  legend->AddEntry("graphSD_dmtpc_surface_bgndsub","DMTPC Surface Run 2009, Bgnd. Subtracted","l");
  //  graphSD_newage_surface->SetName("graphSD_newage_surface");
  //  legend->AddEntry("graphSD_newage_surface","NEWAGE Surface Run 2007","l");
  //  g_naiad->SetName("g_naiad");
  //  legend->AddEntry("g_naiad","NAIAD 2005","l");
  //  g_coupp->SetName("g_coupp");
  //  legend->AddEntry("g_coupp","COUPP 2008","l");

  cout << " Done!" << endl; 
  myCanvas2->Update();



  // compare with published results:

//  cout << "CFR published results: " << endl; 
//  cout << "CDMS 2008 limit on SD_n 2.7 10-38 at 60 GeV " << endl; 
//  cout << "XENON 10 2007 limit on SD_n  6*10-39 at 30 GeV " << endl; 
//  cout << "XENON 10 2007 limit on SD_p  5?*10-37 at 30 GeV " << endl; 
//  cout << "CDMS 2008 limit on SD_n 2.7 10-38 at 60 GeV " << endl; 
//  cout << "NewAge surface run SD_p 2 10-32 at 120? GeV " << endl; 
//  cout << endl;

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
