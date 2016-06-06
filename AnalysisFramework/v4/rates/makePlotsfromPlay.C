#include <vector>
#include <iostream>
#include <fstream>
#include "../../MaxCam/DmtpcSkimEvent.hh"
#include "../../MaxCam/DmtpcSkimDataset.hh"
#include "../../MaxCam/DmtpcAstro.hh"
#include "../../MaxCam/MaxCamSRIM.hh"
#include "../../MaxCam/MaxCamImageTools.hh"
#include "../../MaxCam/waveformtools/include/CspWfVector.hh"
#include "../../MaxCam/waveformtools/include/CspPulse.hh"
#include "../../MaxCam/waveformtools/include/CspWaveform.hh"
#include "src/recoilEnergy.hh"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TTimeStamp.h"
#include "TRandom3.h"

void calculaterayleigh(int n, double* angles, bool axial,		       
		       double& Rbar2, double& phibar, double& Rstar)

{   
   double C=0;   
   double S=0;   
   for(int i=0; i<n; i++)   
   {      
      if(axial)      
      {	 
	 C+=cos(2.0*angles[i]);	 
	 S+=sin(2.0*angles[i]);      
      }      
      else      
      {	 
	 C+=cos(angles[i]);	 
	 S+=sin(angles[i]);      
      }   
   }   
   double R2 = pow(C,2)+pow(S,2);   
   double Cbar = C/double(n);   
   double Sbar = S/double(n);   
   Rbar2 = Cbar*Cbar+Sbar*Sbar;   
   phibar = atan2(Sbar,Cbar);   
   Rstar = (2*n-1)*Rbar2+n/2.0*Rbar2*Rbar2;
//    cout << C << "," << S << "," << R2 << ";" 
// 	<< Cbar << "," << Sbar << "," << Rbar2 << "," << phibar << "," 
// 	<< Rstar << endl;
}
void calculatekupier(int n, double* angles, double& V, bool iszeroto2pi)
{   
   const int npts=n;   
   double Dplus=0, Dminus=0;   
   int iangles[npts];   
   double sangles[npts];   
//double abscissa[npts];   
   TMath::Sort(npts,angles,iangles,kFALSE);   
   for(int i=0; i<npts; i++)   
   {      
// abscissa[i]=(i+1)/double(n+1);      
      sangles[i]=angles[iangles[i]]/(2*TMath::Pi());      
      if(iszeroto2pi==false)	 
	 sangles[i]+=0.5;      
      if((i+1)/double(n)-sangles[i]>Dplus)	 
	 Dplus=(i+1)/double(n)-sangles[i];      
      if(sangles[i]-i/double(n) > Dminus)	 
	 Dminus=sangles[i]-i/double(n);   
   }   
   V=(Dplus+Dminus)*(sqrt(double(npts))+0.155+0.24/sqrt(double(npts)));
}

double mod2pi(double phi)
{   
   while(phi>TMath::Pi() || phi<-TMath::Pi())   
   {      
      if(phi>TMath::Pi())	 
	 phi-=TMath::Pi()*2;      
      if(phi<-TMath::Pi())	 
	 phi+=TMath::Pi()*2;   
   }   
   return phi;
}

double modpi(double phi)
{   
   while(phi>TMath::Pi()/2.0 || phi<-TMath::Pi()/2.0)   
   {      
      if(phi>TMath::Pi()/2.0)	 
	 phi-=TMath::Pi();      
      if(phi<-TMath::Pi()/2.0)	 
	 phi+=TMath::Pi();   
   }   
   return phi;
}

double distToSpacer(DmtpcSkimEvent* ev, int c, int t, DmtpcGainMap* map)
{
   double min=1000;
   
   vector< int > px = ev->cluster(c)->getCluster(t);
      
   int nbinsx = ((TH2F*)ev->image(c))->GetNbinsX()+2;
   int nbinsy = ((TH2F*)ev->image(c))->GetNbinsY()+2;

   for(int i=0; i<px.size(); i++)
   {
      int bx = px[i] % nbinsx; 
      int by = ((px[i]-bx) / nbinsx) % nbinsy; 
      double x = ev->image(c)->GetXaxis()->GetBinCenter(bx);
      double y = ev->image(c)->GetYaxis()->GetBinCenter(by);
      int s;
      double dist = map->distanceToNearestSpacer(x,y,s); 
      if(dist<min) min=dist;
   }

   return min;
   
}


using namespace std;
using namespace RecoilEnergy;

void makePlotsfromPlay(TString playlist, TString tree="skim")
{

   double nang = -175.0*TMath::Pi()/180;
   ifstream play;
   play.open(playlist);

   vector< vector < int > > events;
   events.clear();
   TString temp;
   while(play >> temp)
   {
      vector<int> evdet;
      evdet.clear();
      int run,ev,cam,tr;
      play >> run >> ev >> cam >> tr;
      evdet.push_back(run);
      evdet.push_back(ev);
      evdet.push_back(cam);
      evdet.push_back(tr);

      events.push_back(evdet);
   }
   
   
   TString filepathbase = "/net/zwicky/dmtpc/analysis/10L/skim/dmtpc_10L_";
   TString passpathbase = "pass/dmtpc_10L_";
   TString origpathbase = "/net/zwicky/dmtpc/data/10L/dmtpc_10L_";
//    TString origpathbase = "/net/zwicky/dmtpc/akaboth/projects/DarkMatter/MaxCam/Simulations/v1/output/data/dmtpc_mc_00";
//    TString filepathbase = "skim/dmtpc_mc_00";
//    TString passpathbase = "pass/dmtpc_mc_00";

   DmtpcSkimDataset d;
   DmtpcDataset od;

   TCanvas* c0 = new TCanvas("c0","c0",0,0,1200,1000);
   c0->Divide(2,2);
   TCanvas* c1 = new TCanvas("c1","c1",0,0,1200,1000);
   c1->Divide(2,2);
   TCanvas* c2 = new TCanvas("cTop2","c2",0,0,1200,1000);
   c2->Divide(2,2);
   TCanvas* c3 = new TCanvas("c3","c3",0,0,1200,1000);
   c3->Divide(2,2);
   TCanvas* c4 = new TCanvas("c4","c4",0,0,1600,1000);
   c4->Divide(3,2);

   TCanvas* c5 = new TCanvas("c5","c5",0,0,1400,1000);
   c5->Divide(2,2);
   TCanvas* c6 = new TCanvas("c6","c6",0,0,1600,1000);
   c6->Divide(4,2);

   TTree* pass=0;
   TFile* passfile=0;

   TH2F* hXY[2];
   TH1F* hE[2];
   TH1F* hPhi[2];
   TH1F* hPhiAx[2];
   TH1F* hRedPhi[2];
   TH1F* hRedPhiAx[2];
   TH2F* hMoments[2];
   TH1F* hTMoment[2];
   TH1F* hRuns[2];
   TH1F* hTime[2];
   TH1F* hSpacer[2];
   TH2F* hRvE[2];
   TH1F* hPhiE[2][4];
   
   TH2F* hRvEAll = new TH2F("hRvEAll","hRvEAll",300,0,300,16,0,8);
   hRvEAll->GetXaxis()->SetTitle("Recoil Energy (keV)");
   hRvEAll->GetYaxis()->SetTitle("Range (mm)");
   
   TH1F* hRedPhiAll = new TH1F("hRedPhiAll","hRedPhiAll",20,-3.14159,3.14159);
   hRedPhiAll->GetXaxis()->SetTitle("#phi-#phi_{cygnus}");
   hRedPhiAll->GetYaxis()->SetTitle("N/(#pi/10 radians)");
   
   TH1F* hEAll = new TH1F("hEAll","hEAll",30,60,210);
   hEAll->GetXaxis()->SetTitle("Recoil Energy (keV)");
   hEAll->GetYaxis()->SetTitle("N/(5 keV)");

   TH2F* hPhivEAll = new TH2F("hPhivEAll","hPhivEAll",300,60,210,200,-3.14159,3.14159);
   hPhivEAll->GetXaxis()->SetTitle("Recoil Energy (keV)");
   hPhivEAll->GetYaxis()->SetTitle("#phi-#phi_{cygnus}");
   
   map<TString,int> cameraMap;

   vector< vector <double> > camphi;
   vector< vector <double> > rawphi;
   vector< vector <double> > rawphiax;
   vector< vector <double> > azimuth;
   vector< vector <double> > redphi;
   vector< vector <double> > redphiax;
   
   vector < double > redphiall;
   
   int n80200[2];
   
   TTimeStamp* starttime = new TTimeStamp(2011,05,01,17,32,59);
   TTimeStamp* endtime = new TTimeStamp(2011,06,25,0,47,48);
   

   for(int i=0; i<2; i++)
   {
      
      TString name = "hXY";
      name+=i;
      hXY[i]= new TH2F(name,name,1024,0,1024,1024,0,1024);

      name = "hE";
      name+=i;
      hE[i]=new TH1F(name,name,500,0,1000);
      hE[i]->GetXaxis()->SetTitle("Recoil Energy (keV)");
      hE[i]->GetYaxis()->SetTitle("Counts/3 keV");

      name = "hPhi";
      name+=i;
      hPhi[i]=new TH1F(name,name,20,-3.14159,3.14159);
      hPhi[i]->GetXaxis()->SetTitle("Detector #phi");

      name = "hPhiAx";
      name+=i;
      hPhiAx[i]=new TH1F(name,name,20,-3.14159,3.14159);

      name = "hRedPhi";
      name+=i;
      hRedPhi[i]=new TH1F(name,name,20,-3.14159,3.14159);
      hRedPhi[i]->GetXaxis()->SetTitle("#phi-#phi_{cygnus}");

      name = "hRedPhiAx";
      name+=i;
      hRedPhiAx[i]=new TH1F(name,name,20,-3.14159,3.14159);

      name = "hMoments";
      name+=i;
      hMoments[i] = new TH2F(name,name,100,0,300,100,0,300);

      name="hTMoment";
      name+=i;
      hTMoment[i] = new TH1F(name,name,150,0,15);

      name="hRuns";
      name+=i;
      hRuns[i] = new TH1F(name,name,200,events[0][0],events[events.size()-1][0]+1);

      name="hTime";
      name+=i;
      hTime[i] = new TH1F(name,name,48,starttime->AsDouble(),endtime->AsDouble());
      hTime[i]->GetXaxis()->SetTimeDisplay(1);
      hTime[i]->GetXaxis()->SetTimeFormat("%d-%m-%y%F1970-01-01 00:00:00");
      hTime[i]->GetXaxis()->SetTitle("Date of Event");

      name="hSpacer";
      name+=i;
      hSpacer[i]= new TH1F(name,name,200,0,200);

      name="hRvE";
      name+=i;
      hRvE[i]=new TH2F(name,name,120,0,1000,32,0,8);
      hRvE[i]->GetXaxis()->SetTitle("Recoil Energy (keV)");
      hRvE[i]->GetYaxis()->SetTitle("Projected Range (mm)");
      hRvE[i]->SetMarkerSize(0.4);

      vector <double> temp;
      camphi.push_back(temp);
      rawphi.push_back(temp);
      rawphiax.push_back(temp);
      azimuth.push_back(temp);
      redphi.push_back(temp);
      redphiax.push_back(temp);

      n80200[i]=0;

      for(int j=0; j<4; j++)
      {
	 name="hPhiE";
	 name+=i;
	 name+="_";
	 name+=j;

	 hPhiE[i][j] = new TH1F(name,name,20,-3.14159,3.14159);
      }

   }

   cout << redphi.size() << endl;

   DmtpcGainMap* maps[2];
   int lowtime = 1304607540;
   int hightime= 1308579600;
   TRandom3* rnd = new TRandom3();
//   TTimeStamp* time = new TTimeStamp(lowtime);

 
   for(int i=0; i<int(events.size()); i++)
   {
      if(i%20 == 0) cout << "N: " << i+1 << endl;
      
      TString filepath = filepathbase;
      if(events[i][0]<10000)
	 filepath+="0";
      filepath+=events[i][0];
      filepath+="skim.root";
      TString passpath = passpathbase;
      if(events[i][0]<10000)
	 passpath+="0";
      passpath+=events[i][0];
      passpath+="pass.root";
      TString origfilepath = origpathbase;
      if(events[i][0]<10000)
	 origfilepath+="0";
      origfilepath+=events[i][0];
      origfilepath+=".root";
      
      d.openRootFile(filepath);      
//       d.tree(tree)->SetBranchStatus("*trigger*",0);
//       d.tree(tree)->SetBranchStatus("*clusters*",0);
      d.getEvent(events[i][1]);
      if(i==0)
      {
	 maps[0]=(DmtpcGainMap*)d.getGainMap(0)->Clone("gm0");
	 maps[1]=(DmtpcGainMap*)d.getGainMap(1)->Clone("gm1");

	 c0->cd(1);
//	 maps[0]->drawWithSpacers();
	 c0->cd(2);
//	 maps[1]->drawWithSpacers();
	 c0->Update();
//	 getchar();
	 cameraMap[d.event()->cameraSerialNumber(0)]=0;
	 cameraMap[d.event()->cameraSerialNumber(1)]=1;

      }


      od.openRootFile(origfilepath);
      od.getEvent(events[i][1]);
      

      passfile = new TFile(passpath);

      pass = (TTree*)passfile->Get("cuts");	 
      double _Erecoil[2][15];
      double _calibratedrange[2][15];
      pass->SetBranchAddress("_Erecoil",&_Erecoil);
      pass->SetBranchAddress("_calibratedrange",&_calibratedrange);
      pass->GetEntry(events[i][1]);
      
      TString camSerNum = d.event()->cameraSerialNumber(events[i][2]);

      double distSpac = distToSpacer(d.event(),events[i][2],events[i][3],maps[cameraMap[camSerNum]]);

      hSpacer[cameraMap[camSerNum]]->Fill(distToSpacer(d.event(),events[i][2],events[i][3],maps[cameraMap[camSerNum]]));
      
      if(distSpac>=2 && _Erecoil[events[i][2]][events[i][3]] > 50 &&
		       _Erecoil[events[i][2]][events[i][3]] < 2000 
	 &&  _calibratedrange[events[i][2]][events[i][3]] < 60)
      {

	 if(_Erecoil[events[i][2]][events[i][3]] > 80 && _Erecoil[events[i][2]][events[i][3]] < 200)
	    n80200[cameraMap[camSerNum]]++;
	 
	 hXY[cameraMap[camSerNum]]->Fill(d.event()->x(events[i][2],events[i][3]),
					 d.event()->y(events[i][2],events[i][3]));
	 
	 hE[cameraMap[camSerNum]]->Fill(_Erecoil[events[i][2]][events[i][3]]);
	 
	 hPhi[cameraMap[camSerNum]]->Fill(d.event()->phi(events[i][2],events[i][3]));
	 if(_Erecoil[events[i][2]][events[i][3]]<100)
	    hPhiE[cameraMap[camSerNum]][0]->Fill(d.event()->phi(events[i][2],events[i][3]));
	 if(_Erecoil[events[i][2]][events[i][3]]>=100 && _Erecoil[events[i][2]][events[i][3]]<150)
	    hPhiE[cameraMap[camSerNum]][1]->Fill(d.event()->phi(events[i][2],events[i][3]));
	 if(_Erecoil[events[i][2]][events[i][3]]>=150 && _Erecoil[events[i][2]][events[i][3]]<200)
	    hPhiE[cameraMap[camSerNum]][2]->Fill(d.event()->phi(events[i][2],events[i][3]));
	 if(_Erecoil[events[i][2]][events[i][3]]>=200 && _Erecoil[events[i][2]][events[i][3]]<250)
	    hPhiE[cameraMap[camSerNum]][3]->Fill(d.event()->phi(events[i][2],events[i][3]));
	 
	 hMoments[cameraMap[camSerNum]]->Fill(d.event()->moment(events[i][2],2,events[i][3]),
					      d.event()->transverse_moment(events[i][2],2,events[i][3]));
	 
	 hTMoment[cameraMap[camSerNum]]->Fill(sqrt(d.event()->transverse_moment(events[i][2],2,events[i][3])));
	 
	 hRuns[cameraMap[camSerNum]]->Fill(events[i][0]);
	 
	 TTimeStamp* time = od.event()->UTCtimeStamp();
	 
//	 time->Set(int(rnd->Uniform(lowtime,hightime)),1,0,0);
	 hTime[cameraMap[camSerNum]]->Fill(time->AsDouble());

	 hRvE[cameraMap[camSerNum]]->Fill(_Erecoil[events[i][2]][events[i][3]],_calibratedrange[events[i][2]][events[i][3]]);

	 camphi[cameraMap[camSerNum]].push_back(d.event()->phi(events[i][2],events[i][3]));
	 if(camSerNum=="100439")
	    rawphi[cameraMap[camSerNum]].push_back(mod2pi(TMath::Pi()- camphi[cameraMap[camSerNum]].back()));
	 else
	    rawphi[cameraMap[camSerNum]].push_back(camphi[cameraMap[camSerNum]].back());
	 rawphiax[cameraMap[camSerNum]].push_back(2*modpi(rawphi[cameraMap[camSerNum]].back()));
	 hPhiAx[cameraMap[camSerNum]]->Fill(rawphiax[cameraMap[camSerNum]].back());
	 double eang = rawphi[cameraMap[camSerNum]].back()+nang;
	 azimuth[cameraMap[camSerNum]].push_back(mod2pi(TMath::Pi()/2.0-eang));
	 
	 UInt_t yyu,mmu,ddu,hhu,minu,ssu;
	 time->GetDate(1,0,&yyu,&mmu,&ddu);
	 time->GetTime(1,0,&hhu,&minu,&ssu);

	 int yy = int(yyu);
	 int mm = int(mmu);
	 int dd = int(ddu);
	 double hh = double(hhu);
	 double min = double(minu);
	 double ss = double(ssu);
	 double l,b;	 
	 DmtpcAstro::getWimpGalCoord(yy,mm,dd,hh+min/60.+ss/3600.,l,b);	 
	 b*=-1;	 
	 l+=180.0;
	 double cygnusaz = DmtpcAstro::getAzfromGal(l,b,32.3716667,-103.793611,yy,mm,dd,hh+min/60.+ss/3600.);
	 cygnusaz = mod2pi(cygnusaz*TMath::Pi()/180.0);
	 redphi[cameraMap[camSerNum]].push_back(mod2pi(azimuth[cameraMap[camSerNum]].back()-cygnusaz));
	 hRedPhi[cameraMap[camSerNum]]->Fill(redphi[cameraMap[camSerNum]].back());
	 redphiax[cameraMap[camSerNum]].push_back(modpi(redphi[cameraMap[camSerNum]].back()));
	 hRedPhiAx[cameraMap[camSerNum]]->Fill(redphiax[cameraMap[camSerNum]].back());

	 redphiall.push_back(redphi[cameraMap[camSerNum]].back());

	 hRvEAll->Fill(_Erecoil[events[i][2]][events[i][3]],_calibratedrange[events[i][2]][events[i][3]]);
	 hEAll->Fill(_Erecoil[events[i][2]][events[i][3]]);
	 hRedPhiAll->Fill(redphiall.back());
	 hPhivEAll->Fill(_Erecoil[events[i][2]][events[i][3]],redphiall.back());

      }
      
      passfile->Close();
      
   }
	    

   cout << "Between 80 and 200: " << n80200[0] << "," << n80200[1] << endl;

   c0->cd(1);
   hXY[0]->Draw();
   c0->cd(2);
   hXY[1]->Draw();

   c0->cd(3);
   hTMoment[0]->Draw();
   c0->cd(4);
   hTMoment[1]->Draw();

   c1->cd(1);
   hE[0]->Draw();
   c1->cd(2);
   hE[1]->Draw();
   
   c1->cd(3);
   hPhi[0]->Draw();
   c1->cd(4);
   hPhi[1]->Draw();

   c2->cd(1);
   hTime[0]->Draw();
   c2->cd(2);
   hTime[1]->Draw();

   c2->cd(3);
   hRuns[0]->Draw();
   c2->cd(4);
   hRuns[1]->Draw();

   c3->cd(1);
   hSpacer[0]->Draw();
   c3->cd(2);
   hSpacer[1]->Draw();

   MaxCamSRIM* nuclear = new MaxCamSRIM("SRIM_F_in_CF4_100Torr");
   nuclear->setPressure(60);
   TGraph* gNuclear = nuclear->getRangeVsEnergy();
   gNuclear->SetLineWidth(2);
   gNuclear->SetLineColor(kBlue);

   c3->cd(3);
   hRvE[0]->Draw();
   gNuclear->Draw("SAME");
   c3->cd(4);
   hRvE[1]->Draw();
   gNuclear->Draw("SAME");

   c4->cd(1);
   hPhiAx[0]->Draw();
   c4->cd(4);
   hPhiAx[1]->Draw();
   c4->cd(2);
   hRedPhi[0]->Draw();
   c4->cd(5);
   hRedPhi[1]->Draw();
   c4->cd(3);
   hRedPhiAx[0]->Draw();
   c4->cd(6);
   hRedPhiAx[1]->Draw();

   c5->cd(1);
   hRvEAll->SetMarkerSize(0.3);
   hRvEAll->Draw();
   gNuclear->Draw("SAME L");
   c5->cd(2);
   hEAll->Draw();
   c5->cd(3);
   hPhivEAll->SetMarkerSize(0.3);
   hPhivEAll->Draw();
   c5->cd(4);
   hRedPhiAll->Draw();


   double Rbar2, phibar, Rstar;

   cout << "Camera 0: " << endl;

   calculaterayleigh(camphi[0].size(),&rawphi[0][0],false,Rbar2,phibar,Rstar);      
   cout << "Raw Phi: " << endl << "\t Rbar2: " << Rbar2 << "\t phibar: " << phibar	
	<< "\t Rstar: " << Rstar << "\t Prob: " << TMath::Prob(Rstar,2) << endl;   
   calculaterayleigh(camphi[0].size(),&redphi[0][0],false,Rbar2,phibar,Rstar);      
   cout << "Reduced Phi: " << endl << "\t Rbar2: " << Rbar2 << "\t phibar: " << phibar	
	<< "\t Rstar: " << Rstar << "\t Prob: " << TMath::Prob(Rstar,2) << endl;   

   cout << "Camera 1: " << endl;

   calculaterayleigh(camphi[1].size(),&rawphi[1][0],false,Rbar2,phibar,Rstar);      
   cout << "Raw Phi: " << endl << "\t Rbar2: " << Rbar2 << "\t phibar: " << phibar	
	<< "\t Rstar: " << Rstar << "\t Prob: " << TMath::Prob(Rstar,2) << endl;   
   calculaterayleigh(camphi[1].size(),&redphi[1][0],false,Rbar2,phibar,Rstar);      
   cout << "Reduced Phi: " << endl << "\t Rbar2: " << Rbar2 << "\t phibar: " << phibar	
	<< "\t Rstar: " << Rstar << "\t Prob: " << TMath::Prob(Rstar,2) << endl;   

   cout << "All: " << endl;

   calculaterayleigh(redphiall.size(),&redphiall[0],false,Rbar2,phibar,Rstar);      
   cout << "Reduced Phi: " << endl << "\t Rbar2: " << Rbar2 << "\t phibar: " << phibar	
	<< "\t Rstar: " << Rstar << "\t Prob: " << TMath::Prob(Rstar,2) << endl;  

   c6->cd(1);
   hPhiE[0][0]->Draw();
   c6->cd(2);
   hPhiE[0][1]->Draw();
   c6->cd(3);
   hPhiE[0][2]->Draw();
   c6->cd(4);
   hPhiE[0][3]->Draw();
   c6->cd(5);
   hPhiE[1][0]->Draw();
   c6->cd(6);
   hPhiE[1][1]->Draw();
   c6->cd(7);
   hPhiE[1][2]->Draw();
   c6->cd(8);
   hPhiE[1][3]->Draw();
   c6->Update();
   
}
