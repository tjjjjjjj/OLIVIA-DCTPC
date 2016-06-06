#include "TString.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TApplication.h"
#include "TPluginManager.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TColor.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TTimeStamp.h"
#include "TRandom3.h"

#include "DmtpcSkimEvent.hh"
#include "DmtpcAstro.hh"
#include "MaxCamImageTools.hh"
#include "AnalysisCut.hh"
#include "MaxCamSRIM.hh"
#include "AnalysisConfig.hh"
#include "recoilEnergy.hh"

#include <vector>
#include <iostream>
#include <fstream> 
#include <map>


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

bool checkEdgeSpark(DmtpcSkimEvent* ev, int c, int t)
{
   int nedgepx=0;
   int nbinsx = ev->image(c)->GetNbinsX()+2;
   for(int j=0; j<ev->ntracks(c); j++)
   {
      if(ev->x(c,j)<8)
      {
	 vector< vector< int > > px = ev->cluster(c)->getClusters();
	 for(int k=0; k<px.size(); k++)
	 {
	    for(int m=0; m<px[k].size(); m++)
	    {
	       if(px[k][m]%nbinsx==1)
		  nedgepx++;
	    }
	 }
	 break;
      }
      
   }
   if(nedgepx>40) return 1;
   else return 0;


}


static bool should_collinear_merge(double min_global_rxy,  double min_each_rxy, 
                     double max_join_residual, double weight, 
                     const vector<int> c1, const vector<int> c2, const TH2F * image)
{

     std::cout << "Starting collinear merge!!!" << std::endl; 
     double uxy, ux, uxx,uy,uyy;        
     uxy=0;ux=0;uxx=0; uy=0;uyy=0; 
     double n=0;
     double rxy_each[2];  
     double m[2]; 
     double b[2]; 
     double ctr_x[2]; 
     double ctr_y[2]; 

     for (int clust = 0; clust < 2; clust++)
     {
       vector<int> which = clust == 0 ? c1 : c2; 

       double this_uxy, this_ux, this_uxx,this_uy,this_uyy;        
       this_uxy=0;this_ux=0;this_uxx=0; this_uy=0;this_uyy=0; 
       double this_n= 0; 
       for (unsigned int bin = 0; bin < which.size();bin++)
       {
          double content = image->GetArray()[which[bin]]; 
          content = pow(content,weight); 
          this_n+=content; 
          int x,y,z; 
          image->GetBinXYZ(which[bin],x,y,z); 
          this_ux+=content*x; 
          this_uy+=content*y; 
          this_uxx+=content*x*x; 
          this_uyy+=content*y*y;
          this_uxy+=content*x*y; 
       }

       //TODO: There are some simplications that can be made here
       m[clust] = (this_uxy - 1./this_n*this_ux*this_uy) / (this_uxx - 1./this_n*this_ux*this_ux); 
       rxy_each[clust] = (this_n*this_uxy - this_ux*this_uy) / TMath::Sqrt( (this_n*this_uxx - this_ux*this_ux) * (this_n*this_uyy - this_uy*this_uy) ); 
       b[clust] = this_uy/this_n - m[clust]*this_ux/this_n; 

       std::cout << "m,b,rxy: " << m[clust] << " , " << b[clust] << " , " << rxy_each[clust] << std::endl; 
       ctr_x[clust] = this_ux/this_n; 
       ctr_y[clust] = this_uy/this_n; 

       n += this_n; 
       ux += this_ux/this_n; 
       uy += this_uy/this_n; 
       uyy += this_uyy/this_n; 
       uxx += this_uxx/this_n; 
       uxy += this_uxy/this_n; 

     }

     double global_rxy = (n*uxy - ux*uy) / TMath::Sqrt( (n*uxx - ux*ux) * (n*uyy - uy*uy) ); 


     //There are three cases here:
     //
     //  1) Both clusters have a well defined direction
     //     In this case, see if they match up well. 
     //  2) One cluster has a well defined direction. 
     //     In this case, see if the not so well defined cluster matches up well with the well defined cluster
     //  3) Neither cluster has a well defined direction. In this case, use the global line fit to determine if it's a good match. 
    
    if (fabs(rxy_each[0]) > min_each_rxy && fabs(rxy_each[1]) > min_each_rxy)
    {
        //Test residual for each each track extended to the other 
        //to see if it matches up    
       
        double resid1 = fabs (ctr_y[0] - m[1] * ctr_x[0] - b[1]) / TMath::Sqrt(1 + m[1]*m[1]); 
        double resid2 = fabs (ctr_y[1] - m[0] * ctr_x[1] - b[0]) / TMath::Sqrt(1 + m[0]*m[0]); 

        double resid_mean = TMath::Sqrt(resid1*resid2);  

        std::cout << "BOTH HAVE DIRECTION, RESIDUALS: " << resid2 << " , " << resid1 << "MEAN: " << resid_mean <<  std::endl; 

        return resid_mean < max_join_residual; 
    }


    if (fabs(rxy_each[0]) > min_each_rxy) 
    {
       double resid = fabs (ctr_y[1] - m[0] * ctr_x[1] - b[0]) / TMath::Sqrt(1 + m[0]*m[0]); 
       std::cout << "FIRST HAS DIRECTION, RESIDUAL: " << resid << std::endl; 
       return resid < max_join_residual; 
    } 

    if (fabs(rxy_each[1]) > min_each_rxy) 
    {
       double resid = fabs (ctr_y[0] - m[1] * ctr_x[0] - b[1]) / TMath::Sqrt(1 + m[1]*m[1]); 
       std::cout << "SECOND HAS DIRECTION, RESIDUAL: " << resid << std::endl; 
       return resid < max_join_residual; 
    } 
    
    std::cout << "NEITHER HAS DIRECTION, global_rxy: " << global_rxy << std::endl; 
    return fabs(global_rxy) > min_global_rxy; 
}


bool MergeClusters(const TH2F* image, 
		   vector<vector<int> > & allcluster,
		   int t, 
		   double minDistance, 
		   double minrxy_global, 
		   double minrxy_cluster, 
		   double maxjoin_residual, 
		   double ls_weight, 
		   DmtpcGainMap * map,
		   double spacer_width = 2)
{
   
   if (allcluster.size() <=1) return false; 

  int nbinsx = image->GetNbinsX()+2;
  int nbinsy = image->GetNbinsY()+2;

  vector<vector<double> > dists(allcluster.size(), vector<double>(map ? map->getNSpacers() : 1 ));  
  vector<vector<pair<double,double> > > positions(allcluster.size(), vector<pair<double,double> >(map ? map->getNSpacers() : 1));  
  vector<vector<int> > npositions(allcluster.size(), vector<int>(map ? map->getNSpacers() : 1)); 

  if (map && map->getNSpacers())
  {
    for (unsigned int i = 0; i < allcluster.size(); i++)
    {
      for (int s = 0; s < map->getNSpacers(); s++)
      {
        dists[i][s]=1e10; 
        positions[i][s].first=0; 
        positions[i][s].second=0; 
        npositions[i][s]=0; 
      }

      for (unsigned int b=0; b < allcluster[i].size(); b++)
      {
        int bin = allcluster[i][b]; 
       // std::cout << "LOOKING AT BIN " << b << std::endl; 
        int bx = bin % nbinsx; 
        int by = ((bin -bx) / nbinsx) % nbinsy; 
        double x = image->GetXaxis()->GetBinCenter(bx);
        double y = image->GetYaxis()->GetBinCenter(by);
       // std::cout << "x,y: " << x << "," << y << std::endl; 
        int s;
        double dist = map->distanceToNearestSpacer(x,y,s); 
         positions[i][s].first+=x;  
         positions[i][s].second+=y;  
         npositions[i][s]++; 
       // std::cout << "Nearest Spacer: " << s << " at a distance of " << dist << std::endl; 
        if (dist < dists[i][s])
        {
          dists[i][s] = dist;
        }
      }

      for (int s = 0; s < map->getNSpacers(); s++)
      {
        positions[i][s].first /= npositions[i][s];  
        positions[i][s].second /= npositions[i][s];  
      }

    }
  }


  for (unsigned int i = 0; i < allcluster.size(); i++)
  {
   
     if (int(i) == t) continue; 
     std::cout << " CHECKING " << i << " AGAINST " << t << std::endl; 
     double diag = TMath::Sqrt(2)*image->GetBinWidth(1);
     //Check for the following condition
     // i is close to a spacer
     // j is close to the same spacer
     // i and j are on opposite sides of the spacer
     // i and j are close together
     
     double mind = MaxCamImageTools::minDist(allcluster[i],allcluster[t],const_cast<TH2F*>(image)) ;
     
     cout << "MIN DIST: " << mind << endl;
     
     if (mind == 0 ) 
     {
        continue; 
     }
     
     bool merged = false;
     bool merges = false;
     
     if (mind < minDistance) 
     {
	std::cout << "DISTANCE BETWEEN " << i << " AND " << t << " is, " << mind << ", BELOW THRESH" << std::endl; 
	merged = true; 
     }
     if (map) 
     {
        for (int s = 0; s < map->getNSpacers(); s++)
        {
	   
	   if (dists[i][s] < spacer_width*map->getSpacerWidth(s) + diag 
	       && dists[t][s] < spacer_width*map->getSpacerWidth(s) + diag 
	       && map->crossesSpacer(positions[i][s].first,
				     positions[i][s].second,
				     positions[t][s].first,
				     positions[t][s].second)
	      )
	   {
	      
	      std::cout << "MINIMUM DISTANCE FROM CLUSTER " << i << " TO SPACER " << s << " IS " << dists[i][s] << std::endl;
	      std::cout << "MINIMUM DISTANCE FROM CLUSTER " << t << " TO SPACER " << s << " IS " << dists[t][s] << std::endl;
	      cout << "MERGING ACROSS SPACER" << endl;
	      
	      merges = true;  
	      break; //We've already considered merging these two clusters, so no point in continuing 
	   }
        }
     }
     
     if ((merged && merges) || ((merged || merges) && (mind <= 2 * TMath::Power(image->GetBinWidth(0),2)  || should_collinear_merge(minrxy_global, minrxy_cluster, maxjoin_residual, ls_weight, allcluster[i], allcluster[t],image))))
     {
        std::cout << "MERGING " << i << " WITH " << t << std::endl; 
	return true;
	
     }
  }

  return false;
}



bool cutoffCheck(DmtpcSkimEvent* ev, int c, int t)
{
//    TCanvas* dd = new TCanvas("dd","dd",0,0,1600,1000);
//    dd->Divide(3,2);
   
   TH2F* img = (TH2F*)ev->image(c);

   TH2F* clusthist =  MaxCamImageTools::makeClusterHist(img,ev->cluster(c)->getCluster(t));
   
   TH2F* clusthist_red = (TH2F*)clusthist->Clone("clusthist_red");

   for(int j=1; j<=clusthist->GetNbinsX(); j++)
   {
      for(int k=0; k<=clusthist->GetNbinsY(); k++)
      {
	 float xcent = clusthist->GetXaxis()->GetBinCenter(j);
	 float ycent = clusthist->GetYaxis()->GetBinCenter(k);
	 
	 int xbin = img->GetXaxis()->FindBin(xcent);
	 int ybin = img->GetYaxis()->FindBin(ycent);
	 
	 int bin = img->GetBin(xbin,ybin);
	 
	 if(!ev->cluster(c)->isInCluster(t,bin))
	    clusthist_red->SetBinContent(j,k,0);
	 
      }
   }

//    dd->cd(1);
//    clusthist->Draw("colz");
//    dd->cd(4);
//    clusthist_red->Draw("colz");
   
   TH1D* clusthist_x = clusthist->ProjectionX();
   clusthist_x->Scale(1/double(clusthist_x->GetMaximum()));
   TH1D* clusthist_x_red = clusthist_red->ProjectionX();
   clusthist_x_red->Scale(1/double(clusthist_x_red->GetMaximum()));

//    dd->cd(2);
//    clusthist_x->Draw();
//    dd->cd(5);
//    clusthist_x_red->Draw();
      
   TH1D* clusthist_x_d = (TH1D*)clusthist_x->Clone("clusthist_x_d");
   TH1D* clusthist_x_red_d = (TH1D*)clusthist_x->Clone("clusthist_x_red_d");
   clusthist_x_d->SetBinContent(1,0);
   clusthist_x_red_d->SetBinContent(1,0);
   for(int xb=2; xb<=clusthist_x_d->GetNbinsX()-1; xb++)
   {
      clusthist_x_d->SetBinContent(xb,(clusthist_x->GetBinContent(xb+1)-clusthist_x->GetBinContent(xb))/(clusthist_x->GetBinCenter(xb+1)-clusthist_x->GetBinCenter(xb)));
      clusthist_x_red_d->SetBinContent(xb,(clusthist_x_red->GetBinContent(xb+1)-clusthist_x_red->GetBinContent(xb))/(clusthist_x_red->GetBinCenter(xb+1)-clusthist_x_red->GetBinCenter(xb)));
   }
   
   clusthist_x_d->SetBinContent(clusthist_x_d->GetNbinsX(),0);
   clusthist_x_red_d->SetBinContent(clusthist_x_red_d->GetNbinsX(),0);

//    dd->cd(3);
//    clusthist_x_d->Draw();
//    dd->cd(6);
//    clusthist_x_red_d->Draw();
   
   double derivmax = max(clusthist_x_d->GetMaximum(),fabs(clusthist_x_d->GetMinimum()));
   double derivredmax = max(clusthist_x_red_d->GetMaximum(),fabs(clusthist_x_red_d->GetMinimum()));

//    dd->Update();
//    getchar();

   delete clusthist;

   if(-1*derivmax-derivredmax+0.3<0)
      return 1;
   else
      return 0;

}

int main(int argc, char *argv[])
{
   const Int_t NRGBs = 7;
   const Int_t NCont = 104;
   
   Double_t stops[NRGBs] = { 0.00, 0.10, 0.30, 0.50, 0.65, 0.90, 1.00 };
   Double_t red[NRGBs]   = { 0.31, 0.00, 0.00, 0.10, 0.87, 1.00, 0.51 };
   Double_t green[NRGBs] = { 0.00, 0.00, 0.81, 0.90, 1.00, 0.20, 0.00 };
   Double_t blue[NRGBs]  = { 0.48, 0.51, 1.00, 0.10, 0.12, 0.00, 0.00 };
   
   TStyle *st1 = new TStyle("st1","my style");
   //st1->SetPalette(1);
   TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
   st1->SetNumberContours(NCont);
   st1->SetPadColor(0);
   st1->SetCanvasColor(0);
   st1->SetOptStat(0);

   st1->SetFrameBorderMode(0);
   st1->SetCanvasBorderMode(0);
   st1->SetPadBorderMode(0);
   st1->SetPadColor(0);
   st1->SetCanvasColor(0);
   st1->SetTitleColor(0);
   st1->SetTitleFont(42);
   st1->SetStatColor(0);
   //st1->SetFillColor(0);
   st1->SetTitleColor(1);
   
   // set the paper & margin sizes
   //st1->SetPaperSize(20,26);
   st1->SetPadTopMargin(0.05);
   st1->SetPadRightMargin(0.2); // 0.2
   st1->SetPadBottomMargin(0.2);
   st1->SetPadLeftMargin(0.2);
   
   st1->SetNdivisions(505,"x");
   st1->SetNdivisions(505,"y");
   
   // use large Times-Roman fonts
   st1->SetTextFont(42);
   st1->SetTextSize(0.08);
   
   st1->SetLabelFont(42,"x");
   st1->SetLabelFont(42,"y");
   st1->SetLabelFont(42,"z");
   st1->SetTitleFont(42,"x");
   st1->SetTitleFont(42,"y");
   st1->SetTitleFont(42,"z");
   
   st1->SetTitleOffset(1.25, "x");
   st1->SetTitleOffset(1.25, "y");
   st1->SetTitleOffset(1.25, "z");
   
   st1->SetLabelSize(0.03,"x");
   st1->SetTitleSize(0.04,"x");
   st1->SetLabelSize(0.03,"y");
   st1->SetTitleSize(0.04,"y");
   st1->SetLabelSize(0.03,"z");
   st1->SetTitleSize(0.04,"z");
   
   // stat box
   st1->SetStatBorderSize(1);
   st1->SetStatX(0.95);
   st1->SetStatY(0.95);
   st1->SetStatW(.200);
   st1->SetStatH(.125);
   st1->SetStatColor(0);
   st1->SetStatStyle(0);
   st1->SetStatFont(42);
   
   // title
   st1->SetTitleX(0.3);
   st1->SetTitleW(0.5);
   
   
   // use bold lines and markers
   st1->SetMarkerStyle(20);

   st1->SetHistLineWidth(1.85);
   st1->SetLineStyleString(2,"[12 12]"); // postscript dashes
   
   // do not display any of the standard histogram decorations
   st1->SetOptTitle(0);
   st1->SetOptStat(0);
   st1->SetOptFit(0);
   
   // put tick marks on top and RHS of plots
   st1->SetPadTickX(1);
   st1->SetPadTickY(1);
   
   st1->cd(); 
   
   TApplication theApp("App",&argc,argv);

   gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
					 "*", "TStreamerInfo",
					 "RIO","TStreamerInfo()");

   TString flistfile = theApp.Argv(1);
   vector< TString > runnumbers;
   ifstream flist;
   flist.open(flistfile);
   TString temp;
   while(flist >> temp)
   {
      runnumbers.push_back(temp);
   }

   AnalysisConfig* cfg;
   if(theApp.Argc()>2)
      cfg = new AnalysisConfig(theApp.Argv(2));
   else
      cfg = new AnalysisConfig("cfg/AnalysisDefault.cfg");

   TChain* cuts= new TChain("cuts","cuts");
   TChain* skim= new TChain(cfg->getTree(),cfg->getTree());
   TChain* sim = new TChain("Simulation","Simulation");
   skim->SetBranchStatus("*trigger*",0);
   skim->SetBranchStatus("*clusters*",1);

   for(int i=0; i<int(runnumbers.size()); i++)
   {
//      TString filename = "pass/dmtpc_10L_";
      TString filename = "pass/dmtpc_mc_";
      filename += runnumbers[i];
      filename += "_mixedpass.root";
      int valid = cuts->AddFile(filename, -1);
      cout << runnumbers[i] << "," << valid << "\t";
      if(valid==1)
      {
//	 TString skimfile = "/net/zwicky/dmtpc/analysis/10L/skim/dmtpc_10L_";
//	 TString skimfile = "/net/termite/raid1/dmtpc/analysis/10L/skim/dmtpc_10L_";
	 TString skimfile = "skim/dmtpc_mc_";
	 skimfile+=runnumbers[i];
	 skimfile+="_mixedskim.root";
	 skim->AddFile(skimfile,-1);
	 TString simfile = "/net/zwicky/dmtpc/akaboth/projects/DarkMatter/MaxCam/Simulations/v1/output/data/dmtpc_mc_";
	 simfile+=runnumbers[i];
	 simfile+="_mixed.root";
	 sim->AddFile(simfile,-1);
      }
   }
   cout << endl;
   
   cuts->AddFriend(skim);
   cuts->AddFriend(sim);

   int ncamera=0;
   cuts->SetBranchAddress("ncamera",&ncamera);
   cuts->GetEntry(0);
   const int ncam = ncamera;

   AnalysisCut* spark=0;
   AnalysisCut* nonZero=0;
   AnalysisCut* region=0;
   AnalysisCut* range=0;
   AnalysisCut* worm=0;
   AnalysisCut* edge=0;
   AnalysisCut* basicWorm=0;
   AnalysisCut* rbi=0;
   AnalysisCut* charge=0;
   double _Erecoil[ncam][15];
   double _wormclass[ncam][15];
   DmtpcSkimEvent* ev=0;
   float E;
   float x,y,phi;
   cuts->SetBranchAddress("spark",&spark);
   cuts->SetBranchAddress("nonZero",&nonZero);
   cuts->SetBranchAddress("region",&region);
   cuts->SetBranchAddress("range",&range);
   cuts->SetBranchAddress("worm",&worm);
   cuts->SetBranchAddress("edge",&edge);
   cuts->SetBranchAddress("basicWorm",&basicWorm);
   cuts->SetBranchAddress("rbi",&rbi);
   cuts->SetBranchAddress("charge",&charge);
   cuts->SetBranchAddress("_Erecoil",_Erecoil);
   cuts->SetBranchAddress("_wormclass",_wormclass);
   cuts->SetBranchAddress("event",&ev);
   sim->SetBranchAddress("E",&E);
   sim->SetBranchAddress("x",&x);
   sim->SetBranchAddress("y",&y);
   sim->SetBranchAddress("phi",&phi);

   int nevents = cuts->GetEntries();
   TH1F* hMC[ncam];
   TH1F* hAll[ncam];
   TH1F* hNonZero[ncam];
   TH1F* hRBI[ncam];
   TH1F* hWorm[ncam];
   TH1F* hBasicWorm[ncam];
   TH1F* hRegion[ncam];
   TH1F* hWormClass[ncam];
   TH1F* hCutoff[ncam];
   TH1F* hRejoin[ncam];
   TH1F* hEdgeSpark[ncam];
   TH1F* hSpacer[ncam];
   TH1F* hMoments[ncam];
   TH2F* hPhi[ncam];
   TH2F* hMC2[ncam];
   TH2F* hRedPhi[ncam];
   TH2F* hMC3[ncam];
   double rangecal[ncam];
   double Ecal[ncam];
   cuts->GetEntry(0);
   map<TString,int> cameraMap;
   double mcregion[ncam];
   double wormcutoff[ncam];


   for(int i=0; i<ncam; i++)
   {
      cameraMap[ev->cameraSerialNumber(i)]=i;
      rangecal[i]=cfg->getRangeCal(ev->cameraSerialNumber(i));
      Ecal[i]=cfg->getEnergyCal(ev->cameraSerialNumber(i));
      mcregion[i] = 512*rangecal[i];
      wormcutoff[i]=cfg->getWormCutVal(ev->cameraSerialNumber(i));
      
      hMC[i] = new TH1F("hMC","hMC",120,0,400);
      hAll[i] = new TH1F("hAll","hAll",120,0,400);
      hAll[i]->GetXaxis()->SetTitle("Recoil Energy (keV)");
      hAll[i]->GetYaxis()->SetTitle("Efficiency");
      hNonZero[i] = new TH1F("hNonZero","hNonZero",120,0,400);
      hNonZero[i]->SetLineColor(kRed);
      hRBI[i] = new TH1F("hRBI","hRBI",120,0,400);
      hRBI[i]->SetLineColor(kBlue);
      hWorm[i] = new TH1F("hWorm","hWorm",120,0,400);
      hWorm[i]->SetLineColor(kGreen);
      hBasicWorm[i] = new TH1F("hBasicWorm","hBasicWorm",120,0,400);
      hBasicWorm[i]->SetLineColor(kCyan);
      hRegion[i] = new TH1F("hRegion","hRegion",120,0,400);
      hRegion[i]->SetLineColor(kMagenta);
      hCutoff[i] = new TH1F("hCutoff","hCutoff",120,0,400);
      hCutoff[i]->SetLineColor(kBlue);
      hRejoin[i] = new TH1F("hRejoin","hRejoin",120,0,400);
      hRejoin[i]->SetLineColor(kYellow);
      hEdgeSpark[i] = new TH1F("hEdgeSpark","hEdgeSpark",120,0,400);
      hEdgeSpark[i]->SetLineColor(kBlue-10);
      hSpacer[i] = new TH1F("hSpacer","hSpacer",120,0,400);
      hSpacer[i]->SetLineColor(kCyan);
      hMoments[i] = new TH1F("hEdgeSpark","hEdgeSpark",120,0,400);
      hMoments[i]->SetLineColor(kOrange);
      hWormClass[i] = new TH1F("hWormClass","hWormClass",200,-1,1.5);
      
      hMC2[i] = new TH2F("hMC2","hMC2",20,0,400,10,-3.14159,3.14159);
      hPhi[i] = new TH2F("hPhi","hPhi",20,0,400,10,-3.14159,3.14159);
      hMC3[i] = new TH2F("hMC3","hMC3",20,0,400,10,-3.14159,3.14159);
      hRedPhi[i] = new TH2F("hRedPhi","hRedPhi",20,0,400,10,-3.14159,3.14159);
      
   }

   TFile* initfile = skim->GetCurrentFile();
   TObjArray* gainMaps = (TObjArray*)initfile->Get("gainMaps");
   int lowtime = 1304607540;
   int hightime= 1308579600;
   TRandom3* rnd = new TRandom3();
   TTimeStamp* time = new TTimeStamp(lowtime);

   for(int i=0; i<cuts->GetEntries(); i++)
   {
      if(i%1000==999) cout << "Event: " << i+1 << endl;
      cuts->GetEntry(i); 
      skim->GetEntry(i);
      sim->GetEntry(i);
      for(int c=0; c<ncam; c++)
      {
	 if(x>-1*mcregion[c] && x<mcregion[c] && 
	    y>-1*mcregion[c] && y<mcregion[c])
	 {
	    hMC[c]->Fill(E);
	    time->Set(int(rnd->Uniform(lowtime,hightime)),1,0,0);
	    //for MC, orientation angle is zero
	    double azimuth = TMath::Pi()/2.0-phi;
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
	    double redphi = mod2pi(azimuth-cygnusaz);	   
	    hMC2[c]->Fill(E,phi);
	    hMC3[c]->Fill(E,redphi);
	    for(int t=0; t<15; t++)
	    {
	       double xtrack = (ev->x(c,t)-512)*rangecal[c];
	       double ytrack = (ev->y(c,t)-512)*rangecal[c];
	       if(t<ev->ntracks(c) && fabs(x-xtrack)<10 && fabs(y-ytrack)<10)
	       {
		  hAll[c]->Fill(E);
	       
		  if(nonZero->passes(c,t))
		     hNonZero[c]->Fill(E);
		  if(nonZero->passes(c,t) && rbi->passes(c,t))
		     hRBI[c]->Fill(E);
		  if(nonZero->passes(c,t) && rbi->passes(c,t) 
		     && ev->maxpixel(c,t)<400 && ev->cluster_rms(c,t)>10.5 
		     && ev->maxpixel(c,t)/ev->E(c,t)<0.225)
		  {
		     hBasicWorm[c]->Fill(E);
		     hWormClass[c]->Fill(_wormclass[c][t]);
		  }
		  if(nonZero->passes(c,t) && rbi->passes(c,t)
		     && ev->maxpixel(c,t)<400 && ev->cluster_rms(c,t)>10.5 
		     && ev->maxpixel(c,t)/ev->E(c,t)<0.225 
		     &&_wormclass[c][t] > wormcutoff[c] )
		     hWorm[c]->Fill(E);
		  if(nonZero->passes(c,t) && rbi->passes(c,t) 
		     && ev->maxpixel(c,t)<400 && ev->cluster_rms(c,t)>10.5 
		     && ev->maxpixel(c,t)/ev->E(c,t)<0.225 
		     &&_wormclass[c][t] > wormcutoff[c]
		     && edge->passes(c,t) && region->passes(c,t))
		  {
		     hRegion[c]->Fill(E);

		           
		     if(!cutoffCheck(ev,c,t))
		     {
			hCutoff[c]->Fill(E);
			
			vector< vector < int> > px;
			px.clear();
			px = ev->cluster(c)->getClusters();
			
			if(!MergeClusters(ev->image(c),px,t,4000,0.65,0.8,10,2,(DmtpcGainMap*)gainMaps->At(c)))
			{
			   hRejoin[c]->Fill(E);
			   
			   if(!checkEdgeSpark(ev,c,t))
			   {
			      hEdgeSpark[c]->Fill(E);
			      if(distToSpacer(ev,c,t,(DmtpcGainMap*)gainMaps->At(c))>=0)
			      {
				 hSpacer[c]->Fill(E);
				 hRedPhi[c]->Fill(E,redphi);
				 hPhi[c]->Fill(E,phi);
			      }
			   }
			}
		     }
		  }
		  if(nonZero->passes(c,t) && rbi->passes(c,t) && ev->cluster_rms(c,t) > 10
		     && _wormclass[c][t]> wormcutoff[c] && basicWorm->passes(c,t)
		     && edge->passes(c,t) && region->passes(c,t)
		     && (ev->moment(c,2,t) > 2*ev->transverse_moment(c,2,t)))
		     hMoments[c]->Fill(E);
	       }
	       
	    }
	 }
      }
   }

   for(int i=0; i<ncam; i++)
   {
      hAll[i]->Divide(hMC[i]);
      hAll[i]->SetMaximum(1.1);
      
      hNonZero[i]->Divide(hMC[i]);
      hRBI[i]->Divide(hMC[i]);
      hWorm[i]->Divide(hMC[i]);
      hBasicWorm[i]->Divide(hMC[i]);
      hRegion[i]->Divide(hMC[i]);
      hCutoff[i]->Divide(hMC[i]);
      hRejoin[i]->Divide(hMC[i]);
      hEdgeSpark[i]->Divide(hMC[i]);
      hSpacer[i]->Divide(hMC[i]);
      hMoments[i]->Divide(hMC[i]);
      hPhi[i]->Divide(hMC2[i]);
      hRedPhi[i]->Divide(hMC3[i]);
   }

   TCanvas* c0 = new TCanvas("c0","c0",0,0,1250,800);
//   c0->Divide(2,1);

   c0->cd();
   hAll[0]->Draw();
   hNonZero[0]->Draw("SAME");
//   hRBI[0]->Draw("SAME");
//   hBasicWorm[0]->Draw("SAME");
   hWorm[0]->Draw("SAME");
   hRegion[0]->Draw("SAME");
   hCutoff[0]->Draw("SAME");
//   hRejoin[0]->Draw("SAME");
//   hEdgeSpark[0]->Draw("SAME");
   hSpacer[0]->Draw("SAME");
   c0->Update();
   c0->SaveAs("Eeff.eps");
   
//    c0->cd(2);
//    hAll[1]->Draw();
//    hNonZero[1]->Draw("SAME");
//    hRBI[1]->Draw("SAME");
//    hBasicWorm[1]->Draw("SAME");
//    hWorm[1]->Draw("SAME");
//    hRegion[1]->Draw("SAME");
// //   hMoments[1]->Draw("SAME");
//    hCutoff[1]->Draw("SAME");
//    hRejoin[1]->Draw("SAME");
//    hEdgeSpark[1]->Draw("SAME");

   c0->Update();

   TCanvas* c1 = new TCanvas("c1","c1",0,0,1200,800);
//   c1->Divide(2,1);
   c1->cd();
   hPhi[0]->Draw("colz");
//    c1->cd(2);
//    hPhi[1]->Draw("colz");
    c1->Update();
    c1->SaveAs("PhiEeff.eps");


   TCanvas* c2 = new TCanvas("c2","c2",0,0,1200,800);
//   c2->Divide(2,1);
   c2->cd();
   TH1D* hPy0 =hPhi[0]->ProjectionY("hPy0");
   hPy0->Scale(1/double(hPhi[0]->GetNbinsX()));
   hPy0->SetMinimum(0.3);
   hPy0->Draw();
   TH1D* hRPy0 = hRedPhi[0]->ProjectionY("hRPy0");
   hRPy0->Scale(1/double(hRedPhi[0]->GetNbinsX()));
   hRPy0->SetMinimum(0);
   hRPy0->SetLineColor(kRed);
   hRPy0->Draw("SAME");
//    c2->cd(2);
//    TH1D* hPy1 = hPhi[1]->ProjectionY("hPy1");
// //   hPy1->Scale(1/hPhi[1]->GetNbinsX());
//    hPy1->Draw();
   c2->Update();
   c2->SaveAs("Phieff.eps");

//   TFile* file = new TFile("efficiencies.root","UPDATE");
//   hSpacer[0]->Write(ev->cameraSerialNumber(0));
//   file->Close();
   

   getchar();

   return 0;

}
