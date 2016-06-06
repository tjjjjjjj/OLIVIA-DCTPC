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

#include "DmtpcSkimEvent.hh"
#include "MaxCamImageTools.hh"
#include "AnalysisCut.hh"
#include "MaxCamSRIM.hh"
#include "AnalysisConfig.hh"
#include "recoilEnergy.hh"

#include "DmtpcPulse.hh"
#include "DmtpcWaveformTools.hh"
#include "FastWfVector.hh"
#include "CspWfVector.hh"
#include "FastPulse.hh"
#include "CspPulse.hh"
#include "FastWaveform.hh"
#include "CspWaveform.hh"

#include <vector>
#include <iostream>
#include <fstream> 
#include <map>


bool chargeCuts10L(DmtpcSkimEvent* ev, int cam, int ntrack,double * param)
{  
//Cam == 0: top  
//Cam == 1: bottom  
   CspWfVector* top = (CspWfVector*) ev->waveform_vectors()->At(0);  
   CspWfVector* bot = (CspWfVector*) ev->waveform_vectors()->At(1);  
   int ntrig = top->size();  
   if (ntrig>=25) return false;  
   for (int trig = 0; trig < ntrig; trig++)
   {    
      bool pass = true;    
      if (fabs(cam-param[0])<0.1)
      {      
	 int ncam = (int)cam;          
	 double Eccd = ev->E(ncam,ntrack);  
         //  double x = ev->x(ncam,ntrack);  
         //  double y = ev->y(ncam,ntrack);          
	 double Ea = top->at(trig,0).getPeak()*1000;    
         //Energy Matching:       
	 if (Eccd < param[1]*Ea + param[2])pass=false;      
	 if (Eccd > param[3]*Ea + param[4])pass=false; 
         //     if (pass)cout <<"Energy match: "<<Eccd<<" "<<Ea<<endl;    
         //Check waveform properties:      
	 if (top->at(ntrack).getWfMax()>param[5]||top->at(ntrack).getWfMin()<param[6]) 
	    pass=false;      
	 if (fabs(top->at(ntrack).getBase())>param[7]) 
	    pass=false;      
	 if (top->at(ntrack).getRMS()>param[8]) 
	    pass=false;   
//   if (pass)cout <<"Pass base & rms cuts"<<endl;  
	 if (top->at(ntrack,0).getPeak()<param[12]) 
	    pass=false;     
	 if (top->at(ntrack,0).getPeak()<bot->at(ntrack,0).getPeak()) 
	    pass=false;    
//   if (pass)cout <<"Is highest "<<endl;      
	 if (top->at(ntrack,0).getRise10()-top->at(ntrack,0).getRise90() > param[9]) 
	    pass=false;      
	 if (top->at(ntrack,0).getRise10()-top->at(ntrack,0).getRise90() < param[10]) 
	    pass=false;      
	 if (fabs(top->at(ntrack,0).getFall0())<0.5) 
	    pass = false;// && top->at(ntrack,0).getFall0()<param[11]) pass=false;   
        //   if (pass)cout <<"Pass rise and fall"<<endl;      
      }
      else
      {//bottom         
	 int ncam = (int)cam;      
	 double Eccd = ev->E(ncam,ntrack); 
         //   double x = ev->x(ncam,ntrack); 
         //   double y = ev->y(ncam,ntrack);      
	 double Ea = bot->at(trig,0).getPeak()*1000;    
        //Energy Matching:     //>(40,700), (104,2753)    //<(15,1340), (100,4236)      
	 if (Eccd < param[13]*Ea + param[14])pass=false;      
	 if (Eccd > param[15]*Ea + param[16])pass=false;    
         //  if (pass)cout <<"Energy match: "<<Eccd<<" "<<Ea<<endl;    
         //Check waveform properties:      
	 if (bot->at(ntrack).getWfMax()>param[17]||bot->at(ntrack).getWfMin()<param[18]) 
	    pass=false;    //  if (pass) cout <<"Max/Min"<<endl;      
	 if (fabs(bot->at(ntrack).getBase())>param[19]) 
	    pass=false;    //  if (pass) cout <<"Base"<<endl;      
	 if (bot->at(ntrack).getRMS()>param[20]) 
	    pass=false;    //  if (pass)cout <<"Pass base & rms"<<endl;      
	 if (bot->at(ntrack,0).getPeak()<param[23]) 
	    pass=false;      
	 if (bot->at(ntrack,0).getPeak()<bot->at(ntrack,0).getPeak()) 
	    pass=false;    //   if (pass)cout <<"Is highest"<<endl;      
	 if (bot->at(ntrack,0).getRise10()-bot->at(ntrack,0).getRise90() > param[21]) 
	    pass=false;      
	 if (bot->at(ntrack,0).getRise10()-bot->at(ntrack,0).getRise90() < param[22]) 
	    pass=false;      
	 if (fabs(bot->at(ntrack,0).getFall0())<0.5) pass=false; 
         //     if (pass)cout <<"Pass rise and fall"<<endl;    
      }//camera    
      if (pass) return true;  
   }//triggers  
   return false;
}

bool checkCentroid(DmtpcSkimEvent* ev, int c, int t)
{
   double x = ev->x(c,t);
   double y = ev->y(c,t);

   int bin = ((TH2F*)ev->image(c))->FindBin(x,y);

   return ev->cluster(c)->isInCluster(t,bin);
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

     //cout << "Starting collinear merge!!!" << std::endl; 
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

       //cout << "m,b,rxy: " << m[clust] << " , " << b[clust] << " , " << rxy_each[clust] << std::endl; 
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

        //cout << "BOTH HAVE DIRECTION, RESIDUALS: " << resid2 << " , " << resid1 << "MEAN: " << resid_mean <<  std::endl; 

        return resid_mean < max_join_residual; 
    }


    if (fabs(rxy_each[0]) > min_each_rxy) 
    {
       double resid = fabs (ctr_y[1] - m[0] * ctr_x[1] - b[0]) / TMath::Sqrt(1 + m[0]*m[0]); 
       //cout << "FIRST HAS DIRECTION, RESIDUAL: " << resid << std::endl; 
       return resid < max_join_residual; 
    } 

    if (fabs(rxy_each[1]) > min_each_rxy) 
    {
       double resid = fabs (ctr_y[0] - m[1] * ctr_x[0] - b[1]) / TMath::Sqrt(1 + m[1]*m[1]); 
       //cout << "SECOND HAS DIRECTION, RESIDUAL: " << resid << std::endl; 
       return resid < max_join_residual; 
    } 
    
    //cout << "NEITHER HAS DIRECTION, global_rxy: " << global_rxy << std::endl; 
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
       // //cout << "LOOKING AT BIN " << b << std::endl; 
        int bx = bin % nbinsx; 
        int by = ((bin -bx) / nbinsx) % nbinsy; 
        double x = image->GetXaxis()->GetBinCenter(bx);
        double y = image->GetYaxis()->GetBinCenter(by);
       // //cout << "x,y: " << x << "," << y << std::endl; 
        int s;
        double dist = map->distanceToNearestSpacer(x,y,s); 
         positions[i][s].first+=x;  
         positions[i][s].second+=y;  
         npositions[i][s]++; 
       // //cout << "Nearest Spacer: " << s << " at a distance of " << dist << std::endl; 
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
     //cout << " CHECKING " << i << " AGAINST " << t << std::endl; 
     double diag = TMath::Sqrt(2)*image->GetBinWidth(1);
     //Check for the following condition
     // i is close to a spacer
     // j is close to the same spacer
     // i and j are on opposite sides of the spacer
     // i and j are close together
     
     double mind = MaxCamImageTools::minDist(allcluster[i],allcluster[t],const_cast<TH2F*>(image)) ;
     
     //cout << "MIN DIST: " << mind << endl;
     
     if (mind == 0 ) 
     {
        continue; 
     }
     
     bool merged = false;
     bool merges = false;
     
     if (mind < minDistance) 
     {
	//cout << "DISTANCE BETWEEN " << i << " AND " << t << " is, " << mind << ", BELOW THRESH" << std::endl; 
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
	      
	      //cout << "MINIMUM DISTANCE FROM CLUSTER " << i << " TO SPACER " << s << " IS " << dists[i][s] << std::endl;
	      //cout << "MINIMUM DISTANCE FROM CLUSTER " << t << " TO SPACER " << s << " IS " << dists[t][s] << std::endl;
	      //cout << "MERGING ACROSS SPACER" << endl;
	      
	      merges = true;  
	      break; //We've already considered merging these two clusters, so no point in continuing 
	   }
        }
     }
     
     if ((merged && merges) || ((merged || merges) && (mind <= 2 * TMath::Power(image->GetBinWidth(0),2)  || should_collinear_merge(minrxy_global, minrxy_cluster, maxjoin_residual, ls_weight, allcluster[i], allcluster[t],image))))
     {
        //cout << "MERGING " << i << " WITH " << t << std::endl; 
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
   skim->SetBranchStatus("*trigger*",0);
   skim->SetBranchStatus("*clusters*",0);

   for(int i=0; i<int(runnumbers.size()); i++)
   {
      TString filename = "pass/dmtpc_10L_";
//      TString filename = "pass/dmtpc_mc_";
      filename += runnumbers[i];
//      filename += "_cutoffpass.root";
//      filename += "_mixedpass.root";
      filename += "pass.root";
      int valid = cuts->AddFile(filename, -1);
      cout << runnumbers[i] << "," << valid << "\t";
      if(valid==1)
      {
	 TString skimfile = "/net/zwicky/dmtpc/analysis/10L/skim/dmtpc_10L_";
//	 TString skimfile = "skim/dmtpc_mc_";
//	 TString skimfile = "skim/dmtpc_10L_";
	 skimfile+=runnumbers[i];
//	 skimfile+="_cutoffskim.root";
//	 skimfile+="_mixedskim.root";
	 skimfile+="skim.root";
	 skim->AddFile(skimfile,-1);
      }
   }
   cout << endl;
   
   cuts->AddFriend(skim);

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
   double _calibratedrange[ncam][15];
   double _wormclass[ncam][15];
   DmtpcSkimEvent* ev=0;
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
   cuts->SetBranchAddress("_calibratedrange",_calibratedrange);
   cuts->SetBranchAddress("_wormclass",_wormclass);
   cuts->SetBranchAddress("event",&ev);
   
   int nevents = cuts->GetEntries();
   int nospark[ncam];
   int ntracks[ncam];
   int nonrbi[ncam];
   int nonworm[ncam];
   int nonedge[ncam];
   int inrange[ncam];
   int notcutoff[ncam];
   int notremerge[ncam];
   int notedgespark[ncam];
   int notspacer[ncam];
   int hascharge[ncam];
   int isalpha[ncam];
   int momentratio[ncam];
   double time;
   TH2F* hAll[ncam];
   TH2F* hNonWorm[ncam];
   TH2F* hRegion[ncam];
   TH2F* hNuclear[ncam];
   TH2F* hCutoff[ncam];
   TH2F* hMoment[ncam];
   TH2F* hAlpha[ncam];
   TH2F* hRvEC[ncam];
   TH2F* hPhivE[ncam];
   TH1F* hMax[ncam];
   double rangecal[ncam];
   double Ecal[ncam];
   double wormcutoff[ncam];
   cuts->GetEntry(0);
   map<TString,int> cameraMap;
   skim->GetEntry(0);
   TFile* initfile = skim->GetCurrentFile();
   TObjArray* gainMaps = (TObjArray*)initfile->Get("gainMaps");

   //cout << ((DmtpcGainMap*)gainMaps->At(0))->getSpacerSlope(0) << endl;
   
   for(int i=0; i<ncam; i++)
   {
      cameraMap[ev->cameraSerialNumber(i)]=i;
      rangecal[i]=cfg->getRangeCal(ev->cameraSerialNumber(i));
      Ecal[i]=cfg->getEnergyCal(ev->cameraSerialNumber(i));
      wormcutoff[i]=cfg->getWormCutVal(ev->cameraSerialNumber(i));
      
      nospark[i]=0;ntracks[i]=0;nonrbi[i]=0;nonworm[i]=0;
      nonedge[i]=0;inrange[i]=0;hascharge[i]=0; isalpha[i]=0;
      momentratio[i]=0; notcutoff[i]=0; notremerge[i]=0; notedgespark[i]=0;
      notspacer[i]=0;

      TString recoilname = "hAll_";
      recoilname+=i;
      hAll[i]=new TH2F(recoilname,recoilname,100,0,300,50,0,20);
      hAll[i]->GetXaxis()->SetNdivisions(508);
      hAll[i]->GetXaxis()->SetTitle("Recoil Energy (keV)");
      hAll[i]->GetYaxis()->SetTitle("Range (mm)");
      hAll[i]->SetTitle("");
      

      recoilname = "hNonWorm_";
      recoilname+=i;
      hNonWorm[i]=new TH2F(recoilname,recoilname,100,0,300,50,0,20);
      hNonWorm[i]->GetXaxis()->SetNdivisions(508);
      hNonWorm[i]->GetXaxis()->SetTitle("Recoil Energy (keV)");
      hNonWorm[i]->GetYaxis()->SetTitle("Range (mm)");
      hNonWorm[i]->SetTitle("");
      

      recoilname = "hRegion_";
      recoilname+=i;
      hRegion[i]=new TH2F(recoilname,recoilname,100,0,300,50,0,20);
      hRegion[i]->GetXaxis()->SetNdivisions(508);
      hRegion[i]->GetXaxis()->SetTitle("Recoil Energy (keV)");
      hRegion[i]->GetYaxis()->SetTitle("Range (mm)");
      hRegion[i]->SetTitle("");
      

      recoilname = "hNuclear_";
      recoilname+=i;
      hNuclear[i]=new TH2F(recoilname,recoilname,100,0,300,50,0,20);
      hNuclear[i]->GetXaxis()->SetNdivisions(508);
      hNuclear[i]->GetXaxis()->SetTitle("Recoil Energy (keV)");
      hNuclear[i]->GetYaxis()->SetTitle("Range (mm)");
      hNuclear[i]->SetTitle("");
      
      recoilname = "hCutoff_";
      recoilname+=i;
      hCutoff[i]=new TH2F(recoilname,recoilname,100,0,300,50,0,20);
      hCutoff[i]->GetXaxis()->SetNdivisions(508);
      hCutoff[i]->GetXaxis()->SetTitle("Recoil Energy (keV)");
      hCutoff[i]->GetYaxis()->SetTitle("Range (mm)");
      hCutoff[i]->SetTitle("");
      

      recoilname = "hMoment_";
      recoilname+=i;
      hMoment[i]=new TH2F(recoilname,recoilname,100,0,300,50,0,20);
      hMoment[i]->GetXaxis()->SetNdivisions(508);
      hMoment[i]->GetXaxis()->SetTitle("Recoil Energy (keV)");
      hMoment[i]->GetYaxis()->SetTitle("Range (mm)");
      hMoment[i]->SetTitle("");
      
      recoilname = "hRvE_";
      recoilname+=i;
      hAlpha[i]=new TH2F(recoilname,recoilname,500,0,100000,1024,0,1024);
      hAlpha[i]->GetXaxis()->SetNdivisions(508);
      hAlpha[i]->GetXaxis()->SetTitle("Energy (adu)");
      hAlpha[i]->GetYaxis()->SetTitle("Range (px)");
      hAlpha[i]->SetTitle("");

      recoilname = "hRvEC_";
      recoilname+=i;
      hRvEC[i]=new TH2F(recoilname,recoilname,50,0,500,20,0,20);
      hRvEC[i]->GetXaxis()->SetNdivisions(508);
      hRvEC[i]->GetXaxis()->SetTitle("Energy (keV)");
      hRvEC[i]->GetYaxis()->SetTitle("Range (mm)");
      hRvEC[i]->SetTitle("");

      recoilname = "hPhivE_";
      recoilname+=i;
      hPhivE[i]=new TH2F(recoilname,recoilname,500,0,100000,100,-3.14159,3.14159);
      hPhivE[i]->GetXaxis()->SetNdivisions(508);
      hPhivE[i]->GetXaxis()->SetTitle("Energy (adu)");
      hPhivE[i]->GetYaxis()->SetTitle("Phi");
      hPhivE[i]->SetTitle("");
      

      recoilname = "hMax_";
      recoilname+=i;
      hMax[i]=new TH1F(recoilname,recoilname,200,0,100);
      hMax[i]->GetXaxis()->SetTitle("Trigger Max");
      
      
   }

   //cout << "tag b" << endl;

   TString playname = flistfile;
   playname.ReplaceAll(".txt","narrowsearch.play");
   ofstream play;
   play.open(playname);

   double chargepar[] = 
      {  0,//number of top camera  
	 61.2,-2400,//top charge-light matching lines  
	 64.4,143,  0.097,-0.097,//waveform max/min  
	 0.003,//baseline  
	 0.0012,//rms  
	 5e-6,1e-6,//rise time limits  
	 15e-6,//fall time limit, not used right now  
	 0.008,//min peak height  
	 32.1, -584,//bottom charge-light matching  
	 34.1, 830,  0.097,-0.097,//wf min/max  
	 0.002,0.002,//baseline, rms limits  
	 5e-6,1.5e-6,//rise time limits  
	 0.008 //min peak height  
      };


   for(int i=0; i<cuts->GetEntries(); i++)
   {

      if(i%1000==999) cout << "Event: " << i+1 << endl;
      cuts->GetEntry(i); 
      skim->GetEntry(i);
      
      for(int c=0; c<ncam; c++)
      {
	 if(spark->passes(c,0))
	    nospark[cameraMap[ev->cameraSerialNumber(c)]]++;

	 TString camSerNum = ev->cameraSerialNumber(c);
	 
	 for(int t=0; t<15; t++)
	 {
	    if(nonZero->passes(c,t) && spark->passes(c,t))
	    {
	       ntracks[cameraMap[camSerNum]]++;
	       hAll[cameraMap[camSerNum]]->Fill(_Erecoil[c][t],_calibratedrange[c][t]);
	    }
	    if(nonZero->passes(c,t) && spark->passes(c,t) && rbi->passes(c,t))
	       nonrbi[cameraMap[camSerNum]]++;
	    if(nonZero->passes(c,t) && spark->passes(c,t) && rbi->passes(c,t) &&
	       _wormclass[c][t] > wormcutoff[cameraMap[camSerNum]] &&  
	       ev->maxpixel(c,t)<400 && ev->cluster_rms(c,t)>10.5 
	       && ev->maxpixel(c,t)/ev->E(c,t)<0.225)
	    {
	       nonworm[cameraMap[camSerNum]]++;
	       hNonWorm[cameraMap[camSerNum]]->Fill(_Erecoil[c][t],
						    _calibratedrange[c][t]);
	       hMax[cameraMap[camSerNum]]->Fill(ev->transverse_moment(c,2,t));
				
	    }
	    if(nonZero->passes(c,t) && spark->passes(c,t) && rbi->passes(c,t) &&
	       _wormclass[c][t] > wormcutoff[cameraMap[camSerNum]]  
	       &&  edge->passes(c,t)&& region->passes(c,t)
	       &&  ev->maxpixel(c,t)<400 && ev->cluster_rms(c,t)>10.5  
	       && ev->maxpixel(c,t)/ev->E(c,t)<0.225)
	    {
	       nonedge[cameraMap[camSerNum]]++;
	       hRegion[cameraMap[camSerNum]]->Fill(_Erecoil[c][t],
						    _calibratedrange[c][t]);
	    }
	    if(nonZero->passes(c,t) && spark->passes(c,t) && rbi->passes(c,t) &&
	       _wormclass[c][t] > wormcutoff[cameraMap[camSerNum]]  
	       && edge->passes(c,t) 
	       && region->passes(c,t)  &&  ev->maxpixel(c,t)<400 && ev->cluster_rms(c,t)>10.5
	       && ev->maxpixel(c,t)/ev->E(c,t)<0.225
	       && ev->transverse_moment(c,2,t)<20 && 
	       !(ev->phi(c,t)>TMath::Pi()/2.0-0.1 && ev->phi(c,t)< TMath::Pi()/2.0+0.1) &&
	       !(ev->phi(c,t)>-TMath::Pi()/2.0-0.1 && ev->phi(c,t)< -TMath::Pi()/2.0+0.1) &&
	       ev->range(c,t)>25)
	    {
	       inrange[cameraMap[camSerNum]]++;

	       hNuclear[cameraMap[camSerNum]]->Fill(_Erecoil[c][t],
						    _calibratedrange[c][t]);
	       cuts->SetBranchStatus("*clusters*",1); 
	       cuts->SetBranchStatus("*trigger*",1);
	       ev->Delete(); 
	       ev = 0; 
	       cuts->GetEntry(i); 

	       
	       if((!cutoffCheck(ev,c,t)||1) && checkCentroid(ev,c,t))
	       {
		  notcutoff[cameraMap[camSerNum]]++;

		  vector< vector < int> > px;
		  px.clear();
		  px = ev->cluster(c)->getClusters();

		  if(!MergeClusters(ev->image(c),px,t,4000,0.65,0.8,
				    10,2,(DmtpcGainMap*)gainMaps->At(cameraMap[camSerNum])))
		  {
		     notremerge[cameraMap[camSerNum]]++;
		     if(!checkEdgeSpark(ev,c,t))
		     {
			notedgespark[cameraMap[camSerNum]]++;
			

			if(distToSpacer(ev,c,t,(DmtpcGainMap*)gainMaps->At(cameraMap[camSerNum]))>=0)
			{
			   notspacer[cameraMap[camSerNum]]++;

			   //cout << "10L " << ev->runNumber() << " " << ev->eventNumber() 
//				<< " "  << c << " " << t << endl;
			   play << "10L " << ev->runNumber() << " " << ev->eventNumber() 
				<< " "  << c << " " << t << endl;			   
			   if((camSerNum=="081264" && c==0) || (camSerNum=="100439" && c==1) )
			      chargepar[0]=0;
			   else
			      chargepar[0]=1;
			   
			   if(/*chargeCuts10L(ev,c,t,chargepar)*/ 1)
			   {
			      hCutoff[cameraMap[camSerNum]]->Fill(_Erecoil[c][t],
								  _calibratedrange[c][t]);
			      hascharge[cameraMap[camSerNum]]++;

			   }
			   
			}
		     }
		     
		  }
		  
	       }

	    }
	 }
	    
	      
      }
      cuts->SetBranchStatus("*clusters*",0);
      cuts->SetBranchStatus("*trigger*",0);
      ev->Delete();
      ev = 0;

   }

   cout << "Total: \t\t" << nevents << "\t\t\t" << nevents << "\n"
	<< "Non-spark: \t" << nospark[0] << "\t" << nospark[0]/float(nevents)
	<< "\t" << nospark[1] << "\t" << nospark[1]/float(nevents) << "\n"
	<< "-----------------------------------------------------------------------\n" ;
   nospark[0]*=1.0;
   nospark[1]*=1.0;
   cout << "Tracks: \t" << ntracks[0] << "\t" << ntracks[0]/float(nospark[0]) << "\t"
	<< ntracks[1] << "\t" << ntracks[1]/float(nospark[1]) << "\n"
	<< "Non-RBI: \t" << nonrbi[0] << "\t" << nonrbi[0]/float(nospark[0]) << "\t" 
	<< nonrbi[1] << "\t" << nonrbi[1]/float(nospark[1]) << "\n"
	<< "Non-Worm: \t" << nonworm[0] << "\t" << nonworm[0]/float(nospark[0]) 
	<< "\t" << nonworm[1] << "\t" << nonworm[1]/float(nospark[1]) << "\n"
	<< "Edge/Region: \t" << nonedge[0] << "\t" << nonedge[0]/float(nospark[0]) 
	<< "\t" << nonedge[1] << "\t" << nonedge[1]/float(nospark[1]) << "\n"
	<< "Range: \t\t" << inrange[0] << "\t" << inrange[0]/float(nospark[0]) 
	<< "\t" << inrange[1] << "\t" << inrange[1]/float(nospark[1]) << "\n"
	<< "Not Cutoff: \t\t" << notcutoff[0] << "\t" << notcutoff[0]/float(nospark[0]) 
	<< "\t" << notcutoff[1] << "\t" << notcutoff[1]/float(nospark[1]) << "\n"
	<< "Not Remerge: \t\t" << notremerge[0] << "\t" << notremerge[0]/float(nospark[0]) 
	<< "\t" << notremerge[1] << "\t" << notremerge[1]/float(nospark[1]) << "\n"
	<< "Not Edge Spark: \t\t" << notedgespark[0] << "\t" << notedgespark[0]/float(nospark[0]) 
	<< "\t" << notedgespark[1] << "\t" << notedgespark[1]/float(nospark[1]) << "\n"
	<< "Not Near Spacer: \t\t" << notspacer[0] << "\t" << notspacer[0]/float(nospark[0]) 
	<< "\t" << notspacer[1] << "\t" << notspacer[1]/float(nospark[1]) << "\n"
	<< "Charge: \t" << hascharge[0] << "\t" << hascharge[0]/float(nospark[0]) 
	<< "\t" << hascharge[1] << "\t" << hascharge[1]/float(nospark[1]) << "\n";
   
//    gSystem->Setenv("MCTABLES","/net/zwicky/dmtpc/akaboth/projects/DarkMatter/MaxCam/tables");

//    MaxCamSRIM* alpha = new MaxCamSRIM("SRIM_He_in_CF4_100Torr");
//    alpha->setPressure(60);
//    TGraph* gAlpha = alpha->getRangeVsEnergy();
//    gAlpha->SetLineWidth(3);
//    MaxCamSRIM* prot = new MaxCamSRIM("SRIM_H_in_CF4_100Torr");
//    prot->setPressure(60);
//    TGraph* gProt = prot->getRangeVsEnergy();
//    gProt->SetLineWidth(3);

//    TGraph* gAlphaRec = (TGraph*)gAlpha->Clone("gAlphaRec");
//    TGraph* gProtRec = (TGraph*)gProt->Clone("gProtRec");
//    for(int i=0; i<gAlphaRec->GetN(); i++)
//    {
//       double x,y;
//       gAlphaRec->GetPoint(i,x,y);
//       gAlphaRec->SetPoint(i,RecoilEnergy::visibleToRecoil(x),y);
//    }
//    for(int i=0; i<gProtRec->GetN(); i++)
//    {
//       double x,y;
//       gProtRec->GetPoint(i,x,y);
//       gProtRec->SetPoint(i,RecoilEnergy::visibleToRecoil(x),y);
//    }



//    MaxCamSRIM* nuclear = new MaxCamSRIM("SRIM_F_in_CF4_100Torr");
//    nuclear->setPressure(60);
//    TGraph* gNuclear = nuclear->getRangeVsEnergy();
//    gNuclear->SetLineWidth(3);

//    TCanvas* c0 = new TCanvas("c0","c0",0,0,1600,1000);
//    c0->Divide(5,2);

//    c0->cd(1);
//    hAll[0]->Draw("colz");
//    gNuclear->Draw("SAME");
//    gAlphaRec->Draw("SAME");
//    gProtRec->Draw("SAME");
//    c0->cd(2);
//    hNonWorm[0]->Draw("colz");
//    gNuclear->Draw("SAME");
//    gAlphaRec->Draw("SAME");
//    gProtRec->Draw("SAME");
//    c0->cd(3);
//    hRegion[0]->Draw("colz");
//    gNuclear->Draw("SAME");
//    gAlphaRec->Draw("SAME");
//    gProtRec->Draw("SAME");
//    c0->cd(4);
//    hNuclear[0]->Draw("colz");
//    gNuclear->Draw("SAME");
//    gAlphaRec->Draw("SAME");
//    gProtRec->Draw("SAME");
//    c0->cd(5);
//    hCutoff[0]->Draw("colz");
//    gNuclear->Draw("SAME");
//    gAlphaRec->Draw("SAME");
//    gProtRec->Draw("SAME");
//    c0->cd(6);
//    hAll[1]->Draw("colz");
//    gNuclear->Draw("SAME");
//    gAlphaRec->Draw("SAME");
//    gProtRec->Draw("SAME");
//    c0->cd(7);
//    hNonWorm[1]->Draw("colz");
//    gAlphaRec->Draw("SAME");
//    gNuclear->Draw("SAME");
//    gProtRec->Draw("SAME");
//    c0->cd(8);
//    hRegion[1]->Draw("colz");
//    gNuclear->Draw("SAME");
//    gAlphaRec->Draw("SAME");
//    gProtRec->Draw("SAME");
//    c0->cd(9);
//    hNuclear[1]->Draw("colz");
//    gNuclear->Draw("SAME");
//    gAlphaRec->Draw("SAME");
//    gProtRec->Draw("SAME");
//    c0->cd(10);
//    hCutoff[1]->Draw("colz");
//    gNuclear->Draw("SAME");
//    gAlphaRec->Draw("SAME");
//    gProtRec->Draw("SAME");


//    c0->Update();

//    TCanvas* c1 = new TCanvas("c1","c1",0,0,1600,800);
//    c1->Divide(2,1);

//    c1->cd(1);
// //   hMax[0]->Draw();
//    hRvEC[0]->Draw("colz");
//    gAlphaRec->Draw("SAME");
//    gNuclear->Draw("SAME");
//    // hAlpha[0]->ProjectionY()->Draw();
// //   cout << hAlpha[0]->Integral(0,200,130,200)/hAlpha[0]->Integral(0,200,0,200) << endl;

//    c1->cd(2);
// //   hMax[1]->Draw();
//    hRvEC[1]->Draw("colz");
//    gAlphaRec->Draw("SAME");
//    gNuclear->Draw("SAME");
// //   cout << hAlpha[1]->Integral(0,200,137,200)/hAlpha[1]->Integral(0,200,0,200) << endl;
   
//    c1->Update();

// //    TString outrootname = flistfile;
// //    outrootname.ReplaceAll(".txt",".root");
// //    TFile* outroot = new TFile(outrootname,"RECREATE");

// //    outroot->cd();
// //    hAlpha[0]->Write();
// //    hAlpha[1]->Write();
// //    hRvEC[0]->Write();
// //    hRvEC[1]->Write();
// //    hPhivE[0]->Write();
// //    hPhivE[1]->Write();

// //   getchar();
  
//    TString cname = flistfile;
//    cname.ReplaceAll(".txt",".eps");

//    c0->SaveAs(cname);
    

   return 0;

}
