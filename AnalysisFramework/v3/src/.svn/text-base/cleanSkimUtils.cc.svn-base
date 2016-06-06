#include "../../../MaxCam/DmtpcSkimDataset.hh"
#include <iostream>
#include <strings.h>
#include <stdlib.h>
#include "../../../MaxCam/DmtpcSkimEvent.hh"
#include "../../../MaxCam/MaxCamClusterImage.hh"
#include "TCanvas.h"
#include "TApplication.h" 
#include "style.hh"


using namespace std; 


typedef struct 
{
  int minx; 
  int maxx;
  int miny;
  int maxy; 
} cluster_bounds_t; 


static cluster_bounds_t get_bounds(DmtpcSkimEvent * e, int cam, int t, int radius=20)
{
  MaxCamClusterImage * im = e->cluster(cam); 

  int x,y; 
  cluster_bounds_t bounds; 
  bounds.maxx = 0; 
  bounds.maxy = 0; 
  bounds.miny = 1 << 30; 
  bounds.minx = 1 << 30; 
  
  for (int i = 0; i < im->getCluster(t).size(); i++)
  {
    int bin = im->getCluster(t)[i]; 
    im->getXYFromBinNo(bin,&x,&y,true); 
    bounds.minx = x < bounds.minx ? x : bounds.minx;  
    bounds.miny = y < bounds.miny ? y : bounds.miny;  
    bounds.maxx = x > bounds.maxx ? x : bounds.maxx;  
    bounds.maxy = y > bounds.maxy ? y : bounds.maxy;  

  }

  bounds.minx-=radius;
  bounds.miny-=radius;
  bounds.maxx+=radius;
  bounds.maxy+=radius;

  return bounds; 
} 

static cluster_bounds_t all_bounds()
{
   cluster_bounds_t bounds; 
   bounds.maxx = 1024; 
   bounds.maxy = 1024; 
   bounds.minx = 0; 
   bounds.miny = 0; 
   return bounds; 

}

static int draw_cluster(DmtpcSkimEvent * e, int cam, cluster_bounds_t bounds)
{ 
  MaxCamClusterImage * im = e->cluster(cam); 
  im->getImage()->SetAxisRange(bounds.minx,bounds.maxx,"X");
  im->getImage()->SetAxisRange(bounds.miny,bounds.maxy,"Y");
  im->getImage()->DrawCopy("colz"); 
  return 0; 
}

int help()
{
  cout << "cleanSkimUtils cmd [cmd_args ...]" << endl; 
  cout << "where cmd is one of: " << endl; 
  cout << "\tplotcluster skimfile cam eventnum clusternum [radius]" <<endl; 
  cout << "\tplotframe skimfile cam eventnum " <<endl; 
  cout << "\tplotburnin skimfile cam eventnum clusternum [nburnin] [radius]" <<endl; 
  cout << "\tprintevent skimfile eventumstart [eventnumend]" <<endl; 
  cout << "\tprintcluster skimfile cam eventnum clusternum" <<endl; 
  cout << "\tplotE skimfile cam minrange" <<endl; 
  cout << "\thelp" << endl; 
  return 0; 
}


static int loadEventnum(DmtpcSkimDataset * ds, int n)
{
  int i = n;
  if (i < 0) i = 0; 
  if (i > ds->nevents()) i = ds->nevents()-1;  
  ds->getEvent(i); 
  while(ds->event()->eventNumber()!=n)
  {
    int last = ds->event()->eventNumber(); 
    if (last < n)
    {
      ds->getEvent(++i);  
    } 
    else {
      ds->getEvent(--i); 
    }
    int now = ds->event()->eventNumber(); 
    if ((last <n && now > n) || (last > n && now < n)) return 1; 
  }

  return 0; 
}


int plotframe(int argn, char ** args)
{
  if (argn < 3) return help(); 
  
  TApplication * app = new TApplication("plotcluster",0,0); 
  AnalysisStyle::setStyle();
  DmtpcSkimDataset ds; 
  ds.openRootFile(args[0]); 
  int cam = atoi(args[1]); 
  int eventnum = atoi(args[2]); 
  if(loadEventnum(&ds,eventnum)) return 1; 
  
  TCanvas * canvas = new TCanvas("cluster","Cluster",400,400);
  canvas->Show();  
  canvas->cd(); 
  draw_cluster(ds.event(),cam,all_bounds());  
  canvas->Update(); 
  app->Run(); 
  return 0; 
}

int plotcluster(int argn, char ** args)
{
  if (argn < 4) return help(); 
  
  TApplication * app = new TApplication("plotcluster",0,0); 
  AnalysisStyle::setStyle();
  DmtpcSkimDataset ds; 
  ds.openRootFile(args[0]); 
  int cam = atoi(args[1]); 
  int eventnum = atoi(args[2]); 
  if(loadEventnum(&ds,eventnum)) return 1; 
  int clusternum = atoi(args[3]);
  if (clusternum >= ds.event()->ntracks(cam))
  {
    std::cerr << "Event contains no cluster " << clusternum <<std::endl;
    return 1;
  }
  TCanvas * canvas = new TCanvas("cluster","Cluster",400,400);
  canvas->Show();  
  canvas->cd(); 
  if (argn == 5) 
  {
    int radius = atoi(args[4]);
    draw_cluster(ds.event(),cam,get_bounds(ds.event(),cam,clusternum,radius));  
  }
  else
  {
    draw_cluster(ds.event(),cam,get_bounds(ds.event(),cam,clusternum));  
  }
  canvas->Update(); 
  app->Run(); 
  return 0; 
}

int plotburnin(int argn, char ** args)
{
  if (argn < 4) return help(); 
  DmtpcSkimDataset ds; 
  ds.openRootFile(args[0]);
  int cam = atoi(args[1]);
  int eventnum = atoi(args[2]);
  int clusternum = atoi(args[3]);  
  int nburnin = 20; 
  bool userad = false; 
  int radius; 
  
  if (argn > 4) nburnin = atoi(args[4]);
  if (argn > 5) 
  {
    userad= true; 
    radius = atoi(args[5]); 
  }
  
  TApplication * app = new TApplication("plotburnin",0,0); 
  AnalysisStyle::setStyle(); 
  //TODO: Figure out a "real" heuristic for this  
  int nrows = 2;
  if (nburnin % 2 == 0 && nburnin > 5)
  {
    nrows=4; 
  }
  if (nburnin % 4 == 0 && nburnin > 10)
  {
    nrows= 8; 
  }

  TCanvas * c = new TCanvas("c","Burnin Plots"); 
  c->Divide(2*nburnin/nrows,nrows);
  if(loadEventnum(&ds,eventnum)) return 1;  
 
  
  cluster_bounds_t bounds;
  if(userad) bounds = get_bounds(ds.event(),cam,clusternum,radius);
  else bounds = get_bounds(ds.event(),cam,clusternum); 
  
  c->Show(); 
  int n = 1;
  for (int i = eventnum - nburnin; i <= eventnum + nburnin; i++)
  {
    if (i < 0) continue; 
    if (i >= ds.nevents()) break; 
    if (i==eventnum) continue;
    if(loadEventnum(&ds,i)) continue; 
    c->cd(n++); 
    if (userad) draw_cluster(ds.event(),cam,bounds);
    else        draw_cluster(ds.event(),cam,bounds);
    c->Update(); 
  }
  app->Run(); 

  return 0; 
}


int printcluster(int argn, char ** args)
{
 if (argn <3) return help(); 
 DmtpcSkimDataset ds; 
 ds.openRootFile(args[0]); 
 int c = atoi(args[1]); 
 int eventnum = atoi(args[2]);
 int t = atoi(args[3]); 
 if(loadEventnum(&ds,eventnum))
 {
   cout << "Could not load event" << eventnum << endl; 
   return 1;
 }
 DmtpcSkimEvent * e = ds.event(); 
 e->printOut(c,t);  

 return 0; 
} 

int printevent(int argn, char ** args)
{
  if (argn < 2) return help(); 
  DmtpcSkimDataset ds;
  ds.openRootFile(args[0]);
  int eventstart = atoi(args[1]);
  int eventend = argn > 2 ? atoi(args[2]) : eventstart; 
  
  
  for (int i = eventstart; i <= eventend; i++)
  {
    if (i < 0) continue; 
    if (i >= ds.nevents()) break; 
    if(loadEventnum(&ds,i))
    {
      cout << "Could not load event " << i << endl; 
      continue; 
    }
    ds.event()->printOut(); 
  }
}

int plotE(int argn, char ** args)
{
  if (argn < 3) return help(); 
  TApplication * app = new TApplication("plotE",0,0); 
  AnalysisStyle::setStyle(); 
  DmtpcSkimDataset ds; 
  ds.openRootFile(args[0]); 
  int c = atoi(args[1]); 
  double minrange = atof(args[2]); 

  TH1F * E = new TH1F("E","E",50,0,200000); 
  TH1F * phi = new TH1F("phi","phi",40,-3.14159,3.14159); 
  for (int i = 0; i < ds.nevents(); i++)
  {
    ds.getEvent(i); 
    for (int t = 0; t < ds.event()->ntracks(c); t++)
    {
      if (ds.event()->range(c,t) > minrange) 
      {
        E->Fill(ds.event()->E(c,t)); 
        phi->Fill(ds.event()->phi(c,t)); 
      }
    } 
  }

 // E->Draw(); 
  phi->Draw(); 
  app->Run();
}

int main(int argn, char ** args) 
{
  if (argn==1 || strcmp(args[1],"help")==0)
  {
    return help(); 
  }
  
  if (strcmp(args[1],"plotcluster")==0)
  {
    return plotcluster(argn-2, args+2); 
  }

  if (strcmp(args[1],"plotburnin")==0)
  {
    return plotburnin(argn-2, args+2); 
  }

  if (strcmp(args[1],"printcluster")==0)
  {
    return printcluster(argn-2, args+2); 
  }

  if (strcmp(args[1],"plotframe")==0)
  {
    return plotframe(argn-2, args+2); 
  }
  if (strcmp(args[1],"printevent")==0)
  {
    return printevent(argn-2, args+2); 
  }
  if (strcmp(args[1],"plotE")==0)
  {
    return plotE(argn-2, args+2); 
  }

  
}
