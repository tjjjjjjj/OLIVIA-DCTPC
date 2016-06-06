#include "AnalysisCut.hh"
#include "MultiVariate.hh"
#include <istream>
#include <fstream>
#include "TChain.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "../../../MaxCam/DmtpcSkimPlaylist.hh"
#include "../../../MaxCam/DmtpcRootTools.hh"
#include "TCutG.h"

/** Function Pointer Definitions Here **/

//Range cut
//Parameter 0: minimum length
//Parameter 1: maximum length 
//Parameter n; n>1: lengthcal for camera n-2
bool pass_range_cut( DmtpcSkimEvent * ev, int c , int t,  double * params)
{
  double range = ev->range(c,t) * params[c+2];   
  return (range > params[0] && range < params[1]); 
}; 

//Energy cut, no gainmap
//Parameter 0: minimum E
//Parameter 1: maximum E 
//Parameter n; n>1: Ecal for camera n-2
bool pass_E_cut( DmtpcSkimEvent * ev, int c, int t,  double * params)
{
  double E = ev->E(c,t) / params[c+2];   
  return (E > params[0] && E < params[1]); 
}

//spark cut
// No parameters needed
bool pass_spark_cut( DmtpcSkimEvent * ev, int c, int t __attribute__((unused)),  double * params __attribute__((unused)))
{
  return !ev->spark(c); 
}

//edge cut
//No parameters needed
bool pass_edge_cut( DmtpcSkimEvent * ev, int c, int t,  double * params __attribute__((unused)))
{
  return !ev->edge(c,t); 
}

//region cut
//Parameter 0 Minimum x/y
//Parameter y Maximum x/y

bool pass_region_cut( DmtpcSkimEvent * ev, int c, int t,  double * params)
{
  double x = ev->x(c,t); 
  double y = ev->y(c,t); 

  return (x>params[0] && x < params[1] && y > params[0] && y < params[1]); 
}

//rbi cut 
//Parameter 0 nposthresh
//Parameter 1 sparkrefbins 
bool pass_rbi_cut( DmtpcSkimEvent * ev, int c, int t,  double * params)
{
 return !ev->isRBI(c,t,(int)params[0],(int)params[1]);  
}

bool pass_non_zero_cut( DmtpcSkimEvent *ev, int c, int t, double * params)
{
  return ( ev->E(c,t) > params[0] && 
	   ev->range(c,t) > params[0] && 
	   ev->cluster_rms(c,t) > params[0] && 
	   ev->neighbors(c,t) >= params[0] &&  
	   ev->maxpixel(c,t) > params[0] && 
	   ev->cluster_mean(c,t) > params[0] && 
	   ev->npixel(c,t) > params[0] ); 
}

int main(int nargs, char ** args)
{

  TApplication app("CBA",0,0); 

  //Set up cuts

  vector<AnalysisCut*> cuts; 

  //
  // Make sure output of cleanSkim is non-NULL.  
  //    
  AnalysisCut *nonzerocut =  new AnalysisCut("nonZero",&pass_non_zero_cut,1);
  nonzerocut->setParameter(0,0);
  cuts.push_back(nonzerocut); 

  //Set up i/o
  
  // need to loop over cameras at some point...
  int cam = 0; 
  TString cutfilename = "cba_cuts.root"; 
  if (nargs>2) cam = atoi(args[2]); 
  if (nargs>3) cutfilename = TString(args[3]); 
  
  TChain ch("skim"); 
  
  string temp; 
  ifstream in(args[1]); 

  while (getline(in,temp))
  {
    ch.Add(temp.c_str()); 
  }

  DmtpcSkimEvent * ev = 0; 
  ch.SetBranchAddress("event",&ev); 

  //Set up hists 
  
  TH1F * edist = new TH1F("Edist","Edist",100,0,5000); 
  TH1F * phidist = new TH1F("Phidist","Phidist",10,-3.14159265359,3.14159265359); 
  TH2F * rve = new TH2F("RVE","RVE",100,0,10000,13,0,100);

  DmtpcSkimPlaylist out; 
  //Do the loopie
  
  TFile* cutfile = new TFile(cutfilename,"RECREATE");

  // tree to store analysis-level quantities of tracks that pass all cuts
  TTree* bdat= new TTree("bdat","a really basic analysis tree for DMTPC datasets");
  Double_t theta,phi,E,range,maxpixel,clusterrms,clustermean,x,y;
  Int_t xmin,ymin,xmax,ymax,npixels, neighbors;
  Int_t spark, nburnin, camnum,runnum,evtnum,trknum, isrbi23;

  bdat->Branch("theta",&theta,"theta/D");
  bdat->Branch("phi",&phi,"phi/D");
  bdat->Branch("E",&E,"E/D");
  bdat->Branch("range",&range,"range/D");
  bdat->Branch("maxpixel",&maxpixel,"maxpixel/D");
  bdat->Branch("clusterrms",&clusterrms,"clusterrms/D");
  bdat->Branch("clustermean",&clustermean,"clustermean/D");
  bdat->Branch("x",&x,"x/D");
  bdat->Branch("y",&y,"y/D");
  bdat->Branch("xmin",&xmin,"xmin/I");
  bdat->Branch("ymin",&ymin,"ymin/I");
  bdat->Branch("xmax",&xmax,"xmax/I");
  bdat->Branch("ymax",&ymax,"ymax/I");
  bdat->Branch("npixels",&npixels,"npixels/I");
  bdat->Branch("neighbors",&neighbors,"neighbors/I");
  bdat->Branch("spark",&spark,"spark/I");
  bdat->Branch("isrbi23",&isrbi23,"isrbi23/I");
  bdat->Branch("nburnin",&nburnin,"nburnin/I");
  bdat->Branch("runnum",&runnum,"runnum/I");
  bdat->Branch("evtnum",&evtnum,"evtnum/I");
  bdat->Branch("camnum",&camnum,"camnum/I");
  bdat->Branch("trknum",&trknum,"trknum/I");

  for (int i = 0; i < ch.GetEntries(); i++) {
    ch.GetEntry(i); 
    
    for (size_t cut = 0; cut < cuts.size(); cut++) {
      cuts[cut]->evaluateEvent(ev);       
    }
    
    //    cut_tree->Fill(); 
    
    // loop over the cameras ...
    for(int cam = 0; cam<ev->ncamera(); ++cam) {
      for (int t = 0; t < ev->ntracks(cam);t++) {
	bool fails = false; 
	for (size_t cut = 0; cut < cuts.size(); cut++) {
	  if (!cuts[cut]->passes(cam,t)) {
	    fails = true; 
	    break;
	  }
	}
	
	if (!fails) {
	  cout << "PASSING EVENT!!! Run: " << ev->runNumber() << " Event: " << ev->eventNumber() << " Track: "<<t<< " Cam: " << cam << /* " Classifier: " << result.getClassifier(ev,cam,t) <<*/  endl;  
	  
	  edist->Fill(ev->E(cam,t)); 
	  phidist->Fill(ev->phi(cam,t)); 
	  rve->Fill(ev->E(cam,t),ev->range(cam,t)); 
	  
	  // fill the analysis tree
	  theta=ev->theta(cam,t);
	  phi=ev->phi(cam,t);
	  E=ev->E(cam,t);
	  range=ev->range(cam,t);
	  maxpixel=ev->maxpixel(cam,t);
	  clusterrms=ev->cluster_rms(cam,t);
	  neighbors=ev->neighbors(cam,t);
	  clustermean=ev->cluster_mean(cam,t);
	  x=ev->x(cam,t);
	  y=ev->y(cam,t);
	  spark=( ev->spark(cam) ? 1 : 0 );
	  isrbi23=ev->isRBI(cam,t,2,3);
	  nburnin=ev->nburnin(cam,t);
	  npixels=ev->npixel(cam,t);
	  // extract the cluster bounds
	  //	MaxCamClusterImage* thisCluster = (MaxCamClusterImage*)sd.event()->cluster(icam);
	  ((MaxCamClusterImage*)ev->cluster(cam))->clusterBounds(t,&xmin,&xmax,&ymin,&ymax,0);
	  
	  camnum=cam;
	  runnum=ev->runNumber();
	  evtnum=ev->eventNumber();
	  trknum=t;
	  
	  bdat->Fill();
	  
	  out.add("Raytheon",ev->runNumber(),ev->eventNumber(),cam,t); 
	}
      }
    }
  }
  
  out.save("cba_out.play");

  DmtpcRootTools::setColorStandard1();
  TCanvas c1("EDistC","Energy Distribution",400,400); 
  edist->Draw(); 
  TCanvas c2("PhiDistC","Phi Distribution",400,400); 
  phidist->Draw(); 
  TCanvas c3("RvE","Range vs. Energy",400,400); 
  gStyle->SetPalette(1); 
  rve->Draw("colz"); 
  
  //  cut_tree->Write(); 
  bdat->Write();
  //  app.Run(); 
}
