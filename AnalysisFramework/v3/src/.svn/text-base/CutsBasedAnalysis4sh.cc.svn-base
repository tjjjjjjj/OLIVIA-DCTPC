#include "AnalysisCut.hh"
#include "MultiVariate.hh"
#include <istream>
#include <fstream>
#include "TChain.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "../../../MaxCam/DmtpcSkimPlaylist.hh"
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

//active_region cut
//No parameters needed
bool pass_active_region_cut( DmtpcSkimEvent * ev, int c, int t,  double * params __attribute__((unused)))
{
  if(c==0){ // camera 0 ...

    // geometry->pixels
    Double_t sideboundarywidth=50;
    Double_t fieldcageinnerradius=955;
    Double_t fieldcageouterradius=fieldcageinnerradius+135;
    Double_t fieldcagexcenter=0;
    Double_t fieldcageycenter=1000;

    Bool_t leftedge=false;
    Bool_t topedge=false;
    Bool_t fieldcageedge=false;
    Bool_t insidefieldcage=false;
    Bool_t outsidefieldcage=false;
    
    // get the pixels corresponding to this track cluster ..
    vector<int> pix = ev->cluster(c)->getCluster(t);

    // loop over the track pixels
    for (int ipix=0; ipix<pix.size(); ipix++) {
      
      Int_t binx,biny;
      ev->cluster(c)->getXYFromBinNo(pix[ipix],&binx,&biny,true);
      //      cout << "binx: " << binx << " biny: " << biny << endl;

      if(binx>0&&biny>0){      
	if(binx<sideboundarywidth) leftedge=true;
	if(biny>1024-sideboundarywidth) topedge=true;
	
	Double_t pixelDistanceFromCenterOfFieldCage=TMath::Sqrt( TMath::Power( ( binx-fieldcagexcenter ) , 2 ) +
	 							 TMath::Power( ( biny-fieldcageycenter ) , 2 ) );
	if(pixelDistanceFromCenterOfFieldCage>fieldcageinnerradius) outsidefieldcage=true;
	if(pixelDistanceFromCenterOfFieldCage<fieldcageinnerradius) insidefieldcage=true;
	if(pixelDistanceFromCenterOfFieldCage<fieldcageinnerradius && pixelDistanceFromCenterOfFieldCage>(fieldcageinnerradius-sideboundarywidth)) fieldcageedge=true;
      }
    }
    
    // clean up (do I need this?  is this a pointer or data in memory that I created?)
    pix.erase(pix.begin(),pix.end());
    
    return !leftedge && !topedge && !outsidefieldcage && !fieldcageedge;
  }
  
  return false;
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

/*
//MVA worm cut 
//Parameter 0 pointer to correctly initialized MultiVariateResult
//Parameter 1 pointer to double cutoff value
bool pass_worm_cut( DmtpcSkimEvent *ev, int c, int t, vector<void *> params)
{
  MultiVariate::MultiVariateResult * r = (MultiVariate::MultiVariateResult*) params[0]; 
  double cutoff = *((double*) params[1]); 
  if (ev->range(c,t) <= 0) return false; 
  //cout << ev->eventNumber() << " " << cutoff << " " << r->getClassifier(ev,c,t) << endl; 
  return r->getClassifier(ev,c,t) > cutoff; 
}
*/

//cut around where we know the alpha source was in the data
//Parameter 0 pointer to correctly initialized TCutG 
bool pass_known_am241_source_position_cut( DmtpcSkimEvent *ev, int c, int t, vector<void *> params)
{
  TCutG * tcg = (TCutG*) params[0]; 

  Bool_t in_known_am241_source_position=false;

  // get the pixels corresponding to this track cluster ..
  vector<int> pix = ev->cluster(c)->getCluster(t);
  // loop over the track pixels
  for (int ipix=0; ipix<pix.size(); ipix++) {
    
    Int_t binx,biny;
    
    ev->cluster(c)->getXYFromBinNo(pix[ipix],&binx,&biny,true);
    //      cout << "binx: " << binx << " biny: " << biny << endl;
    
    if(binx>0&&biny>0){      
      if( tcg->IsInside( binx,biny ) ){
	in_known_am241_source_position=true;
	break;
      }
    }

  }
    
  // clean up (do I need this?  is this a pointer or data in memory that I created?)
  pix.erase(pix.begin(),pix.end());
  
  return !in_known_am241_source_position;
} 

//cut around where we know the nuclear recoils should be in range VS energy
//Parameter 0 pointer to correctly initialized TCutG 
bool pass_nuclear_range_energy_ROI_cut( DmtpcSkimEvent *ev, int c, int t, vector<void *> params)
{
  TCutG * tcg = (TCutG*) params[0]; 
  return tcg->IsInside( ev->E(c,t),ev->range(c,t) );
} 

//simple worm cut
//Parameter 0 maximum maxpixel/range
//Parameter 1 maximum maxpixel
bool pass_simple_worm_cut( DmtpcSkimEvent * ev, int c, int t,  double * params)
{
  return ( ( ev->range(c,t) > 0 ? ev->maxpixel(c,t)/ev->range(c,t)<params[0] : false ) && //maxpixel/range cut
	   ( ev->maxpixel(c,t)<params[1] ) );                                             //maxpixel cut
}

bool pass_non_zero_cut( DmtpcSkimEvent *ev, int c, int t, double * params)
{
  return ( ev->E(c,t)>params[0] && 
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
  // is in the active region of the tpc, and away from the rings/walls?
  //
  AnalysisCut* activeregioncut = new AnalysisCut("activeregion",&pass_active_region_cut,0); 
  cuts.push_back(activeregioncut); 

  //
  // spark cut
  //
  AnalysisCut* sparkcut = new AnalysisCut("spark",&pass_spark_cut,0); 
  cuts.push_back(sparkcut); 
  
  //
  // rbi cut
  //
  // DEFAULT 10L SETUP:WORK TO UNDERSTAND WHAT THIS IS DOING
  AnalysisCut* rbicut = new AnalysisCut("rbi",&pass_rbi_cut,2); 
  rbicut->setParameter(0,2); 
  rbicut->setParameter(1,3); 
  cuts.push_back(rbicut); 
  
  //
  // worm cut
  //
  // NEED THE FISHER DISCRIMINANT-BASED CUT, BUT USE THIS FOR NOW
  AnalysisCut* simplewormcut = new AnalysisCut("simpleworm",&pass_simple_worm_cut,2); 
  simplewormcut->setParameter(0,100); 
  simplewormcut->setParameter(1,3000); 
  cuts.push_back(simplewormcut); 

  /*
    AnalysisCut* wormcut = new AnalysisCut("worm",&pass_worm_cut); 
    MultiVariate::MultiVariateResult result;
    result.setResult("wr4bottom"); 
    double cutoff = 3; 
    vector<void*>worm_parameters; 
    worm_parameters.push_back((void*) &result); 
    worm_parameters.push_back((void*) &cutoff); 
    wormcut->setParameters(worm_parameters); 
    cuts.push_back(wormcut); 
  */ 

  //
  // known americium source position
  //
  TFile alphasourcecutfile("../../4shooter/analysis/cf252/positioncuts.root");
  TCutG *alphasourcecut=(TCutG*)alphasourcecutfile.Get("alphasourcecut");
  cout << "alphasourcecut=" << alphasourcecut << endl;

  AnalysisCut* knownam241sourcepositioncut = new AnalysisCut("knownam241sourcepositioncut",&pass_known_am241_source_position_cut); 
  vector<void*>knownam241sourcepositioncut_parameters;
  knownam241sourcepositioncut_parameters.push_back((void*) alphasourcecut); 
  knownam241sourcepositioncut->setParameters(knownam241sourcepositioncut_parameters); 
  cuts.push_back(knownam241sourcepositioncut); 

  //
  // nuclear recoil range-energy ROI
  //
  TFile nrrvseroifile("../../4shooter/analysis/cf252/NuclearRangeVSECut.root");
  TCutG *nrrvseroicut=(TCutG*)nrrvseroifile.Get("NuclearRangeVSECut");
  
  AnalysisCut* nuclearrangeenergyroicut = new AnalysisCut("nuclearrangeenergyroicut",&pass_nuclear_range_energy_ROI_cut);
  vector<void*>nuclearrangeenergyroicut_parameters;
  nuclearrangeenergyroicut_parameters.push_back((void*) nrrvseroicut);
  nuclearrangeenergyroicut->setParameters(nuclearrangeenergyroicut_parameters); 
  cuts.push_back(nuclearrangeenergyroicut); 

  //
  // time since last spark (can keep track of the spark times, and use them
  //  to actively cut)
  //
  /*
  AnalysisCut* nuclearrangeenergyroicut = new AnalysisCut("nuclearrangeenergyroicut",&pass_nuclear_range_energy_ROI_cut);
  vector<void*>nuclearrangeenergyroicut_parameters;
  nuclearrangeenergyroicut_parameters.push_back((void*) nrrvseroicut);
  nuclearrangeenergyroicut->setParameters(nuclearrangeenergyroicut_parameters); 
  cuts.push_back(nuclearrangeenergyroicut); 
  */
  
  //
  // Make sure output of cleanSkim is non-NULL.  
  //    INVESTIGATE THE RANGE=0 TRACKS; IT SEEMS LIKE THE REASON WHY MOST TRACKS
  //    FAIL THIS CUT IS THAT RANGE=0.
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
  ch.SetBranchStatus("*trigger*",0);
  //  ch.SetBranchStatus("*clusters*",0);


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
  TTree* cut_tree = new TTree("cuts","cuts");

  for (size_t cut = 0; cut < cuts.size(); cut++)
  {
    cut_tree->Branch(cuts[cut]->GetName(),"AnalysisCut",&(cuts[cut]),128000,1);
  }
  
  for (int i = 0; i < ch.GetEntries(); i++)
  {
    ch.GetEntry(i); 

    for (size_t cut = 0; cut < cuts.size(); cut++)
    {

      cuts[cut]->evaluateEvent(ev);       
    }

    cut_tree->Fill(); 

    for (int t = 0; t < ev->ntracks(cam);t++)
    {
      bool fails = false; 
      for (size_t cut = 0; cut < cuts.size(); cut++)
      {
        if (!cuts[cut]->passes(cam,t))
        {
          fails = true; 
          break;
        }
      }

      if (!fails)
      {
	cout << "PASSING EVENT!!! Run: " << ev->runNumber() << " Event: " << ev->eventNumber() << " Track: "<<t<< " Cam: " << cam << /* " Classifier: " << result.getClassifier(ev,cam,t) <<*/  endl;  

        edist->Fill(ev->E(cam,t)); 
        phidist->Fill(ev->phi(cam,t)); 
        rve->Fill(ev->E(cam,t),ev->range(cam,t)); 

        out.add("4sh",ev->runNumber(),ev->eventNumber(),cam,t); 

      }
    }
  }

  out.save("cba_out.play");
  TCanvas c1("EDistC","Energy Distribution",400,400); 
  edist->Draw(); 
  TCanvas c2("PhiDistC","Phi Distribution",400,400); 
  phidist->Draw(); 
  TCanvas c3("RvE","Range vs. Energy",400,400); 
  gStyle->SetPalette(1); 
  rve->Draw("colz"); 
  nrrvseroicut->Draw("lp");
  
  
  cut_tree->Write(); 
  app.Run(); 
}
