// waveformtools includes
//  for analysis of generic waveforms
#include "../../../MaxCam/waveformtools/include/WaveformVector.hh"
#include "../../../MaxCam/waveformtools/include/SkimWaveform.hh"
#include "../../../MaxCam/waveformtools/include/DmtpcPulse.hh"
//  for analysis of voltage sensitive amplitifier waveforms
#include "../../../MaxCam/waveformtools/include/FastWfVector.hh"
#include "../../../MaxCam/waveformtools/include/FastWaveform.hh"
#include "../../../MaxCam/waveformtools/include/FastPulse.hh"
//  for analysis of charge sensitive pre-amp waveforms
#include "../../../MaxCam/waveformtools/include/CspWfVector.hh"
#include "../../../MaxCam/waveformtools/include/CspWaveform.hh"
#include "../../../MaxCam/waveformtools/include/CspPulse.hh"
//  for analysis of charge sensitive pmt waveforms
#include "../../../MaxCam/waveformtools/include/PMTWfVector.hh"
#include "../../../MaxCam/waveformtools/include/PMTWaveform.hh"
#include "../../../MaxCam/waveformtools/include/PMTPulse.hh"
//  general functions for waveform analysis
#include "../../../MaxCam/waveformtools/include/WaveformAnalysis.hh"
#include "../../../MaxCam/waveformtools/include/WaveformTools.hh"

#include "AnalysisCut.hh"
#include "MultiVariate.hh"
#include <istream>
#include <fstream>
#include "TChain.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "../../../MaxCam/DmtpcSkimPlaylist.hh"
#include "../../../MaxCam/Dmtpc4ShooterStitcher.hh"
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
bool pass_active_region_position_cut( DmtpcSkimEvent *ev, int c, int t, vector<void *> params)
{
  // assume this track is fully contained within the active region
  bool track_is_contained_in_active_region=true;
  
  Dmtpc4ShooterStitcher * d4ss = (Dmtpc4ShooterStitcher*) params[0]; 
  double innerRadius_clearance_in_pixels = *((double*) params[1]); 
  
  double innerRadius=d4ss->innerRadius(c);
  double outerRadius=d4ss->outerRadius(c);
  double xCenter=d4ss->xCenter(c);
  double yCenter=d4ss->yCenter(c);
  

  // why is this necessary?!?!?
  if( t < ev->ntracks(c) ){
    
    //    cout << "pass_active_region_position_cut( DmtpcSkimEvent *ev, int c, int t, vector<void *> params): c=" << c
    //	       << " t=" << t 
    //	       << endl;


    // get the pixels corresponding to this track cluster ..
    vector<int> pix = ev->cluster(c)->getCluster(t);

    // loop over the track pixels
    for (int ipix=0; ipix<pix.size(); ipix++) {
      
      Int_t binx,biny;
      ev->cluster(c)->getXYFromBinNo(pix[ipix],&binx,&biny,true);
      //      cout << "binx: " << binx << " biny: " << biny << endl;
      
      // find the distance from this pixel of the track to the center of the anode
      double pixel_distance_from_this_track_pixel_to_anode_center=TMath::Sqrt( TMath::Power( ( binx-xCenter ) , 2 ) +
									       TMath::Power( ( biny-yCenter ) , 2 ) );
      
      // see if this pixel is contained inside the inner-radius and the inner-radius
      // tolerance specified by the user
      if( pixel_distance_from_this_track_pixel_to_anode_center < ( innerRadius-innerRadius_clearance_in_pixels ) ); // do nothing
      else track_is_contained_in_active_region=false;
    }
    
    // clean up (do I need this?  is this a pointer or data in memory that I created?)
    pix.erase(pix.begin(),pix.end());
  }
  
  return track_is_contained_in_active_region;
} 

/*
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
*/

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
  /*
  cout << "pass_non_zero_cut( DmtpcSkimEvent *ev, int c, int t, double * params): c=" << c
       << " t=" << t 
       << endl;
  cout << "\t"
       << ( ev->E(c,t) > params[0] && 
	    ev->range(c,t) > params[0] && 
	    ev->cluster_rms(c,t) > params[0] && 
	    ev->neighbors(c,t) >= params[0] &&  
	    ev->maxpixel(c,t) > params[0] && 
	    ev->cluster_mean(c,t) > params[0] && 
	    ev->npixel(c,t) > params[0] )
       << endl;
  */

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

  /* 
  //
  // is in the active region of the tpc, and away from the rings/walls?
  //
  //  AnalysisCut* activeregioncut = new AnalysisCut("activeregion",&pass_active_region_cut,0); 
  //  cuts.push_back(activeregioncut); 
  
  //
  // spark cut
  //
  //  AnalysisCut* sparkcut = new AnalysisCut("spark",&pass_spark_cut,0); 
  //  cuts.push_back(sparkcut); 
  
  //
  // rbi cut
  //
  // DEFAULT 10L SETUP:WORK TO UNDERSTAND WHAT THIS IS DOING
  //  AnalysisCut* rbicut = new AnalysisCut("rbi",&pass_rbi_cut,2); 
  //  rbicut->setParameter(0,2); 
  //  rbicut->setParameter(1,3); 
  //  cuts.push_back(rbicut); 
  
  //
  // worm cut
  //
  // NEED THE FISHER DISCRIMINANT-BASED CUT, BUT USE THIS FOR NOW
  //  AnalysisCut* simplewormcut = new AnalysisCut("simpleworm",&pass_simple_worm_cut,2); 
  //  simplewormcut->setParameter(0,100); 
  //  simplewormcut->setParameter(1,3000); 
  //  cuts.push_back(simplewormcut); 
  
  AnalysisCut* wormcut = new AnalysisCut("worm",&pass_worm_cut); 
  MultiVariate::MultiVariateResult result;
  result.setResult("wr4bottom"); 
  double cutoff = 3; 
  vector<void*>worm_parameters; 
  worm_parameters.push_back((void*) &result); 
  worm_parameters.push_back((void*) &cutoff); 
  wormcut->setParameters(worm_parameters); 
  cuts.push_back(wormcut); 
  
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
  AnalysisCut* nuclearrangeenergyroicut = new AnalysisCut("nuclearrangeenergyroicut",&pass_nuclear_range_energy_ROI_cut);
  vector<void*>nuclearrangeenergyroicut_parameters;
  nuclearrangeenergyroicut_parameters.push_back((void*) nrrvseroicut);
  nuclearrangeenergyroicut->setParameters(nuclearrangeenergyroicut_parameters); 
  cuts.push_back(nuclearrangeenergyroicut); 
  */
  
  TFile fstitch("stitch.root");
  Dmtpc4ShooterStitcher *stitcher=(Dmtpc4ShooterStitcher*)fstitch.Get("stitch");

  /*
  double innerRadius_clearance_in_pixels=50.;  

  AnalysisCut* activeregionpositioncut = new AnalysisCut("activeregionposition",&pass_active_region_position_cut); 
  vector<void*>activeregionpositioncut_parameters;
  activeregionpositioncut_parameters.push_back((void*) stitcher);  
  activeregionpositioncut_parameters.push_back((void*) &innerRadius_clearance_in_pixels);
  activeregionpositioncut->setParameters(activeregionpositioncut_parameters);
  cuts.push_back(activeregionpositioncut); 
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
  //  ch.SetBranchStatus("*trigger*",0);
  // ch.SetBranchStatus("*clusters*",0);

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
  // I DON'T KNOW WHAT THIS IS FOR, SO I'VE GOT IT OFF UNTIL
  // I FIGURE OUT WHAT IT'S SUPPOSED TO DO ...
  //  TTree* cut_tree = new TTree("cuts","cuts");

  /*
  for (size_t cut = 0; cut < cuts.size(); cut++)
  {
    cut_tree->Branch(cuts[cut]->GetName(),"AnalysisCut",&(cuts[cut]),128000,1);
  }
  */
  
  // tree to store analysis-level quantities of tracks that pass all cuts
  TTree* bdat= new TTree("bdat","a really basic analysis tree for DMTPC datasets");
  Double_t theta,phi,E,range,maxpixel,clusterrms,clustermean,x,y,mindist2innerR,mindist2outerR;
  Double_t pmtintegral,anodeq,vetoq,fastv;
  Int_t xmin,ymin,xmax,ymax,npixels, neighbors;
  Int_t spark, nburnin, camnum,runnum,evtnum,trknum, isrbi23;

  bdat->Branch("theta",&theta,"theta/D");
  bdat->Branch("phi",&phi,"phi/D");
  bdat->Branch("E",&E,"E/D");
  bdat->Branch("range",&range,"range/D");
  bdat->Branch("maxpixel",&maxpixel,"maxpixel/D");
  bdat->Branch("clusterrms",&clusterrms,"clusterrms/D");
  bdat->Branch("clustermean",&clustermean,"clustermean/D");
  bdat->Branch("mindist2innerR",&mindist2innerR,"mindist2innerR/D");
  bdat->Branch("mindist2outerR",&mindist2outerR,"mindist2outerR/D");
  bdat->Branch("x",&x,"x/D");
  bdat->Branch("y",&y,"y/D");
  bdat->Branch("pmtintegral",&pmtintegral,"pmtintegral/D");
  bdat->Branch("anodeq",&anodeq,"anodeq/D");
  bdat->Branch("vetoq",&vetoq,"vetoq/D");
  bdat->Branch("fastv",&fastv,"fastv/D");
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
	   
	  /////////////////////////////////////////////////////////////////////////////////////////////
	  //
	  // Find the minimum track-inner radius and the minimum track-outer radius distances ...
	  //
	  // get the pixels corresponding to this track cluster ..
	  vector<int> pix = ev->cluster(cam)->getCluster(t);

	  double innerRadius=stitcher->innerRadius(cam);
	  double outerRadius=stitcher->outerRadius(cam);
	  double xCenter=stitcher->xCenter(cam);
	  double yCenter=stitcher->yCenter(cam);
	    
	  mindist2innerR=0.; mindist2outerR=0.;
	  
	  // loop over the track pixels
	  for (int ipix=0; ipix<pix.size(); ipix++) {
	    
	    Int_t binx,biny;
	    ev->cluster(cam)->getXYFromBinNo(pix[ipix],&binx,&biny,true);
	    //      cout << "binx: " << binx << " biny: " << biny << endl;
	    
	    // find the distance from this pixel of the track to the center of the anode
	    double pixel_distance_from_this_track_pixel_to_anode_center=TMath::Sqrt( TMath::Power( ( binx-xCenter ) , 2 ) +
		 								     TMath::Power( ( biny-yCenter ) , 2 ) );

	    // negative if the track is outside the radius in question
	    if(!mindist2innerR) mindist2innerR=(innerRadius-pixel_distance_from_this_track_pixel_to_anode_center);
	    if(!mindist2outerR) mindist2outerR=(outerRadius-pixel_distance_from_this_track_pixel_to_anode_center);

	    if( TMath::Abs( innerRadius-pixel_distance_from_this_track_pixel_to_anode_center ) < TMath::Abs(mindist2innerR) ){
	      mindist2innerR=(innerRadius-pixel_distance_from_this_track_pixel_to_anode_center);
	    }
	    if( TMath::Abs( outerRadius-pixel_distance_from_this_track_pixel_to_anode_center ) < TMath::Abs(mindist2outerR) ){
	      mindist2outerR=(outerRadius-pixel_distance_from_this_track_pixel_to_anode_center);
	    }
	  }
	  
	  // clean up (do I need this?  is this a pointer or data in memory that I created?)
	  pix.erase(pix.begin(),pix.end());
	  //
	  // End of finding the minimum track-inner radius and the minimum track-outer radius distances ...
	  //
	  /////////////////////////////////////////////////////////////////////////////////////////////

	  // get the PMT summed integral
	  pmtintegral=0.;
	  anodeq=0.; vetoq=0.; fastv=0.;

	  TObjArray* waveform_vectors = ev->waveform_vectors();
	  for(Int_t ich=0; ich<waveform_vectors->GetEntries(); ++ich){
	    
	    if( ((TObject*)(*waveform_vectors)[ich])->GetName()==TString("mesh") ){
	      FastWfVector *mesh=((FastWfVector*)(*waveform_vectors)[ich]);  
	      for(Int_t itrg=0; itrg<mesh->size(); ++itrg) {
		if(!itrg) fastv=mesh->at(itrg).at(0).getPeak();
		if(mesh->at(itrg).at(0).getPeak()>fastv) fastv=mesh->at(itrg).at(0).getPeak();
	      }
	    }

	    if( ((TObject*)(*waveform_vectors)[ich])->GetName()==TString("anode") ){
	      CspWfVector *anode=((CspWfVector*)(*waveform_vectors)[ich]);  
	      for(Int_t itrg=0; itrg<anode->size(); ++itrg) {
		if(!itrg) anodeq=anode->at(itrg).at(0).getPeak();
		if(anode->at(itrg).at(0).getPeak()>anodeq) anodeq=anode->at(itrg).at(0).getPeak();
	      }
	    }

	    if( ((TObject*)(*waveform_vectors)[ich])->GetName()==TString("veto") ){
	      CspWfVector *veto=((CspWfVector*)(*waveform_vectors)[ich]);  
	      for(Int_t itrg=0; itrg<veto->size(); ++itrg) {
		if(!itrg) vetoq=veto->at(itrg).at(0).getPeak();
		if(veto->at(itrg).at(0).getPeak()>vetoq) vetoq=veto->at(itrg).at(0).getPeak();
	      }
	    }

	    if( ((TObject*)(*waveform_vectors)[ich])->IsA()->GetName()==TString("PMTWfVector") ){
	      /*
	      cout << "ich="
		   << ich 
		   << " channel="
		   << ((TObject*)(*waveform_vectors)[ich])->GetName()
		   << " class="
		   << ((TObject*)(*waveform_vectors)[ich])->IsA()->GetName()
		   << endl;

	      cout << "((PMTWfVector*)(*waveform_vectors)[ich])->size()=" 
		   << ((PMTWfVector*)(*waveform_vectors)[ich])->size()
		   << endl;
	      */

	      PMTWfVector *pmt=((PMTWfVector*)(*waveform_vectors)[ich]);  
	      //	      cout << "pmt->size()=" << pmt->size() << endl;
	      for(Int_t itrg=0; itrg<pmt->size(); ++itrg) {
		pmtintegral+=pmt->at(itrg).at(0).getIntegral();
	      }
	    }

	  }

	  bdat->Fill();
	  
	  out.add("4sh",ev->runNumber(),ev->eventNumber(),cam,t); 
	}
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
  //  nrrvseroicut->Draw("lp");
  
  
  //  cut_tree->Write(); 
  bdat->Write();
  //  app.Run(); 
}
