#include "DmtpcSkimDataset.hh"
#include "TDatime.h"
#include "TTimeStamp.h"
#include "TChain.h"
#include "TMath.h"
#include "TH2F.h"
#include "Dmtpc4ShooterStitcher.hh"
#include "TH1F.h"
#include "TSystem.h"
#include "TObjString.h"
#include "TKey.h"
#include <iostream>
#include "MaxCamClusterImage.hh"
#include "DmtpcStringTools.hh"
#include <assert.h>
using std::cout;
using std::endl;
using std::cerr;

using namespace std;

ClassImp(DmtpcSkimDataset)
DmtpcSkimDataset::DmtpcSkimDataset(unsigned max) 
{ 
  _biasFrames = 0; 
  _file = 0; 
  _max = max; 

  _events = (DmtpcSkimEvent**) malloc(max * sizeof(DmtpcSkimEvent*)); 
  _trees = (TChain**) malloc(max * sizeof(TChain*)); 

  for (unsigned i = 0; i < max; i++)
  {
    _events[i] = 0; 
    _trees[i] = 0; 
  }

  _gainMaps =0; 
  bias_loaded = false; 
  orig_loaded = false; 
  d = new DmtpcDataset;
  _config = 0; 

  ntrees = 0; 

}

DmtpcSkimDataset::DmtpcSkimDataset(const DmtpcSkimDataset &other) : TObject() {
  cout << GetName() << "Copy Constructor not done" << endl;
  _file = other._file;

}

DmtpcSkimDataset::~DmtpcSkimDataset() 
{

  //std::cout << "Inside Destructor!" << std::endl; 



    if (_file){
        closeRootFile(); 
    }
    delete d; 

    free(_events); 
    free(_trees); 

}

void DmtpcSkimDataset::loadDmtpcEvent(bool load , const char * orig)
{
 
  if (!_file || !_file->IsOpen())
  {
    cerr << "Trying to load DmtpcDataset before DmtpcSkimDataset. Nothing done." <<endl;
    return; 
  }

  //If told to unload but nothing's open, just return
  if (!load and !orig_loaded) return; 

  //Unload if told to unload and already loaded
  else if (!load && orig_loaded)
  {
    d->closeRootFile(); 
    orig_loaded = false; 
    orig = 0; 
    return; 
  }

  if (orig == NULL && !orig_loaded)
  {
    //Figure out the file name 
    TString n; 
    n+=TString(DEFAULT_DATA_DIR); 
    TString f = TString(
          DmtpcStringTools::basename(getFileName())
        ).ReplaceAll("skim","");

    //Figure out if the detector has a det tag 
    //if it does
    //the file name will be of the form
    //dmtpc_tag_xxxxx.root
    if (f.Data()[5] == '_')
    {
      for (int i = 6; i < 9; i++)
      {
        n+=f.Data()[i]; 
      }
      n+= "/"; 
    }

    n+=f;

    d->openRootFile(n.Data()); 
  }

  else if (!orig_loaded) d->openRootFile(orig); 

  orig_loaded = true; 
  d->getEvent(current_index); 

  return; 
}

void
DmtpcSkimDataset::getEvent(int i) {

  for (unsigned j = 0; j < ntrees; j++)
  {
    TChain * t =  _trees[j]; 
    if (i<0 || i>=t->GetEntries()) {
        cerr << GetName() << ": event number incorrect! " << i << endl;
        return; 
    }
     
    t->GetEntry(i);
  }

  if (orig_loaded) d->getEvent(i); 
  current_index = i; 
}
  


int 
DmtpcSkimDataset::getEventNumber(int n, int startingGuess)
{
  int i;
  if (startingGuess > 0) i = startingGuess;
  else i = n;

  if (i < 0) i = 0; 
  if (i >= nevents()) i = nevents() - 1; 
  getEvent(i);
  while (event()->eventNumber()!=n)
  {
    int last = event()->eventNumber();
    if (last < n)
    {
      getEvent(++i); 
    }
    else
    {
      getEvent(--i); 
    }
      
    int now = event()->eventNumber(); 
    //cout << "last, now = " << last << ", " << now << endl;
    if ( (last < n && now > n) || (last > n && now < n) || (now==last)) return -1; 
  }
  return i; 
}

void 
DmtpcSkimDataset::clearEventMemory() {
}


void DmtpcSkimDataset::mergeTempTreesv5(TTree * tmpskim, list<vector<vector<vector<BurninEncoded_t> > > * > * burnin, list<vector<pair<int,int> >*> * sparkref, int run_n)
{
  if(!_file->IsOpen())
  {
    cerr << "No out file is open." << endl;  
    return; 
  }
  if (tmpskim->GetEntries() != burnin->size())
  {
    cerr << "Two trees have non-matching dimensions. The merging is not complete" << endl;
    return; 
  } 

  int n = tmpskim->GetEntries(); 
  cout << "n=" << n << endl;
 
  DmtpcSkimEvent * _event = new DmtpcSkimEvent;

  tmpskim->SetBranchAddress("tmp_event",&_event); 

  _file->cd();

  TTree * imageTree = new TTree(tmpskim->GetName(), "Skimmed DMTPC events");
  TObjArray* waveform_vectors = NULL;

  imageTree->Bronch("event", "DmtpcSkimEvent", &_event, 128000, 1);  

  for (int i = 0; i < n; i++)
  {
    cout << "Merging Event: " << i << endl;

    tmpskim->GetEvent(i);

    //copypasta from Shawn 
    // waveform_vectors appears to get saved twice... I guess that's fine. 
    if(!waveform_vectors)
    {
      cout << "_event->_waveform_vectors=" << _event->_waveform_vectors << endl;
      cout << "_event->_waveform_vectors->GetEntries()=" << _event->_waveform_vectors->GetEntries() << endl; 
      cout << "generating the waveform_vectors TObjArray in the skim tree" << endl;
      waveform_vectors = new TObjArray(_event->_waveform_vectors->GetEntries());
      for(Int_t ich=0; ich<_event->_waveform_vectors->GetEntries(); ++ich)
      {
        
        waveform::tools::addWaveformVectorToTObjArray( waveform_vectors ,
                                                       ich , 
                                                       (*_event->_waveform_vectors)[ich]->GetName() ,
                                                       ((TObject*)(*_event->_waveform_vectors)[ich])->IsA()->GetName() );
        
      }
      //      waveform_vectors->AddAt( new WaveformVector(0,0,0,(*_event->_waveform_vectors)[ii]->GetName()), ii);
      imageTree->Branch(waveform_vectors,32000,1);
    }
    
    for(Int_t ii=0; ii<_event->_waveform_vectors->GetEntries(); ++ii) (*waveform_vectors)[ii]=(*_event->_waveform_vectors)[ii];


    //end copypasta

    _event->copyBurnin(burnin->front()); 
    delete burnin->front();
    burnin->pop_front(); 

    _event->copySparkRef(sparkref->front()); 
    delete[] sparkref->front();
    sparkref->pop_front(); 

    imageTree->Fill();
  }
 

  imageTree->Write(); 

  delete imageTree; 

}


void 
DmtpcSkimDataset::mergeTempTrees( TTree * tmpskim, list<vector<vector<vector<BurninEncoded_t> > > * > * burnin, list<vector<pair<int,int> >*> * sparkref, int run_n)
{
  if(!_file->IsOpen())
  {
    cerr << "No out file is open." << endl;  
    return; 
  }
  if (tmpskim->GetEntries() != burnin->size())
  {
    cerr << "Two trees have non-matching dimensions. The merging is not complete" << endl;
    return; 
  } 

  int n = tmpskim->GetEntries(); 
  cout << "n=" << n << endl;
 
  //We need these first
  int ncam; 
  tmpskim->SetBranchAddress("ncamera",&ncam);
  tmpskim->GetEvent(0); 
  
  DmtpcSkimEvent * _event = new DmtpcSkimEvent;
  //We need pointer to pointers to memory, so define these here 
  double  theta[ncam][15] ; bool thetaBr=0;
  double timenow[ncam]; bool timenowBr=0;
  double  range[ncam][15] ; bool rangeBr=0;
  double  diffusedrange[ncam][15] ; bool diffusedrangeBr=0;
  double  majoraxis[ncam][15] ; bool majoraxisBr=0;
  double  minoraxis[ncam][15] ; bool minoraxisBr=0;
  double  phi[ncam][15] ; bool phiBr=0;
  double  E[ncam][15] ; bool EBr=0;
  double  EGainMap[ncam][15] ; bool EGainMapBr=0;
  double  x[ncam][15] ; bool xBr=0;
  double  y[ncam][15] ; bool yBr=0;
  int  xbegin[ncam][15] ; bool xbeginBr=0;
  int  xend[ncam][15] ; bool xendBr=0;
  int  ybegin[ncam][15] ; bool ybeginBr=0;
  int  yend[ncam][15] ; bool yendBr=0;
  double  skewness[ncam][15] ; bool skewnessBr=0;
  double  cluster_rms[ncam][15] ; bool cluster_rmsBr=0;
  double  cluster_mean[ncam][15] ; bool cluster_meanBr=0;
  double  maxpixel[ncam][15] ; bool maxpixelBr=0;
  double  cygnus_angle[ncam][15] ; bool cygnus_angleBr=0;
  bool    edge[ncam][15] ; bool edgeBr=0;
  double  integral[ncam] ; bool integralBr=0;
  bool    spark[ncam] ; bool sparkBr=0;
  int     lastspark[ncam]; bool lastsparkBr=0;
  int     ntracks[ncam] ; bool ntracksBr=0;
  int     neighbors[ncam][15] ; bool neighborsBr=0;
  int     npixel[ncam][15] ; bool npixelBr=0;
  int     npixel_red[ncam][15] ; bool npixel_redBr=0;
  double  energy_density[ncam][15] ; bool energy_densityBr=0;
  double  ra[ncam][15] ; bool raBr=0;
  double  dec[ncam][15] ; bool decBr=0;
  double  glon[ncam][15] ; bool glonBr=0;
  double  glat[ncam][15] ; bool glatBr=0;
  double  moments[ncam][4][15] ; bool momentsBr=0;
  double  transverse_moments[ncam][4][15] ; bool transverse_momentsBr=0;
  double  image_mean[ncam] ; bool image_meanBr=0;
  double  image_rms[ncam] ; bool image_rmsBr=0;
  int  pixels_killed[ncam] ; bool pixels_killedBr=0;
  vector < string >* cameraSerialNumber=0; bool cameraSerialNumberBr=0;

  /* waveform stuff because idk how else to accomplish this
  int leftint[ncam] ; bool leftintBr=0;
  int leftint2[ncam] ; bool leftint2Br=0;
  int rightint[ncam] ; bool rightintBr=0;
  int rightint2[ncam] ; bool rightint2Br=0;
  int jmin[ncam] ; bool jminBr=0;
  int jbragg[ncam] ; bool jbraggBr=0;
  int jterm[ncam] ; bool jtermBr=0;
  int jterm2[ncam] ; bool jterm2Br=0;
  int jterm3[ncam] ; bool jterm3Br=0;
  int origin[ncam] ; bool originBr=0;
  vector<int> termmatch[ncam] ; bool termmatchBr=0;
  double rec_SD[ncam] ; bool rec_SDBr=0;
  int wfd_delta[ncam] ; bool wfd_deltaBr=0;
  int peak1[ncam] ; bool peak1Br=0;
  int peak2[ncam] ; bool peak2Br=0;
  double peak1val[ncam] ; bool peak1valBr=0;
  double peak2val[ncam] ; bool peak2valBr=0;
  int half1[ncam] ; bool half1Br=0;
  int half2[ncam] ; bool half2Br=0;
  double termdist[ncam] ; bool termdistBr=0;
  int maxdevloc[ncam] ; bool maxdevlocBr=0;
  double rms_left[ncam] ; bool rms_leftBr=0;
  double rms_right[ncam] ; bool rms_rightBr=0;
  double rms_outer[ncam] ; bool rms_outerBr=0;
  double rms_full[ncam] ; bool rms_fullBr=0;
  double dt[ncam] ; bool dtBr=0;

  if(tmpskim->FindBranch("leftint"))
    {tmpskim->SetBranchAddress("leftint", leftint); leftintBr=true;}
  if(tmpskim->FindBranch("leftint2"))
    {tmpskim->SetBranchAddress("leftint2",leftint2);leftint2Br=true;}
  if(tmpskim->FindBranch("rightint"))
    {tmpskim->SetBranchAddress("rightint",rightint);rightintBr=true;}
  if(tmpskim->FindBranch("rightint2"))
    {tmpskim->SetBranchAddress("rightint2",rightint2);rightint2Br=true;}
  if(tmpskim->FindBranch("jmin"))
    {tmpskim->SetBranchAddress("jmin",jmin);jminBr=true;}
  if(tmpskim->FindBranch("jbragg"))
    {tmpskim->SetBranchAddress("jbragg",jbragg);jbraggBr=true;}
  if(tmpskim->FindBranch("jterm"))
    {tmpskim->SetBranchAddress("jterm",jterm);jtermBr=true;}

  //Not yet finished... do I even need this stuff here? Not sure what I needed this for in the first place, actually.
  */

  /* end waveform stuff */

  if(tmpskim->FindBranch("theta"))
    {tmpskim->SetBranchAddress("theta",  theta); thetaBr=true;}

  if(tmpskim->FindBranch("timenow"))
  {tmpskim->SetBranchAddress("timenow",  timenow); timenowBr=true;}

  if(tmpskim->FindBranch("range"))
  {     tmpskim->SetBranchAddress("range",  range); rangeBr=true;}
  if(tmpskim->FindBranch("diffusedrange"))
  { tmpskim->SetBranchAddress("diffusedrange",  diffusedrange); diffusedrangeBr=true;}
  if(tmpskim->FindBranch("majoraxis"))
  {tmpskim->SetBranchAddress("majoraxis", majoraxis); majoraxisBr=true;}
  if(tmpskim->FindBranch("minoraxis"))
  {tmpskim->SetBranchAddress("minoraxis", minoraxis); minoraxisBr=true;}
  if(tmpskim->FindBranch("phi"))
  {tmpskim->SetBranchAddress("phi",  phi); phiBr=true;}
  if(tmpskim->FindBranch("E"))
  {tmpskim->SetBranchAddress("E",  E); EBr=true;}
  if(tmpskim->FindBranch("EGainMap"))
  {tmpskim->SetBranchAddress("EGainMap",  EGainMap); EGainMapBr=true;}
  if(tmpskim->FindBranch("x"))
  {tmpskim->SetBranchAddress("x",  x); xBr=true;}
  if(tmpskim->FindBranch("y"))
  {tmpskim->SetBranchAddress("y",  y); yBr=true;}

   if(tmpskim->FindBranch("xbegin"))
  {tmpskim->SetBranchAddress("xbegin",  xbegin); xbeginBr=true;}

   if(tmpskim->FindBranch("xend"))
  {tmpskim->SetBranchAddress("xend",  xend); xendBr=true;}

   if(tmpskim->FindBranch("ybegin"))
  {tmpskim->SetBranchAddress("ybegin",  ybegin); ybeginBr=true;}

   if(tmpskim->FindBranch("yend"))
  {tmpskim->SetBranchAddress("yend",  yend); yendBr=true;}

  if(tmpskim->FindBranch("skewness"))
  {tmpskim->SetBranchAddress("skewness",  skewness); skewnessBr=true;}
  if(tmpskim->FindBranch("cluster_rms"))
  {tmpskim->SetBranchAddress("cluster_rms", cluster_rms); cluster_rmsBr=true;}
  if(tmpskim->FindBranch("cluster_mean"))
  {tmpskim->SetBranchAddress("cluster_mean", cluster_mean); cluster_meanBr=true;}
  if(tmpskim->FindBranch("maxpixel"))
  {tmpskim->SetBranchAddress("maxpixel", maxpixel); maxpixelBr=true;}
  if(tmpskim->FindBranch("cygnus_angle"))
  {tmpskim->SetBranchAddress("cygnus_angle", cygnus_angle); cygnus_angleBr=true;}
  if(tmpskim->FindBranch("edge"))
  {tmpskim->SetBranchAddress("edge", edge); edgeBr=true;}
  if(tmpskim->FindBranch("integral"))
  {tmpskim->SetBranchAddress("integral", integral); integralBr=true;}
  if(tmpskim->FindBranch("spark"))
  {tmpskim->SetBranchAddress("spark", spark); sparkBr=true;}
  if(tmpskim->FindBranch("lastspark"))
  {tmpskim->SetBranchAddress("lastspark",  lastspark); lastsparkBr=true;}
  if(tmpskim->FindBranch("ntracks"))
  {tmpskim->SetBranchAddress("ntracks", ntracks); ntracksBr=true;}
  if(tmpskim->FindBranch("neighbors"))
  {tmpskim->SetBranchAddress("neighbors", neighbors); neighborsBr=true;}
  if(tmpskim->FindBranch("npixel"))
  {tmpskim->SetBranchAddress("npixel", npixel); npixelBr=true;}
  if(tmpskim->FindBranch("ncamera"))
  {tmpskim->SetBranchAddress("ncamera", &_event->_ncamera);}
  if(tmpskim->FindBranch("clusters"))
  {tmpskim->SetBranchAddress("clusters", &_event->_clusters); }
  if(tmpskim->FindBranch("energy_density"))
  {tmpskim->SetBranchAddress("energy_density",energy_density); energy_densityBr=true;}
  if(tmpskim->FindBranch("eventnum"))
  {tmpskim->SetBranchAddress("eventnum",&_event->_eventNumber); }
  if(tmpskim->FindBranch("ra"))
  {tmpskim->SetBranchAddress("ra",ra);   raBr=true;}
  if(tmpskim->FindBranch("dec"))
  {tmpskim->SetBranchAddress("dec",dec);   decBr=true;}
  if(tmpskim->FindBranch("image_mean"))
  {tmpskim->SetBranchAddress("image_mean",image_mean); image_meanBr=true;}
  if(tmpskim->FindBranch("image_rms"))
  {tmpskim->SetBranchAddress("image_rms",image_rms); image_rmsBr=true;}
  if(tmpskim->FindBranch("pixels_killed"))
  {tmpskim->SetBranchAddress("pixels_killed",pixels_killed); pixels_killedBr=true;}
  if(tmpskim->FindBranch("npixel_red"))
  {tmpskim->SetBranchAddress("npixel_red",npixel_red); npixel_redBr=true;}
  if(tmpskim->FindBranch("glon"))
  {tmpskim->SetBranchAddress("glon",glon);   glonBr=true;}
  if(tmpskim->FindBranch("glat"))
  {tmpskim->SetBranchAddress("glat",glat);   glatBr=true;}
  if(tmpskim->FindBranch("moments"))
  {tmpskim->SetBranchAddress("moments",moments);   momentsBr=true;}
  if(tmpskim->FindBranch("transverse_moments"))
  {tmpskim->SetBranchAddress("transverse_moments",transverse_moments);   transverse_momentsBr=true;}
  if(tmpskim->FindBranch("trigger_groups"))
  {tmpskim->SetBranchAddress("trigger_groups", &_event->_trigger_groups); }
  if(tmpskim->FindBranch("waveform_vectors"))
  {tmpskim->SetBranchAddress("waveform_vectors", &_event->_waveform_vectors); }
  if(tmpskim->FindBranch("cameraSerialNumber"))
  {tmpskim->SetBranchAddress("cameraSerialNumber",&cameraSerialNumber); cameraSerialNumberBr=true;}
  
  _file->cd();
  
  TTree * imageTree = new TTree(tmpskim->GetName(), "Skimmed DMTPC events");
  TObjArray* waveform_vectors = NULL;

  imageTree->Branch("event", "DmtpcSkimEvent", &_event, 128000, 1);  

  _event->initVectors(ncam);
  for (int i = 0; i < n; i++)
  {

    cout << "Merging Event: " << i << endl;

    tmpskim->GetEvent(i);

    if(tmpskim->FindBranch("waveform_vectors")){
      if(!waveform_vectors){
        cout << "_event->_waveform_vectors=" << _event->_waveform_vectors << endl;
        cout << "_event->_waveform_vectors->GetEntries()=" << _event->_waveform_vectors->GetEntries() << endl; 
        cout << "generating the waveform_vectors TObjArray in the skim tree" << endl;
        waveform_vectors = new TObjArray(_event->_waveform_vectors->GetEntries());
        for(Int_t ich=0; ich<_event->_waveform_vectors->GetEntries(); ++ich){
          
          waveform::tools::addWaveformVectorToTObjArray( waveform_vectors ,
                                                         ich , 
                                                         (*_event->_waveform_vectors)[ich]->GetName() ,
                                                         ((TObject*)(*_event->_waveform_vectors)[ich])->IsA()->GetName() );
          
        }
        //      waveform_vectors->AddAt( new WaveformVector(0,0,0,(*_event->_waveform_vectors)[ii]->GetName()), ii);
        imageTree->Branch(waveform_vectors,32000,1);
      }
      
      for(Int_t ii=0; ii<_event->_waveform_vectors->GetEntries(); ++ii) (*waveform_vectors)[ii]=(*_event->_waveform_vectors)[ii];
    }
    
    _event->_runNumber = run_n;
    _event->_index = i; 

    //Copy from array to vectors
    
    for(int c = 0; c< ncam; c++)
    {

      if(integralBr)
        _event->_integral[c] = integral[c]; 
      if(sparkBr)
        _event->_spark[c] = spark[c];
      if(lastsparkBr)
        _event->_lastspark[c] = lastspark[c];
      if(ntracksBr)
        _event->_ntracks[c] = ntracks[c]; 
      if(image_rmsBr)
        _event->_image_rms[c] = image_rms[c]; 
       
      if(timenowBr)
        _event->_timenow[c] = timenow[c]; 
 
      if(image_meanBr)
        _event->_image_mean[c] = image_mean[c]; 
      if(pixels_killedBr)
        _event->_pixels_killed[c] = pixels_killed[c]; 
      if(cameraSerialNumberBr)
        _event->_cameraSerialNumber[c] = (*cameraSerialNumber)[c];
      for (int t = 0; t< 15 ; t++)
        {
          if(thetaBr)
            _event->_theta[c][t] = theta[c][t]; 
          if(rangeBr)
            _event->_range[c][t] = range[c][t]; 
          if(diffusedrangeBr)
            _event->_diffusedRange[c][t] = diffusedrange[c][t];
          if(majoraxisBr)
            _event->_majoraxis[c][t] = majoraxis[c][t];
          if(minoraxisBr)
            _event->_minoraxis[c][t] = minoraxis[c][t];
         if(phiBr)
           _event->_phi[c][t] = phi[c][t]; 
         if(EBr)
            _event->_E[c][t] = E[c][t]; 
         if(EGainMapBr)
            _event->_EGainMap[c][t] = EGainMap[c][t]; 
         if(xBr)
            _event->_x[c][t] = x[c][t]; 
         if(yBr)
            _event->_y[c][t] = y[c][t]; 
	 if(xbeginBr)
	   _event->_xbegin[c][t] = xbegin[c][t];
         if(ybeginBr)
	   _event->_ybegin[c][t] = ybegin[c][t];
	 if(xbeginBr)
	   _event->_xend[c][t] = xend[c][t];
	 if(yendBr)
	   _event->_yend[c][t] = yend[c][t];
	 if(skewnessBr)
            _event->_skewness[c][t] = skewness[c][t]; 
         if(cluster_rmsBr)
            _event->_cluster_rms[c][t] = cluster_rms[c][t]; 
         if(cluster_meanBr)
            _event->_cluster_mean[c][t] = cluster_mean[c][t]; 
         if(maxpixelBr)
            _event->_maxpixel[c][t] = maxpixel[c][t]; 
         if(cygnus_angleBr)
            _event->_cygnus_angle[c][t] = cygnus_angle[c][t]; 
         if(edgeBr)
            _event->_edge[c][t] = edge[c][t]; 
         if(neighborsBr)
            _event->_neighbors[c][t] = neighbors[c][t]; 
         if(npixelBr)
            _event->_npixel[c][t] = npixel[c][t]; 
         if(npixel_redBr)
            _event->_npixel_red[c][t] = npixel_red[c][t]; 
         if(energy_densityBr)
            _event->_energy_density[c][t] = energy_density[c][t]; 
         if(raBr)
            _event->_ra[c][t] = ra[c][t]; 
         if(decBr)
            _event->_dec[c][t] = dec[c][t]; 
         if(glonBr)
            _event->_glon[c][t] = glon[c][t]; 
         if(glatBr)
            _event->_glat[c][t] = glat[c][t]; 

        for (int m = 0; m < 4; m++)
        {
           if(momentsBr)
              _event->_moments[c][m][t] = moments[c][m][t]; 
           if(transverse_momentsBr)
              _event->_transverse_moments[c][m][t] = transverse_moments[c][m][t]; 
        }
      }
   }

    _event->copyBurnin(burnin->front()); 
    delete burnin->front();
    burnin->pop_front(); 

    _event->copySparkRef(sparkref->front()); 
    delete[] sparkref->front();
    sparkref->pop_front(); 

    imageTree->Fill();
  }

  imageTree->Write(); 

  delete imageTree;
 
}


void DmtpcSkimDataset::writeStitch(Dmtpc4ShooterStitcher * stitch)
{

  //only write once 

  if (_file->Get("stitch"))
  {
     cerr << "Dataset already contains stitch object!!! Not Writing! "<< endl;  
     return; 
  }

  _file->cd(); 
  stitch->Write("stitch",1); 
  gROOT->cd(); 
}

void DmtpcSkimDataset::writeGainMaps()
{
  //Only write once
  if (_file->Get("gainMaps"))
    return; 
  _file->cd();
  _gainMaps->Write("gainMaps",1);
  gROOT->cd();
}

void 
DmtpcSkimDataset::createRootFile( const char *fname, TString foption)
{
    foption.ToLower();    

    if (_file){
      //if (_file->IsOpen()) _file->Close();
      //delete _file;

      closeRootFile(); 

    }
    //
    // write
    //
    if (fname && foption=="create") {
      _file = TFile::Open(fname, foption);
      if(!_file->IsOpen())
      {
        cerr << "Error opening file! Files cannot be overwritten.\n";
      }
      _gainMaps = new TObjArray;  

    }
    //
    // read
    //
    else if (fname && foption=="") {
       
       _file = TFile::Open(fname);
       _file->cd(); 
       TIter next(_file->GetListOfKeys());
       while(TKey * key = (TKey*) next.Next())
       {
          TObject * obj = _file->Get(key->GetName());   
          if (!strcmp(obj->IsA()->GetName(),"TTree"))
          {
             if (_indices.count(key->GetName())) continue; 

         //    std::cout << " event address of " << key->GetName() << _events[ntrees] << std::endl; 
         //    std::cout << " chain address of " << key->GetName() << _trees[ntrees] << std::endl; 
             gROOT->cd(); 
             TChain * ch = new  TChain(key->GetName(),"Reconstructed Tree"); 
             ch->SetDirectory(0); 
             _trees[ntrees] = ch; 
         //    std::cout << " chain address of " << key->GetName() << ch << std::endl; 
         //    std::cout << " event address of " << key->GetName() << _events[ntrees] << std::endl; 
             ch->SetBranchAddress("event",&(_events[ntrees])); 
         //    std::cout << " event address of " << key->GetName() << _events[ntrees] << std::endl; 
             ch->Add(fname); 
             _nevents = ch->GetEntries(); 
             _indices[std::string(key->GetName())] = ntrees++; 
          }

          obj->Delete(); 
       }
        
       _gainMaps = (TObjArray*)_file->Get("gainMaps");
       if (_file->Get("config"))
       {
         _config = (TObjString*) _file->Get("config"); 
       }
    }
    else if (fname) assert(0);
}

DmtpcGainMap* DmtpcSkimDataset::getGainMap(TString serialNumber)
{
   for(int i=0; i<_gainMaps->GetEntries(); i++)
   {
      if(((DmtpcGainMap*)_gainMaps->At(i))->GetName()==serialNumber)
      {
         return (DmtpcGainMap*)_gainMaps->At(i);
      }
   }
   cout << "No gain map found with serial number " << serialNumber
        << ". \n";
   return 0;

}

const char* DmtpcSkimDataset::getFileName() { return _file->GetName(); }

void 
DmtpcSkimDataset::openRootFile( const char *fname) {
  createRootFile(fname,"");
}

void 
DmtpcSkimDataset::loadBiasFrames(bool load, const char *fname)
{
  
  /* If bias not loaded and asked to not load, do nothing */
  if (!bias_loaded && !load) return; 
  
  /* Otherwise, we either want to delete because we no longer 
   * need the bias frames or because we want a different file */
  if (bias_loaded)
  {
    bias_loaded = false;
    for (int i = 0; i < nbias; i++)
    {
      _biasFrames[i]->Delete(); 
    } 
    
    delete _biasFrames; 
  }
  
  
  TFile * biasfile; 
  if (strcmp(fname,"")==0)
  {
    //figure out bias file from file name
    TString fn = TString(getFileName()).ReplaceAll("skim.root","bias.root");
    biasfile = new TFile(fn.Data()); 
  }
  else
  {
      biasfile = new TFile(fname);  
  }
  
  if (biasfile == 0)
  {
    cerr << fname << "could not be loaded" << endl; 
    return;
  }
  
  gROOT->cd(); 
  
  TTree * bias = (TTree *) biasfile->Get("bias"); 
  TH2F * frame = 0; 
  bias->SetBranchAddress("biasframe",&frame);
  
  nbias = bias->GetEntries(); 
  _biasFrames = new TH2F*[nbias]; 
  
  for (int i = 0; i < nbias; i++)
  {
    bias->GetEvent(i);
    _biasFrames[i] = frame;
    //cout << frame->GetName() << " " << frame << endl; 
    frame = 0; 
  }
   
   bias->Delete();
   biasfile->Close(); 
   biasfile->Delete(); 
   bias_loaded = true; 
}


void
DmtpcSkimDataset::closeRootFile()
{
  if (!_file) return;
  if (!_file->IsOpen()) return;
  if (orig_loaded) d->closeRootFile(); 

    for (unsigned i = 0; i <ntrees;i++)
  {
    //std::cout << "Deleting " << _trees[i] << std::endl; 
    _trees[i]->Delete();
    _trees[i] = 0; 
    delete _events[i]; 
    _events[i] = 0; 
  }
  _file->Close();
  delete _file;
  _indices.clear(); 
  ntrees = 0; 

 if (_gainMaps)
  {
    _gainMaps->Delete(); 
    delete _gainMaps; 
    _gainMaps = 0; 
  }
 _file = 0; 

 if (_config)
 {
    _config->Delete(); 
 } 

}


void
DmtpcSkimDataset::newRootFile(const char *fname)
{
  //spitz 
  createRootFile(fname,"create");
}

