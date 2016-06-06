#include "../../MaxCam/MaxCamMC.hh"
#include "../../MaxCam/MaxCamImageTools.hh"
#include "../../MaxCam/MaxCamCluster.hh"
#include "../../MaxCam/DmtpcEvent.hh"
#include "../../MaxCam/MaxCamConfig.hh"
#include "../../MaxCam/DmtpcDataset.hh"
#include "../../MaxCam/DmtpcKeys.hh"
#include "../../MaxCam/DmtpcSkimDataset.hh"
#include "../../MaxCam/DmtpcSkimEvent.hh"
#include "../../MaxCam/DmtpcSkimRunSummary.hh"
#include "../../MaxCam/MaxCamTriggerGroup.hh"

#include "TBenchmark.h"
#include "TF1.h"
#include "TSystem.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TProfile.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TList.h"
#include "TCutG.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <list>
#include <stdlib.h>
#include <sstream>
#include <stdio.h>
#include "cleanSkimConfig.hh"
#include "cleanSkimProcessor.hh"
#include "TGraphPolar.h"
#include "TText.h"

using namespace std;

/********* Global Variables ******/ 

//Printout levels: | together 

static const u_int32_t dbg_error     =  1;     //Fatal unexpected condition printouts (these go to cerr)
static const u_int32_t dbg_warning   =  2;     //Non-fatal unexpected condition printouts
static const u_int32_t dbg_progress  =  4;     //Progress printouts
static const u_int32_t dbg_bias      =  8;     //Bias algorithm printouts
static const u_int32_t dbg_cluster   =  16;    //Cluster finding printouts
static const u_int32_t dbg_burnin    =  32;    //Burnin printouts
static const u_int32_t dbg_config    =  64;    //Config printout

static u_int32_t dbg_level = dbg_progress | dbg_error | dbg_warning; 


/******** Prototypes **********/
static int parse_args( int nargs, char ** args, TString * files, TString * keys, TString * out, char ** );
static bool burnin_test( double xdelta, double ydelta, CleanSkimConfig * conf);
int Convert(TString,int,int,DmtpcDataset & d );
void initGUI(TCanvas*);
//static void refreshTime();
void refreshPhiHistogram(vector<Double_t>, vector<Double_t>,TGraphPolar*, TCanvas*);

/******** end prototypes ******/

/******** global variables for mturek ************/
/*TString folder = "/skim";
Bool_t skipInitialContents = false;
int maxTracks = 50;
// Data structures:
vector<TString> processedFiles;
DmtpcSkimDataset sd;
Int_t EventsRead = 0;
vector<Double_t> phi;
vector<Double_t> radii;
// GUI elements:
TCanvas * canvas;
TText * currentTimeLabel;
TH1D * phiHistogram;
TGraphPolar * polarGraph;
TGraphPolar * polarArrow;
TText * numberDetected;
TText * averagePhi;
TText * rms;*/
TText currentTimeLabel;
TText numberDetected;
TText averagePhi;
TText rms;

/******************* end global variables ******/

/******** cleanSkim *********/
int cleanSkim(DmtpcDataset & d, TString & key, TTree * rawtree, int runnum, 
              DmtpcSkimDataset & sd, bool delete_intermediate, CleanSkimConfig * conf, unsigned algo_index)
{

  if(dbg_level & dbg_error){
  std::cout <<"hasOverscan(): " << conf->hasOverscan() << std::endl; 
  }

  TBenchmark* myBench=new TBenchmark();

  TString tmpoutfilename = TString(sd.getFileName()).ReplaceAll(".root","_tmp.root");
  if(dbg_level & dbg_error) cout << tmpoutfilename << endl;

  /* Define Temporary out files and trees */
  TFile* tmpoutfile = new TFile(tmpoutfilename,"RECREATE");
  TTree* tmpskimtree;
  TH2F* tempimg=0;
  TH2F* os_tempimg=0;

  /* Get some eventwide information */
  d.getEvent(0);
  const int ncamera = d.event()->ccdData()->GetEntries();
  const int nbinsx =  d.event()->ccdData(0)->GetNbinsX();
  const int nbinsy =  d.event()->ccdData(0)->GetNbinsY();
  const int eventVersion = d.event()->IsA()->GetClassVersion();
  int ncam = ncamera;

  if(dbg_level & dbg_error) cout << "ncamera " << ncamera << endl;

  /* Support processing on datasets with no actual images */ 
  bool fake_image = false; 
  if (ncamera ==0 || d.event()->ccdData(0)==NULL || d.event()->ccdData(0)->Integral() ==0)
  {
    fake_image = true; 
    tempimg = new TH2F ("tempimg","tempimg",256,0,1024,256,0,1024); 
  }

  if(dbg_level & dbg_error) cout << "Getting camera id" << endl;
  	
 //Get camera id's 
 string camera_ids[ncamera];  
 for (int u = 0; u < ncamera; u++)
 {
    camera_ids[u] = string( ((MaxCamConfig*)d.event()->ccdConfig(u))->serialNumber.Data()); 

    if(dbg_level & dbg_error) cout << "camera_ids " << camera_ids << endl;

    if (camera_ids[u] == "")
    {
      camera_ids[u] = std::string(conf->getFallBackCameraID(u)); 
      if (camera_ids[u]!="")
      {
        if(dbg_level & dbg_error) std::cout << "Fell back to camera id: " << camera_ids[u] << std::endl; 
      }
    }
 }

 if(dbg_level & dbg_error) cout << "Loading complete" << endl;

  /* Create the temporary skim tree without burnin data.  */
  tmpskimtree = new TTree(key,"Skimmed Events");
  double theta[ncamera][15], phi[ncamera][15], E[ncamera][15];
  double EGainMap[ncamera][15],diffusedrange[ncamera][15];
  double range[ncamera][15], x[ncamera][15], y[ncamera][15];
  double skewness[ncamera][15], integral[ncamera];
  bool edge[ncamera][15], spark[ncamera];
  int ntracks[ncamera], lastspark[ncamera],date, time, eventnum;
  int pixels_killed[ncamera], npixel_red[ncamera][15];
  double image_mean[ncamera],image_rms[ncamera];
  double os_mean[ncamera],os_rms[ncamera];
  MaxCamClusterImage* clusti=0;
  double cluster_rms[ncamera][15], cluster_mean[ncamera][15], energy_density[ncamera][15];
  int neighbors [ncamera][15], npixel[ncamera][15]; 
  double maxpixel[ncamera][15], cygnus_angle[ncamera][15]; 
  double moments[ncamera][4][15]; 
  double transverse_moments[ncamera][4][15]; 
  double ra[ncamera][15], dec[ncamera][15]; 
  double glat[ncamera][15],  glon[ncamera][15]; 
  DmtpcEvent* skimevent = new DmtpcEvent();
  TObjArray* clust = new TObjArray(ncamera);
  clust->SetOwner(kTRUE);  
  TObjArray* trigger_groups = new TObjArray;
  vector< string > * cameraSerialNumber = new vector< string >;
  double majoraxis[ncamera][15],minoraxis[ncamera][15];
  // make a fixed length array of waveform vectors, with only as many entries as there are scope channels
  if(dbg_level & dbg_error) cout << "Temporary tree created."  << endl;
  
  tmpskimtree->Branch("ncamera",&ncam,"ncamera/I");
  tmpskimtree->Branch("spark",&spark,"spark[ncamera]/O");
  tmpskimtree->Branch("lastspark",&lastspark,"lastspark[ncamera]/I");
  tmpskimtree->Branch("integral",&integral,"integral[ncamera]/D");
  tmpskimtree->Branch("eventnum",&eventnum,"eventnum/I");
  tmpskimtree->Branch("theta",theta,"theta[ncamera][15]/D");
  tmpskimtree->Branch("phi",phi,"phi[ncamera][15]/D");
  tmpskimtree->Branch("E",E,"E[ncamera][15]/D");
  tmpskimtree->Branch("EGainMap",EGainMap,"EGainMap[ncamera][15]/D");
  tmpskimtree->Branch("energy_density", energy_density, "energy_density[ncamera][15]/D");
  tmpskimtree->Branch("range",range,"range[ncamera][15]/D");
  tmpskimtree->Branch("majoraxis",majoraxis,"majoraxis[ncamera][15]/D");
  tmpskimtree->Branch("minoraxis",minoraxis,"minoraxis[ncamera][15]/D");
  tmpskimtree->Branch("diffusedrange",diffusedrange,"diffusedrange[ncamera][15]/D");
  tmpskimtree->Branch("x",x,"x[ncamera][15]/D");
  tmpskimtree->Branch("y",y,"y[ncamera][15]/D");
  tmpskimtree->Branch("ntracks",ntracks,"ntracks[ncamera]/I");
  tmpskimtree->Branch("edge",edge,"edge[ncamera][15]/O");
  tmpskimtree->Branch("skewness",skewness,"skewness[ncamera][15]/D");
  tmpskimtree->Branch("ra",ra,"ra[ncamera][15]/D");
  tmpskimtree->Branch("dec",dec,"dec[ncamera][15]/D");
  tmpskimtree->Branch("glat",glat,"glat[ncamera][15]/D");
  tmpskimtree->Branch("glon",glon,"glon[ncamera][15]/D");
  tmpskimtree->Branch("moments",moments,"moments[ncamera][4][15]/D");
  tmpskimtree->Branch("transverse_moments",transverse_moments,"transverse_moments[ncamera][4][15]/D");
  tmpskimtree->Branch("date",&date,"date/I");
  tmpskimtree->Branch("time",&time,"time/I");
  tmpskimtree->Branch("clusters","TObjArray",&clust,128000,0);
  tmpskimtree->Branch("cluster_rms",cluster_rms,"cluster_rms[ncamera][15]/D");
  tmpskimtree->Branch("cluster_mean",cluster_mean,"cluster_mean[ncamera][15]/D");
  tmpskimtree->Branch("maxpixel",maxpixel,"maxpixel[ncamera][15]/D");
  tmpskimtree->Branch("neighbors",neighbors,"neighbors[ncamera][15]/I");
  tmpskimtree->Branch("cygnus_angle", cygnus_angle, "cygnus_angle[ncamera][15]/D"); 
  tmpskimtree->Branch("npixel", npixel, "npixel[ncamera][15]/I"); 
  tmpskimtree->Branch("npixel_red",npixel_red,"npixel_red[ncamera][15]/I");
  tmpskimtree->Branch("pixels_killed",pixels_killed,"pixels_killed[ncamera]/I");
  tmpskimtree->Branch("image_mean",image_mean,"image_mean[ncamera]/D");
  tmpskimtree->Branch("image_rms",image_rms,"image_rms[ncamera]/D");
  //tmpskimtree->Branch("trigger_groups","TObjArray",&trigger_groups,128000,0);
  tmpskimtree->Branch("cameraSerialNumber",&cameraSerialNumber);
  if(dbg_level & dbg_error) cout << "Branches created."  << endl;

  gROOT->cd();


  /* Open up file for bias frame saving */
  TFile * biasoutfile;
  if (algo_index == 0)
  {
    TString biasoutfilename=TString(sd.getFileName()).ReplaceAll(key+".root","bias.root");
    //biasoutfilename.ReplaceAll(key+".root","bias.root");
    biasoutfile = new TFile(biasoutfilename,"RECREATE");
    if(dbg_level & dbg_error) cout << "Bias file created: " << biasoutfilename << endl;

  }
  

  TH2F* biasframe = d.event()->ccdData(0);
  TTree* biastree = new TTree("bias","Bias Information");
  biastree->Branch("biasframe","TH2F",&biasframe,128000,0);
  gROOT->cd();

  double sigmathr[ncamera];
  int nframes[ncamera];
  for(int u=0; u<ncamera; u++){nframes[u]=0;}
  TH2F* secondarybias[ncamera];
  TH2F* os_secondarybias[ncamera];
  double bias_mean[ncamera], bias_rms[ncamera],os_bias_mean[ncamera],os_bias_rms[ncamera];
  double cleaned_os_bias_mean[ncamera];
  double last_mean[ncamera]; // for the spark cut



  // handle the bias frame
  //      define the bias threshold for runs<500 and MC
  for(int u=0; u<ncamera; u++)
    {
      if (fake_image)
      {
         break; 
      }
      if (dbg_level & dbg_bias) std::cout << "Camera " << u << " image cleaning:" << std::endl;
      secondarybias[u] = (TH2F*)d.getBiasFrame(u+1)->Clone();
      
      // the cluster finding threshold is very sensitive to this number!
      // using 1.5 sigma reproduces cluster finding efficiency from AF v1
      sigmathr[u]=1.5*MaxCamImageTools::getRMS(secondarybias[u]);
      
      
      if (dbg_level & dbg_bias)
      {  
        if(dbg_level & dbg_error) std::cout << "     camera bias MEAN: " << 
             MaxCamImageTools::getMean(secondarybias[u]) << std::endl;
        if(dbg_level & dbg_error) std::cout << "     camera bias RMS: " << 
             MaxCamImageTools::getRMS(secondarybias[u]) << std::endl;
      }
    }

  const int nev = rawtree->GetEntries();

  if (dbg_level & dbg_error) cout << "Number of events: " << nev << endl;

  /* Bias frame stuff */
  for(int u=0; u<ncamera; u++)
    {
      if (fake_image) 
      {
        secondarybias[u] = new TH2F("bias","bias",256,0,1024,256,0,1024); 
      }
      else
      {
        if(nframes[u] !=0 && runnum < 500) {
          if (dbg_level & dbg_bias) std::cout << "run number " << runnum << " using averaged bias " << std::endl;
          secondarybias[u]->Scale(1/double(nframes[u]));
        } else {
          secondarybias[u] = (TH2F*)d.getBiasFrame(u+1)->Clone();
          if (conf->hasOverscan())
          {
            os_secondarybias[u]= (TH2F*)d.getBiasFrameOverscan(u+1)->Clone();
          }
          if (dbg_level & dbg_bias) std::cout << "run number " << runnum << " using camera bias " << std::endl;
        }
        // kill 5 sigma outliers in the bias frame, iterate 3 times
        for (int it=0; it<3; it++) {
           double mean, rms,osmean,osrms;
          mean=MaxCamImageTools::getMean(secondarybias[u]);
          rms=MaxCamImageTools::getRMS(secondarybias[u]);
          if (conf->hasOverscan())
          {
            osmean=MaxCamImageTools::getMean(os_secondarybias[u]);
            osrms=MaxCamImageTools::getRMS(os_secondarybias[u]);
          }
          if (conf->isMC()) { // MC
            MaxCamImageTools::killPixels(secondarybias[u],mean+10.*rms);
          } else {
            MaxCamImageTools::killPixels(secondarybias[u],mean*conf->getOutlierFactor());
            if (conf->hasOverscan())
            {
              MaxCamImageTools::killPixels(os_secondarybias[u],osmean*conf->getOutlierFactor());
            }
          }          if (dbg_level & dbg_bias) 
          {
            if(dbg_level & dbg_error) std::cout << "     iteration " << it << " averaged bias image mean (no outliers): " << mean << std::endl;
            if(dbg_level & dbg_error) std::cout << "     iteration " << it << " averaged bias image rms (no outliers): " << rms << std::endl;
          }
          bias_mean[u]=mean;
          bias_rms[u]=rms;
 
         if (conf->hasOverscan())
          {
	    os_bias_mean[u]=osmean;
            os_bias_rms[u]=osrms;
          }
        }

	if (conf->hasOverscan())
          {
	    // we use this later to correct how much images have drifted in overall mean 
	    // since the bias frame was taken.  It's important that we use the cleaned image
	    // for this mean because we might otherwise pick up a hot pixel
	    cleaned_os_bias_mean[u]=MaxCamImageTools::getMean(os_secondarybias[u]);
          }
      }
      biasframe = secondarybias[u];
      if (algo_index==0)
      {
	biasoutfile->cd();
        biastree->Fill();
      }
      gROOT->cd();   
    }
    
  if (dbg_level & dbg_progress) std::cout << "All preclean activities done" << std::endl;
  int nempty=0;

  /* Running list of positions so we don't have to read back through skim file
     all the time. It is a List of a Vector[2] of a double[2].*/
  list<vector<double *>**> positions; 
  

  /* List of spark ref excluded pixels */
  vector<pair<int,int> > * sparkref_running = new vector<pair<int,int> >[ncamera]; 
  list<vector<pair<int,int> >*>  sparkref; 

  list<vector<vector<vector<BurninEncoded_t> > >*> burnin_temp; //Temporary list that is updated a lot
  list<vector<vector<vector<BurninEncoded_t> > >*> burnin; //permanent list 
  
  int position_offset  = 0; 

  //define cut for sparks with overscan
  TCutG* os_cut;
  if (conf->hasOverscan())
  {
    os_cut = new TCutG("os_cut",5);
    os_cut->SetPoint(0,0,1028);
    os_cut->SetPoint(1,1024,1028);
    os_cut->SetPoint(2,1024,1032);
    os_cut->SetPoint(3,0,1032);
    os_cut->SetPoint(4,0,1028);
  }

  //Load gain maps if we are using them; 

  DmtpcGainMap * gainmaps[ncamera]; 

  if (conf->useGainMap())
  {
    TFile gfile(conf->getGainMapFile()); 
    gfile.ls(); 
    for (int u = 0; u < ncamera; u++)
    {
      std::string ser =  camera_ids[u]; 
      gainmaps[u]=  ser == "" ? 0 : (DmtpcGainMap*) gfile.Get(ser.c_str())->Clone();
      gainmaps[u]->getGainMap()->Scale(gainmaps[u]->getGainMap()->GetNbinsX()*gainmaps[u]->getGainMap()->GetNbinsY()/gainmaps[u]->getGainMap()->Integral());
      sd.addGainMap((DmtpcGainMap*)gainmaps[u]->Clone());
    }
    sd.writeGainMaps();
    gSystem->Sleep(500); 
    gfile.Close(); 
  }
  else
  {
    for (int u = 0; u < ncamera; u++)  
    {
      gainmaps[u] = 0; 
    }
  }

  // loop over events


  for(int i = 0; i<nev; i++)
  {
    myBench->Start("Event");

    if (dbg_level & dbg_progress) std::cout << "Processing event " << i << std::endl;
    d.getEvent(i);
    int sumtracks=0;
    int sparks=0;
    clust->Clear();
    trigger_groups->Clear(); 
    
    vector<vector<vector<BurninEncoded_t> > > * this_event_burnin= new vector<vector<vector<BurninEncoded_t> > >;
    myBench->Start("Saving scope data");
            
    myBench->Start("Loop over cameras");
    
    for(int u=0; u<ncamera; u++)
      {
        //save the camera serial numbers
        cameraSerialNumber->push_back(camera_ids[u]); 
	
        vector< vector<BurninEncoded_t> >  this_camera_burnin ;
        MaxCamClusterImage::CAMERA_ORIENTATION orientation = conf->getCameraOrientation(camera_ids[u].c_str());         
	
        if (dbg_level & dbg_progress) std::cout << "\t" << u << ":" << std::endl;
        // keep track of the last frame for the spark cut
        if (i==0) {
          last_mean[u] = bias_mean[u];
          lastspark[u]=0;
        }
	
        if(!fake_image)
	  {
	    tempimg = (TH2F*)d.event()->ccdData(u)->Clone("tempimg");
	    os_tempimg = conf->hasOverscan() ? (TH2F*)d.event()->overscan(u)->Clone("tempimg") : 0;
	  }
	
        double mean, rms, osmean, osrms;
        mean= fake_image ? 0 : MaxCamImageTools::getMean(tempimg);
        rms= fake_image ? 0 : MaxCamImageTools::getRMS(tempimg);
        osmean= fake_image  || !conf->hasOverscan() ? 0 : MaxCamImageTools::getMean(os_tempimg);
        osrms= fake_image || !conf->hasOverscan() ? 0 : MaxCamImageTools::getRMS(os_tempimg);
	
        double mean_nooutliers = 0; 
        double rms_nooutliers = 0;
	
        double os_mean_nooutliers = 0; 
        double os_rms_nooutliers = 0;
	myBench->Start("meanRMSNoOutliers");
        if (!fake_image)
	  {
	    MaxCamImageTools::meanRMSNoOutliers(tempimg, mean_nooutliers, rms_nooutliers);          
	    if (conf->hasOverscan())
	      {
		MaxCamImageTools::meanRMSNoOutliers(os_tempimg, os_mean_nooutliers, os_rms_nooutliers);          
	      }
	  }
	myBench->Show("meanRMSNoOutliers");

        if (dbg_level & dbg_cluster)
	  {
	    std::cout << "     event: " << i << std::endl;
	    std::cout << "     image raw mean: " << mean << std::endl;
	    std::cout << "     image raw rms: " << rms << std::endl;
	    std::cout << "     image n.o. mean: " << mean_nooutliers << std::endl;
	    std::cout << "     image n.o. rms: " << rms_nooutliers << std::endl;
	    std::cout << "     bias mean: " << bias_mean[u] << std::endl;
	    std::cout << "     bias rms: " << bias_rms[u] << std::endl;
	    if (runnum > 0) 
	      std::cout << "     spark test: " << (mean/last_mean[u]) << std::endl;
	    
	  }
        image_mean[u] = mean_nooutliers;
	image_rms[u] = rms_nooutliers;
        os_mean[u] = os_mean_nooutliers;
        os_rms[u] = os_rms_nooutliers;
        // check for sparks
	
        double os_spark_mean = conf->hasOverscan() ? MaxCamImageTools::getMean(os_tempimg,os_cut): 0 ;
	
	if( fake_image || 
	    conf->isMC() ||
	    ( mean/last_mean[u]<conf->getSparkCut(const_cast<char*>(camera_ids[u].c_str())) && !conf->hasOverscan() ) || 
	    ( (mean/last_mean[u])<conf->getSparkCut(const_cast<char*>(camera_ids[u].c_str())) && 
	      ( conf->hasOverscan() && mean-os_spark_mean< conf->getOverscanSparkCut(const_cast<char*>(camera_ids[u].c_str())) ) ) )
	  {
	    spark[u]=0;
	    last_mean[u]=mean;
	  } else {
          spark[u]=1;
          //std::cout << "overscan spark cut for " << ser << ": " << conf->getOverscanSparkCut(const_cast<char*>(ser.Data())) << std::endl; 
	  
          //Populate sparkref ----------  Taken out for live running; wlk
          TH2F* sparkrefimg = (TH2F*) tempimg->Clone("sparkrefimg"); 
          sparkrefimg->Add(secondarybias[u],-1); 
	  
          int nsat = 0;  
	  myBench->Start("Loading sparkrefimg");
          for (int m = 1; m <= sparkrefimg->GetNbinsX(); m++)
	    {
	      for (int n = 1; n <= sparkrefimg->GetNbinsY(); n++)
		{
		  if (sparkrefimg->GetBinContent(m,n) > conf->getSatThresh())
		    {
		      nsat++; 
		    }
		}
	    }
	  myBench->Show("Loading sparkrefimg");
	  
	  myBench->Start("Thresholding");
          if (nsat > conf->getNSatThresh())
	    {
	      for (int m = 1; m <= sparkrefimg->GetNbinsX(); m++)
		{
		  for (int n = 1; n <= sparkrefimg->GetNbinsY(); n++)
		    {
		      if (sparkrefimg->GetBinContent(m,n) > conf->getSatThresh())
			{
			  if(dbg_level & dbg_error) std::cout << "Adding sparkref pixels " << m << "," << n << std::endl; 
			  sparkref_running[u].push_back(pair<int,int>(m,n)); 
			}
		    }
		}
	    }
	  myBench->Show("Thresholding");
	  
          sparkrefimg->Delete(); 
	  
        }
        
        if (dbg_level & dbg_cluster) 
	  {
	    if(spark[u] == 1) std::cout << "     SPARK!" << std::endl;
	  }
	
        // clean up image before cluster finding, same way as bias frame
        if (conf->isMC()) { // MC
          if (!conf->isNoKill())
            pixels_killed[u] = MaxCamImageTools::killLonePixels2(tempimg, mean + 10*rms);
        } else if (fake_image)
	  {
            pixels_killed[u] = 0; 
	  } else {
          pixels_killed[u] = MaxCamImageTools::killLonePixels2(tempimg,conf->getOutlierFactor()*mean);
          if (conf->hasOverscan())
	    {
	      MaxCamImageTools::killLonePixels2(os_tempimg,conf->getOutlierFactor()*osmean);
	    }
        }
	
	// we use this later to correct how much images have drifted in overall mean 
	// since the bias frame was taken.  It's important that we use the cleaned image
	// for this mean because we might otherwise pick up a hot pixel
	//double cleaned_osmean=conf->hasOverscan() ? MaxCamImageTools::getMean(os_tempimg) : 0;;
	
        if (!fake_image)
	  {
	    tempimg->Add(secondarybias[u],-1);
	    if (conf->hasOverscan())
	      {
		os_tempimg->Add(os_secondarybias[u],-1);
	      }
	  }
        //      subtract off remaining pedestal
        //
	//        double perpx = fake_image ? 0 : tempimg->Integral()/(nbinsx*nbinsy);
	
	
        if(dbg_level & dbg_error) cout << i << "," << u << ": " << bias_mean[u] << "," << mean << ";" 
		  << os_bias_mean[u] << "," << osmean << endl;
	
	// this is the correction meant to account for a possible drift from the mean since the bias frame was taken
	// it's extremely important to use the hot pixel *cleaned* overscan means for this drift correction, 
	//double perpx = fake_image ? 0 : conf->hasOverscan() ? cleaned_osmean-cleaned_os_bias_mean[u] : tempimg->Integral()/(nbinsx*nbinsy) ;
	
	//cout << perpx << endl;
	//if (dbg_level & dbg_cluster) cout << "     image pedestal correction: " << perpx << endl;

	//look for partials here
	myBench->Start("precluster");

	if(spark[u]==0 && !fake_image)
	  {
	    //MaxCamImageTools::subtractPedestal(tempimg,perpx);
	    
	    //look for partials here
	    TH2C* edges = MaxCamImageTools::edgeDetect(tempimg,
						       conf->getPartialBlur(),
						       conf->getPartialLowThresh(),
						       conf->getPartialHighThresh());
	    
	    TH1D* edges_p = edges->ProjectionX("edges_p",2,edges->GetNbinsX()-1);
	    
	    const int validbins=tempimg->GetNbinsX();
	    double edgesval[validbins];
	    int edgesind[validbins];
	    for(int j=2; j<edges_p->GetNbinsX(); j++)
	      {
		edgesval[j-1]=edges_p->GetBinContent(j);
	      }
	    edgesval[0]=0; edgesval[255]=0;
	    TMath::Sort(validbins,edgesval,edgesind);
	    
	    if(edgesval[edgesind[0]]>conf->getPartialNPrimaryThresh() )
	      {
		int p=1;
		bool secedges=false;
		while(edgesval[edgesind[p]] > conf->getPartialNSecondaryThresh())
		  {
		    if(abs(edgesind[p]-edgesind[0])>=conf->getPartialDistanceLow()&& 
		       abs(edgesind[p]-edgesind[0])<=conf->getPartialDistanceHigh())
		      {
			secedges=true;
		      }
		    p++;
		  }
		spark[u]=!secedges;
	      }
	    
	    delete edges; delete edges_p;
	  }
	
	integral[u]=rms; // this variable is written to the output tree
	                 // useful to use rms here since the spark cut
	// now cuts on rms (instead of integral) 
	
	
        // save spark images            
        if(spark[u]==1) 
	  {
	    sparks++;
	    ntracks[u] = 0; 
	    lastspark[u]=0;
	    TH2F * sparkimg = (TH2F*) tempimg->Clone("sparkimg");
	    if(eventVersion==1)
	      {
		clusti= new MaxCamClusterImage(sparkimg,d.event()->timeStamp());
	      }
	    else
	      {
		clusti= new MaxCamClusterImage(sparkimg,d.event()->UTCtimeStamp());
	      }
	    
	    // Find oversaturated pixels          
	    
	    
	  }  

	myBench->Show("precluster");

	myBench->Start("Cluster finding");
        if (fake_image)
	  {
	    clusti = new MaxCamClusterImage((TH2F*)tempimg->Clone("baseimg"),d.event()->UTCtimeStamp());  
	    ntracks[u] = 0; 
	  }
        // find clusters
        else if(spark[u]==0)
          {
            if(i>0)
              lastspark[u]++;
            if(dbg_level & dbg_error) cout << u << "," << lastspark[u] << endl;
	    
            // look for clusters
	    
            TH2F* clustimg = (TH2F*)tempimg->Clone("clustimg");
	    
            //old cluster finding algorithm 
            if (!strcmp(conf->getClusterFindingAlgorithm(algo_index),"ci")) 
	      {
		TH2F* baseimage = (TH2F*)tempimg->Clone("baseimage");
		if(eventVersion==1)
		  {
		    clusti = new MaxCamClusterImage(baseimage,d.event()->timeStamp());
		  }
		else
		  {
		    clusti = new MaxCamClusterImage(baseimage,d.event()->UTCtimeStamp());
		  }
		// rebin
		baseimage->Rebin2D(2,2);
		// blur
		baseimage = (TH2F*)MaxCamImageTools::blur(baseimage,1,conf->getBlurAmount());
		
		if (conf->useGainMap() && conf->getClusterFindCIUseGainMapToMerge())
		  {
		    ntracks[u]=  MaxCamImageTools::findClustersGM(baseimage,
								  clusti,
								  conf->getClusterMinSigma(),
								  conf->getClusterMaxSigma(),
								  conf->getClusterMinSize(),
								  conf->getClusterMinDist(),
								  conf->getClusterJoinMinRxyGlobal(),
								  conf->getClusterJoinMinRxyCluster(),
								  conf->getClusterJoinMaxJoinResidual(),
								  conf->getClusterJoinLeastSquaresWeight(),
								  gainmaps[u],
								  conf->getClusterJoinSpacerWidth()
								  );
		  }
		else
		  {
		    ntracks[u]=  MaxCamImageTools::findClustersCI(baseimage,
								  clusti,
                                                          conf->getClusterMinSigma(),
								  conf->getClusterMaxSigma(),
								  conf->getClusterMinSize(),
								  conf->getClusterMinDist());
		  }
		
		if (dbg_level & dbg_cluster) cout << "     unsmeared cluster pixel threshold: " << conf->getClusterMinSigma()*rms_nooutliers << endl;
		clusti->changeImageWithThreshold(clustimg,conf->getClusterReducedThreshold()*rms_nooutliers);  

		// now we must delete baseimage
		gROOT->Delete("baseimage"); 
		
            }
            
            //Otherwise use a new one
            else
	      {
		clusti = new MaxCamClusterImage(clustimg,d.event()->UTCtimeStamp());
		if(!strcmp(conf->getClusterFindingAlgorithm(algo_index),"seed"))
		  {
		    ntracks[u] = MaxCamImageTools::findClustersGMSeed(clustimg, 
								      clusti,
								      conf->getClusterSeedThreshold(), 
								      conf->getClusterThreshRingRatio(),
								      conf->getClusterMaxWrongProb(),
								      conf->getClusterMinThreshold(),
								      conf->getClusterBlurRadius(),
								      conf->getGaussianBlurAmount(),
								      conf->getClusterNeighborsThresholdForFilling(),
								      conf->getClusterMinNeighborsToKeepPixel(),
								      conf->getClusterMinSizeUnbinned(),
								      conf->getClusterMinDist(),
								      conf->getClusterJoinMinRxyCluster(),
								      conf->getClusterJoinMinRxyGlobal(),
								      conf->getClusterJoinMaxJoinResidual(),
								      conf->getClusterJoinLeastSquaresWeight(),
								      gainmaps[u],
								      conf->getClusterJoinSpacerWidth()); 
		    
		  }
		else if (!strcmp(conf->getClusterFindingAlgorithm(algo_index),"ad"))
		  {
		    ntracks[u] = MaxCamImageTools::findClustersADHysteresisGM(clustimg, clusti, 
									      conf->getClusterFindADK(), 
									      conf->getClusterFindADLambda(), 
									      MaxCamImageTools::TUKEY, 
									      conf->getClusterFindADNIter(), 
									      conf->getClusterFindADGradientBlurAmount(),
									      MaxCamImageTools::SOBEL, 
									      conf->getClusterFindADHighThresh(),
									      conf->getClusterFindADLowThresh(),
									      conf->getClusterNeighborsThresholdForFilling(),
									      conf->getClusterMinSizeUnbinned(),
									      conf->getClusterMinDist(),
									      conf->getClusterJoinMinRxyCluster(),
									      conf->getClusterJoinMinRxyGlobal(),
									      conf->getClusterJoinMaxJoinResidual(),
									      conf->getClusterJoinLeastSquaresWeight(),
									      gainmaps[u],
									      conf->getClusterJoinSpacerWidth());
		  }
		else
		  {
		    ntracks[u] = 0; 
		    std::cout << "NO CLUSTER FINDING ALGORITHM!!!!" << std::endl; 
		  }
		
		clusti->applyRedThreshold(conf->getClusterReducedThreshold() * rms_nooutliers); 
	      }
	    myBench->Show("Cluster finding");
	    
            sumtracks+=ntracks[u];
	    
            cout << "Ntracks: " << ntracks[u] << endl;
            //Loop over the tracks here

	    myBench->Start("Loop over tracks");
            for(int v=0; v<ntracks[u]; v++)
              {
                if (v>=15) break; 
                vector<BurninEncoded_t> this_track_burnin;
                
                vector<int> pxs_smear = clusti->getCluster(v);
                if (dbg_level & dbg_cluster) std::cout << "     smeared cluster size: " << int(pxs_smear.size()) << std::endl;
                vector<int> pxs_nosmear = clusti->getClusterRed(v);
                if (dbg_level & dbg_cluster) std::cout << "     unsmeared cluster size: " <<  int(pxs_nosmear.size()) << std::endl;
		
                if(!conf->useGainMap())
		  {
		    E[u][v] = clusti->getIntegral(v);
		    EGainMap[u][v] = 0;
		  }
                else
		  {
		    E[u][v] = clusti->getIntegral(v);
		    EGainMap[u][v] = clusti->getIntegralWithGainMap(v,gainmaps[u]);
		  }
		edge[u][v] = clusti->hitsEdge(v);
                skewness[u][v] = clusti->getSkewness(v,phi[u][v]);
                theta[u][v] = 0; //TODO: When we figure out how to get theta, fix this
                switch(conf->getPhiAlgorithm())
		  {
                  case 1: phi[u][v] = clusti->getPhi(v); break;
                  case 2: phi[u][v] = clusti->getPhi2(v); break;
                  case 3: phi[u][v] = clusti->getPhi3(v); break; 
                  case 4: phi[u][v] = clusti->getPhi4(v); break;
                  default: phi[u][v] = 0; 
		    if (dbg_level & dbg_warning) std::cout << "Warning: invalid phi algorithm selected" << std::endl; 
		  }
                double xb,yb,xe,ye;
                switch(conf->getRangeAlgorithm())
		  {                  
                  case 1: range[u][v] = clusti->getLength(v,xb,yb,xe,ye); break;
                  case 2: range[u][v] = clusti->getLength2(v,phi[u][v],1); break;
                  default: range[u][v] = 0;
		    if (dbg_level & dbg_warning) std::cout << "Warning: invalid range algorithm selected" << std::endl; 
		  }
		diffusedrange[u][v]=clusti->getDiffusedLength(v,xb,yb,xe,ye);
		clusti->getEllipseAxes(v,majoraxis[u][v],minoraxis[u][v]);
		
                clusti->getXY(v,x[u][v],y[u][v]);
                //Get cluster mean, rms, maxpixel
                
                cluster_mean[u][v] = clusti->getMean(v);
                cluster_rms[u][v] = clusti->getRMS(v,cluster_mean[u][v]);
                energy_density[u][v] = clusti->getEnergyDensity(v); 
                npixel[u][v] = clusti->getCluster(v).size(); 
                npixel_red[u][v] = clusti->getClusterRed(v).size(); 
                int maxBin; 
                maxpixel[u][v] = clusti->getMax(v,&maxBin);
                neighbors[u][v] = clusti->getNumNeighbors(v,conf->getClusterMinSigma(),maxBin);
                cygnus_angle[u][v] = clusti->getCygnusAngle(v, conf->getNorthAngle(),
                                                            orientation,
                                                            conf->getLatitude(), conf->getLongitude(),
                                                            phi[u][v],
                                                            theta[u][v]+ TMath::Pi()/2. //TODO: when we figure out theta, remove this argument
                                                            );
		
		
                clusti->getRADec( phi[u][v], theta[u][v] + TMath::Pi()/2., //TODO: ditto as above
				  d.event()->timeStamp(), 
				  conf->getLatitude(), conf->getLongitude(), 
				  conf->getNorthAngle(), orientation, 
				  ra[u][v],dec[u][v],glat[u][v],glon[u][v]);
		
                for (int m = 0; m <4; m++)
		  { 
		    moments[u][m][v] = clusti->getMoment(v,m+1,phi[u][v], 4,"pixelPerBin"); 
		    transverse_moments[u][m][v] = clusti->getMoment(v,m+1,phi[u][v]+TMath::Pi()/2., 4,"pixelPerBin"); 
		  }
		
                //Loop over last up to N_LOOPBACK events in search of things at the same position
                int pos_i = 0; //position index (since we're using iterators)

		
                list<vector<double *>**>::iterator pos_iter; //iterator for previous positions
                list<vector<vector<vector<BurninEncoded_t> > >* >::iterator burnin_iter; 
                
                for (pos_iter = positions.begin(), burnin_iter = burnin_temp.begin(); 
                     pos_iter != positions.end() ; 
                     pos_iter++, burnin_iter++)
		  {
		    vector<double *> ** test_positions = *pos_iter; 
		    vector<vector<vector<BurninEncoded_t> > > * test_burnin = *burnin_iter;
                    
                    
		    //Loop over the tracks in these events
		    for(unsigned int z = 0; z<test_positions[u]->size(); z++)
		      {
			double * position = (double *)(test_positions[u]->at(z));
			if (position[0]==0.0 && position[1] == 0.0) continue; //ignore events with no tracks
			double x_delta = x[u][v] - position[0];
			double y_delta = y[u][v] - position[1];

                         if(burnin_test(x_delta, y_delta, conf))
                         {
                           int test_event_index = position_offset + pos_i; 
                           this_track_burnin.push_back(encodeBurnin(z,test_event_index));

                           //Mark this event on the other events record...
                           int our_index = position_offset + positions.size();
                           assert(test_event_index != our_index);

                           test_burnin->at(u)[z].push_back(encodeBurnin(v,our_index));
                         }

		      }    
		    pos_i++; 
		  }//End Loop over last up to N_LOOKBACK events
		
                this_camera_burnin.push_back(this_track_burnin); 
                
              }//End Loop over tracks
	    myBench->Show("Loop over tracks");
	    
            /* initialize the rest of the event array */
            for(int v=ntracks[u]; v<15; v++)
	      {
                theta[u][v] = 0;
                phi[u][v]=0;
                E[u][v]=0;
                range[u][v]=0;
                x[u][v]=0;
                y[u][v]=0;
                edge[u][v]=0;
                skewness[u][v] = 0;
                cluster_rms[u][v]=0;
                cluster_mean[u][v]=0;
                maxpixel[u][v]=0;
                cygnus_angle[u][v]=0;
                neighbors[u][v]=-1;
                energy_density[u][v]=0;
                npixel[u][v]=0;
                dec[u][v] = 0; 
                ra[u][v] = 0; 
                glat[u][v] = 0;  
                glon[u][v] = 0; 
                npixel_red[u][v] = 0;
                for (int m = 0; m< 4; m++) 
		  {
		    moments[u][m][v] = 0; 
		    transverse_moments[u][m][v] = 0; 
		  }
		
	      }
	    
	  } // End if spark[u] == 0
        
        if (algo_index > 0)
	  {
	    const_cast<TH2F*>((TH2F*)clusti->getImage())->Delete(); 
	    clusti->forgetImage(); 
	  }
        clust->Add(clusti);
        this_event_burnin->push_back(this_camera_burnin);
	
	
	
      } //End loop over cameras
    
    myBench->Show("Loop over cameras");
    
    /* Copy sparkref to list */
    vector<pair<int,int> > * this_sparkref = new vector<pair<int,int> >[ncamera]; 
    for (int u = 0; u < ncamera; u++)
      {
	for (unsigned int v = 0; v < sparkref_running[u].size(); v++)
	  {
	    this_sparkref[u].push_back(sparkref_running[u][v]); 
	  }
	
	assert (this_sparkref[u].size()==sparkref_running[u].size()); 
      } 
    sparkref.push_back(this_sparkref); 
    

    //Clear temporary images
    if (!fake_image)
      {
	gROOT->Delete("tempimg"); 
	gROOT->Delete("copy");
      }
    
    // fill tree
    myBench->Start("Filling tree");
    
    if(sumtracks == 0 && sparks==0) nempty++;
    
    // if(sumtracks>0 || nempty%100 ==0 || sparks>0)
    // Now just write out every track... 
    { 
      for(unsigned int csn=0; csn<cameraSerialNumber->size(); csn++)
        {
          //std::cout << (*cameraSerialNumber)[csn] << std::endl;
        }
      eventnum=i;
      skimevent = d.event();
      tmpoutfile->cd();
      tmpskimtree->Fill();
      gROOT->cd();
      cameraSerialNumber->clear();
      
      trigger_groups->SetOwner(kTRUE); 
      trigger_groups->Clear(); 
      clust->SetOwner(kTRUE); 
      clust->Clear(); 
      //Add self to position list, removing first entry if list is already N_LOOKBACK
      vector<double *> ** these_positions = new vector<double *>*[ncamera];
      for (int u = 0; u < ncamera; u++)
        {
          these_positions[u] = new vector<double *>; 
          for (int v = 0; v < ntracks[u]; v++)
	    {
	      if (v >=15) continue; 
	      double * this_pos = new double[2];
	      this_pos[0] = x[u][v];
	      this_pos[1] = y[u][v];
	      these_positions[u]->push_back(this_pos); 
	    }
        } 
      
      myBench->Show("Filling tree");
      
      positions.push_back(these_positions); 
      
      //Store the burnin_info in the buffer until we go past it enough times
      burnin_temp.push_back(this_event_burnin); 
      
      
      //After we start filling the buffer, 
      // it is time to clear the first entry
      // and to write out the first event
      
      if (positions.size() >(unsigned int) conf->getBurninNumEvents() && !fake_image)
	{
	  
	  //Delete the stored positions. 
	  for (int u = 0; u < ncamera; u++)                 
	    {
	      for (unsigned int z = 0; z < positions.front()[u]->size(); z++)
		{
		  delete positions.front()[u]->at(z);
		}
	      delete positions.front()[u];
	    }
	  delete positions.front();
	  positions.pop_front(); 
	  
	  
	  burnin.push_back(burnin_temp.front());
	  burnin_temp.pop_front(); 
	  
	  position_offset++;
	  gROOT->cd(); 
	}
      myBench->Show("Event");
      myBench->Reset();
    }
  }//End loop over events
  
  
  if(dbg_level & dbg_error) std::cout << "End of loop over events" << std::endl;
  
  //Empty out remaining event buffer. 
  while (burnin_temp.size() > 0)
    {
      burnin.push_back(burnin_temp.front());
      burnin_temp.pop_front(); 
    }
  
  
  //Write out and clean up
  tmpoutfile->cd();
  tmpskimtree->Write();
  
  if (algo_index == 0)
    {
      biasoutfile->cd();
      biastree->Write();
    }
  // std::cout << sd.getGainMaps()->GetEntries() << std::endl;
  //  std::cout << sd.getGainMap(0)->GetName() << std::endl;
  
  
  //Save configuration
  if (algo_index == 0)
    {
      std::stringstream str; 
      conf->print(str); 
      sd.setConfig(str.str().c_str()); 
      sd.writeConfig();
    }
  
  sd.mergeTempTrees(tmpskimtree,&burnin, &sparkref, runnum);

  if(dbg_level & dbg_error)  cout << "Merge complete " << endl;

  if (conf->useGainMap())
    {
      for (int i = 0; i < ncamera; i++) 
	{
	  if (gainmaps[i]) gainmaps[i]->Delete(); 
	}
    }
  
  
  
  while(burnin.size() > 0)
    {
      delete burnin.front(); 
      burnin.pop_front(); 
    }
  
  /*while(sparkref.size() > 0)
    {
      delete sparkref.front(); 
      sparkref.pop_front(); 
      }*/
  
  delete tmpskimtree;

  if(dbg_level & dbg_error) cout << "Closing file: " << tmpoutfilename << endl;
  tmpoutfile->Close();
  gSystem->Sleep(500);
  if(dbg_level & dbg_error) cout << "File closed" << endl;

  
  delete biastree;
  if (algo_index == 0) 
    {
      biasoutfile->Close();
    }
  
  if (delete_intermediate)
    {
      if (dbg_level & dbg_progress) std::cout << "Deleting " << tmpoutfilename << std::endl;
      unlink(tmpoutfilename.Data());
    }
/*TString rmcomm("rm ");
  rmcomm+=biasoutfilename;
  gSystem->Exec(rmcomm);*/
  return 0;
}

static const int nreqmod = 0; 

/* Main function */ 

int main(int argn, char ** argv)
{
   int return_val = 0; 
   vector<Double_t> phi;
   vector<Double_t> radii;
   TGraphPolar* phiHistogram;
   TCanvas* canvas;
   //initGUI(canvas);

  /* Argument Parsing */ 
   TString rawdatafiles = "files.txt";
   TString keyfilename = "keys.txt";
   TString outdir = "/temp/"; 
   TString indir = "/short/";
   TString finaldir = " /skim/";
   char * config = 0;
   if(parse_args(argn, argv, &rawdatafiles, &keyfilename, &outdir,&config)) return -1; 
   //if (dbg_level & dbg_config) config->print(); 
   
   
  /* Handle keys */
   TString key = "skim"; 
   TString reqmod[nreqmod];
   bool pass; 
   //DmtpcKeys k(keyfilename, rawdatafiles, key, outdir, nreqmod, reqmod, pass); 
      //if (!pass) return -1; 
    
   

   /* Loop through input files */

   //std::cout << "k.getNFiles " << k.getNFiles() << std::endl;


   const char* file1;
   Int_t loop_counter = 0;
     while(true){
       cout << "Cycling... " << loop_counter << endl;
     loop_counter+=1;
     if(loop_counter>25) return -1;  ////  Quick fix to overcome the trouble with cycling...
     void *dir1 = gSystem->OpenDirectory(indir);
     while(file1 = gSystem->GetDirEntry(dir1)){
       if(strcmp(file1,".")!=0 && strcmp(file1,"..")!=0){
	 TString filename1(file1);
	 TString fullname(indir);
	 fullname += file1;
	 cout << "Skimming file:  " << fullname << endl;
       
      ///Need to add snippit to get base tree without using DmtpcKeys
      //May speed up process by getting bastree once for entire while loop
      DmtpcDataset d; 

      int i=Convert(filename1,0,0,d);

      //cout << "Root file: "  << fullname << " opened." << endl;

      //Use sample file to get the raw data tree....
      TTree* basetree=0;
      TFile* file2 = new TFile("/logs/sample.root");
      basetree = (TTree*)file2->Get("dmtpc");

      
      //Extract run number from file       
      d.getEvent(0); 
      int runnum = d.event()->runNumber();
      if(dbg_level & dbg_error) cout << "Run Number:  " << runnum << endl; 
      
      //Add friends
      //basetree->AddFriend("skim",fullname); 
      //cout << "Friend added??" << endl;
      
      key =  "skim"; 
      TString outfilename = filename1.ReplaceAll(".root","skim.root"); 
      DmtpcSkimDataset sd; 
      
      sd.newRootFile(outdir + outfilename); 
      if(dbg_level & dbg_error) cout << "Second root file oppened: "  << outdir + outfilename << endl;

      /* Now ready to process this file */
      CleanSkimConfig * conf = skim::preprocess(&d,config); 
      
      for (unsigned int nprocess = 0; nprocess < conf->getNClusterFindingAlgorithms(); nprocess++)
      {
        if(nprocess > 0) key = TString(conf->getClusterFindingAlgorithm(nprocess)); 

        return_val += cleanSkim(d,key,basetree,runnum,sd,true,conf,nprocess); 
      }
      /************* Move clean file for GUI, remove temp files *****************/
      TString mvcom = "mv ";
      TString space = " ";
      mvcom += outdir + outfilename + space + finaldir + outfilename;
      if(dbg_level & dbg_error) cout << "Moving final file: " << mvcom << endl;
      gSystem->Exec(mvcom);
      TString rmcom = "rm ";
      rmcom += outdir + outfilename.ReplaceAll("skim.root","bias.root");
      if(dbg_level & dbg_error) cout << rmcom << endl;
      gSystem->Exec(rmcom);
      TString rmcom_two = rmcom.ReplaceAll("bias.root",".root");
      if(dbg_level & dbg_error) cout << rmcom_two << endl;
      gSystem->Exec(rmcom_two);
       }
     }
     }
     gSystem->Sleep(500);
   return return_val;
}

int Convert(TString filename, int runno, int eventno, DmtpcDataset & a)
{
  //
  // Input file
  if(dbg_level & dbg_error) std::cout << "filename " << filename << std::endl;
  TString infile="/short/";
  infile+=filename;
  TFile inFile(infile);
  if(dbg_level & dbg_error) cout << "After TFile" << endl;
  while(inFile.IsZombie()){
    ///Return to top of function???
    if(dbg_level & dbg_error) cout << "Error opening file" << endl;
    return -1;
  }

  if(dbg_level & dbg_error) inFile.ls();
  //
  // Load 
  TH2S* ccd=(TH2S*) inFile.Get("ccd_0;1");
  ccd->SetDirectory(0);
  TH2F* biasFrame1=(TH2F*) inFile.Get("biasFrame1;1");
  biasFrame1->SetDirectory(0);
  TFile out("junk.root","recreate");
  ccd->Write();
  biasFrame1->Write();
  inFile.Close();
  //out.ls();
  //
  // Setup output file
  TString outfile="/temp/";
  outfile+=filename;

  if(dbg_level & dbg_error) cout << "Outfile name: " << outfile << endl;
  //outfile+=".root";
  //
  // Setup dataset
  a.createRootFile(outfile,"recreate");
  if(dbg_level & dbg_error) cout << "File created" << endl;
  a.setComment("Single event");
  a.setLocation("44-013");
  a.comment()->Write();
  a.keyword()->Write();
  a.location()->Write();
  if(dbg_level & dbg_error)   cout << "Comments written" << endl;
  //
  // Setup configuration
  MaxCamConfig* ccdconfig=new MaxCamConfig("ccdconfig","Andor camera configuration");
  ccdconfig->cameraID=0;
  ccdconfig->row_width=1024;
  ccdconfig->img_rows=1024;
  ccdconfig->hbin=1024;
  ccdconfig->vbin=1024;
  ccdconfig->ul_x=0;
  ccdconfig->ul_y=0;
  ccdconfig->lr_x=1024;
  ccdconfig->lr_y=1024;
  ccdconfig->CCDTemp=-20;
  ccdconfig->CCDTempSet=-20;
  ccdconfig->exposureTime=1000;
  ccdconfig->frameType=0;
  ccdconfig->nFlushes=-1;
  ccdconfig->bitDepth=65535;
  ccdconfig->daqTime=-1;
  ccdconfig->serialNumber="0";

  if(dbg_level & dbg_error) cout << "ccdconfig stuff done" << endl;
  //
  // Write bias frame
  biasFrame1->Write("biasFrame1");
  if(dbg_level & dbg_error) cout << "Biasframe 1 written" << endl;
  //
  // Write event
  DmtpcEvent* ev=a.event();
  ev->setEventNumber(eventno);
  ev->setRunNumber(runno);
  new ((*a.event()->rawCcdData())[0]) TH2S(*ccd);
  new((*a.event()->ccdConfig())[0]) MaxCamConfig(*ccdconfig);
  a.fill();
  a.write();
      TString rmcom = "rm ";
      rmcom += infile;
      if(dbg_level & dbg_error) cout << rmcom << endl;
      gSystem->Exec(rmcom);
}

static bool burnin_test( double x_delta, double y_delta, CleanSkimConfig * conf)
{
  double xd = TMath::Abs(x_delta);
  double yd = TMath::Abs(y_delta); 

  double th = conf->getBurninDistanceThresh(); 
  burnin_method_t method= conf->getBurninMethod();

  switch(method)
  {
    case SQUARE:
      return (xd <= th && yd <= th);
    case CIRCLE:
      return (xd*xd <= th*th && yd*yd <= th*th); 
    case CROSS:
      return (xd <= th || yd <= th);
  }
  
  //if we got here... something's wrong
  if (dbg_level & dbg_warning)
  {
    cerr << "Invalid Burnin Test: " << method << endl; 
  }
  return false; 
}


static int parse_args( int nargs, char ** args, TString * files, TString * keys, TString * out, char** cfg)
{
  int c = 0; 
  for (int i = 1; i < nargs; i++)
  {
    /* Handle switches */
    if (args[i][0]=='-')
    {
      //Check for debug level
      if(strcmp(args[i],"-d")==0) 
      {
        char * endptr; 
        dbg_level = (u_int32_t) strtol(args[++i], &endptr, 10);
        if (strlen(endptr) > 0) // This means a digit wasn't passed
        {
           if (dbg_level & dbg_error) cerr << "Non number was passed to debug level " << std::endl; 
           return 1; 
        }
      }

      if(strcmp(args[i],"-c")==0) 
      {
        //free(*cfg); 
        *cfg = args[++i]; 
      }
                  
    }
    else
    {
      switch(c++)
      {
        case 0: *files = TString(args[i]); break; 
        case 1: *keys = TString(args[i]); break;
        case 2: *out = TString(args[i]); break; 
        default: 
          if (dbg_level & dbg_error) cerr << "Unexpected Extra Argument. Aborting!" << std::endl;
          return 1; 
      }
    }
  }
  
  return 0; 
}

void initGUI(TCanvas* canvas)
{
	// Create a 3 by 3 grid with default margins
        canvas = new TCanvas("canvas", "Neutron detector GUI");
	canvas->Divide(3,3);
	
	// Enter time pad (top-left)
	canvas->cd(1);
	
	TTimeStamp *startTime = new TTimeStamp();
	TString startTimeText = startTime->AsString("s");
	
	TText *startTimeLabel = new TText(0.0, 0.9, startTimeText);
	startTimeLabel->SetTextAlign(13); // left-top adjusted
	startTimeLabel->SetTextSize(0.10);
	startTimeLabel->Draw();
	
	/*currentTimeLabel = new TText(0.0, 1.0, startTimeText);
	currentTimeLabel.SetTextAlign(13); // left-top adjusted
	currentTimeLabel.SetTextSize(0.10);
	currentTimeLabel.Draw();
	
	// Enter summary pad (top-right)
	canvas->cd(3);
	numberDetected = new TText(1.0, 1.0, "Tracks detected: "); 
	numberDetected.SetTextAlign(33); // right-top adjusted
	numberDetected.SetTextSize(0.10);
	numberDetected.Draw();
	
	averagePhi = new TText(1.0, 0.9, "Average direction: ");
	averagePhi.SetTextAlign(33); // right-top adjusted
	averagePhi.SetTextSize(0.10);
	averagePhi.Draw();
	
	rms = new TText(1.0, 0.8, "RMS:");
	rms.SetTextAlign(33); // right-top adjusted
	rms.SetTextSize(0.10);
	rms.Draw();*/

		
}

void refreshPhiHistogram(vector<Double_t> phi, vector<Double_t> radii, TGraphPolar* phiHistogram, TCanvas* canvas)
{
	delete phiHistogram;

	canvas->cd(7);	

	phiHistogram = new TGraphPolar();
	
	for(Int_t i = 0; i < phi.size(); i++)
	{
	  phiHistogram->SetPoint(i,phi.at(i),radii.at(i));
	}
	
	phiHistogram->Draw();
}
