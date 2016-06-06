#include "../../../MaxCam/MaxCamMC.hh"
#include "../../../MaxCam/MaxCamImageTools.hh"
#include "../../../MaxCam/MaxCamCluster.hh"
#include "../../../MaxCam/DmtpcEvent.hh"
#include "../../../MaxCam/DmtpcDataset.hh"
#include "../../../MaxCam/DmtpcKeys.hh"
#include "../../../MaxCam/DmtpcSkimDataset.hh"
#include "../../../MaxCam/DmtpcSkimEvent.hh"
#include "../../../MaxCam/DmtpcSkimRunSummary.hh"
#include "../../../MaxCam/MaxCamTriggerGroup.hh"
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
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <list>
#include <stdlib.h>
#include <sstream>
#include <stdio.h>
#include "cleanSkimConfig.hh"


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
static void samepos_copy_next(int * burnin, int * nburnin, list<int *> * same_positions, int ncamera, int );
static int parse_args( int nargs, char ** args, TString * files, TString * keys, TString * out, CleanSkimConfig ** );
static void handle_legacy_events(int,int,DmtpcDataset &, TH2F *, double *, double *, TH2F **, int *, double *, int, CleanSkimConfig *);
static bool burnin_test( double xdelta, double ydelta, CleanSkimConfig * conf);
/******** end prototypes ******/

/******** cleanSkim *********/
int cleanSkim(DmtpcDataset & d, TString & key, TTree * rawtree, int runnum, 
              DmtpcSkimDataset & sd, bool delete_intermediate, CleanSkimConfig * conf)
{

  TString tmpoutfilename = TString(sd.getFileName()).ReplaceAll(".root","_tmp.root");

  /* Define Temporary out files and trees */
  TFile* tmpoutfile = new TFile(tmpoutfilename,"RECREATE");
  TTree* tmpskimtree;
  TTree* burnintree; 
  TH2F* tempimg;

  /* Get some eventwide information */
  d.getEvent(0);
  const int ncamera = d.event()->ccdData()->GetEntries();
  const int nbinsx =  d.event()->ccdData(0)->GetNbinsX();
  const int nbinsy =  d.event()->ccdData(0)->GetNbinsY();
  const int eventVersion = d.event()->IsA()->GetClassVersion();
  int ncam = ncamera;

  /* Support processing on datasets with no actual images */ 
  bool fake_image = false; 
  if (ncamera ==0 || d.event()->ccdData(0)==NULL || d.event()->ccdData(0)->Integral() ==0)
  {
    fake_image = true; 
    tempimg = new TH2F ("tempimg","tempimg",256,0,1024,256,0,1024); 
  }

	


  /* Create the temporary skim tree without burnin data.  */
  tmpskimtree = new TTree(key,"Skimmed Events");
  double theta[ncamera][15], phi[ncamera][15], E[ncamera][15];
  double range[ncamera][15], x[ncamera][15], y[ncamera][15];
  double skewness[ncamera][15], integral[ncamera];
  bool edge[ncamera][15], spark[ncamera];
  int ntracks[ncamera], date, time, eventnum;
  int pixels_killed[ncamera], npixel_red[ncamera][15];
  double image_mean[ncamera],image_rms[ncamera];
  MaxCamClusterImage* clusti;
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
  
  tmpskimtree->Branch("ncamera",&ncam,"ncamera/I");
  tmpskimtree->Branch("spark",&spark,"spark[ncamera]/O");
  tmpskimtree->Branch("integral",&integral,"integral[ncamera]/D");
  tmpskimtree->Branch("eventnum",&eventnum,"eventnum/I");
  tmpskimtree->Branch("theta",theta,"theta[ncamera][15]/D");
  tmpskimtree->Branch("phi",phi,"phi[ncamera][15]/D");
  tmpskimtree->Branch("E",E,"E[ncamera][15]/D");
  tmpskimtree->Branch("energy_density", energy_density, "energy_density[ncamera][15]/D");
  tmpskimtree->Branch("range",range,"range[ncamera][15]/D");
  tmpskimtree->Branch("x",x,"x[ncamera][15]/D");
  tmpskimtree->Branch("y",y,"y[ncamera][15]/D");
  tmpskimtree->Branch("ntracks",ntracks,"ntracks[ncamera]/I");
  tmpskimtree->Branch("edge",edge,"edge[ncamera][15]/O");
  tmpskimtree->Branch("skewness",skewness,"skewness[ncamera][15]/D");
  tmpskimtree->Branch("ra",ra,"ra[ncamera][15]/D");
  tmpskimtree->Branch("dec",dec,"dec[ncamera][15]/D");
  tmpskimtree->Branch("glat",glat,"glat[ncamera][15]/D");
  tmpskimtree->Branch("glon",glon,"glon[ncamera][15]/D");
  tmpskimtree->Branch("moments",moments,"moments[ncamera][15][4]/D");
  tmpskimtree->Branch("transverse_moments",transverse_moments,"transverse_moments[ncamera][15][4]/D");
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
  tmpskimtree->Branch("trigger_groups","TObjArray",&trigger_groups,128000,0);

  gROOT->cd();

  /* Open up file for bias frame saving */
  TString biasoutfilename = TString(sd.getFileName()).ReplaceAll(key+".root","bias.root");
  biasoutfilename.ReplaceAll(key+".root","bias.root");
  TFile* biasoutfile = new TFile(biasoutfilename,"RECREATE");
  TH2F* biasframe = d.event()->ccdData(0);
  TTree* biastree = new TTree("bias","Bias Information");
  biastree->Branch("biasframe","TH2F",&biasframe,128000,0);
  gROOT->cd();

  double sigmathr[ncamera];
  int nframes[ncamera];
  for(int u=0; u<ncamera; u++){nframes[u]=0;}
  TH2F* secondarybias[ncamera];
  double bias_mean[ncamera], bias_rms[ncamera];
  double last_mean[ncamera]; // for the spark cut

  // handle the bias frame
  //      define the bias threshold for runs<500 and MC
  for(int u=0; u<ncamera; u++)
    {
      if (fake_image)
      {
         break; 
      }
      if (dbg_level & dbg_bias) cout << "Camera " << u << " image cleaning:" << endl;
      secondarybias[u] = (TH2F*)d.getBiasFrame(u+1)->Clone();
      
      // the cluster finding threshold is very sensitive to this number!
      // using 1.5 sigma reproduces cluster finding efficiency from AF v1
      sigmathr[u]=1.5*MaxCamImageTools::getRMS(secondarybias[u]);
      
      
      if (dbg_level & dbg_bias)
      {  
        cout << "     camera bias MEAN: " << 
             MaxCamImageTools::getMean(secondarybias[u]) << endl;
        cout << "     camera bias RMS: " << 
             MaxCamImageTools::getRMS(secondarybias[u]) << endl;
      }
    }

  
  const int nev = rawtree->GetEntries();
  
  /* Special bias handling for events with run number < 500 */
  if ( runnum < 500 &&  strstr(sd.getFileName(),"dmtpc_run")!=NULL) {
    handle_legacy_events(nev, ncamera, d, tempimg, bias_mean, bias_rms, secondarybias, nframes, sigmathr, runnum, conf);
  } 
  
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
          if (dbg_level & dbg_bias) cout << "run number " << runnum << " using averaged bias " << endl;
          secondarybias[u]->Scale(1/double(nframes[u]));
        } else {
          secondarybias[u] = (TH2F*)d.getBiasFrame(u+1)->Clone();
          if (dbg_level & dbg_bias) cout << "run number " << runnum << " using camera bias " << endl;
        }
        // kill 5 sigma outliers in the bias frame, iterate 3 times
        for (int it=0; it<3; it++) {
          double mean, rms;
          mean=MaxCamImageTools::getMean(secondarybias[u]);
          rms=MaxCamImageTools::getRMS(secondarybias[u]);
          if (conf->isMC()) { // MC
            MaxCamImageTools::killPixels(secondarybias[u],mean+10.*rms);
          } else {
            MaxCamImageTools::killPixels(secondarybias[u],mean*conf->getOutlierFactor());
          }
          if (dbg_level & dbg_bias) 
          {
            cout << "     iteration " << it << " averaged bias image mean (no outliers): " << mean << endl;
            cout << "     iteration " << it << " averaged bias image rms (no outliers): " << rms << endl;
          }
          bias_mean[u]=mean;
          bias_rms[u]=rms;
        }
      }
      biasframe = secondarybias[u];
      biasoutfile->cd();
      biastree->Fill();
      gROOT->cd();   
    }
    
  if (dbg_level & dbg_progress) cout << "All preclean activities done" << endl;
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

  // loop over events
  for(int i = 0; i<nev; i++)
  {
     if (dbg_level & dbg_progress) cout << i << endl;
     d.getEvent(i);
     int sumtracks=0;
     int sparks=0;
     clust->Clear();
     trigger_groups->Clear(); 
     
     vector<vector<vector<BurninEncoded_t> > > * this_event_burnin= new vector<vector<vector<BurninEncoded_t> > >;

     //Save Scope Data
     int nch = conf->getNChannelsPerTrigger();  
     int nent = d.event()->scopeData()->GetEntries(); 
     int ntr = nent/nch; 
     //cout << "Ntriggers: " << ntr << endl; 
     //cout << "Nentries: " << nent << endl; 
     sumtracks+= ntr; 
     for (int g = 0; g < ntr; g++)
     {
       MaxCamTriggerGroup * grp =  new MaxCamTriggerGroup(g); 

       for (int c = g; c < nent; c+=ntr)
       {
      //   cout << "c: " << c << " g: " << g << endl;
         MaxCamWaveform * wf = new MaxCamWaveform(d.event()->scopeData(c), (unsigned int) g, 
                          conf->getChannelId(c/ntr), conf->getInvertWaveform(c/ntr),
                          conf->getWaveformNBaseline(c/ntr), conf->getWaveformNBlur(c/ntr));
      //   cout << "wf pointer: " << wf << endl;  
      //   cout << "wf->getWaveForm" << wf->getWaveform() << endl; 
      //   cout << "grp pointer: " << grp << endl;  
      //   cout << "wf rise time: " << wf->getRiseTime() << endl; 
         grp->addWaveform(wf); 
       }
       trigger_groups->Add(grp); 
     }
       
     for(int u=0; u<ncamera; u++)
     {
        vector< vector<BurninEncoded_t> >  this_camera_burnin ;
        MaxCamClusterImage::CAMERA_ORIENTATION orientation = conf->getCameraOrientation(u);         

        if (dbg_level & dbg_progress) cout << "\t" << u << ":" << endl;
        // keep track of the last frame for the spark cut
        if (i==0) {
          last_mean[u] = bias_mean[u];
        }

        if(!fake_image)
          tempimg = (TH2F*)d.event()->ccdData(u)->Clone("tempimg");

        double mean, rms;
        mean= fake_image ? 0 : MaxCamImageTools::getMean(tempimg);
        rms= fake_image ? 0 : MaxCamImageTools::getRMS(tempimg);
        double mean_nooutliers = 0; 
        double rms_nooutliers = 0;
        if (!fake_image)
          MaxCamImageTools::meanRMSNoOutliers(tempimg, mean_nooutliers, rms_nooutliers);          
        
        if (dbg_level & dbg_cluster)
        {
          cout << "     image raw mean: " << mean << endl;
          cout << "     image raw rms: " << rms << endl;
          cout << "     image n.o. mean: " << mean_nooutliers << endl;
          cout << "     image n.o. rms: " << rms_nooutliers << endl;
          cout << "     bias mean: " << bias_mean[u] << endl;
          cout << "     bias rms: " << bias_rms[u] << endl;
          if (runnum > 0) 
            cout << "     spark test: " << (mean/last_mean[u]) << endl;

        }
        image_mean[u] = mean_nooutliers;
        image_rms[u] = rms_nooutliers;
        // check for sparks
        if(fake_image || (mean/last_mean[u])<1.01) {
          spark[u]=0;
          last_mean[u]=mean;
        } else {
          spark[u]=1;

          //Populate sparkref 
          TH2F* sparkrefimg = (TH2F*) tempimg->Clone("sparkrefimg"); 
          sparkrefimg->Add(secondarybias[u],-1); 
         
          int nsat = 0;  
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
          if (nsat > conf->getNSatThresh())
          {
            for (int m = 1; m <= sparkrefimg->GetNbinsX(); m++)
            {
              for (int n = 1; n <= sparkrefimg->GetNbinsY(); n++)
              {
                if (sparkrefimg->GetBinContent(m,n) > conf->getSatThresh())
                {
                  cout << "Adding sparkref pixels " << m << "," << n << endl; 
                  sparkref_running[u].push_back(pair<int,int>(m,n)); 
                }
              }
            }
          }

          sparkrefimg->Delete(); 
        }
        if (conf->isMC()) spark[u] = 0; // MC, for now
        
        if (dbg_level & dbg_cluster) 
        {
          if(spark[u] == 1) cout << "     SPARK!" << endl;
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
        }

        if (!fake_image)
          tempimg->Add(secondarybias[u],-1);
        //      subtract off remaining pedestal
        //
        double perpx = fake_image ? 0 : tempimg->Integral()/(nbinsx*nbinsy);
        if (dbg_level & dbg_cluster) cout << "     image pedestal correction: " << perpx << endl;
        if(spark[u]==0 && !fake_image)
           MaxCamImageTools::subtractPedestal(tempimg,perpx);
        integral[u]=rms; // this variable is written to the output tree
                         // useful to use rms here since the spark cut
                         // now cuts on rms (instead of integral) 

        // save spark images            
        if(spark[u]==1) 
        {
           sparks++;
           ntracks[u] = 0; 
           TH2F * sparkimg = (TH2F*) tempimg->Clone("sparkimg");
           if(eventVersion==1)
           {
             clusti= new MaxCamClusterImage(sparkimg,d.event()->timeStamp());
           }
           else
           {
             clusti= new MaxCamClusterImage(sparkimg,d.event()->UTCtimeStamp());
           }
          
           /* Find oversaturated pixels */           


        }


        if (fake_image)
        {
          clusti = new MaxCamClusterImage((TH2F*)tempimg->Clone("baseimg"),d.event()->UTCtimeStamp());  
          ntracks[u] = 0; 
        }
        // find clusters
        else if(spark[u]==0)
          {
            TH2F* baseimage = (TH2F*)tempimg->Clone("baseimage");
            // rebin
            baseimage->Rebin2D(2,2);
            // blur
            baseimage = MaxCamImageTools::blur(baseimage,1,conf->getBlurAmount());
            

            // look for clusters
            if(eventVersion==1)
            {
              clusti = new MaxCamClusterImage(baseimage,d.event()->timeStamp());
            }
            else
            {
               clusti = new MaxCamClusterImage(baseimage,d.event()->UTCtimeStamp());
            }
              
            ntracks[u]=  MaxCamImageTools::findClustersCI(baseimage,
                                                        clusti,
                                                        conf->getClusterMinSigma(),
                                                        conf->getClusterMaxSigma(),
                                                        conf->getClusterMinSize(),
                                                        conf->getClusterMinDist());
            sumtracks+=ntracks[u];


            // get reconstructed quantities
            TH2F* clustimg = (TH2F*)tempimg->Clone("clustimg");

            if (dbg_level & dbg_cluster) cout << "     unsmeared cluster pixel threshold: " << conf->getClusterMinSigma()*rms_nooutliers << endl;
            clusti->changeImageWithThreshold(clustimg,conf->getClusterMinSigma()*rms_nooutliers);  
            // now we must delete baseimage
            gROOT->Delete("baseimage"); 


            //Loop over the tracks here
            for(int v=0; v<ntracks[u]; v++)
              {
                if (v>=15) break; 
                vector<BurninEncoded_t> this_track_burnin;
                
                vector<int> pxs_smear = clusti->getCluster(v);
                if (dbg_level & dbg_cluster) cout << "     smeared cluster size: " << int(pxs_smear.size()) << endl;
                vector<int> pxs_nosmear = clusti->getClusterRed(v);
                if (dbg_level & dbg_cluster) cout << "     unsmeared cluster size: " <<  int(pxs_nosmear.size()) << endl;

                E[u][v] = clusti->getIntegral(v);
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
                  if (dbg_level & dbg_warning) cout << "Warning: invalid phi algorithm selected" << endl; 
                }
                double xb,yb,xe,ye;
                switch(conf->getRangeAlgorithm())
                {                  
                  case 1: range[u][v] = clusti->getLength(v,xb,yb,xe,ye); break;
                  case 2: range[u][v] = clusti->getLength2(v,phi[u][v],1); break;
                  default: range[u][v] = 0;
                  if (dbg_level & dbg_warning) cout << "Warning: invalid range algorithm selected" << endl; 
                }

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
                       for(int z = 0; z<test_positions[u]->size(); z++)
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
        
        clust->Add(clusti);
        this_event_burnin->push_back(this_camera_burnin);

       

     } //End loop over cameras


     /* Copy sparkref to list */
     vector<pair<int,int> > * this_sparkref = new vector<pair<int,int> >[ncamera]; 
     for (int u = 0; u < ncamera; u++)
     {
       for (int v = 0; v < sparkref_running[u].size(); v++)
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
     if(sumtracks == 0 && sparks==0) nempty++;

    // if(sumtracks>0 || nempty%100 ==0 || sparks>0)
    // Now just write out every track... 
     { eventnum=i;
        skimevent = d.event();
        tmpoutfile->cd();
        tmpskimtree->Fill();
        gROOT->cd();

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

        positions.push_back(these_positions); 

        //Store the burnin_info in the buffer until we go past it enough times
        burnin_temp.push_back(this_event_burnin); 
        

        //After we start filling the buffer, 
        // it is time to clear the first entry
        // and to write out the first event
      
       if (positions.size() > conf->getBurninNumEvents() && !fake_image)
       {

         //Delete the stored positions. 
         for (int u = 0; u < ncamera; u++)                 
         {
            for (int z = 0; z < positions.front()[u]->size(); z++)
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
     }
  }//End loop over events

  //Empty out remaining event buffer. 
  while (burnin_temp.size() > 0)
  {
     burnin.push_back(burnin_temp.front());
     burnin_temp.pop_front(); 
  }


  //Write out and clean up
  tmpoutfile->cd();
  tmpskimtree->Write();
  biasoutfile->cd();
  biastree->Write();
  
  sd.mergeTempTrees(tmpskimtree,&burnin, &sparkref, runnum);    
  
  
  while(burnin.size() > 0)
  {
    delete burnin.front(); 
    burnin.pop_front(); 
  }

  while(sparkref.size() > 0)
  {
    delete sparkref.front(); 
    sparkref.pop_front(); 
  }
  
  delete tmpskimtree;
  tmpoutfile->Close();
  
  delete biastree;
  biasoutfile->Close();
  


  if (delete_intermediate)
  {
    if (dbg_level & dbg_progress) cout << "Deleting " << tmpoutfilename << endl;
    unlink(tmpoutfilename.Data());
  }

  return 0;
}

static const int nreqmod = 0; 

/* Main function */ 

int main(int argn, char ** argv)
{
   int return_val = 0; 

  /* Argument Parsing */ 
   TString rawdatafiles = "files.txt";
   TString keyfilename = "keys.txt";
   TString outdir = "./skim/"; 
   //TString sumdir = "./sum/";
   CleanSkimConfig * config = new CleanSkimConfig; 
   if(parse_args(argn, argv, &rawdatafiles, &keyfilename, &outdir,&config)) return -1; 
   TString sumdir = outdir;//summary file directory
   sumdir+="sum/";
   if (dbg_level & dbg_config) config->print(); 
   
   
  /* Handle keys */
   TString key = "skim"; 
   TString reqmod[nreqmod];
   bool pass; 
   
   DmtpcKeys k(keyfilename, rawdatafiles, key, outdir, nreqmod, reqmod, pass); 
   
   if (!pass) return -1; 
    
   /* Loop through input files */
   for (int f = 0; f < k.getNFiles(); f++)
   {
        if (dbg_level & dbg_progress) cout << k.getFile(f) << endl;
                  
        /* Create DMTPC Database and draw out tree */
        DmtpcDataset d; 
        d.openRootFile(k.getRootDirName() + k.getFile(f)); 
        
        TTree * rawtree = k.getBaseTree(f); 
        
        //Extract run number from file name 
	d.getEvent(0); 
        int runnum = d.event()->runNumber(); 
      
      //Add friends
      k.addFriends(rawtree, f); 
      
      TString outfilename = k.getFile(f).ReplaceAll(".root",key+".root"); 
      DmtpcSkimDataset sd; 
      
      sd.newRootFile(outdir + outfilename); 

      /* Now ready to process this file */
      return_val += cleanSkim(d,key,rawtree,runnum,sd,true,config); 
   }
   return return_val;
}


/************* Static Void Functions *****************/
static void handle_legacy_events(int nev, int ncamera, DmtpcDataset &d, TH2F * tempimg, double * bias_mean, 
                                 double * bias_rms, TH2F ** secondarybias, int * nframes, double * sigmathr, int runnum,
                                 CleanSkimConfig * conf
                                 )
{                                 
  for(int i=0; i<nev; i=i+100)
  {
    if (dbg_level & dbg_bias) cout << "making bias frame: event " << i << endl;
    d.getEvent(i);
    for(int u=0; u<ncamera; u++)
    {
      tempimg = (TH2F*)d.event()->ccdData(u)->Clone("backimg");
      // clean up the bias frame
      MaxCamImageTools::meanRMSNoOutliers(tempimg, bias_mean[u], bias_rms[u]);
      if (dbg_level & dbg_bias)
      {
        cout << "     bias image mean (no outliers): " << bias_mean[u] << endl;
        cout << "     bias image rms (no outliers): " << bias_rms[u] << endl;
      }
      // make sure the bias is reasonable
      if( (MaxCamImageTools::getRMS(tempimg) < sigmathr[u]) )
      {
        // kill 5 sigma outliers in the bias frame, iterate 3 times
        for (int it=0; it<3; it++) {
          double mean, rms;
          mean=MaxCamImageTools::getMean(tempimg);
          rms=MaxCamImageTools::getRMS(tempimg);
          if (conf->isMC()) { // MC
            MaxCamImageTools::killPixels(tempimg,mean+10.*rms);
          } else {
            MaxCamImageTools::killPixels(tempimg,mean*1.25);
          }
          bias_mean[u]=mean;
          bias_rms[u]=rms;
        }
        if(nframes[u]==0) 
          secondarybias[u]=(TH2F*)tempimg->Clone("average");
        else
          secondarybias[u]->Add(tempimg,1);
        nframes[u]++;
        }
    }
  }
  gROOT->Delete("backimg;*");
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

static int parse_args( int nargs, char ** args, TString * files, TString * keys, TString * out, CleanSkimConfig ** cfg)
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
           if (dbg_level & dbg_error) cerr << "Non number was passed to debug level " << endl; 
           return 1; 
        }
      }

      if(strcmp(args[i],"-c")==0) 
      {
        free(*cfg); 
        *cfg = new CleanSkimConfig(args[++i]); 
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
          if (dbg_level & dbg_error) cerr << "Unexpected Extra Argument. Aborting!" << endl;
          return 1; 
      }
    }
  }
  
  return 0; 
}



