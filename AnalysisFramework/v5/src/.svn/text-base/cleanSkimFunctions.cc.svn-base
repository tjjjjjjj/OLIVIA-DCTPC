#include "cleanSkimFunctions.hh"
#include "MaxCamImageTools.hh"
#include <stdio.h>
#include <dirent.h>
#include "WaveformTools.hh"
#include "WaveformAnalysis.hh"
#include "CspWaveform.hh"
#include "PMTWaveform.hh"
#include "FastWaveform.hh"
#include "PMTPulse.hh"
#include "DmtpcIterators.hh"
#include "DmtpcDataConverter.hh"
#include "CspPulse.hh"
#include "FastPulse.hh"

#define DEBUGOUT 

TH2 * cleanSkimFunctions::cleanImage(const TH2 * raw_img, const TH2 * bias, int & nkilled, double raw_rms,  double raw_mean, const TH2 * raw_overscan, double raw_os_mean, double os_bias_mean, const CleanSkimConfig * conf)
{

  TH2F * tempimg = (TH2F*) raw_img->Clone("tempimg");

  image_clean_method_t meth = conf->getImageCleanMethod(); 

  if (meth == IMAGE_CLEAN_BIAS_SUBTRACT_FIRST || meth == IMAGE_CLEAN_MEDIAN_DIFFERENCE)
  {
    tempimg->Add(bias,-1);
  }



  double perpx = 0;


  //Left in for legacy purposes, but this seems silly
  if (conf->isMC() && !conf->isNoKill())
  {
     nkilled = MaxCamImageTools::killLonePixels2(tempimg, raw_mean + 10*raw_rms);
  }
  else if (!conf->isNoKill() && meth == IMAGE_CLEAN_BIAS_SUBTRACT_FIRST || meth == IMAGE_CLEAN_TRADITIONAL)
  {
     nkilled = MaxCamImageTools::killLonePixels2(tempimg,conf->getOutlierFactor()*raw_mean);
     if (conf->hasOverscan()) 
     {
        TH2F * ostmp =  (TH2F*) raw_overscan->Clone("tempos"); 
        MaxCamImageTools::killLonePixels2(ostmp, conf->getOutlierFactor() * raw_os_mean); 
        perpx = MaxCamImageTools::getMean(ostmp) - os_bias_mean; 
        delete ostmp; 
     }
  }
  else if (!conf->isNoKill() && meth == IMAGE_CLEAN_MEDIAN_DIFFERENCE)
  {
    nkilled =  MaxCamImageTools::killLonePixelsMedianDifference(tempimg, conf->getOutlierFactor() * raw_rms); 
    if (conf->hasOverscan()) 
    {
       TH2F * ostmp =  (TH2F*) raw_overscan->Clone("tempos"); 
       MaxCamImageTools::killLonePixelsMedianDifference(ostmp, conf->getOutlierFactor() * raw_rms); 
       perpx = MaxCamImageTools::getMean(ostmp) - os_bias_mean; 
       delete ostmp; 
    }
  }

  if (meth == IMAGE_CLEAN_TRADITIONAL)
  {
    tempimg->Add(bias,-1);
  }

  if (!conf->hasOverscan())
  {
    perpx = tempimg->Integral() / (tempimg->GetNbinsX() * tempimg->GetNbinsY()); 
  }

  MaxCamImageTools::subtractPedestal(tempimg, perpx); 

  return tempimg; 
}


int cleanSkimFunctions::loadGainMaps(DmtpcGainMap ** gainmaps, DmtpcSkimDataset * sd,  int n, const CleanSkimConfig * conf, const Dmtpc4ShooterStitcher * stitch, string * camera_ids)
{

  if (!conf->useGainMap()) 
  {
      return 0; 
  }


    TFile gfile(conf->getGainMapFile()); 
    vector<const TH2*> imgs(n); 
    vector<const DmtpcGainMap*> gms(n); 


    for (int u = 0; u < n; u++)
    {
      std::string ser =  camera_ids[u]; 

      if (!stitch)
      {
        gainmaps[u]=  ser == "" ? 0 : (DmtpcGainMap*) gfile.Get(ser.c_str())->Clone();
        gainmaps[u]->getGainMap()->Scale(gainmaps[u]->getGainMap()->GetNbinsX()*gainmaps[u]->getGainMap()->GetNbinsY()/gainmaps[u]->getGainMap()->Integral());
        if (sd) sd->addGainMap((DmtpcGainMap*)gainmaps[u]->Clone());
      }
      else
      {
        int idx = stitch->getIndex(ser.c_str());
        gms[idx] =   (DmtpcGainMap*) gfile.Get(ser.c_str());
        imgs[idx] =  gms[idx]->getGainMap(); 
      }
    }

    gROOT->cd(); 

    if (stitch)
    {
      TH2F* stitched_map =  (TH2F*) stitch->stitch(&imgs); 
      DmtpcGainMap * newgm = new DmtpcGainMap("stitched_gainmap"); 
      newgm->setGainMap(stitched_map); 

      //We need to add spacers now

      vector<vector<TVector3> > spacers;  //r theta w

      double rthresh = conf->getSpacerJoinR(); 
      double theta_thresh = conf->getSpacerJoinTheta(); 

      for (int i = 0; i < n; i++)
      {
        for (int j = 0; j < stitch->getNSpacers(i); j++)
        {
          double m = stitch->getSpacerSlope(i,j); 
          double b = stitch->getSpacerIntercept(i,j); 
          double w = 32 * stitch->getScale(i); 

          double th = stitch->getRotation(i); 
          double x0 = stitch->xCenter(i); 
          double y0 = stitch->yCenter(i); 
          double c = stitch->getScale(i); 

          double theta = atan(m); 
          double global_theta = theta - th; 
          double global_m = tan(global_theta); 
          double global_b = cos(theta) > 0.1 ?  c * (b - y0 + m*x0) / (cos(th) + m*sin(th)) : c * ( b/m-y0/m + x0)/(cos(th)/m + sin(th)); 

          double global_r = global_b / (global_m * global_m + 1); 

          if (isnan(global_b)) continue; 


          int merge = false; 
          for (unsigned k = 0; k < spacers.size(); k++)
          {

            for (unsigned l = 0; l < spacers[k].size(); l++)
            {

               double dtheta = fabs(sin(global_theta - spacers[k][l].y())); 
 //              cout << dtheta << endl; 
               if (dtheta > theta_thresh) continue; 

               double dr = fabs( global_b - spacers[k][l].x()); 
               if (dr <= rthresh)
               {
                  merge = true; 
//                  cout << "MERGE" << endl; 
                  spacers[k].push_back(TVector3(global_b,global_theta,w)); 
                  break; 
               } 
            }

            if (merge) break; 
          }

          if (!merge)
          {
            vector<TVector3> newv; 
            newv.push_back(TVector3(global_b,global_theta,w)); 
            spacers.push_back(newv); 
          }
        }
      }


      for (unsigned k = 0; k < spacers.size(); k++)
      {

        double bsum = 0; 
        double wsum = 0; 
        double sinsum = 0;
        double cossum = 0; 

        unsigned ns = spacers[k].size(); 
        for (unsigned l = 0; l < ns; l++)
        {
          bsum += spacers[k][l].x(); 
          sinsum +=(sin  (spacers[k][l].y())); 
          cossum +=fabs(cos  (spacers[k][l].y())); 
          wsum += spacers[k][l].z(); 
        }
        double b = bsum/ns; 
        double m = sinsum/cossum; 
        double w = wsum/ns; 
//        cout << m << " " << b << endl; 
//        newgm->addSpacer(m,b,w); 
        newgm->addSpacer(0.000001,b,w); 
      }

      gainmaps[0] = newgm; 
      if (sd) sd->addGainMap(newgm); 

    }

    if (sd)  sd->writeGainMaps(); 
    gfile.Close(); 

}


TH2 * cleanSkimFunctions::cleanBiasFrame(const TH2 * in, const CleanSkimConfig * conf)
{

  TH2F * out = (TH2F*) in->Clone(); 
  //kill 5 sigma outliers in the bias frame, iterate 3 times

  if (conf->isMC() || conf->getImageCleanMethod() == IMAGE_CLEAN_BIAS_SUBTRACT_FIRST || conf->getImageCleanMethod() == IMAGE_CLEAN_MEDIAN_DIFFERENCE) 
  {
    return out; 
  }

  for (int it = 0; it < 3; it++)
  {
    double mean = MaxCamImageTools::getMean(out); 
    double rms = MaxCamImageTools::getRMS(out); 

    MaxCamImageTools::killPixels(out, mean * conf->getOutlierFactor()); 

  }

  return out; 
}

TH2 * cleanSkimFunctions::medianBiasFrameStack(TTree * tree, int frame, const CleanSkimConfig * conf, bool os) 
{
  char branchname[32]; 
  sprintf(branchname, "%s%d", os ? "os_bias" : "bias", frame); 
  TH2S* bias = 0; 
  int n = tree->GetEntries(); 
  tree->SetBranchAddress(branchname, &bias); 
  TreeIterator<TH2S*> beg(tree,&bias,0); 
  TreeIterator<TH2S*> end(tree,&bias,n); 
  TH2 * median = MaxCamImageTools::histStackNthElement<unsigned short> (beg,end,n/2, n, 2<<16,true); 

  //median will be a TH2S but we want a TH2F
  

  sprintf(branchname,"%s%d",os ? "biasFrameOverscan":"biasFrame",frame+1); 

  TH2F * answer = DmtpcDataConverter::ccdExpand((TH2S*) median,branchname); 
  delete median; 
  return answer; 
}


Dmtpc4ShooterStitcher * cleanSkimFunctions::loadStitch(int run, const CleanSkimConfig * conf)
{

  const char * stitchFile = 0; 
  string buf; 

  
  if (strcasecmp(conf->getStitchFile(), "auto"))
  {
    stitchFile = conf->getStitchFile(); 
  }
  else
  {

    const char * stitch_dir = conf->getStitchDir();  

    DIR * dp; 
    struct dirent * dirp; 

    if ( (dp = opendir(stitch_dir)) == NULL)
    {
      cerr << "Could not open stitch_dir: " << stitch_dir << endl; 
      return NULL; 
    }

    vector<int> stitch_runs; 

    while (dirp = readdir(dp))
    {
         
        int i;
        if (sscanf(dirp->d_name,"dmtpc_4sh_%05dstitch.root", &i))
        {
          stitch_runs.push_back(i); 
        }
    }

    closedir(dp); 


    sort(stitch_runs.begin(), stitch_runs.end()); 

     
    int r = -1; 
    for (vector<int>::iterator it = stitch_runs.begin(); it!= stitch_runs.end(); it++)
    {
      if (*it > run) break; 
      r = *it; 
    }

    if (r == -1)
    {
      cerr << "Could not find stitch file for run " << run << endl; 
    }
    
    buf.resize(strlen(stitch_dir) + 64); 
    sprintf(&buf[0], "%s/dmtpc_4sh_%05dstitch.root", stitch_dir, r); 
    stitchFile = buf.c_str(); 
  }


  //Try to open file

  TFile f(stitchFile); 

  if (!f.IsOpen())
  {
    cerr << "Could not open file " << stitchFile << endl; 
    return NULL; 
  }

  gROOT->cd(); 
  Dmtpc4ShooterStitcher * stitch = (Dmtpc4ShooterStitcher*) (f.Get("stitch")->Clone("cs_stitch")); 
  

  return stitch; 
}


bool cleanSkimFunctions::checkPartialSpark(const TH2 * img, const CleanSkimConfig * conf)
{
   TH2C * edges = MaxCamImageTools::edgeDetect(img, conf->getPartialBlur(), 
                                               conf->getPartialLowThresh(),
                                               conf->getPartialHighThresh());  

   int nbinsx = edges->GetNbinsX(); 
   TH1D * edges_proj = edges->ProjectionX("edges_p",2,nbinsx-1); 


   double *  edgesval= new double[nbinsx]; 
   int  * edgesind= new int[nbinsx]; 

   for (int j = 2; j  < edges_proj->GetNbinsX(); j++)
   {
     edgesval[j-1] = edges_proj->GetBinContent(j); 
   }

   edgesval[0] = 0; 
   edgesval[nbinsx-1] = 0; 


   TMath::Sort(img->GetNbinsX(), edgesval, edgesind); 

   bool spark = false; 
   if (edgesval[edgesind[0]] > conf->getPartialNPrimaryThresh())
   {
     int p = 1; 
     bool secedges = false; 

     while (edgesval[edgesind[p]] > conf->getPartialNSecondaryThresh())
     {
       if (abs(edgesind[p] - edgesind[0]) >= conf->getPartialDistanceLow() && 
           abs(edgesind[p] - edgesind[0]) <= conf->getPartialDistanceHigh()) 
       {
          secedges = true; 
       }

       p++; 
     }

     spark = !secedges; 

 }

 delete edges; 
 delete edges_proj; 
 delete edgesval; 
 delete edgesind; 


 return spark; 
}



void cleanSkimFunctions::populateSparkRef(TH2 * raw_img, TH2 * biases, const CleanSkimConfig * conf, vector<pair<int,int> > * sparkref_running) 
{


  TH2 * sparkimg = (TH2*) raw_img->Clone("sparkrefimg"); 
  raw_img->Add(biases,-1); 
  int nsat = 0; 

  for (int x = 1; x <= sparkimg->GetNbinsX(); x++)
  {
    for (int y = 1; y <= sparkimg->GetNbinsY(); y++)
    {
      if (sparkimg->GetBinContent(x,y) > conf->getSatThresh())
      {
        nsat++; 
      }
    }
  }


  if (nsat > conf->getNSatThresh())
  {
    for (int x = 1; x <= sparkimg->GetNbinsX(); x++)
    {
      for (int y = 1; y <= sparkimg->GetNbinsY(); y++)
      {
         if (sparkimg->GetBinContent(x,y) > conf->getSatThresh())          
         {
            sparkref_running->push_back(pair<int,int>(x,y)); 
         }
      }
    }

  }

  sparkimg->Delete(); 

}


bool cleanSkimFunctions::checkSpark(const char * cam_id, double mean, double last_mean,  double os_spark_mean, const  CleanSkimConfig * conf)
{
  
   bool spark =  mean / last_mean >  ((CleanSkimConfig*) conf)->getSparkCut(cam_id); 

   if (conf->hasOverscan())
   {
      spark = spark || mean - os_spark_mean > ((CleanSkimConfig*)conf)->getOverscanSparkCut(cam_id); 
   }

   return spark;
}



int cleanSkimFunctions::fillWaveformVectorsInTObjArray( DmtpcDataset & d, TObjArray * wfvlist) 
{

  //  cout << "inside fillWaveformVectorsInTObjArray( "
  //       << &d
  //       << " , " 
  //       << wfvlist
  //       << " ) ..." 
  //       << endl;


  int nch=wfvlist->GetEntries();
  for (int ich=0; ich<nch; ++ich){
    
    TString tsTempWaveformVectorClassType(((TObject*)(*wfvlist)[ich])->IsA()->GetName());
    //    cout << "tsTempWaveformVectorClassType=" << tsTempWaveformVectorClassType << endl;
    if( tsTempWaveformVectorClassType == (TString)"CspWfVector" ){
      
      //      cout << "it's a CspWfVector!" << endl;
      CspWfVector *wv=static_cast<CspWfVector*> ((*wfvlist)[ich]);
      wv->clear();

      fillCspWfVector(d, wv, ich, nch);
      
    } else if ( tsTempWaveformVectorClassType == (TString)"FastWfVector" ) {
      
      //      cout << "it's a FastWfVector!" << endl;
      FastWfVector *wv=static_cast<FastWfVector*> ((*wfvlist)[ich]);
      wv->clear();
      
      fillFastWfVector(d, wv, ich, nch);
      
    } else if ( tsTempWaveformVectorClassType == (TString)"PMTWfVector" ) {
      
      //      cout << "it's a PMTWfVector!" << endl;
      PMTWfVector *wv=static_cast<PMTWfVector*> ((*wfvlist)[ich]);
      wv->clear();
      
      fillPMTWfVector(d, wv, ich, nch);
 
    } else {
      
      //      cout << "it's an else!" << endl;
      WaveformVector *wv=static_cast<WaveformVector*> ((*wfvlist)[ich]);
      wv->clear();
      
      //      fillWaveformVector(d, wv, ich, conf);
      
    }
  }

  //  cout << "leaving fillWaveformVectorsInTObjArray( "
  //       << &d
  //       << " , " 
  //       << wfvlist
  //       << " ) ..." 
  //       << endl;

  return 0;

}


int cleanSkimFunctions::fillCspWfVector(DmtpcDataset & d, CspWfVector *cspwv, int ich, int nch)
{
  //  cout << "inside fillCspWfVector( ... ) ..." << endl;

  // general scope info for this event
  int nent = d.event()->scopeData()->GetEntries(); 
  int ntr = nent/nch; 
  //  cout << "Nchannels: " << nch << endl;
  //  cout << "Ntriggers: " << ntr << endl; 
  //  cout << "Nentries: " << nent << endl; 

  // loop over the triggers for this specific channel, and pull the waveform and use it to 
  // fill a SkimWaveform object for each trigger
  for(int itr=0; itr<ntr; ++itr){ 

    TH1F *hwf=d.event()->scopeData(ich,itr);
    
    CspWaveform cspwf;
    waveform::analysis::analyzeCSP( hwf , cspwf , 0);
    
    // push this SkimWaveform object onto the WaveformVector for this channel
    cspwv->add(cspwf);
  }
  
  //  cout << "cspwv->size()=" << cspwv->size() << endl;
  //  cout << "exiting fillCspWfVector( ... ) ..." << endl;
  
  return 0;
}

int cleanSkimFunctions::fillPMTWfVector(DmtpcDataset & d, PMTWfVector *pmtwv, int ich, int nch)
{
  //  cout << "inside fillPMTWfVector( ... ) ..." << endl;

  // general scope info for this event
  int nent = d.event()->scopeData()->GetEntries(); 
  int ntr = nent/nch; 
  //  cout << "Nchannels: " << nch << endl;
  //  cout << "Ntriggers: " << ntr << endl; 
  //  cout << "Nentries: " << nent << endl; 

  // loop over the triggers for this specific channel, and pull the waveform and use it to 
  // fill a SkimWaveform object for each trigger
  for(int itr=0; itr<ntr; ++itr){ 

    TH1F *hwf=d.event()->scopeData(ich,itr);
    
    PMTWaveform pmtwf;
    waveform::analysis::analyzePMT( hwf , pmtwf , 1.0);
    
    // push this SkimWaveform object onto the WaveformVector for this channel
    pmtwv->add(pmtwf);
  }
  
  //  cout << "pmtwv->size()=" << pmtwv->size() << endl;
  //  cout << "exiting fillPMTWfVector( ... ) ..." << endl;
  
  return 0;
}

int cleanSkimFunctions::fillFastWfVector(DmtpcDataset & d, FastWfVector *fastwv, int ich, int nch)
{
  //  cout << "inside fillFastWfVector( ... ) ..." << endl;
  
  // general scope info for this event
  int nent = d.event()->scopeData()->GetEntries(); 
  int ntr = nent/nch; 
  //  cout << "Nchannels: " << nch << endl;
  //  cout << "Ntriggers: " << ntr << endl; 
  //  cout << "Nentries: " << nent << endl; 
  
  // loop over the triggers for this specific channel, and pull the waveform and use it to 
  // fill a SkimWaveform object for each trigger
  for(int itr=0; itr<ntr; ++itr){ 
    
    TH1F *hwf=d.event()->scopeData(ich,itr);
    
    FastWaveform fastwf;
    //    cout << "hwf=" << hwf << endl;
    //    cout << "ich=" << ich << endl;
    //    cout << "itr=" << itr << endl;
    waveform::analysis::analyzeFast( hwf , fastwf );
    //    cout << "just got done with analyzeFast( ... ) ... " << endl;
    
    // push this SkimWaveform object onto the WaveformVector for this channel
    fastwv->add(fastwf);
  } 
  
  //  cout << "fastwv->size()=" << fastwv->size() << endl;
  //  cout << "exiting fillFastWfVector( ... ) ..." << endl;
  
  return 0;
}

bool cleanSkimFunctions::burnin_test( double x_delta, double y_delta, const CleanSkimConfig * conf)
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
  cerr << "Invalid Burnin Test: " << method << endl; 
  return false; 
}

int cleanSkimFunctions::parse_args( int nargs, char ** args, TString * files, TString * keys, TString * out, char** cfg)
{
  int c = 0; 
  for (int i = 1; i < nargs; i++)
  {
    /* Handle switches */
    if (args[i][0]=='-')
    {
      
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
           cerr << "Unexpected Extra Argument. Aborting!" << endl;
          return 1; 
      }
    }
  }
  
  return 0; 
}


vector<BurninEncoded_t> cleanSkimFunctions::computeTrackBurnin(list<vector<double*> **> * positions,
                                                               list<vector<vector<vector<BurninEncoded_t> > >* > * burnin_temp, 
                                                               int position_offset, double x, double y, int cam, int track, const CleanSkimConfig * conf)
{
    vector<BurninEncoded_t> this_track_burnin; 
    int pos_i = 0; 
    int u = cam; 
    int v = track; 

    //Loop over last up to N_LOOPBACK events in search of things at the same position
    list<vector<double *>**>::iterator pos_iter; //iterator for previous positions
    list<vector<vector<vector<BurninEncoded_t> > >* >::iterator burnin_iter; 
    
    for (pos_iter = positions->begin(), burnin_iter = burnin_temp->begin(); 
         pos_iter != positions->end() ; 
         pos_iter++, burnin_iter++)
        {
           vector<double *> ** test_positions = *pos_iter; 
           vector<vector<vector<BurninEncoded_t> > > * test_burnin = *burnin_iter;
           
           
           //Loop over the tracks in these events
           for(unsigned int z = 0; z<test_positions[u]->size(); z++)
           {
             double * position = (double *)(test_positions[u]->at(z));
             if (position[0]==0.0 && position[1] == 0.0) continue; //ignore events with no tracks
             double x_delta = x - position[0];
             double y_delta = y - position[1];

             if(burnin_test(x_delta, y_delta, conf))
             {
               int test_event_index = position_offset + pos_i; 
               this_track_burnin.push_back(encodeBurnin(z,test_event_index));

               //Mark this event on the other events record...
               int our_index = position_offset + positions->size();
               assert(test_event_index != our_index);

               test_burnin->at(u)[z].push_back(encodeBurnin(v,our_index));
             }
           }    
           pos_i++; 
        }//End Loop over last up to N_LOOKBACK events
      

    return this_track_burnin; 

}


int cleanSkimFunctions::clusterFind(const TH2F * img, const DmtpcEvent * ev, MaxCamClusterImage ** clusti, const DmtpcGainMap * gainmap, const CleanSkimConfig * conf, const DmtpcStitchInfo * stitch,
                                    const double * image_rms, const double * blurred_image_rms,  const double * image_mean, const double * blurred_image_mean,  int indx, int u)
{


  static int nth = 0; 

  TH2F * clustimg = (TH2F*) img->Clone("clustimg"); 

  int ntracks = 0;
  char * algo = conf->getClusterFindingAlgorithm(indx); 

  if (!strcasecmp(algo,"ci"))
  {
    TH2F * baseimage = (TH2F*) img->Clone("baseimage"); 
    const int version = ev->IsA()->GetClassVersion(); 

    (*clusti) = version == 1 ? new MaxCamClusterImage(baseimage,  ev->timeStamp()) 
                          : new MaxCamClusterImage(baseimage,  ev->UTCtimeStamp()) ;

    baseimage->Rebin2D(2,2); 
    baseimage = (TH2F*) MaxCamImageTools::blur(baseimage,1,conf->getBlurAmount()); 

    if (conf->useGainMap() && conf->getClusterFindingAlgorithm())
    {
      ntracks= MaxCamImageTools::findClustersGM(baseimage, 
                                                *clusti, 
                                                conf->getClusterMinSigma(), 
                                                conf->getClusterMaxSigma(), 
                                                conf->getClusterMinSize(),
                                                conf->getClusterMinDist(), 
                                                conf->getClusterJoinMinRxyGlobal(),
                                                conf->getClusterJoinMinRxyCluster(), 
                                                conf->getClusterJoinMaxJoinResidual(),
                                                conf->getClusterJoinLeastSquaresWeight(), 
                                                gainmap,
                                                conf->getClusterJoinSpacerWidth()
                                                );

    }
    else 
    {
      ntracks=  MaxCamImageTools::findClustersCI(baseimage, 
                                                *clusti, 
                                                conf->getClusterMinSigma(),
                                                conf->getClusterMaxSigma(),
                                                conf->getClusterMinSize(),
                                                conf->getClusterMinDist());

    }

    (*clusti)->changeImageWithThreshold(clustimg,conf->getClusterReducedThreshold()*image_rms[u]);  
    // now we must delete baseimage
    gROOT->Delete("baseimage"); 
    gROOT->Delete("copy"); 

  }

  else
  {
    (*clusti) = new MaxCamClusterImage(clustimg, ev->UTCtimeStamp()); 
    if (!strcasecmp(algo,"seed"))
    {


      if (!stitch)
      {
        ntracks = MaxCamImageTools::findClustersGMSeed(clustimg, 
                                                       *clusti, 
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
                                                       gainmap,
                                                       conf->getClusterJoinSpacerWidth(), 
                                                       conf->getSeedClusterFindReproduceV4Bug()
                                                       ); 

      }
      else 
      {
#ifdef DEBUGOUT
  char buf[64]; 
  sprintf(buf,"debug/%d.root",nth); 
#endif
        
        ntracks = MaxCamImageTools::findClustersGMSeedStitch(clustimg, 
                                                             *clusti, 
                                                             stitch, 
                                                             image_mean, 
                                                             blurred_image_mean, 
                                                             image_rms, 
                                                             blurred_image_rms, 
                                                             conf->getClusterSeedThreshold(), 
                                                             conf->getClusterReducedThreshold(), 
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
                                                             gainmap,
                                                             conf->getClusterJoinSpacerWidth()
#ifdef DEBUGOUT
                                                             , buf, 
#else
                                                             , 0, 
#endif                                                              
                                                             conf->getSeedClusterFindReproduceV4Bug()
                                                             ); 





      }

    }

    else if (!strcasecmp(algo,"ad"))
    {
        ntracks = MaxCamImageTools::findClustersADHysteresisGM(clustimg, *clusti, 
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
                                                               gainmap,
                                                               conf->getClusterJoinSpacerWidth());
    }
    else if (!strcasecmp(algo,"ring"))
    {
      char * debug = 0; 
#ifdef DEBUGOUT
      char buf[64]; 
      debug = buf; 
      sprintf(debug,"debug/ring%d.root",nth); 
#endif
      ntracks = MaxCamImageTools::findClustersGMRing(clustimg, *clusti, stitch, image_rms,image_mean, 
                                                     conf->getClusterFindRingSpaceSigma(), 
                                                     conf->getClusterFindRingRMSSigma(), 
                                                     conf->getClusterFindRingCoreThreshHigh(), 
                                                     conf->getClusterFindRingCoreThreshLow(), 
                                                     conf->getClusterFindRingRingThresh(), 
                                                     conf->getClusterFindRingNSigma(), 
                                                     conf->getClusterMinSizeUnbinned(), 
                                                     MaxCamImageTools::BILATERAL_GAUSSIAN,  
                                                     conf->getClusterFindRingNCleanup(), 
                                                     conf->getClusterMinDist(), 
                                                     conf->getClusterJoinMinRxyCluster(),
                                                     conf->getClusterJoinMinRxyGlobal(), 
                                                     conf->getClusterJoinMaxJoinResidual(),
                                                     conf->getClusterJoinLeastSquaresWeight(), 
                                                     gainmap,
                                                     conf->getClusterJoinSpacerWidth(),
                                                     debug, 12

          ); 

    }
 
    else 
    {
       ntracks = 0; 
       std::cout << "NO CLUSTER FINDING ALGORITHM!!!!" << std::endl; 

    }

    if (!conf->stitch() && strcasecmp(algo,"ring")) {(*clusti)->applyRedThreshold(conf->getClusterReducedThreshold() * image_rms[u]); }

  }


  nth++; 
  return ntracks; 
}


/** Reconstructed Track Parameters Go Here **/ 
void cleanSkimFunctions::fillTrackInfo(DmtpcSkimEvent * tmp_event, int u, int v, MaxCamClusterImage * clusti, const CleanSkimConfig * conf, const DmtpcGainMap * gm, const Dmtpc4ShooterStitcher * stitch, const char * cam_id, DmtpcEvent * orig_ev)
{
  

  //Energy
  tmp_event->_E[u][v] = clusti->getIntegral(v); 
  if (conf->useGainMap()) tmp_event->_EGainMap[u][v] = clusti->getIntegralWithGainMap(v,gm); 

  //Phi
  switch(conf->getPhiAlgorithm())
  {
    case 1: tmp_event->_phi[u][v] = clusti->getPhi(v); break;
    case 2: tmp_event->_phi[u][v] = clusti->getPhi2(v); break;
    case 3: tmp_event->_phi[u][v] = clusti->getPhi3(v); break; 
    case 4: tmp_event->_phi[u][v] = clusti->getPhi4(v); break;
    default: tmp_event->_phi[u][v] = 0; 
    cout << "Warning: invalid phi algorithm selected" << endl; 
  }

  //Skewness
  tmp_event->_skewness[u][v] = clusti->getSkewness(v,tmp_event->_phi[u][v]); 

  //Theta (not really)
  tmp_event->_theta[u][v] = 0; 

  //Range
  double xb,yb,xe,ye;
  switch(conf->getRangeAlgorithm())
  {                  
      case 1: tmp_event->_range[u][v] = clusti->getLength(v,xb,yb,xe,ye); break;
      case 2: tmp_event->_range[u][v] = clusti->getLength2(v,tmp_event->_phi[u][v],1); break;
      default: tmp_event->_range[u][v] = 0;
      cout << "Warning: invalid range algorithm selected" << endl; 
  }

  //Diffused Range
	tmp_event->_diffusedRange[u][v]=clusti->getDiffusedLength(v,xb,yb,xe,ye);

  //Ellipse Axes
	clusti->getEllipseAxes(v,tmp_event->_majoraxis[u][v],tmp_event->_minoraxis[u][v]);
		
  //Position
  clusti->getXY(v,tmp_event->_x[u][v],tmp_event->_y[u][v]);
  double r_offset = stitch ? 0 : 512; //TODO: unhardcode this
  tmp_event->_r[u][v] = sqrt(pow(tmp_event->_x[u][v] - r_offset,2) + pow(tmp_event->_y[u][v] - r_offset,2)); 

  //Cluster Mean
  tmp_event->_cluster_mean[u][v] = clusti->getMean(v);
  
  //Cluster RMS
  tmp_event->_cluster_rms[u][v] = clusti->getRMS(v,tmp_event->_cluster_mean[u][v]);

  //Energy Density (is this different from mean?) 
  tmp_event->_energy_density[u][v] = clusti->getEnergyDensity(v); 

  //Npixel (red) 
  tmp_event->_npixel[u][v] = clusti->getCluster(v).size(); 
  tmp_event->_npixel_red[u][v] = clusti->getClusterRed(v).size(); 

  //Max Pixel
  int maxBin; 
  tmp_event->_maxpixel[u][v] = clusti->getMax(v,&maxBin);

  //Neighbors
  tmp_event->_neighbors[u][v] = clusti->getNumNeighbors(v,conf->getClusterMinSigma(),maxBin);

  //Cygnus Angle
  MaxCamClusterImage::CAMERA_ORIENTATION orientation = conf->getCameraOrientation(cam_id);         
  tmp_event->_cygnus_angle[u][v] = clusti->getCygnusAngle(v, conf->getNorthAngle(),
                                              orientation,
                                              conf->getLatitude(), conf->getLongitude(),
                                              tmp_event->_phi[u][v],
                                              tmp_event->_theta[u][v]+ TMath::Pi()/2. //TODO: when we figure out theta, remove this argument
                                              );
                                              

  //Right Ascension, Declination, Galactic Coords
  clusti->getRADec( tmp_event->_phi[u][v], tmp_event->_theta[u][v] + TMath::Pi()/2., //TODO: ditto as above
                   orig_ev->timeStamp(), 
                   conf->getLatitude(), conf->getLongitude(), 
                   conf->getNorthAngle(), orientation, 
                   tmp_event->_ra[u][v],tmp_event->_dec[u][v],tmp_event->_glat[u][v],tmp_event->_glon[u][v]); 


  //Moments
  for (int m = 0; m <4; m++)
  { 
    tmp_event->_moments[u][m][v] = clusti->getMoment(v,m+1,tmp_event->_phi[u][v], 4,"pixelPerBin"); 
    tmp_event->_transverse_moments[u][m][v] = clusti->getMoment(v,m+1,tmp_event->_phi[u][v]+TMath::Pi()/2., 4,"pixelPerBin"); 
  }


  //edge / veto
  tmp_event->_edge[u][v] = !stitch ? clusti->hitsEdge(v) :  clusti->hitsVeto(v, stitch->innerRadius(0), stitch->outerRadius(0));  

  //Stitch only stuff
  if (stitch)
  {
    tmp_event->_inactive[u][v] = clusti->hitsInactive(v, stitch->outerRadius(0)); 
    tmp_event->_crossing[u][v] = clusti->crossesCameras(v, orig_ev->ccdData(0)->GetNbinsX() == 256 ? 
                                                           stitch->getStitchInfo256() : stitch->getStitchInfo1024()); 
  }


}

void cleanSkimFunctions::updateGlobalSparkRef(vector<pair<int,int> > **sparkref_running, list<vector<pair<int,int> > * > * sparkref, int ncamera)
{

   vector<pair<int,int> > * this_sparkref = new vector<pair<int,int> >[ncamera]; 
   for (int u = 0; u < ncamera; u++)
   {
     for (unsigned int v = 0; v < (*sparkref_running)[u].size(); v++)
     {
       this_sparkref[u].push_back((*sparkref_running)[u][v]); 
     }
     
     assert (this_sparkref[u].size()==(*sparkref_running)[u].size()); 
   } 
   sparkref->push_back(this_sparkref); 

}

void cleanSkimFunctions::updateGlobalBurnin( list<vector<double*>**> * positions, list<vector<vector<vector<BurninEncoded_t> > > *> * burnin_temp, 
                                            list<vector<vector<vector<BurninEncoded_t> > > *> * burnin, DmtpcSkimEvent * tmp_event, const CleanSkimConfig * conf, int ncamera, vector<vector<vector<BurninEncoded_t> > > * this_event_burnin) 
{

   //Add self to position list, removing first entry if list is already N_LOOKBACK
    vector<double *> ** these_positions = new vector<double *>*[ncamera];
    for (int u = 0; u < ncamera; u++)
    {
      these_positions[u] = new vector<double *>; 
      for (int v = 0; v < tmp_event->ntracks(u); v++)
      {
         double * this_pos = new double[2];
         this_pos[0] = tmp_event->x(u,v); 
         this_pos[1] = tmp_event->y(u,v); 
         these_positions[u]->push_back(this_pos); 
      }
    } 

    positions->push_back(these_positions); 

    //Store the burnin_info in the buffer until we go past it enough times
    burnin_temp->push_back(this_event_burnin); 
    

    //After we start filling the buffer, 
    // it is time to clear the first entry
    // and to write out the first event
  
   if (positions->size() >(unsigned int) conf->getBurninNumEvents())
   {

     //Delete the stored positions. 
     for (int u = 0; u < ncamera; u++)                 
     {
        for (unsigned int z = 0; z < positions->front()[u]->size(); z++)
        {
           delete positions->front()[u]->at(z);
        }
        delete positions->front()[u];
     }
     delete positions->front();
     positions->pop_front(); 
    
     
     burnin->push_back(burnin_temp->front());
     burnin_temp->pop_front(); 
     
    gROOT->cd(); 
   }
}




