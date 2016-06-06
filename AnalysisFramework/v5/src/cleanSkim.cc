#include "cleanSkim.hh"

using namespace std;

int cleanSkim(DmtpcDataset & d, int runnum, DmtpcSkimDataset & sd, bool delete_intermediate, 
              const CleanSkimConfig * conf, unsigned algo_index)
{


  /************************/
  /**** Initialization  ***/ 
  /************************/

  /** Temporary tree without temporal information **/ 
  TString tmpoutfilename = TString(sd.getFileName()).ReplaceAll(".root","_tmp.root");
  TFile* tmpoutfile = new TFile(tmpoutfilename,"RECREATE");
  tmpoutfile->SetCompressionLevel(conf->getCompressionSettings()); 
  sd.file()->SetCompressionLevel(conf->getCompressionSettings()); 
  TTree * tmpskimtree; 
  tmpskimtree = new TTree(algo_index == 0 ? "skim" : conf->getClusterFindingAlgorithm(algo_index),"Temp Skim Tree"); 

  DmtpcSkimEvent * tmp_event = new DmtpcSkimEvent; 

  /* Load first event to get information */ 
  d.getEvent(0); 
  const int ncam_in = d.event()->ccdData()->GetEntries(); 
  const int nbinsx =  d.event()->ccdData(0)->GetNbinsX();
  const int nbinsy =  d.event()->ccdData(0)->GetNbinsY();
  const double xmin =  d.event()->ccdData(0)->GetXaxis()->GetXmin();
  const double xmax =  d.event()->ccdData(0)->GetXaxis()->GetXmax();
  const double ymin =  d.event()->ccdData(0)->GetYaxis()->GetXmin();
  const double ymax =  d.event()->ccdData(0)->GetYaxis()->GetXmax();

  TH2F * blank = new TH2F("blank","blank", nbinsx, xmin, xmax, nbinsy, ymin,ymax); 
  
  const int eventVersion = d.event()->IsA()->GetClassVersion();

  double last_spark_mean[ncam_in] ;
  MaxCamClusterImage * clusti = 0; 

  /* Support processing on datasets with no actual images */ 
  bool fake_image = false; 
  if (ncam_in ==0 || d.event()->ccdData(0)==NULL || d.event()->ccdData(0)->Integral() ==0)
  {
    fake_image = true; 
  }

 //Get camera id's 
 string camera_ids[ncam_in];  
 for (int u = 0; u < ncam_in; u++)
 {
    camera_ids[u] = string( ((MaxCamConfig*)d.event()->ccdConfig(u))->serialNumber.Data()); 
    if (camera_ids[u] == "")
    {
      camera_ids[u] = std::string(conf->getFallBackCameraID(u)); 
      if (camera_ids[u]!="")
      {
        std::cout << "Fell back to camera id: " << camera_ids[u] << std::endl; 
      }
    }
 }

 int ncam_out = conf->stitch() ? 1 : ncam_in; 
 tmp_event->initCamVectors(ncam_out); 
// tmp_event->initTrackVectors(ncam_out, 15); //backwards compatibility
 tmp_event->_clusters->Expand(ncam_out); 
 tmp_event->_ncamera = ncam_out; 
 tmpskimtree->Bronch("tmp_event","DmtpcSkimEvent",&tmp_event,32000,1); 
 gROOT->cd();

 tmp_event->_clusters->SetOwner(kTRUE); 

 /* Open up file for bias frame saving */
 TFile * biasoutfile;
 vector<TH2*> biases(ncam_in); 
 vector<TH2*> os_biases(ncam_in); 
 vector<double> bias_mean(ncam_in); 
 vector<double> os_bias_mean(ncam_in); 
 TTree * biastree; 

 if (algo_index == 0) // We only handle bias frame once
 {
    TString biasoutfilename = TString(sd.getFileName()).ReplaceAll("skim.root","bias.root");
    biasoutfile = new TFile(biasoutfilename,"RECREATE");

    TH2F * biasframe = 0;
    TH2F * osframe;
    biastree = new TTree("bias","Bias Information");
    biastree->Branch("biasframe","TH2F",&biasframe,128000,0);
    gROOT->cd();

    for (int u = 0; u < ncam_in; u++)
    {
      //Get cleaned bias frame

      if (fake_image)
      {
        biasframe = new TH2F;   
      }
      else
      {
        TH2F * rawframe = conf->getBiasCleanMethod() == BIAS_CLEAN_MEDIAN_STACK  
                          ? (TH2F*) cleanSkimFunctions::medianBiasFrameStack((TTree*) d.file()->Get("allbias"),u,conf,false)
                          : (TH2F*) d.getBiasFrame(u+1);                              

        biasframe =  (TH2F*)cleanSkimFunctions::cleanBiasFrame(rawframe, conf); 

        if (conf->getBiasCleanMethod() == BIAS_CLEAN_MEDIAN_STACK)
        {
          delete rawframe; 
        }
      }
      
      bias_mean[u] = MaxCamImageTools::getMean(biasframe);     
      last_spark_mean[u] = bias_mean[u]; 

      if (conf->hasOverscan())
      {

         TH2F * rawframe = conf->getBiasCleanMethod() == BIAS_CLEAN_MEDIAN_STACK  
                          ? (TH2F*) cleanSkimFunctions::medianBiasFrameStack((TTree*) d.file()->Get("allbias"),u,conf,true)
                          : (TH2F*) d.getBiasFrameOverscan(u+1);                              
         osframe = (TH2F*) cleanSkimFunctions::cleanBiasFrame(rawframe, conf); 

         if (conf->getBiasCleanMethod() == BIAS_CLEAN_MEDIAN_STACK)
         {
            delete rawframe;   
         }

         os_bias_mean[u] = MaxCamImageTools::getMean(osframe); 
         os_biases[u] = osframe; 
      }
       
      biasoutfile->cd(); 
      biastree->Fill(); 
      gROOT->cd(); 
      biases[u] = biasframe; 
    }
  }

   cout << "All preclean activities done" << endl;


  /* Running list of positions so we don't have to read back through skim file
     all the time. It is a List of a Vector[2] of a double[2].*/
  list<vector<double *>**> positions; 
  

  /* List of spark ref excluded pixels */
  vector<pair<int,int> > * sparkref_running = new vector<pair<int,int> >[ncam_in]; 
  list<vector<pair<int,int> >*>  sparkref; 

  list<vector<vector<vector<BurninEncoded_t> > >*> burnin_temp; //Temporary list that is updated a lot
  list<vector<vector<vector<BurninEncoded_t> > >*> burnin; //permanent list 
  

  const int nev = d.tree()->GetEntries();

  Dmtpc4ShooterStitcher * stitch = conf->stitch() ? cleanSkimFunctions::loadStitch(runnum, conf) : 0 ; 
  if (stitch) sd.writeStitch(stitch); 

  //Load Gain Maps; 
  vector<DmtpcGainMap *> gainmaps (ncam_out);
  cleanSkimFunctions::loadGainMaps(&(gainmaps[0]), &sd,  ncam_in, conf, stitch, camera_ids); 
  
  vector<int> last_spark(ncam_in,0); 

  //initialize waveform stuff
  //
  for(Int_t ich=0; ich<conf->getNChannelsPerTrigger(); ++ich){
    waveform::tools::addWaveformVectorToTObjArray( tmp_event->_waveform_vectors,
						   ich,
						   conf->getChannelId(ich),
						   conf->getChannelType(ich) );
  }


  tmp_event->_waveform_vectors->SetOwner(true); 
  tmp_event->_runNumber = runnum; 

  /************************/
  /**** Main Event Loop ***/ 
  /************************/
  for (int i = 0; i < nev; i++)
  {
     cout << i << endl;
     d.getEvent(i); 

     tmp_event->_clusters->Clear(); 

     vector<vector<vector<BurninEncoded_t> > > * this_event_burnin= new vector<vector<vector<BurninEncoded_t> > >;

     cleanSkimFunctions::fillWaveformVectorsInTObjArray( d , tmp_event->_waveform_vectors );

    //cout << "waveform_vectors->GetEntries()=" << tmp_event->_waveform_vectors->GetEntries() << endl;
    //cout << "waveform_vectors->GetSize()=" << tmp_event->_waveform_vectors->GetSize() << endl;

    TCutG os_cut("os_cut",5);
    if (conf->hasOverscan())
    {
      os_cut.SetPoint(0,0,1028);
      os_cut.SetPoint(1,1024,1028);
      os_cut.SetPoint(2,1024,1032);
      os_cut.SetPoint(3,0,1032);
      os_cut.SetPoint(4,0,1028);
    }

     int nspark = 0; 
     vector<const TH2*> cleaned_images(ncam_in); 
     vector<bool> spark(ncam_in,false); 
     vector<double> image_rms(ncam_in,0); 
     vector<double> cleaned_image_rms(ncam_in,0); 
     vector<double> cleaned_blurred_image_rms(ncam_in,0); 
     vector<double> cleaned_blurred_image_mean(ncam_in,0); 
     vector<double> cleaned_image_mean(ncam_in,0); 
     vector<double> raw_rms(ncam_in,0); 
     vector<double> image_mean(ncam_in,0); 
     vector<double> os_mean(ncam_in,0); 
     vector<double> os_rms(ncam_in,0); 
     vector<int> nkilled(ncam_in,0); 

     /**************************************/
     /**** Cleaning and Spark Detection  ***/ 
     /**************************************/
     for(int u=0; u<ncam_in; u++)
     {
        TH2F * raw_img = d.event()->ccdData(u); 
        TH2F * raw_overscan =  conf->hasOverscan() ? d.event()->overscan(u) : 0 ; 

        cout << "\t CLEANING "<< u << "" << endl;
        raw_rms[u] = MaxCamImageTools::getRMS(raw_img); 
        // keep track of the last frame for the spark cut


        //Check for spark
        
        double spark_mean = fake_image ? 0 : MaxCamImageTools::getMean(raw_img); 
        double os_spark_mean = conf->hasOverscan() ? MaxCamImageTools::getMean(raw_overscan, &os_cut) : 0; 

        if (!fake_image && !conf->isMC()) 
        {
          spark[u] = cleanSkimFunctions::checkSpark(camera_ids[u].c_str(), spark_mean, last_spark_mean[u], os_spark_mean, conf); 
          last_spark_mean[u] = spark_mean; 


         //Populate Sparkref if necessary 
          if (spark[u]) cleanSkimFunctions::populateSparkRef(raw_img, biases[u], conf, &sparkref_running[u]); 
          if (spark[u]) nspark++; 
        }
         
        //Clean Images
        int clean_idx = stitch ? stitch->getIndex(camera_ids[u].c_str()) : u; 
        cleaned_images[clean_idx] = fake_image ? (TH2F*) blank->Clone() : 
                                         spark[u] ? (TH2F*) raw_img->Clone("sparkimg") : 
                                                    ((TH2F*) cleanSkimFunctions::cleanImage(raw_img, biases[u], nkilled[u],
                                                    raw_rms[u], spark_mean, raw_overscan, os_spark_mean, os_bias_mean[u], conf)); 



        //Get image, overscan, mean and rms
        MaxCamImageTools::meanRMSNoOutliers(raw_img, image_mean[u], image_rms[u]); 
        MaxCamImageTools::meanRMSNoOutliers((TH2F*)cleaned_images[clean_idx], cleaned_image_mean[clean_idx], cleaned_image_rms[clean_idx]); 

        // we need blurred rms and mean for stitched seedClusterFind
        if (stitch && !strcmp(conf->getClusterFindingAlgorithm(algo_index),"seed"))
        {
          TH2 * blurred = MaxCamImageTools::gaussianBlur(cleaned_images[clean_idx], conf->getClusterBlurRadius(), conf->getGaussianBlurAmount()); 
          MaxCamImageTools::meanRMSNoOutliers(blurred,cleaned_blurred_image_mean[clean_idx],cleaned_blurred_image_rms[clean_idx]); 
          delete blurred; 
        }

        if (conf->hasOverscan()) MaxCamImageTools::meanRMSNoOutliers(raw_overscan, os_mean[u], os_rms[u]); 


        
        //Check for partial spark
        if (!fake_image && !conf->isMC() && !spark[u])
        {
//          cout << "Checking partial " << endl; 
          spark[u] = cleanSkimFunctions::checkPartialSpark(cleaned_images[clean_idx], conf); 
        }

        if (spark[u]) last_spark[u] = 0; 
        else last_spark[u]++; 
     }

     cout << "\t DONE CLEANING" << endl;

     TH2F * stitched = stitch ? (TH2F*) stitch->stitch(&cleaned_images, conf->getInterpolationMethod()) : 0; 
    
     tmp_event->_stitched = stitch ? true : false; 

     tmp_event->clearTrackVectors(); 
     for (int u = 0; u < ncam_out; u++)
     {

       TH2F * img = stitch ? stitched : (TH2F*)  cleaned_images[u]; 
       vector< vector<BurninEncoded_t> >  this_camera_burnin ;

       //Average camera for stitched images for output tree for now
       tmp_event->_integral[u] = stitch ? std::accumulate(raw_rms.begin(), raw_rms.end(),0.) / ncam_in : raw_rms[u]; 
       tmp_event->_image_rms[u] = stitch ? std::accumulate(image_rms.begin(), image_rms.end(),0.) / ncam_in : image_rms[u]; 
       tmp_event->_image_mean[u] = stitch ? std::accumulate(image_mean.begin(), image_mean.end(),0.) / ncam_in : image_mean[u]; 
       //tmp_event->_os_mean[u] = stitch ? std::accumulate(os_mean.begin(), os_mean.end(),0) / ncam_in : os_mean[u]; 
       //tmp_event->_os_rms[u] = stitch ? std::accumulate(os_rms.begin(), os_rms.end(),0) / ncam_in : os_rms[u]; 
       tmp_event->_pixels_killed[u] = stitch ? std::accumulate(nkilled.begin(), nkilled.end(),0)  : nkilled[u]; 
       tmp_event->_lastspark[u] = stitch ? *(std::min_element(last_spark.begin(), last_spark.end())) : last_spark[u]; 
       if (!stitch)
       {
         tmp_event->_cameraSerialNumber[u] = camera_ids[u]; 
       }
       else
       {
         tmp_event->_cameraSerialNumber[u] = "Stitched"; 
//         for (int c = 0; c < ncam_in; c++)
//         {
//           if (c > 0) 
//           {
//             tmp_event->_cameraSerialNumber[u] += "&"; 
//           }
//           tmp_event->_cameraSerialNumber[u] += camera_ids[c]; 
//         }
       }


       if (spark[u] || (stitched && nspark > 0))
       {
          cout << "Spark!" << endl; 
          //Just write out the spark 
          tmp_event->_ntracks[u] = 0; 
          tmp_event->_spark[u] = 1; 
          TH2F * sparkimg = (TH2F*) img->Clone("sparkimg"); 
          clusti = eventVersion == 1 ? new MaxCamClusterImage(sparkimg,d.event()->timeStamp()) 
                                     : new MaxCamClusterImage(sparkimg, d.event()->UTCtimeStamp()); 
          tmp_event->addTrackVector(15); //backwards compatibility
       }
       else if (fake_image)
       {
          tmp_event->_spark[u] = 0; 
          tmp_event->_ntracks[u] = 0;
          tmp_event->addTrackVector( 15); //backwards compatibility
       }
       else
       {
  
         tmp_event->_spark[u] = 0; 
         const DmtpcStitchInfo * sinfo = 0; 
         if (stitch) sinfo = stitch->getStitchInfo(cleaned_images[0]->GetNbinsX());
         tmp_event->_ntracks[u] = cleanSkimFunctions::clusterFind(img, d.event(), &clusti, gainmaps[u], conf, sinfo,  &(cleaned_image_rms[0]),  &(cleaned_blurred_image_rms[0]) ,  &(cleaned_image_mean[0]), &(cleaned_blurred_image_mean[0]), algo_index, u);

         tmp_event->addTrackVector( tmp_event->_ntracks[u] > 15 ? tmp_event->_ntracks[u] : 15 ); //backwards compatibility
//         cout << tmp_event->_range.size() << endl; 
//         cout << tmp_event->_range[0].size() << endl; 
         /************************/
         /*** Loop over tracks ***/ 
         /************************/

         for (int v = 0; v < tmp_event->_ntracks[u]; v++)
         {
 
            cleanSkimFunctions::fillTrackInfo(tmp_event, u,v,clusti,conf,gainmaps[u], stitch, camera_ids[u].c_str(),d.event());  

      
            this_camera_burnin.push_back(cleanSkimFunctions::computeTrackBurnin(&positions, &burnin_temp, 
                                                                               i, tmp_event->x(u,v), 
                                                                               tmp_event->y(u,v), u, v,conf)); 
         }

     }


     if (algo_index > 0) 
     {
       clusti->forgetImage(); 
     }
     else if(conf->getHistSaveType() != HIST_SAVE_TYPE_FLOAT) 
     {
       switch (conf->getHistSaveType())
       {
         case HIST_SAVE_TYPE_SHORT: 
            clusti->changeHistType('S');  
            break; 
         case HIST_SAVE_TYPE_INT: 
            clusti->changeHistType('I');  
            break; 
        default: 
            cerr << "Unrecognized histogram type "<< endl; 
       }

     }
     else if (conf->getRoundAmount() != 0)
     {
      clusti->roundValues(conf->getRoundInClusters(), conf->getRoundAmount()); 
     }   


     tmp_event->_clusters->AddAt(clusti,u);
     this_event_burnin->push_back(this_camera_burnin);


   } /*End loop over cameras */


   if (!fake_image) cleanSkimFunctions::updateGlobalSparkRef(&sparkref_running, &sparkref, ncam_out); 


   //Fill Temp Tree 
   tmp_event->_eventNumber = i; 
   tmpoutfile->cd(); 
   tmpskimtree->Fill(); 
   gROOT->cd(); 

   if (stitched) delete stitched; 
   for (int u = 0; u < cleaned_images.size(); u++) delete cleaned_images[u]; 

   //Update global burnin list 
   if (!fake_image) cleanSkimFunctions::updateGlobalBurnin(&positions,&burnin_temp, &burnin, tmp_event, conf, ncam_out, this_event_burnin); 

  } //End loop over events


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
 

  //Save configuration
  if (algo_index == 0)
  {
    std::stringstream str; 
    conf->print(str); 
    sd.setConfig(str.str().c_str()); 
    sd.writeConfig(); 
  }

  sd.mergeTempTreesv5(tmpskimtree,&burnin, &sparkref, runnum);    

  if (conf->useGainMap())
  {
    for (int i = 0; i < ncam_out; i++) 
    {
      if (gainmaps[i]) gainmaps[i]->Delete(); 
    }
  }
  
  
  
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
  if (algo_index == 0) 
  {
    biasoutfile->Close();
  }
  
  if (delete_intermediate)
  {
     cout << "Deleting " << tmpoutfilename << endl;
    unlink(tmpoutfilename.Data());
  }

}

/* Main function */ 

int main(int argn, char ** argv)
{
   int return_val = 0; 

  /* Argument Parsing */ 
   TString rawdatafiles = "files.txt";
   TString keyfilename = "keys.txt";
   TString outdir = "./skim/"; 
   //TString sumdir = "./sum/";
   char * config = 0;
   if(cleanSkimFunctions::parse_args(argn, argv, &rawdatafiles, &keyfilename, &outdir,&config)) return -1; 
   TString sumdir = outdir;//summary file directory
   sumdir+="sum/";
  // if (dbg_level & dbg_config) config->print(); 
   
   

  /* Handle keys */
   TString key = "skim"; 
   bool pass; 
   DmtpcKeys k(keyfilename, rawdatafiles, key, outdir, 0, 0, pass); 
   
   if (!pass) return -1; 
    
   /* Loop through input files */
   for (int f = 0; f < k.getNFiles(); f++)
   {
       cout << k.getFile(f) << endl;
                 
      /* Create DMTPC Database and draw out tree */
      DmtpcDataset d; 
      d.openRootFile(k.getRootDirName() + k.getFile(f)); 
        
      TTree * rawtree = k.getBaseTree(f); 
        
      //Extract run number from file 
      d.getEvent(0); 
      int runnum = d.event()->runNumber(); 
      
      //Add friends
      k.addFriends(rawtree, f); 
      
      key =  "skim"; 
      TString outfilename = k.getFile(f).ReplaceAll(".root",key+".root"); 
      DmtpcSkimDataset sd; 
      
      sd.newRootFile(outdir + outfilename); 

      /* Now ready to process this file */
      CleanSkimConfig * conf = skim::preprocess(&d,config); 
      conf->print(); 
      
      for (unsigned int nprocess = 0; nprocess < conf->getNClusterFindingAlgorithms(); nprocess++)
      {
        if(nprocess > 0) key = TString(conf->getClusterFindingAlgorithm(nprocess)); 

        return_val += cleanSkim(d,runnum,sd,true,conf,nprocess); 
      }
   }
   return return_val;
}


