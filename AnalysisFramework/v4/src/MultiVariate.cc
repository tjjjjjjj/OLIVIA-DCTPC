#include <iostream>
#include "MultiVariate.hh"
#include "TChain.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "../../../MaxCam/DmtpcSkimEvent.hh"
#include "../../../MaxCam/DmtpcMCDataset.hh"
#include "../../../MaxCam/MaxCamClusterImage.hh"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"
#include "TMVA/Factory.h"

ClassImp(MultiVariate::MultiVariateResult)

static char * opt_E = "log(E)"; 
static char * opt_range = "log(range)"; 
static char * opt_cluster_rms = "log(cluster_rms)"; 
static char * opt_maxpixel = "log(maxpixel)";
static char * opt_npixel = "log(npixel)";
static char * opt_npixel_red = "log(npixel_red)";

void add_variables(TMVA::Factory * ptr)    
{                                             
  ptr->AddVariable("E",'F');                
  ptr->AddVariable("range",'F');            
  ptr->AddVariable("cluster_rms",'F');        
  ptr->AddVariable("maxpixel",'F');         
  ptr->AddVariable("neighbors",'F');       
  ptr->AddVariable("npixel",'F');            
  ptr->AddVariable("npixel_red",'F');            
  ptr->AddVariable(opt_E,'F'); 
  ptr->AddVariable(opt_range,'F'); 
  ptr->AddVariable(opt_cluster_rms,'F'); 
  ptr->AddVariable(opt_maxpixel,'F'); 
  ptr->AddVariable(opt_npixel,'F'); 
  ptr->AddVariable(opt_npixel_red,'F'); 
}


void MultiVariate::runMVSimple( const char * job_name,
                         TFile * f,
                         TTree * signal,
                         TTree * background,
                         const char * sig_cut, 
                         const char * bg_cut, 
                         TMVA::Types::EMVA method,
                         const char * method_options
)
{

  const char * method_name = 0;
  switch(method)
  {
    case TMVA::Types::kVariable:     method_name = "Variable"; break;
    case TMVA::Types::kCuts:         method_name = "Cuts"; break;
    case TMVA::Types::kLikelihood:   method_name = "Likelihood"; break; 
    case TMVA::Types::kPDERS:        method_name = "PDERS"; break;         
    case TMVA::Types::kHMatrix:      method_name = "HMatrix"; break;         
    case TMVA::Types::kFisher:       method_name = "Fisher"; break;  
    case TMVA::Types::kKNN:          method_name = "KNN"; break;  
    case TMVA::Types::kCFMlpANN:     method_name = "CFMlpANN"; break;   
    case TMVA::Types::kTMlpANN:      method_name = "TMlpANN"; break;    
    case TMVA::Types::kBDT:          method_name = "BDT"; break;   
    case TMVA::Types::kDT:           method_name = "DT"; break;   
    case TMVA::Types::kRuleFit:      method_name = "RuleFit"; break;    
    case TMVA::Types::kSVM:          method_name = "SVM"; break;    
    case TMVA::Types::kMLP:          method_name = "MLP"; break;      
    case TMVA::Types::kBayesClassifier: 
                                     method_name = "BayesClassifier"; break;    
    case TMVA::Types::kFDA:          method_name = "FDA"; break;    
    case TMVA::Types::kCommittee:    method_name = "Committee"; break;    
    case TMVA::Types::kBoost:        method_name = "Boost"; break;    
    case TMVA::Types::kPDEFoam:      method_name = "PDEFoam"; break;    
    case TMVA::Types::kLD:           method_name = "LD"; break;    
    case TMVA::Types::kPlugins:      method_name = "Plugins"; break;    
    case TMVA::Types::kMaxMethod:    method_name = "MaxMethod"; break;    
    default: method_name = "unknown"; 
  }

  TMVA::Tools::Instance(); 	

  TMVA::Factory  * factory = new TMVA::Factory(job_name,f,"V"); 

  add_variables(factory); 

  factory->AddSignalTree(signal); 
  factory->AddBackgroundTree(background,1.); 
  
 
  //Prepare trees
  factory->PrepareTrainingAndTestTree(sig_cut, bg_cut,""); 

  //Select methods
  factory->BookMethod(method,method_name,method_options); 


  unsigned canary [512]; 
  memset(canary,0xcd,sizeof(canary)); 
  factory->TrainAllMethods(); 

  printf("Start of canary: %p\n",canary); 
  for (uint32_t i = 0; i < 512; i++)
  {
    if (canary[i] != 0xcdcdcdcd)
    { 
      printf("Corruption detected at %p!. Expected: %x, found %x\n", canary + i,0xcdcdcdcd,canary[i]);  
    }
  }
  printf("end of canary: %p\n",canary + 512); 

  factory->TestAllMethods(); 
  factory->EvaluateAllMethods(); 

  f->Close();
}

void MultiVariate::runMVMultiple(const char * name,
                         TFile * f, 
                         TTree * signal,
                         TTree ** background,
                         int nbackground, 
                         const char * sig_cut,
                         const char * bg_cut,
                         const double * background_weights,
                         const TMVA::Types::EMVA * methods,
                         const char ** method_names,
                         const char ** method_options,
                         int nmethods
                         )
{
  TMVA::Tools::Instance(); 	
  TMVA::Factory  * factory = new TMVA::Factory(TString(name),f,"V"); 

  factory->AddSignalTree(signal); 

  for (int i = 0; i < nbackground; i++)
  {
    double weight = 1; 
    if (background_weights!=NULL) weight = background_weights[i]; 
    factory->AddBackgroundTree(background[i],weight); 
  }

  //Add variables
  add_variables(factory); 

  //Prepare trees
  factory->PrepareTrainingAndTestTree(sig_cut,bg_cut,""); 

  //Select methods
  for (int i = 0 ; i < nmethods; i++)
  {
    factory->BookMethod(methods[i],TString(method_names[i]),TString(method_options[i])); 
  }

  factory->TrainAllMethods(); 
  factory->TestAllMethods(); 
  factory->EvaluateAllMethods(); 
  f->Close(); 
  f->Delete(); 
}

TTree * MultiVariate::buildMVTree(const char * name, TFile * f, const std::vector<std::string>  *files, int nburnin_thresh, const char * cam_id, const char * tree_name, double true_cluster_distance) 
{

  DmtpcSkimEvent * ev = 0; 
  TChain * skim = new TChain(tree_name);
 
  DmtpcMCDataset * mc = 0;  
  for (unsigned int i = 0; i < files->size(); i++)
  {
    skim->Add(files->at(i).c_str()); 
  } 

  skim->SetBranchAddress("event",&ev);
  skim->SetBranchStatus("*clusters",0);

  double  E, range,skewness,cluster_rms,maxpixel,npixel,npixel_red,maxpixel_over_E, cluster_rms_over_E, neighbors;
  f->cd(); 
  TTree * mv = new TTree(name,name); 
  gROOT->cd(); 
  mv->Branch("E",&E); 
  mv->Branch("range",&range); 
  mv->Branch("skewness",&skewness); 
  mv->Branch("cluster_rms",&cluster_rms); 
  mv->Branch("maxpixel",&maxpixel); 
  mv->Branch("neighbors",&neighbors);
  mv->Branch("npixel",&npixel);
  mv->Branch("npixel_red",&npixel_red); 


  int mc_ev = 0; 
  char * last_file = 0; 

  for (unsigned i = 0; i < skim->GetEntries();i++)
  {
    skim->GetEntry(i);
    if (true_cluster_distance > 0 )
    {
      if (!last_file || strcmp(last_file, skim->GetCurrentFile()->GetName()))
      {
        if (mc) mc->Delete(); 
        mc = new DmtpcMCDataset(); 
        if (last_file) free (last_file); 
        last_file = strdup(skim->GetCurrentFile()->GetName()); 
        mc->loadFile(TString(last_file).ReplaceAll("AnalysisFramework/v4/skim/","MaxCam/Simulations/v1/output/data/").ReplaceAll("skim",""), true); 
        mc_ev = 0; 
      }
      mc->getEvent(mc_ev++); 
    }

    int cam_min, cam_max; 
    if (cam_id && !strcasecmp("all",cam_id))
    {
      cam_min = 0; 
      cam_max = ev->ncamera(); 
    }
    else
    {
      cam_min = cam_id ? ev->findSerialNumber(cam_id) : 0; 
      if (cam_min < 0) continue; 

      cam_max = cam_min+1; 
    }

    for (int c = cam_min; c < cam_max; c++) 
    {
      if (ev->spark(c)) continue;  //skip sparks
    
      const MaxCamClusterImage  * clust = true_cluster_distance > 0 ?   mc->getCluster(c) : 0; 
      double x,y; 
      if (true_cluster_distance > 0 && clust->getNCluster() > 0)
      {
        clust->getXY(0,x,y); 
      }
      for (int t = 0; t < ev->ntracks(c); t++)
      {
        if (ev->edge(c,t)) continue; 
        if (nburnin_thresh > 0 && ev->nburnin(c,t) > nburnin_thresh) continue; 

        E = ev->E(c,t); 
        if (E<=0) continue; 
        range = ev->range(c,t); 
        if (range<=0) continue; 

        if (true_cluster_distance > 0)
        {
          double dist = sqrt(pow(ev->x(c,t) - x,2) + pow(ev->y(c,t) - y,2)); 
          std::cout << dist << std::endl; 
          if (dist > true_cluster_distance || dist != dist) continue; 
        }
       
        skewness = ev->skewness(c,t); 
        cluster_rms = ev->cluster_rms(c,t); 
        if (cluster_rms <= 0) continue; 
        maxpixel = ev->maxpixel(c,t); 
        if (maxpixel <= 0) continue; 
        neighbors = (double)ev->neighbors(c,t);
        npixel = (double)ev->npixel(c,t);
        npixel_red = (double)ev->npixel_red(c,t);

        f->cd(); 
        mv->Fill(); 
        gROOT->cd();
      }
    }
  }

  f->cd(); 
  mv->Write(); 
  gROOT->cd(); 

  return mv;
}

MultiVariate::MultiVariateResult::~MultiVariateResult() 
{
  if (current!=NULL) free(current); 
  reader->Delete(); 
}

MultiVariate::MultiVariateResult::MultiVariateResult() 
{

  reader = new TMVA::Reader(""); 

  reader->AddVariable("E",&E);
  reader->AddVariable("range",&range);            
  reader->AddVariable("cluster_rms",&cluster_rms);        
  reader->AddVariable("maxpixel",&maxpixel);         
  reader->AddVariable("neighbors",&neighbors);       
  reader->AddVariable("npixel",&npixel);            
  reader->AddVariable("npixel_red",&npixel_red);            
  reader->AddVariable(opt_E,&logE); 
  reader->AddVariable(opt_range,&logrange); 
  reader->AddVariable(opt_cluster_rms,&logcluster_rms); 
  reader->AddVariable(opt_maxpixel,&logmaxpixel); 
  reader->AddVariable(opt_npixel,&lognpixel); 
  reader->AddVariable(opt_npixel_red,&lognpixel_red); 
  current = NULL;

}

bool MultiVariate::MultiVariateResult::setResult(const char * job_name,
           const char * method_name, const char * path_to_weights_dir)
{
  TString path_string = TString(path_to_weights_dir) + 
                        TString("/") + TString(job_name)
                        +TString("_") + TString(method_name)
                        +TString(".weights.xml");
  if (current!=NULL) free(current); 
  current=strdup(method_name);
  return reader->BookMVA(method_name,path_string.Data()); 
}

void MultiVariate::MultiVariateResult::setrdvars(const DmtpcSkimEvent * ev, int cam, int track) 
{

  E = ev->E(cam,track); 
  range = ev->range(cam,track); 
  cluster_rms = ev->cluster_rms(cam,track);
  maxpixel = ev->maxpixel(cam,track);
  neighbors = ev->neighbors(cam,track); 
  npixel = ev->npixel(cam,track); 
  npixel_red = ev->npixel_red(cam,track); 
  logE = TMath::Log(E); 
  logcluster_rms = TMath::Log(cluster_rms); 
  logmaxpixel = TMath::Log(maxpixel); 
  logrange = TMath::Log(range); 
  lognpixel =   TMath::Log(npixel);
  lognpixel_red =   TMath::Log(npixel_red);
}

double MultiVariate::MultiVariateResult::getClassifier(const DmtpcSkimEvent * ev, int cam,
    int track, const char * method_name) 
{
  setrdvars(ev,cam,track); 
  if (method_name==NULL)
  {
    return reader->EvaluateMVA(current); 
  }
  return reader->EvaluateMVA(method_name); 
}

bool MultiVariate::MultiVariateResult::getCutsClassifier(const DmtpcSkimEvent * ev, int cam,
    int track, double sigeff,  const char * method_name) 
{
  setrdvars(ev,cam,track); 
  if (method_name==NULL)
  {
    return reader->EvaluateMVA(current,sigeff); 
  }
  return reader->EvaluateMVA(method_name,sigeff); 
}

double MultiVariate::MultiVariateResult::getSignalProbability(const DmtpcSkimEvent * ev, int cam,
    int track, double sigfrac, const char * method_name) 
{
  setrdvars(ev,cam,track); 
  if (method_name==NULL)
  {
    return reader->GetProba(current,sigfrac); 
  }
  return reader->GetProba(current,sigfrac); 
}


void MultiVariate::evaluate(MultiVariate::MultiVariateResult * r, MVTree * s, MVTree * b, const char * out_file,
                            const char * sig_cam, const char * bg_cam, int maxnburn, double true_cluster_distance)
{
  s->SetBranchStatus("*clusters*",0); 
  s->SetBranchStatus("*trigger*",0); 
  b->SetBranchStatus("*clusters*",0); 
  b->SetBranchStatus("*trigger*",0); 

  DmtpcSkimEvent * ev = 0; 

  s->SetBranchAddress("event",&ev);
  b->SetBranchAddress("event",&ev);

  std::cout << "GET_ENTRIES: "<<  s->GetEntries() << " " << b->GetEntries() << std::endl;

  TFile f(out_file,"RECREATE"); 

  //Find minimum and maximum
  double min = 1e300; 
  double max = -1e300; 
  
  double Emin = 1e300; 
  double Emax = -1e300; 

  std::vector<double> signal;
  std::vector<double> bg;
  std::vector<double> Energy;

  TH1F * sig_range_dist = new TH1F("sig_range","sig_range",100,0,200); 
  TH1F * bg_range_dist = new TH1F("bg_range","bg_range",100,0,200); 

  bool isSignal; 
  double classifier,range,E,maxpixel,cluster_rms;
  int neighbors,npixel,npixel_red, run, event, track,cam; 

  TTree * eval = new TTree("eval","eval"); 

  eval->Branch("isSignal",&isSignal); 
  eval->Branch("classifier",&classifier); 
  eval->Branch("npixel",&npixel); 
  eval->Branch("npixel_red",&npixel_red); 
  eval->Branch("cluster_rms",&cluster_rms); 
  eval->Branch("neighbors",&neighbors); 
  eval->Branch("maxpixel",&maxpixel); 
  eval->Branch("run",&run); 
  eval->Branch("event",&event); 
  eval->Branch("track",&track); 
  eval->Branch("cam",&cam); 
  eval->Branch("E",&E); 
  eval->Branch("range",&range); 

  DmtpcMCDataset mc; 
  int mc_ev = 0; 

  char * last_file = 0; 
  for (int i = 0; i < s->GetEntries(); i++)
  {
    s->GetEntry(i);         
    
    int cam_min, cam_max; 
    if (sig_cam && !strcasecmp("all",sig_cam))
    {
      cam_min = 0;  
      cam_max = ev->ncamera(); 
    }
    else
    {
      cam_min = sig_cam ? ev->findSerialNumber(sig_cam) : 0; 
      if (cam_min == -1) continue; 
      cam_max = cam_min + 1; 
    }
    for (cam = cam_min; cam < cam_max; cam++)
    {
      for (int t = 0; t < ev->ntracks(cam); t++) 
      {
        sig_range_dist->Fill(ev->range(cam,t));
        if ((maxnburn > 0 && ev->nburnin(cam,t) > maxnburn)
            || ev->edge(cam,t) || ev->range(cam,t)==0)
        {
          continue;
        }

        double c = r->getClassifier(ev,cam,t); 
        if (c < min) min = c; 
        if (c > max) max = c; 

        E = ev->E(cam,t); 
        if (E < Emin) Emin = E; 
        if (E > Emax) Emax = E; 
        Energy.push_back(E); 
        signal.push_back(c); 

        classifier = c; 
        isSignal = true; 
        npixel = ev->npixel(cam,t); 
        npixel_red = ev->npixel_red(cam,t); 
        neighbors = ev->neighbors(cam,t); 
        cluster_rms = ev->cluster_rms(cam,t); 
        maxpixel = ev->maxpixel(cam,t); 
 //       maxpixel_over_E = maxpixel/E; 
 //       cluster_rms_over_E = cluster_rms_over_E/E; 
        range = ev->range(cam,t); 
        track = t; 
        event = ev->eventNumber(); 
        run = ev->runNumber(); 
        eval->Fill(); 

      }
    }
    ev->Delete(); 
    ev = 0; 
  }
  

  std::cout << "Signal Size: "  << signal.size() << std::endl;
  for (int i = 0; i < b->GetEntries(); i++)
  {
    b->GetEntry(i);         
    int min_cam, max_cam; 
    if (bg_cam && !strcasecmp("all",bg_cam)) 
    {
      min_cam = 0; 
      max_cam = ev->ncamera(); 
    }
    else
    {
      min_cam = bg_cam ? ev->findSerialNumber(bg_cam): 0; 
      if (min_cam == -1) continue; 
      max_cam = min_cam+ 1; 
    }
    for (cam = min_cam; cam < max_cam; cam++)
    {
      for (int t = 0; t < ev->ntracks(cam); t++) 
      {
        bg_range_dist->Fill(ev->range(cam,t));
        if ((maxnburn > 0 && ev->nburnin(cam,t) > maxnburn)
            || ev->edge(cam,t) || ev->range(cam,t)==0)
        {
          continue;
        }

        double c = r->getClassifier(ev,cam,t); 
        if (c < min) min = c; 
        if (c > max) max = c; 
        bg.push_back(c); 

        classifier = c; 
        isSignal = false; 
        npixel = ev->npixel(cam,t); 
        npixel_red = ev->npixel_red(cam,t); 
        neighbors = ev->neighbors(cam,t); 
        cluster_rms = ev->cluster_rms(cam,t); 
//        cluster_rms_over_E = cluster_rms/E; 
        maxpixel = ev->maxpixel(cam,t); 
//        maxpixel_over_E = maxpixel/E; 
        range = ev->range(cam,t); 
        track = t; 
        event = ev->eventNumber(); 
        run = ev->runNumber(); 
        E = ev->E(cam,t); 
        eval->Fill(); 
      }
    }
    ev->Delete(); 
    ev = 0; 
  }


  eval->Write(); 
  std::cout << "Background Size: "  << bg.size() << std::endl;
  int nbins = 100; 

  TH1F * sighist = new TH1F("sig_value","sig_value",nbins,min,max);
  TH1F * bghist = new TH1F("bg_value","bg_value",nbins,min,max);
  TH2F * valve = new TH2F("value_v_E","value_v_E",nbins,Emin,Emax,nbins,min,max); 

  for (unsigned int i = 0; i < signal.size();i++)
  {
    sighist->Fill(signal[i]); 
    valve->Fill(Energy[i],signal[i]); 
  }

  for (unsigned int i = 0; i < bg.size();i++)
  {
    bghist->Fill(bg[i]); 
  }

  double * sig_eff = new double[nbins]; 
  double * bg_eff = new double[nbins]; 
  double * purity = new double[nbins];
  double * bins = new double[nbins];

  double cutoff; 

  for (int i = 1; i <=nbins; i++)
  {
     sig_eff[i-1] = sighist->Integral(1,i)/sighist->Integral(1,nbins); 
     bg_eff[i-1] =1 - bghist->Integral(1,i)/bghist->Integral(1,nbins);  

     if (bg_eff[i-1] == 0 && i>1 && bg_eff[i-2] > 0)
     {
       cutoff = bghist->GetBinCenter(i); 
     }
     bins[i-1] =sighist->GetBinCenter(i); 
     purity[i-1] = (sighist->Integral(i,nbins)/sighist->Integral(1,nbins)) / 
                    (sighist->Integral(i,nbins)/sighist->Integral(1,nbins) + bghist->Integral(i,nbins)/bghist->Integral(1,nbins)); 
  } 

  TGraph * sig_eff_g = new TGraph(nbins,bins,sig_eff); 
  sig_eff_g->Write("signal_efficiency"); 
  TGraph * bg_eff_g = new TGraph(nbins,bins,bg_eff); 
  bg_eff_g->Write("background_eff_cmp"); 

  TGraph * purity_g = new TGraph(nbins,bins,purity); 
  purity_g->Write("purity_g"); 
  
  sig_range_dist->Write(); 
  bg_range_dist->Write(); 

  bghist->Write();
  sighist->Write();
  valve->Write(); 

  delete sig_eff;
  delete bg_eff; 
  delete purity; 
  delete bins;

  f.Close(); 
}




