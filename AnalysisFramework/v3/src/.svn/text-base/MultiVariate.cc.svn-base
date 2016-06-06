#include <iostream>
#include "MultiVariate.hh"
#include "TChain.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "../../../MaxCam/DmtpcSkimEvent.hh"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"
#include "TMVA/Factory.h"

ClassImp(MultiVariate::MultiVariateResult)
ClassImp(MultiVariate::MultiVariateEvaluation)

static char * opt_E = "log(E)"; 
static char * opt_range = "log(range)"; 
static char * opt_cluster_rms = "log(cluster_rms)"; 
static char * opt_maxpixel = "log(maxpixel)"; 
static char * opt_npixel = "10*log(npixel)";
static char * opt_neighbors = "neighbors"; 

void add_variables(TMVA::Factory * ptr, bool lda_optimized = true)    
{                                             
  if (!lda_optimized)                         
  {                                           
    ptr->AddVariable("E",'F');                
    ptr->AddVariable("range",'F');            
    ptr->AddVariable("cluster_rms",'F');        
    ptr->AddVariable("maxpixel",'F');         
    ptr->AddVariable("neighbors",'I');       
    ptr->AddVariable("npixel",'I');            
  }   
  else 
  { 
    ptr->AddVariable(opt_E,'F'); 
    ptr->AddVariable(opt_range,'F'); 
    ptr->AddVariable(opt_cluster_rms,'F'); 
    ptr->AddVariable(opt_maxpixel,'F'); 
    ptr->AddVariable(opt_neighbors,'I');
    ptr->AddVariable(opt_npixel,'I'); 
  }
}


void MultiVariate::runMVSimple( const char * job_name,
                         const char * outputfile,
                         TTree * signal,
                         TTree * background,
                         TMVA::Types::EMVA method,
                         const char * method_options,
                         bool optimize_for_lda
)
{

  char * method_name;
  switch(method)
  {
    case TMVA::Types::kVariable:     method_name = "Variable"; break;
    case TMVA::Types::kCuts:         method_name = "Cuts"; break;
    case TMVA::Types::kSeedDistance: method_name = "SeedDistance"; break;         
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
  TFile* f = new TFile(TString(outputfile),"RECREATE");
  TMVA::Factory  * factory = new TMVA::Factory(job_name,f); 

  factory->AddSignalTree(signal); 
  factory->AddBackgroundTree(background,1.); 
  
  add_variables(factory,optimize_for_lda); 
 
  //Prepare trees
  factory->PrepareTrainingAndTestTree("",""); 

  //Select methods
 
  factory->BookMethod(method,method_name,method_options); 

  factory->TrainAllMethods(); 
  factory->TestAllMethods(); 
  factory->EvaluateAllMethods(); 


  f->Close();
}

void MultiVariate::runMVMultiple(const char * name,
                         const char * outputfile, 
                         TTree * signal,
                         TTree ** background,
                         int nbackground, 
                         const double * background_weights,
                         const TMVA::Types::EMVA * methods,
                         const char ** method_names,
                         const char ** method_options,
                         int nmethods,
                         bool optimize_for_lda
                         )
{
  TMVA::Tools::Instance(); 	
  TFile * f = new TFile(outputfile,"RECREATE");
  TMVA::Factory  * factory = new TMVA::Factory(TString(name),f,"V"); 

  factory->AddSignalTree(signal); 

  for (int i = 0; i < nbackground; i++)
  {
    double weight = 1; 
    if (background_weights!=NULL) weight = background_weights[i]; 
    factory->AddBackgroundTree(background[i],weight); 
  }

  //Add variables
  add_variables(factory,optimize_for_lda); 

  //Prepare trees
  factory->PrepareTrainingAndTestTree("",""); 

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


TTree * MultiVariate::buildMVTree(const char * name, const std::vector<std::string>  *files, int nburnin_thresh,  uint32_t cams)
{

  DmtpcSkimEvent * ev = 0; 
  TChain * skim = new TChain("skim");
 
  for (unsigned int i = 0; i < files->size(); i++)
  {
    skim->Add(files->at(i).c_str()); 
  } 

  skim->SetBranchAddress("event",&ev);
  
  double  E, range,skewness,cluster_rms,maxpixel;
  int neighbors,npixel; 


  skim->SetBranchStatus("*clusters",0);
  skim->SetBranchStatus("*trigger_groups",0);

  TTree * mv = new TTree(name,name); 
  mv->Branch("E",&E,"E/D"); 
  mv->Branch("range",&range,"range/D"); 
  mv->Branch("skewness",&skewness,"skewness/D"); 
  mv->Branch("cluster_rms",&cluster_rms,"cluster_rms/D"); 
  mv->Branch("maxpixel",&maxpixel,"maxpixel/D"); 
  mv->Branch("neighbors",&neighbors,"neighbors/I");
  mv->Branch("npixel",&npixel,"npixel/I");

  for (int i = 0; i < skim->GetEntries();i++)
  {
    skim->GetEntry(i);
    for (int c = 0; c < ev->ncamera();c++)
    {
      if (!((c+1)&cams)) continue;  
      if (ev->spark(c)) continue;  //skip sparks

      for (int t = 0; t < ev->ntracks(c); t++)
      {
        if (ev->edge(c,t)) continue; 
        if (nburnin_thresh > 0 && ev->nburnin(c,t) > nburnin_thresh) continue; 

        E = ev->E(c,t); 
        if (E<=0) continue; 
        range = ev->range(c,t); 
        if (range<=0) continue; 
        skewness = ev->skewness(c,t); 
        cluster_rms = ev->cluster_rms(c,t); 
        maxpixel = ev->maxpixel(c,t); 
        neighbors = ev->neighbors(c,t);
        npixel = ev->npixel(c,t);
        mv->Fill(); 

      }
    }
  }

  return mv;
}

MultiVariate::MultiVariateResult::~MultiVariateResult() 
{
  if (current!=NULL) free(current); 
  reader->Delete(); 
}

MultiVariate::MultiVariateResult::MultiVariateResult( bool lda_optimized ) 
{

  reader = new TMVA::Reader(""); 

  optimized = lda_optimized; 
  if (!lda_optimized)                         
  {                                           
    reader->AddVariable("E",&E);
    reader->AddVariable("range",&range);            
    reader->AddVariable("cluster_rms",&cluster_rms);        
    reader->AddVariable("maxpixel",&maxpixel);         
   reader->AddVariable("neighbors",&neighbors);       
    reader->AddVariable("npixel",&npixel);            
  }   
  else 
  { 
    reader->AddVariable(opt_E,&E); 
    reader->AddVariable(opt_range,&range); 
    reader->AddVariable(opt_cluster_rms,&cluster_rms); 
    reader->AddVariable(opt_maxpixel,&maxpixel); 
    reader->AddVariable(opt_neighbors,&neighbors);
    reader->AddVariable(opt_npixel,&npixel); 
  }
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

  if (optimized)
  {
    E = TMath::Log(E); 
    cluster_rms = TMath::Log(cluster_rms); 
    maxpixel = TMath::Log(maxpixel); 
    range = TMath::Log(range); 
    npixel =  (int) (10 * TMath::Log(npixel)); //Multiply by 10 to separate integers 
  }
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


const MultiVariate::MultiVariateEvaluation 
MultiVariate::evaluate(MultiVariate::MultiVariateResult * r, MVTree * s, 
                       MVTree * b, int sig_cam, int bg_cam, int maxnburn)
{
  s->SetBranchStatus("*clusters*",0); 
  s->SetBranchStatus("*trigger*",0); 
  b->SetBranchStatus("*clusters*",0); 
  b->SetBranchStatus("*trigger*",0); 

  DmtpcSkimEvent * ev = 0; 

  s->SetBranchAddress("event",&ev);
  b->SetBranchAddress("event",&ev);

  std::cout << "GET_ENTRIES: "<<  s->GetEntries() << " " << b->GetEntries() << std::endl;

  //Find minimum and maximum
  double min = 1e300; 
  double max = -1e300; 
  

  std::vector<double> signal;
  std::vector<double> bg;

  int cam = sig_cam;
  for (int i = 0; i < s->GetEntries(); i++)
  {
    s->GetEntry(i);         
    for (int t = 0; t < ev->ntracks(cam); t++) 
    {
      if ((maxnburn > 0 && ev->nburnin(cam,t) > maxnburn)
          || ev->edge(cam,t) || ev->range(cam,t)==0)
      {
        continue;
      }

      double c = r->getClassifier(ev,cam,t); 
      if (c < min) min = c; 
      if (c > max) max = c; 

      signal.push_back(c); 
    }
    ev->Delete(); 
    ev = 0; 
  }

  std::cout << "Signal Size: "  << signal.size() << std::endl;
  cam = bg_cam; 
  for (int i = 0; i < b->GetEntries(); i++)
  {
    b->GetEntry(i);         
    for (int t = 0; t < ev->ntracks(cam); t++) 
    {
      if ((maxnburn > 0 && ev->nburnin(cam,t) > maxnburn)
          || ev->edge(cam,t) || ev->range(cam,t)==0)
      {
        continue;
      }

      double c = r->getClassifier(ev,cam,t); 
      if (c < min) min = c; 
      if (c > max) max = c; 
      bg.push_back(c); 
    }
    ev->Delete(); 
    ev = 0; 
  }


  std::cout << "Background Size: "  << bg.size() << std::endl;
  int nbins = 100; 

  TH1F * sighist = new TH1F("sighist","sighist",nbins,min,max);
  TH1F * bghist = new TH1F("bghist","bghist",nbins,min,max);

  for (unsigned int i = 0; i < signal.size();i++)
  {
    sighist->Fill(signal[i]); 
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
  TGraph * bg_eff_g = new TGraph(nbins,bins,bg_eff); 
  TGraph * purity_g = new TGraph(nbins,bins,purity); 

  bghist->Delete();
  sighist->Delete();

  delete sig_eff;
  delete bg_eff; 
  delete purity; 
  delete bins;

  return MultiVariateEvaluation(purity_g,sig_eff_g,bg_eff_g,cutoff); 
}


MultiVariate::MultiVariateEvaluation::~MultiVariateEvaluation()
{
  _purity->Delete(); 
  _efficiency->Delete(); 
  _background->Delete(); 
}

MultiVariate::MultiVariateEvaluation::MultiVariateEvaluation(TGraph * purity, TGraph * efficiency, TGraph * background, double cutoff)
{
  _purity = purity; 
  _efficiency = efficiency; 
  _background = background; 
  _cutoff = cutoff; 

}


