//--*-C++-*--

#include "../../MaxCam/MaxCamMC.hh"
#include "../../MaxCam/MaxCamImageTools.hh"
#include "../../MaxCam/MaxCamCluster.hh"
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
#include "TNtuple.h"
#include "TH2F.h"
#include "TPolyMarker.h"
#include "TClonesArray.h"

#include <iostream>
#include <vector>
#include <math.h>
using std::vector;

#include "../../MaxCam/DmtpcEvent.hh"

using namespace std;

TCanvas c;


// Adjust these as necessary.  Used by Image(), sum() and sumImages2()

int ccdN=0;         // ccd #: switch between 0 or 1 to access top or bottom CCDs
double cLim=-1.5;   // cold pixel threshold
double sLim=1.5;    // threshold for burned in spots from sparks
double hLim=2.5;    // hot pixel threshold
double pFac=13;     
int sprkB=1;
int blurW=4;        // size of bluring window
double blurF=0.8;   // weight for surrounding pixels during bluring


// Use for looking at single images from a run
// runNum is run number
// Img is event number
// if hot = 1, leave hot pixels in final image
// if biasB = 1, subtract bias frame
// if verbose = 1, print information about event

void Image(TString runNum, int Img, int hot=0, int biasB=1, int verbose=0)
  // If bias==0, display just image
  // If bias==1, display image with bias subtraction
  // If Img<0, display bias image instead of data; Img number ignored
  // If hot==1, show hot pixels as is
{  
  TFile* file1;
  TH2F* ImgSum;
  TString fname;
  int ntracks[2];
  bool sprk[2];
  TTree* skim;
  TTree* bias;

  if (verbose==1)
    {
      fname = "skim/dmtpc_run" + runNum + "skim.root"; // open skim file
      file1 = new TFile(fname,"READ");
      skim = (TTree*) file1->Get("skim");
      skim->SetBranchAddress("ntracks",ntracks);
      skim->SetBranchAddress("spark",sprk);
      skim->GetEntry(Img);
      cout << "ntracks: " << ntracks[ccdN] << "  spark: " << sprk[ccdN] << endl;
    }
  if (biasB==0)
    {
      fname = "skim/dmtpc_run" + runNum + "skim.root"; // open skim file
      file1 = new TFile(fname,"READ");
      DmtpcEvent *dm1;
      skim = (TTree*) file1->Get("skim");
      skim->SetBranchAddress("event",&dm1);
      skim->GetEntry(Img);
      ImgSum = (TH2F*)dm1->ccdData(ccdN);
    }
  elseif (Img<0)
    {
      fname = "skim/dmtpc_run" + runNum + "bias.root";
      file1 = new TFile(fname,"READ");
      bias = (TTree*) file1->Get("bias");
      bias->SetBranchAddress("biasframe",&ImgSum);
      bias->GetEntry(ccdN);
    }
  else
    {
      TH2F* bImg;

      fname = "skim/dmtpc_run" + runNum + "skim.root"; // open skim file
      file1 = new TFile(fname,"READ");
      DmtpcEvent *dm1;
      skim = (TTree*) file1->Get("skim");
      skim->SetBranchAddress("event",&dm1);
      skim->GetEntry(Img);
      ImgSum = (TH2F*)dm1->ccdData(ccdN);

      fname = "skim/dmtpc_run" + runNum + "bias.root";
      file1 = new TFile(fname,"READ");
      bias = (TTree*) file1->Get("bias");
      bias->SetBranchAddress("biasframe",&bImg);
      bias->GetEntry(ccdN);
      ImgSum->Add(ImgSum,bImg,1,-1);
    }
    
  TH2F* hotImg = findHot(ImgSum,hLim);

  c.Clear();
  c.Divide(2,1);
  c.cd(1);
  hotImg->Draw("colz");

  c.cd(2);
  if (hot==0)
    {
      for(int i=1; i<=ImgSum->GetNbinsX(); i++)
	{
	  for(int j=1; j<=ImgSum->GetNbinsY(); j++)
	    {
	      if (hotImg->GetBinContent(i,j)>=1)
		ImgSum->SetBinContent(i,j,0.0);
	    }
	}
      //   hotImg->Add(hotImg,1);
      }
  ImgSum->Draw("colz");
}



// Use to create integrated image of frames from a run.
// runNum is the string for the run number
// sImg is the event number to start integration
// eImg is the event number to end on
// if blurB = 1, blur image after integration
// if biasB = 1, carry out bias subtraction

TH2F* sum(TString runNum, int sImg, int eImg, int blurB=1, int biasB=1)
{  
  TFile* file1;
  TFile* file2;
  DmtpcEvent* dm1;
  TH2F* ImgSum;
  TH2F* bImg;
 
  fname = "skim/dmtpc_run" + runNum + "bias.root";
  file2 = new TFile(fname,"READ");
  bias = (TTree*) file2->Get("bias");
  bias->SetBranchAddress("biasframe",&bImg);
  bias->GetEntry(ccdN);
  //file2->Close();

  TString fname = "skim/dmtpc_run" + runNum + "skim.root"; // open skim file
  file1 = new TFile(fname,"READ");
  skim = (TTree*) file1->Get("skim");
  skim->SetBranchAddress("event",&dm1);
  skim->GetEntry(sImg);
  ImgSum = dm1->ccdData(ccdN);
  ImgSum->Scale(0.000000001);
  ImgSum = sumImages2(ImgSum,bImg,runNum,sImg,eImg,blurB,biasB);
  
  return ImgSum;
}



// Used by TH2F* sum()
TH2F* sumImages2(TH2F* ImgSum, TH2F* bImg, TString runNum, int sImg, int eImg, int blurB=1, int biasB=1)//
{
  TFile* file1;
  TFile* file2;
  DmtpcEvent* dm1;
  TH2F* tempImg;
  TH2F* tempImg2;
  TH2F* tempImg3;
  TH2F* wormImg;
  TH2F* sprkFill;
  int ntracks[2];
  bool sprk[2];
  TTree* skim;
  TTree* sksort;
  TObjArray* TBdeleted = new TObjArray();

  int nbinsx = ImgSum->GetNbinsX();
  int nbinsy = ImgSum->GetNbinsY();

  TH1D* pixCount = new TH1D("","",eImg-sImg,sImg,eImg);

  // Get params from skim file
  TString fname = "skim/dmtpc_run" + runNum + "skim.root"; // open skim file
  file1 = new TFile(fname,"READ");
  skim = (TTree*) file1->Get("skim");
  skim->SetBranchAddress("event",&dm1);
  skim->SetBranchAddress("ntracks",ntracks);
  skim->SetBranchAddress("spark",sprk);
  skim->GetEntry(sImg);
  skim->SetBranchAddress("clusters",&TBdeleted);
  //ImgSum = dm1->ccdData(ccdN);

  double alpha[2][15], frecoil[2][15], worm[2][15], other[2][15];

  // Get params from sksort file
  fname = "sort/dmtpc_run" + runNum + "sksort.root";
  file2 = new TFile(fname,"READ");
  sksort = (TTree*) file2->Get("sksort");
  sksort->SetBranchAddress("alpha",alpha);
  sksort->SetBranchAddress("frecoil",frecoil);
  sksort->SetBranchAddress("worm",worm);
  sksort->SetBranchAddress("other",other);
  sksort->GetEntry(sImg);

  // Initialize 2-D histograms
  TH2F* hotImg = findHot(ImgSum);
  sprkFill = (TH2F*) ImgSum->Clone("sprkFill");
  sprkFill->Scale(0.000000001);
 
  //cout<<skim->GetEntries()<<endl;
  double counter=0;     // number of frames used (excludes spark/alpha frames)
  double firstBlank=0;  // number frames used in integrated image for filling in burned in spots
  int track_check;      // boolean for whether a track is found by AnalysisFramework
  int dImg = (eImg-sImg)/10;
  double yieldSumT=0  ; // total pixel yield sum
  double yieldSum=0;

  double yS=0;
  double byS = bImg->Integral();

  // Loop over every frame and add non-spark/alpha frames
  for(int i=sImg; i<eImg; i++)
    {
      skim->GetEntry(i);
      sksort->GetEntry(i);

      TBdeleted->Delete(); 

      if (!(sprk[ccdN]))
	{
	  track_check=0;
	  for(int j=0; j<ntracks[ccdN]; j++)
	    {
	      if (alpha[0][j]==1) track_check=1;
	      if (frecoil[0][j]==1) track_check=1;
	      // if (worm[0][j]==1) track_check=1;  // worms cut using hot pixel threshold
	      // if (other[0][j]==1) track_check=1;
	    }
	  if (track_check==1) continue;

	  tempImg = dm1->ccdData(ccdN);
	  counter++;

	  yS = tempImg->Integral();
	  pixCount->SetBinContent(i,yS-byS);
	  yieldSumT += yS;
	  ImgSum->Add(ImgSum,tempImg,1,1);
	  if (i%dImg == 0)
	    {
	      tempImg = findHot2(tempImg,5);
	      hotImg->Multiply(tempImg,hotImg,1,1);
	    }
	}

      // If spark encountered, drop frame
      // If first spark, save integrated image to correct burned in spots
       else 
 	{
	  if (firstBlank==0)
	    {
	      tempImg3 = (TH2F*) ImgSum->Clone("tempImg3");
	      firstBlank = counter;
	      yieldSum = yieldSumT;
	    }
	}
    }
	
  if (firstBlank==0)
    {
      tempImg3 = (TH2F*) ImgSum->Clone("tempImg3");
      firstBlank = counter;
      yieldSum = yieldSumT;
    }

  c.Clear();
  c.Divide(2,2);

  if (biasB) ImgSum->Add(ImgSum,bImg,1,-counter);
  double yieldSumT2 = ImgSum->Integral();
  ImgSum->Scale(1/counter);

  tempImg = findHot(ImgSum,3);
  tempImg2 = findCold(ImgSum,-3);
  tempImg2->Add(tempImg,tempImg2,1,1);
  AddTH2(tempImg,-1,1);
  ImgSum->Multiply(ImgSum,tempImg,1,1);

  c.cd(2);
  //hotImg->Draw("colz");
  pixCount->Draw();

  c.cd(3);
  wormImg=findHot(bImg,.11);
  wormImg->Add(wormImg,hotImg,1,-1);
  wormImg->Draw("colz");

  hotImg->Scale(0.0);

  c.cd(4);

  tempImg3->Add(tempImg3,bImg,1,-firstBlank);
  double yieldSum2 = tempImg3->Integral();
  tempImg3->Scale(1/firstBlank);
  sprkFill->Multiply(sprkFill,tempImg3,1,1);
  sprkFill->Draw("colz");

  c.cd(1);

  if (blurB) tempImg = findCold(ImgSum,cLim);
  tempImg2->Add(tempImg,tempImg2,1,1);

  if (sprkB)
    {
      tempImg3 = findHot(ImgSum,sLim);
      AddTH2(tempImg3,-1,1);
      ImgSum->Multiply(ImgSum,tempImg3,1,1);
      ImgSum->Add(ImgSum,sprkFill,1,firstBlank/counter*yieldSumT2/yieldSum2);
    }

  if (blurB)
    {
      tempImg = findHot(ImgSum,hLim);
      tempImg->Add(tempImg,tempImg2,1,1);
      AddTH2(tempImg,1,0.01);
      ImgSum = blur(ImgSum,tempImg,blurW,blurF);
    }

  ImgSum->Draw("colz");  

  //cout << counter << " " << firstBlank << endl;

  return ImgSum;
}


void sumImages(TH2F* ImgSum, TString runNum)
{
  TFile* file1;
  DmtpcEvent* dm1;
  TH2F* tempImg;
  int ntracks;

  TString fname = "skim/dmtpc_run" + runNum + "skim.root"; // open skim file
  file1 = new TFile(fname,"READ");
  skim->SetBranchAddress("event",&dm1);
  skim->SetBranchAddress("ntracks",&ntracks);
  skim->GetEntry(0);
  //ImgSum = dm1->ccdData(ccdN);
 
  cout<<skim->GetEntries()<<endl;
  for(int i=1; i<skim->GetEntries(); i++) // Sum skimmed images
    {
      //cout<<"loop "<<endl;
      skim->GetEntry(i);
      if (ntracks>=0)
	{
	  tempImg = dm1->ccdData(ccdN);
	  ImgSum->Add(ImgSum,tempImg,1,1);
	}
    }

}



// Use to monitor the activity of up to three different pixels over the course of a run.
// if biasB=1, bias frame is subtracted
void pixels(TString runNum, int x1, int y1, int x2, int y2, int x3, int y3, int biasB=1)
{  
  TFile* file1;
  TFile* file2;
  DmtpcEvent* dm1;
  TH2F* Img;
  TH2F* bImg;

  int ntracks[2];
  bool sprk[2];

  fname = "skim/dmtpc_run" + runNum + "bias.root";
  file2 = new TFile(fname,"READ");
  bias = (TTree*) file2->Get("bias");
  bias->SetBranchAddress("biasframe",&bImg);
  bias->GetEntry(0);

  TString fname = "skim/dmtpc_run" + runNum + "skim.root"; // open skim file
  file1 = new TFile(fname,"READ");
  skim = (TTree*) file1->Get("skim");
  skim->SetBranchAddress("event",&dm1);
  skim->SetBranchAddress("ntracks",ntracks);
  skim->SetBranchAddress("spark",sprk);
  skim->GetEntry(0);
  Img = dm1->ccdData(ccdN);

  int numE = skim->GetEntries();
  TH1D* h1 = new TH1D("h1","",numE,0,numE);
  TH1D* h2 = new TH1D("h2","",numE,0,numE);
  TH1D* h3 = new TH1D("h3","",numE,0,numE);
  TH1D* h4 = new TH1D("h4","",numE,0,numE);
 
  double nbinsx = Img->GetNbinsX();
  double nbinsy = Img->GetNbinsY();
  double bavg = bImg->Integral()/nbinsx/nbinsy;

  double b1 = bImg->GetBinContent(x1,y1);
  double b2 = bImg->GetBinContent(x2,y2);
  double b3 = bImg->GetBinContent(x3,y3);

  double alpha[2][15], frecoil[2][15], worm[2][15], other[2][15];

  // Get params from sksort file
  fname = "sort/dmtpc_run" + runNum + "sksort.root";
  file2 = new TFile(fname,"READ");
  sksort = (TTree*) file2->Get("sksort");
  sksort->SetBranchAddress("alpha",alpha);
  sksort->SetBranchAddress("frecoil",frecoil);
  sksort->SetBranchAddress("worm",worm);
  sksort->SetBranchAddress("other",other);
  int track_check;

  double sum;
  double sqsum;
  int counter=0;
  for (int i=0; i<numE; i++)
    {
 
      skim->GetEntry(i);
      sksort->GetEntry(i);
      if (sprk[ccdN]) continue;

	  track_check=0;
	  for(int j=0; j<ntracks[ccdN]; j++)
	    {
	      if (alpha[0][j]==1) track_check=1;
	      if (frecoil[0][j]==1) track_check=1;
	      //if (worm[0][j]==1) track_check=1;
	    }
	  if (track_check==1) continue;

	  counter++;

      Img = dm1->ccdData(ccdN);
      h1->SetBinContent(i,(Img->GetBinContent(x1,y1)-b1));
      h2->SetBinContent(i,(Img->GetBinContent(x2,y2)-b2));
      h3->SetBinContent(i,(Img->GetBinContent(x3,y3)-b3));
      h4->SetBinContent(i,(Img->Integral()/nbinsx/nbinsy)-bavg);

      sum+=(Img->GetBinContent(x1,y1)-b1);
      sqsum+=pow((Img->GetBinContent(x1,y1)-b1),2);
    }

  c.Clear();
  c.Divide(2,2);
  c.cd(1);
  h1->Draw();
  c.cd(2);
  h2->Draw();
  c.cd(3);
  h3->Draw();
  c.cd(4);
  h4->Draw();
  cout << h1->Integral() <<" "<< h2->Integral() <<" "<< h3->Integral() <<endl;

  cout << sqsum/counter-pow(sum/counter,2) << endl;

}


///////////////



// Identifies pixels in Img that are below threshold (see thresh in code below)

TH2F* findCold(TH2F* Img, double sigma=-1.5)
{
  TH2F* hotImg = (TH2F*) Img->Clone("hotImg");
  hotImg->Scale(0.0);
  int nbinsx = Img->GetNbinsX();
  int nbinsy = Img->GetNbinsY();

  double avg = Img->Integral()/nbinsx/nbinsy;
  double rms = GetRMS2(Img);
  double thresh = sigma*rms+avg; 

  for(int i=1; i<=nbinsx; i++)
    {
      for(int j=1; j<=nbinsy; j++)
	{
	  if (Img->GetBinContent(i,j)<thresh)
	    {hotImg->SetBinContent(i,j,1.0);}
	}
    }

  return hotImg;
}



// Identifies pixels in Img that are above threshold 


TH2F* findHot(TH2F* Img, double sigma=1.5)
{
  TH2F* hotImg = (TH2F*) Img->Clone("hotImg");
  hotImg->Scale(0.0);
  int nbinsx = Img->GetNbinsX();
  int nbinsy = Img->GetNbinsY();

  double avg = Img->Integral()/nbinsx/nbinsy;
  double rms = GetRMS2(Img);
  double thresh = sigma*rms+avg;

  for(int i=1; i<=nbinsx; i++)
    {
      for(int j=1; j<=nbinsy; j++)
	{
	  if (Img->GetBinContent(i,j)>thresh)
	    {hotImg->SetBinContent(i,j,1.0);}
	}
    }

  return hotImg;
}



TH2F* findHot2(TH2F* Img, double weight=0.1)
{
  TH2F* hotImg = (TH2F*) Img->Clone("hotImg");
  hotImg->Scale(0.0);
  int nbinsx = Img->GetNbinsX();
  int nbinsy = Img->GetNbinsY();

  TH1D* profile1 = (TH1D*) Img->ProjectionY("profile1");
  double avg = profile1->Integral()/nbinsx/nbinsx;
  double max = Img->GetMaximum();
  double thresh = weight*max+(1-weight)*avg;

  for(int i=1; i<=nbinsx; i++)
    {
      for(int j=1; j<=nbinsy; j++)
	{
	  if (Img->GetBinContent(i,j)>thresh)
	    {hotImg->SetBinContent(i,j,1.0);}
	  //else
	  //  hotImg->SetBinContent(i,j,0.0);
	}
    }

  return hotImg;
}



// Scale Img by Iw and add constant w to all pixels

void AddTH2(TH2F* Img, double Iw, double w)
{
  int a = w;
  int nbinsx = Img->GetNbinsX();
  int nbinsy = Img->GetNbinsY();

  double val;
  for(int i=1; i<=nbinsx; i++)
    {
      for(int j=1; j<=nbinsy; j++)
	{
	  val = Img->GetBinContent(i,j);
	  Img->SetBinContent(i,j,(val*Iw+w));
	}
    }
}



// Scale Img1 by c1 and Img2 by c2
// Divide pixel values of Img1 by those of Img2
// Yes, c2 is redundant

TH2F*  DivTH2(TH2F* Img1, TH2F* Img2, double c1, double c2)
{
  int nbinsx = Img1->GetNbinsX();
  int nbinsy = Img1->GetNbinsY();
  TH2F* Img = (TH2F*) Img1->Clone("Img");
  Img->Scale(0.0);

  double val1,val2;
  for(int i=1; i<=nbinsx; i++)
    {
      for(int j=1; j<=nbinsy; j++)
	{
	  val1 = Img1->GetBinContent(i,j);
	  val2 = Img2->GetBinContent(i,j);
	  if !((val2*c2)==0) //(val1*c1)/(val2*c2)>0)
	      Img->SetBinContent(i,j,(val1*c1)/(val2*c2));	    
	  else
	      Img->SetBinContent(i,j,0.5);
	}
    }
  return Img;
}



// Invert bin values

TH2F*  Invert(TH2F* Img1)
{
  int nbinsx = Img1->GetNbinsX();
  int nbinsy = Img1->GetNbinsY();
  TH2F* Img = (TH2F*) Img1->Clone("Img");
  Img->Scale(0.0);

  double val;
  for(int i=1; i<=nbinsx; i++)
    {
      for(int j=1; j<=nbinsy; j++)
	{
	  val = Img1->GetBinContent(i,j);
	  if ((val)!=0) 
	      Img->SetBinContent(i,j,(1.0/val));	    
	  else
	      Img->SetBinContent(i,j,1.0);
	}
    }
  return Img;
}



// Get variance of pixels in Img

double GetRMS2(TH2F* Img)
{
  int nbinsx = Img->GetNbinsX();
  int nbinsy = Img->GetNbinsY();

  double val;
  double sqSum;
  double sum;
  for(int i=1; i<=nbinsx; i++)
    {
      for(int j=1; j<=nbinsy; j++)
	{
	  val = Img->GetBinContent(i,j);
	  sqSum+=(val*val);
	  sum+=val;
	}
    }
  double total=nbinsx*nbinsy;
  sqSum=sqSum/total-sum*sum/total/total;
  return sqrt(fabs(sqSum));
}
   


// Blur pixels in image while excluding pixels identified to be out of range
// pixels of value 1 in hotImg should correspond to pixels beyond threshold in image
// 2*blurn+1 is the size of the blurring window
// blurfrac is the weight given to neighboring pixels during blurring

TH2F* blur(TH2F* image, TH2F* hotImg, int blurn, double blurfrac)
{
   TH2F* copy = (TH2F*)image->Clone("copy");
   int nx = copy->GetNbinsX();
   int ny = copy->GetNbinsY();

   for(int i=1; i<=nx; i++)
   {
      for(int j=1; j<=ny; j++)
      {	 
	 double sum=0;
	 int nin=0;

	 for(int k=i-blurn; k<= i+blurn; k++)
	 {
	    for(int l=j-blurn; l<=j+blurn; l++)
	    {
	      if(k>0 && l>0 && k<=nx && l<=ny && hotImg->GetBinContent(k,l)<1)
	       {
		  nin++;
		  if(k==i && l==j) sum+=image->GetBinContent(k,l);
		  else sum+=image->GetBinContent(k,l)*blurfrac;
	       } 

	    }
	 }
	 sum = sum/(1+(nin-1)*blurfrac);

	 copy->SetBinContent(i,j,sum);

      }
   }
   
   return copy;

}



// Produces image showing theoretical lens effect

TH2F* lensEffect(TH2F* img)
{

  TH2F* lImg = (TH2F*) img->Clone("lImg");
  lImg->Scale(0);

  int nx = lImg->GetNbinsX();
  int ny = lImg->GetNbinsY();
  
  // all cm's
  double vertDis = 40.0; 
  double vertDisSq = pow(vertDis,2);
  double dFac = 15.0/256.0;
  double effFac = 1.0;

  double distance;
  double cosDis;
  double x;
  double y;
  
  for(int i=1; i<=nx; i++)
    {
      for(int j=1; j<=ny; j++)
	{
	  x = double(i-nx/2);
	  y = double(j-ny/2);
	  distance = dFac*sqrt(double(x*x+y*y));
	  cosDis = vertDis/sqrt(vertDisSq + distance*distance);
	  lImg->SetBinContent(i,j,(effFac*pow(cosDis,4)));
	}
    }

  cout << cosDis;
  return lImg;
}



// Use to subtract one image from another.
// If images are shifted (not rotated) relative to each other, use dx and dy
// Currently only implemented for POSITIVE dx and dy
// dx is horizontal shift (in pixels) applied to Img2
// dx is vertical shift applied to Img2
// Img2 is subtracted from Img1
// The result is truncated to show just the region overlapped by Img1 and Img2

TH2F* contrast(TH2F* Img1, TH2F* Img2, int dx, int dy)
{
  
  int nbinsx = Img1->GetNbinsX();
  int nbinsy = Img1->GetNbinsY();

  int nx = abs(dx-nbinsx);
  int ny = abs(dy-nbinsy);
  TH2F* result = new TH2F("result","",nx,0,nx,ny,0,ny);

  double val1,val2;
  for(int i=1; i<=nx; i++)
    {
      for(int j=1; j<=ny; j++)
	{
	  val1 = Img1->GetBinContent(dx+i,dy+j);
	  val2 = Img2->GetBinContent(i,j);
	  result->SetBinContent(i,j,(val1-val2));
	}
    }

  return result;
}


// Same as contrast() except the result is not truncated
TH2F* contrast2(TH2F* Img1, TH2F* Img2, int dx, int dy)
{
  int nbinsx = Img1->GetNbinsX();
  int nbinsy = Img1->GetNbinsY();

  int nx = abs(dx-nbinsx);
  int ny = abs(dy-nbinsy);
  TH2F* result = new TH2F("result","",nbinsx,0,nbinsx,nbinsy,0,nbinsy);
  AddTH2(result,0.0,1.0);

  double val1,val2;
  for(int i=1; i<=nx; i++)
    {
      for(int j=1; j<=ny; j++)
	{
	  val1 = Img1->GetBinContent(dx+i,dy+j);
	  val2 = Img2->GetBinContent(i,j);
	  result->SetBinContent(dx+i,dy+j,(val1-val2));
	}
    }

  return result;
}




TH2F* Offset(TH2F* Img1, int dx, int dy)
{
  int nbinsx = Img1->GetNbinsX();
  int nbinsy = Img1->GetNbinsY();

  int nx = abs(dx-nbinsx);
  int ny = abs(dy-nbinsy);
  TH2F* result = new TH2F("result","",nbinsx,0,nbinsx,nbinsy,0,nbinsy);
  AddTH2(result,0.0,1.0);

  double val1;
  for(int i=1; i<=nbinsx; i++)
    {
      for(int j=1; j<=nbinsy; j++)
	{
	  if (dx+i>0 && dx+i<nbinsx && dy+j>0 && dy+j<nbinsy)
	    {
	      val1 = Img1->GetBinContent(i,j);
	      result->SetBinContent(dx+i,dy+j,val1);
	    }
	}
    }

  return result;
}



// Change pixel values to their absolute values

TH2F* AbsTH2(TH2F* IImg)
{
  TH2F* Img = (TH2F*) IImg->Clone("Img");  

  int nbinsx = Img->GetNbinsX();
  int nbinsy = Img->GetNbinsY();

  double val;
  for(int i=1; i<=nbinsx; i++)
    {
      for(int j=1; j<=nbinsy; j++)
	{
	  val = Img->GetBinContent(i,j);
	  Img->SetBinContent(i,j,fabs(val));
	}
    }

  return Img;
}



// Compares every pair of bin values for IImg1 and IImg2.
// Returns image where every bin value is the lower of each pair.
// cc1 and cc2 scale the images before extraction - not currently implemented.

TH2F* extract(TH2F* IImg1, TH2F* IImg2) //, double cc1=1.0, double cc2=1.0)
{
  TH2F* ssubImg = contrast(IImg1,IImg2,0,0);
  TH2F* absImg = AbsTH2(ssubImg);

  IImg2->Scale(-1.0);
  TH2F* ssmImg = contrast(IImg1, IImg2,0,0);
  IImg2->Scale(-1.0);
  
  TH2F* finImg = contrast(ssmImg, absImg,0,0);

  finImg->Scale(0.5);
  finImg->Draw("colz");

  return finImg;
}



// rotate image by 90 degrees

TH2F* rotate(TH2F* IImg)
{
  TH2F* Img = (TH2F*) IImg->Clone("Img");  

  int nbinsx = Img->GetNbinsX();
  int nbinsy = Img->GetNbinsY();

  double val;
  for(int i=1; i<=nbinsx; i++)
    {
      for(int j=1; j<=nbinsy; j++)
	{
	  val = IImg->GetBinContent(i,j);
	  Img->SetBinContent(nbinsy-j,i,val);
	}
    }

  return Img;
}
