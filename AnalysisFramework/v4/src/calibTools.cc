#include "calibTools.hh"
#include "TH2.h"
#include "TCanvas.h"
#include "TString.h"
#include "TColor.h"
#include "TStyle.h"
#include "TObjArray.h"
#include "../../../MaxCam/DmtpcSkimEvent.hh"
#include "../../../MaxCam/MaxCamImageTools.hh"
#include <vector>

TObjArray* calibTools::sumImagesWithN(TObjArray* datasets)
{
   TCanvas* c = new TCanvas("c","c",0,0,1500,500);
   c->Divide(3,1);
   
   DmtpcSkimDataset* d0 = (DmtpcSkimDataset*)datasets->At(0);
   d0->getEvent(0);
   const int ncam = d0->event()->ncamera();

   TH2F* sumOfImages[ncam];
   TH2F* hN[ncam];
   TH2F* lastSparkRef[ncam];
   for(int i=0; i<ncam; i++) 
   {
      hN[i]=(TH2F*)d0->event()->cluster(i)->getImage()->Clone("hN");hN[i]->Reset();
      sumOfImages[i]=(TH2F*)hN[i]->Clone("sumOfImages");
      lastSparkRef[i]=(TH2F*)hN[i]->Clone("lastSparkRef");
      for(int j=1; j<=lastSparkRef[i]->GetNbinsX(); j++)
      {
	 for(int k=1; k<=lastSparkRef[i]->GetNbinsY(); k++)
	 {
	    //initalize first spark ref mask to 1
	    //assumes first run has no abnormality
	    lastSparkRef[i]->SetBinContent(j,k,1);
	 }
      }
   }
   
   TH2F* hMask = (TH2F*)d0->event()->cluster(0)->getImage()->Clone("hMask");
   
   for(int n=0; n<datasets->GetEntries(); n++)
   {
      DmtpcSkimDataset* d = (DmtpcSkimDataset*)datasets->At(n);
      d->getEvent(0);
      cout << d->event()->runNumber() << endl;

      for(int i=0; i<d->nevents(); i++)
      {
	 if (i%50==0) cout << "N: " << i << endl;
	 d->getEvent(i);
	 for(int j=0; j<ncam; j++)
	 {
	    double mean = MaxCamImageTools::getMean(d->event()->image(j));
	    if(!d->event()->spark(j) && d->event()->lastspark(j)>2 && mean>0 && mean <100)
	    {  
	       hMask->Reset();
	       for(int m=1; m<=hMask->GetNbinsX(); m++)
	       {
		  for(int n=1; n<=hMask->GetNbinsY(); n++)
		  {
		     hMask->SetBinContent(m,n,1);
		  }
	       }
	       hMask->Multiply(lastSparkRef[j]);
	       for(int k=0; k<d->event()->cluster(j)->getNCluster(); k++)
	       {
		  vector<int> px = d->event()->cluster(j)->getCluster(k);
		  for(int m=0; m<int(px.size());m++)
		     hMask->SetBinContent(px[m],0);
	       }
	       for(int k=0; k<d->event()->nsparkref(j); k++)
	       {
		  hMask->SetBinContent(d->event()->sparkrefX(j,k),
				       d->event()->sparkrefY(j,k),0);
	       }
	       hN[j]->Add(hMask);
	       hMask->Multiply(d->event()->cluster(j)->getImage());
 
	       sumOfImages[j]->Add(hMask);
	    }
	 }
	 if(i==d->nevents()-1)
	 {
	    for(int j=0; j<ncam; j++)
	    {
	       for(int k=1; k<=lastSparkRef[j]->GetNbinsX(); k++)
	       {
		  for(int l=1; l<=lastSparkRef[j]->GetNbinsY(); l++)
		  {
		     lastSparkRef[j]->SetBinContent(k,l,1);
		  }
	       }
	       for(int k=0; k<d->event()->nsparkref(j); k++)
	       {
		  lastSparkRef[j]->SetBinContent(d->event()->sparkrefX(j,k),
						 d->event()->sparkrefY(j,k),0);
	       }
	    }
	 }
	 
      }

   }

   
   TObjArray* sumArray= new TObjArray(ncam);
   sumArray->SetOwner(kTRUE);

   for(int n=0; n<ncam; n++)
   {
      sumOfImages[n]->Divide(hN[n]);
      for(int i=1; i<=sumOfImages[n]->GetNbinsX(); i++)
      {
	 for(int j=1; j<=sumOfImages[n]->GetNbinsY();j++)
	 {

	    if(sumOfImages[n]->GetBinContent(i,j) <-20 || sumOfImages[n]->GetBinContent(i,j) >20)
	       sumOfImages[n]->SetBinContent(i,j,0);
	 }
      }
      sumArray->Add(sumOfImages[n]);

   }
   
   
   return sumArray;

}


TObjArray* calibTools::sumImages(DmtpcSkimDataset& ds,TString type)
{

   type.ToLower();

   ds.getEvent(0);
   const int ncam = ds.event()->ncamera();

   ds.loadDmtpcEvent(false);

   double min, max;
   if(type=="gamma")
   {min=-10;max=10;}
   else if(type=="alpha")
   {min=-1000;max=1000;}
   else
   {cout << "Type not found; quitting!" << endl; return 0;}

   TH2F* sumOfImages[ncam];
   TH2F* hN[ncam];
   int nimages[ncam];
   for(int i=0; i<ncam; i++) 
   {
      hN[i]=(TH2F*)ds.event()->cluster(i)->getImage()->Clone("hN");hN[i]->Reset();
      sumOfImages[i]=(TH2F*)hN[i]->Clone("sumOfImages");
   }
   for(int i=0; i<ds.nevents(); i++)
      //for(int i=0; i<100; i++)
   {

      if (i%50==0) cout << "N: " << i << endl;
      ds.getEvent(i);
      for(int j=0; j<ncam; j++)
      {
	 
	 if(type=="gamma")
	 {
	    if(!ds.event()->spark(j) && ds.event()->lastspark(j)>2)
	   {  
	      TH2F* hMask = (TH2F*)ds.event()->cluster(j)->getImage()->Clone("hMask");
	      hMask->Reset();
	      for(int m=1; m<=hMask->GetNbinsX(); m++)
	      {
		 for(int n=1; n<=hMask->GetNbinsY(); n++)
		 {
		    hMask->SetBinContent(m,n,1);
		 }
	      }
	      for(int k=0; k<ds.event()->cluster(j)->getNCluster(); k++)
	      {
		 vector<int> px = ds.event()->cluster(j)->getCluster(k);
		 for(int m=0; m<int(px.size());m++)
		    hMask->SetBinContent(px[m],0);
	      }
	      for(int k=0; k<ds.event()->nsparkref(j); k++)
	      {
//		 hMask->SetBinContent(ds.event()->sparkrefX(j,k),
//				      ds.event()->sparkrefY(j,k),0);
	      }
	      hN[j]->Add(hMask);
	      hMask->Multiply(ds.event()->cluster(j)->getImage());
	      sumOfImages[j]->Add(hMask);

	   }


	 }
	 else if(type=="alpha")
	 {
	    if(!ds.event()->spark(j))
	    {  
	       TH2F* hMask = (TH2F*)ds.event()->cluster(j)->getImage()->Clone("hMask");
	       hMask->Reset();
	       for(int m=1; m<=hMask->GetNbinsX(); m++)
	       {
		  for(int n=1; n<=hMask->GetNbinsY(); n++)
		  {
		     hMask->SetBinContent(m,n,1);
		  }
	       }
	       
	       hN[j]->Add(hMask);
	       hMask->Multiply(ds.event()->cluster(j)->getImage());
	       sumOfImages[j]->Add(hMask);
	       
	    }
	    
	 }



      }
   }



   TObjArray* sumArray= new TObjArray(ncam);
   sumArray->SetOwner(kTRUE);

   for(int n=0; n<ncam; n++)
   {
      sumOfImages[n]->Divide(hN[n]);
      for(int i=1; i<=sumOfImages[n]->GetNbinsX(); i++)
      {
	 for(int j=1; j<=sumOfImages[n]->GetNbinsY();j++)
	 {

	    if(sumOfImages[n]->GetBinContent(i,j) <-20 || sumOfImages[n]->GetBinContent(i,j) >20)
	       sumOfImages[n]->SetBinContent(i,j,0);
	 }
      }

      sumArray->Add(sumOfImages[n]);
   }
   
   
   return sumArray;


}

TObjArray* calibTools::defineRegions(TObjArray* sumImages)
{
   TObjArray* regions = new TObjArray(sumImages->GetEntries());
   regions->SetOwner(kTRUE);
   
   for(int i=0; i<sumImages->GetEntries(); i++)
   {
      TH2F* sumClust=(TH2F*)sumImages->At(i)->Clone("sumClust");
      sumClust->Rebin2D(2,2);
      sumClust = MaxCamImageTools::blur(sumClust,1,0.5);
      TTimeStamp* t = new TTimeStamp(0);
      MaxCamClusterImage* clust = new MaxCamClusterImage(sumClust,t);
      int ncl = MaxCamImageTools::findClustersCI(sumClust,clust,2.5,10000,128,1600);
      clust->changeImageWithThreshold((TH2F*)sumImages->At(i),sumClust->GetMaximum()*0);


      
      int size=0;
      int nth=0;
      for(int j=0; j<ncl; j++)
      {
	 vector<int> px = clust->getCluster(j);
	 if(int(px.size())>size) {size=px.size(); nth=j;}
      }

      cout << "N px in cluster: " << size << endl;
      
      vector<int> px = clust->getClusterRed(nth);
      TH2F* mask= (TH2F*)sumImages->At(i)->Clone("mask"+i);
      for(int i=0; i<px.size(); i++)
      {
	 mask->SetBinContent(px[i],5000);
      }
      for(int k=1; k<=mask->GetNbinsX(); k++)
      {
	 for(int l=1; l<=mask->GetNbinsY(); l++)
	 {
	    if(mask->GetBinContent(k,l)!=5000)
	       mask->SetBinContent(k,l,0);
	    else
	       mask->SetBinContent(k,l,1);
	 }
      }

      regions->Add(mask);
   }
   
   return regions;
}

bool calibTools::isInRegion(TH2F* region, vector<int> cluster)
{
   for(int i=0; i<cluster.size(); i++)
   {
      if(region->GetBinContent(cluster[i])==0) return false;
   }
   
   return true;
}




TH2F* calibTools::makeSpacerMask(TH2F* image, int nspacers, double* m, double* b, double* width)
{
   TH2F* mask = (TH2F*)image->Clone("mask");
   mask->Reset();
   for(int i=1; i<=mask->GetNbinsX(); i++)
   {
      for(int j=1; j<=mask->GetNbinsY(); j++)
      {
	 mask->SetBinContent(i,j,1);
	 double x0 = mask->GetXaxis()->GetBinCenter(i);
	 double y0 = mask->GetYaxis()->GetBinCenter(j);
	 for(int k=0; k<nspacers; k++)
	 {
	    double d = fabs(m[k]*x0+b[k]-y0)/sqrt(m[k]*m[k]+1);
	    if(d/width[k] < 1.3) {mask->SetBinContent(i,j,0); break;}
	 }
	 
      }
   }

   return mask;

}

TH2F* calibTools::blurWithoutSpacers(TH2F* image, TH2F* spacermask, double blurfrac)
{
   TH2F* copy = (TH2F*)image->Clone("clone");
   copy->Multiply(spacermask);
   
   TH2F* blurred = (TH2F*)image->Clone("blurred");

   int nx = blurred->GetNbinsX();
   int ny = blurred->GetNbinsY();

   for(int i=1; i<=nx; i++)
   {
      for(int j=1; j<=ny; j++)
      {	 
	 double sum=0;
	 int nin=0;

	 if(copy->GetBinContent(i,j) != 0)
	 {
	    for(int k=i-1; k<= i+1; k++)
	    {
	       for(int l=j-1; l<=j+1; l++)
	       {
		  if(k > 0 && l > 0 && k <=nx && l <= ny
		     && copy->GetBinContent(k,l)!=0)
		  {
		     nin++;
		     if(k==i && l==j) sum+=image->GetBinContent(k,l);
		     else sum+=image->GetBinContent(k,l)*blurfrac;
		  } 
	       }
	    }
	    sum = sum/(1+(nin-1)*blurfrac);
	    
	    blurred->SetBinContent(i,j,sum);
	 }


      }
   }

   return blurred;

}
