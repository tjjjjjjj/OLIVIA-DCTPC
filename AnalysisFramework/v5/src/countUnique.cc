#include "DmtpcSkimDataset.hh" 
#include "DmtpcDataset.hh" 
#include <set> 
#include <cstring>
#include <cstdio>

int main (int nargs, char ** args) 
{
  int cam = 0; 
  int event = 0; 
  if (nargs < 2) 
  {
    printf("Usage: ./countUnique  file.root [cam=0] [event=0]\n"); 
  }

  if (nargs>2)
  {
    cam = atoi(args[2]); 
  }
  if (nargs>3)
  {
    event = atoi(args[3]); 
  }



  std::set<double> unique; 

  int size; 

  if (!strstr(args[1],"skim.root"))
  {

    DmtpcDataset d; 
    d.openRootFile(args[1]); 
    d.getEvent(event); 
    TH2 * img = d.event()->ccdData(cam); 

    size = img->GetNbinsX() * img->GetNbinsY(); 
    for (int i = 1; i < img->GetNbinsX(); i++)
    {
      for (int j = 1; j < img->GetNbinsY(); j++)
      {
        unique.insert(img->GetBinContent(i,j)); 
      }
    }
  }
  else
  {

    DmtpcSkimDataset d; 
    d.openRootFile(args[1]); 
    d.getEvent(event); 
    TH2 * img = (TH2*) d.event()->cluster(cam)->getImage(); 
    size = img->GetNbinsX() * img->GetNbinsY(); 

    for (int i = 1; i < img->GetNbinsX(); i++)
    {
      for (int j = 1; j < img->GetNbinsY(); j++)
      {
        unique.insert(img->GetBinContent(i,j)); 
      }
    }
  }


  printf("Unique: %d out of %d (%f)\n",unique.size(), size, unique.size()/double(size)); 




}
