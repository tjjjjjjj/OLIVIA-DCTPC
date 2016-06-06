#include "../../../MaxCam/DmtpcSkimDataset.hh"
#include "../../../MaxCam/MaxCamClusterImage.hh"
#include <iostream>


int main(int nargs, char ** args)
{

  DmtpcSkimDataset d1; 
  DmtpcSkimDataset d2; 

  d1.openRootFile(args[1]);
  d2.openRootFile(args[2]);

  if (d1.nevents() != d2.nevents())
  {
     cout << "nevents differs!!" << endl; 
     return 1; 
  }


  for (int i = 0; i < d1.nevents(); i++)
  {

    d1.getEvent(i);
    d2.getEvent(i);
    for (int c = 0; c < d1.event()->ncamera(); c++)
    {
      for (int x = 1; x <= d1.event()->cluster(c)->getImage()->GetNbinsX(); x++)
      {
        for (int y = 1; y <= d1.event()->cluster(c)->getImage()->GetNbinsY(); y++)
        {
          if (d1.event()->cluster(c)->getImage()->GetBinContent(x,y) != d2.event()->cluster(c)->getImage()->GetBinContent(x,y))
          {
            cout << "FAIL!!!" << endl; 
          }
        }
      }
    }
   

   for (int t = 0; t < d1.event()->ntriggers(); t++)
   {
      for (int n = 0; n < d1.event()->trigger_group(t)->nWaveForms(); n++)  
      {
        for (int b = 0; b < d1.event()->trigger_group(t)->getWaveform(n)->getWaveform()->GetNbinsX(); b++)
        {
           if (d1.event()->trigger_group(t)->getWaveform(n)->getWaveform()->GetBinContent(b)
               !=d2.event()->trigger_group(t)->getWaveform(n)->getWaveform()->GetBinContent(b))
               {

                  cout << "SCOPE FAILL!!!" << endl; 
               }
        }
      }
    }
  }
}
