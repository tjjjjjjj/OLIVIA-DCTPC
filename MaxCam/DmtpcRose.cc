#include "DmtpcRose.hh"

#include <iostream>
#include "TCrown.h"
#include "TPaveText.h"
#include "TPad.h"
#include "TH1.h"
#include "TString.h"
#include "TObjArray.h"
#include "TMath.h"

using namespace std;



void DmtpcRose::DrawRose(TObjArray* slices, TString opt)
{
  opt.ToUpper();
  if(opt.Contains("SAME"))
    {
      for(int i=0; i<slices->GetEntries(); i++)
        {
          slices->At(i)->Draw("SAME");
        }

    }
  else
    {
      double maxval=0;
      for(int i=0; i<slices->GetEntries(); i++)
        {
          if(((TCrown*)slices->At(i))->GetR2()>maxval)
            maxval=((TCrown*)slices->At(i))->GetR2();
        }
      maxval=maxval*1.01;

      TObjArray* Grid = new TObjArray(80);
      TObjArray* Labels = new TObjArray(14);
      TCrown* GridE;
      TPaveText* LabelT;
      for(int i=1; i<=10; i++)
        {
          float radius = float(int(maxval/10.*i));
          cout << radius << endl;
          double centerx = 0.0;
          double centery = radius*0.95;
          double deltax = radius*0.01;
          double deltay = maxval*0.03;
          LabelT = new TPaveText(centerx-deltax,radius-deltay,
                                 centerx+deltax,radius-maxval*0.005);
          TString radtext = "";
          radtext+=int(maxval/10.*i);
          LabelT->AddText(radtext);
          LabelT->SetTextSize(0.02);
          Labels->Add(LabelT);

          for(int j=0; j<8; j++)
            {
              GridE = new TCrown(0,0,0,radius,45.0*j,45.0*(j+1));
              GridE->SetLineStyle(3);
              GridE->SetFillStyle(0);
              if(i<10)
                GridE->SetNoEdges();
              else
                {
                  TString anglabel = "";
                  anglabel+=int(45*j);
                  cout << anglabel << "," << radius*cos(TMath::Pi()/8.0*j) << endl;
                  centerx = radius*cos(TMath::Pi()/4.0*j)*1.05;
                  centery = radius*sin(TMath::Pi()/4.0*j)*1.05;
                  deltax = radius*0.01;
                  deltay = radius*0.01;
                  LabelT = new TPaveText(centerx-deltax,centery-deltay,
                                         centerx+deltax,centery+deltay);
                  LabelT->AddText(anglabel);
                  LabelT->SetTextSize(0.02);
                  Labels->Add(LabelT);
                }
              Grid->Add(GridE);
            }
        }

      maxval=maxval*1.2;
      gPad->Range(-maxval,-maxval,maxval,maxval);
      for(int i=0; i<Grid->GetEntries(); i++)
        {
          Grid->At(i)->Draw();
        }
      for(int i=0; i<Labels->GetEntries(); i++)
        {
          Labels->At(i)->Draw("SAME");
        }
      for(int i=0; i<slices->GetEntries(); i++)
        {
          slices->At(i)->Draw("SAME");
        }

    }
}


// converts a TH1 w/ x-axis in RADIANS into a set of TCrowns to be plotted later by e.g. DrawRose().                
TObjArray* DmtpcRose::makeRoseSlices(TH1* h)
{

  TObjArray* Rose = new TObjArray(h->GetNbinsX());
  Rose->SetOwner(kTRUE);
  TCrown* slice;
  double rad2deg = 180.0/TMath::Pi();
  //cout << "h->GetNbinsX() = " << h->GetNbinsX() << endl;
  for(int i=1; i<=h->GetNbinsX(); i++)
    {
      //cout << "h->GetBinContent(" << i << ") = " << h->GetBinContent(i) << endl;
      //cout << "h->GetBinLowEdge(" << i << ") = " << h->GetBinLowEdge(i) << endl;
      //cout << "h->GetBinWidth(" << i << ")   = " << h->GetBinWidth(i) << endl;
      slice = new TCrown(0,0,0,h->GetBinContent(i),
                         rad2deg*h->GetBinLowEdge(i),
                         rad2deg*(h->GetBinLowEdge(i)+h->GetBinWidth(i)));
      slice->SetFillStyle(h->GetFillStyle());
      slice->SetFillColor(h->GetFillColor());
      slice->SetLineColor(h->GetLineColor());
      slice->SetLineStyle(h->GetLineStyle());
      slice->SetLineWidth(h->GetLineWidth());
      Rose->Add(slice);
    }

  return Rose;
}
