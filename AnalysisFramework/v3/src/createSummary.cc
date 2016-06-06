#include "../../../MaxCam/DmtpcSkimDataset.hh"
#include "../../../MaxCam/DmtpcSkimEvent.hh"
#include "../../../MaxCam/DmtpcDataset.hh"
#include "../../../MaxCam/DmtpcEvent.hh"
#include "../../../MaxCam/MaxCamConfig.hh"
#include "../../../MaxCam/MaxCamImageTools.hh"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TVectorD.h"
#include "TObject.h"
#include "TString.h"
#include "TCanvas.h"

#include <cmath>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
using namespace std;

int createSummary(TString runFile, TString skimdir = "./", TString outdir="./");
void setStyle();


int 
main(int argc, char** argv){
  if (argc < 4) {
    cerr << "Usage: not enough arguments"<<endl;
    return -1;
  }
  else if (argc > 4) {
    cerr <<"Usage: Too many arguments"<<endl;
    return -2;
  }
  TString file(argv[1]);
  TString outdir(argv[3]);
  TString skimdir(argv[2]);
  int ret = createSummary(file,skimdir,outdir);
  return ret;
}

int
createSummary(TString runfile, TString skimdir, TString outdir){

  setStyle();
  DmtpcSkimDataset d;
  d.openRootFile(skimdir+runfile+"skim.root");
  d.loadBiasFrames(true);
  d.tree()->SetBranchStatus("_clusters",0);
  d.getEvent(0);
  static int ncamera = d.event()->ncamera();
 

  //Create Plots and save data
  double biasMean[ncamera], biasRMS[ncamera];
  int nspark[ncamera], totalRBI[ncamera], totalTrack[ncamera];
  double pixelsKilled[ncamera], pkSigma[ncamera];
  double imageMean[ncamera], meanSigma[ncamera], imageRMS[ncamera], rmsSigma[ncamera]; 
  for (int i = 0; i < ncamera; i++){

    //Clear values:
    nspark[i] = 0;
    totalRBI[i] = 0;
    totalTrack[i] = 0;
    pixelsKilled[i] = 0;
    pkSigma[i] = 0;
    imageMean[i] = 0;
    meanSigma[i] = 0;
    imageRMS[i] = 0;
    rmsSigma[i] = 0;
    biasMean[i] = 0;
    biasRMS[i] = 0;

    TCanvas *c1 = new TCanvas("c1","c1");
    TString num("");
    num+=i;

    TH2F* bias = (TH2F*) d.getBiasFrame(i)->Clone("biasClone"+num);
    double bmn, brms;
    MaxCamImageTools::meanRMSNoOutliers(bias,bmn,brms);
    biasMean[i] = bmn;
    biasRMS[i] = brms;
    nspark[i] = d.tree()->Draw("_spark["+num+"]","_spark["+num+"]>0","goff");
    d.tree()->Draw("_nburnin["+num+"]","!_spark["+num+"] &&_nburnin["+num+"]>=5 && _ntracks["+num+"]>0","goff");
    d.tree()->Draw("_ntracks>>trackPlot","!_spark","goff");




    //Plot distribution of image means
    d.tree()->Draw("_image_mean["+num+"]>>meanPlot","!_spark[" + num+"]","goff");
    TH1F* meanPlot = (TH1F*) gROOT->FindObject("meanPlot");
    imageMean[i] = meanPlot->GetMean();
    meanSigma[i] = meanPlot->GetRMS();
    meanPlot->SetXTitle("Image Mean");
    meanPlot->SetTitle("Image Mean Distribution");
    meanPlot->GetXaxis()->SetRangeUser(max(imageMean[i]-3*meanSigma[i],meanPlot->GetXaxis()->GetBinLowEdge(1)),min(imageMean[i]+3*meanSigma[i],meanPlot->GetXaxis()->GetBinUpEdge(meanPlot->GetNbinsX())));
    meanPlot->Draw();
    c1->Print(outdir+runfile+"_mean"+num+".png");

    //Plot distribution of image rms
    d.tree()->Draw("_image_rms["+num+"]>>rmsPlot","!_spark[" + num+"]","goff");
    TH1F* rmsPlot = (TH1F*) gROOT->FindObject("rmsPlot");
    imageRMS[i] = rmsPlot->GetMean();
    rmsSigma[i] = rmsPlot->GetRMS();
    rmsPlot->SetXTitle("Image RMS");
    rmsPlot->SetTitle("Image RMS Distribution");
    rmsPlot->GetXaxis()->SetRangeUser(max(imageRMS[i]-3*rmsSigma[i],rmsPlot->GetXaxis()->GetBinLowEdge(1)),min(imageRMS[i]+3*rmsSigma[i],meanPlot->GetXaxis()->GetBinUpEdge(rmsPlot->GetNbinsX())));
    rmsPlot->Draw();
    c1->Print(outdir+runfile+"_rms"+num+".png");

    //Plot distribution of number of pixels killed
    d.tree()->Draw("_pixels_killed["+num+"]>>pkPlot","!_spark[" + num+"]","goff");
    TH1F* pkPlot = (TH1F*) gROOT->FindObject("pkPlot");
    pixelsKilled[i] = pkPlot->GetMean();
    pkSigma[i] = pkPlot->GetRMS();
    pkPlot->SetXTitle("Pixels Killed");
    pkPlot->SetTitle("Distribution of Number of Pixels Killed");
    pkPlot->Draw();
    c1->Print(outdir+runfile+"_pk"+num+".png");

    bias->Reset();
    bias->Rebin2D(4,4);
    bias->SetXTitle("x [pixel]");
    bias->SetYTitle("y [pixel]");
    TH2F* rbipos = (TH2F*) bias->Clone("rbipos");
    bias->SetTitle("Positions of Non-RBI tracks");
    rbipos->SetTitle("Positions of RBI");
    TH1F* trackPlot = (TH1F*) gROOT->FindObject("trackPlot");
    trackPlot->Reset();
    trackPlot->SetXTitle("Number of Tracks");
    trackPlot->SetTitle("Distribution of Tracks");
    int trackevent=0;
    int rbievent=0;

    //Make plots of things with respect to time
    TProfile trVsTime("trVsTime","Non-RBI Tracks", d.nevents()/ 40,0,d.nevents(),"s");
    TProfile rbiVsTime("rbiVsTime","RBI Events",d.nevents()/40,0,d.nevents(),"s");

    //Image mean stability
    TProfile meanVsTime("meanVsTime","Image Mean",d.nevents()/40,0,d.nevents(),"s");
    meanVsTime.GetXaxis()->SetTitle("Event Number");
    meanVsTime.GetYaxis()->SetTitle("Image Mean");
    d.tree()->Draw("_image_mean["+num+"]:_eventNumber["+num+"]>>meanVsTime","!_spark[" + num+"]","goff");
    meanVsTime.Draw("e1");
    c1->Print(outdir+runfile+"_meanVsTime"+num+".png");

    //Image rms stability
    TProfile rmsVsTime("rmsVsTime","Image RMS",d.nevents()/40,0,d.nevents(),"s");
    rmsVsTime.GetXaxis()->SetTitle("Event Number");
    rmsVsTime.GetYaxis()->SetTitle("Image RMS");
    d.tree()->Draw("_image_rms["+num+"]:_eventNumber["+num+"]>>rmsVsTime","!_spark[" + num+"]","goff"); 
    rmsVsTime.Draw("e1"); 
    c1->Print(outdir+runfile+"_rmsVsTime"+num+".png");

    for (int j = 0; j < d.nevents(); j++){
      d.getEvent(j);
      int ntrack = d.event()->ntracks(i);
      if (ntrack==0 || d.event()->spark(i)) continue;
      int tracks = 0, rbi = 0;
      bool hasTrack = 0, hasRBI=0;
      for (int k = 0; k < ntrack; k++){
        if (d.event()->nburnin(i,k) > 5){
          hasRBI=1;
          rbi++;
          totalRBI[i]++;
          rbipos->Fill(d.event()->x(i,k),d.event()->y(i,k));
        }else{
          hasTrack=1;
          tracks++;
	  totalTrack[i]++;
          bias->Fill(d.event()->x(i,k),d.event()->y(i,k));
        }
      }//Loop over tracks
      if (hasTrack) trackevent++;
      if (hasRBI) rbievent++;
      if (tracks>0) trackPlot->Fill(tracks);
      trVsTime.Fill(d.event()->eventNumber(),tracks);
      rbiVsTime.Fill(d.event()->eventNumber(),rbi);
    }//Loop over events
    trackPlot->Draw();
    c1->Print(outdir+runfile+"_" + "tr"+num+".png");

    bias->Draw("colz");
    c1->Print(outdir+runfile+"_trxy"+num+".png");
   
    rbipos->Draw("colz");
    c1->Print(outdir+runfile+"_rbixy"+num+".png");
    
    trVsTime.SetXTitle("Event Number");
    trVsTime.SetYTitle("Number of Tracks");
    trVsTime.Draw("e1"); 
    c1->Print(outdir+runfile+"_trVsTime"+num+".png");

    rbiVsTime.SetXTitle("Event Number");
    rbiVsTime.SetYTitle("Number of RBI");
    rbiVsTime.Draw("e1"); 
    c1->Print(outdir+runfile+"_rbiVsTime"+num+".png");

    //Plot track positions
    /*
    bias->Rebin2D(4,4);
    d.tree()->Draw("_x["+num+"]:_y["+num+"]>>biasClone"+num,"!_spark["+num+"]+_nburnin["+num+"]<5 && _ntracks["+num+"] && _E["+num+"]>0","colz");
    bias->SetTitle("Positions of Non-RBI Tracks");
    c1->Print(outdir+runfile+"_trxy"+num+".png");
 
    //Plot RBI positions
    d.tree()->Draw("_x["+num+"]:_y["+num+"]>>biasClone"+num,"!_spark["+num+"]+_nburnin["+num+"]>=5 && _ntracks["+num+"] && _E["+num+"]>0","colz");
    bias->SetTitle("Positions of Non-RBI Tracks");
    c1->Print(outdir+runfile+"_rbixy"+num+".png");
    */
    delete c1;
    delete bias;
    delete rbipos;
    delete trackPlot;
    delete meanPlot;
    delete rmsPlot;
    delete pkPlot;
   }//Loop over cameras

  d.loadDmtpcEvent();
  d.getEvent(0);
  TString beginTime = d.orig_event()->UTCtimeStamp()->AsString();
  beginTime.Remove(0,5);
  beginTime.ReplaceAll("+0000 (GMT) ","");
  d.getEvent(d.nevents()-1);
  TString endTime = d.orig_event()->UTCtimeStamp()->AsString();
  endTime.Remove(0,5); 
  endTime.ReplaceAll("+0000 (GMT) ","");

  //Open the tex file for output.
  FILE* outfile = fopen(outdir+runfile+".tex","w");
  if (!outfile) {
    cerr <<"Could not open TeX file for writing"<<endl;
    return 1;
  }

  //Begin printing latex headers.
  fprintf(outfile,"\\documentclass{article}\n\\usepackage{graphicx}\n");
  TString correctedName(runfile);
  correctedName.ReplaceAll("_","\\_");
  fprintf(outfile,"\\title{Summary of Run %s}\n",correctedName.Data());
  fprintf(outfile,"\\author{DMTPC}\n");

  //Start the document
  fprintf(outfile,"\\begin{document}\n");
  fprintf(outfile,"\\maketitle\n");

  //Print number of cameras and number of events in run
  fprintf(outfile,"Number of Cameras: %i\n",ncamera);
  fprintf(outfile,"\nNumber of Events: %i\n\n",d.nevents());
  fprintf(outfile,"Start: %s\n\n",beginTime.Data());
  fprintf(outfile,"End: %s\n\n",endTime.Data());

  //Print some info about the bias frames
  fprintf(outfile,"\\begin{tabular}{|c||c|c|c|}\n");
  fprintf(outfile,"\\hline Camera & Serial & Bias Mean & Bias RMS \\\\\n");
  for (int i = 0; i < ncamera; i++)
    fprintf(outfile,"\\hline %i & %s & %4.2f & %4.2f \\\\\n",i,d.orig_event()->ccdConfig(i)->serialNumber.Data(),biasMean[i],biasRMS[i]);
  fprintf(outfile,"\\hline\n\\end{tabular}\n\n");
  //Print some tables for each camera
 
  fprintf(outfile,"\n\\begin{tabular}{ |c||c|c|c|c|c|c|c|c|c|}\n");
  fprintf(outfile,"\\hline Camera & Tracks & RBI  & Sparks & Mean & $\\sigma$ \\\\\n");
  for (int i = 0; i < ncamera; i++){
    fprintf(outfile,"\\hline%i & %i & %i & %i & %4.2f & %4.2f \\\\\n",
	    i,totalTrack[i],totalRBI[i],nspark[i],imageMean[i],meanSigma[i]);
  } //Loop over cameras
  fprintf(outfile,"\\hline\\end{tabular}\n\\newline\\newline\n\n");
  fprintf(outfile,"\\begin{tabular}{|c||c|c|c|c|}\n"); 
  fprintf(outfile,"\\hline Camera & RMS & $\\sigma$ & Pixels Killed & $\\sigma$ \\\\\n");
  for (int i = 0; i < ncamera; i++){;
    fprintf(outfile,"\\hline %i & %4.2f & %4.2f & %4.2f & %4.2f \\\\\n",
	    i,imageRMS[i],rmsSigma[i],pixelsKilled[i],pkSigma[i]);
  } //Loop over cameras
  fprintf(outfile,"\\hline\\end{tabular}\n");


  //Print the plots for each camera:
  for (int i = 0; i < ncamera; i++){
    fprintf(outfile,"\\clearpage\n");
    fprintf(outfile,"\\section*{Camera %i Plots}\n",i);
    //    fprintf(outfile,"\\newpage\n");
    fprintf(outfile,"\\begin{figure}[b]\n");
    fprintf(outfile,"\\includegraphics[scale=0.30]{%s%s_tr%i.png}\n",outdir.Data(),runfile.Data(),i);
    fprintf(outfile,"\\includegraphics[scale=0.30]{%s%s_mean%i.png}\n",outdir.Data(),runfile.Data(),i);
    fprintf(outfile,"\\includegraphics[scale=0.30]{%s%s_rms%i.png}\n",outdir.Data(),runfile.Data(),i);
    fprintf(outfile,"\\includegraphics[scale=0.30]{%s%s_pk%i.png}\n",outdir.Data(),runfile.Data(),i);
    fprintf(outfile,"\\end{figure}\n");
    fprintf(outfile,"\\begin{figure}[p]\n");
    fprintf(outfile,"\\includegraphics[scale=0.30]{%s%s_trVsTime%i.png}\n",outdir.Data(),runfile.Data(),i);
    fprintf(outfile,"\\includegraphics[scale=0.30]{%s%s_rbiVsTime%i.png}\n",outdir.Data(),runfile.Data(),i);
    fprintf(outfile,"\\includegraphics[scale=0.30]{%s%s_trxy%i.png}\n",outdir.Data(),runfile.Data(),i);
    fprintf(outfile,"\\includegraphics[scale=0.30]{%s%s_rbixy%i.png}\n",outdir.Data(),runfile.Data(),i);
    fprintf(outfile,"\\end{figure}\n");
    fprintf(outfile,"\\begin{figure}[p]\n");
    fprintf(outfile,"\\includegraphics[scale=0.30]{%s%s_meanVsTime%i.png}\n",outdir.Data(),runfile.Data(),i);
    fprintf(outfile,"\\includegraphics[scale=0.30]{%s%s_rmsVsTime%i.png}\n",outdir.Data(),runfile.Data(),i);
    fprintf(outfile,"\\end{figure}\n");

  }
  //End the document
  fprintf(outfile,"\\end{document}\n");
  system("cat "+outdir+runfile+".tex");

  //Run pdflatex
  fclose(outfile);
  TString pdfcall("pdflatex -interaction batchmode -jobname ");
  pdfcall += outdir+runfile+" "+outdir+runfile+".tex";

  cout << pdfcall<<endl;
  system(pdfcall);
  
  //Tar everything and leave just the pdf:
    TString files(outdir+runfile+".tex "+ outdir+runfile+".pdf " + outdir+runfile+"*.png");
    system("tar -cvzf " + outdir+runfile+".tar.gz " +files+"; wait;");

  //Cleanup remaining files:
  TString filePath(outdir+runfile);
  TString delFiles(filePath);
  delFiles += ".tex " + filePath+".aux " + filePath+".log " + filePath+".nav "+filePath+".toc "+filePath+".snm "; 
  delFiles += filePath+"*.png";
  system("rm -f " + delFiles);
  
  return 0;
}

void setStyle(){

  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetTitleColor(0);
  gStyle->SetTitleFont(42);
  gStyle->SetStatColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleColor(1);

  gStyle->SetLegendBorderSize(1);

  // set the paper & margin sizes
  //gStyle->SetPaperSize(20,26);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.2); // 0.2
  gStyle->SetPadBottomMargin(0.2);
  gStyle->SetPadLeftMargin(0.2);

  gStyle->SetNdivisions(505,"x");
  gStyle->SetNdivisions(505,"y");

  // use large Times-Roman fonts
  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.08);

  gStyle->SetLabelFont(42,"x");
  gStyle->SetLabelFont(42,"y");
  gStyle->SetLabelFont(42,"z");
  gStyle->SetTitleFont(42,"x");
  gStyle->SetTitleFont(42,"y");
  gStyle->SetTitleFont(42,"z");

  gStyle->SetTitleOffset(1, "x");
  gStyle->SetTitleOffset(1, "y");
  gStyle->SetTitleOffset(1, "z");

  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetTitleSize(0.06,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleSize(0.06,"y");
  gStyle->SetLabelSize(0.05,"z");
  gStyle->SetTitleSize(0.06,"z");

  /*
  double stops[5] = {0,0.34,0.61,0.84,1.0};
  double red[5] = {0.0,0.0,0.87,1.0,0.51};
  double green[5] = {0.0,0.81,1.00,0.2,0.0};
  double blue[5] = {0.51,1.0,0.12,0.0,0.0};
  gStyle->CreateGradientColorTable(5,stops,red,green,blue,255);
  // gStyle->SetNumberContours(255);
  */
  // stat box
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatX(0.95);
  gStyle->SetStatY(0.95);
  gStyle->SetStatW(.200);
  gStyle->SetStatH(.125);
  gStyle->SetStatColor(0);
  gStyle->SetStatStyle(0);
  gStyle->SetStatFont(42);

  // title
  gStyle->SetTitleX(0.3);
  gStyle->SetTitleW(0.5);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);

  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(1);
  //  gStyle->SetOptStat(1111);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPalette(1);
  gStyle->SetNumberContours(255);


}
