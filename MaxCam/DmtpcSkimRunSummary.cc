#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TObject.h"
#include "TVectorD.h"
#include "TObject.h"
#include "DmtpcSkimRunSummary.hh"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

ClassImp(DmtpcSkimRunSummary)


DmtpcSkimRunSummary::DmtpcSkimRunSummary()
{
  outfile = NULL;
}


DmtpcSkimRunSummary::DmtpcSkimRunSummary(int cam, int nev)
{
  setNEvents(cam,nev);
}

DmtpcSkimRunSummary::~DmtpcSkimRunSummary()
{
//  if (outfile) fclose(outfile);
}


int
DmtpcSkimRunSummary::openFile()
{
  outfile = fopen(outdir + outname,"w");
  if (!outfile) return 1;
  return 0; 
}

void
DmtpcSkimRunSummary::setNEvents(int ncamera, int nev)
{
  nevents = nev;
  ncam = ncamera;
  totalSparks.resize(ncamera);
  totalBurnin.resize(ncamera);
  totalTracks.resize(ncamera);
  aveMean.resize(2*ncamera);
  aveRMS.resize(2*ncamera);
  aveMaxPixel.resize(2*ncamera);
  for (int i = 0; i < ncamera; i++){
    totalSparks[i] = 0;
    totalBurnin[i] = 0;
    totalTracks[i] = 0;
    aveMean[i] = 0;
    aveRMS[i] = 0;
    aveMaxPixel[i] = 0;
    aveMean[i+ncamera] = 0;
    aveRMS[i+ncamera] = 0;
    aveMaxPixel[i+ncamera] = 0;
  }
//  imageRMS.resize(nev*ncamera);
//  imageMean.resize(nev*ncamera);
//  totalSparks.resize(ncamera);
//  totalBurnin.resize(ncamera);
//  ntracks.resize(nev*ncamera);
//  nburnin.resize(nev*ncamera);
}

int
DmtpcSkimRunSummary::openFile(TString dir, TString file)
{
  outdir = dir;
  outname = file;
  return openFile();
}

void
DmtpcSkimRunSummary::closeFile()
{

  fclose(outfile);
}


int
DmtpcSkimRunSummary::setImageRMS(int cam, int n, double rms)
{
  if (n >= nevents || cam >=ncam){
    cout <<"Invalid event number"<<endl;
    return 1;
  }
  if (imageRMS.size() == 0) imageRMS.resize(ncam*nevents);
  imageRMS[nevents*cam + n] = rms;
  aveRMS[cam] += rms;
  aveRMS[ncam+cam] += rms*rms;
  return 0;
}

int
DmtpcSkimRunSummary::setImageMean(int cam, int n, double mean)
{
  if (n >= nevents || cam >= ncam){
    cout <<"Invalid event number"<<endl;
    return 1;
  }
  if (imageMean.size() == 0) imageMean.resize(ncam*nevents);
  imageMean[nevents*cam+n] = mean;
  aveMean[cam] += mean;
  aveMean[ncam+cam] += mean*mean;
  return 0;
}

int
DmtpcSkimRunSummary::setNTracks(int cam, int n, int nt)
{
  if (n >= nevents || cam >= ncam){
    cout <<"Invalid event number"<<endl;
    return 1;
  }
  if (ntracks.size() == 0) ntracks.resize(ncam*nevents);
  ntracks[nevents*cam+n] = nt;
  totalTracks[cam] += nt;
  return 0;
}

int
DmtpcSkimRunSummary::setNburnin(int cam, int n, int nb)
{
  if (n >= nevents || cam >= ncam){
    cout <<"Invalid event number"<<endl;
    return 1;
  }
  if (nburnin.size() == 0) nburnin.resize(ncam*nevents);
  nburnin[nevents*cam+n] = nb;
  totalBurnin[cam] += nb;
  return 0;
}

int
DmtpcSkimRunSummary::setMaxPixel(int cam, int n, double mp)
{
  if (n >= nevents || cam >= ncam){
    cout <<"Invalid event number"<<endl;
    return 1;
  }
  if (maxPixel.size() == 0) maxPixel.resize(ncam*nevents);
  maxPixel[nevents*cam+n] = mp;
  aveMaxPixel[cam] += mp;
  aveMaxPixel[cam+ncam] += mp*mp;
  return 0;
}

void
DmtpcSkimRunSummary::createImages(int cam)
{

  vector<double>::iterator rmsIt = imageRMS.begin() + nevents*cam;
  vector<double>::iterator meanIt = imageMean.begin() + nevents*cam;
  vector<int>::iterator ntrackIt = ntracks.begin() + nevents*cam;
  vector<int>::iterator burninIt = nburnin.begin() + nevents*cam;
  vector<double>::iterator maxpixIt = maxPixel.begin() + nevents*cam;

  TVectorD rmsVect(nevents);
  TVectorD meanVect(nevents);
  TVectorD trackVect(nevents);
  TVectorD eventVect(nevents);
  TVectorD burninVect(nevents);
  TVectorD maxpixVect(nevents);

  {
    TH1F rmsHist("rmshist","Camera " + (TString) cam + " RMS",40,
		 aveRMS[cam]-0.5*aveRMS[cam+ncam],aveRMS[cam]+0.5*aveRMS[cam+ncam]);
    rmsHist.SetXTitle("Image RMS");

    TH1F meanHist("meanhist","Camera " + (TString) cam + " Mean",40,
		  aveMean[cam]-0.5*aveMean[cam+ncam],aveMean[cam]+0.5*aveMean[cam+ncam]);
    meanHist.SetXTitle("Image Mean");

    int minTrack = *min_element(ntrackIt,ntrackIt+nevents);
    int maxTrack = *max_element(ntrackIt,ntrackIt+nevents);
    TH1F trackHist("trackhist","Camera " + (TString) cam + " Tracks",40,minTrack,maxTrack);
    trackHist.SetXTitle("Number of Tracks");

    TH1F maxpixHist("maxpixHist","Camera " + (TString) cam + " Max Pixel",40,
		    aveMaxPixel[cam]-0.5*aveMaxPixel[cam+ncam],aveMaxPixel[cam]+0.5*aveMaxPixel[cam+ncam]);
    maxpixHist.SetXTitle("Max Pixel Value");
//    double minBurnin = min_element(burninIt,burninIt+nevents);
//    double maxBurnin = max_element(burninIt,burninIt+nevents);
//    TH1F burninHist("burninhist","Camera " + (TString) ncam + " RBI",40,minBurnin,maxBurnin);
//  burninHist.SetXTitle("Number of RBI");  

    for (int i = 0; i < nevents; i++){
      eventVect[i] = i;
      rmsVect[i] = *rmsIt;
      meanVect[i] = *meanIt;
      trackVect[i] = *ntrackIt;
      burninVect[i] = *burninIt;
      maxpixVect[i] = *maxpixIt;      

      rmsHist.Fill(*rmsIt);
      meanHist.Fill(*meanIt);
      trackHist.Fill(*ntrackIt);
      maxpixHist.Fill(*maxpixIt);
      rmsIt++;
      meanIt++;
      ntrackIt++;
      burninIt++;
      maxpixIt++;
    }

    TCanvas c1("c1","c1",600,600);
    c1.Divide(2,2);
    c1.cd(1);
    meanHist.Draw();
    c1.cd(2);
    rmsHist.Draw();
    c1.cd(3);
    trackHist.Draw();
    c1.cd(4);
    maxpixHist.Draw();
 //   c1.cd(4);
 //   burninHist.Draw();
    TString can1name(outdir);
    can1name += "can1Cam";
    can1name += cam;
    can1name +="_" +key+".png";
    c1.Print(can1name);

  }

  {
    TString cameraName("Cam. ");
    cameraName += cam;
    TGraph rmsGr(eventVect,rmsVect);
    rmsGr.GetXaxis()->SetTitle("Event");
    rmsGr.GetYaxis()->SetTitle("Image RMS");
    rmsGr.SetTitle("Image RMS "+cameraName);
    TGraph meanGr(eventVect,meanVect);
    meanGr.GetXaxis()->SetTitle("Event");
    meanGr.GetYaxis()->SetTitle("Image Mean");
    meanGr.SetTitle("Image Mean " + cameraName);
    TGraph trackGr(eventVect,trackVect);
    trackGr.GetXaxis()->SetTitle("Event");
    trackGr.GetYaxis()->SetTitle("Number of Tracks");
    trackGr.SetTitle("Number of Tracks " + cameraName);
    TGraph burninGr(eventVect,burninVect);
    burninGr.GetXaxis()->SetTitle("Event");
    burninGr.GetYaxis()->SetTitle("Number of RBI");
    burninGr.SetTitle("Number of burnin tracks " + cameraName );   
    TGraph maxpixGr(eventVect,maxpixVect);
    maxpixGr.GetXaxis()->SetTitle("Event");
    maxpixGr.GetYaxis()->SetTitle("Max Pixel Value");
    maxpixGr.SetTitle("Max Pixel "+ cameraName);
    TCanvas c2("c2","c2",600,900);
    c2.Divide(2,3);
    c2.cd(1);
    meanGr.Draw("al*");
    c2.cd(2);
    rmsGr.Draw("al*");
    c2.cd(3);
    trackGr.Draw("al*");
    c2.cd(4);
    burninGr.Draw("al*");
    c2.cd(5);
    maxpixGr.Draw("al*");
    gPad->SetLogy(1);
    TString can2name(outdir);
    can2name += "can2Cam";
    can2name += cam;
    can2name +="_"+ key+".png";
    c2.Print(can2name);
  }

}



void
DmtpcSkimRunSummary::setStyle()
{
 gStyle->SetPalette(1);

 gStyle->SetFrameBorderMode(0);
 gStyle->SetCanvasBorderMode(0);
 gStyle->SetPadBorderMode(0);
 gStyle->SetPadColor(0);
 gStyle->SetCanvasColor(0);
 gStyle->SetTitleColor(0);
 gStyle->SetTitleBorderSize(0);
 gStyle->SetTitleFont(42);
 gStyle->SetStatColor(0);
 gStyle->SetTitleFillColor(0);
 gStyle->SetTitleColor(1);

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

  gStyle->SetTitleOffset(.7, "x");
  gStyle->SetTitleOffset(.7, "y");
  gStyle->SetTitleOffset(.7, "z");

  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetTitleSize(0.06,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleSize(0.06,"y");
  gStyle->SetLabelSize(0.05,"z");
  gStyle->SetTitleSize(0.06,"z");

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


  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  gStyle->SetLineWidth(2);
  gStyle->SetMarkerSize(0.7);
  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPalette(1);
}
void
DmtpcSkimRunSummary::beginFile()
{


//Set the average and RMS values for output into the file
  for (int cam = 0; cam < ncam; cam++){
    aveMean[cam] /= nevents;
    aveRMS[cam] /= nevents;
    aveMaxPixel[cam] /= nevents;
    aveMean[ncam + cam] = sqrt(aveMean[ncam+cam] / nevents - pow(aveMean[cam],2));
    aveRMS[ncam + cam] = sqrt(aveRMS[ncam+cam] / nevents - pow(aveRMS[cam],2));
    aveMaxPixel[ncam + cam] = sqrt(aveMaxPixel[ncam+cam] / nevents - pow(aveMaxPixel[cam],2));
  }
//Begin printing latex file.
  fprintf(outfile,"\\documentclass{article}\n\\usepackage{graphicx}\n");
  TString correctedName = outname;
  correctedName.ReplaceAll("_","\\_");
  fprintf(outfile,"\\title{Summary of Run %s}\n",correctedName.Data());
  fprintf(outfile,"\\author{DMTPC}\n");
  time_t currentTime;
  time(&currentTime);
  fprintf(outfile,"\\date{%s}\n",ctime(&currentTime));
  fprintf(outfile,"\\begin{document}\n");
  fprintf(outfile,"\\maketitle\n");
//Print number of cameras and number of events in run
  fprintf(outfile,"Number of Cameras: %i\n",ncam);
  fprintf(outfile,"\nNumber of Events: %i\n",nevents);
//Tables including total tracks, RBI, sparks and also basic image properties
  fprintf(outfile,"\\begin{table}\n");
  fprintf(outfile,"\\centering\n");
  fprintf(outfile,"\n\\begin{tabular}{ |c||c|c|c|c|c|}\n");
  fprintf(outfile,"\\hline Camera & Tracks & RBI  & Sparks & Mean & \\sigma \\\\\n");
  for (int i = 0; i < ncam; i++){;
    fprintf(outfile,"\\hline%i & %i & %i & %i & %4.2f & %4.2f  \\\\\n",
	    i,totalTracks[i],totalBurnin[i],totalSparks[i],aveMean[i],aveMean[i+ncam]);
  } 
  fprintf(outfile,"\\hline\\end{tabular}\\newline\\newline\n");

  fprintf(outfile,"\\begin{tabular}{|c||c|c|c|c|}\n"); 
  fprintf(outfile,"\\hline Camera & RMS & \\sigma & MaxPixel & \\sigma \\\\\n");
  for (int i = 0; i < ncam; i++){;
    fprintf(outfile,"\\hline %i & %4.2f & %4.2f & %4.2f & %4.2f \\\\\n",
	    i,aveRMS[i],aveRMS[i+ncam],aveMaxPixel[i],aveMaxPixel[i+ncam]);
  } 
  fprintf(outfile,"\\hline\\end{tabular}\n");

  fprintf(outfile,"\\end{table}");
}

void 
DmtpcSkimRunSummary::endFile()
{
  fprintf(outfile,"\\end{document}\n");
}

void
DmtpcSkimRunSummary::outputCamera(int cam)
{
  fprintf(outfile,"\\newpage\n");
//  fprintf(outfile,"\\section{Camera %i}\n",cam);
  fprintf(outfile,"\\begin{figure}[ht]");
  fprintf(outfile,"\\centering\n");
  fprintf(outfile,"\\includegraphics[scale=0.7]{%scan1Cam%i_%s.png}\n",outdir.Data(),cam,key.Data());
  fprintf(outfile,"\\end{figure}\n");
  fprintf(outfile,"\\begin{figure}[ht]");
  fprintf(outfile,"\\centering\n");
  fprintf(outfile,"\\includegraphics[scale=0.7]{%scan2Cam%i_%s.png}\n",outdir.Data(),cam,key.Data());
  fprintf(outfile,"\\end{figure}\n");

}

void
DmtpcSkimRunSummary::outputAll()
{
  setStyle();
  //cout <<"Style set"<<endl;
  beginFile();
  //cout <<"File begun"<<endl;
  for (int i = 0; i < ncam; i++) {
    createImages(i);
    outputCamera(i);
    //cout <<"Camera "<<i <<" output"<<endl;
  }
  endFile();
  //cout <<"File ended"<<endl;
}

int
DmtpcSkimRunSummary::makeTeXFile()
{
  int opened = openFile();
  //cout <<"File opened"<<endl;
  if (opened) {
    cerr << "Error opening file"<<endl;
    return 1;
  }
  outputAll();
  //cout <<"Output to tex file complete"<<endl;
  closeFile();
  //cout <<"File closed"<<endl;
  return 0;
}

int
DmtpcSkimRunSummary::makeTeXFile(TString dir, TString file)
{
  outdir = dir;
  outname = file;
  return makeTeXFile();
}

int
DmtpcSkimRunSummary::pdfLatex()
{
  if (system("which pdflatex")) {
    cout <<"pdflatex not found!" <<endl;
    return 1;
  }
  TString out(outdir);
  TString pdfname(outdir);
  pdfname += outname;
  pdfname.ReplaceAll(".tex","");
  out +=  outname;
  TString pdfcall("pdflatex -interaction nonstopmode -jobname ");
  pdfcall += pdfname;
  pdfcall += " ";
  pdfcall += outdir+outname;
  pdfcall += " ; wait;";
  cout <<pdfcall<<endl;
  system(pdfcall);
  return 0;
}

int
DmtpcSkimRunSummary::tar()
{
  if (system("which tar")){
    cout <<"tar not found!!"<<endl;
    return 1;
  }
  TString tarname(outdir);
  tarname += outname+".tar.gz";
  tarname.ReplaceAll(".tex","");
  TString files(outdir);
  files += "*"+key+".png";
  TString out(outdir);
  out += outname;
  out.ReplaceAll(".tex","");
  files += " " + out+".tex" + " " + out +".pdf";
  system("tar -cvzf " + tarname+" "+files+"; wait;");
  return 0;
}

void
DmtpcSkimRunSummary::cleanup(bool deletePics, bool deleteTeX, bool deletePDF)
{

  TString files;
  TString outname0 = outname;
  outname0.ReplaceAll(".tex","");
  cout <<"Preparing list of files to delete"<<endl;
  if (deletePics){
    TString pics(outdir);
    pics += "*"+key+".png ";
    files += pics;
  }
  if (deleteTeX){
    files += " "+outdir +outname0 +".tex ";
  }
  if (deletePDF){
    files += outdir+outname0+".pdf ";
  }
  TString texStuff(outdir);
  texStuff += outname0+".aux " + outdir+outname0+".log "+outdir+outname0+".nav "+outdir+outname0+".toc " + outdir + outname0+ ".snm";
  files += texStuff;
  cout <<"Deleting "<<files<<endl;
  system("rm -f "+files);
}

int 
DmtpcSkimRunSummary::runAll(bool deletePics, bool deleteTeX,bool deletePDF){
  int errorNum = 0;

  if (!makeTeXFile())  {
    if (!pdfLatex()){
      if (!tar()){
        cout <<"Everything Complete"<<endl;
      }else errorNum = 3;
    }else errorNum = 2;
  }else errorNum = 1;

  switch (errorNum){
    case 1:
	cout <<"TeX file creation failed!"<<endl;
	break;
    case 2: 
        cout <<"pdflatex failed!"<<endl;
	break;
    case 3:
        cout <<"tar failed!"<<endl;
	break;
    default:
	break;
  } 
  cout <<"Cleaning up"<<endl;
  cleanup(deletePics,deleteTeX, deletePDF);
  return errorNum;
}
