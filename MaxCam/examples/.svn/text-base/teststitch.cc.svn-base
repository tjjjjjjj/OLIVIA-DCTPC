#include "TFile.h"
#include "../Dmtpc4ShooterStitcher.hh"
#include "../DmtpcLensCorrection.hh"
#include "../DmtpcCameraMap.hh"
#include <vector>
#include <TH2.h>
#include "TApplication.h"


int main(int nargs, char **args)
{
  TApplication app("app",0,0); 
  TFile out ("stitch.root","RECREATE"); 
  Dmtpc4ShooterStitcher * stitch = new Dmtpc4ShooterStitcher("stitch"); 

  DmtpcCameraMap map; 
  map.loadMap("camera.map"); 
  std::vector<const TH2*> images; 

  TString rotation_guess[4]; 
  TString serials[4]; 

  double weights[4]; 
  double lens_params[3] = {1,0,3.2e-8}; 
  DmtpcLensCorrection * l = new DmtpcLensCorrection("lens",2,lens_params); 
  char *  suffix = "_2896.root"; 
  char *  dir = "../../4shooter/analysis/ImageAnode/masterdir/AnodeImage_4sh_"; 

//TFile f("/net/zwicky/dmtpc/data/4sh/ImageAnode/AnodeImage_4sh_081013_20110628_192904.root"); 
//  TFile f("/net/zwicky/dmtpc/cozzyd/projects/DarkMatter/4shooter/analysis/ImageAnode/masterdir/AnodeImage_4sh_A80333_20120125_235022.root"); 
  TFile f(TString::Format("%sA80333%s",dir,suffix)); 
  
//  images.push_back((const TH2*) l.correctDistortion((TH2*) f.Get("MasterAnode_A80333"))); 
  images.push_back((const TH2*)  f.Get("MasterAnode_A80333")); 
  rotation_guess[0] = map.getRot("A80333"); 
  serials[0] = TString("A80333"); 

  weights[0] = 4.42; 
 // TFile f2("/net/zwicky/dmtpc/data/4sh/ImageAnode/AnodeImage_4sh_100534_20110628_192904.root"); 
  //TFile f2("/net/zwicky/dmtpc/cozzyd/projects/DarkMatter/4shooter/analysis/ImageAnode/masterdir/AnodeImage_4sh_100534_20120125_235022.root"); 
  TFile f2(TString::Format("%s100534%s",dir,suffix)); 
//  images.push_back((const TH2*) l.correctDistortion((TH2*) f2.Get("MasterAnode_100534"))); 
  images.push_back((const TH2*)  f2.Get("MasterAnode_100534")); 
  rotation_guess[1] = map.getRot("100534"); 
  serials[1] = TString("100534"); 
  weights[1] = 6.26; 

  //TFile f3("/net/zwicky/dmtpc/data/4sh/ImageAnode/AnodeImage_4sh_110121_20110628_192904.root"); 
//  TFile f3("/net/zwicky/dmtpc/cozzyd/projects/DarkMatter/4shooter/analysis/ImageAnode/masterdir/AnodeImage_4sh_110121_20120125_235022.root"); 
  TFile f3(TString::Format("%s110121%s",dir,suffix)); 
  images.push_back((const TH2*) f3.Get("MasterAnode_110121")); 
  rotation_guess[2] = map.getRot("110121"); 
  serials[2] = TString("110121"); 
  weights[2] = 5.15; 

  //TFile f4("/net/zwicky/dmtpc/cozzyd/projects/DarkMatter/4shooter/analysis/ImageAnode/masterdir/AnodeImage_4sh_A80334_20120125_235022.root"); 
  TFile f4(TString::Format("%sA80334%s",dir,suffix)); 
  images.push_back((const TH2*) f4.Get("MasterAnode_A80334")); 
  rotation_guess[3] = map.getRot("A80334"); 
  serials[3] = TString("A80334"); 
  weights[3] = 6; 

  out.cd(); 
//  stitch->setWeights(weights); 
  stitch->setLensCorrection(l); 
  stitch->train(&images, rotation_guess, serials, true);
  stitch->writeOverlayFile("A80334.overlay",3); 
  stitch->writeOverlayFile("110121.overlay",2); 
  stitch->writeOverlayFile("100534.overlay",1); 
  stitch->writeOverlayFile("A80333.overlay",0); 
  stitch->Write(); 
  std::cout << "Done Training!" << std::endl; 
//  TFile out ("stitch.root"); 
// Dmtpc4ShooterStitcher * stitch = (Dmtpc4ShooterStitcher*) out.Get("stitch"); 
  TFile out2 ("stitched.root","RECREATE"); 
  TH2 * result = stitch->stitch(&images); 
  std::cout << result << std::endl; 
  std::cout << "Done Stitching!" << std::endl; 
  out.Close(); 
  out2.cd(); 
  result->Write(); 
  out2.Close(); 
  std::cout << "Done" <<std::endl; 
}
