{

  TFile f("/net/zwicky/dmtpc/cozzyd/projects/DarkMatter/4shooter/analysis/ImageAnode/masterdir/AnodeImage_4sh_100534_20120125_235022.root"); 
  TH2 * img = (TH2*) f.Get("MasterAnode_100534"); 
  
  DmtpcLensCorrection c("lens",2); 
  c.setParameter(0,1); 
  c.setParameter(1,0); 
  c.setParameter(2,1e-7); 

  TH2 * distorted = c.correctDistortion(img,0,"bicubic",0,false); 
  TH2 * corrected = c.correctDistortion(distorted,0,"bicubic",0,true); 
  DmtpcRootTools::setColorStandard1(); 

  TCanvas * canv  = new TCanvas("c","c",1000,300); 
  canv->Divide(3,1); 
  canv->cd(1); 
  img->DrawCopy("colz"); 
  canv->cd(2); 
  distorted->DrawCopy("colz"); 
  canv->cd(3); 
  corrected->DrawCopy("colz"); 
}
