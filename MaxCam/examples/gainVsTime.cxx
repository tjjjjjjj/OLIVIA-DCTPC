/************************************************************
//
//  A script for analyzing data.
//  Run the script on any machine from cint:
//
//  .L gainVsTime.cxx
//  all()
//
*************************************************************/


MaxCamRead *ana;
int i, n, nparam;
TF1 *fun;
float x[9];
long time0;
int day0;
TGraphErrors *gr=0;
float sqrt2pi=sqrt(2*TMath::Pi());
float ccdGain=1;
TTree *nt;

// 0 = single track (cut around wire along x and make projection to y)
// 1 = wire integral (project image to x axis and integrate the whole peak).
// 2 = same as 0, but cut&count
int type=0;  


struct wireData {
  float y, yerr, mean, width, chi2;
  int ndof;
} wire0, wire1, wire2, wire3, wire4, wire5, *wire;

struct expConfig {
  int time;
  float wirehv;
  float meshhv;
  float pressure, setpress;
  int   expotime;
} expconf;

void init() {
  ana = new MaxCamRead("data/ccdrun_00073.root");
  ana->readBiasFrame();
  ana->addRun("data/ccdrun_00074.root");
  ana->addRun("data/ccdrun_00075.root");
  ana->addRun("data/ccdrun_00076.root");
  ana->addRun("data/ccdrun_00077.root");

  i=0; 
  nparam=5;

  fun = new TF1("fun","[0]+[1]*exp(-0.5*(x-[2])**2/[3]**2)"); 

  n=ana->tree()->GetEntries(); 
  cout << "TOTAL = " << n << endl;

  ana->getEvent(0);
  time0=ana->timeStamp()->Get();
  day0 =ana->timeStamp()->GetDay();

  nt = new TTree("res","");
  nt->Branch("conf",  &expconf,  "time/I:wirehv/F:meshhv:pressure:setpress:expotime/I");
  nt->Branch("wire0", &wire0, "y/F:yerr:mean:width:chi2:ndof/I");
  nt->Branch("wire1", &wire1, "y/F:yerr:mean:width:chi2:ndof/I");
  nt->Branch("wire2", &wire2, "y/F:yerr:mean:width:chi2:ndof/I");
  nt->Branch("wire3", &wire3, "y/F:yerr:mean:width:chi2:ndof/I");
  nt->Branch("wire4", &wire4, "y/F:yerr:mean:width:chi2:ndof/I");
  nt->Branch("wire5", &wire5, "y/F:yerr:mean:width:chi2:ndof/I");
}



void eventConfig() {
  expconf.time = ana->timeStamp()->Get()-time0 - (ana->timeStamp()->GetDay()-day0)*3600*9;
  expconf.wirehv = ana->wire()->currentValue;
  expconf.meshhv = ana->mesh()->currentValue;
  expconf.pressure = ana->pressure()->currentValue;  
  expconf.setpress = ana->pressure()->setValue;  
  expconf.expotime = ana->ccdConfig()->exposureTime;
}


void wireYield(int iwire) { 

  int wireMinPixel=11; int wireMaxPixel=16; // to be fixed! 
  wire=&wire0;
  switch(iwire) {
  case 1: wireMinPixel=5; wireMaxPixel=19; wire=&wire1; break;
    //case 1: wireMinPixel=22; wireMaxPixel=29; wire=&wire1; break;
  case 2: wireMinPixel=29; wireMaxPixel=34; wire=&wire2; break;
    //case 2: wireMinPixel=10; wireMaxPixel=40; wire=&wire2; break;
  case 3: wireMinPixel=38; wireMaxPixel=43; wire=&wire3; break;
  case 4: wireMinPixel=54; wireMaxPixel=59; wire=&wire4; break;
  case 5: wireMinPixel=73; wireMaxPixel=78; wire=&wire5; break;
  }
  

  // find intensity
  TH1D *hy = 0;

  switch (type) {
  case 0: hy=ana->ccdImage()->ProjectionY("_py",wireMinPixel,wireMaxPixel,"e"); break; // use only one wire, one track
  case 1: hy=ana->ccdImage()->ProjectionX("_px");  break; // use only one wire, all tracks
  default: assert(0);
  }

  // find bin density
  double brho = hy->GetNbinsX()/(hy->GetXaxis()->GetXmax()-hy->GetXaxis()->GetXmin());
  

  fun->SetParameters( hy->GetMinimum(),
		      hy->GetMaximum()-hy->GetMinimum(), 
		      hy->GetBinCenter(hy->GetMaximumBin()), // educated guesses
		      5);
  if (type==1) {
    hy->GetXaxis()->SetRange(wireMinPixel,wireMaxPixel);
    fun->SetParLimits(2,210,210);
    fun->SetParLimits(3,4.45,4.45);
  }

  hy->Fit("fun","Q0L"); // Q0L
  //c1->Update(); gSystem->Sleep(1000); 


  float I = hy->GetFunction("fun")->GetParameter(1)
    * hy->GetFunction("fun")->GetParameter(3)
    * sqrt2pi*brho
    / ccdGain;

  float IErr = I* sqrt(
    pow(hy->GetFunction("fun")->GetParError(1)/hy->GetFunction("fun")->GetParameter(1),2) +
    pow(hy->GetFunction("fun")->GetParError(3)/hy->GetFunction("fun")->GetParameter(3),2) );

  // time stamp
  wire->y     = I; 
  wire->yerr  = IErr;
  wire->mean  = hy->GetFunction("fun")->GetParameter(2);
  wire->width = hy->GetFunction("fun")->GetParameter(3);
  wire->chi2  = hy->GetFunction("fun")->GetChisquare();
  wire->ndof  = hy->GetFunction("fun")->GetNDF();

}


void makeNtuple() {
  while (i<=n-1){
    ana->getEvent(i);
    eventConfig();
    //for (int iw=0;iw<6;iw++) wireYield(iw);
    //type=0; wireYield(3);
    //type=0; wireYield(4);
    //type=0; wireYield(5);
    type=1; wireYield(1);
    type=1; wireYield(2);
    nt->Fill();
    cout <<"event="<<i<<endl;
    i++;
  }  
}


TProfile *hprof25,*hprof50,*hprof100, *hprof150, *hprof170;
void plotFrame() {
  TH2F *h=new TH2F("h","",100,1.3, 2.3, 100,1000, 10000);
  h->Draw();
  plot25Torr("same");
  plot50Torr("same");
  plot100Torr("same");
  plot150Torr("same");
  plot170Torr("same");
  TLegend *leg = new TLegend(0.3,0.6,0.5,0.8);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hprof25,"25 Torr","P");
  leg->AddEntry(hprof50,"50 Torr","P");
  leg->AddEntry(hprof100,"100 Torr","P");
  leg->AddEntry(hprof150,"150 Torr","P");
  leg->AddEntry(hprof170,"170 Torr","P");
  leg->Draw();
}
void plot25Torr(TString opt="") {
  hprof25=new TProfile("hprof25","",4,1.425,1.625);
  nt->Draw("wire1.y:conf.wirehv","wire1.y<4000&&setpress>20&&setpress<30","box"+opt);
  nt->Draw("wire1.y:conf.wirehv>>hprof25","wire1.y<4000&&setpress>20&&setpress<30","prof same");
  hprof25->SetMarkerColor(3);
}
void plot50Torr(TString opt="") {
  hprof50=new TProfile("hprof50","",3,1.525,1.675);
  nt->Draw("wire1.y:conf.wirehv","wire1.y<4000&&setpress>40&&setpress<60&&conf.wirehv<1.675","box"+opt);
  nt->Draw("wire1.y:conf.wirehv>>hprof50","wire1.y<4000&&setpress>40&&setpress<60","prof same");
  hprof50->SetMarkerColor(2);
}
void plot100Torr(TString opt="") {
  hprof100=new TProfile("hprof100","",4,1.575,1.775);
  nt->Draw("wire1.y:conf.wirehv","wire1.y<5000&&setpress>90&&setpress<120","box"+opt);
  nt->Draw("wire1.y:conf.wirehv>>hprof100","wire1.y<5000&&setpress>90&&setpress<120","prof same");
}
void plot150Torr(TString opt="") {
  hprof150=new TProfile("hprof150","",4,1.625,1.825);
  nt->Draw("wire1.y:conf.wirehv","wire1.y<5000&&setpress>140&&setpress<160","box"+opt);
  nt->Draw("wire1.y:conf.wirehv>>hprof150","wire1.y<5000&&setpress>140&&setpress<160","prof same");
  hprof150->SetMarkerColor(6);
}
void plot170Torr(TString opt="") {
  hprof170=new TProfile("hprof170","",7,1.875,2.225);
  nt->Draw("wire2.y:conf.wirehv","wire2.y<1e4&&wire2.y>0&&setpress>160&&setpress<180","box"+opt);
  nt->Draw("wire2.y:conf.wirehv>>hprof170","wire2.y<1e4&&wire2.y>0&&setpress>160&&setpress<180","prof same");
  hprof170->SetMarkerColor(4);
}


void all() {

  init();

  makeNtuple();
  
}
