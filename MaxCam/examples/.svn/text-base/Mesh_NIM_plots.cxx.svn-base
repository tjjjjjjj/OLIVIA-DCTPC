
TString what="angleres"; // gain, width, resolution, spark, angleres


TH1F*
makeProfile(int n, float min, float max) {
    TF1 ftmp("ftmp","gaus");
    //TF1 ftmp("ftmp","[0]*exp(0.5*((x-[1])/[2])**2)+[3]");
    TGraph *graph = (TGraph*)gPad->GetPrimitive("Graph")->Clone("graph"); 
    TH1F *hpr = new TH1F("hpr","",n, min, max);
    float tmin=0, tmax=60; // ADU/keV
    if (what=="width") tmax=1000;
    else if (what=="spark") { tmin=0.0; tmax=110; }
    else if (what=="angleres") { tmin=-2; tmax=2; }
    for (int i=0; i<n; i++) {
        TH1F htmp("htmp","",100, tmin, tmax); 
        for (int j=0; j<graph->GetN(); j++) {
            if ( graph->GetX()[j]>hpr->GetXaxis()->GetBinLowEdge(i+1) &&
                 graph->GetX()[j]<hpr->GetXaxis()->GetBinUpEdge(i+1) ) {
                htmp.Fill( graph->GetY()[j] );
                //cout << "fill " << graph->GetY()[j] << "  for " <<  graph->GetX()[j] << endl;
            }
        }
        if (htmp.GetEntries()<1) continue;
        if (what!="spark") htmp.GetXaxis()->SetRange( htmp.GetMaximumBin()-20,  htmp.GetMaximumBin()+20);
        ftmp.SetParameters( htmp.GetMaximum(), htmp.GetBinCenter( htmp.GetMaximumBin() ),  1, 0);
        htmp.Fit("ftmp","LL");
        htmp.DrawCopy();
        if (what=="resolution") {
            float res=  htmp.GetFunction("ftmp")->GetParameter(2)/htmp.GetFunction("ftmp")->GetParameter(1);
             hpr->SetBinContent(i+1, res);
             hpr->SetBinError( i+1, res*sqrt( pow(htmp.GetFunction("ftmp")->GetParError(1)/htmp.GetFunction("ftmp")->GetParameter(1),2) +
                                              pow(htmp.GetFunction("ftmp")->GetParError(2)/htmp.GetFunction("ftmp")->GetParameter(2),2) )  );
        }
        else if (what=="spark") {
            float rate= htmp.GetMean();
            hpr->SetBinContent(i+1, rate );
            hpr->SetBinError( i+1,  0 );
        }
        else if (what=="angleres") {
            hpr->SetBinContent(i+1, htmp.GetFunction("ftmp")->GetParameter(2)*180./3.14 );
            //hpr->SetBinError( i+1,  htmp.GetFunction("ftmp")->GetParError(2)*180./3.14 );
        } else {
            hpr->SetBinContent(i+1, htmp.GetFunction("ftmp")->GetParameter(1) );
            hpr->SetBinError( i+1,  htmp.GetFunction("ftmp")->GetParError(1) );
        }
        //c1->Update(); getchar(); 
    }
    delete graph;
    return hpr;
}



void cu_gain() {
    
    // Segment intensity from bins 87->77 which are
    // approximately 1cm away from the source (around 110).
    // Assume upper, lower sources have the same energy.
    //
    double d1=(110-87)*8/16.4;
    double d2=(87-77+1)*8/16.4+d1;
    MaxCamSRIM srim("SRIM_He_in_CF4_100Torr");
    srim.setPressure(50);   double dE050=srim.calcEnergyLoss(5350, d1,d2); 
    srim.setPressure(75);   double dE075=srim.calcEnergyLoss(5350, d1,d2); 
    srim.setPressure(100);  double dE100=srim.calcEnergyLoss(5350, d1,d2);
    srim.setPressure(150);  double dE150=srim.calcEnergyLoss(5350, d1,d2); 
    srim.setPressure(200);  double dE200=srim.calcEnergyLoss(5350, d1,d2); 
    cout << "dE( 50Torr)="<< dE050 << endl;    
    cout << "dE( 75Torr)="<< dE075 << endl;    
    cout << "dE(100Torr)="<< dE100 << endl;
    cout << "dE(150Torr)="<< dE150 << endl;
    cout << "dE(200Torr)="<< dE200 << endl;
    TString comm="segment.E()/";

    gStyle->SetPadRightMargin(0.05);

    if (what=="width") { comm="segment.W()/"; dE050=dE075=dE100=dE150=dE200=0.0164; }
    else if (what=="spark") { comm="segment.isDischarge()/"; dE050=dE075=dE100=dE150=dE200=0.01; }
    
    leg= new TLegend(0.75,0.65,0.90,0.85); leg->SetBorderSize(0); leg->SetFillStyle(0);
    
    TFile *f200u= new TFile("RESULTS/results_run378.root");
    TString com200=comm; com200+=dE200; com200+=":setup.wirehv";
    TString cut200u="segment.W()>4&&segment.W()<10&&segment.Mean()>280&&segment.Mean()<360&&segment.nseg()==1";
    if (what=="spark") cut200u="";
    res->Draw(com200, cut200u);
    TH1F* hg200=makeProfile(10, 0.96-0.005, 1.07+0.005); hg200->SetMarkerStyle(20);
    
    
    TFile *f150u= new TFile("RESULTS/results_run377.root");
    TString com150=comm; com150+=dE150; com150+=":setup.wirehv";
    TString cut150u="segment.W()>4&&segment.W()<10&&segment.Mean()>280&&segment.Mean()<360&&segment.nseg()==1";
    res->Draw(com150, cut150u);
    TH1F* hg150=makeProfile(11, 0.83-0.005, 0.93+0.005);  hg150->SetMarkerStyle(24); 


    TFile *f100u= new TFile("RESULTS/results_run379.root");
    TString com100=comm; com100+=dE100; com100+=":setup.wirehv";
    TString cut100u="segment.W()>4&&segment.W()<10&&segment.Mean()>280&&segment.Mean()<360&&segment.nseg()==1";
    if (what=="spark") cut100u="";
    res->Draw(com100, cut100u);
    TH1F* hg100=makeProfile(10, 0.72-0.005, 0.81+0.005);  hg100->SetMarkerStyle(21);

    
    TFile *f075u= new TFile("RESULTS/results_run380.root");
    TString com075=comm; com075+=dE075; com075+=":setup.wirehv";
    TString cut075u="segment.W()>4&&segment.W()<10&&segment.Mean()>280&&segment.Mean()<360&&segment.nseg()==1&&!segment.isDischarge()";
    if (what=="spark") cut075u="";
    res->Draw(com075, cut075u);
    TH1F* hg075=makeProfile(9, 0.66-0.005, 0.74+0.005);  hg075->SetMarkerStyle(25); 


    TFile *f050u= new TFile("RESULTS/results_run381.root");
    TString com050=comm; com050+=dE050; com050+=":setup.wirehv";
    TString cut050u="segment.W()>4&&segment.W()<10&&segment.Mean()>280&&segment.Mean()<360&&segment.nseg()==1&&!segment.isDischarge()";
    if (what=="spark") cut050u="";
    res->Draw(com050, cut050u);
    TH1F* hg050=makeProfile(7, 0.63-0.005, 0.69+0.005); hg050->SetMarkerStyle(22);

    
    float max=60; // ADU/keV
    if (what=="width") max=1000;
    else if (what=="resolution") max=0.5;
    else if (what=="spark")  max=100.0;
    TH2F* frame=new TH2F("frame","",100,0.5, 1.2, 100, 0, max); frame->SetXTitle("V_{anode} (kV)");
    if (what=="gain") frame->SetYTitle("I_{CCD} (ADC/keV)");
    else if (what=="resolution") frame->SetYTitle("#DeltaE/E");
    else if (what=="width") frame->SetYTitle("Width (#mum)");
    else if (what=="spark") frame->SetYTitle("rate (%)");
    if (what=="gain") {  hg200->Fit("expo"); hg150->Fit("expo"); hg100->Fit("expo"); hg075->Fit("expo");  hg050->Fit("expo"); }
    frame->Draw();
    hg200->Draw("same"); hg150->Draw("same"); hg100->Draw("same"); hg075->Draw("same"); hg050->Draw("same");
    TLegend* leg=new TLegend(0.6, 0.5, 0.7, 0.9); leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->AddEntry( hg050, "50 Torr",  "p");
    leg->AddEntry( hg075, "75 Torr",  "p");
    leg->AddEntry( hg100, "100 Torr", "p");
    leg->AddEntry( hg150, "150 Torr", "p");
    leg->AddEntry( hg200, "200 Torr", "p");
    leg->Draw(); 
}






void ITO_gain() {
    
    // Segment intensity from bins 80->70 which are
    // approximately 1cm away from the source (around 100).
    // Assume upper, lower sources have the same energy.
    //
    double d1=(100-77)*8/16.4;
    double d2=(77-67+1)*8/16.4+d1;
    MaxCamSRIM srim("SRIM_He_in_CF4_100Torr");
    srim.setPressure(50);   double dE050=srim.calcEnergyLoss(5350, d1,d2); 
    srim.setPressure(75);   double dE075=srim.calcEnergyLoss(5350, d1,d2); 
    srim.setPressure(100);  double dE100=srim.calcEnergyLoss(5350, d1,d2);
    srim.setPressure(150);  double dE150=srim.calcEnergyLoss(5350, d1,d2); 
    srim.setPressure(200);  double dE200=srim.calcEnergyLoss(5350, d1,d2); 
    cout << "dE( 50Torr)="<< dE050 << endl;    
    cout << "dE( 75Torr)="<< dE075 << endl;    
    cout << "dE(100Torr)="<< dE100 << endl;
    cout << "dE(150Torr)="<< dE150 << endl;
    cout << "dE(200Torr)="<< dE200 << endl;
    TString comm="segment.E()/";

    gStyle->SetPadRightMargin(0.05);
    
    if (what=="width") { comm="segment.W()/"; dE050=dE075=dE100=dE150=dE200=0.0164; }
    else if (what=="spark") { comm="segment.isDischarge()/"; dE050=dE075=dE100=dE150=dE200=0.01; }
    
    leg= new TLegend(0.75,0.65,0.90,0.85); leg->SetBorderSize(0); leg->SetFillStyle(0);
    
    TFile *f200u= new TFile("RESULTS/results_run366_upper.root");
    TString com200=comm; com200+=dE200; com200+=":setup.wirehv";
    TString cut200u="segment.W()>4&&segment.W()<12&&segment.Mean()>120&&segment.Mean()<170&&segment.nseg()==1";
    if (what=="spark") cut200u="";
    res->Draw(com200, cut200u);
    TH1F* hg200=makeProfile(9, 1-0.005, 1.08+0.005); hg200->SetMarkerStyle(20);

    
    TFile *f150u= new TFile("RESULTS/results_run367_upper.root");
    TString com150=comm; com150+=dE150; com150+=":setup.wirehv";
    TString cut150u="segment.W()>4&&segment.W()<12&&segment.Mean()>120&&segment.Mean()<170&&segment.nseg()==1";
    res->Draw(com150, cut150u);
    TH1F* hg150=makeProfile(10, 0.86-0.005, 0.95+0.005);  hg150->SetMarkerStyle(24); 


    TFile *f100u= new TFile("RESULTS/results_run368_upper.root");
    TString com100=comm; com100+=dE100; com100+=":setup.wirehv";
    TString cut100u="segment.W()>4&&segment.W()<12&&segment.Mean()>120&&segment.Mean()<170&&segment.nseg()==1";
    if (what=="spark") cut100u="";
    res->Draw(com100, cut100u);
    TH1F* hg100=makeProfile(10, 0.75-0.005, 0.84+0.005);  hg100->SetMarkerStyle(21);  

    
    TFile *f075u= new TFile("RESULTS/results_run375_upper.root");
    TString com075=comm; com075+=dE075; com075+=":setup.wirehv";
    TString cut075u="segment.W()>4&&segment.W()<12&&segment.Mean()>80&&segment.Mean()<130&&segment.nseg()==1&&!segment.isDischarge()";
    if (what=="spark") cut075u="";
    res->Draw(com075, cut075u);
    TH1F* hg075=makeProfile(8, 0.70-0.005, 0.77+0.005);  hg075->SetMarkerStyle(25); 


    TFile *f050u= new TFile("RESULTS/results_run376_upper.root");
    TString com050=comm; com050+=dE050; com050+=":setup.wirehv";
    TString cut050u="segment.W()>4&&segment.W()<12&&segment.Mean()>80&&segment.Mean()<130&&!segment.isDischarge()";
    if (what=="spark") cut050u="";
    res->Draw(com050, cut050u);
    TH1F* hg050=makeProfile(8, 0.65-0.005, 0.72+0.005); hg050->SetMarkerStyle(22);

    
    float max=50; // ADU/keV
    if (what=="width") max=1000;
    else if (what=="resolution") max=0.5;
    else if (what=="spark")  max=100.0;
    TH2F* frame=new TH2F("frame","",100,0.5, 1.2, 100, 0, max); frame->SetXTitle("V_{anode} (kV)");
    if (what=="gain") frame->SetYTitle("I_{CCD} (ADC/keV)");
    else if (what=="resolution") frame->SetYTitle("#DeltaE/E");
    else if (what=="width") frame->SetYTitle("Width (#mum)");
    else if (what=="spark") frame->SetYTitle("rate (%)");
    if (what=="gain") {  hg200->Fit("expo"); hg150->Fit("expo"); hg100->Fit("expo"); hg075->Fit("expo");  hg050->Fit("expo"); }
    frame->Draw();
    hg200->Draw("same");  hg150->Draw("same");  hg100->Draw("same");  hg075->Draw("same"); hg050->Draw("same");
    TLegend* leg=new TLegend(0.5, 0.7, 0.7, 0.85); leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->AddEntry( hg050, "50 Torr",  "p");
    leg->AddEntry( hg075, "75 Torr",  "p");
    leg->AddEntry( hg100, "100 Torr", "p");
    leg->AddEntry( hg150, "150 Torr", "p");
    leg->AddEntry( hg200, "200 Torr", "p");
    leg->Draw();
}




void
ITO_drift_voltage() {
   // Segment intensity from bins 80->70 which are
    // approximately 1cm away from the source (around 100).
    // Assume upper, lower sources have the same energy.
    //
    double d1=60*8/16.4+18;
    double d2=(70+1)*8/16.4+18;
    MaxCamSRIM srim("SRIM_He_in_CF4_100Torr");
    srim.setPressure(100);  double dE100=srim.calcEnergyLoss(5350, d1,d2);
    cout << "dE(100Torr)="<< dE100 << endl;
    TString comm="segment.E()/";

    if (what=="width") { comm="segment.W()/"; dE100=0.0164; }
    
    leg= new TLegend(0.75,0.65,0.90,0.85); leg->SetBorderSize(0); leg->SetFillStyle(0);
    
    TFile *f100u= new TFile("RESULTS/results_run371_upper.root");
    TString com100=comm; com100+=dE100; com100+=":setup.meshhv/2.5";
    TString cut100u="segment.W()>4&&segment.W()<12&&segment.Mean()>110&&segment.Mean()<190&&segment.nseg()==1";
    //TString cut100u="segment.W()>4&&segment.W()<8&&segment.Mean()>320&&segment.Mean()<380&&segment.nseg()==1";
    res->Draw(com100, cut100u);
    TH1F* hg100=makeProfile(14, (0.2-0.005)/2.5, (1.5+0.005)/2.5 ); hg100->SetMarkerStyle(20);
    hg100->Draw(); hg100->SetXTitle("Drift field (kV/cm)");
    if (what=="gain") hg100->SetYTitle("I_{CCD} (ADC/keV)");
    else if (what=="width")  hg100->SetYTitle("Width (#mum)");
}




void ITO_depthOfField() {
    // based on runs 373-374
    

    float x[]={3, 2.75, 2.5, 3.5, 3.25, 2.25};
    float xe[]={0,0,0,0,0,0};
    float y[] ={ 5.93, 5.66, 9.40, 11.40, 8.66, 12.7};
    float ye[]={ 0.24, 0.18, 0.32, 0.45,  0.38, 0.43};

    int n=sizeof(x)/sizeof(float);
    for (int i=0;i<n; i++) {
        x[i]*=25.6; // inch -> mm
        y[i]/=16.4; // pixel -> mm
        ye[i]/=16.4; // pixel -> mm
    }

    gStyle->SetOptFit();
    TGraphErrors* g= new TGraphErrors(n,x,y,xe,ye);
    g->Draw("AP");
    TF1* fun=new TF1("fun","[0]*(x-[1])**2+[2]");
    fun->SetParNames("a","x_{0}", "#sigma_{0}");
    fun->SetParameters(0.002, 75, 0.34);
    g->Fit("fun");
    g->GetHistogram()->SetXTitle("#Delta d (mm)");
    g->GetHistogram()->SetYTitle("Segment width (mm)");
}




void
ITO_drift_voltage() {
   // Segment intensity from bins 80->70 which are
    // approximately 1cm away from the source (around 100).
    // Assume upper, lower sources have the same energy.
    //
    double d1=60*8/16.4+18;
    double d2=(70+1)*8/16.4+18;
    MaxCamSRIM srim("SRIM_He_in_CF4_100Torr");
    srim.setPressure(100);  double dE100=srim.calcEnergyLoss(5350, d1,d2);
    cout << "dE(100Torr)="<< dE100 << endl;
    TString comm="segment.E()/";

    if (what=="width") { comm="segment.W()/"; dE100=0.0164; }
    
    leg= new TLegend(0.75,0.65,0.90,0.85); leg->SetBorderSize(0); leg->SetFillStyle(0);
    
    TFile *f100u= new TFile("RESULTS/results_run371_upper.root");
    TString com100=comm; com100+=dE100; com100+=":setup.meshhv/2.5";
    TString cut100u="segment.W()>4&&segment.W()<12&&segment.Mean()>110&&segment.Mean()<190&&segment.nseg()==1";
    //TString cut100u="segment.W()>4&&segment.W()<8&&segment.Mean()>320&&segment.Mean()<380&&segment.nseg()==1";
    res->Draw(com100, cut100u);
    TH1F* hg100=makeProfile(14, (0.2-0.005)/2.5, (1.5+0.005)/2.5 ); hg100->SetMarkerStyle(20);
    hg100->Draw(); hg100->SetXTitle("Drift field (kV/cm)");
    if (what=="gain") hg100->SetYTitle("I_{CCD} (ADC/keV)");
    else if (what=="width")  hg100->SetYTitle("Width (#mum)");
}



void
ITO_eventImage() {
    int col[]={10, 19,18,17,16,11,20,15,14,13,12, 1};
    gStyle->SetPalette(12,col);
    MaxCamRead* ana= new MaxCamRead("data/ccdrun_00373.root");
    ana->readBiasFrame();
    ana->findHotPixels("hotpixels.dat");
    ana->getEvent(6);
    TH2F* hccd =MaxCamImageTools::resizeImage( ana->ccdImage(), 0, 768./16.4, 0, 512./16.4);
    hccd->SetXTitle("X (mm)");
    hccd->SetYTitle("Y (mm)");
    gStyle->SetPadRightMargin(0.120);
    gStyle->SetPadLeftMargin(0.15);
    TPaveText *alpha1=new TPaveText(44.5, 20.0, 50, 23.5);
    alpha1->SetBorderSize(1);
    alpha1->SetFillStyle(1001);
    alpha1->SetFillColor(10)
    alpha1->AddText("Lower");
    alpha1->AddText("^{241}Am");
    TPaveText *alpha2=new TPaveText(44.5, 8.0, 50, 11.5);
    alpha2->SetBorderSize(1);
    alpha2->SetFillStyle(1001);
    alpha2->SetFillColor(10)
    alpha2->AddText("Upper");
    alpha2->AddText("^{241}Am");
    hccd->Draw("colz");
    alpha1->Draw();
    alpha2->Draw();
}


void mesh_gain() {
    
    // Segment intensity from bins 85->75 which are
    // approximately 1cm away from the source (around 100).
    // Assume upper, lower sources have the same energy.
    //
    double d1=(100-85)/16.4;
    double d2=(85-75+1)*8/16.4+d1;
    MaxCamSRIM srim("SRIM_He_in_CF4_100Torr");
    srim.setPressure(50);   double dE050=srim.calcEnergyLoss(5350, d1,d2); 
    srim.setPressure(75);   double dE075=srim.calcEnergyLoss(5350, d1,d2); 
    srim.setPressure(100);  double dE100=srim.calcEnergyLoss(5350, d1,d2);
    srim.setPressure(150);  double dE150=srim.calcEnergyLoss(5350, d1,d2); 
    srim.setPressure(200);  double dE200=srim.calcEnergyLoss(5350, d1,d2); 
    cout << "dE( 50Torr)="<< dE050 << endl;    
    cout << "dE( 75Torr)="<< dE075 << endl;    
    cout << "dE(100Torr)="<< dE100 << endl;
    cout << "dE(150Torr)="<< dE150 << endl;
    cout << "dE(200Torr)="<< dE200 << endl;
    TString comm="segment.E()/";

    gStyle->SetPadRightMargin(0.05);

    if (what=="width") { comm="segment.W()/"; dE050=dE075=dE100=dE150=dE200=0.0164; }
    else if (what=="spark") { comm="segment.isDischarge()/"; dE050=dE075=dE100=dE150=dE200=0.01; }
    
    leg= new TLegend(0.75,0.65,0.90,0.85); leg->SetBorderSize(0); leg->SetFillStyle(0);
    
    TFile *f200u= new TFile("RESULTS/results_run392.root");
    TString com200=comm; com200+=dE200; com200+=":setup.wirehv";
    TString cut200u="segment.W()>4&&segment.W()<10&&segment.Mean()>200&&segment.Mean()<250&&segment.nseg()==1";
    if (what=="spark") cut200u="";
    res->Draw(com200, cut200u);
    TH1F* hg200=makeProfile(10, 1.05-0.005, 1.14+0.005); hg200->SetMarkerStyle(20);
    
    
    TFile *f150u= new TFile("RESULTS/results_run391.root");
    TString com150=comm; com150+=dE150; com150+=":setup.wirehv";
    TString cut150u="segment.W()>4&&segment.W()<10&&segment.Mean()>200&&segment.Mean()<250&&segment.nseg()==1";
    res->Draw(com150, cut150u);
    TH1F* hg150=makeProfile(7, 0.93-0.005, 0.98+0.005);  hg150->SetMarkerStyle(24); 


    TFile *f100u= new TFile("RESULTS/results_run390.root");
    TString com100=comm; com100+=dE100; com100+=":setup.wirehv";
    TString cut100u="segment.W()>4&&segment.W()<10&&segment.Mean()>200&&segment.Mean()<250&&segment.nseg()==1";
    if (what=="spark") cut100u="";
    res->Draw(com100, cut100u);
    TH1F* hg100=makeProfile(7, 0.79-0.005, 0.85+0.005);  hg100->SetMarkerStyle(21);

    
    TFile *f075u= new TFile("RESULTS/results_run394.root");
    TString com075=comm; com075+=dE075; com075+=":setup.wirehv";
    TString cut075u="segment.W()>4&&segment.W()<10&&segment.Mean()>200&&segment.Mean()<250&&segment.nseg()>=0&&!segment.isDischarge()";
    if (what=="spark") cut075u="";
    res->Draw(com075, cut075u);
    TH1F* hg075=makeProfile(6, 0.75-0.005, 0.80+0.005);  hg075->SetMarkerStyle(25);


    TFile *f050u= new TFile("RESULTS/results_run393.root");
    TString com050=comm; com050+=dE050; com050+=":setup.wirehv";
    TString cut050u="segment.W()>4&&segment.W()<10&&segment.Mean()>200&&segment.Mean()<250&&segment.nseg()>=0&&!segment.isDischarge()";
    if (what=="spark") cut050u="";
    res->Draw(com050, cut050u);
    TH1F* hg050=makeProfile(7, 0.68-0.005, 0.74+0.005); hg050->SetMarkerStyle(22);
    
    
    float max=50; // ADU/keV
    if (what=="width") max=1000;
    else if (what=="resolution") max=0.5;
    else if (what=="spark")  max=100.0;
    TH2F* frame=new TH2F("frame","",100,0.5, 1.2, 100, 0, max); frame->SetXTitle("V_{anode} (kV)");
    if (what=="gain") frame->SetYTitle("I_{CCD} (ADC/keV)");
    else if (what=="resolution") frame->SetYTitle("#DeltaE/E");
    else if (what=="width") frame->SetYTitle("Width (#mum)");
    else if (what=="spark") frame->SetYTitle("rate (%)");
    if (what=="gain") {  hg200->Fit("expo"); hg150->Fit("expo"); hg100->Fit("expo"); hg075->Fit("expo");  hg050->Fit("expo");
    }
    frame->Draw();
    hg200->Draw("same"); hg150->Draw("same"); hg100->Draw("same"); hg075->Draw("same"); hg050->Draw("same");
    TLegend* leg=new TLegend(0.6, 0.5, 0.7, 0.9); leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->AddEntry( hg050, "50 Torr",  "p");
    leg->AddEntry( hg075, "75 Torr",  "p");
    leg->AddEntry( hg100, "100 Torr", "p");
    leg->AddEntry( hg150, "150 Torr", "p");
    leg->AddEntry( hg200, "200 Torr", "p");
    leg->Draw(); 
}





void
mesh_uniformity() {
    int col[]={10, 19,18,17,16,11,20,15,14,13,12, 1};
    //gStyle->SetPalette(12,col);
    MaxCamRead* ana= new MaxCamRead("data/ccdrun_00395.root");
    ana->readBiasFrame();
    ana->findHotPixels("hotpixels.dat");
    ana->addHotPixel(47,47);
    ana->addHotPixel(34,45);
    TH2F* hall=ana->accumulatePressure(-1, 0, 100); hall->SetMinimum(0); hall->SetMaximum(3e4);
    //hall=MaxCamImageTools::resizeImage( hall, 0,768./16.4, 0,512./16.4);
    hall->Draw("colz");
    TLine *l=new TLine(4,152, 760,181); l->SetLineWidth(4); l->SetLineStyle(2); l->Draw();
    l=new TLine(0,500, 400,512); l->SetLineWidth(4); l->SetLineStyle(2); l->Draw();
}



void tpc(TString opt="ito") {
    
    TCanvas *c1=0;
    if (opt.Contains("tpc")) {
        c1 = new TCanvas("c1", "c1",700,630);
        c1->Range(-120,-55,120,75);
    } else {
        c1 = new TCanvas("c1", "c1",700,180);
        c1->Range(-120,-15,140,15);
    }
    
    gStyle->SetOptStat(0);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetTickx();
    c1->SetTicky();
    c1->SetLeftMargin(0.2);
    c1->SetRightMargin(0.2);
    c1->SetTopMargin(0.05);
    c1->SetBottomMargin(0.2);
    c1->SetFrameBorderMode(0);

    TLatex *text= new TLatex; text->SetTextSize(0.15);
    
    if (opt.Contains("ito")) {
        // acrylic
       TPave *acrylic = new TPave(-105,-3.2,105,3.2,1,"br");
       if (!opt.Contains("tpc")) { text->DrawLatex(117,-3,"Acrylic"); TArrow *aacr=new TArrow(116, -1, 107, 0, 0.012, "|>"); aacr->Draw(); }
       acrylic->Draw();
       // ITO-mylar films
       TPave *itoU = new TPave(-105,-4.0,105,-3.6,0,"br"); itoU->SetFillColor(16); itoU->Draw();
       if (!opt.Contains("tpc")) { text->DrawLatex(118,2.6,"ITO"); TArrow *aito=new TArrow(116, 3.6, 107, 3.6, 0.012, "|>"); aito->Draw(); }
       TPave *itoD = new TPave(-105, 4.0,105, 3.6,0,"br"); itoD->SetFillColor(16); itoD->Draw();
    } else if (opt.Contains("mesh")) {
        TLine *anodeU = new TLine(-105,  4, 105,  4); anodeU->SetLineStyle(2);  anodeU->SetLineWidth(2);  anodeU->Draw();
        TPave *anodeL = new TPave(-105, -3.9,-92.2, 3.9, 0,"br"); anodeL->SetFillColor(12); anodeL->Draw();
        TPave *anodeR = new TPave( 105, -3.9, 92.2, 3.9, 0,"br"); anodeR->SetFillColor(12); anodeR->Draw();
        //TLine *anodeD = new TLine(-105, -4, 105, -4); anodeD->SetLineStyle(2);  anodeD->SetLineWidth(2);  anodeD->Draw();
        text->DrawLatex(118,2.6,"Mesh"); TArrow *amesh=new TArrow(116, 3.6, 107, 3.6, 0.012, "|>"); amesh->Draw();
        text->DrawLatex(118,-3,"Al"); TArrow *aal=new TArrow(116, -1, 107, 0, 0.012, "|>"); aal->Draw();
    } else if (opt.Contains("cu")) {
        TPave *g10 = new TPave(-105,-3.2,105,3.2,1,"br"); g10->SetFillStyle(3021); g10->SetFillColor(1);
       text->DrawLatex(118,-3,"G10"); TArrow *ag10=new TArrow(116, -1, 107, 0, 0.012, "|>"); ag10->Draw();
       g10->Draw();
       // copper
       TPave *copperU = new TPave(-105,4.0,105,3.6,0,"br"); copperU->SetFillColor(16); copperU->Draw();
       text->DrawLatex(118,2.6,"Cu"); TArrow *acu=new TArrow(116, 3.6, 107, 3.6, 0.012, "|>"); acu->Draw();
    }

       
   // fishing lines
   for (float x=-100; x<110; x+= 20) {
       TEllipse *ellipseU = new TEllipse(x,  5, 0.5,0.5,0,360,0); ellipseU->Draw();
       if (opt.Contains("ito")) {
           TEllipse *ellipseD = new TEllipse(x, -5, 0.5,0.5,0,360,0); ellipseD->Draw();
       }
   }
   
   // mesh
   TLine *groundU = new TLine(-105, 6, 105, 6); groundU->SetLineStyle(2);  groundU->SetLineWidth(2); groundU->Draw();
   if (!opt.Contains("tpc")) { text->DrawLatex(114, 10,"Mesh"); TArrow *amesh=new TArrow(116, 9, 107, 6, 0.012, "|>"); amesh->Draw();}
   if (opt.Contains("ito")) {
       TLine *groundD = new TLine(-105, -6, 105, -6); groundD->SetLineStyle(2);  groundD->SetLineWidth(2);  groundD->Draw();
   }

   
   // resistors
   if (opt.Contains("tpc")) {
       // drift mesh
       TLine *driftU = new TLine(-105, 51, 105, 51); driftU->SetLineStyle(2); driftU->SetLineWidth(2);  driftU->Draw();
       TLine *driftD = new TLine(-105, -51, 105, -51); driftD->SetLineStyle(2); driftD->SetLineWidth(2); driftD->Draw();

       // aluminum rings
       for (int i=0; i<5; i++) {
           TPave *meshUL = new TPave(-105, 6+i*10,-92.2, 11+i*10, 0,"br"); meshUL->SetFillColor(12); meshUL->Draw();
           TPave *meshLL = new TPave(-105,-6-i*10,-92.2,-11-i*10, 0,"br"); meshLL->SetFillColor(12); meshLL->Draw();
           TPave *meshUR = new TPave( 105, 6+i*10, 92.2, 11+i*10, 0,"br"); meshUR->SetFillColor(12); meshUR->Draw();
           TPave *meshLR = new TPave( 105,-6-i*10, 92.2,-11-i*10, 0,"br"); meshLR->SetFillColor(12); meshLR->Draw();
       }   
       
       for (int i=0; i<5; i++) {
           TLine* wireLR=new TLine(-112.5, 9+i*10, -105,   9+i*10); wireLR->Draw();
           TLine* wireUD=new TLine(-112.5, i==0 ? 0 : 8+i*10, -112.5, i==4 ? 49 : 10+i*10); wireUD->Draw();
           wireLR=new TLine(-112.5, -9-i*10, -105,-9-i*10); wireLR->Draw();
           wireUD=new TLine(-112.5, i==0 ? 0 : -8-i*10, -112.5, i==4 ? -49 : -10-i*10); wireUD->Draw();
           if (i>3) continue;
           TPave* rU=new TPave(-115, 18+i*10, -110, 10+i*10); rU->SetBorderSize(1); rU->Draw();
           TPave* rD=new TPave(-115,-18-i*10, -110,-10-i*10); rD->SetBorderSize(1); rD->Draw(); 
       }
       TLine * wireG=new TLine(-117, 0, -112.5,0); wireG->Draw();
       wireG=new TLine(-117, 0, -117, -3); wireG->Draw();
       wireG=new TLine(-120, -3, -114, -3); wireG->Draw();
       wireG=new TLine(-119, -4, -115, -4); wireG->Draw();
       wireG=new TLine(-118, -5, -116, -5); wireG->Draw();
       wireG=new TLine(-117.5, -6, -116.5, -6); wireG->Draw();
   

       // HV
       //TEllipse *anodeHV = new TEllipse(112, 0, 6, 4); anodeHV->Draw();
       TText *anodeHV_text = new TText(107, -1, "+V"); anodeHV_text->SetTextSize(0.06); anodeHV_text->Draw();
       TText *driftHV_text = new TText(107, -50, "-V"); driftHV_text->SetTextSize(0.06); driftHV_text->Draw();
       driftHV_text = new TText(107, 44, "-V"); driftHV_text->SetTextSize(0.06); driftHV_text->Draw();
       

       // Am241 sources
       TPaveText *alpha1=new TPaveText(60, 25.0, 80, 30);
       alpha1->SetBorderSize(1);
       alpha1->SetFillStyle(1001);
       alpha1->SetFillColor(10);
       alpha1->AddText("^{241}Am");
       alpha1->Draw();
       TArrow *a1=new TArrow(50, 27.5, 60, 27.5, 0.03, "<|");
       a1->Draw();
       TPaveText *alpha2=new TPaveText(60,-25.0, 80, -30);
       alpha2->SetBorderSize(1);
       alpha2->SetFillStyle(1001);
       alpha2->SetFillColor(10);
       alpha2->AddText("^{241}Am");
       alpha2->Draw();
       TArrow *a2=new TArrow(50, -27.5, 60, -27.5, 0.03, "<|");
       a2->Draw();   

       // Fe55 source
       TBox *fe55=new TBox(20, 52, 40, 54); fe55->SetFillStyle(1001); fe55->SetFillColor(14); fe55->Draw();
       text->SetTextSize(0.035); text->DrawLatex(20,55,"^{55}Fe");

       // lens
       TArc *lens=new TArc(0, 0, 65, 80, 100); lens->SetNoEdges(); lens->SetFillStyle(0); lens->Draw();
       lens=new TArc(0, 128, 65, 260, 280); lens->SetNoEdges(); lens->SetFillStyle(0); lens->Draw();

       // Camera
       TBox *ccd=new TBox(-20, 66, 20, 75); ccd->SetFillStyle(1001); ccd->SetFillColor(14); ccd->Draw();
       text->DrawLatex(-50,72,"Camera");
       text->DrawLatex(-50,67,"+ Lens");
      
       
   }

   
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
   c1->ToggleToolBar();
}



//
//    sensitivity
//
TGraph *plotTRIM(TString file) {
    ifstream fin(file);
    if (!fin.is_open()) { cout << "File not found" << endl; return; }
    char buf[256];
    for (int i=0; i<15; i++) { fin.getline(buf, 256); }//cout << buf << endl; }
    int ion, i=0;
    float e0, x0, y0, z0, eloss0;
    float e1, x1, y1, z1, eloss1;
    float x[100000], y[100000];
    while (!fin.eof()) {
        if (!(i%2)) {
            fin >> ion >>e0 >> x0 >> y0 >> z0 >> eloss0;
            //cout << i << "  " << e0 << "  " << x0 << endl;
            x[i/2]=e0;
        } else {
            fin >> ion >> e1 >> x1 >> y1 >> z1 >> eloss1;
            y[i/2]=sqrt( pow(x1-x0,2) + pow(y1-y0,2) + pow(z1-z0,2) ) * 1e-7;
            //cout << i << "  " << e1 << "  " << x1 << endl;
        }
        if (fin.eof()) break;
        i++;
        if (i>2000) break;
    }
    cout << "Total number of tracks = " << i/2 << "   for " << file << endl;
    g= new TGraph(i/2, x, y);
    return g;
}

float sigmaCCD=26; // sigma ADC
float gainCCD=10; // ADC/keV
float Lbin=0.5; //8./16.4; // bin in mm
void convertToSignalSignificance(TGraph *g) {
    for (int i=0; i<g->GetN(); i++) {
        g->GetY()[i] = g->GetX()[i]/(sigmaCCD/gainCCD*sqrt(g->GetY()[i]/Lbin));
    }
}

void convertToRangeOverEnergy(TGraph *g) {
    for (int i=0; i<g->GetN(); i++) {
        g->GetY()[i] = log( g->GetY()[i] / g->GetX()[i] );
    }
}

void smearRange(TGraph *g) {
    // add 1mm resolution
    for (int i=0; i<g->GetN(); i++) {
        g->GetY()[i] = gRandom->Gaus(g->GetY()[i], 1);
    }
}


void mcrange(TString what) {

    float press=100; // Torr

    if (what=="srim") {
        MaxCamSRIM *F= new MaxCamSRIM("../tables/SRIM_F_in_CF4_100Torr_lowE");
        F->setPressure(press);
        MaxCamSRIM *alpha= new MaxCamSRIM("../tables/SRIM_He_in_CF4_100Torr");
        alpha->setPressure(press);
        TGraph *g=F->getRangeVsEnergy(); g->SetLineWidth(3);
        for (int i=0; i<g->GetN(); i++) g->GetY()[i] = g->GetX()[i]/(sigmaCCD/gainCCD*sqrt(g->GetY()[i]/Lbin));
        g->Draw("L");   
        TGraph *ga=alpha->getRangeVsEnergy(); ga->SetLineStyle(2); ga->SetLineWidth(3);
        for (int i=0; i<ga->GetN(); i++) ga->GetY()[i] = ga->GetX()[i]/(sigmaCCD/gainCCD*sqrt(ga->GetY()[i]/Lbin));
        ga->Draw("L");
    }
    else if (what=="geant") {
        TString plotcom="sqrt((xf-x)**2+(yf-y)**2+(zf-z)**2):Ekin*1e3";
        fel = new TFile("MonteCarloSamples/mcdark_el.root");
        trackList->Draw(plotcom,"trackid==1");
        TGraph *gel = (TGraph*)gPad->GetPrimitive("Graph")->Clone("gel"); gel->SetMarkerStyle(7); gel->SetMarkerColor(2);
        smearRange(gel);
        
        falpha = new TFile("MonteCarloSamples/mcdark_alpha.root");
        trackList->Draw(plotcom,"trackid==1");
        TGraph *galpha = (TGraph*)gPad->GetPrimitive("Graph")->Clone("galpha"); galpha->SetMarkerStyle(7); galpha->SetMarkerColor(3);
        smearRange(galpha);
        
        fC = new TFile("MonteCarloSamples/mcdark_C.root");
        trackList->Draw(plotcom,"trackid==1");
        TGraph *gC = (TGraph*)gPad->GetPrimitive("Graph")->Clone("gC"); gC->SetMarkerStyle(7); gC->SetMarkerColor(5);
        smearRange(gC);

        fF = new TFile("MonteCarloSamples/mcdark_F.root");
        trackList->Draw(plotcom,"trackid==1");
        TGraph *gF = (TGraph*)gPad->GetPrimitive("Graph")->Clone("gF"); gF->SetMarkerStyle(7); gF->SetMarkerColor(4);
        smearRange(gF);

        gPad->Clear();
        galpha->GetHistogram()->SetXTitle("Energy (keV)");
        //galpha->GetHistogram()->SetYTitle("E/#sigma_{E,CCD}"); galpha->SetMaximum(150);
        //galpha->GetHistogram()->SetYTitle("log(Range/Energy)"); galpha->SetMaximum(150);
        galpha->GetHistogram()->SetYTitle("Range (mm)"); galpha->SetMaximum(150);
        galpha->Draw("AP");
        gel->Draw("P");
        gC->Draw("P");
        gF->Draw("P");
    }
    else if (what=="trim") {
        TGraph *gF=plotTRIM("MonteCarloSamples/EXYZ_F.txt");  gF->SetMarkerStyle(7);gF->SetMarkerColor(17);
        //for (int i=0; i<gF->GetN(); i++) gF->GetY()[i] = gF->GetX()[i]/(sigmaCCD/gainCCD*sqrt(gF->GetY()[i]/Lbin));
	gF->Draw("P");

        TGraph *gC=plotTRIM("MonteCarloSamples/EXYZ_C.txt");  gC->SetMarkerStyle(7); gC->SetMarkerColor(13);
        //for (int i=0; i<gC->GetN(); i++) gC->GetY()[i] = gC->GetX()[i]/(sigmaCCD/gainCCD*sqrt(gC->GetY()[i]/Lbin));
	gC->Draw("P");

        TGraph *gAlpha=plotTRIM("MonteCarloSamples/EXYZ_Alpha.txt"); gAlpha->SetMarkerStyle(7);
	gAlpha->Draw("P");
    }

}




//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//            Cf-252 run
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

float dataMinE=100;
float dataMaxE=950;
int   dataNbin=20;

void Cf_gain_calibration() {
    double d1=(17-0)*8/16.4;
    double d2=(17-7+1)*8/16.4+d1;
    MaxCamSRIM srim("SRIM_He_in_CF4_100Torr");
    srim.setPressure(75);   double dE075=srim.calcEnergyLoss(5350, d1,d2); 
    cout << "dE( 75Torr)="<< dE075 << endl;    
    TString comm="segment.E()/sqrt(1+segment.correlation()**2)/";

    gStyle->SetPadRightMargin(0.05);
    
    TChain res("res");
    res.Add("RESULTS/results_run382.root");
    res.Add("RESULTS/results_run385.root");
    res.Add("RESULTS/results_run387.root");
    TString com075=comm; com075+=dE075; com075+=":(setup.time-883.385e6)/3600";
    TString cut075 ="(segment.Mean()>260&&segment.Mean()<320&&setup.wirehv>0.735&&setup.run==382)";
    cut075 += "||(segment.Mean()>180&&segment.Mean()<250&&setup.run==385)";
    cut075 += "||(segment.Mean()>200&&segment.Mean()<270&&setup.run==387)";
    res.Draw(com075, cut075 );
    TH1F* hg075=makeProfile(7, 0-0.5, 6+0.5); hg075->SetMarkerStyle(20); hg075->SetXTitle("hours"); hg075->SetYTitle("Gain (ADU/keV)");
    hg075->Draw();
}



TString CfDataCut="!segment.isDischarge()&&segment.n()>5&&!segment.isBorderline()&&(segment.Ix()/segment.Iy()>5||segment.Iy()/segment.Ix()>5)&&segment.L()*8/16.4>-2+segment.E()/2700&&segment.L()*8/16.4<12&&!(setup.run==384&&setup.event==1474||setup.run==383&&setup.event==1634)";
//TString CfMcCut="segment.n()>5&&(segment.Ix()/segment.Iy()>5||segment.Iy()/segment.Ix()>5)";
TString CfMcCut="segment.n()>5&&(segment.Ix()/segment.Iy()>5||segment.Iy()/segment.Ix()>5)";
void Cf_energy_vs_length(TString opt="") {
    TChain res("res");
    gStyle->SetMarkerSize(1.2);
    if (opt.Contains("mc")) {
        frame=new TH2F("frame","", dataNbin,dataMinE,dataMaxE, 20, 0,14); frame->SetXTitle("Energy (keV)"); frame->SetYTitle("Range (mm)");  frame->SetLineColor(12);
        res.Add("RESULTS/results_mcrun_00057.root");
        res.Draw("segment.L()*8/16.4:segment.E()/27>>frame",CfMcCut, "box");
    } else if (opt.Contains("alpha")) {
        res.Add("RESULTS/results_mcrun_00058.root");
        res.Draw("segment.L()*8/16.4:segment.E()/27",CfMcCut, "box same");
        
    } else {
        res.Add("RESULTS/results_run383.root");
        res.Add("RESULTS/results_run384.root");
        res.Add("RESULTS/results_run386.root");
        res.Draw("segment.L()*8/16.4:segment.E()/27",CfDataCut,"same");
    }
}

void Cf_skewness(TString opt="") {
    TChain res("res");
    gStyle->SetMarkerSize(1.2);
    if (opt.Contains("mc")) {
        frame=new TH2F("frame","",  dataNbin,dataMinE,dataMaxE, 20, -1,1); frame->SetXTitle("Energy (keV)"); frame->SetYTitle("Skewness"); frame->SetLineColor(12);
        res.Add("RESULTS/results_mcrun_00057.root");
        res.Draw("segment.skewness():segment.E()/27>>frame",CfMcCut,"box");
        TLine *l=new TLine(dataMinE,0, dataMaxE, 0); l->Draw();
    } else if (opt.Contains("alpha")) {
        res.Add("RESULTS/results_mcrun_00058.root");
        res.Draw("segment.skewness():segment.E()/27",CfMcCut, "box same");
    } else {
        res.Add("RESULTS/results_run383.root");
        res.Add("RESULTS/results_run384.root");
        res.Add("RESULTS/results_run386.root");
        res.Draw("segment.skewness():segment.E()/27",CfDataCut,"same");
    }  
}


void Cf_cosRecoil(TString opt="") {
    TChain res("res");
    gStyle->SetMarkerSize(1.2);
    if (opt.Contains("mc")) {
        res.Add("RESULTS/results_mcrun_00057.root"); 
        //res.Add("RESULTS/results_mcrun_000151.root"); 
        if (opt.Contains("2d")) {
            frame=new TH2F("frame","", dataNbin,dataMinE,dataMaxE, 20, -1,1); frame->SetXTitle("Energy (keV)"); frame->SetYTitle("cos(#theta_{Recoil})"); frame->SetLineColor(12);
            res.Draw("-segment.cosRecoil():segment.E()/27>>frame",CfMcCut+"&&segment.skewness()>0","box");
            res.Draw("segment.cosRecoil():segment.E()/27>>+frame",CfMcCut+"&&segment.skewness()<=0","box");
            TLine *l=new TLine(dataMinE,0, dataMaxE, 0); l->Draw();
        } else {
            hmc= new TH1F("hmc","",12,-1,1); hmc->SetXTitle("cos(#theta_{Recoil})");
            res.Draw("-segment.cosRecoil()>>hmc",CfMcCut+"&&segment.skewness()>0","goff");
            res.Draw("segment.cosRecoil()>>+hmc",CfMcCut+"&&segment.skewness()<=0","goff");
        }
    } else if (opt.Contains("alpha")) {
        res.Add("RESULTS/results_mcrun_00058.root");
        res.Draw("-segment.cosRecoil():segment.E()/27",CfMcCut+"&&segment.skewness()>0", "box same");
        res.Draw("segment.cosRecoil():segment.E()/27",CfMcCut+"&&segment.skewness()<=0", "box same");
    } else {
        res.Add("RESULTS/results_run383.root");
        res.Add("RESULTS/results_run384.root");
        res.Add("RESULTS/results_run386.root");
        if (opt.Contains("2d")) {
            res.Draw("-segment.cosRecoil():segment.E()/27",CfDataCut+"&&segment.skewness()>0","same");
            res.Draw("segment.cosRecoil():segment.E()/27",CfDataCut+"&&segment.skewness()<=0","same");
        } else {
            if (opt.Contains("nht")) {
                hdata= new TH1F("hdata","",6,0,1);hdata->SetXTitle("cos(#theta_{Recoil})");
                res.Draw("segment.cosRecoil()>>hdata",CfDataCut,"goff");
            } else {
                res.Draw("-segment.cosRecoil()>>hdata",CfDataCut+"&&segment.skewness()>0", "goff");
                res.Draw("segment.cosRecoil()>>+hdata",CfDataCut+"&&segment.skewness()<=0", "goff");
                hdata->Draw(" E");
            }
        }
    }  
}
///res.Scan("setup.event:setup.run:segment.E()/27:segment.L()*8/16.4:segment.skewness()",CfDataCut)    


void Cf_Er(TString opt="") {
    TChain res("res");
    if (opt.Contains("mc")) {
        hmc=new TH1F("hmc","",5, dataMinE, dataMaxE);
        res.Add("RESULTS/results_mcrun_00057.root"); 
        res.Draw("segment.E()/27>>hmc",CfMcCut);
        if (hdata) hmc->SetNormFactor( hdata->Integral() );
        hmc->Draw(opt);
    } else {
        hdata=new TH1F("hdata","",5,dataMinE, dataMaxE);
        res.Add("RESULTS/results_run383.root");
        res.Add("RESULTS/results_run384.root");
        res.Add("RESULTS/results_run386.root");
        res.Draw("segment.E()/27>>hdata",CfDataCut,"goff");
        hdata->Draw(opt);
        hdata->SetXTitle("Energy (keV)");
    }
    
}



Cf_energy_efficiency() {
    TChain res("res");
    gStyle->SetMarkerSize(1.2);
    res.Add("RESULTS/results_mcrun_00057.root");
    hmc= new TH1F("hmc","",47,10,950);
    hmctruth= new TH1F("hmctruth","",47,10,950);  hmctruth->SetXTitle("E (keV)"); hmctruth->SetMinimum(0);
    res.Draw("segment.E()/27>>hmc",CfMcCut,"goff");
    res.Draw("segment.E()/27>>hmctruth","","goff");
    hmctruth->Draw();
    hmc->Draw("same");
}

Cf_energy_Q(TString what="", TString opt="") {
    TChain res("res");
    gStyle->SetMarkerSize(1.2);
    res.Add("RESULTS/results_mcrun_00057.root");
    //res.Add("RESULTS/results_mcrun_000152.root");
    float min=0, max=410;
    int n=10; //10
    hmc= new TH1F("hmc","",n,min,max); hmc->SetXTitle("E (keV)"); hmc->SetYTitle("Q_{HT}"); hmc->SetLineWidth(3); hmc->SetMinimum(0);
    /* hmc->SetMarkerStyle(1);
    hmc->SetFillStyle(4100);
    hmc->SetFillColor(4);
    hmc->SetMarkerColor(4);
    hmc->SetLineColor(4);
    hmc->GetXaxis()->SetAxisColor(1);
    hmc->GetXaxis()->SetLabelColor(10);
    hmc->GetXaxis()->SetTitleColor(10);
    hmc->GetYaxis()->SetAxisColor(1);
    hmc->GetYaxis()->SetLabelColor(10);
    hmc->GetYaxis()->SetTitleColor(10);*/


    hmctruth= new TH1F("hmctruth","",n,min,max);  
    hmcwrong= new TH1F("hmcwrong","",n,min,max);
    if (!what.Contains("truth")) {
        res.Draw("segment.E()/27>>hmc",CfMcCut,"goff");
        res.Draw("segment.E()/27>>hmctruth","","goff");
        res.Draw("segment.E()/27>>hmcwrong",CfMcCut+"&&segment.skewness()>0","goff");
    } else {
        res.Draw("setup.emc>>hmc",CfMcCut,"goff");
        res.Draw("setup.emc>>hmctruth","","goff");
        res.Draw("setup.emc>>hmcwrong",CfMcCut+"&&segment.skewness()>0","goff");
    }
    for (int i=1; i<=hmc->GetNbinsX(); i++) {
        float eff = hmc->GetBinContent(i)/hmctruth->GetBinContent(i);
        float w = hmcwrong->GetBinContent(i)/hmc->GetBinContent(i);
        hmc->SetBinContent( i, eff*pow(1-2*w,2) );
    }
    hmc->DrawCopy(opt);
}




Cf_event_images() {

    int col[]={10, 19,18,17,16,11,20,15,14,13,12, 1};
    //gStyle->SetPalette(12,col);

    MaxCamRead* ana= new MaxCamRead("data/ccdrun_00383.root");
    ana->readBiasFrame();
    ana->findHotPixels("hotpixels.dat");

    ana->getEvent(1248);
    TH2F *hall=ana->ccdImage();
    hall->GetXaxis()->SetRange(1,35); hall->GetYaxis()->SetRange(1,22); hall->SetMinimum(-100); hall->SetMaximum(850);
    hall->Draw("colz");
    c1->Print("run00383_ev1249_recoilCandidate.pdf");

    ana->getEvent(762);
    hall=ana->ccdImage();
    hall->GetXaxis()->SetRange(5,40); hall->GetYaxis()->SetRange(35,57); hall->SetMinimum(-100); hall->SetMaximum(650);
    hall->Draw("colz");
    c1->Print("run00383_ev0763_recoilCandidate.pdf");

    ana->getEvent(948);
    hall=ana->ccdImage();
    hall->GetXaxis()->SetRange(20,36); hall->GetYaxis()->SetRange(55,62); hall->SetMinimum(-100); hall->SetMaximum(750);
    hall->Draw("colz");
    c1->Print("run00383_ev0949_recoilCandidate.pdf");

    ana->getEvent(418);
    hall=ana->ccdImage();
    hall->GetXaxis()->SetRange(10,45); hall->GetYaxis()->SetRange(18,40); hall->SetMinimum(-100); hall->SetMaximum(750);
    hall->Draw("colz");
    c1->Print("run00383_ev0419_recoilCandidate.pdf");

    MaxCamRead* ana= new MaxCamRead("data/ccdrun_00384.root");
    ana->readBiasFrame();
    ana->findHotPixels("hotpixels.dat");

     
    ana->getEvent(238);
    hall=ana->ccdImage();
    hall->GetXaxis()->SetRange(20,60); hall->GetYaxis()->SetRange(5,21); hall->SetMinimum(-100); hall->SetMaximum(850);
    hall->Draw("colz");
    c1->Print("run00384_ev0239_recoilCandidate.pdf");

    ana->getEvent(1202);
    hall=ana->ccdImage();
    hall->GetXaxis()->SetRange(55,90); hall->GetYaxis()->SetRange(28,50); hall->SetMinimum(-100); hall->SetMaximum(550);
    hall->Draw("colz");
    c1->Print("run00384_ev1203_recoilCandidate.pdf");


    MaxCamRead* ana= new MaxCamRead("data/ccdrun_00386.root");
    ana->readBiasFrame();
    ana->findHotPixels("hotpixels.dat");
    ana->getEvent(1462);
    hall=ana->ccdImage();
    hall->GetXaxis()->SetRange(20,60); hall->GetYaxis()->SetRange(34,56); hall->SetMinimum(-100); hall->SetMaximum(500);
    hall->Draw("colz");
    c1->Print("run00386_ev1462_recoilCandidate.pdf");
   
}





void
Cf_recoilAngle_resolution() {
    
    TChain res("res");
    gStyle->SetMarkerSize(1.2);
    float min=20, max=500;
    hangle= new TH1F("hangle","",50,0,TMath::Pi());
    hmcangle= new TH1F("hmcangle","",50,0,TMath::Pi());
    hmcangle= new TH1F("hmcangle","",50,0,TMath::Pi());
    hresid= new TH1F("hresid","",50,-1,1);
    //res.Add("RESULTS/results_mcrun_00057.root");
    res.Add("RESULTS/results_mcrun_000152.root");
    res.Draw("acos(segment.cosRecoil())>>hangle",CfMcCut,"goff");
    //res.Draw("3.14159-abs(mcsegment.Vect().Phi())>>hmcangle",CfMcCut,"goff");
    res.Draw("3.14159-abs(mcsegment.Vect().Phi())-acos(segment.cosRecoil()):setup.emc");//,CfMcCut);
    //res.Draw("3.14159-abs(mcsegment.Vect().Phi())-acos(segment.cosRecoil())>>hresid",CfMcCut,"goff");
    //hangle->Draw(); hmcangle->Draw("same");
    //hresid->Draw();
    TH1F *h=makeProfile(8, min, max);
    h->SetXTitle("Energy (keV)");
    h->SetYTitle("Resolution (degrees)");
    h->SetLineWidth(3);
    h->Draw();
}







//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//           WIMPs
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////


float eMin=50;
float eMax=150;
// NIM = 0.01
// Soudan 0.003 bg/kg/y/keV,
// DUSEL(upper) 1e-5 bg/kg/y/keV
// NEWAGE: 1686 in (100-400keV) or 2/3*1686 in (e100-200keV) and 0.15kg*day so <rate>=2/3*1686/0.151*365/100keV=2.7e4 for (100-200keV)
float backgroundRatePerKgYearKeV=0.01;
float years=1.;
float kilos=100.;
float wimpDensity=0.3; // GeV/cm3
float wimpMass=100; // GeV
bool  hasHeadTail=true;
bool  is2D=true;
float targetIonMass=19; // GeV
bool  plotFit=true;

TString WimpCut="segment.n()>4";

void
WIMP_cosRecoil() {

    float cosMin= hasHeadTail ? -1 : 0;
    
    TChain res("res");
    gStyle->SetMarkerSize(1.2);
    if (wimpMass==800) res.Add("RESULTS/results_mcrun_00063.root");
    else if (wimpMass==400) res.Add("RESULTS/results_mcrun_00062.root");
    else if (wimpMass==200) res.Add("RESULTS/results_mcrun_00061.root");
    else if (wimpMass==100) res.Add("RESULTS/results_mcrun_00059.root");
    else if (wimpMass==80) res.Add("RESULTS/results_mcrun_00068.root");
    else if (wimpMass==50) res.Add("RESULTS/results_mcrun_00064.root");
    else if (wimpMass==40) res.Add("RESULTS/results_mcrun_00067.root");
    else if (wimpMass==30) res.Add("RESULTS/results_mcrun_00066.root");
    else if (wimpMass==25) res.Add("RESULTS/results_mcrun_00065.root");
    else { cout << "Simulation missing for WIMP mass = " << wimpMass << endl; exit(1); }
    if (is2D) {
        hmc =new TH2F("hmc","", int(eMax-eMin)/20,eMin,eMax, 10, cosMin,1); hmc->SetXTitle("Energy (keV)"); hmc->SetLineColor(12);
        hmcbkg=new TH2F("hmcbkg","", 1,eMin,eMax,  1, cosMin,1); hmcbkg->SetBinContent(1,1,1);
        if (hasHeadTail) {
            hmc->SetYTitle("cos(#theta_{Recoil})"); 
            res.Draw("-segment.cosRecoil():segment.E()/27>>hmc",WimpCut+"&&segment.skewness()>0","box");
            res.Draw("segment.cosRecoil():segment.E()/27>>+hmc",WimpCut+"&&segment.skewness()<=0","box");
            TLine *l=new TLine(eMin,0, eMax, 0); l->Draw();
        } else {
            hmc->SetYTitle("|cos(#theta_{Recoil})|"); 
            res.Draw("segment.cosRecoil():segment.E()/27>>hmc",WimpCut,"box");
        }
    } else {
        hmc= new TH1F("hmc","",6,cosMin,1);
        if (!hasHeadTail) {
            hmc->SetXTitle("|cos(#theta_{Recoil})|");
            res.Draw("segment.cosRecoil()>>hmc",WimpCut,"goff");
        } else {
            hmc->SetXTitle("cos(#theta_{Recoil})");
            res.Draw("-segment.cosRecoil()>>hmc",WimpCut+"&&segment.skewness()>0","goff");
            res.Draw("segment.cosRecoil()>>+hmc",WimpCut+"&&segment.skewness()<=0","goff");
        }
        hmc->Draw("same");
    }
}

TH1F*
WIMP_energy_Q(TString what="pdf truth", TString opt="") {
    TChain res("res");
    gStyle->SetMarkerSize(1.2);
    res.Add("RESULTS/results_mcrun_000132.root");
    float min=eMin, max=eMax;
    int n=10;
    hmc= new TH1F("hmc","",n,min,max); hmc->SetXTitle("E (keV)"); hmc->SetYTitle("Q_{HT}"); hmc->SetLineWidth(3);
    hmctruth= new TH1F("hmctruth","",n,min,max);  
    hmcwrong= new TH1F("hmcwrong","",n,min,max);
    if (!what.Contains("truth")) {
        res.Draw("segment.E()/27*1.3>>hmc",WimpCut,"goff");
        res.Draw("segment.E()/27*1.3>>hmctruth","","goff");
        res.Draw("segment.E()/27*1.3>>hmcwrong",WimpCut+"&&segment.skewness()>0","goff");
    } else {
        res.Draw("setup.emc>>hmc",WimpCut,"goff");
        res.Draw("setup.emc>>hmctruth","","goff");
        res.Draw("setup.emc>>hmcwrong",WimpCut+"&&segment.skewness()>0","goff");
    }
     if (!what.Contains("pdf")) {
         for (int i=1; i<=hmc->GetNbinsX(); i++) {
             float truth = hmctruth->GetBinContent(i);
             float reco  = hmc->GetBinContent(i);
             float wrong = hmcwrong->GetBinContent(i);
             float eff    = truth>0 ? reco/truth : 0;
             float efferr = truth>0 ? sqrt(eff*(1-eff)/truth) : 0;
             float w    = reco>0 ? wrong/reco : 0;
             float werr = reco>0 ? sqrt(w*(1-w)/reco) : 0;
             y = eff*pow(1-2*w,2);
             yerr = sqrt( pow(efferr*pow(1-2*w,2), 2) + pow(eff*2*(1-2*w)*2*werr, 2) );
             hmc->SetBinContent( i, y );
             hmc->SetBinError( i, yerr );
         }
    }
    hmc->Draw(opt);
    return hmc;
}

TF1*
WIMP_surface_bg() {
    // from NEWAGE
    float x[]={100,110,120, 200, 275, 400};
    float y[]={250,200,150, 50,   30, 20};
    gbg=new TGraph(6,x,y);
    gbg->Draw("APL");
    TF1 * fbg=new TF1("fbg","[0]*exp([1]*x)+[2]");
    gbg->Fit("fbg");
    return gbg->GetFunction("fbg");
}



void
WIMP_significance() {
    // take cosRecoil only, for now
    // add energy later
    int ngen= plotFit ? 1 : 1000;
    double nBackground=plotFit ? 1000 : backgroundRatePerKgYearKeV*(eMax-eMin)*years*kilos;
    double nSignal= plotFit ? 1000 : 0;
    float cosMin= hasHeadTail ? -1 : 0;
    bool keepData=true;
    
    RooRealVar *cosRecoil=new RooRealVar("cosRecoil", "cos#theta_{recoil}", 0, cosMin, 1);
    RooRealVar *eRecoil=new RooRealVar("eRecoil", "E_{recoil}", 100, eMin, eMax);

    // signal PDF
    WIMP_cosRecoil(); RooDataHist *hsig= new RooDataHist("hsig","signal mc", RooArgList(*eRecoil, *cosRecoil), hmc); 
    RooHistPdf *sigPdf=new RooHistPdf("sigPdf","signal PDF", RooArgSet(*eRecoil, *cosRecoil), *hsig);
    //WIMP_energy_Q("pdf"); RooDataHist *hsig= new RooDataHist("hsig","signal mc", RooArgList(*eRecoil), hmc); 
    //RooHistPdf *sigPdf=new RooHistPdf("sigPdf","signal PDF", RooArgSet(*eRecoil), *hsig);
    
    // background PDF
    //RooRealVar bgexpo("bgexpo","", -0.0166);
    //RooExponential *eRecoilBkgPdf=new RooExponential("eRecoilBkgPdf","", *eRecoil, bgexpo);    
    RooArgList bkgpars;
    RooPolynomial *eRecoilBkgPdf=new RooPolynomial("eRecoilBkgPdf","", *eRecoil, bkgpars, 0);
    RooPolynomial *cosRecoilBkgPdf=new RooPolynomial("cosRecoilBkgPdf","", *cosRecoil, bkgpars, 0);
    RooProdPdf *bkgPdf=new RooProdPdf("bkgPdf","Background PDF", *eRecoilBkgPdf, *cosRecoilBkgPdf );
    
    // model
    RooRealVar *nsig=new RooRealVar("nsig","Number of signal events", nSignal, -100, 100000);
    RooRealVar *nbkg=new RooRealVar("nbkg","Number of background events", nBackground, -100, 100000);

    RooAddPdf *model=new RooAddPdf ("model",  "dm search model", RooArgSet (*sigPdf, *bkgPdf), RooArgList(*nsig, *nbkg) );
    //RooAddPdf *model=new RooAddPdf ("model",  "dm search model", RooArgSet (*sigPdf, *eRecoilBkgPdf), RooArgList(*nsig, *nbkg) );
    
    RooMCStudy *toy=new RooMCStudy(*model, RooArgSet(*eRecoil, *cosRecoil), RooFit::Extended(), RooFit::FitOptions("rhem") );
    toy->generateAndFit(ngen, nSignal+nBackground, keepData);

    // plot 1st sample
    if (plotFit) {
        cp=new TCanvas; cp->Divide(2,1); 
        cp->cd(1); RooPlot *p1=toy->genData(0)->plotOn( cosRecoil->frame(10) ); model->plotOn( p1, RooFit::ProjWData(*toy->genData(0)) );  p1->Draw();
        cp->cd(2); RooPlot *p2=toy->genData(0)->plotOn( eRecoil->frame(10) ); model->plotOn( p2, RooFit::ProjWData(*toy->genData(0)) );  p2->Draw();
        return;
    }

    TH1F *hres=new TH1F("hres","",1000, -50, 50);
    for (int i=0; i<ngen; i++) {
        if (toy->fitResult(i)->covQual() != 3) continue;
        RooRealVar *fsig=toy->fitResult(i)->floatParsFinal().find("nsig");
        hres->Fill( fsig->getVal() );
    }
    hres->Draw();
    float nsig90CL=0;
    for (int i=1; i<= hres->GetNbinsX(); i++) {
        if (hres->Integral(1,i)/hres->Integral()<0.9) continue;
        nsig90CL=hres->GetBinCenter(i);
        cout << "Nsig (90% CL) = " << nsig90CL << " for " << hres->Integral() << " samples"<< endl;
        break;
    }
    TFeldmanCousins fc;
    cout << "Feldman-Cousins Nsig (90% CL) = " << fc.CalculateUpperLimit(nBackground, nBackground) << endl;
    cout << "Estimated total background = " << nBackground
         << " for energy range " << eMin<< "-"<< eMax
         << "  kilos="<<kilos
         << "  years="<<years
         << "  Mwimp="<<wimpMass
         << endl;
    
    // turn into cross section:
    double R = MaxCamWIMP::dRdE_ERange(eMin, eMax, wimpMass, targetIonMass); // for sigma=1pb, mass=1kg, days=1
    R *= kilos * (years*365);

    // form factor?
    // efficiency?
    
    
    // divide 90% CL with
    cout << "Estimated total WIMP rate for 1pb = " << R << endl;
    cout << "CS limit (90% CL) = " <<  nsig90CL/R << " pb" << endl;
    
}



void
WIMP_CS() {
    c= new TCanvas();
    gPad->SetLogy();
    gPad->SetLogx();
    TLatex l; l.SetTextSize(0.035);

    // needed to reach MSSM limit
    // 100 kg, ERecoil, cosRecoil, head-tail, 50keV threshold, 300m3 @ 75Torr, 0.01 bg/kg/y/keV
    float x_100kgy[]={  25,  30,   40,  50,   80,  100,  200,  400,  800};
    float y_100kgy[]={ 4.1, 2.8,  2.5, 2.4,  2.7,  3.1,  5.5,  10,  19.3};
    g_100kgy= new TGraph(9, x_100kgy, y_100kgy);
    for (int i=0;i<g_100kgy->GetN();i++) {
        double spinFactor = 0.75/0.674;
        double protonFactor = 1./0.46;
        double muSqFactor = pow( 1. / 19.* (19+x_100kgy[i])/(1+x_100kgy[i]) ,  2);// 
        g_100kgy->GetY()[i]*= 1e-3 * spinFactor * protonFactor * muSqFactor;
        cout << 1./ ( spinFactor * protonFactor * muSqFactor ) << endl;
    }
    g_100kgy->SetLineWidth(3);
    g_100kgy->GetHistogram()->SetXTitle("WIMP mass (GeV)");
    g_100kgy->GetHistogram()->SetYTitle("#sigma_{p} (pb)");
    g_100kgy->SetMinimum(1e-6);
    g_100kgy->SetMaximum(3);//000);
    g_100kgy->Draw("AL");
    l.DrawLatex(200,2.8e-3,"DMTPC, 100kg-y");



    
    // MSMM theory
    // Ellis, Ferstl, Olive, PRD63, Fig5c tanbeta=10
    float x_th[]={ 53,  60,  70,  79,   81, 100,  130,  135, 141,  156,  161,  171, 186,  194,  200};
    float y_th[]={0.1, 0.6, 2.5, 3.9, 4.02, 0.6, 0.36, 0.29, 0.3, 0.25, 0.25, 0.22, 0.12, 0.12, 0.1};
    g_th= new TGraph(15, x_th, y_th);
    for (int i=0;i<g_th->GetN();i++) g_th->GetY()[i]*=1e-3; 
    g_th->SetFillColor(7);
    //g_th->Draw("LF");
    //g_100kgy->Draw("L");
    //l.DrawLatex(70,2e-4,"MSSM");

    // 0.1kg-y, scale from 100kg-y
    g_1kgy= new TGraph(9, x_100kgy, y_100kgy);
    for (int i=0;i<g_1kgy->GetN();i++) {
        double spinFactor = 0.75/0.674;
        double protonFactor = 1./0.46;
        double muSqFactor = pow( 1. / 19.* (19+x_100kgy[i])/(1+x_100kgy[i]) ,  2);// 
        g_1kgy->GetY()[i] *= sqrt(100./0.1) * 1e-3 * spinFactor * protonFactor * muSqFactor;
    }
    g_1kgy->Draw("L");    g_1kgy->SetLineWidth(3); 
    l.DrawLatex(200,1e-1,"DMTPC, 0.1kg-y");


    // 0.001kg-y, scale from 100kg-y
    //g_10liter= new TGraph(9, x_100kgy, y_100kgy);
    //for (int i=0;i<g_1kgy->GetN();i++) g_10liter->GetY()[i] *= 1e-3*sqrt(100./0.004); 
    //g_10liter->Draw("L");    g_10liter->SetLineWidth(3); 
    //l.DrawLatex(200,1e-1,"DMTPC, 0.004kg-y");


    
    // needed for best directional search (surface)
    // 0.15 kg*day, ERecoil, cosRecoil, head-tail, 100-200keV threshold, 300m3 @ 75Torr, 2.7e4 bg/y/kg/keV
    //float x_150gday[]={/*    25,*/    30,     40,    50,     80,    100,    200,    400,    800};
    //float y_150gday[]={/* 55355,*/ 24411,  21183, 16832,  14269,  13113,  17531,  29328,  51410};
    //g_150gday= new TGraph(8, x_150gday, y_150gday);
    //g_150gday->SetLineWidth(3);  g_150gday->SetLineWidth(3); 
    //g_150gday->Draw("L");
    //l.DrawLatex(250,5e3,"DMTPC, 0.15kg-d");
    //l.DrawLatex(250,2e3,"(surface)");


    // COUPP
    float x_coupp[]={  4,   5,   6,   10,   20,   40,   60,   100,  200,  400,  700,  900, 1000};
    float y_coupp[]={ 50,   6, 2.3,  0.6, 0.36, 0.28, 0.28,  0.36, 0.59,    1,   1.7,   2,   2.3};
    g_coupp= new TGraph(12, x_coupp, y_coupp);
    g_coupp->SetLineWidth(2);  g_coupp->SetLineStyle(2); 
    g_coupp->Draw("C");
    l.SetTextAngle(15); l.SetTextSize(0.025);
    l.DrawLatex(400,1.10,"COUPP (2007)");


    // KIMS
    float x_kims[]={   20,   30,   40,   50,   60,   100,  200,  300,   400,  800,  1000};
    float y_kims[]={ 13.7, 0.36, 0.20, 0.17, 0.16,  0.17, 0.22, 0.28,  0.36,  0.7,   0.8};
    g_kims= new TGraph(11, x_kims, y_kims);
    g_kims->SetLineWidth(2);  g_kims->SetLineStyle(2); 
    g_kims->Draw("C");
    l.DrawLatex(400,0.38,"KIMS (2007)");
    
    // NAIAD
    float x_naiad[]={   20,   30,   40,   50,   60,   70,   80,   90,  100,  200,  300,    400,   800,  1000};
    float y_naiad[]={  0.5, 0.36, 0.30, 0.28, 0.25, 0.243, 0.25, 0.268, 0.28, 0.40, 0.52,  0.67,  1.30,  1.49};
    g_naiad= new TGraph(14, x_naiad, y_naiad);
    g_naiad->SetLineWidth(2);  g_naiad->SetLineStyle(2); 
    g_naiad->Draw("C");
    l.DrawLatex(400,0.69,"NAIAD (2005)");
    
}


//  20       30        40       50       80        110       
//  0.0011   0.00151   0.00224  0.0029   0.0079    0.0236
