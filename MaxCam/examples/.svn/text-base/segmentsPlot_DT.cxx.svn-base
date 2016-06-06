TString cut0_data="wire0.y>wire1.y&&wire0.widthR-wire0.widthL>0&&wire0.y<1e5&&wire0.widthL>12&&wire0.widthR<83"; // 12 83
TString cut1_data="wire0.y<wire1.y&&wire1.widthR-wire1.widthL>0&&wire1.y<1e5&&wire1.widthL>12&&wire1.widthR<83";
TString cut0_mc  ="wire0.y>wire1.y&&wire0.widthR-wire0.widthL>0&&wire0.y<1e5&&wire0.widthL>10&&wire0.widthR<85";
TString cut1_mc  ="wire0.y<wire1.y&&wire1.widthR-wire1.widthL>0&&wire1.y<1e5&&wire1.widthL>10&&wire1.widthR<85";

TChain data("res");       
TChain data180("res");       
TChain mcF("res");        
TChain mcHe("res");       

TProfile* hdata, *hmcHe, *hmcF;
TH2F *sframe;
TH1F *hstat, *hstat_neg;
TLegend *leg;
void loadRuns() {
    data.SetMarkerStyle(24); data.SetMarkerColor(1);
    data180.SetMarkerStyle(20); data.SetMarkerColor(1);
    mcF.SetMarkerStyle(7);
    mcHe.SetMarkerStyle(7); mcHe.SetMarkerColor(3);

    // problems with pressure measurement
    // h: eps, taup
    // i: 250Torr
    // j: 200Torr

    //data.Add("RESULTS/results_run139h.root");
    //data.Add("RESULTS/results_run146h.root");
    //data180.Add("RESULTS/results_run149h.root");

    //data.Add("segout.root");

    // 90deg runs
    //data.Add("RESULTS/results_run185g.root");
    //data.Add("RESULTS/results_run186g.root");
    //data.Add("RESULTS/results_run187g.root");
    //data.Add("RESULTS/results_run188g.root");
   
    // zero aperture in concrete wall
    //data.Add("RESULTS/results_run147.root");
    
    
    //mcF.Add("RESULTS/results_mcrun_00044.root");
    mcF.Add("RESULTS/results_mcrun_00054.root");
    

    
    //float x[]={10,30, 50, 70, 90, 110, 130, 150}; int nx=7;
    //float x[]={300,600,900, 1500, 2100}; int nx=4;
    //float x[]={100,400,700, 1000, 1300, 1600, 1900}; int nx=6;
    
    //float x[]={0,5000}; int nx=1;
    float x[]={5,10,20,30,40,50}; int nx=5;
    for (int i=0; i<=nx; i++) x[i]/=8;
    //float x[]={0,15,30,45}; int nx=3;
    hdata=new TProfile("hdata","", nx,   x, "S"); hdata->SetMarkerStyle(20);  hdata->SetMarkerColor(1); hdata->SetFillColor(5); hdata->SetFillStyle(4100);
    hmcHe=new TProfile("hmcHe","", nx-2, x, "S"); hmcHe->SetMarkerStyle(7);  hmcHe->SetMarkerColor(3);  hmcHe->SetLineColor(3); 
    hmcHe->SetFillStyle(4100); hmcHe->SetFillColor(3);
    hmcF=new TProfile("hmcF","", nx, x, "S");   hmcF->SetMarkerStyle(1); hmcF->SetFillColor(11); hmcF->SetLineColor(1); hmcF->SetFillStyle(4100);
    sframe= new TH2F("sframe","", 2,x[0]-50,x[nx-1]+50, 2,-1,1);

    leg=new TLegend(0.5,0.6,0.8,0.8);
    leg->AddEntry(hdata, "data", "p");
    leg->AddEntry(hmcF, "MC", "lf");
    leg->SetBorderSize(0);

    hstat = new TH1F("hstat","", nx,x);
    hstat_neg = new TH1F("hstat_neg","", nx,x);
}


void stopping(TString what="", TString opt="same") {
        loadRuns();

        TString wire0_mc  ="wire0.y:(wire0.widthR-wire0.widthL)";
        TString wire1_mc  ="wire1.y:(wire1.widthR-wire1.widthL)";
        //TString wire0_data  ="wire0.y/8:wire0.esrim";
        //TString wire1_data  ="wire1.y/8:wire1.esrim";
        TString wire0_data="wire0.y/8:(wire0.widthR-wire0.widthL)/8.25";
        TString wire1_data="wire1.y/8:(wire1.widthR-wire1.widthL)/8.25";
 
        if (what.Contains("prof")){
            wire0_data+=">>+hdata";
            wire1_data+=">>+hdata";
            wire0_mc+=">>+hmcF"; 
            wire1_mc+=">>+hmcF"; 
        }

        
        h= new TH2F("h","",50,0,50, 50,0,20000); h->Draw();

        
        data.Draw(wire0_data, cut0_data, what);
        data.Draw(wire1_data, cut1_data, opt);
        data180.Draw(wire0_data, cut0_data, opt);
        data180.Draw(wire1_data, cut1_data, opt);

        mcF.Draw(wire0_mc, cut0_mc, what);
        mcF.Draw(wire1_mc, cut1_mc, opt);

        //mcHe.Draw(wire0_mc, cut0_mc, opt);
        //mcHe.Draw(wire1_mc, cut1_mc, opt);

        //data.Draw(wire0_data, cut0_data, opt);
        //data.Draw(wire1_data, cut1_data, opt);
}


void skewness(TString what="prof", TString opt="same") {

        loadRuns();
        
        TString wire0_data="wire0.skewness:(wire0.widthR-wire0.widthL)/8"; if (what.Contains("prof")) wire0_data+=">>+hdata";
        TString wire1_data="wire1.skewness:(wire1.widthR-wire1.widthL)/8"; if (what.Contains("prof")) wire1_data+=">>+hdata";
        //TString wire0_data="wire0.skewness:wire0.esrim"; if (what.Contains("prof")) wire0_data+=">>+hdata";
        //TString wire1_data="wire1.skewness:wire1.esrim"; if (what.Contains("prof")) wire1_data+=">>+hdata";

        TString wire0_mcF="wire0.skewness:(wire0.widthR-wire0.widthL)/8";  if (what.Contains("prof")) wire0_mcF+=">>+hmcF";
        TString wire1_mcF="wire1.skewness:(wire1.widthR-wire1.widthL)/8";  if (what.Contains("prof")) wire1_mcF+=">>+hmcF";
        //TString wire0_mcF="wire0.skewness:wires.emc";  if (what.Contains("prof")) wire0_mcF+=">>+hmcF";
        //TString wire1_mcF="wire1.skewness:wires.emc";  if (what.Contains("prof")) wire1_mcF+=">>+hmcF";

        TString wire0_mcHe="wire0.skewness:wire0.widthR-wire0.widthL"; if (what.Contains("prof")) wire0_mcHe+=">>+hmcHe";
        TString wire1_mcHe="wire1.skewness:wire1.widthR-wire1.widthL"; if (what.Contains("prof")) wire1_mcHe+=">>+hmcHe";

        h= new TH2F("h","",50,0,50/8., 50,-1,1); h->Draw(); h->SetYTitle("Skewness"); h->SetXTitle("Projected range (mm)");
        
        data.Draw(wire0_data, cut0_data, what); 
        data.Draw(wire1_data, cut1_data, opt);
        data180.Draw(wire0_data, cut0_data, opt); 
        data180.Draw(wire1_data, cut1_data, opt); 

        mcF.Draw(wire0_mcF, cut0_mc, what);
        mcF.Draw(wire1_mcF, cut1_mc, opt);

        //mcHe.Draw(wire0_mcHe, cut0_mc, opt);
        //mcHe.Draw(wire1_mcHe, cut1_mc, opt);

        if (what.Contains("prof")) {
            h->Draw();
            hmcF->Draw("e2 h same");
            hdata->Draw("same");
        }
}

void emc() {
    loadRuns();
    TString wire0_mcF="wire0.y*0.04:wires.emc"; 
    TString wire1_mcF="wire1.y*0.04:wires.emc";  
    mcF.Draw(wire0_mcF, cut0_mc, "box");
    mcF.Draw(wire1_mcF, cut1_mc, "box same"); 
}

void ecalib() {
    loadRuns();
    TString wire0_data="wire0.y/8:wire0.esrim"; 
    TString wire1_data="wire1.y/8:wire1.esrim";  
    data.Draw(wire0_data, cut0_data,"");
    data.Draw(wire1_data, cut1_data,"same");
    data180.Draw(wire0_data, cut0_data,"same");
    data180.Draw(wire1_data, cut1_data,"same");
}


TGraphAsymmErrors *g, *gmc;
void statSkew() {
    loadRuns();

     TString wire0_data="(wire0.widthR-wire0.widthL)/8";
     TString wire1_data="(wire1.widthR-wire1.widthL)/8";
//      TString wire0_data="wire0.y/20";
//      TString wire1_data="wire1.y/20";
    
     data.Draw(wire0_data+">>+hstat_neg", cut0_data+"&&wire0.skewness<0");
     data.Draw(wire1_data+">>+hstat_neg", cut1_data+"&&wire1.skewness<0");
    data.Draw(wire0_data+">>+hstat", cut0_data);
    data.Draw(wire1_data+">>+hstat", cut1_data);
    
    TGraphAsymmErrors *g=MaxCamStat::fractionError(hstat_neg, hstat);
    g->Draw("AP");
    g->GetHistogram()->SetYTitle("Fraction with #gamma<0");
    g->GetHistogram()->SetXTitle("Projected range (mm)");
}

void mcstatSkew(TString what="asym") {
        loadRuns();
  
        hstat = new TH1F("hstat","", 10,10,210);
        hstat_neg = new TH1F("hstat_neg","", 10, 10, 210);

        mcF.Draw("wires.emc>>+hstat_neg", cut0_mc+"&&wire0.skewness<0", "goff");
        mcF.Draw("wires.emc>>+hstat", cut0_mc, "goff");
        
        if (what=="asym") {
            gmc=MaxCamStat::asymmetryError(hstat_neg, hstat);
            gmc->Draw("AP");
            //gmc->GetHistogram()->SetYTitle("Fraction with #gamma<0");
            gmc->GetHistogram()->SetYTitle("Head-tail asymmetry");
            gmc->GetHistogram()->SetXTitle("Recoil energy (keV)");
        } else {
            gmc=MaxCamStat::fractionError(hstat_neg, hstat);
            gmc->Draw("AP");
            gmc->GetHistogram()->SetYTitle("Fraction with #gamma<0");
            gmc->GetHistogram()->SetXTitle("Recoil energy (keV)");
        }
            
}
