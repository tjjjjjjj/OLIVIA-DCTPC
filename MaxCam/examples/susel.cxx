void tpc(TString opt="ito") {
    
    TCanvas *c1=0;
    if (opt.Contains("tpc")) {
        c1 = new TCanvas("c1", "c1",840,800);
        c1->Range(-220,-400, 320,400);
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

    TBox *shield=new TBox(-165,-320, 165,320); shield->SetLineWidth(3); shield->SetFillStyle(3247); shield->SetLineColor(1); shield->SetFillColor(1); shield->Draw();
    TBox *shield0=new TBox(-165,-320, 165,320); shield0->SetLineWidth(2); shield0->SetFillStyle(0);  shield0->Draw();
    TBox *vessel=new TBox(-150,-300, 150,300); vessel->SetLineWidth(3); vessel->SetFillStyle(1001); vessel->SetFillColor(10); vessel->Draw();
    TBox *vessel0=new TBox(-150,-300, 150,300); vessel0->SetLineWidth(2); vessel0->SetFillStyle(0); vessel0->Draw();

    TLine *feed1a=new TLine(165,270, 195, 270); feed1a->SetLineWidth(2); feed1a->Draw();
    TLine *feed1b=new TLine(200,270, 250, 270); feed1b->SetLineWidth(2); feed1b->Draw();
    TLine *feed1c=new TLine(165,277, 255, 277); feed1c->SetLineWidth(2); feed1c->Draw();
    TLine *feed2a=new TLine(200,270, 200, -270); feed2a->SetLineWidth(2); feed2a->Draw();
    TLine *feed2b=new TLine(195,270, 195, -270); feed2b->SetLineWidth(2); feed2b->Draw();
    TLine *feed3a=new TLine(165,-270, 195, -270); feed3a->SetLineWidth(2); feed3a->Draw();
    TLine *feed3b=new TLine(200,-270, 220, -270); feed3b->SetLineWidth(2); feed3b->Draw();
    TLine *feed3c=new TLine(165,-277, 220, -277); feed3c->SetLineWidth(2); feed3c->Draw();
    TLine *feed4a=new TLine(250,270, 250, 50); feed4a->SetLineWidth(2); feed4a->Draw();
    TLine *feed4b=new TLine(255,277, 255, 50); feed4b->SetLineWidth(2); feed4b->Draw();
    TPaveText *bottle=new TPaveText(238,-50, 268,50); bottle->SetLineWidth(2); bottle->SetFillStyle(0); bottle->SetBorderSize(1);
    bottle->AddText("CF_{4}"); bottle->AddText("gas"); bottle->AddText("bottle"); bottle->AddText("(24 lb)"); bottle->Draw();
    TPaveText *pump=new TPaveText(220,-295, 260,-255); pump->SetLineWidth(2); pump->SetFillStyle(0); pump->SetBorderSize(1);
    pump->AddText("vacuum"); pump->AddText("pump");  pump->Draw();
    
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

   
   if (opt.Contains("tpc")) {
       text->SetTextSize(0.02);
       
       // drift mesh
       TLine *driftU = new TLine(-105,  (11+24*10), 105,  (11+24*10)); driftU->SetLineStyle(2); driftU->SetLineWidth(2);  driftU->Draw();
       TLine *driftD = new TLine(-105, -(11+24*10), 105, -(11+24*10)); driftD->SetLineStyle(2); driftD->SetLineWidth(2); driftD->Draw();

       // scale
       TArrow *driftscale1 = new TArrow(120, 5, 120, 250, 0.01, "<||>"); driftscale1->SetLineWidth(1); driftscale1->Draw();
       text->DrawLatex(-35, 59, "23 cm");
       TArrow *driftscale2 = new TArrow(120, -5, 120, -250, 0.01, "<||>"); driftscale2->SetLineWidth(1); driftscale2->Draw();   
       text->DrawLatex(122, 115, "25 cm");
       TArrow *ring = new TArrow(-90, 55, 90, 55, 0.01, "<||>"); ring->SetLineWidth(1); ring->Draw();   
       text->DrawLatex(122, -155, "25 cm");

       
       // aluminum rings
       for (int i=0; i<25; i++) {
           TPave *meshUL = new TPave(-105, 6+i*10,-92.2, 11+i*10, 0,"br"); meshUL->SetFillColor(12); meshUL->Draw();
           TPave *meshLL = new TPave(-105,-6-i*10,-92.2,-11-i*10, 0,"br"); meshLL->SetFillColor(12); meshLL->Draw();
           TPave *meshUR = new TPave( 105, 6+i*10, 92.2, 11+i*10, 0,"br"); meshUR->SetFillColor(12); meshUR->Draw();
           TPave *meshLR = new TPave( 105,-6-i*10, 92.2,-11-i*10, 0,"br"); meshLR->SetFillColor(12); meshLR->Draw();
       }   
       
       // resistors
       for (int i=0; i<25; i++) {
           TLine* wireLR=new TLine(-112.5, 9+i*10, -105,   9+i*10); wireLR->Draw();
           TLine* wireUD=new TLine(-112.5, i==0 ? 0 : 8+i*10, -112.5, i==24 ? 49 : 10+i*10); wireUD->Draw();
           wireLR=new TLine(-112.5, -9-i*10, -105,-9-i*10); wireLR->Draw();
           wireUD=new TLine(-112.5, i==0 ? 0 : -8-i*10, -112.5, i==24 ? -49 : -10-i*10); wireUD->Draw();
           if (i>23) continue;
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
       //TText *anodeHV_text = new TText(107, -1, "+V"); anodeHV_text->SetTextSize(0.06); anodeHV_text->Draw();
       //TText *driftHV_text = new TText(107, -50, "-V"); driftHV_text->SetTextSize(0.06); driftHV_text->Draw();
       //driftHV_text = new TText(107, 44, "-V"); driftHV_text->SetTextSize(0.06); driftHV_text->Draw();
       


       // lens
       //TArc *lens=new TArc(0, 0, 65, 80, 100); lens->SetNoEdges(); lens->SetFillStyle(0); lens->Draw();
       //lens=new TArc(0, 128, 65, 260, 280); lens->SetNoEdges(); lens->SetFillStyle(0); lens->Draw();

       // Camera
       //TBox *ccd=new TBox(-20, 66, 20, 75); ccd->SetFillStyle(1001); ccd->SetFillColor(14); ccd->Draw();
       //text->DrawLatex(-50,72,"Camera");
       //text->DrawLatex(-50,67,"+ Lens");
      


       //TLine *wimpin = new TLine(-200, 250, -20, 140); wimpin->SetLineWidth(3); wimpin->Draw();   
       //TLine *wimpout = new TLine(-20, 140, 200, 150); wimpout->SetLineWidth(3); wimpout->Draw();   
       
   }


   
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
   c1->ToggleToolBar();
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
    for (int i=0;i<g_100kgy->GetN();i++) g_100kgy->GetY()[i]*=1e-3; 
    g_100kgy->SetLineWidth(3);
    g_100kgy->GetHistogram()->SetXTitle("WIMP mass (GeV)");
    g_100kgy->GetHistogram()->SetYTitle("#sigma_{p} (pb)");
    g_100kgy->SetMinimum(1e-3);
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
    for (int i=0;i<g_1kgy->GetN();i++) g_1kgy->GetY()[i] *= 1e-3*sqrt(100./0.1); 
    g_1kgy->Draw("L");    g_1kgy->SetLineWidth(3); 
    l.DrawLatex(200,1e-1,"DMTPC, 0.1kg-y");


    // 0.001kg-y, scale from 100kg-y
    g_10liter= new TGraph(9, x_100kgy, y_100kgy);
    for (int i=0;i<g_1kgy->GetN();i++) g_10liter->GetY()[i] *= 1e-3*sqrt(100./0.004); 
    g_10liter->Draw("L");    g_10liter->SetLineWidth(3); 
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





float eMin=50;
float eMax=250;
float years=1.;
float kilos=100;
float wimpDensity=0.3; // GeV/cm3
float wimpMass=100; // GeV
float targetIonMass=19; // GeV
float pWIMPCS=0.1;


TString WimpCut="segment.n()>4";




TH1F*
WIMP_QHT_total(TString what="", TString opt="") {
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

    float min=eMin, max=eMax;
    int n=10;
    hmc= new TH1F("hmc","",n,min,max); hmc->SetXTitle("E (keV)"); hmc->SetYTitle("N_{R}"); hmc->SetLineWidth(3);
    hmctruth= new TH1F("hmctruth","",n,min,max);  
    hmcwrong= new TH1F("hmcwrong","",n,min,max);
    heff= new TH1F("hceff","",n,min,max);
    hQ= new TH1F("hQ","",n,min,max);hQ->SetYTitle("Q_{HT}"); 

    res.Draw("segment.E()/27*1.3>>hmc",WimpCut,"goff");
    res.Draw("segment.E()/27*1.3>>hmctruth","","goff");
    res.Draw("segment.E()/27*1.3>>hmcwrong",WimpCut+"&&segment.skewness()>0","goff");

    float totalQ=0;
    for (int i=1; i<=hmc->GetNbinsX(); i++) {
        float truth = hmctruth->GetBinContent(i);
        float reco  = hmc->GetBinContent(i);
        float wrong = hmcwrong->GetBinContent(i);
        float eff    = truth>0 ? reco/truth : 0;
        float efferr = truth>0 ? sqrt(eff*(1-eff)/truth) : 0;
        float w    = reco>0 ? wrong/reco : 0;
        float werr = reco>0 ? sqrt(w*(1-w)/reco) : 0;
        float Q = eff*pow(1-2*w,2);
        float Qerr = sqrt( pow(efferr*pow(1-2*w,2), 2) + pow(eff*2*(1-2*w)*2*werr, 2) );
        hQ->SetBinContent( i, Q );
        hQ->SetBinError( i, Qerr );
        heff->SetBinContent( i, eff );
        heff->SetBinError( i, efferr  );
        totalQ += Q*reco;
    }
    double R = MaxCamWIMP::dRdE_ERange(eMin, eMax, wimpMass, targetIonMass); // for sigma=1pb, mass=1kg, days=1
    R *= kilos * (years*365) * pWIMPCS;
    cout << "Estimated total WIMP recoil rate = " << R << endl
         << "Energy range: " << eMin << "-" << eMax << "keV" << endl
         << "CF4 mass = " << kilos << "kg" << endl
         << "Duration = " << years << "years" << endl
         << "p-cross section = " << pWIMPCS << "pb" << endl;
    float scale=R/hmctruth->Integral();
    cout << "WIMP rate scaling = " << scale << endl;
    hmc->Scale(scale);
    totalQ *= scale;
    cout << "Total Q = " << totalQ << " significance = " << sqrt(totalQ) << "  P = " << 1-TMath::Prob(totalQ,1) << endl;
    
    
    hmc->Draw(opt);
    return hmc;
}



void plotCrossSections() {

   MaxCamENDF *he = new MaxCamENDF("ENDF_CS_n_on_4He", "CS");
   MaxCamENDF *c = new MaxCamENDF("ENDF_CS_n_on_12C", "CS");
   MaxCamENDF *f = new MaxCamENDF("ENDF_CS_n_on_19F", "CS");

   for (int i=0; i<he->getCrossSection()->GetN(); i++) he->getCrossSection()->GetX()[i]*=0.001; //->MeV
   for (int i=0; i<c->getCrossSection()->GetN(); i++) c->getCrossSection()->GetX()[i]*=0.001; //->MeV
   for (int i=0; i<f->getCrossSection()->GetN(); i++) f->getCrossSection()->GetX()[i]*=0.001; //->MeV

   frame= new TH2F("frame","", 100, 0.1, 4000, 100, 0.05, 150); frame->Draw(); gPad->SetLogy(); gPad->SetLogx();
   frame->SetXTitle("Neutron energy (MeV)");
   frame->SetYTitle("Cross section (b)");
   f->getCrossSection()->Draw("L");   f->getCrossSection()->SetLineColor(1); f->getCrossSection()->SetLineWidth(3); f->getCrossSection()->SetLineStyle(1);
   he->getCrossSection()->Draw("L"); he->getCrossSection()->SetLineColor(2); he->getCrossSection()->SetLineWidth(3); he->getCrossSection()->SetLineStyle(2);
   c->getCrossSection()->Draw("L");   c->getCrossSection()->SetLineColor(4); c->getCrossSection()->SetLineWidth(3); c->getCrossSection()->SetLineStyle(6);



   leg= new TLegend(0.5, 0.5, 0.8,0.7);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->AddEntry( he->getCrossSection(), "He(n,n)He", "L");  
   leg->AddEntry( c->getCrossSection(),  "C(n,n)C",   "L");  
   leg->AddEntry( f->getCrossSection(),  "F(n,n)F",   "L");  
   leg->Draw();
}


void plotAngles( float energy=1e3, TString opt="AL") {
   MaxCamENDF *he = new MaxCamENDF("ENDF_DCS_n_on_4He", "scatt");
   TGraph2D *g=he->getDifferentialCrossSection();
   float x[1000], y[1000];
   int n=0;
   for (float c=0.99999; c>-0.99999; c-=0.01) {
       x[n]=acos(c)*180./3.14;
       y[n]=g->Interpolate(energy+1e-5, c);
       //if (y[n]<0) break;
       //cout << x[n]<<"  "<<y[n]<<endl;
       n++;
   }
   TGraph *gr=new TGraph(n,x,y);
   gr->Draw(opt);
}



TH1D* calcNeutronSpectrum(
    TString fname="MonteCarloSamples/mcrun_000154.root",
    float ionMass=4,
    int   nevents=3875) {
    
    
    TChain ntp("ntp");
    ntp.Add(fname);
    
    int n=16;
    float neutronMaxE=8000;
    float r=4*1*ionMass/pow(ionMass+1,2);
    TH2F *h2=new TH2F("h2","", n,0,neutronMaxE, n,0,neutronMaxE*r);
    ntp.Draw("er:en>>h2","","goff");
    TH1D *en=h2->ProjectionX("en");
    TH1D *er=h2->ProjectionY("er");
    TH1D *enreco = new TH1D("enreco","", n,0,neutronMaxE);
    float x[100], y[100];
    for (int i=n; i>0; i--) {
        x[i-1]=er->GetBinCenter( i );
        y[i-1]=h2->GetBinContent(i, i)>0 ? en->GetBinContent( i )/h2->GetBinContent(i, i) : 0;        
    }

    
    ntp.Draw("er:en>>h2","er>50","goff", nevents);
    TH1D *er=h2->ProjectionY("er");
    for (int i=n; i>0; i--) {
        double nreco = er->GetBinContent( i ) * y[i-1];
        double nrecoErr = er->GetBinError( i ) * y[i-1]; // assume scaling error is small (no systematics)
        enreco->SetBinContent( i, nreco );
        enreco->SetBinError( i, nrecoErr );

        //cout << x[i-1] << "  " << y[i-1] << "  " <<  er->GetBinContent( i ) << " -> " << nreco <<endl;

        for (int j=1; j<n; j++) {
            double newer=er->GetBinContent(j)-h2->GetBinContent(i, j);
            er->SetBinContent( j, newer );
            //cout << newer <<", ";
        }
        //cout << endl;
    }

    
    //g= new TGraph(n,x,y); g->Draw("APL");
    //enreco->Draw();

    return enreco;
}

void neutronSpectrum() {
    TH1D *he=calcNeutronSpectrum("MonteCarloSamples/mcrun_000154.root", 4, 4325); he->SetXTitle("Neutron energy (keV)");
    cout << "he="<<he->Integral() << endl;
    TH1D *f=calcNeutronSpectrum("MonteCarloSamples/mcrun_000155.root", 19, 1767); f->SetMarkerStyle(25);
    cout << "f="<<f->Integral() << endl;
    
    he->Draw(); f->Draw("same"); 
    leg= new TLegend(0.55, 0.5, 0.85,0.65);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry( he, "He(n,n)He", "P");  
    leg->AddEntry( f,  "F(n,n)F",   "P");  
    leg->Draw();
}


TH1D* calcNeutronHT(TString fname, int nevents) {
    
    
    TChain ntp("ntp");
    ntp.Add(fname);
    
    int n=5;
    TH1F* hpos=new TH1F("hpos","",n,-TMath::PiOver2(),TMath::PiOver2());
    TH1F *hneg=new TH1F("hneg","",n,-TMath::PiOver2(),TMath::PiOver2());
    TH1D *hasy=new TH1D("hasy","",n,-TMath::PiOver2(),TMath::PiOver2());
    ntp.Draw("atan(recoilY/recoilX)>>hpos","recoilX>0&&er>100","goff", nevents);
    ntp.Draw("atan(recoilY/recoilX)>>hneg","recoilX<0&&er>100","goff", nevents);
    for (int i=1; i<=n; i++) {
        double tot= hpos->GetBinContent(i)+hneg->GetBinContent(i);
        double asym= tot>0 ? (( hpos->GetBinContent(i)-hneg->GetBinContent(i) )/tot) : 0;
        double asymerr=tot ? (sqrt(  hpos->GetBinContent(i)+hneg->GetBinContent(i) )/tot) : 0;
        hasy->SetBinContent(i,  asym);
        hasy->SetBinError(i, asymerr);   
    }

    return hasy;
}


void htasymmetry() {

    TH1D *he=calcNeutronHT("MonteCarloSamples/mcrun_000156.root", 4325);
    
    TH1D *f=calcNeutronHT("MonteCarloSamples/mcrun_000157.root", 1767); f->SetMarkerStyle(25);

    TH2F *hframe=new TH2F("hframe","",100,-TMath::PiOver2(),TMath::PiOver2(),  100, -0.3,1);
    hframe->SetXTitle("Neutron direction (rad)"); hframe->SetYTitle("Head-Tail Asymmetry");
    hframe->Draw();
    he->Draw("same");
    f->Draw("same");
    leg= new TLegend(0.55, 0.5, 0.85,0.65);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry( he, "He(n,n)He", "P");  
    leg->AddEntry( f,  "F(n,n)F",   "P");  
    leg->Draw();
}


TH1D* calcRecoilEnergy(TString fname, int nevents) {
    TChain ntp("ntp");
    ntp.Add(fname);
    int n=20;
    TH1D *her = new TH1D("her","er>50",n,0,5000);
    ntp.Draw("er>>her","","goff", nevents);
    return her;
}


void erecoil() {
   TH1D *he=calcRecoilEnergy("MonteCarloSamples/mcrun_000156.root", 4325);
   he->SetXTitle("Recoil energy (keV)");
   TH1D *f =calcRecoilEnergy("MonteCarloSamples/mcrun_000157.root", 1767); f->SetMarkerStyle(25);
   he->Draw("e");
   f->Draw("e same");
   leg= new TLegend(0.55, 0.5, 0.85,0.65);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry( he, "He(n,n)He", "P");  
    leg->AddEntry( f,  "F(n,n)F",   "P");
    leg->Draw();
}


double muonRate(float h) {
    // differential 1/s/cm2/sr
//     double I1=8.6e-6;  // 1/s/cm2/sr
//     double I2=0.44e-6; // 1/s/cm2/sr
//     double l1=450; // m.w.e.
//     double l2=870; // m.w.e

    // flat Earth -> 1/cm2/s
    double I1=67.97e-6;  // 1/s/cm2/sr
    double I2=2.071e-6; // 1/s/cm2/sr
    double l1=285; // m.w.e.
    double l2=698; // m.w.e

    return I1*exp(-h/l1) + I2*exp(-h/l2);
}

double neutronsFromMuonRate(float h) {
    // rate in n/cm2/s
    double P0=0.4e-6;
    double P1=860;//mwe
    return P0*P1/h*exp(-h/P1);
}
