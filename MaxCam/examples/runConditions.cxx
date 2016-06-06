
float get_max_dEdX_times_multiplication() {
	MaxCamSRIM srimHe("SRIM_He_in_CF4_100Torr");
	srimHe.setPressure(75);
	TGraph* dEdXele = srimHe.getStoppingVsEnergy(0);
	float maxdEdX = 0;
	for (int i = 0; i < dEdXele->GetN(); i++) {
		float dEdX = dEdXele->GetY()[i];
		if (dEdX > maxdEdX) maxdEdX = dEdX;
	}
	//cout << maxdEdX;
	// value 38 obtained from paper 0804.4827v1
	return maxdEdX * 38;
}

//float max_dEdX_times_multiplication=70*20; // keV/mm * ADU/keV 
float max_dEdX_times_multiplication=get_max_dEdX_times_multiplication();

TH2F *F_dedx, *F_npho, *F_sgnf;
TH1F *F_ADU;

float vixel=0.596;
float pixelNoise=6.6;

std::vector<TString> string_vector = 0;
TObjArray* dedx_vector = new TObjArray(4);
TObjArray* adu_vector = new TObjArray(4);
TObjArray* npho_vector = new TObjArray(4);
TObjArray* sgnf_vector = new TObjArray(4);

/*
void drawdEdX() {
	MaxCamSRIM srimHe("SRIM_He_in_CF4_100Torr");
	srimHe.setPressure(75);
	TGraph* dEdX = srimHe.getStoppingVsEnergy(false);
	TCanvas* c = new TCanvas("c", "Stopping Power");
	c->cd();
	dEdX->SetMarkerStyle(7);
	dEdX->Draw("AP");
}*/

void fillHisto( TH2F *h, TH1F* hap, TH2F *hr, TH2F* hsgf, TString fname) {
    MaxCamSRIM stop(fname);
    int nx=h->GetXaxis()->GetNbins();
    int ny=h->GetYaxis()->GetNbins();
    for (int i=1; i<=nx; i++) {
	float dummy_val = h->GetXaxis()->GetBinCenter(i);
        stop.setPressure( h->GetXaxis()->GetBinCenter(i) );
        
        // stopping power vs pressure
        float max_dedx_elec=-1;
        for (int j=1; j<=ny; j++) {
            float dedx_elec = stop.getStoppingVsEnergy(false)->Eval( h->GetYaxis()->GetBinCenter(j) );
            h->SetBinContent(i,j, dedx_elec );
            if (dedx_elec>max_dedx_elec) max_dedx_elec=dedx_elec;            
        }
        
        // maximum charge gain vs pressure
        float maxADU = max_dEdX_times_multiplication/max_dedx_elec;
        //cout << max_dEdX_times_multiplication << "  " << max_dedx_elec << "  " << maxADU << endl;
        hap->SetBinContent(i, maxADU);
   
        // photon yield vs pressure
        for (int j=1; j<=ny; j++) {
            float dedx_elec = stop.getStoppingVsEnergy(false)->Eval( h->GetYaxis()->GetBinCenter(j) );
            float range = stop.getRangeVsEnergy()->Eval( h->GetYaxis()->GetBinCenter(j) );
            float npho=range*dedx_elec*maxADU;
            hr->SetBinContent(i,j, npho );

            float noise=pixelNoise*sqrt(range/vixel);
            hsgf->SetBinContent(i,j, npho/noise);
        }
        

    }
}

void build() {
	string_vector.push_back("F");
	string_vector.push_back("C");
	string_vector.push_back("He");
	string_vector.push_back("H");

	TString dedx_name("");
	TString adu_name("");
	TString npho_name("");
	TString sgnf_name("");

	for (int i = 0; i < string_vector.size(); i++) {
		TString prefix = string_vector.at(i);

		cout << "Building " << prefix << "\n";
	
		dedx_name = prefix + "_dedx";
		adu_name = prefix + "_ADU";
		npho_name = prefix + "_npho";
		sgnf_name = prefix + "_sgnf";

		h_dedx = new TH2F(dedx_name, "dE/dx (keV/mm)", 10, 0, 200, 10, 0, 500);
		h_dedx->SetXTitle("Pressure (Torr)");
		h_dedx->SetYTitle("Energy (keV)");

		TString hname="Spark threshold = ";
		hname += max_dEdX_times_multiplication;
		hname += " ADU/mm";

		h_ADU = new TH1F(adu_name, hname, 10, 0, 200);
		h_ADU->SetXTitle("Pressure (Torr)");
		h_ADU->SetYTitle("Max Gain (ADU/keV)");

		h_npho = new TH2F(npho_name, "Radiated photons", 10, 0, 200, 10, 0, 500);
		h_npho->SetXTitle("Pressure (Torr)");
		h_npho->SetYTitle("Energy (keV)");

		h_sgnf = new TH2F(sgnf_name, "Significance", 10, 0, 200, 10, 0, 500);
		h_sgnf->SetXTitle("Pressure (Torr)");
		h_sgnf->SetYTitle("Energy (keV)");

		TString fileName("SRIM_");
		fileName += prefix;
		fileName += "_in_CF4_100Torr";

//		cout << "Filling histograms for " << prefix << " with " << fileName << "\n";
		fillHisto(h_dedx, h_ADU, h_npho, h_sgnf, fileName);
//		cout << "Done filling histogram.\n";

		dedx_vector->Add(h_dedx);
		adu_vector->Add(h_ADU);
		npho_vector->Add(h_npho);
		sgnf_vector->Add(h_sgnf);
	}
}

void plot() {
	build();
	plot("F");
}

TH1* plotTest() {
	cout << "Size: " << dedx_vector->GetEntries() << "\n";
	return (TH1*) dedx_vector->At(0);
}

int getIndexOf(TString what) {
	int index = -1;
	for (int i = 0; i < string_vector.size(); i++) {
		if (what == string_vector.at(i)) {
			index = i;
		}
	} 

	if (index >= string_vector.size()) {
		cout << what << " not built.";
		return -1;
	} else return index;
}

void plot(TString what) {

	int index = getIndexOf(what);
	if (index == -1) return;
	cout << "Found " << what << " at index " << index << "\n";

	TH2F* dedx = dedx_vector->At(index);
	TH1F* ADU = adu_vector->At(index);
	TH2F* npho = npho_vector->At(index);
	TH2F* sgnf = sgnf_vector->At(index);

	cout << "Assigned histograms.\n";

	gStyle->SetPalette(1);
	gStyle->SetOptTitle(1);

	c = new TCanvas("c", what);
	c->Divide(2,2);
	c->cd(1);  dedx->SetStats(kFALSE); dedx->Draw("colz");
	c->cd(2);  ADU->SetStats(kFALSE); ADU->Draw();
	c->cd(3);  npho->SetStats(kFALSE); npho->Draw("colz");
	c->cd(4);  sgnf->SetStats(kFALSE); sgnf->Draw("colz");

}

void plot2(TString what) {
	if (what == "F") {
		F_dedx = new TH2F("F_dedx", "dE/dx (keV/mm)", 10, 0, 200,   10, 0, 500);
		F_dedx->SetXTitle("Pressure (Torr)");
		F_dedx->SetYTitle("Energy (keV)");

		TString hname="Spark threshold =";
		hname += max_dEdX_times_multiplication;
		hname += " ADU/mm";
		F_ADU = new TH1F("F_ADU", hname, 10, 0, 200);
		F_ADU->SetXTitle("Pressure (Torr)");	
		F_ADU->SetYTitle("Max Gain (ADU/keV)");
    
		F_npho = new TH2F("F_npho", "Radiated photons", 10, 0, 200,   10, 0, 500);
		F_npho->SetXTitle("Pressure (Torr)");
		F_npho->SetYTitle("Energy (keV)");    
    
		F_sgnf = new TH2F("F_sgnf", "Significance", 10, 0, 200,   10, 0, 500);
		F_sgnf->SetXTitle("Pressure (Torr)");
		F_sgnf->SetYTitle("Energy (keV)");

		fillHisto( F_dedx, F_ADU, F_npho, F_sgnf, "SRIM_F_in_CF4_100Torr");

 		sgnf_vector->Add(F_sgnf);        

		gStyle->SetPalette(1);
		gStyle->SetOptTitle(1);
		c=new TCanvas("c","Fluorine");
		c->Divide(2,2);
		c->cd(1);  F_dedx->SetStats(kFALSE); F_dedx->Draw("colz");
		c->cd(2);  F_ADU->SetStats(kFALSE); F_ADU->Draw();
		c->cd(3);  F_npho->SetStats(kFALSE); F_npho->Draw("colz");
		c->cd(4);  F_sgnf->SetStats(kFALSE); F_sgnf->Draw("colz");
	}
	else if (what == "C") {
                C_dedx = new TH2F("C_dedx", "dE/dx (keV/mm)", 10, 0, 200,   10, 0, 500);
                C_dedx->SetXTitle("Pressure (Torr)");
                C_dedx->SetYTitle("Energy (keV)");

                TString hname="Spark threshold =";
                hname += max_dEdX_times_multiplication;
                hname += " ADU/mm";
                C_ADU = new TH1F("C_ADU", hname, 10, 0, 200);
                C_ADU->SetXTitle("Pressure (Torr)");
                C_ADU->SetYTitle("Max Gain (ADU/keV)");

                C_npho = new TH2F("C_npho", "Radiated photons", 10, 0, 200,   10, 0, 500);
                C_npho->SetXTitle("Pressure (Torr)");
                C_npho->SetYTitle("Energy (keV)");

                C_sgnf = new TH2F("C_sgnf", "Significance", 10, 0, 200,   10, 0, 500);
                C_sgnf->SetXTitle("Pressure (Torr)");
                C_sgnf->SetYTitle("Energy (keV)");

                fillHisto( C_dedx, C_ADU, C_npho, C_sgnf, "SRIM_C_in_CF4_100Torr");

		sgnf_vector->Add(C_sgnf);

		gStyle->SetPalette(1);
                gStyle->SetOptTitle(1);
                c=new TCanvas("c","Carbon");
                c->Divide(2,2);
                c->cd(1);  C_dedx->SetStats(kFALSE); C_dedx->Draw("colz");
                c->cd(2);  C_ADU->SetStats(kFALSE); C_ADU->Draw();
                c->cd(3);  C_npho->SetStats(kFALSE); C_npho->Draw("colz");
                c->cd(4);  C_sgnf->SetStats(kFALSE); C_sgnf->Draw("colz");
	}

	else if (what == "H" || what == "p") {
		H_dedx = new TH2F("H_dedx", "dE/dx (keV/mm)", 10, 0, 200,   10, 0, 500);
                H_dedx->SetXTitle("Pressure (Torr)");
                H_dedx->SetYTitle("Energy (keV)");

                TString hname="Spark threshold =";
                hname += max_dEdX_times_multiplication;
                hname += " ADU/mm";
                H_ADU = new TH1F("H_ADU", hname, 10, 0, 200);
                H_ADU->SetXTitle("Pressure (Torr)");
                H_ADU->SetYTitle("Max Gain (ADU/keV)");

                H_npho = new TH2F("H_npho", "Radiated photons", 10, 0, 200,   10, 0, 500);
                H_npho->SetXTitle("Pressure (Torr)");
                H_npho->SetYTitle("Energy (keV)");

                H_sgnf = new TH2F("H_sgnf", "Significance", 10, 0, 200,   10, 0, 500);
                H_sgnf->SetXTitle("Pressure (Torr)");
                H_sgnf->SetYTitle("Energy (keV)");

                fillHisto( H_dedx, H_ADU, H_npho, H_sgnf, "SRIM_H_in_CF4_100Torr");

		sgnf_vector->Add(H_sgnf);

		gStyle->SetPalette(1);
                gStyle->SetOptTitle(1);
                c=new TCanvas("c","Proton");
                c->Divide(2,2);
                c->cd(1);  H_dedx->SetStats(kFALSE); H_dedx->Draw("colz");
                c->cd(2);  H_ADU->SetStats(kFALSE); H_ADU->Draw();
                c->cd(3);  H_npho->SetStats(kFALSE); H_npho->Draw("colz");
                c->cd(4);  H_sgnf->SetStats(kFALSE); H_sgnf->Draw("colz");
	
	}

	else if (what == "He") {
		He_dedx = new TH2F("C_dedx", "dE/dx (keV/mm)", 10, 0, 200,   10, 0, 500);
                He_dedx->SetXTitle("Pressure (Torr)");
                He_dedx->SetYTitle("Energy (keV)");

                TString hname="Spark threshold =";
                hname += max_dEdX_times_multiplication;
                hname += " ADU/mm";
                He_ADU = new TH1F("C_ADU", hname, 10, 0, 200);
                He_ADU->SetXTitle("Pressure (Torr)");
                He_ADU->SetYTitle("Max Gain (ADU/keV)");

                He_npho = new TH2F("C_npho", "Radiated photons", 10, 0, 200,   10, 0, 500);
                He_npho->SetXTitle("Pressure (Torr)");
                He_npho->SetYTitle("Energy (keV)");

                He_sgnf = new TH2F("C_sgnf", "Significance", 10, 0, 200,   10, 0, 500);
                He_sgnf->SetXTitle("Pressure (Torr)");
                He_sgnf->SetYTitle("Energy (keV)");

                fillHisto( He_dedx, He_ADU, He_npho, He_sgnf, "SRIM_He_in_CF4_100Torr");

		sgnf_vector->Add(He_sgnf);

		gStyle->SetPalette(1);
                gStyle->SetOptTitle(1);
                c=new TCanvas("c","Helium");
                c->Divide(2,2);
                c->cd(1);  He_dedx->SetStats(kFALSE); He_dedx->Draw("colz");
                c->cd(2);  He_ADU->SetStats(kFALSE); He_ADU->Draw();
                c->cd(3);  He_npho->SetStats(kFALSE); He_npho->Draw("colz");
                c->cd(4);  He_sgnf->SetStats(kFALSE); He_sgnf->Draw("colz");
	
	}

}

TH1D* getSignificanceVsPressure(TString what) {
	//TH2F* sigvpres = new TH2F("sigvpres", "Significance v. Pressure", 10, 0, 200, 10, 0, 500);

	if (what == "all") return getSigVPressureAll();

	int index = getIndexOf(what);
	if (index == -1) return NULL;
	else {
		TH2F* hsgnf = sgnf_vector->At(index);
		TString pname(what);
		pname+= "_sgnfvP";
		return hsgnf->ProjectionX(pname, 0, 10, "");
	}
}

TH1D* getSignificanceVsEnergy(TString what) {
	//TH2F* sigvpres = new TH2F("sigvpres", "Significance v. Pressure", 10, 0, 200, 10, 0, 500);

	if (what == "all") return getSigVEnergyAll();

	int index = getIndexOf(what);
	if (index == -1) return NULL;
	else {
		TH2F* hsgnf = sgnf_vector->At(index);
		TString pname(what);
		pname+= "_sgnfvE";
		return hsgnf->ProjectionY(pname, 0, 10, "");
	}
}

void drawSigVEnergyAll() {
	TString prefix("");
	Double_t max = 0.0;
	TObjArray* arr2 = new TObjArray(string_vector.size());
	for (int i = 0; i < string_vector.size(); i++) {
		prefix = string_vector.at(i);
//		int index = getIndexOf(prefix);
		int color = getColorOf(prefix);
		TH1D* sigve = getSignificanceVsEnergy(prefix);
		if (sigve->GetMaximum() > max) max = sigve->GetMaximum();
		sigve->SetLineColor(color);
		arr2->Add(sigve);
		//if (i == 0) sigvp->Draw();
		//else sigvp->Draw("same");
	}
	TCanvas* cSvE = new TCanvas("cSvE", "Significance vs. Energy");
	for (int i = 0; i < arr2->GetEntries(); i++) {
		if (i == 0) {
			((TH1D*) (arr2->At(i)))->SetMaximum(max);
			((TH1D*) (arr2->At(i)))->Draw();
		}
		else {
			((TH1D*) (arr2->At(i)))->Draw("same");
		}
	}
}

void drawSigVPressureAll() {
	TString prefix("");
	Double_t max = 0.0;
	TObjArray* arr = new TObjArray(string_vector.size());
	for (int i = 0; i < string_vector.size(); i++) {
		prefix = string_vector.at(i);
//		int index = getIndexOf(prefix);
		int color = getColorOf(prefix);
		TH1D* sigvp = getSignificanceVsPressure(prefix);
		if (sigvp->GetMaximum() > max) max = sigvp->GetMaximum();
		sigvp->SetLineColor(color);
		arr->Add(sigvp);
		//if (i == 0) sigvp->Draw();
		//else sigvp->Draw("same");
	}
	TCanvas* cSvP = new TCanvas("cSvP", "Significance vs. Pressure");
	for (int i = 0; i < arr->GetEntries(); i++) {
		if (i == 0) {
			((TH1D*) (arr->At(i)))->SetMaximum(max);
			((TH1D*) (arr->At(i)))->Draw();
		}
		else {
			((TH1D*) (arr->At(i)))->Draw("same");
		}
	}
}

TH1D* getDiffSignificance(TString compareWhat, TString what1, TString what2) {
	Bool_t compareEnergy = false;
	if (compareWhat == "energy" || compareWhat == "E" || compareWhat == "Energy") compareEnergy = true;

	TH1D* sigvX_1 = 0;
	TH1D* sigvX_2 = 0;
	if (compareEnergy) {
		sigvX_1 = getSignificanceVsEnergy(what1);
		sigvX_2 = getSignificanceVsEnergy(what2);
	} else {
		sigvX_1 = getSignificanceVsPressure(what1);
		sigvX_2 = getSignificanceVsPressure(what2);
	}

/*
	sigvX_1->Fit("pol4", "R", "goff", 0, 200);
	TList *fns_1 = sigvX_1->GetListOfFunctions();
	TF1* fit_1 = (TF1*) fns_1->FindObject("pol4");
	fit_1->SetName("fit_1");

	sigvX_2->Fit("pol4", "R", "goff", 0, 200);
	TList *fns_2 = sigvX_2->GetListOfFunctions();
	TF1* fit_2 = (TF1*) fns_2->FindObject("pol4");
	fit_2->SetName("fit_2");

	//Double_t diff0 = fit_2->GetParameter(0) - fit_1->GetParameter(0);
	//Double_t diff1 = fit_2->GetParameter(1) - fit_1->GetParameter(1);
	//Double_t diff2 = fit_2->GetParameter(2) - fit_1->GetParameter(2);
	//Double_t diff3 = fit_2->GetParameter(3) - fit_1->GetParameter(3);

	TF1* ratFcn = new TF1("ratFcn", "(pol4(0)) / (pol4(4))", 0, 200);
	ratFcn->SetParameter(0, fit_2->GetParameter(0));
	ratFcn->SetParameter(1, fit_2->GetParameter(1));
	ratFcn->SetParameter(2, fit_2->GetParameter(2));
	ratFcn->SetParameter(3, fit_2->GetParameter(3));
	ratFcn->SetParameter(4, fit_1->GetParameter(0));
	ratFcn->SetParameter(5, fit_1->GetParameter(1));
	ratFcn->SetParameter(6, fit_1->GetParameter(2));
	ratFcn->SetParameter(7, fit_1->GetParameter(3));

	return ratFcn;
*/
	TString title("Significance of ");
	title += what2;
	title += " over significance of";
	title += what1;
	TH1D* hRatio = new TH1D("hRatio", title, 10, 0, 200);
	for (int i = 0; i <= hRatio->GetNbinsX(); i++) {
		Double_t first = sigvX_1->GetBinContent(i);
		Double_t second = sigvX_2->GetBinContent(i);
		if (first != 0) {
			Double_t ratio = second/first;
			hRatio->SetBinContent(i, ratio);
		}
	}

	return hRatio;

}

Int_t getColorOf(TString what) {
	if (what == "He") return 3;
	else if (what == "C") return 5;
	else if (what == "F") return 4;
	else if (what == "H") return 6;
}
