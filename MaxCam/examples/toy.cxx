void toyRecoil() {

        MaxCamMC mc;
        mc.fillSrimTable("/Users/ddujmic/Desktop/DCH/SRIM/tables/F_in_CF4_100Torr"); // SRIM table for F -> CF4
        mc.setStopping(180); // 180 Torr

        c1= new TCanvas("c1");

        TFile f("toy_tmp.root","RECREATE"); // file to save output
        TH2F *hs=new TH2F("hs","",10,0,1000, 50,-1,1);

        for (double e=100; e<1000; e+=100) { // recoil energy

        for (int i=0; i<100; i++) { // numer of recoils per energy

                double phi=2*TMath::Pi()*gRandom->Rndm(); // random phi angle
                mc.makeRecoil(e,phi); // compute kinematics of recoil, set recoil track
                mc.setRecoilCoord(0,0,0); // define recoil vertex
                mc.event(); // track recoil through detector applying diffusion and ccd resolution
                mc.applyBiasADC(mc.getWireImage()); // smear ccd image due to adc bias
                TH1D *hwire1=mc.getWireImage()->ProjectionX("x1",161,171 ); // these are know positions of wires
                TH1D *hwire0=mc.getWireImage()->ProjectionX("x0",120,130 ); //        - ii -
                hwire0->Add(hwire1); // sum signal

                //hwire0->Draw();gSystem->Sleep(1000);c1->Update();mc.getWireImage()->Draw("colz"); gSystem->Sleep(1000); c1->Update();
                double skewness=hwire0->GetSkewness();//mc.getTrackImage()->ProjectionX("x")->GetSkewness();
                cout <<skewness<<endl;

                hs->Fill(e,skewness);

        }
        }
        hs->Write();
        f.Write();

}


void toyAlpha() {

        MaxCamMC mc;
        mc.fillSrimTable("/Users/ddujmic/Desktop/DCH/SRIM/tables/He_in_CF4_100Torr"); // SRIM table for He -> CF4
        mc.setStopping(180); // 180 Torr
        mc.setWireImage(96,0,12,  64,-1.4,6.4);

        // 'recoil' track will be propagted through the chamber so        
        // set alpha mass, energy 
        mc.setRecoil(-500,0,0,4e6); 
        mc.setRecoilCoord(10,4,0); // define starting point for alpha
        mc.event(); // propagate through detector
        mc.applyBiasADC(mc.getWireImage()); // smear ccd image due to adc bias
 
        MaxCamImageTools::resizeImage(mc.getWireImage(), 0,768, 0,512)->Draw("colz");
}



void compareAlphaFluorine() {

        float theta=64.4;
        theta*=TMath::Pi()/180.;

        MaxCamMC mc;
        mc.fillSrimTable("/Users/ddujmic/Desktop/DCH/SRIM/tables/He_in_CF4_100Torr"); // SRIM table for He -> CF4
        mc.setStopping(180); // 180 Torr
        mc.setWireImage(96,0,12,  64,-1.4,6.4);
        mc.setRecoil(-500*cos(theta),500*sin(theta),0,4e6); 
        mc.setRecoilCoord(10,5,0); // define starting point for alpha
        mc.event(); // propagate through detector

        mc.fillSrimTable("/Users/ddujmic/Desktop/DCH/SRIM/tables/F_in_CF4_100Torr"); // SRIM table for F -> CF4
        mc.setStopping(180); // 180 Torr
        mc.setRecoil(-500*cos(theta),500*sin(theta),0,19e6); 
        mc.setRecoilCoord(10,0,0); // define starting point for F
        mc.event(false); // propagate through detector

        mc.applyBiasADC(mc.getWireImage()); // smear ccd image due to adc bias 
        mc.getWireImage()->Draw("colz"); 
}