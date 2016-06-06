#include "MaxCamSRIM.hh"
#include "TMath.h"
#include "TFile.h"
#include "TSystem.h"
#include "TH2.h"
#include <iostream>
#include <fstream>
#include <vector>
using std::cout;
using std::endl;
using std::cerr;
using std::vector;

#include <string>
using std::string;
#include <sstream>
using std::istringstream;



ClassImp(MaxCamSRIM)

//____________________
//
MaxCamSRIM::MaxCamSRIM( const char *fname, TString energyUnitString )  : _targetDensityAt600Torr(0), _projectileMass(0) {
  // Constructor uses name of file with MC and a pointer to 
  // the random number generator.
  
   TFile *df=new TFile("current.root","RECREATE"); df->GetName();
   _srimTable = new TNtuple("srimTable","","E:dEdx_ele:dEdx_nuc:range:lstrag:tstrag:dR:dLStraggle:dTStraggle");

   // default energy units are keV, default length units are mm
   SetEnergyUnits(energyUnitString);

   fillSrimTable( fname ); 
   setStopping(600); 
}


MaxCamSRIM::MaxCamSRIM(const MaxCamSRIM &other) {
        _stopping=other._stopping;
        cout << GetName() << "Copy Constructor not done" << endl;
}

MaxCamSRIM::~MaxCamSRIM() {
  delete _stopping;
  delete _stoppingElec;
  delete _range;
  delete _rangeElec;
  delete _rangenergy;
  delete _rangenergyElec;
  delete _enerange;
  delete _enerangeElec;
  
  // added 4/10/2010 by Shawn; deletes TGraphs of the raw numbers in the SRIM file
  delete _rawtotalstopping;
  delete _rawelectricstopping;
  delete _rawnuclearstopping;
  delete _rawrangenergy;
}



void
MaxCamSRIM::setStopping(float p) {
        // Convert stopping power in MeV/(mg/cm^2) to keV/mm for given
        // pressure. Member getStopping()->Eval(energy) returns stopping power for
        // this this pressure at given energy.
        //

    f=new TFile("$DCTPC_TOP_DIR/MaxCam/tables/10MeV_He_recoil_in_600torr_125CF4_875He.root");
    hist = (TH2D*)f->Get("hr");

    _press = p;
    //float factor = density(p)*1e5; // convert to keV/mm
    float fac = (600./_press); //spitz original code: float fac = (100./_press) changed to 600 for 600 torr
    //spitz: it looks like this expects a SRIM input file which is made with 100 torr. Divide by 6 because Big and little DCTPC nominally have 600 torr. This is going to lead to confusion at some point!
    float factor = (_targetDensityAt600Torr/6.)*1e5*_press/100; // convert to keV/mm
        
        cout<<"factor "<<(_targetDensityAt600Torr/6.)*1e5*_press/100<<endl;
        
    vector<float> t,u,v,w,x,y,z,s,ry,rr,rw,rq;
    int n=_srimTable->GetEntries();
    for (int i=0; i<n; i++) {
        _srimTable->GetEvent(i);
        float *vars=_srimTable->GetArgs();
        x.push_back( vars[0] ); // energy
        y.push_back( (vars[1]+vars[2])*factor ); // total dE/dx
	ry.push_back( (vars[1]+vars[2]) ); // raw total dE/dx
        w.push_back( vars[1]*factor ); // elec dE/dx
        rw.push_back( vars[1] ); // raw elec dE/dx
        rr.push_back( vars[2] ); // raw nucl dE/dx
	rq.push_back( vars[3] ); // raw expected range
        u.push_back( vars[4]*fac ); // long. straggling  
        v.push_back( vars[5]*fac ); // tran. straggling
        s.push_back( vars[6] );
//         u.push_back( vars[7]*fac ); // long. straggling  
//         v.push_back( vars[8]*fac ); // tran. straggling
//         cout << vars[0] << "  " << (vars[1]+vars[2])*factor << endl;
    }
    _stopping = new TGraph(n, &x[0], &y[0]);
    assert(_stopping);
    _stoppingElec = new TGraph(n, &x[0], &w[0]);
    assert(_stoppingElec);
    _lstraggle = new TGraph(n, &x[0], &u[0]);
    _tstraggle = new TGraph(n, &x[0], &v[0]);
    _stepSize  = new TGraph(n, &x[0], &s[0]);

    // added 4/10/2010 by Shawn; fills the TGraphs of the raw numbers in the SRIM file
    _rawtotalstopping = new TGraph(n, &x[0], &ry[0]);
    _rawelectricstopping = new TGraph(n, &x[0], &rw[0]);
    _rawnuclearstopping = new TGraph(n, &x[0], &rr[0]);

    _rawrangenergy = new TGraph(n, &x[0], &rq[0]);

    t.clear();
    w.clear();
    x.clear();
    y.clear();
    z.clear();
    double dstep=0.05; // 50um step
    double range=0;
    double dEdx=0, dEdx_elec=0;
    double ionEnergy=0.009; // start from 9eV
    double ionEnergy_elec = 0.5*ionEnergy; // assume 1/2 of energy goes to scintillation
    double maxEnergy=_stopping->GetX()[_stopping->GetN()-1];
    while ( ionEnergy<maxEnergy ) {
        dEdx = _stopping->Eval( ionEnergy );
        dEdx_elec = _stoppingElec->Eval( ionEnergy );
        ionEnergy += dEdx*dstep;
        ionEnergy_elec += dEdx_elec*dstep;
        range += dstep;
        t.push_back( ionEnergy_elec );
        y.push_back( dEdx );
        x.push_back( range );
        z.push_back( ionEnergy );
        w.push_back( dEdx_elec );
    }        
    _range = new TGraph(x.size(), &x[0], &y[0]);
    assert(_range);
    _rangenergy = new TGraph(z.size(), &z[0], &x[0]);
    assert(_rangenergy);
    _enerange = new TGraph(z.size(), &x[0], &z[0]);
    assert(_enerange);
    
    _rangeElec = new TGraph(x.size(), &x[0], &w[0]);
    assert(_rangeElec);
    _rangenergyElec = new TGraph(t.size(), &t[0], &x[0]);
    assert(_rangenergyElec);
    _enerangeElec = new TGraph(t.size(), &x[0], &t[0]);
    assert(_enerangeElec);
}



double
MaxCamSRIM::calcEnergyLoss(double ionEnergy, double distance0, double distance1, bool isTotal, double stepsize) {
        // Compute the energy loss for given ion starting energy
        // between two spacial points.

        double dstep=stepsize/1000.;
        double range=0;
        double sumE=0, dE=0, dElec;
        while (range<distance1 && ionEnergy>0) {
                dE = _stopping->Eval( ionEnergy )*dstep;
		dElec = _stoppingElec->Eval( ionEnergy )*dstep;
                if ( range>distance0 && range<distance1 ) sumE += isTotal ? dE : dElec ;
                range += dstep;
                ionEnergy -= dE;
                //cout << dE << " -> " << ionEnergy << "  " << range << endl;
        }
        return sumE;
}








TGraph*
MaxCamSRIM::readSRIM(const char* fileName, TString what, float minE, float maxE, float dEdxFactor) {
        // OBSOLETE FUNCTION.Reads SRIM tables and makes a graph object for given relation.
        // Inputs are file name with SRIM table, qualtities to be related, energy range and
        // a conversion factor to desired units. 
        //
        cout <<"readSRIM is OBSOLETE, please use fillSrimTable !!!"<<endl; assert(0);

    vector<double> x;
    vector<double> y;
    double Energy, Range, dEdxElec, dEdxNucl, zz;
    double ElecSum=0;
    double NuclSum=0;
    double elecFrac=1;
    double oldRange=-1;
    ifstream fin(fileName);
    char buff[256];
    TString unit;	
    // read header
    for (int i=0; i<26; i++) fin.getline(buff, 256);
    while (!fin.eof()) {
    	fin >> Energy >> unit;
	if (unit=="eV")  Energy*=0.001;
	if (unit=="MeV") Energy*=1000;
	if (maxE>0 && Energy>maxE) break;
	if (minE>0 && Energy<minE) continue;
	fin >> dEdxElec >> dEdxNucl >> Range >> unit >> zz >> buff >> zz >> buff ;
	if (unit=="um") Range*=0.001;

	//cout << dEdxElec <<"  "<< dEdxNucl <<"  "<< Range <<"  "<< unit <<"  "<< zz <<"  "<< buff <<"  "<< zz <<"  "<< buff <<endl;;
        // integrate Elec and Nucl contributions to extract fraction of energy lost due to Elec processes
        elecFrac=1;
        if (what.Contains("elec") && oldRange>0) {
                ElecSum += dEdxElec * (Range-oldRange) * dEdxFactor * 1000; // keV
                NuclSum += dEdxNucl * (Range-oldRange) * dEdxFactor * 1000; // keV
                elecFrac = ElecSum/(ElecSum+NuclSum);
        }
        oldRange=Range;
        //cout << ElecSum << "  " << NuclSum << "   " << ElecSum+NuclSum << "   " << Energy << endl;

//        if (what.Contains("xproject")) {
//                Range *= calcCosRecoil( Energy ); // account for finite recoil angle
//        }


	if (what.Contains("range") && what.Contains("energy")) {
   		x.push_back(Energy*elecFrac);
    		y.push_back(Range);
		cout << Energy << "  " << Range << endl;
	}
	else if (what.Contains("energy") && what.Contains("dEdx")) {
                x.push_back(Energy);
                double dEdx= what.Contains("elec") ? dEdxElec : (dEdxElec + dEdxNucl);
                dEdx *= dEdxFactor;
    		y.push_back( dEdx );
		cout << Energy << "  " << dEdx << endl;
	}
	else if (what.Contains("range") && what.Contains("dEdx")) {
                x.push_back(Range);
                double dEdx=  what.Contains("elec") ? dEdxElec : (dEdxElec + dEdxNucl);
                dEdx *= dEdxFactor;
    		y.push_back( dEdx );
		cout << Range << "  " << dEdx << endl;
	}

	
    }
    assert( x.size()==y.size());
    if (x.size()<1) return 0;
    return new TGraph( x.size(), &x[0], &y[0]);
}




void
MaxCamSRIM::fillSrimTable(const char* fileName) {
        // Load SRIM table into memory. dE/dx is given in units of density so
        // before using this table, one has to specify working pressure - see
        // function setStopping(pressure) for details.

    const char* fileDir=gSystem->Getenv("MCTABLES");
    if (!fileDir) {
        cout << GetName() <<": Environment variable MCTABLES not defined." << endl; 
        assert(0);
    }    
    TString fullName(fileDir);
    fullName += "/";
    fullName += fileName;

    _srimTable->Reset();
    //cout << fullName << endl;

    double Energy, Range, dEdxElec, dEdxNucl, LStraggle, TStraggle;
    double delRange, delLStraggle, delTStraggle;
    double oldRange=0, oldLStraggle=0, oldTStraggle=0;
    ifstream fin((const char*)fullName);
    assert( fin.is_open() );
    TString unit;	
    string line;
    // read header
    while (1) {
        getline( fin, line);
        if (line.find("Ion =")<line.size()) {
            istringstream slineion(line);
            slineion >> unit >> unit >> unit >> unit >> unit >> unit >> unit >> _projectileMass;
            _projectileMass *= 1e6*0.93*GetEnergyConversionFactor(); // keV if GetEnergyConversionFactor() returns 1
            getline( fin, line);
            getline( fin, line);
            istringstream slinedensity(line);
            slinedensity >> unit >> unit >> unit >> _targetDensityAt600Torr; // all normalized to 100 Torr
	    break;
        }
    }
    while (1) {
        getline( fin, line);
        if (line.find("-----------  ----------")<line.size()) break;
    }   
    while (!fin.eof()) {
        getline( fin, line);
        if (line.find("---")<line.size()) break;
        istringstream sline(line);
    	sline >> Energy >> unit;
	if (unit=="eV")  Energy*=0.001*GetEnergyConversionFactor();
	if (unit=="keV")  Energy*=1.*GetEnergyConversionFactor();
	if (unit=="MeV") Energy*=1e3*GetEnergyConversionFactor();
	if (unit=="GeV") Energy*=1e6*GetEnergyConversionFactor();
	sline >> dEdxElec >> dEdxNucl >> Range >> unit;
	if (unit=="m")  Range*=1e3;
	if (unit=="um") Range*=1e-3;
	if (unit=="A") Range*=1.e-7;
        sline >> LStraggle >> unit;
	if (unit=="m")  LStraggle*=1e3;        
	if (unit=="um") LStraggle*=1e-3;        
	if (unit=="A")  LStraggle*=1e-7;        
        sline  >> TStraggle >> unit ;
	if (unit=="m")  TStraggle*=1e3;        
	if (unit=="um") TStraggle*=1e-3;        
	if (unit=="A")  TStraggle*=1e-7;        

        delRange = Range - oldRange;
        delLStraggle = sqrt( LStraggle*LStraggle - oldLStraggle*oldLStraggle );
        delTStraggle = sqrt( TStraggle*TStraggle - oldTStraggle*oldTStraggle );
        oldRange=Range;
        oldLStraggle=delLStraggle;
        oldTStraggle=delTStraggle;
        
        _srimTable->Fill( Energy, dEdxElec, dEdxNucl, Range, LStraggle, TStraggle, delRange, delLStraggle, delTStraggle);
    }
    // load the conversions

    _srimTable->SetTitle( fileName );
}


void MaxCamSRIM::Print() {

    cout << "Table Name ....... " << _srimTable->GetTitle() << endl;
    cout << "Pressure (Torr)... " << _press << endl;
}


double
MaxCamSRIM::density(double pTorr, double A, double T) {
        // Calculate density from pressure, molecular mass, temparature.
        

    double amu=1.66e-24;
    return amu*A*numberDensity(pTorr,T); //  g/cm^3
}

double
MaxCamSRIM::getdedx(double energy) {
        // Calculate density from pressure, molecular mass, temparature.

 
    int bin=hist->GetXaxis()->FindBin(energy);    
    proj = hist->ProjectionY("",bin,bin);
    double dedx=proj->GetRandom();
    return dedx; //  
}


double
MaxCamSRIM::numberDensity(double pTorr, double T) {
        // Calculate density from pressure, molecular mass, temparature.
        
    double k = 1.38e-23; //  J/K
    double P = pTorr*101325./760;
    double rho = P/(k*T);
    return rho*1e-6; //  1/cm^3
}


double
MaxCamSRIM::X0(int A, int Z) {
    // Helper function to return radiation length for gammas
    // using Dahl approximation (in g/cm2)
    
    return 716.4 * A / ( Z*(Z+1)*log(287/sqrt(Z)) );
}

void
MaxCamSRIM::SetEnergyUnits(TString input_string)
{ 
  // right now only MeV and keV are supported

  // we're already returning energy in the requested units!
  if( ( input_string=="MeV" ) ||
      ( input_string=="keV" ) ){

    _energyUnitString=input_string;

  } else { 
    cout << "requested energy unit \"" << input_string << "\" not supported!  Giving you energies in keV!" << endl;
    _energyUnitString="keV";
  }
  return;
}

float
MaxCamSRIM::GetEnergyConversionFactor(void){
  // code is pre-tooled for returning keV energies
  if( _energyUnitString=="keV" ) return 1.;
  if( _energyUnitString=="MeV" ) return 1./1000.; // convert all energies to MeV!
  // we shouldn't ever get to this return, but just do energies in keV
  // if we do...
  return 1.;
}

Double_t
MaxCamSRIM::ROOTv24TGraphEval(TGraph* gr,Double_t x){
  // MERCILESSLY STOLEN FROM ROOT v24'S TGraph MEMBER FUNCTION Eval(...)
  //linear interpolation
  //find point in graph immediatly below x
  //In case x is < fX[0] or > fX[fNpoints-1] return the extrapolated point
  Int_t low = TMath::BinarySearch(gr->GetN(),gr->GetX(),x);
  Int_t up = low+1;
  if (low == gr->GetN()-1) {up=low; low = up-1;}
  if (low == -1) {low=0; up=1;}
  if (gr->GetX()[low] == gr->GetX()[up]) return gr->GetY()[low];
  Double_t yn = x*(gr->GetY()[low]-gr->GetY()[up]) +gr->GetX()[low]*gr->GetY()[up] - gr->GetX()[up]*gr->GetY()[low];
  return yn/(gr->GetX()[low]-gr->GetX()[up]);
}


