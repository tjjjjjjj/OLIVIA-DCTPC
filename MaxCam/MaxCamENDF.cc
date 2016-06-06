#include "MaxCamENDF.hh"
#include "MaxCamTwoBodyKinematics.hh"
#include "TMath.h"
#include "TFile.h"
#include "TSystem.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1D.h"
#include "TRandom.h"

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

#define legmax 14

double fLegendre(double *x, double *var) {
        double cosTheta=x[0];
        return  MaxCamTwoBodyKinematics::Legendre(cosTheta, var, legmax); 
}


ClassImp(MaxCamENDF)

//____________________
//
    MaxCamENDF::MaxCamENDF( const char *fileName, TString opt)  {
  // Constructor uses name of file with an ENDF table.
  //
  
  _fDCS = new TF1("_fDCS",fLegendre, -1,1,legmax);
  opt.ToLower();
  if (opt.Contains("scat"))      fillScatteringTable( fileName );
  else if (opt.Contains("fis"))  fillFissionTable( fileName );
  else if (opt.Contains("cs"))  fillCSTable( fileName );
  else assert(0);
}



MaxCamENDF::MaxCamENDF(const MaxCamENDF &other) {
        _fDCS=other._fDCS;
        _DCS=other._DCS;
        cout << GetName() << "Copy Constructor not done" << endl;
}


MaxCamENDF::~MaxCamENDF() {
}



void
MaxCamENDF::fillScatteringTable(const char* fileName) {
        // Load ENDF table into memory. 

        _DCS = new TGraph2D();
        
        const char* fileDir=gSystem->Getenv("MCTABLES");
        if (!fileDir) {
                cout << GetName() <<": Environment variable MCTABLES not defined." << endl; 
                assert(0);
        }    
        TString fullName(fileDir);
        fullName += "/";
        fullName += fileName;

        ifstream fin( fullName );

        string line;
        double E;
        int NW, n=0;
        int  exponent;
        for (int i=0; i<5; i++) { 
                getline( fin, line); //cout << line << endl; 
        }
        while (!fin.eof()) {
                getline( fin, line);
                istringstream sline(line);
                sline >> E >> exponent >> E >> exponent;
                sline >> NW >> NW >> NW;
                if (!NW) break;
                E *= pow(10, exponent-3); // keV

                //cout << "E=" << E << "   NW="<<NW<<endl;
                for (int i=0; i<_fDCS->GetNpar(); i++) _fDCS->SetParameter(i, 0);
                double P;
                istringstream *pline=0;
                _fDCS->SetParameter(0, 1);
                for (int i=0; i<NW; i++) {
                    if (i<NW) {
                        if (i%6==0) { if (pline) delete pline; getline(fin, line); pline = new istringstream(line);  }
                        (*pline) >> P >> exponent;
                        P *= pow(10, exponent);
                        _fDCS->SetParameter(i+1, P);
                    }
                }
                //for (int i=0; i<=NW; i++) cout << _fDCS->GetParameter(i) <<", "; cout << endl;
                //_fDCS->Print();
                
                // build points for tgraph2d object 
                double calcDCS;
                _maxCS=0;
                for (double cosTheta=-1; cosTheta<1.01; cosTheta+=0.1) {
                        calcDCS = _fDCS->Eval( cosTheta );
                        if (calcDCS>_maxCS) _maxCS=calcDCS;
                        //cout << calcDCS << "  " ;    
                        _DCS->SetPoint(n++, E, cosTheta, calcDCS);
                }
                //cout << endl;
        }

}




void
MaxCamENDF::fillCSTable(const char* fileName) {
        // Load ENDF table into memory. 

        _CS = new TGraph();

        const char* fileDir=gSystem->Getenv("MCTABLES");
        if (!fileDir) {
                cout << GetName() <<": Environment variable MCTABLES not defined." << endl; 
                assert(0);
        }    
        TString fullName(fileDir);
        fullName += "/";
        fullName += fileName;

        ifstream fin( fullName );

        string line;
        double E, CS;
        int  exponent, n=0, ntot=0;
        for (int i=0; i<2; i++) { 
            getline( fin, line); //cout << line << endl;
                
        }
        getline( fin, line); 
        istringstream nline(line);
        nline >> ntot;
        _maxCS=-1;
        while (!fin.eof()) {
            getline( fin, line);
            if (line.find("99999")<line.size()) break;
            istringstream sline(line);
            for (int i=0; i<3; i++) {
                sline >> E >> exponent;
                E *= pow(10, exponent-3); // keV
                sline >> CS >> exponent;
                CS *= pow(10, exponent);
                //cout << "E=" << E << "   CS="<<CS<<endl;
                _CS->SetPoint(n++, E, CS);
                if (CS>_maxCS) _maxCS=CS;
                if (n>=ntot) break;
            }
        }
}


double
MaxCamENDF::generateCosAngleCMS(double energy) {
    // for given neutron energy, generate a cosine
    // of scattering angle in CMS

    double cosAngle, tmprnd, tmpcs;
    //int ncount=0;
    while (1) {
        tmprnd = gRandom->Rndm();
        cosAngle = (gRandom->Rndm()-0.5)*2;
        tmpcs = _DCS->Interpolate(energy+1e-5, cosAngle);
        assert(tmpcs>=0);// cout << "Negative DCS"<<endl;;
        if ( tmpcs > tmprnd * _maxCS) break;
        //if (ncount++>10000) return -2;
    }
    return cosAngle;
}



void
MaxCamENDF::fillFissionTable(const char* fileName) {
    // Load ENDF table into memory. 
    
    _DCS = new TGraph2D();
    
    const char* fileDir=gSystem->Getenv("MCTABLES");
    if (!fileDir) {
        cout << GetName() <<": Environment variable MCTABLES not defined." << endl; 
        assert(0);
    }    
    TString fullName(fileDir);
    fullName += "/";
    fullName += fileName;
        
    ifstream fin( fullName );

    string line;
    int NW;
    int  exponent;
    for (int i=0; i<9; i++) { 
        getline( fin, line); //cout << line << endl; 
    }
    istringstream sline(line);
    sline >> NW;
    assert (NW);
    //cout << "   NW="<<NW<<endl;
    vector<double> x,y;
    istringstream *pline=0;
    double tmp;
    _maxCS=-1;
    while (!fin.eof()) {
        
        for (int i=0; i<NW; i+=2) {                        
            if (i%6==0) { if (pline) delete pline; getline(fin, line); pline = new istringstream(line);  }
            (*pline) >> tmp >> exponent;
            tmp *= pow(10, exponent-3); // keV
            if (x.size() && tmp<x[x.size()-1]) { i=NW; break; }
            x.push_back( tmp );
            (*pline) >> tmp >> exponent;
            tmp *= pow(10, exponent);
            y.push_back( tmp );
            if (_maxCS<tmp) _maxCS=tmp;
            //cout<<tmp<<endl;
        }
    }
    _CS = new TGraph( x.size(), &x[0], &y[0]);
}


double
MaxCamENDF::generateEnergy(double emin, double emax) {
    // generate energy based on the cross section curve
    // dsigma/dE.
    //

    assert(_CS);
    double e=0;
    while (1) {
        e=gRandom->Rndm()*(emax-emin)+emin;
        if (gRandom->Rndm()<_CS->Eval(e)/_maxCS) break;
    }
    return e;
}


bool
MaxCamENDF::acceptEnergy(double Eneutron) {
    // scale energy relative CS; maximum CS has acceptance 100%
    if (getCrossSection()->Eval(Eneutron)/maxCS()<gRandom->Rndm()) return false;
    return true;
}

double
MaxCamENDF::generateRockNeutronEnergy(TString opt) {
    // this is a rough approximation to Mei,Hime simulation of
    // neutron spectrum in generic rock
    // Units are  1e-6/cm2/s/keV
    
    opt.ToLower();
    static TF1 *fRockNeutrons=0;
    if (!fRockNeutrons) {
        if (opt=="gransasso") 
            fRockNeutrons = new TF1("fRockNeutronsGranSasso","0.06*exp(-x/2)+0.07*exp(-0.5*(x-2.3)**2/0.5**2)", 0,8);
        else
            fRockNeutrons = new TF1("fRockNeutrons",
                                    "0.105*exp(-x/1.5)+0.15*exp(-0.5*(x-0.8)**2/0.2**2)+0.1*exp(-0.5*(x-2.5)**2/0.4**2)",0,8);
    }
        
    return fRockNeutrons->GetRandom(0,8)*1e3; //keV
}


double
MaxCamENDF::generateCosmicNeutronEnergy() {
  // spectrum based on Nakamura et al,J.Nucl.Sci.Tech,42,p843-853 (2005) 
  // Values scanned by Kazu Terao

  static TH1D *hNaka=0;

  if (!hNaka) {
    double y[]={ 0.000201645, 0.000100147, 7.056e-05, 6.00765e-05, 6.6612e-05, 8.1643e-05, 0.00010086, 0.000151868, 0.000215576,
		 0.000272886, 0.000527356, 0.000737288, 0.000843325, 0.000953602, 0.00100452, 0.000923981, 0.000714119, 0.000442771,
		 0.000536077, 0.000627261, 0.000898657, 0.000907154, 0.000928373, 0.00106408, 0.00107258, 0.001062, 0.00105566,
		 0.0010684, 0.000767358, 0.000661377, 0.000614756, 0.000451518, 0.000364606, 0.000222582, 0.0001272, 0.000127208};
    double x[]={ 2.01624e-07, 7.32594e-07, 2.93609e-06, 1.36784e-05, 0.000237753, 0.00188933, 0.0091776, 0.0511658, 0.154604,
		 0.354705, 0.690443, 1.1076, 1.63986, 2.42986, 3.60585, 5.36473, 8.17613, 11.8422, 15.3887, 19.4166, 24.389,
		 30.2095, 38.194, 49.0356, 61.7644, 78.1854, 100.549, 132.201, 173.345, 243.019, 319.627, 376.126, 452.811,
		 638.313, 887.314, 1014.29};
  
    int nbounds=sizeof(x)/sizeof(double)+1;
    double *xbounds=new double[ nbounds ];
    xbounds[0]=0;
    for (int i=0; i<nbounds-1; i++) x[i]*=1e3; // keV

    for (int i=1; i<=nbounds; i++) {
      xbounds[i]=2*x[i-1]-xbounds[i-1]; 
    }

    hNaka=new TH1D("hNaka","", nbounds-1, xbounds);
    for (int i=1; i<=nbounds; i++) hNaka->SetBinContent(i, y[i-1]);
  }

  return hNaka->GetRandom();

}
