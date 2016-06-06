#include "MaxCamMC.hh"
//#include "MaxCam.hh"
#include "MaxCamTwoBodyKinematics.hh"
#include "MaxCamImageTools.hh"
#include "MaxCamChannel.hh"
#include "TMath.h"
#include "TRotation.h"
#include "TChain.h"
#include "TF1.h"

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

#include "MaxCamUnits.hh"
using MaxCamUnits::mm;
using MaxCamUnits::cm;

ClassImp(MaxCamMC)

//____________________
//
//  Class for MC simulation of scintillation of tracks in drift
//  chamber. Wires are defined along x-direction.
//
    MaxCamMC::MaxCamMC(TString fname) : MaxCamSRIM(), MaxCamDataset(fname, "recreate") {
        // Constructor uses name of file with MC and a
        // pointer to the random number generator.
  
    _rnd=new TRandom2;
    setPhotonsPerkeV( 20 ); 
    _minEnergy=0.2; 
    _stepSize=0.0001;

    
   _rotation = new TRotation;
   setProjectile(-14.1e3, 0,0, 1e6); // default is 14.1MeV neutron
   setRecoil(100,0,0,19e6); // default is 100keV fluorine
   _recoilCoord=0;

   setElecScintillation(1.0); // set all dEdx_elec creates light
   setNuclScintillation(1.0); // set all dEdx_nucl creates light

   // set signal width
   setAvalancheWidth(0.1);
   setDiffusionConstTerm(0.1046);
   setDiffusionDzTerm(0.01325);

   setStraggling(true);

   setTrackImage(); //  2x2cm^2    
   setWireImage();  //  2x2cm^2 
   setCCDImage();   //  768 x 512, 8 x 8 bins

   setNoiseADC(25); // noise at 1sec exposure

   /////for (int i=0; i<20; i++) addWire( -50+i*5 ); // wires in mm


   _elf= new MaxCamElectricField();   
   
   _imageTree->Branch("calibration", &_calibration,
                      "photons_keV/D:pixel_mm:elecScint:nuclScint:noiseADC:avalWidth:diff0:diffDz");
   _imageTree->Branch("recoil","TLorentzVector",&_recoil, 32000, 0);
   _imageTree->Branch("projectile","TLorentzVector",&_projectile, 32000, 0);

    }


MaxCamMC::MaxCamMC(const MaxCamMC &other) :
    MaxCamSRIM(other), MaxCamDataset(other)  {
    _rnd=other._rnd;
    _projectile=other._projectile;
    _recoil=other._recoil;
    cout << GetName() << "Copy Constructor not done" << endl;
}

MaxCamMC::~MaxCamMC() {
  delete _rnd;
  delete _projectile;
  delete _recoil;
}


void 
MaxCamMC::setPressure(double p) {
    // Set gas pressure
    
    _pressure->currentValue=p;
    _pressure->setValue=p;
    setStopping(p);
}



void 
MaxCamMC::setProjectile(double px, double py, double pz, double M) {
    // Set projectile momentum vector, mass.

    TVector3 p1(px,py,pz);
    TVector3 p2=p1.Orthogonal();
    _projectile = new TLorentzVector(p1, sqrt(px*px+py*py+pz*pz+M*M));
    _rotation->SetXAxis(p1,p2);
}

void
MaxCamMC::setIsotropicProjectile(double P, double M) {
    double phi=gRandom->Rndm()*TMath::TwoPi();
    double costheta=(gRandom->Rndm()-0.5)*2;
    double sintheta=sqrt(1-costheta*costheta);
    double px=P*sintheta*cos(phi);
    double py=P*sintheta*sin(phi);
    double pz=P*costheta;
    setProjectile(px, py, pz, M);
}



void
MaxCamMC::setRandomRecoilCoord( double x1, double x2, double y1, double y2, double z1, double z2) {
    // Set random interaction point withing given volume.

    float x = x1 + (x2-x1)*_rnd->Rndm();
    float y = y1 + (y2-y1)*_rnd->Rndm();        
    float z = z1 + (z2-z1)*_rnd->Rndm();
    setRecoilCoord(x,y,z);
}


double
MaxCamMC::calcCosRecoil(double p) {
    // Elastic, non-relativistic two body scattering.
    // Computes cosine of recoil angle for given recoil 
    // energy and mass:
    // 
    //  cos(Theta) = 1/2 * sqrt(E/E0) * (1+M/m) / sqrt(M/m)
    //

    return MaxCamTwoBodyKinematics::calcCosRecoilFromRecoilEnergy( p, _projectile->Vect().Mag(), _recoil->M(), _projectile->M());
}


int
MaxCamMC::makeRecoil(double p, double phi) {
        // Set recoil parameters (NR) using momentum and
        // azimuthal (phi) angle as input. The recoil angle
        // is computed from the energy so this function allows
        // study of tracks with certain energy.

        double cosRecoil = calcCosRecoil(p);
        if (cosRecoil<-1 || cosRecoil>1) { cout << "cosRecoil="<<cosRecoil<<endl; return -1; }
        double sinRecoil = sqrt(1-cosRecoil*cosRecoil);
            
        double px = p*cosRecoil;
        double py = p*sinRecoil*cos(phi);
        double pz = p*sinRecoil*sin(phi);
        double E = sqrt(px*px+py*py+pz*pz+_recoil->M()*_recoil->M());       

        // set recoil in projectile frame
        TVector3 pRecoil(px,py,pz);
        // rotate to lab frame
        pRecoil = (*_rotation) * pRecoil;        
                
        if (_recoil) delete _recoil;
        _recoil = new TLorentzVector(pRecoil,E);
        return 0;
} 

int
MaxCamMC::setRecoilEnergy(double p) {
        // Set recoil energy
        //

        double phi=TMath::TwoPi()*_rnd->Rndm(); 
        return makeRecoil( p, phi );
}


int
MaxCamMC::setRecoilEnergyAngle(double p, double phi) {
        // Set recoil energy and azimuthal angle (in projectile frame)
        //

        return makeRecoil( p, phi );
}

void
MaxCamMC::setRecoilMass(double M) {
    setRecoil(0,0,0,M);
}


double
MaxCamMC::applyDiffusion(double x) {
    // Apply measured diffusion to the coordinate along wire.
    // tmp solution: assume the end of drift region is z=0.

    double dz=getRecoilCoord()->z()*0.1; // in cm
    if (dz<0) dz=0;
    // sigma in mm, z coord in cm
    double sigma=sqrt( getDiffusionConstTerm() + getDiffusionDzTerm()*dz ); 
    return _rnd->Gaus(x, sigma);
}


void
MaxCamMC::applyBiasADC(TH2F* h) {
    // tmp solution is hardcoded width width of 10 and no bias
    
    int nb=(h->GetNbinsX()+1)*(h->GetNbinsY()+1);
    for (int i=1; i<=nb; i++) h->SetBinContent( i, _rnd->Gaus(h->GetBinContent(i), getNoiseADC() ) );
}

void
MaxCamMC::applyRadialEffect(TH2F* h)
{
   int xbins=h->GetNbinsX()+2;
   int ybins=h->GetNbinsY()+2;

   double xlow = h->GetXaxis()->GetXmin();
   double xhigh = h->GetXaxis()->GetXmax();
   double ylow = h->GetYaxis()->GetXmin();
   double yhigh = h->GetYaxis()->GetXmax();
   double xcenter = (xlow+xhigh)/2;
   double ycenter = (ylow+yhigh)/2;

   //we want to scale all things by the maximum radius, 
   //which is the center to any corner.
   double maxradius = sqrt(pow((xlow-xcenter),2)+pow(ylow-ycenter,2));

   TF1* radialfn = new TF1("radialfn","[0]+[1]*x",0,1);
   radialfn->SetParameters(1,0);
   
   for(int i=1;i<xbins;i++)
   {
      for(int j=1;j<ybins;j++)
      {
	 double x = h->GetXaxis()->GetBinCenter(i);
	 double y = h->GetYaxis()->GetBinCenter(j);
	 double radius = sqrt(pow(x-xcenter,2)+pow(y-ycenter,2));
	 radius = radius/maxradius;
	 double scalefactor = radialfn->Eval(radius);
	 h->SetBinContent(i,j,scalefactor*h->GetBinContent(i,j));
      }
   }

   
}

double
MaxCamMC::applyAvalancheWidth(double r) {
    // Apply measured width of the coordinate perpendicular to wires.
    // tmp solution: hardcode width to 100um

    return _rnd->Gaus(r, getAvalancheWidth() );
}


bool
MaxCamMC::propagateRecoil() {
        // Propagate a recoil product by a step-size.

        double dR=getStepSize();

        assert( getRecoil() );
        assert( getStoppingVsEnergy() );
        assert( getTrackImage() );
        getTrackImage()->Clear();


        // current energy, step size
        //double dR = 0.0001;//getRangeStep()->Eval(E);
        double E  = getRecoil()->Vect().Mag();
        double scaleR = 1;//sqrt( dR / getRangeVsEnergy()->Eval( E ) );
        double lstraggle = _rnd->Gaus(  1, getLStraggleVsEnergy()->Eval(E)/getRangeVsEnergy()->Eval(E) * scaleR );
        double tstraggle = _rnd->Gaus(  0, getTStraggleVsEnergy()->Eval(E)/getRangeVsEnergy()->Eval(E) * scaleR );
        double nstraggle = _rnd->Gaus(  0, getTStraggleVsEnergy()->Eval(E)/getRangeVsEnergy()->Eval(E) * scaleR );
        TVector3 longVec = getRecoil()->Vect().Unit();
        TVector3 tranVec = longVec.Orthogonal();
        TVector3 normVec = longVec.Cross(tranVec);
        TVector3 stepVec = longVec*lstraggle + tranVec*tstraggle + normVec*nstraggle;
        stepVec=stepVec.Unit();
        
        //getRecoil()->Vect().Print();
        
         // get stopping power
        double dE = getStoppingVsEnergy()->Eval( E ) * dR; 
        //cout << E<<"  "<<dE<<"  "<< stepVec.Z()<<"  " <<getLStraggleVsEnergy()->Eval(E) <<endl;
        double dE_elec = getStoppingVsEnergy(false)->Eval( getRecoil()->Vect().Mag() ) * dR; // elec dE/dx

        // create poisson n photons based on dE_elec*f_elec + dE_nucl*f_nucl
        // where f_elec (f_nucl) is the fraction of electronic (nuclear) stopping power 
        // that scintillates
        double dE_scint = (dE-dE_elec)*getNuclScintillation() + dE_elec*getElecScintillation();
        double npho = _rnd->PoissonD( getPhotonsPerkeV()*dE_scint );
        _totalPhotons += npho;
        // fill XY histogram
        getTrackImage()->Fill( getRecoilCoord()->x(), getRecoilCoord()->y(), npho );
       
        // drift electrons
        driftElectrons(npho);

        // move track
        //TVector3 stepVec=getRecoil()->Vect().Unit();
        *_recoilCoord += stepVec*dR;

        
        // decrease energy
        TVector3 recoilVec = getRecoil()->Vect() - stepVec*(dE);
        //TVector3 recoilVec = stepVec*(E-dE);
        getRecoil()->SetVectM( recoilVec, getRecoil()->M() );
        
        return E-dE>_minEnergy && E-dE>dE;
}

void
MaxCamMC::driftElectrons(double npho) {

// old code:
    double xdiff, ydiff;        
    xdiff = applyDiffusion(getRecoilCoord()->x());
    ydiff = applyDiffusion(getRecoilCoord()->y());
    if (getAnodeWireListX().size() || getAnodeWireListY().size()) {  // if using wires...
        bool isX = findAnodeWire( xdiff, ydiff );
        if (isX)  xdiff = applyAvalancheWidth( xdiff );
        else      ydiff = applyAvalancheWidth( ydiff );
    }
    getWireImage()->Fill( xdiff, ydiff, npho);
    return;
    
// electron-drift code: (NOT DONE)
    double xdrift = getRecoilCoord()->x();
    double ydrift = getRecoilCoord()->y();
    double zdrift = getRecoilCoord()->z();
    double x0,y0,z0;
    //double D  = getDiffusionDzTerm();
    double step = 0.1*mm;

    
    bool isClose=false;
    xdrift*=mm;
    ydrift*=mm;
    zdrift*=mm;
    while (!isClose) {
        x0=xdrift;
        y0=ydrift;
        z0=zdrift;
        _elf->driftStep(ydrift, zdrift, step); // convert to mm        
        // ask if close enough
        for (unsigned int i=0; i<getAnodeWireListY().size(); i++) {
            if ( fabs(getAnodeWireListY()[i]*mm-ydrift)<3*step &&
                 fabs(zdrift)<3*step ) { isClose=true; break; } // assumes anode at z=0!!!!!!!
        }
    }
    xdrift/=mm;
    ydrift/=mm;
    zdrift/=mm;
    
    findAnodeWire( xdrift, ydrift);
    xdrift = applyDiffusion(xdrift);
    ydrift = applyAvalancheWidth( ydrift );
    
    getWireImage()->Fill( xdrift, ydrift, npho);
}



 bool
 MaxCamMC::findAnodeWire(double &x, double &y) {
        int ix=findClosestWire(x, _anodeWireListX ); 
        int iy=findClosestWire(y, _anodeWireListY );
        assert (ix>-1 || iy>-1);
        bool isX=true;
        //cout << "anode wire for " << x <<","<<y<<"  -> " << ix << "  " << iy;;
        if (ix>-1 && iy<0)  isX=true; // only x-wires
        else if (ix<0  && iy>-1) isX=false; // only y-wires
        else isX = fabs(getAnodeWireListX()[ix]-x)<fabs(getAnodeWireListY()[iy]<y) ? true : false; // both x&y
        if (isX) {
            x=getAnodeWireListX()[ix]; 
        }
        else {
            y=getAnodeWireListY()[iy];
        }
        //cout << "   new: " << x << ", " << y << endl;
        return isX;
    }



void
MaxCamMC::event(bool resetImage) {

        if (resetImage) {
            getTrackImage()->Reset();
            getWireImage()->Reset();
            getCCDImage()->Reset();
            _totalPhotons=0;
        }

        //cout <<"Anode wires="<< getAnodeWireListX().size() << "   " << getAnodeWireListY().size() << endl;
        
        //   Propagate recoil
        //
        //TVector3 x0=*getRecoilCoord(); // initial position
        TLorentzVector recoil0( *_recoil );

        // straggling
        double E  = getRecoil()->Vect().Mag();
        double R  = getRangeVsEnergy()->Eval(E);
	double lstraggle,tstraggle,nstraggle;
	if(getStraggling())
	{

	   lstraggle = _rnd->Gaus(  R, getLStraggleVsEnergy()->Eval(E) );
	   tstraggle = _rnd->Gaus(  0, getTStraggleVsEnergy()->Eval(E) );
	   nstraggle = _rnd->Gaus(  0, getTStraggleVsEnergy()->Eval(E) );
	}
	else
	{
	   lstraggle = R;
	   tstraggle = 0;
	   nstraggle = 0;
	}
        TVector3 longVec = getRecoil()->Vect().Unit();
        TVector3 tranVec = longVec.Orthogonal();
        TVector3 normVec = longVec.Cross(tranVec);
        TVector3 stepVec = longVec*lstraggle + tranVec*tstraggle + normVec*nstraggle;
        TVector3 recoilVec  = stepVec.Unit()*(E*lstraggle/R);
        getRecoil()->SetVectM( recoilVec, getRecoil()->M() );
        
        
        while( propagateRecoil() );
        //TVector3 range=*getRecoilCoord()-x0;
        //cout << "Range="<< range.Mag()<<"  " << range.x()<<","<<range.y()<<","<<range.z()<<endl;

	if(getSpacerListX().size() || getSpacerListY().size()) 
	{applySpacers(getWireImage());}

	//apply radial effect
	applyRadialEffect(getWireImage());

        //   Make CCD Image
        // 
        MaxCamImageTools::resizeImage( getWireImage(), getCCDImage() ); // convert to CCD coordinates

        //  Add ADC noise
//         static int m1=0; if (m1++<10) cout <<"no ADC bias"<<endl;
        applyBiasADC( getCCDImage() );

        // save to ntuple
        *_recoil=recoil0;
        _img_histo=getCCDImage();
        _imageTree->Fill();
}



TGraph*
MaxCamMC::getRangeVsEnergyProject(TString opt) {
        // Compute projection of range to longitudinal (opt=long) or transverse (opt=tran)
        // axis of the initial neutron direction.
        // 
        
        int n=getRangeVsEnergy()->GetN();
        vector<double> RcosRecoil;
        double E,R,cosRec;
        bool isLong = opt.Contains("long") ? true : false;
        for (int i=0; i<n; i++) {
                getRangeVsEnergy()->GetPoint(i,E,R);
                cosRec = calcCosRecoil(E);
                if (fabs(cosRec)>1) cosRec=0;
                R *= isLong ? cosRec : sqrt(1-cosRec*cosRec);
                RcosRecoil.push_back( R );
        }
        return new TGraph(n, getRangeVsEnergy()->GetX(), &RcosRecoil[0]);
}



TGraph*
MaxCamMC::getEnergyVsRangeProject(TString opt) {
    // Compute total energy from projection of range to longitudinal (opt=long)
    // or transverse (opt=tran) axis of the initial neutron direction.
    // 
    
    int n=getEnergyVsRange()->GetN();
    vector<double> RcosRecoil;
    double E,R,cosRec;
    bool isLong = opt.Contains("long") ? true : false;
    for (int i=0; i<n; i++) {
        getEnergyVsRange()->GetPoint(i,R,E);
        cosRec = calcCosRecoil(E);
        if (fabs(cosRec)>1) break;
        R *= isLong ? cosRec : sqrt(1-cosRec*cosRec);
        RcosRecoil.push_back( R );
    }
    return new TGraph(RcosRecoil.size(), &RcosRecoil[0], getEnergyVsRange()->GetY());
}




void
MaxCamMC::applySpacers(TH2* image) {
   bool isX=false;
   if(getSpacerListX().size()) {isX=true;}
   if(getSpacerListY().size()) {isX=false;}

   vector<double> spacerlist;
   vector<double> spacerlow;
   vector<double> spacerup;
   int nbins;
   int nopp;
   TAxis* axis;
   TAxis* opaxis;

   if(isX)
   {
      spacerlist = getSpacerListX();
      nbins = image->GetNbinsX();
      nopp = image->GetNbinsY();
      axis = image->GetXaxis();
      opaxis = image->GetYaxis();
   }

   if(!isX)
   {
      spacerlist = getSpacerListY();
      nbins = image->GetNbinsY();
      nopp = image->GetNbinsX();
      axis = image->GetYaxis();
      opaxis = image->GetXaxis();
   }


   for(unsigned int i = 0; i<spacerlist.size(); i++)
   {
      spacerlow.push_back(spacerlist[i]-getSpacerDiameter()/2);
      spacerup.push_back(spacerlist[i]+getSpacerDiameter()/2);

      for(int j = 0; j <= nbins; j++)
      {
	 double binlow = axis->GetBinLowEdge(j);
	 double binup = axis->GetBinUpEdge(j);
	 double redfac=1;
	 
	 // fitted spacers
	 // fitted fxn is 1-0.85*exp(-abs(x)^1.5/18) where x is in mm
	 double spacermid = (spacerup[i] + spacerlow[i])/2;

	 // If you know actual spacer diameter, modify scaling factor in line below
	 double sFactor = .3784; // empirical scaling factor
	 double spacerdiam = sFactor*getSpacerDiameter();

	 if(spacerlow[i] > binlow && spacerup[i] < binup)
	   {
	     redfac = 1-0.85*0.3018*spacerdiam/(binup-binlow);
	   }

	 // If pixel overlaps only part of spacer, integrate over overlapping region
	 if ((spacerlow[i] < binlow && spacerup[i] > binup) ||
	     (spacerlow[i] > binlow && spacerlow[i] < binup && spacerup[i] > binup) ||
	     (spacerup[i] > binlow && spacerup[i] < binup && spacerlow[i] < binlow))
	   {
	     double tempSum=0;
	     int counter = 0;
	     // Stepsize is currently set relative to bin size.
	     // Used counter so that absolute stepsize can be chosen.
	     double stepsize = (binup-binlow)/10;
	     for(double mm=binlow; mm<=binup; mm=mm+stepsize)
	       {
		 tempSum = tempSum + exp(-pow(fabs((mm-spacermid)/spacerdiam),1.5)/18);
		 counter = counter+1;
	       }
	     double Sum = tempSum/counter;
	     redfac = 1 - 0.85*Sum;
	   }


	 // The cosine and square spacers aren't as accurate
	 // cosine spacers
// 	 double pi=3.141592654;
// 	 double spacermid = (spacerup[i] + spacerlow[i])/2;
// 	 double spacerdiam = getSpacerDiameter();
// 	 if(spacerlow[i] > binlow && spacerup[i] < binup)
// 	    redfac = 1-spacerdiam/(2*(binup-binlow));
// 	 if(spacerlow[i] < binlow && spacerup[i] > binup)
// 	    redfac = 1-((sin(2*pi*(binup-spacermid)/spacerdiam)
// 			 - sin(2*pi*(binlow-spacermid)/spacerdiam))/2/pi
// 			+ (binup-binlow)/spacerdiam);
// 	 if(spacerlow[i] > binlow && spacerlow[i] < binup && spacerup[i] > binup)
// 	    redfac = 1-(sin(2*pi*(binup-spacermid)/spacerdiam)/2/pi
// 			+ (binup-spacerlow[i])/spacerdiam);
// 	 if(spacerup[i] > binlow && spacerup[i] < binup && spacerlow[i] < binlow)
// 	    redfac = 1-(-sin(2*pi*(binlow-spacermid)/spacerdiam)/2/pi
// 			+ (spacerup[i]-binlow)/spacerdiam);

	 // square spacers
// 	 if(spacerlow[i] > binlow && spacerup[i] < binup)
// 	    redfac = 1-getSpacerDiameter()/(binup-binlow);
// 	 if(spacerlow[i] < binlow && spacerup[i] > binup)
// 	    redfac = 0;
// 	 if(spacerlow[i] > binlow && spacerlow[i] < binup && spacerup[i] > binup)
// 	    redfac = (spacerlow[i]-binlow)/(binup-binlow);
// 	 if(spacerup[i] > binlow && spacerup[i] < binup && spacerlow[i] < binlow)
// 	    redfac = (binup-spacerup[i])/(binup-binlow);

	 if(redfac < 1)
	 {
	    for(int k = 0;  k <= nopp; k++)
	    {
	       if(isX)
	       {
		  double contents = image->GetBinContent(j,k)*redfac;
		  image->SetBinContent(j,k,contents);
	       }
	       if(!isX)
	       {
		  double contents = image->GetBinContent(k,j)*redfac;
		  image->SetBinContent(k,j,contents);
	       }
	    } 
	 }
      }
   }
}
