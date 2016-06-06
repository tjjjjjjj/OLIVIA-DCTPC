#ifndef MAXCAM_ELECTRICFIELD_HH
#define MAXCAM_ELECTRICFIELD_HH

#include <assert.h>

#include "TROOT.h"
#include "TString.h"
#include "MaxCamUnits.hh"
#include "TVectorD.h"

#include "MaxCamUnits.hh"

#include <vector>
using std::vector;

class TGraph;




class MaxCamElectricField {

public:

  // Ctors
    MaxCamElectricField();

    //MaxCamElectricField(const MaxCamElectricField &other) {}

    virtual ~MaxCamElectricField() {}

    TString GetName() { return "MaxCamElectricField"; }
    

    double V(double *var);
    double V(double x, double y);
    double Ex(double x, double y);
    double Ey(double x, double y);
    TGraph* fieldLine(double x, double y, double pressure=200, double step=2e-6, int charge=-1);
    void driftStep(double &x, double &y, double step);

    void resetCharges();


    TGraph* getFieldLine() { return _fieldLine; }
    double  getElectricFieldSum() { return _electricFieldSum; }
    double  getAlphaSum() { return _alphaSum; }
    double  alpha(double E, double P); // compute ionization for given field, pressure
    TGraph* getAlphaTable() { return _ap_vs_ep; }
    TGraph* getAttachmentTable() { return _at_vs_ep; }
    double  attachment(double E, double P); // compute attachemnt for given field, pressure
    
    int getStatFlag() { return _statFlag; }

    struct Wire {
        double V;
        double x, y;
        double radius;
    };
    vector<Wire> wireList;
    
    void setDriftPlane(double dr) { _dr=dr; resetCharges(); }
    void setFieldPlane(double h2, double voltage) {
        _h2=h2; if (_h2<0) _h2=-_h2; _Vp=voltage; resetCharges();
    }
    void setWirePlane(double y, double voltage, double radius, double pitch, int n, double x0) {
        for (int i=0; i<n; i++) {
            Wire wire;
            wire.V=voltage;
            wire.x=x0+i*pitch;;
            wire.y=y;
            wire.radius=radius;
            wireList.push_back(wire);
        }
        resetCharges();
    }

    void Print();


    
private:

    double calcL(double x, double y, double xn, double yn, double r);
    
    TVectorD* calcCharge();

    TGraph *_fieldLine;
    double _electricFieldSum;
    double _alphaSum;
    int _statFlag;

    TGraph *_ap_vs_ep; // alpha/P vs. E/P
    TGraph *_at_vs_ep; // attachment vs. E/P
    
    
// ------ voltages:
    
    double _Vc;
    double _Va;
    double _Vp;
    
// ------ geometry:
    
    double _dr; //  26*mm;
    double _h1; //  3*mm;
    double _h2; //  3*mm;

    
// ------- derived:
    
    double _d; //  h1-h2+dr;
    double _D; //  h1+h2+dr;
    double _a; //  pi/D;
    double _VR; //  Vc/Va;
    double _VP; //  Vp/Va;
    TVectorD *_Qa;  //
    
    ClassDef(MaxCamElectricField,0)

        };

#endif

