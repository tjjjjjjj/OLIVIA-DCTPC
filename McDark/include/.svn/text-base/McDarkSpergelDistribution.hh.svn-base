#ifndef McDarkSpergelDistribution_H
#define McDarkSpergelDistribution_H 1

#include "TF2.h"

class McDarkSpergelDistribution {

public:

    McDarkSpergelDistribution();

    void GetRandomValues(Double_t[]);

    void setDarkMatterMass(double m) { differentialRateSpergel->SetParameter(0,m); }

    void setTargetAtomicNumber(int A) { differentialRateSpergel->SetParameter(1,A); }

    void setDarkMatterAvgVelocity(double v) { differentialRateSpergel->SetParameter(2,v); }

    void setEarthVelocityAroundSun(double v) { differentialRateSpergel->SetParameter(3,v); }
    
    void setEarthVelocityVariation(double dv) { differentialRateSpergel->SetParameter(4,dv); }

    void seTimeOfYear(double y) {  differentialRateSpergel->SetParameter(5,y); }

    void setDarkMatterCrossSection(double cs) {  differentialRateSpergel->SetParameter(6,cs); }

    void setDarkMatterDensity(double d) {  differentialRateSpergel->SetParameter(7,d); }
   
    
    TF2* getFunction() { return differentialRateSpergel; }
    

    
private:
    
    TF2* differentialRateSpergel;

};
    
#endif
