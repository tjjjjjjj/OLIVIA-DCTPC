

#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include "TF1.h"

double minu=1;
double sec=1./60;
double hour=60;
double day=hour*24;
double year=day*365;


double T12_U238  = 5e9*year;
double T12_Th234 = 0;
double T12_Pa234 = 0;
double T12_U234  = 0.2e6*year;
double T12_Th230 = 75e3*year;
double T12_Ra226 = 1.6e3*year;

double T12_Rn222 = 3.82*day;
double T12_Po218 = 3.1*minu;
double T12_Pb214 = 27*minu;
double T12_Bi214 = 20*minu;
double T12_Po214 = 160e-6*sec;


TMatrixD *Lambda=0;
TMatrixD *D=0;
TMatrixD *V=0;
TMatrixD *VInv=0;
TMatrixD *N0=0;

bool doCheck=true;



void initDecay() {

  int n=5;
  Lambda = new TMatrixD( n, n );
  (*Lambda)(0,0) = - log(2)/T12_Rn222;
  (*Lambda)(1,0) =   log(2)/T12_Rn222;   (*Lambda)(1,1) = - log(2)/T12_Po218;
  (*Lambda)(2,1) =   log(2)/T12_Po218;   (*Lambda)(2,2) = - log(2)/T12_Pb214;
  (*Lambda)(3,2) =   log(2)/T12_Pb214;   (*Lambda)(3,3) = - log(2)/T12_Bi214;
  (*Lambda)(4,3) =   log(2)/T12_Bi214;   (*Lambda)(4,4) = - log(2)/T12_Po214;

  N0=new TMatrixD(n,1);
  (*N0)(0,0)=1e6;

  TMatrixDEigen eig( *Lambda );
  D = new TMatrixD( eig.GetEigenValues() );
  //D->Print();

  V=new TMatrixD ( eig.GetEigenVectors() );
  //V->Print();

  VInv = new TMatrixD( *V );
  VInv->Invert();
  //VInv->Print();

  if (doCheck) {
    Lambda->Print();
    TMatrixD check1(n,n), check2(n,n);
    check1.Mult(*V, *D);
    check2.Mult(check1, *VInv);
    check2.Print();
  }

}



double 
Bateman(double *var, double *par) {
  
  double t=var[0];
  double daughterType=par[0];

  int n=Lambda->GetNrows();
  TMatrixD Dt(n,n);
  for (int i=0; i<n; i++) Dt(i,i)=exp( (*D)(i,i)*t);
  //D->Print();
  //Dt.Print();

  TMatrixD tmp1(n,n), tmp2(n,n), tmp3(n,1);
  tmp1.Mult( *V, Dt); //tmp1.Print();
  tmp2.Mult( tmp1, *VInv); //tmp2.Print();
  tmp3.Mult( tmp2, *N0 ); //tmp3.Print();
  
  int type=int( daughterType+0.5 );

  return tmp3(type,0);

}
  


void init() {

  initDecay();
  TF1 *fDecay= new TF1("fDecay",Bateman, 0, 10*60, 1);
  fDecay->SetParameter(0,  2);  // 0=Rn222, 1=Po218, 2=Pb214, 3=Bi214, 4=Po214
  fDecay->Draw();


  (*Lambda)(0,0) = - log(2)/T12_Rn222;
  (*Lambda)(1,0) =   log(2)/T12_Rn222;   (*Lambda)(1,1) = - log(2)/T12_Po218;
  (*Lambda)(2,1) =   log(2)/T12_Po218;   (*Lambda)(2,2) = - log(2)/T12_Pb214;
  (*Lambda)(3,2) =   log(2)/T12_Pb214;   (*Lambda)(3,3) = - log(2)/T12_Bi214;
  (*Lambda)(4,3) =   log(2)/T12_Bi214;   (*Lambda)(4,4) = - log(2)/T12_Po214;


}

