#include "assert.h"
#include "RobinHoodTriangleMaker.hh"
#include "TVector3.h"
#include "TRotation.h"
#include "TMath.h"

#include <iostream>
using std::cout;
using std::endl;
using std::ios;

#include <fstream>
using std::ofstream;

ClassImp(RobinHoodTriangleMaker)


void 
RobinHoodTriangleMaker::cylinder( TVector3 *center, TVector3 *norm, 
				  double ID, double OD, double thick, int nTriangles,
				  TString opt) {

  assert( ID>0 && OD>0 && thick>0 && nTriangles>0 && center && norm );

  if (!opt.Contains("add")) _list.Clear();



  // Distribute triangles over surfaces
  double area_flat  = TMath::Pi() * 0.25 * ( OD*OD - ID*ID );
  double area_outer = TMath::Pi() * OD * thick;
  double area_inner = TMath::Pi() * ID * thick;
  double area = area_flat * 2 + area_outer + area_inner;

  int nt_flat  = nTriangles * area_flat  / area;
  int nt_outer = nTriangles * area_outer / area;
  int nt_inner = nTriangles * area_inner / area;


  cout << nt_flat<< "  "<<  nt_outer<<"  " <<  nt_inner<<endl;


  //
  //
  // Split outer surface:
  //
  //
  int nt_outer_edge = sqrt(nt_outer * TMath::Pi()*OD / thick);
  int nt_outer_thick = nt_outer / nt_outer_edge;
  if (nt_outer_edge<1) nt_outer_edge=1;
  if (nt_outer_thick<1) nt_outer_thick=1;

  double dphi = TMath::TwoPi() / nt_outer_edge;
  
  TVector3 tCenter = (*center)+thick*0.5*(*norm).Unit();
  TVector3 dthick = (*norm).Unit();
  dthick *= thick / nt_outer_thick;
  TVector3 nphi = (*norm).Orthogonal().Unit()*0.5*OD;

  TRotation *rotation = new TRotation;
  rotation->SetZAxis( (*norm), nphi );

  if (!opt.Contains("notout")) for (int jt=0; jt<nt_outer_thick; jt++) {

    TVector3 t1 = tCenter - dthick*jt;
    TVector3 t2 = t1 - dthick;

    for (int it=0; it<nt_outer_edge; it++) {

      TVector3* a=new TVector3( t1+rotation->RotateZ(dphi*0.5)*nphi );
      TVector3* b=new TVector3( t2+rotation->RotateZ(dphi*0.5)*nphi );

      TVector3* c=new TVector3( t1+rotation->RotateZ(dphi*0.5)*nphi );
      TVector3* d=new TVector3( t2+rotation->RotateZ(dphi*0.5)*nphi );
      
      rotation->RotateZ(-dphi);

      RobinHoodTriangle *T1=new RobinHoodTriangle( a, b, c, 0);
      RobinHoodTriangle *T2=new RobinHoodTriangle( b, c, d, 0);

      _list.Add(T1);
      _list.Add(T2);

    }
  }

  cout << "Total triangles from outer surface = " << _list.GetSize() << endl;

  //
  //
  //     Top and bottom surface
  //
  //
  TVector3 dr=nphi.Unit();
  int nt_flat_phi   = sqrt( nt_flat*(OD+ID)*0.5*TMath::Pi()/(OD-ID) );
  int nt_flat_r = nt_flat/nt_flat_phi; 
  if (nt_flat_r<1) nt_flat_r=1;
  if (nt_flat_phi<1) nt_flat_phi=1;
  cout << nt_flat_r << "   " << nt_flat_phi << endl;
  dr *= 0.5*(OD-ID)/nt_flat_r;
  dphi = TMath::TwoPi() / nt_flat_phi; 
  TVector3 bCenter = (*center) - thick*0.5*(*norm).Unit();
  if (!opt.Contains("nottop")) for (int jt=0; jt<nt_flat_r; jt++) {

    for (int it=0; it<nt_flat_phi; it++) {
    
      TVector3 da = rotation->RotateZ(dphi*0.5)*nphi;
      TVector3 db = rotation->RotateZ(dphi*0.5)*(nphi-dr);
      TVector3 dc = rotation->RotateZ(dphi*0.5)*nphi;
      TVector3 dd = rotation->RotateZ(dphi*0.5)*(nphi-dr);

      TVector3* a=new TVector3( tCenter+da );
      TVector3* b=new TVector3( tCenter+db );
      TVector3* c=new TVector3( tCenter+dc );
      TVector3* d=new TVector3( tCenter+dd );
      
      TVector3* A=new TVector3( bCenter+da );
      TVector3* B=new TVector3( bCenter+db );
      TVector3* C=new TVector3( bCenter+dc );
      TVector3* D=new TVector3( bCenter+dd );
      
      rotation->RotateZ(-dphi);
      
      RobinHoodTriangle *T1=new RobinHoodTriangle( a, b, c, 0);
      RobinHoodTriangle *T2=new RobinHoodTriangle( b, c, d, 0);
      RobinHoodTriangle *T3=new RobinHoodTriangle( A, B, C, 0);
      RobinHoodTriangle *T4=new RobinHoodTriangle( B, C, D, 0);

      _list.Add(T1);
      _list.Add(T2);
      _list.Add(T3);
      _list.Add(T4);
    }
    nphi -= dr;
  }

  cout << "Total triangles from top+bottom+outer surface = " << _list.GetSize() << endl;



  //
  //
  // Split inner surface:
  //
  //
  int nt_inner_edge = sqrt(nt_inner * TMath::Pi()*ID / thick);
  int nt_inner_thick = nt_inner / nt_inner_edge;
  if (nt_inner_edge<1) nt_inner_edge=1;
  if (nt_inner_thick<1) nt_inner_thick=1;

  dphi = TMath::TwoPi() / nt_inner_edge;
  dthick = (*norm).Unit();
  dthick *= thick / nt_inner_thick;
  nphi = (*norm).Orthogonal().Unit()*0.5*ID;


  if (!opt.Contains("notin")) for (int jt=0; jt<nt_inner_thick; jt++) {
    TVector3 t1 = tCenter - dthick*jt;
    TVector3 t2 = t1 - dthick;

    for (int it=0; it<nt_outer_edge; it++) {


      TVector3* a=new TVector3( t1+rotation->RotateZ(dphi*0.5)*nphi );
      TVector3* b=new TVector3( t2+rotation->RotateZ(dphi*0.5)*nphi );

      TVector3* c=new TVector3( t1+rotation->RotateZ(dphi*0.5)*nphi );
      TVector3* d=new TVector3( t2+rotation->RotateZ(dphi*0.5)*nphi );
      
      rotation->RotateZ(-dphi);

      RobinHoodTriangle *T1=new RobinHoodTriangle( a, b, c, 0);
      RobinHoodTriangle *T2=new RobinHoodTriangle( b, c, d, 0);

      _list.Add(T1);
      _list.Add(T2);

    }
  }
  cout << "Total triangles = " << _list.GetSize() << endl;

  


}
  

void 
RobinHoodTriangleMaker::disk( TVector3 *center, TVector3 *norm, 
			      double Diameter, int nTriangles,
			      TString opt) {

  assert( Diameter>0 && nTriangles>0 && center && norm );

  if (!opt.Contains("add")) _list.Clear();


  TVector3 nphi = (*norm).Orthogonal().Unit()*0.5*Diameter;
  TRotation *rotation = new TRotation;
  rotation->SetZAxis( (*norm), nphi );


  TVector3 dr=nphi.Unit();
  int nt_flat_phi   = sqrt( nTriangles );
  int nt_flat_r = nt_flat_phi; 
  if (nt_flat_r<1) nt_flat_r=1;
  if (nt_flat_phi<1) nt_flat_phi=1;
  cout << nt_flat_r << "   " << nt_flat_phi << endl;
  dr *= 0.5*Diameter/nt_flat_r;
  double dphi = TMath::TwoPi() / nt_flat_phi; 

  for (int jt=0; jt<nt_flat_r; jt++) {

    for (int it=0; it<nt_flat_phi; it++) {
    
      TVector3 da = rotation->RotateZ(dphi*0.5)*nphi;
      TVector3 db = rotation->RotateZ(dphi*0.5)*(nphi-dr);
      TVector3 dc = rotation->RotateZ(dphi*0.5)*nphi;
      TVector3 dd = rotation->RotateZ(dphi*0.5)*(nphi-dr);

      TVector3* a=new TVector3( (*center)+da );
      TVector3* b=new TVector3( (*center)+db );
      TVector3* c=new TVector3( (*center)+dc );
      TVector3* d=new TVector3( (*center)+dd );
      
      
      rotation->RotateZ(-dphi);
      
      RobinHoodTriangle *T1=new RobinHoodTriangle( a, b, c, 0);
      RobinHoodTriangle *T2=new RobinHoodTriangle( b, c, d, 0);
 
      _list.Add(T1);
      _list.Add(T2);
    }
    nphi -= dr;
    rotation->RotateZ(-dphi*0.5);
  }

  cout << "Total triangles from disk = " << _list.GetSize() << endl;

}
  


void 
RobinHoodTriangleMaker::sphere( TVector3 *center, TVector3 *norm, double radius, double height, int nTriangles) {}
  

  
void 
RobinHoodTriangleMaker::rectangle( TVector3 *center, TVector3 *side1, TVector3 *side2, int nTriangles, TString opt) {

  assert(  nTriangles>0 && center && side1 && side2 );

  if (!opt.Contains("add")) _list.Clear();

  int nt_side1  = sqrt( nTriangles * side1->Mag() / side2->Mag() );
  int nt_side2  = nTriangles/nt_side1;
  if (nt_side1<1) nt_side1=1;
  if (nt_side2<1) nt_side2=1;
  cout << nt_side1 << "  " <<  nt_side2 << endl;

  TVector3 ds1 = (*side1)*(1./nt_side1);
  TVector3 ds2 = (*side2)*(1./nt_side2);

  for (int it=0; it<nt_side1; it++) {

    for (int jt=0; jt<nt_side2; jt++) {

      TVector3* a=new TVector3( (*center) + it*ds1     + jt*ds2 );
      TVector3* b=new TVector3( (*center) + (it+1)*ds1 + jt*ds2 );
      TVector3* c=new TVector3( (*center) + it*ds1     + (jt+1)*ds2 );
      TVector3* d=new TVector3( (*center) + (it+1)*ds1 + (jt+1)*ds2 );
      
      RobinHoodTriangle *T1=new RobinHoodTriangle( a, b, d, 0);
      RobinHoodTriangle *T2=new RobinHoodTriangle( a, c, d, 0);

      _list.Add(T1);
      _list.Add(T2);

    }
  }

  cout << "Total triangles from rectangle = " << _list.GetSize() << endl;
}





void
RobinHoodTriangleMaker::dumpTriangles(const char *fname) {
  ofstream fout;
  fout.open(fname);

  TIter next(&_list);
  while ( RobinHoodTriangle *tri = (RobinHoodTriangle*) next()  ) {
    for (int i=0; i<3; i++) {
      fout << tri->getVertex(i)->X() << "  " << tri->getVertex(i)->Y() << "  " << tri->getVertex(i)->Z() << "     ";
    }
    fout << endl;
  }


  fout.close();

}


void 
RobinHoodTriangleMaker::addConductorToRhc(const char* rhcFile, double volts, const char * geoFile) {

  
  ofstream fout;
  fout.open(rhcFile, ios::app);

  fout << endl;
  fout << "o1" << endl; // o0=point charge   o1=metal    o2=dielectric
  fout << 1 << endl; // object is active
  fout << 1 << endl; // double accuracy
  fout << 1 << endl; // 0=fixed charge    1=fixed potential
  fout << volts << endl; // value of potential
  fout << 1 << endl; // dummy
  fout << geoFile << endl; // geometry file with triangles
  fout << -1 << endl;
  fout << -1 << endl;
  fout <<  0 << endl;
  fout << (const char*)(TString(geoFile)+".tmp" ) << endl;

  fout.close();


}
