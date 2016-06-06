#include "AnalysisCut.hh"
#include <stdlib.h>
#include <assert.h>
#include <cmath>

ClassImp(AnalysisCut);

AnalysisCut::AnalysisCut()
{

}


AnalysisCut::AnalysisCut(TString name, bool(*fcn)(DmtpcSkimEvent*, int, int, double*), int nparameters)
{
   SetName(name);
   _type = 0;
   _func = fcn;
   _npar = nparameters;
   _params = new double[_npar];

}

AnalysisCut::AnalysisCut(TString name, bool(*fcn)(DmtpcSkimEvent*, int, int, vector<void*> ))
{
   SetName(name);
   _type = 1;
   _vfunc = fcn;
   
}


AnalysisCut::~AnalysisCut()
{

   if(_params) delete [] _params;
}

void AnalysisCut::setParameters(vector<void*> params)
{
   for(int i=0; i<int(params.size()); i++)
   {
      _vparams.push_back(params[i]);
   }
}

void AnalysisCut::evaluateEvent(DmtpcSkimEvent* ev)
{
   _pass.clear();
   vector< bool > tracks;
   for(int i=0; i<ev->ncamera(); i++)
   {
      tracks.clear();
      for(int j=0; j<15; j++)
      {
	 if(_type==0)
	    tracks.push_back((_func)(ev,i,j,_params));
	 if(_type==1)
	    tracks.push_back((_vfunc)(ev,i,j,_vparams));
      }
      _pass.push_back(tracks);
   }

}
