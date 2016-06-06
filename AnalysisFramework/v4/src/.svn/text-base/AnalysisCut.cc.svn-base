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
   _once = false;

}

AnalysisCut::AnalysisCut(TString name, bool(*fcn)(DmtpcSkimEvent*, int, int, vector<void*> ))
{
   SetName(name);
   _type = 1;
   _vfunc = fcn;
   _once = false;
   
}


AnalysisCut::~AnalysisCut()
{

   //if(_params) delete [] _params;
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
   vector< int > tracks;
   for(int i=0; i<ev->ncamera(); i++)
   {
      tracks.clear();
      bool onceval=false;
      for(int j=0; j<15; j++)
      {
	 if(_once)
	 {
	    if(!j)
	    {
	       if(_type==0)
	       {
		  onceval = (_func)(ev,i,j,_params);
		  tracks.push_back(onceval);
	       }
	       if(_type==1)
	       {
		  onceval = (_vfunc)(ev,i,j,_vparams);
		  tracks.push_back(onceval);
	       }
	    }
	    else
	    {
	       tracks.push_back(onceval);
	    }
	 }
	 else
	 {
	    if(_type==0)
	       tracks.push_back((_func)(ev,i,j,_params));
	    if(_type==1)
	       tracks.push_back((_vfunc)(ev,i,j,_vparams));
	 }
      }
      _pass.push_back(tracks);
   }

}
