#ifndef ANALYSIS_CUT_HH
#define ANALYSIS_CUT_HH

#include "TClonesArray.h"
#include "TObject.h"
#include "TMath.h"
#include "Math/ParamFunctor.h"
#include "TNamed.h"
#include "../../../MaxCam/DmtpcSkimEvent.hh"
#include <ostream>
#include <iostream>
#include <vector>

class AnalysisCut : public TNamed
{

   public:

      /**
	 Default constructor
      */
      AnalysisCut();
      /**
	 Constructor for any cut which only depends on numbers; type "0"
	 \param name Name of the cut
	 \param fcn Function taking a DmtpcSkimEvent, the camera index, the track index, 
	 and an array of doubles for cut parameters. The function MUST be in this format
	 \param the number of parameters
      */
      AnalysisCut(TString name, bool (*fcn)(DmtpcSkimEvent* , int , int, double*), int nparameters);
      /**
	 Constructor for any cut which depends on numbers, or objects; type "1"
	 \param name Name of the cut
	 \param fcn Function taking a DmtpcSkimEvent, the camera index, the track index, 
	 and a vector of pointers. The pointers can be to anything, but the user is 
	 responsible for keeping track of the types of the pointers and using them 
	 appropriately, including deleting them, as required.
      */
      AnalysisCut(TString name, bool (*fcn)(DmtpcSkimEvent* , int , int, vector<void*>));

      /**
	Default destructor
      */
      ~AnalysisCut();

      /**
	 Returns the boolean pass value of a subevent
	 \param c camera index
	 \param t track index
      */
      bool passes(int c, int t) {return _pass[c][t];}
      /**
	 Returns the boolean pass value of a subevent
	 \param c camera index
	 \param t track index
      */
      vector < vector <bool> > passes() {return _pass;}

      /** 
	  Sets a parameter if all parameters are doubles
	  \param i index of parameter to set
	  \param val value of parameter
       */
      void setParameter(int i, double val) {if(_type==0)_params[i]=val;}
      /**
	 Gets the ith parameter if all parameters are doubles
	 \param i index of parameter to get
	 \return parameter[i]
      */
      double getParameter(int i) {if(_type==0)return _params[i];}
      /**
	 Sets the parameters in the case of passing pointers
	 \param params the vector<void*> of parameters
      */
      void setParameters(vector<void*> params);
      /**
	 evaluates the cut on an event by looping over the number of cameras in the event, 
	 and up to 15 tracks. 
	 \param ev the DmtpcSkimEvent to evaluate.
      */
      void evaluateEvent(DmtpcSkimEvent* ev);

   private:
      bool (* _func)(DmtpcSkimEvent*,int,int,double*); //! pointer to funtion
      bool (* _vfunc)(DmtpcSkimEvent*,int,int,vector<void*>); //! pointer to funtion
      int _npar; //! number of parameters
      double* _params; //! array for parameters
      vector<void*> _vparams; //! vector for parameters
      vector < vector <bool> > _pass; //nested vector of passes in cam, track
      int _type; //! type of cut (double* = 0, vector< void* >=1)
      
   ClassDef(AnalysisCut,1);

};


#endif
