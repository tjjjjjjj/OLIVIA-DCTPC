#include <malloc.h>
#include "AnalysisResult.hh"
#include <string.h>
#include <iostream>
ClassImp(AnalysisResult);

AnalysisResult::AnalysisResult()
{
  _data = 0; 
}

AnalysisResult::~AnalysisResult()
{
  if (_data!=0)
    free(_data); 
}

AnalysisResult::AnalysisResult(int ncam, int nev, UInt_t * data)
{
  _ncam = ncam; 
  _nev = nev;
  _arr_size = 15*_nev*_ncam; 
  
  _data = (UInt_t *) malloc(sizeof(UInt_t) * _arr_size); 
  std::cout << "arr_size: " << _arr_size << std::endl;
  if (data!=NULL)
  {
    for (int i = 0; i < _arr_size; i++)
    {
      _data[i] = data[i];
    }
  } 
}

AnalysisResult::AnalysisResult(const AnalysisResult & other)
{
  _ncam = other._ncam;
  _nev = other._nev; 
  _arr_size = other._arr_size;
  _data = (UInt_t *) malloc(sizeof(UInt_t) * _arr_size); 
  memcpy(_data,other._data,sizeof(UInt_t) * _arr_size); 
}
