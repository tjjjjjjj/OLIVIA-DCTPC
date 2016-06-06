#ifndef DMTPC_ITERATORS_HH
#define DMTPC_ITERATORS_HH

#include <iterator>
#include "DmtpcDataset.hh" 
#include "DmtpcSkimDataset.hh" 
#include "TH2.h" 
#include "TTree.h"

using namespace std; 

class DmtpcDatasetImageIterator : public iterator<input_iterator_tag, const TH2*> 
{
  private: 
    DmtpcDataset * d; 
    int c; 
    int iev; 

  public: 
    DmtpcDatasetImageIterator(DmtpcDataset * dataset, int camera = 0, int init = 0) 
      { d = dataset; c = camera; iev = init; } 
    DmtpcDatasetImageIterator(const DmtpcDatasetImageIterator & oth) 
      { d = oth.d; c = oth.c; iev = oth.iev; } 
    DmtpcDatasetImageIterator & operator++() 
      { iev++; return *this;} 
    DmtpcDatasetImageIterator operator++(int) 
      { DmtpcDatasetImageIterator tmp(*this); operator++(); return tmp;} 
    bool operator==(const DmtpcDatasetImageIterator & oth)
      { return d == oth.d && c == oth.c && iev == oth.iev;} 
    bool operator!=(const DmtpcDatasetImageIterator & oth)
      { return d != oth.d || c != oth.c || iev != oth.iev;} 
    const TH2 * operator*() 
      { d->getEvent(iev); return d->event()->ccdData(c); } 
}; 

class DmtpcSkimDatasetImageIterator : public iterator<input_iterator_tag, const TH2*> 
{
  private: 
    DmtpcSkimDataset * d; 
    int c; 
    int iev; 

  public: 
    DmtpcSkimDatasetImageIterator(DmtpcSkimDataset * dataset, int camera = 0, int init = 0) { d = dataset; c = camera; iev = init; } 
    DmtpcSkimDatasetImageIterator(const DmtpcSkimDatasetImageIterator & oth) { d = oth.d; c = oth.c; iev = oth.iev; } 
    DmtpcSkimDatasetImageIterator & operator++()  { iev++; return *this;} 
    DmtpcSkimDatasetImageIterator operator++(int)  { DmtpcSkimDatasetImageIterator tmp(*this); operator++(); return tmp;} 
    bool operator==(const DmtpcSkimDatasetImageIterator & oth) { return d == oth.d && c == oth.c && iev == oth.iev;} 
    bool operator!=(const DmtpcSkimDatasetImageIterator & oth) { return d != oth.d || c != oth.c || iev != oth.iev;} 
    const TH2 * operator*() { d->getEvent(iev); return d->event()->cluster(c)->getImage(); } 
}; 

template <class T = double> 
class TreeIterator : public iterator<input_iterator_tag, T> 
{
  private:  
    TTree *t; 
    T * ptr; 
    int ev; 
  public: 
    TreeIterator(TTree * tree, T *branch_pointer, int i = 0) { t = tree; ptr = branch_pointer; ev = i; } 
    TreeIterator(const TreeIterator & oth) { t = oth.t; ptr = oth.ptr; ev = oth.ev; } 
    TreeIterator & operator++() { ev++; return *this;}
    TreeIterator operator++(int) { TreeIterator tmp(*this); ev++; return tmp;}
    bool operator==(const TreeIterator & oth) { return t == oth.t && ptr == oth.ptr && ev == oth.ev; } 
    bool operator!=(const TreeIterator & oth) { return !(t == oth.t && ptr == oth.ptr && ev == oth.ev); } 
    T operator *() { t->GetEntry(ev); return *ptr; } 
};



#endif

