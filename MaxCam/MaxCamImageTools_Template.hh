#ifndef MAXCAM_IMAGE_TOOLS_TEMPLATE_HH
#define MAXCAM_IMAGE_TOOLS_TEMPLATE_HH

#include "DmtpcRootTools.hh"
#include "TF1.h" 
#include "TH1.h" 
#include <vector>
/** Templated functions go here since they take up so much space  **/ 

namespace MaxCamImageTools
{


   /** returns a vector of the fit parameters, as well as chi square*/ 
   template <class T, class InputIterator>
   vector<TH2*> * valueFit(InputIterator first, InputIterator last, int total = -1, TF1 * fn = NULL, const char * fitpars = "NQ", int maxalloc = 1<<16, bool verbose = false) 
   {
     const TH2 * hist = *first; 
     if (!hist) 
     {
       return 0; 
     }

     int npar = fn ? fn->GetNpar() : 3; 
     vector<TH2*> *answer = new vector<TH2*>(npar+1); 

     static char * default_names[] = {"A","mu","sigma"}; 
     for (int i = 0; i < npar; i++)
     {
       TString name = fn ? fn->GetParName(i) : default_names[i]; 
       (*answer)[i] = DmtpcRootTools::newTH2StealSize(hist,'D',name,name); 
     }

     (*answer)[npar] = DmtpcRootTools::newTH2StealSize(hist,'D',"chisq","chisq"); 

     int nx = (*answer)[0]->GetNbinsX(); 
     int ny = (*answer)[0]->GetNbinsX(); 
     int nbins = nx * ny;

     if (total < 0)
     {
       total = 0; 
       for (InputIterator it(first); it!=last; it++)
       {
          total++; 
       }
     }

     if (verbose)
     {
        cout << "Allocating " << sizeof(T) * maxalloc * total << " bytes " << endl; 
     }


     vector<vector<T> > values (maxalloc, vector<T>(total)); 
     bool leftover = nbins % maxalloc != 0;
    
     if (verbose)
       cout << "leftover: " << leftover<< endl;

     int ntimes = nbins/maxalloc + leftover; 
    
     if (verbose)
       cout << "ntimes: " <<ntimes << endl;

     for (int i = 0; i < ntimes ; i++)
     {
      if (verbose) cout << i << endl; 
      int idx = 0; 
      int max = (leftover && i == ntimes - 1) ? nbins % maxalloc : maxalloc; 

      vector<T> maxes(maxalloc); 
      vector<T> mins(maxalloc); 

      for (InputIterator it(first); it!=last; it++)
      {
        hist = *it; 
        if (verbose) cout << i << endl; 
        for (int j = 0; j < max ; j++)
        {
            int bin = i * maxalloc + j; 
            int xbin = bin % nx + 1; 
            int ybin = bin / nx + 1; 
            values[j][idx] = (T) hist->GetBinContent(xbin,ybin); 

            if (!idx)
            {
                maxes[j] = values[j][idx];
                mins[j] = values[j][idx];
            }
            else
            {
              if (values[j][idx] > maxes[j]) maxes[j] = values[j][idx]; 
              if (values[j][idx] < mins[j]) mins[j] = values[j][idx]; 
            }
//            cout << j << " " <<  idx << endl; 
        }
        idx++; 
      }

      for (int j = 0; j < max; j++)
      {
        int bin = i * maxalloc + j; 
        int xbin = bin % nx + 1; 
        int ybin = bin / nx + 1; 

        TH1 * fithist = DmtpcRootTools::newTH1StealType(hist,"fithist","fithist",int(maxes[j]-mins[j] + 1), mins[j],maxes[j]); 
        for (unsigned ii = 0; ii < values[j].size(); ii++)
        {
          fithist->Fill(values[j][ii]); 
        }

        bool deletefn = fn; 
        if (!fn)
        {
           fn = new TF1("tmpfn","gaus",mins[j],maxes[j]); 
        }
        fithist->Fit(fn,fitpars); 
        for (int par = 0; par < npar; par++)
        {
          (*answer)[par]->SetBinContent(xbin,ybin,fn->GetParameter(par)); 
        }

        (*answer)[npar]->SetBinContent(xbin,ybin,fn->GetChisquare()); 

        delete fithist; 
        if (deletefn) 
        {
          delete fn; 
          fn = 0; 
        }
      }
    }

    if (verbose)
    {
      cout << "Done!" <<endl; 
    }


    return answer; 


 }; 

   template <class T, class InputIterator>
   TH2 * histStackNthElement(InputIterator first, InputIterator last, int nth, int total = -1, int maxalloc = 1024, bool verbose = false) 
   {
    const TH2 * hist = *first; 

    if (!hist)
    {
        return 0; 
    }

    TH2 * answer = DmtpcRootTools::newTH2StealType(hist, "hnth","hnth", hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(), 
                                                                        hist->GetNbinsY(), hist->GetYaxis()->GetXmin(), hist->GetYaxis()->GetXmax()) ; 
    int nx = answer->GetNbinsX(); 
    int ny = answer->GetNbinsY(); 
    int nbins = answer->GetNbinsX() * answer->GetNbinsY(); 


    if (total < 0)
    {
       total = 0; 
       for (InputIterator it(first); it!=last; it++)
       {
          total++; 
       }
    }


    if (nth > total) return 0; 

    if (verbose)
    {
      cout << "Allocating " << sizeof(T) * maxalloc * total << " bytes " << endl; 
    }

    vector<vector<T> > values ( maxalloc, vector<T>(total)); 

    bool leftover = nbins % maxalloc != 0;
    
    if (verbose)
      cout << "leftover: " << leftover<< endl;
    int ntimes = nbins/maxalloc + leftover; 
    
    if (verbose)
      cout << "ntimes: " <<ntimes << endl;


    for (int i = 0; i < ntimes ; i++)
    {
      if (verbose) cout << i << endl; 
      int idx = 0; 
      int max = (leftover && i == ntimes - 1) ? nbins % maxalloc : maxalloc; 
      for (InputIterator it(first); it!=last; it++)
      {
        hist = *it; 
        for (int j = 0; j < max ; j++)
        {
            int bin = i * maxalloc + j; 
            int xbin = bin % nx + 1; 
            int ybin = bin / nx + 1; 
            values[j][idx] = (T) hist->GetBinContent(xbin,ybin); 
//            cout << j << " " <<  idx << endl; 
        }
        idx++; 
      }

      for (int j = 0; j < max; j++)
      {
        int bin = i * maxalloc + j; 
        int xbin = bin % nx + 1; 
        int ybin = bin / nx + 1; 
        std::nth_element(values[j].begin() , values[j].begin() + nth, values[j].end() ) ; 
        T elem = values[j][nth]; 
        answer->SetBinContent(xbin,ybin,elem); 
      }
    }

    if (verbose)
    {
      cout << "Done!" <<endl; 
    }


    return answer; 
  }

};
#endif
