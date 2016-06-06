#ifndef DMTPC_ROOT_TOOLS_HH
#define DMTPC_ROOT_TOOLS_HH

class TGraph;
class TChain; 
class TH2; class TH1; class TCanvas; 
class TArray; 

//typedef unsigned int size_t; 
#include <unistd.h>
#include <iostream>

/** Namespace of functions for ROOT related tasks that are absent from older ROOT releases but that are useful...
e.g. TGraph.Integral() is not available in 5.24.
*/


namespace DmtpcRootTools {

  /** Adds files from a text file list to chain 
   * */ 
  size_t addFilesToChain(TChain * ch, const char * file, 
                        const char * prefix = 0, const char * det_tag = 0, const char * suffix = 0); 

  double integral(TGraph *gr, int first, int last);

  //Monochromatic Schemes
  void setColorGrayscale(double p = 1);
  void setColorGreen(double p = 1);
  void setColorBlue(double p = 1);
  void setColorRed(double p = 1);
  void setColorCyan(double p = 1);
  void setColorMagenta(double p = 1);
  void setColorYellow(double p = 1);
  //Similar to ROOT palette 1
  void setColorStandard1();
  //Found in ROOT documentation somewhere
  void setColorPalette1();//Blue
  void setColorPalette2();//Pink/green
  void setColorPalette3();//green/pink
  //Matlab Color Schemes
  void setColorHot();
  void setColorCool(double p =1);
  void setColorJet();
  void setColorCopper();
  void setColorBone();
  void setColorWinter(double p = 1);
  void setColorSpring(double p = 1);
  void setColorSummer(double p = 1);
  void setColorAutumn(double p = 1);
  //Random other schemes (school colors, whatever else we want)
  void setColorMIT(double p = 1);
  void setColorColumbia(double p = 1);
  void setColorBU(double p = 1);

  //Random colors
  void setColorRandom(int n = 5);

  /** 
   * Checks if a file contains a given tree. Returns tree if it does, false if it doesn't or if the
   * file can't be opened. 
   *
   * \param file the name of the file to check 
   * \param treename the name of the TTree to check for
   * \return true if the file contains the tree
   */ 

  bool checkIfFileContainsTree(const char * file, const char * treename); 

  /** Create a new TH2 of the desired type (C,S,I,F,D) */ 
  TH2 * newTH2(char type, const char * name, const char * title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax);
  TH1 * newTH1(char type, const char * name, const char * title, int nbinsx, double xmin, double xmax, TH1 * ptr = 0); 
  TH2 * newTH2StealType(const TH2* type_giver, const char * name, const char * title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax);
  TH2 * newTH2StealTypeAndSize(const TH2* type_giver, const char * name, const char * title); 
  TH2 * newTH2StealSize(const TH2* size_giver, char type, const char * name, const char * title); 
  TH1 * newTH1StealSize(const TH1* size_giver, char type, const char * name, const char * title); 
  TH1 * newTH1StealType(const TH1* type_giver, const char * name, const char * title, int nbinsx, double xmin, double xmax, TH1 * placement_ptr = 0); 
  unsigned int TH1GetNbins(const TH1* hist); 

  template<class T>
  bool checkNaN(size_t N, const T * v) 
  {
    for (size_t i = 0; i < N; i++)
    {
      if (isnan(v[i]))
      {
        return true; 
      }
    }
    return false; 
  }

  template <class T>
  void printVec(size_t N, const T * v) 
  {
    std::cout << "[ "; 
    for (size_t i = 0; i < N; i++)
    {
      std::cout << v[i];   
      std::cout << (i < N-1 ? " , " : " ]"); 
    }
    std::cout << std::endl; 
  }


  /** Stalls until the given canvas is closed */ 
  void waitForCanvasClose(TCanvas * c); 

  char getType(const TH1* in); 
  char getType(const TArray* in); 


};  
#endif
