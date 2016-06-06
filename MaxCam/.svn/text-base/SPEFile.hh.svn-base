
#include <iostream>
#include <fstream>
#include "TH2F.h"

using namespace std;


class SPEFile {

  struct header {
    short xpix;
    char junk1[64];
    short data_type;
    char junk2[546];
    short ypix;
    char junk3[788];
    long nimages;
  };
  
  struct header myheader;
  
  int nbytesperim;
  
  fstream myfile;

public:

  SPEFile (const char* filename);

  int Getxpix();

  int Getypix();

  int Getdata_type();
  
  long int Getnimages();

  TH2F* GetImagetoROOT(int a);


};




