#include "DmtpcDataConverter.hh"
#include <iostream>

/** Front end for DmtpcDataConverter::convertDataFile */ 

int main(int nargs, char ** args)
{

  if (nargs ==  3) 
    DmtpcDataConverter::convertDataFile(args[1],args[2]); 
  else
    std::cout << "Usage: DataConverter infile outfile " << std::endl; 

  return 0; 
}

