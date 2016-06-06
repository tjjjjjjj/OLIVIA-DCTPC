#ifndef McDarkUIcmdWithAnIntegerAndDouble_H
#define McDarkUIcmdWitHAnIntegerAndDoubleCmd_H 1

#include "G4UIcommand.hh"

// class description:
//  A concrete class of G4UIcommand. The command defined by this class
// takes a double value and a unit string.
//  General information of G4UIcommand is given in G4UIcommand.hh.

class McDarkUIcmdWithAnIntegerAndDouble : public G4UIcommand
{
  public: // with description
    McDarkUIcmdWithAnIntegerAndDouble
    (const char * theCommandPath,G4UImessenger * theMessenger);
};

#endif
