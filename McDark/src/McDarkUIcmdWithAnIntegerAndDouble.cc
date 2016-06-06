
#include "McDarkUIcmdWithAnIntegerAndDouble.hh"
#include "G4Tokenizer.hh"
#include "G4UnitsTable.hh"
#include <sstream>

McDarkUIcmdWithAnIntegerAndDouble::McDarkUIcmdWithAnIntegerAndDouble
(const char * theCommandPath,G4UImessenger * theMessenger)
:G4UIcommand(theCommandPath,theMessenger)
{
  G4UIparameter * iCameraNumber = new G4UIparameter('i');
  SetParameter(iCameraNumber);
  iCameraNumber->SetParameterName("Integer_Paramter");
 // iCameraNumber->SetParameterRange("CameraNumber >= 0");

  G4UIparameter * doubleParam = new G4UIparameter('d');
  SetParameter(doubleParam);
  doubleParam->SetParameterName("Double_Parameter");
 // noise_sigma->SetParameterRange("Noise_sigma > 0");
}

/*
G4int getCameraNumber(const char* newValue) {
	G4Tokenizer tok(newValue);
	G4int num = tok();
	return num;
}

G4int getPixelsPerSide(const char* newValue) {
	//std::istringstream is(newValue);
	G4Tokenizer tok(newValue);
	tok();
	G4int num;
	is >> num >> num;
	return num;
}

G4int getPixelsPerBinLength(const char* newValue) {
	std::istringstream is(newValue);
	G4int num;
	is >> num >> num >> num;
	return num;
}

G4double getViewFieldLengthX(const char* newValue) {
	std::istringstream is(newValue);
	G4double dump;
	is >> dump >> dump >> dump;
	G4double lengthX;
	is >> lengthX;
	//G4double dump;
	char[30] units;
	is >> dump >> dump >> dump;
	is >> units;
	G4String dUnit = units;
	return (lengthX * valueOf(dUnit));
}

G4double getViewFieldLengthY(const char* newValue) {
	std::istringstream is(newValue);
	G4double dump;
	is >> dump >> dump >> dump >> dump;
	G4double lengthY;
	is >> lengthY;
	char[30] units;
	is >> dump >> dump;
	is >> units;
	G4String dUnit = units;
	return (lengthY * valueOf(dUnit));
}

G4double getOffsetX(const char* newValue) {
	std::istringstream is(newValue);
	G4double dump;
	is >> dump >> dump >> dump >> dump >> dump;
	G4double offsetX;
	is >> offsetX;
	char[30] units;
	is >> dump;
	is >> units;
	G4String dUnit = units;
	return (offsetX * valueOf(dUnit));
}

G4double getOffsetY(const char* newValue) {
	std::istringstream is(newValue);
	G4double dump;
	is >> dump >> dump >> dump >> dump >> dump >> dump;
	G4double offsetY;
	is >> offsetY;
	char[30] units;
	is >> units;
	G4String dUnit = units;
	return (lengthY * valueOf(dUnit));
} */
