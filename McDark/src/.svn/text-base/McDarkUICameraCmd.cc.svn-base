
#include "McDarkUICameraCmd.hh"
#include "G4Tokenizer.hh"
#include "G4UnitsTable.hh"
#include <sstream>

McDarkUICameraCmd::McDarkUICameraCmd
(const char * theCommandPath,G4UImessenger * theMessenger)
:G4UIcommand(theCommandPath,theMessenger)
{
  G4UIparameter * iCameraNumber = new G4UIparameter('i');
  SetParameter(iCameraNumber);
  iCameraNumber->SetParameterName("CameraNumber");
  iCameraNumber->SetParameterRange("CameraNumber >= 0");

  G4UIparameter * iPixelsPerSide = new G4UIparameter('i');
  SetParameter(iPixelsPerSide);
  iPixelsPerSide->SetParameterName("PixelsPerSide");
  iPixelsPerSide->SetParameterRange("PixelsPerSide >= 1");

  G4UIparameter * iPixelsPerBinSide = new G4UIparameter('i');
  SetParameter(iPixelsPerBinSide);
  iPixelsPerBinSide->SetParameterName("PixelsPerBinSide");
  iPixelsPerBinSide->SetParameterRange("PixelsPerBinSide >= 1");

  G4UIparameter * x_length = new G4UIparameter('d');
  SetParameter(x_length);
  x_length->SetParameterName("SizeOfViewField_x");
  x_length->SetParameterRange("SizeOfViewField_x > 0");

  G4UIparameter * y_length = new G4UIparameter('d');
  SetParameter(y_length);
  y_length->SetParameterName("SizeOfViewField_y");
  y_length->SetParameterRange("SizeOfViewField_y > 0");

  G4UIparameter * x_offset = new G4UIparameter('d');
  SetParameter(x_offset);
  x_offset->SetParameterName("X_offset");

  G4UIparameter * y_offset = new G4UIparameter('d');
  SetParameter(y_offset); 
  y_offset->SetParameterName("Y_offset");

  G4UIparameter * untParam = new G4UIparameter('s');
  SetParameter(untParam);
  untParam->SetParameterName("Unit");
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
