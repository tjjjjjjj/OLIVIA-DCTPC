//Root files
#include "TSystem.h"
#include "TROOT.h"
#include "TString.h"
//MC files
#ifdef EXEC
#include "RunGenerator.hh"
#else
#include "code/RunGenerator.hh"
#endif
using std::cout;
using std::endl;

#ifdef EXEC
int main(int nargs, char* argv[])
{
  if (nargs < 2) {
    cout <<"Usage error: Need filename"<<endl;
    return 1;
  }
  TString filename(argv[1]);
#else
int runFromFile(TString filename)
{
#endif

  //gSystem->Setenv("MCTABLES","../../tables");

  cout <<"Starting run"<<endl;

  RunGenerator rungen;
  std::system("rm waveforms/*");
  rungen.readParametersFromFile(filename);
  rungen.findFileNames();
  rungen.initialize();
  rungen.run();
  rungen.outputTextFile(true,filename);
  cout <<"Finished!"<<endl;
  return 0;
}

