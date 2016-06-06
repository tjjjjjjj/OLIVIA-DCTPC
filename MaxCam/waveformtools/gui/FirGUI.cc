#include "GuiFirFrame.hh"
#include "TApplication.h"
#include "TGWindow.h"
#include "TStyle.h"
#include "TROOT.h"

int main(int nargs, char** args)
{
  TApplication app("app",0,0);
  GuiFirFrame* fr = new GuiFirFrame(gClient->GetRoot());
  app.Run();
  return 0;
}
