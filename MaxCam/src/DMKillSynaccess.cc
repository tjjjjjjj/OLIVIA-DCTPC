#include <iostream>
using std::cout;
using std::endl;

#include "TROOT.h"
#include "TSystem.h"
#include "TSocket.h"
#include "MaxCamChannel.hh"
#include <cmath>

static char * default_host = "aragorn.lns.mit.edu "; 

#define NSOCKETS 8 
#define SLEEPTIME 300

bool debug=true;


int main(int argc, char* argv[]) {
  
  cout << "Open socket" << endl;
  //TSocket *redsock = new TSocket("powmgmt.lns.mit.edu", 23);
  char * host; 
  if (argc>1) host = argv[1];
  else host = default_host; 
  TSocket *redsock = new TSocket(host, 23);
  if (!redsock->IsValid()) {
      cout << "Socket not valid " << endl;
      return 1;
  } else {
    if (debug) {cout<<"Successfully opened socket"<<endl;}
  }

  // check initial message
  char rmsg[2048];
  cout << "Check initial message"<<endl;
  int nbytes = redsock->RecvRaw( rmsg, 2048, kDontBlock);
  cout << rmsg << endl;
  nbytes = redsock->RecvRaw( rmsg, 2048, kDontBlock);
  cout << rmsg << endl;

  /*
  //Close all valves
  char msg[] = "pset 0 0\n";
  for (int i = 1; i <= NSOCKETS; i++)
  {
    *(msg+5) = (char) (48+i); 
    cout << msg << endl; 
    redsock->SendRaw(msg, 10); // make sure gas is closed
    gSystem->Sleep(SLEEPTIME);
  }
  */

  cout << "Closing socket"<<endl;
  redsock->Close();
  delete redsock;
  return 0;
}
