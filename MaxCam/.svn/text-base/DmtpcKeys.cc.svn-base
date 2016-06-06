#include "DmtpcKeys.hh"
#include "MaxCamImage.hh"
#include "TDatime.h"
#include "MaxCamConfig.hh"
#include "ScopeDataInfo.hh"
//#include "ScopeWaveform.hh"

#include "TMath.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"

#include <iostream>
#include <fstream>
using namespace std;

DmtpcKeys::DmtpcKeys(TString keyfilename, TString rawdatafiles, 
		     TString key,TString newdir, 
		     int nreqmod, TString* reqmod, bool& pass)
{

   _keyfilename = keyfilename;
   _key = key;
   _keynewdir = newdir;
   pass = 1;
   ifstream infile(rawdatafiles,ios::in);
   ifstream keyfile(keyfilename,ios::in);
   keyfile >> _rootdirname;
   
   TString keydirtemp, keytemp;
   int nk=0;
   _isSkim=false;
   while(keyfile >> keydirtemp)
   {
      keyfile >> keytemp;
      if(keytemp == "skim") {_isSkim=true; _skimindex=nk;}
      _keydir.push_back(keydirtemp); 
      _keys.push_back(keytemp);
      nk++;
   }
   
   if(_keys.size() == 1) cout << "There is 1 key. \n";
   else cout << "There are " << _keys.size() << " keys. \n";
   
   TString filenametemp;
   while(infile >> filenametemp)
   {
      _files.push_back(filenametemp);
   }
   
   if(_files.size() == 1) cout << "There is 1 file. \n";
   else cout << "There are " << _files.size() << " files. \n";
   keyfile.close();
   infile.close();

   if(contains(_key))
   {
      cout << "Key already exists in key file; Aborting \n";
      pass = 0;
   }
   if(!(contains(nreqmod,reqmod)))
   { 
      cout << "Some required keys not found; Aborting! \n";
      pass = 0;
   }
   
   if(pass) writeKey();
   
}

DmtpcKeys::~DmtpcKeys()
{

}

bool DmtpcKeys::contains(TString testkey)
{
   for(int i=0; i<int(_keys.size()); i++)
   {
      if(_keys[i] == testkey) return true;
   }
   
   return false;
}

bool DmtpcKeys::contains(int n, TString* testkeys)
{
   bool foundall = true;
   for(int i=0; i<n; i++)
   {
      bool found = contains(testkeys[i]);
      if(!found)
      {
	 cout << "Required key " << testkeys[i] << " missing. \n";
	 foundall = false;
      }
   }
   
   return foundall;
}

void DmtpcKeys::writeKey()
{
   ofstream keyfileo(_keyfilename,ios::app);
   keyfileo << _keynewdir << "\t" << _key << "\n";
   keyfileo.close();
}

void DmtpcKeys::addFriends(TTree* tree, int f)
{
   if(!_isSkim)
   {
      for(int i=0; i<int(_keys.size()); i++)
      {
	 if(!(_keys[i].BeginsWith("sk")))
	 {
	    TString tempfname = _files[f];
	    tempfname.ReplaceAll(".",_keys[i]+".");
	    cout << _keys[i] << "fname: " << _keydir[i]+tempfname << endl;
	    tree->AddFriend(_keys[i],_keydir[i]+tempfname);
	 }
      }
   }
   if(_isSkim)
   {
      for(int i=0; i<int(_keys.size()); i++)
      {
	 if(_keys[i].BeginsWith("sk") && _keys[i]!="skim")
	 {
	    TString tempfname = _files[f];
	    tempfname.ReplaceAll(".",_keys[i]+".");
	    cout << _keys[i] << "fname: " << _keydir[i]+tempfname << endl;
	    tree->AddFriend(_keys[i],_keydir[i]+tempfname);
	 }
      }
   }
}

bool DmtpcKeys::openOutputFile(TFile* routfile, int i)
{
   TString routfilename = _files[i];
   routfilename.ReplaceAll(".root",_key+".root");
   routfile = new TFile(_keynewdir+routfilename,"CREATE");

   if(routfile->IsOpen())
   {
      return true;
   }
   else
   {
      cout << "Output file already exists; Aborting!" << endl;
      return false;
   }
   
}

TTree* DmtpcKeys::getBaseTree(int f)
{  
   TTree* basetree = 0;

   cout << "Skim:" << _isSkim << endl;

   if(!_isSkim)
   {
      cout << _rootdirname+_files[f] << endl;
      
      TFile* file = new TFile(_rootdirname+_files[f]);
      basetree = (TTree*)file->Get("dmtpc");
   }
   if(_isSkim)
   {
      TString skimfile = _files[f];
      skimfile.ReplaceAll(".root","skim.root");
      cout << _keydir[_skimindex]+skimfile << endl;

      TFile* file = new TFile(_keydir[_skimindex]+skimfile);
      basetree = (TTree*)file->Get("skim");
   }

   return basetree;
   
}
