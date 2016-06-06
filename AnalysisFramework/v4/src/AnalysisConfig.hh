#ifndef ANALYSIS_CONFIG_HH
#define ANALYSIS_CONFIG_HH

/* AnalysisConfig.hh 
  
   author: ACK 
   
   Parsing for analysis configuration
*/

#include <string>
#include <strings.h>
#include <ostream>
#include <iostream>
#include <map>
#include "../../../MaxCam/SimpleConfig.hh"
#include "TString.h"

class AnalysisConfig :  public SimpleConfig
{
   public:
      AnalysisConfig();
      AnalysisConfig(const char *f);
      
      /* Prints out the current config to the given output 
	 stream. */
      void print(std::ostream & out = std::cout); 
      
      /* Getter/setter */
      char* getTree() const {return _tree;}
      void setTree(char* tree) {_tree=tree;}

      char* getRangeCalFileName() const {return _rangefilename;}
      void setRangeCalFileName(char* rangefilename) {_rangefilename=rangefilename;}
      
      double getRangeCal(const char* serial) {return _rangecal.count(std::string(serial)) ? _rangecal[std::string(serial)] : -1;}

      double getEnergyCal(const char* serial) {return _energycal.count(std::string(serial)) ? _energycal[std::string(serial)] : -1;}
      void setEnergyCal(const char* serial, double value) {_energycal[serial]=value;}

      TString getWormFileName(const char* serial) {return _wormfilename.count(std::string(serial)) ? _wormfilename[std::string(serial)] : "";}

      double getWormCutVal(const char* serial) {return _wormcutval.count(std::string(serial)) ? _wormcutval[std::string(serial)] : -1;}


   private:
      void parsePair(std::string key, std::string value, int n);
      
      char* _tree;
      char* _rangefilename;
      char** _cameras;
      std::map<std::string, double> _rangecal;
      std::map<std::string, double> _energycal;
      std::map<std::string, TString> _wormfilename;
      std::map<std::string, double> _wormcutval;

      void setDefaults();
      
      
      
};

#endif
