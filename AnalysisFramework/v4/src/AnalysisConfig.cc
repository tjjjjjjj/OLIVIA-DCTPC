#include "AnalysisConfig.hh" 
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>

AnalysisConfig::AnalysisConfig()
{
  setDefaults(); 
}

AnalysisConfig::AnalysisConfig(const char * file)
{
  setDefaults();
  parseFile(file); 
}

void AnalysisConfig::setDefaults()
{
   _tree="skim";
   registerNewString("tree", &_tree);

   _rangecal["081264"]=0.137;
   _rangecal["100439"]=0.168;
   _rangecal["A80334"]=0.171;
   _rangecal["A80333"]=0.137;

   _energycal["081264"]=22;
   _energycal["100439"]=18;
   _energycal["A80334"]=12;
   _energycal["A80333"]=15;

   _wormfilename["081264"]="wr5_081264";
   _wormfilename["100439"]="wr5_100439";
   _wormfilename["A80334"]="wr5_A80334";
   _wormfilename["A80333"]="wr5_A80333";

   _wormcutval["081264"]=0;
   _wormcutval["100439"]=0;
   _wormcutval["A80334"]=0;
   _wormcutval["A80333"]=0;

   _rangefilename="cfg/range.root";
   registerNewString("rangefilename",&_rangefilename);

   
}


void AnalysisConfig::parsePair(std::string key, std::string value, int n) 
{
  if (true) 
  {
    std::cout << "key: \"" << key << "\"\tvalue: \"" << value << "\"" << std::endl; 
  }

  if (key == "rangecal") 
  {
    _rangecal.clear(); 

    if (!strcasecmp(value.c_str(),"none") || !strcasecmp(value.c_str(),""))
    {
       return;   
    }

    int ntokens = 1; 

    for (unsigned int i = 0; i < value.length() ; i++)
    {
      if (value[i] == ',') ntokens++; 
    }

    char buf[80]; 
    double f; 
    if (ntokens == 1)
    {
      sscanf(value.c_str(),"%[^:]:%lf",buf,&f);  
      _rangecal[buf] = f; 
      return;
    }
   
   char * ch; 
   char * strcp = strdup(value.c_str()); 
   ch = strtok(strcp,","); 
   while (ch!=NULL)
   {
     sscanf(ch,"%[^:]:%lf",buf,&f);  
     _rangecal[buf] = f; 
     std::cout << buf << ": " << f << std::endl; 
     ch = strtok(NULL,","); 
   }

   free (strcp); 
   return;
  }

  if (key == "energycal") 
  {
    _energycal.clear(); 

    if (!strcasecmp(value.c_str(),"none") || !strcasecmp(value.c_str(),""))
    {
       return;   
    }

    int ntokens = 1; 

    for (unsigned int i = 0; i < value.length() ; i++)
    {
      if (value[i] == ',') ntokens++; 
    }

    char buf[80]; 
    double f; 
    if (ntokens == 1)
    {
      sscanf(value.c_str(),"%[^:]:%lf",buf,&f);  
      _energycal[buf] = f; 
      return;
    }
   
   char * ch; 
   char * strcp = strdup(value.c_str()); 
   ch = strtok(strcp,","); 
   while (ch!=NULL)
   {
     sscanf(ch,"%[^:]:%lf",buf,&f);  
     _energycal[buf] = f; 
     std::cout << buf << ": " << f << std::endl; 
     ch = strtok(NULL,","); 
   }

   free (strcp); 
   return;
  }

  if (key == "wormfilename") 
  {
    _wormfilename.clear(); 

    if (!strcasecmp(value.c_str(),"none") || !strcasecmp(value.c_str(),""))
    {
       return;   
    }

    int ntokens = 1; 

    for (unsigned int i = 0; i < value.length() ; i++)
    {
      if (value[i] == ',') ntokens++; 
    }

    char buf[80]; 
    char f[25]; 
    if (ntokens == 1)
    {
      sscanf(value.c_str(),"%[^:]:%s",buf,f);  
      _wormfilename[buf] = f; 
      return;
    }
   
   char * ch; 
   char * strcp = strdup(value.c_str()); 
   ch = strtok(strcp,","); 
   while (ch!=NULL)
   {
     sscanf(ch,"%[^:]:%s",buf,f);  
     _wormfilename[buf] = f; 
     std::cout << buf << ": " << f << std::endl; 
     ch = strtok(NULL,","); 
   }

   free (strcp); 
   return;
  }

  if (key == "wormcutval") 
  {
    _wormcutval.clear(); 

    if (!strcasecmp(value.c_str(),"none") || !strcasecmp(value.c_str(),""))
    {
       return;   
    }

    int ntokens = 1; 

    for (unsigned int i = 0; i < value.length() ; i++)
    {
      if (value[i] == ',') ntokens++; 
    }

    char buf[80]; 
    double f; 
    if (ntokens == 1)
    {
      sscanf(value.c_str(),"%[^:]:%lf",buf,&f);  
      _wormcutval[buf] = f; 
      return;
    }
   
   char * ch; 
   char * strcp = strdup(value.c_str()); 
   ch = strtok(strcp,","); 
   while (ch!=NULL)
   {
     sscanf(ch,"%[^:]:%lf",buf,&f);  
     _wormcutval[buf] = f; 
     std::cout << buf << ": " << f << std::endl; 
     ch = strtok(NULL,","); 
   }

   free (strcp); 
   return;
  }


  if (parseRegistered(key,value,n)) 
  {
    return; 
  }
      
//If we get here, that means there was an unrecognized key
  std::cerr << "Unrecognized key \""<< key <<"\" on line " 
            << n << ". Skipping." << std::endl;

}

