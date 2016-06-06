#include "SimpleConfig.hh"
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <strings.h>
#include <cstring>



int SimpleConfig::writeOutConfig(char * file)
{
  std::ofstream f; 
  f.open(file); 
  if (!f.is_open())
  {
    std::cerr << "Could not open file " << file << " for some reason." << std::endl;
    return 1; 
  }
  print(f); 
  f.close(); 
  return 0; 
}

static std::string removeWhitespace(std::string in)
{
  std::stringstream out; 

  for (std::string::size_type i = 0; i < in.size(); i++)
  {
    if (in[i]!=' ' && in[i]!='\t')
        out << in[i]; 
  }

  return out.str(); 

}


void SimpleConfig::parseFile(const char * file)
{

  std::ifstream f(file);  
  if (f.is_open()) parseStream(f); 
  else std::cerr << "Could not open file " << file << std::endl; 
}

void SimpleConfig::parseStream(std::istream & in)
{
  std::string line; 
  int n = 0; 
  while(! in.eof())
  {
    std::getline(in,line); 
    n++;

    //Get rid of all white space
    line = removeWhitespace(line); 

    if (line.size() == 0) continue; 
    
    //Get rid of comments
    std::string::size_type comment_loc = line.find("#"); 
    if (comment_loc != std::string::npos)
    {
      line.replace(comment_loc,line.size() - comment_loc, ""); 
    }

    if (line.size() ==0) continue;  
    
    std::string::size_type equal_sign_loc = line.find("="); 
    
    //Make sure there is an equal sign
    if (equal_sign_loc == std::string::npos )
    {
      std::cerr << "No equal sign found in non-empty line " << n
                << ". Skipping." << std::endl;
      continue; 
    }
    //Make sure only one equal sign
    if (equal_sign_loc != line.rfind("="))
    {
      std::cerr << "More than one = sign found in config file line " << n
                <<". Skipping." << std::endl;
      continue; 
    }
    
    //Make sure there is something both before and after the equal sign
    if (equal_sign_loc == 0 || equal_sign_loc == line.size() - 1)
    {
      std::cerr << "No key or no value on line " << n << ". Skipping." 
                << std:: endl;

      continue; 
    }
    
    if (CFG_DBG)
    {
      std::cout << "line: \"" << line << "\"" << std::endl;
    }

    parsePair(line.substr(0,equal_sign_loc), line.substr(equal_sign_loc+1),n);

  }
   
}

static const int TYPE_DBL = 0; 
static const int TYPE_INT = 1; 
static const int TYPE_UINT = 2; 
static const int TYPE_STR = 3; 
static const int TYPE_ENUM = 4; 
static const int TYPE_BOOL = 5; 
static const int TYPE_INVALID=-1; 
void SimpleConfig::registerNewDouble(char * name, double * ptr)
{
  SimpleConfigStore_t store; 
  store.name = name; 
  store.type = TYPE_DBL; 
  store.ptr = (void *) ptr;
  registered.push_back(store);  
}

void SimpleConfig::registerNewUInt(char * name, unsigned int * ptr)
{
  SimpleConfigStore_t store; 
  store.name = name; 
  store.type = TYPE_INT; 
  store.ptr = (void *) ptr; 
  registered.push_back(store); 
}

void SimpleConfig::registerNewInt(char * name, int * ptr)
{
  SimpleConfigStore_t store; 
  store.name = name; 
  store.type = TYPE_INT; 
  store.ptr = (void *) ptr; 
  registered.push_back(store); 
}

void SimpleConfig::registerNewBool(char * name, bool * ptr)
{
  SimpleConfigStore_t store; 
  store.name = name; 
  store.type = TYPE_BOOL; 
  store.ptr = (void *) ptr; 
  registered.push_back(store); 
}

void SimpleConfig::registerNewString(char * name, char ** ptr)
{
  SimpleConfigStore_t store; 
  store.name = name; 
  store.type = TYPE_STR; 
  store.ptr = ptr; 
  registered.push_back(store); 
}

void SimpleConfig::registerNewEnum(char * name, int n, char ** keys, int * indices, void * ptr)
{
  SimpleConfigStore_t store; 
  store.name = name; 
  store.type = TYPE_ENUM; 
  store.ptr = ptr; 
  store.n = n; 
  store.keys = keys; 
  store.indices = indices; 
  registered.push_back(store); 
}

void SimpleConfig::parsePair(std::string key, std::string value, int n)
{
  parseRegistered(key,value,n); 
}

bool SimpleConfig::parseRegistered(std::string key, std::string value, int n)
{
  
  
  for (std::string::size_type i = 0; i < registered.size(); i++)
  {
      if (strcasecmp(key.c_str(), registered[i].name) == 0)
      {
          switch(registered[i].type)
          {
            case TYPE_DBL:
              *((double*)registered[i].ptr) = atof(value.c_str()); break; 
            
            case TYPE_INT:
              *((int*) registered[i].ptr) = atoi(value.c_str()); break; 

            case TYPE_BOOL:
              *((bool*) registered[i].ptr) = strcasecmp(value.c_str(),"true") == 0 || atoi(value.c_str()) !=0; break; 

            case TYPE_UINT:
              *((unsigned int*) registered[i].ptr) = atoi(value.c_str()); break; 
            
            case TYPE_STR:
              *((char**) registered[i].ptr) = strdup(value.c_str()); break; 
            
            case TYPE_ENUM:
              for (int j = 0; j < registered[i].n; j++)
              {
                if ( strcasecmp(value.c_str(), registered[i].keys[j]) ==0)
                {
                  if (registered[i].indices==NULL)
                    *((int*)registered[i].ptr) = j; 
                  else 
                    *((int*) registered[i].ptr) = registered[i].indices[j]; 
                  break; 
                }
                if (j == registered[i].n - 1)
                {
                  std::cerr << "Warning, don't understand value \""<<value<<"\" on line " << n << ". Skipping." << std::endl; 
                  registered[i].type = TYPE_INVALID; 
                }
              }
          }
          return true; 
      }
  }

  return false; 
}

void SimpleConfig::print(std::ostream & out)
{
  printRegistered(out); 
}

void SimpleConfig::printRegistered(std::ostream & out) const
{
  for (std::string::size_type i = 0; i < registered.size(); i++)
  {
    char * printme; 
    int val; 
    switch(registered[i].type)
    {
      case TYPE_DBL: 
        out << registered[i].name << " = " <<*( (double*) registered[i].ptr) << std::endl; 
        break; 
      case TYPE_INT: 
        out << registered[i].name << " = " <<*( (int*) registered[i].ptr) << std::endl; 
        break; 
      case TYPE_BOOL: 
        printme = *( (bool*) registered[i].ptr) ? strdup("true") : strdup("false");
        out << registered[i].name << " = " << printme << std::endl; 
        break; 
      case TYPE_UINT: 
        out << registered[i].name << " = " <<*( (unsigned int*) registered[i].ptr) << std::endl; 
        break; 
      case TYPE_STR:
        out << registered[i].name << " = ";
        printme =  *((char**) registered[i].ptr); 
        if (printme == NULL) printme= ""; 
        out << printme << std::endl; 
        break; 
      case TYPE_ENUM:
        out << registered[i].name << " = "; 

        val = *( (int*) registered[i].ptr); 
        if (registered[i].indices == NULL) 
        {
          out << registered[i].keys[val] << std::endl; 
        }
        else 
        {
          for (int j = 0; j < registered[i].n ; j++)
          {
            if (registered[i].indices[j] == val)
            {
              out << registered[i].keys[j] << std::endl; 
              break; 
            }
          } 
        }

      case TYPE_INVALID:
        break; 

    }
  }
  
}


