#include "DmtpcStringTools.hh"
 
using namespace std;

void DmtpcStringTools::ltrim(string& str) {
  string::size_type pos = 0;
  while (pos < str.size() && isspace(str[pos])) pos++;
  str.erase(0, pos);
}
 
void DmtpcStringTools:: rtrim(string& str) {
  string::size_type pos = str.size();
  while (pos > 0 && isspace(str[pos - 1])) pos--;
  str.erase(pos);
}
 
void DmtpcStringTools:: trim(string& str) {
  ltrim(str);
  rtrim(str);
}

string DmtpcStringTools:: basename(string str)
{
  int pos = str.rfind("/"); 

  if (pos == string::npos) 
    pos = 0; 
  else
    pos++; 

  return str.substr(pos,string::npos); 
}
