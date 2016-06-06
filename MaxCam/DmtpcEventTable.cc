#include "DmtpcEventTable.hh"
#include "TObjString.h"
#include <iomanip>
ClassImp(DmtpcEventTableEntry); 
DmtpcEventTableEntry::DmtpcEventTableEntry(const char * r,const char * c)
{
  row =  TString(r);
  column = TString(c);
  TString hashstring = row +"_"+column; 
  hashnum = hashstring.Hash(); 
  val = 0; 
}

DmtpcEventTableEntry::DmtpcEventTableEntry(DmtpcEventTableEntry & other)
{
  row = TString(other.row);
  column = TString(other.column);
  val = other.val; 
  hashnum = other.hashnum; 
}

Bool_t DmtpcEventTableEntry::IsEqual(const TObject * other) const
{
  DmtpcEventTableEntry * o = (DmtpcEventTableEntry *) other; 
  return  (Hash() == o->Hash() && rowName() == o->rowName()
           && columnName() == o->columnName());
}

DmtpcEventTableEntry::~DmtpcEventTableEntry()
{
  //Nothing to see here, please move along
}

ClassImp(DmtpcEventTable); 




DmtpcEventTable::DmtpcEventTable(const char * name, const char * title)
{
  fName = TString(name);
  fTitle = TString(title); 
  hash= new THashTable(); 
  rows = new TObjArray();
  columns = new TObjArray() ;
}

DmtpcEventTable::DmtpcEventTable(const TString & name, const TString & title)
{
  fName = name;
  fTitle = title; 
  hash= new THashTable(); 
  rows = new TObjArray();
  columns = new TObjArray() ;
}

DmtpcEventTable::DmtpcEventTable(const DmtpcEventTable & other)
{

  hash = (THashTable *) other.hash->Clone();
  rows = (TObjArray *) other.rows->Clone();
  columns = (TObjArray *) other.columns->Clone();
}

DmtpcEventTable::~DmtpcEventTable()
{
  hash->Delete();
  rows->Delete();
  columns->Delete(); 
}

int DmtpcEventTable::add(const char * row, const char * column, int amt)
{
  TObjString *r = new TObjString(row);
  TObjString *c = new TObjString(column);   
  
  if (!rows->Contains(r))
  {
    rows->Add(r); 
  }
  else
  {
    delete r; 
  }

  if (!columns->Contains(c))
  {
    columns->Add(c); 
  }
  else
  {
    delete c; 
  }

  DmtpcEventTableEntry entry(row,column);
  DmtpcEventTableEntry * e = (DmtpcEventTableEntry *) hash->FindObject(&entry);

  if (e!=NULL) {
    e->add(amt); 
    return e->value(); 
  }
  else
  {
    DmtpcEventTableEntry * copy = new DmtpcEventTableEntry(entry);
    copy->add(amt); 
    hash->Add(copy);  
    return amt; 
  }

  return 0; 
}

void DmtpcEventTable::combineWith(DmtpcEventTable * o)
{
  //Merge Column Names
  for (int i = 0; i < o->getColumns()->GetEntries(); i++)
  {
    TObject * c = o->getColumns()->At(i);
    if (!columns->Contains(c))
      columns->Add(c); 
  } 
  
  //Merge Row Names
  for (int i = 0; i < o->getRows()->GetEntries(); i++)
  {
    TObject * r = o->getRows()->At(i);
    if (!rows->Contains(r))
      rows->Add(r); 
  } 
  

  //Add all values
  for (int i = 0; i < rows->GetEntries(); i++)
  {
    for (int j = 0; j < columns->GetEntries(); j++)
    {
      const char * r = ((TObjString*)rows->At(i))->String().Data(); 
      const char * c = ((TObjString*)columns->At(j))->String().Data(); 
      add(r,c,o->value(r,c)); 
    }
  }

}

int DmtpcEventTable::value(const char * row, const char * column)
{
  DmtpcEventTableEntry entry(row,column); 

  if (hash->Contains(&entry))
  {
    return ((DmtpcEventTableEntry *) hash->FindObject(&entry))->value(); 
  }
 
  else return 0; 
}

void DmtpcEventTable::print(ostream & out)
{

 //Print out column Names:
 out << "            | " << std::left;  
 for (int i = 0; i < columns->GetEntries(); i++)
 {
   out << std::setw(12) << ((TObjString *)columns->At(i))->GetString() << "| "; 
 }

 out << "\n";
 for (int i = 0; i <=columns->GetEntries(); i++)
    out << "---------------";
 out << "\n";
 for (int i = 0; i < rows->GetEntries(); i++)
 {
    out << std::setw(12) <<((TObjString*)rows->At(i))->GetString() << "| "; 
    
    for (int j = 0; j < columns->GetEntries(); j++)
    {
       out << std::setw(12); 
       out << value( ((TObjString *)rows->At(i))->GetString().Data(),
                     ((TObjString *)columns->At(j))->GetString().Data()
                     ); 

       out << "| ";
    }
    out << "\n";
 } 
}

void DmtpcEventTable::printLatex(ostream & out)
{
  out << "Not yet implemented\n";
}

