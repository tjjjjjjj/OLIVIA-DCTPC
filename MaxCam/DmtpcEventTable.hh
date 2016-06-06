#ifndef DMTPC_EVENT_TABLE_HH
#define DMTPC_EVENT_TABLE_HH
#include "THashTable.h"
#include "TNamed.h"
#include "TClonesArray.h"
#include "DmtpcEventTableEntry.hh"
#include <ostream>
#include <iostream>

/** 
 *   A DmtpcEvent table provides a convenient way to store information 
 *   about the number of events in a certain category (e.g. passing 
 *   certain cuts). 
 *
 *  
 *   
 */
class DmtpcEventTable : public TNamed
{
  public:
    /**
     *  Copy Constructor
     *  @param other Another DmtpcEventTable to make a copy of
     */
    DmtpcEventTable(const DmtpcEventTable & other); 

    /** 
     * Create DmtpcEventTable with a given name or title
     * @param name the name to use for the table
     * @param title the title to use for the table
     */ 
    DmtpcEventTable(const char * name = "", const char * title = ""); 

    /** 
     * Create DmtpcEventTable with a given name or title
     * @param name the name to use for the table
     * @param title the title to use for the table
     */ 
    DmtpcEventTable(const TString & name , const TString & title);

    /** 
     * Class destructor
     */
    ~DmtpcEventTable();

    /** 
     * Merges the contents of another table with this table. 
     * @param other A pointer to another DmtpcEventTable
     */
    void combineWith(DmtpcEventTable * other); 

    /** 
     * Increments an entry. If it does not exist already it is created.
     * @sa add()
     * @param row the name of the row to increment
     * @param column the name of the column to increment
     * @return the new value of the entry
     */ 
    int increment(const char * row, const char * column){ return add(row,column,1);} 

    /** 
     * Adds a value to an entry. If it does not exist already it is created.
     * @sa add()
     * @param row the name of the row to increment
     * @param column the name of the column to increment
     * @param amt the amount to add (can be negative)
     * @return the new value of the entry
     */ 
    int add(const char * row, const char * column, int amt); 

    /** 
     * Returns the value of value to an entry. If it does not exist already it is created.
     * @param row the name of the row to check
     * @param column the name of the column to check
     * @return the value of the entry
     */ 
    int value(const char * row, const char * column); 

    /** 
     * Printout the table in fancy ascii format
     * @param out the name of the out stream to print to
     */
    void print(ostream & out = std::cout);
    /** 
     * Printout the table in Latex format
     * @param out the name of the out stream to print to
     * \warning Not yet implemented
     */
    void printLatex(ostream & out = std::cout); 

    /** 
     * Returns a TObjArray of the rows. 
     * members should be TObjStrings
     * @return TObjArray of rows
     */
    TObjArray * getRows() { return rows; }

    /** 
     * Returns a TObjArray of the columns
     * members should be TObjStrings
     * @return TObjArray of columns
     */
    TObjArray * getColumns() { return columns; }

  private:
    /** hash table that lets you find the entry from a row and column name */
    THashTable * hash; 

    /** row names **/
    TObjArray * rows; 

    /** column names */ 
    TObjArray * columns; 
   
    ClassDef(DmtpcEventTable,1); 
}; 
#endif

