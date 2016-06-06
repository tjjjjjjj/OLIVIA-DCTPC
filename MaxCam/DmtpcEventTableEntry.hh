#ifndef DMTPC_EVENT_TABLE_ENTRY_HH
#define DMTPC_EVENT_TABLE_ENTRY_HH
#include "TString.h"
#include "TObject.h"
#include <iostream>
/** 
 * This class represents an entry in a DmtpcEventTable
 */
class DmtpcEventTableEntry : public TObject
{
  public:
    /** 
     * Creates a DmtpcEventTable with a given row and column title. 
     * \param row The name of the row
     * \param column The name of the column
     */
    DmtpcEventTableEntry(const char * row="", const char * column="");
    
    /**
     * Copy Constructor
     * \param other the DmtpcEventTableEntry to clone
     */
    DmtpcEventTableEntry(DmtpcEventTableEntry & other); 

    /** 
     * Destructor
     */
    ~DmtpcEventTableEntry(); 

    /** 
     * Increments the entry by one
     * @return the updated value
     */
    int increment() {return ++val;} 

    /** 
     * Adds an amount to the entry
     * @param amt the amount to add
     * @return the updated value
     */
    void add(int amt) {val+=amt;}

    /** 
     * Returns the value of the entry
     * @return the value of the entry
     */
    int value(){return val;} 

    /** 
     * Returns the hash number of the entry. This is
     * calculated using the concatenation of the row and column names
     * @return the hash numbe
     */
    ULong_t Hash() const{ return hashnum;} 
    /**
     * Checks if an entry is equal to another entry. 
     * @param other a TObject to compare against. Will crash if not of right type. 
     * @return a boolean indicating the equality
     */
    Bool_t IsEqual(const TObject * other) const; 

    /**
     * Gets the row name of the entry. 
     * @return row name as a TString
     */
    TString  rowName() const{ return row;}

    /**
     * Gets the column name of the entry. 
     * @return column name as a TString
     */
    TString  columnName() const{ return column;}
  private:

    TString row; ///<row name
    TString column;///< column name
    ULong_t hashnum; ///< hashnumber calculated from concatenated row and column
    int val; ///< value in this entry

    ClassDef(DmtpcEventTableEntry,1);
};

#endif
