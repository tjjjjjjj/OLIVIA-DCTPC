#ifndef SIMPLE_CONFIG_HH
#define SIMPLE_CONFIG_HH

#include <string>
#include <ostream>
#include <iostream>
#include <vector>

#define CFG_DBG 1 

/** This class is used to represent a field in memory */
class SimpleConfigStore_t 
{
  public: 
   char * name;  ///<Field name
   char type;    ///<Field type

   void * ptr;   ///<Pointer to data 
   int n; ///< If an enum, the number of choices
   char **  keys; ///< If an enum, a string array of the keys
   int * indices;  ///< If an enum, a string array of the indices
};


/** This class is used to store and parse config files of a config files of the following format. 
 *  To use this, you should inherit from this class.
 *  
   <pre>
      #comment 
      key = value #comment
      key = value
   </pre>

   
   <p>For common types, you may use the registerNewXXX protected methods in the subclass. This will allow you to parse
   any strings, numbers and enum types. 

   <p>For example, suppose you have a config with a value called foo;
   In your subclass you would have something like:

   <pre>
    class FooConfig : public SimpleConfig
    {
      public:
        FooConfig(); 
      private:
        double foo; 
    };

    ...

    FooConfig::FooConfig()
    {
      registerNewDouble("foo",&foo); 
    }
   </pre>

   If you want something more complicated, you may override the parsePair and print( ) methods to parse 
   a particular field. 

   See MaxCamChannelConfig for a good example of how to use this 
*/


class SimpleConfig
{
  public:
    /** Parses the given file and replaces any matching
       existing values with values from file 
      
      @param file filename of the file 
    */
    void parseFile(const char * file);  

    /** Parse a stream 
     * @param in the stream to parse
     * */ 
    void parseStream(std::istream & in); 

    /** Write out the current config to the named file
     @param file the filename of the file
     */
    int writeOutConfig(char * file); 

    /** Prints out the current config to the given output stream 
     *
     * @param out the output stream to print to (default is stdout)
     * */
    void print(std::ostream & out = std::cout) ; 

    /** Destructor. Does nothing */
    virtual ~SimpleConfig(){;}

  protected:
    /** This method is attempts to parse a key value pair from a file 
     * using the registered fields. 
     *
     * @param key the key to parse (the thing on the left hand side of the = )
     * @param value the key to value (the thing on the right hand side of the = ) 
     * @param line the line number of the line to be parsed.
     * @return true if the key was parsed, false otherwise
     *
     */
    bool parseRegistered(std::string key, std::string value, int line); 

    /** prints out all registered fields to the given output stream 
     * in the config format
     * @param out the output stream to print to (default is stdout)
     */
    void printRegistered(std::ostream & out = std::cout) const; 

    /** register a new double to parse
     *  @param name the name of the field to look for
     *  @param ptr the pointer to store the parsed value
     */
    void registerNewDouble(char * name, double * ptr); 

    /** register a new int to parse
     *  @param name the name of the field to look for
     *  @param ptr the pointer to store the parsed value
     */
    void registerNewInt(char * name, int * ptr); 

    /** register a new unsigned int to parse
     *  @param name the name of the field to look for
     *  @param ptr the pointer to store the parsed value
     */
    void registerNewUInt(char * name, unsigned int * ptr); 

   /** register a new string to parse. You need to prepare
     *  @param name the name of the field to look for
     *  @param ptr the pointer to store the parsed value
     */
    void registerNewString(char * name, char ** ptr); 

    /** register a new enum to parse. You 
     * must prepare a char * array of the valid keys
     * and an int array of the indices if they are
     * not ordinal. 
     *
     *  @param name the name of the field to look for
     *  @param n the number of choices in the enum
     *  @param keys char * array of valid keys
     *  @param indices array of the indices of the enum or NULL to use normal ordinal order (0,1,2,3...)
     *  @param ptr pointer to store the parsed value
     */
    void registerNewEnum( char * name, int n, char ** keys, int * indices, void * ptr); 

    /** register a new bool to parse
     *  @param name the name of the field to look for
     *  @param ptr the pointer to store the parsed value
     */
    void registerNewBool(char * name, bool * ptr); 

    std::string asString(); 
    
  private:
    /** This method attempts to parse a key value pair from a file. 
     * By default it just calls parseRegistered to attempt to parse a registered
     * field but can be overriden to do something more complicated
     *
     * @param key the key to parse (the thing on the left hand side of the = )
     * @param value the key to value (the thing on the right hand side of the = ) 
     * @param line the line number of the line to be parsed.
     *
     */
    virtual void parsePair(std::string key, std::string value, int line) ; 
    
    /** Stores the registered fields */
    std::vector<SimpleConfigStore_t> registered; 
};

#endif
