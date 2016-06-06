#ifndef DMTPC_STRING_TOOLS_HH
#define DMTPC_STRING_TOOLS_HH

#include <string>

namespace DmtpcStringTools {
  
  /** 
   * Strip leading whitespace from a string.
   *
   * @param[in] str Reference to the string to be trimmed. Modified in-place.
   * @return        void
   */
  void ltrim(std::string& str);
  
  /**
   * Strip trailing whitespace from a string.
   *
   * @param[in] str Reference to the string to be trimmed. Modified in-place.
   * @return        void
   */
  void rtrim(std::string& str);
  
  /**
   * Strip both leading and trailing whitespace from a string.
   *
   * @param[in] str Reference to the string to be trimmed. Modified in-place.
   * @return        void
   */
  void trim(std::string& str);
  

  /** 
   * Takes the part of a filepath after the last /
   *
   * @param str Reference to the string to get basename of
   * @return a new string corresponding to the basename
   */

   std::string basename(std::string str); 
};

#endif // DMTPC_STRING_TOOLS_HH
