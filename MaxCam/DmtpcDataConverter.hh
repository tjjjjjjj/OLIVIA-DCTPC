#ifndef DMTPC_DATA_CONVERTER_HH
#define DMTPC_DATA_CONVERTER_HH

#include "DmtpcDataset.hh"
#include "DmtpcEvent.hh"
#include "TH1C.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2S.h"
#include "ScopeDataInfo.hh"
#include "ScopeWaveformData.hh"

/** This namespace contains functions to facilitate conversion between  
 *  old-style TH2F and TH1F containers to smaller TH2S and ScopeWaveformData (which is a type 
 *  of TH1C with additional information in it. 
 *
 *  \author Cosmin Deaconu 
 */

namespace DmtpcDataConverter
{

  /** Converts an old format (<5) DmtpcEvent to the new (>=5) format which uses 
   * smaller data types 
   *
   * @param old The DmtpcEvent to convert
   * @param placement_ptr if you want to overwrite the contents of an existing DmtpcEvent, pass it here, otherwise use NULL
   * */
  DmtpcEvent * convertHistDataTypes(DmtpcEvent * old, DmtpcEvent * placement_ptr = NULL); 

  /** Converts a data file from the old format (<5) to the new (>5) format which uses smaller data types. 
   *  @param infile The path of the file to convert
   *  @param outfile The path of the converted file 
   *  @param overwrite Allow writing over an existing file for the outfile 
   */
  void convertDataFile(const char * infile, const char * outfile, bool overwrite=true); 

  /** Convert a scope waveform from a TH1F to a TH1C to reduce size. The TH1C is stored
   * using the arbitrary digitizer units which are unsigned chars. However, root TH1C's store signed 
   * chars, so attempting to plot the raw data is not particular useful (since the zero point
   * is generaly 128). Use scopeExpand to convert the ScopeWaveformData objects back to TH1F's in physical units. 
   *
   * @param old Pointer to the waveform to convert
   * @param sinfo Pointer to the associated ScopeDataInfo 
   * @param new_name The name of the new waveform, or NULL to use the old name 
   * @param placement_ptr Use existing ScopeWaveformData memory, or NULL to allocate new memory 
   */
  ScopeWaveformData * scopeReduce(const TH1F * old, const ScopeDataInfo * sinfo, const char * new_name = NULL, TH1F * placement_ptr = NULL); 

  /** Convert a CCD tiomage from a TH2F to a TH2S to reduce size. Both will be in ccd units, although because TH2S's are
   * signed, the TH2S may have wrapping around to negative values (in practice this can only happen with sparks). Use 
   * ccdExpand to correctly convert the TH2S back to a TH2F without wrapover issue. 
   * @param old Pointer to image to convert 
   * @param new_name The name of the new image, or NULL to use the old name 
   * @param placement_ptr Use the memory of an existing TH2S, or NULL to allocate new memory
   */
  TH2S * ccdReduce(const TH2F * old, const char * new_name = NULL, TH2S * placement_ptr = NULL); 

  /** Converts a reduced ScopeWaveformData object to a TH1F object using physical units. As with legacy TH1F's,
   *  the timestamp is copied into the underflow bin of the TH1F. 
   *
   *  @param reduced The ScopeWaveformData to expand
   *  @param new_name The new name of the histogram or NULL to use the existing name
   *  @param placement_ptr Use the memory of an existing TH1F, or NULL to allocate new memory
   */
  TH1F * scopeExpand(const ScopeWaveformData * reduced, const char * new_name = NULL, TH1F * placement_ptr = NULL); 

  /** Converts a reduced TH2S * back to a TH2F * while doing the necessary casting to ensure values over 32767 are preserved properly
   *
   * @param reduced The TH2S* to expand  
   * @param new_name The new name of the image or NULL to use the same name 
   * @param placement_ptr Use the memory of an existing TH2F, or NULL to allocate new memory
   */
  TH2F * ccdExpand(const TH2S * reduced, const char * new_name = NULL, TH2F * placement_ptr = NULL); 

}
#endif 
