#ifndef DMTPC_CAMERAMAP_HH
#define DMTPC_CAMERAMAP_HH

#include <vector>
#include "TString.h"
#include "TObject.h"

using std::vector;

/** 
 * Class to specify the relative rotations between cameras.  Useful
 * for the 4shooter where all cameras point in the same direction.
 * Initially generated to help visualize 4shooter camera data without
 * requiring image stitching.  Uses information about camera serial
 * number so that the scheme is robust against power-cycling the
 * cameras (which changes the "camera number" in the DmtpcDataset).
 *
 */
class DmtpcCameraMap : public TObject {

public:
    
  //
  // Ctors
  DmtpcCameraMap();
  
  virtual ~DmtpcCameraMap();
  
  virtual const char* GetName() const { return "DmtpcCameraMap"; }  

  /**
   *
   * @param[in] filename Name of file that contains information about the 
   *                     desired camera map.  File is ASCII and should
   *                     have 3 columns (serialNumber, pad, rotationType)
   *   For example:
   *   <pre>
   *   # serialNumber pad rotatePerfect
   *   081013   1   HALF_TURN
   *   100534   2   QUARTER_LEFT
   *   100439   3   QUARTER_RIGHT
   *   081264   4   NONE
   *   </pre>
   *
   * For information about rotationType, see MaxCamImageTools::rotatePerfect()
   *
   */
  void loadMap(TString filename);
  
  /**
   * Add a new camera to the list
   *
   * This is called internally by loadMap().  The user should not 
   * call this directly
   *
   * @param[in] serialNumber CCD serial number (see MaxCamConfig.serialNumber)
   * @param[in] pad Display this camera on this pad divided ROOT TCanvas 
   *                (see image below)
   * @param[in] rotationType How to rotate the image before displaying it
   *              for value options, see MaxCamImageTools::rotatePerfect()
   * 
   * Pad Numbering details:
   *       For the 4shooter, this will be either 1, 2, 3 or 4.
   *       the numbering is like so:
   * 
   * <pre>
   *           ----------------------
   *           |          |         |
   *           |          |         |
   *           |    1     |    2    |
   *           |          |         |
   *           ----------------------
   *           |          |         |
   *           |          |         |
   *           |    3     |    4    |
   *           |          |         |
   *           ----------------------
   * </pre>
   */
  void addCamera(TString serialNumber, Int_t pad, TString rotationType);

  Int_t getPad(TString sn); // return the pad # for a given serial number
  TString getRot(TString sn); // return the required rotation for a given sn
  double getRotRadians(TString sn); // return the required rotation for a given sn in radians
  Int_t getIdOfSN(TString sn);

  void  setVerbose(Int_t v) { _verbose = v;}
  Int_t getVerbose() { return _verbose; }
  
private:  
  
  vector< TString > _serialNumbers;
  vector< Int_t > _pads;
  vector< TString > _rotations;
  Int_t _verbose;

  ClassDef(DmtpcCameraMap,1)
        
};

#endif
