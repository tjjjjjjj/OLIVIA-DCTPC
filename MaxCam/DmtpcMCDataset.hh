#ifndef DMTPC_MC_DATASET_HH
#define DMTPC_MC_DATASET_HH

#include "TROOT.h"
#include "TFile.h"
#include "TObject.h"
#include "MaxCamClusterImage.hh"
#include "TH2F.h"
#include "TTree.h"
#include <vector>
#include <list>

class DmtpcMCDataset : public TObject
{

  public:
    DmtpcMCDataset(); 
    virtual ~DmtpcMCDataset(); 

    /** Load a MC Data file. 
     * \param filename the name of the file to open
     * \param load_true_clusters attempt to load the true cluster information (must be enabled when generating MC) 
     */
    bool loadFile(const char * filename, bool load_true_clusters = false);

    /** Get the specified event */ 
    void getEvent(int ev) { if (tree) tree->GetEntry(ev); } 

    /** Returns the number of cameras */ 
    int ncamera() const { return ncam;} 
 
   /** Returns the original cluster (if you didn't load true clusters, this probably will crash **/
    const MaxCamClusterImage * getCluster(int cam = 0) const { return  clusters ? (MaxCamClusterImage*) clusters->At(cam) : 0 ; } 

    /** Get the MC track phi */
    float getPhi() const {return phi;} 

    /** Get the MC track starting  x
     *
     *  @param physical return in mm world units rather than pixel units
     * */
    float getX(bool physical=false) const {return physical ? x : (x-campos[0]) / lengthcal + pix[0]/2.;} 

    /** Get the MC track starting  y
     *
     *  @param physical return in mm world units rather than pixel units
     * */
    float getY(bool physical=false) const {return physical ? y : (y-campos[1]) / lengthcal + pix[1]/2.;} 

   /** Get the MC track starting  z
     *
     *  @return z in mm
     * */
    float getZ() const {return  z;} 

    int getSequence() const {return sequence;}

    /** Get the MC track theta */
    float getTheta() const {return theta;} 

    /** Get the MC deposited energy in ccd units */
    float getE(int cam = 0) const {return integ[cam];}  

    /** Get the MC Particle energy in kev */
    float getParticleE(bool scint = false) const {return scint ? Escint :  E; }  

    /** Get the MC range in ccd units */ 
    float getRange() const { return length/lengthcal; } 

    float getRangeZ() const { return zlength; }

    /** Get lengthcal **/ 
    float getLengthcal() const { return lengthcal; } 


    /** Get pressure **/ 
    float getPressure() const { return pressure; }   

    /** Get gain **/ 
    float getGain() const { return gain; }   
    
    /** Get noise **/ 
    float getNoise() const { return noise; }   

    /** Get bias **/ 
    float getBias() const { return bias; } 
  
  private:
    float phi;
    float theta; 
    float lengthcal; 
    float length;
    float zlength;
    float E; 
    float Escint; 
    int sequence;
    float x; 
    float y; 
    float pressure; 
    float noise; 
    float bias;
    float gain; 
    float z; 
    float * integ;
    float campos[2];
    int ncam;
    int pix[2]; 
    float width[2]; 
    TTree * tree; 
    TFile * f; 
    TObjArray * clusters; 

  ClassDef(DmtpcMCDataset,0)

};

#endif
