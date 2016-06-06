//calibration constants for Little DCTPC

//order of calibrations
//position of detector center (changes every set)
//energy calibration (changes every sequence)
//radial calibration (never changes...but be aware that detector center does change)

const int TOTAL_SET_NUMBER = 3;
const int TOTAL_SEQ_NUMBER = 18;


//set-by-set position calibration. Each number denotes the pixel corresponding to the center of the detector. When considering position in the analysis, the center of the detector (rather than the center of the image) shall be defined as (0,0).
double POSITION_X_CALIB[TOTAL_SET_NUMBER] = {0.,0.,0.};
double POSITION_Y_CALIB[TOTAL_SET_NUMBER] = {0.,0.,0.};
double POSITION_RADIUS_CALIB = -1.;//radius in pixels (radius is zzz cm)


//sequence-by-sequence energy calibration
//note that the waveform and CCD calibration constants can be a bit different because the charge/light ratio is dependent on gain
double ENERGY_CALIB_wf[TOTAL_SEQ_NUMBER] =   {1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00};
double ENERGY_CALIB_ccd[TOTAL_SEQ_NUMBER] =   {1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00};


//sequence by sequence energy calibration                             
double ENERGY_CALIB[TOTAL_SEQ_NUMBER]= {0.955114,1.06241,1.07105,1.07669,1.02827,1.15574,1.18828,1.08779,0.998239,1.05118,1.03128,1.02639,0.992746,0.925486,1.1664,1.05791,0.989169};


//radial calibration, according to the equation: [0] / ([1] + radius)^[2]
double RADIAL_PARAM0_CALIB = 132.449;
double RADIAL_PARAM1_CALIB = 147.474;
double RADIAL_PARAM2_CALIB = 0.99942;

//declaring variables used in macros
double r2;
double weight;
double rho_start;
double rho_end;
double verticallength;
double horizontallength;
double tracklength;
double lengthenergy;
double meshcalib;
double ccdcalib;
double anodecalib;
double lengthcalib;
double delta_anodeccd;
double delta_meshccd;
double delta_anodemesh;
double delta_meshlength;
double delta_anodelength;
double delta_ccdlength;
double mesh_over_ccd;
double ccdenergy;

//drift velocity
double driftspeed =3.1e06;//2.6!!2.5 worked 3.7e06; //in units of mm/sec//DRIFT_VELOCITY
double secondspersample=1;//4.e-0.9


//diffusion constant
double DIFFUSION_CONSTANT = -1.; 


//pixel size; the CCD is 1024x1024 pixels. We bin at 4x4 but a pixel is defined as 1 of the 1024 (rather than as 1 of 256)
double MMPERPIXEL = 166.6/1024.;


//detector orientation (angle) relative to ??
double DETECTOR_ORIENTATION = -1.;


//detector height
double DETECTOR_HEIGHT = -1.;


//digitizer sample rate
double DIGITIZER_SAMPLERATE= 250.; //250 MegaSamples/second


//The exposure in seconds for each sequence. This can be found with reduced_analysis/LittleDCTPC_physics/run_exposure.C                      
double EXPOSURE_SEQ[TOTAL_SEQ_NUMBER]={0., 345130., 87757., 248539., 367407., 364018., 89584., 251442., 324605., 358127., 359194., 357272., 365276., 354349., 370056., 125739., 103122., 1111687.};


//Anode, Mesh and CCD Calibration Constants
double ANODE_CALIB=0.91705;
double MESH_CALIB=0.620001;
double CCD_CALIB=1.31299;//option 2 -> 1.10665, option 3-> 0.6554649
double Eccd_0=-135.834;
double Eccd_1=0.760093;
double Eccd_2=-3.05039e-05;
double LENGTH_CALIB=1.04766;

//Voltage to Energy
double VOLT_TO_ENERGY=6227.9;

//Calibration Fit Parameters

//anode calibration gaus fit
double anode_gausfit_par0=152.401;
double anode_gausfit_par1=4589.25;
double anode_gausfit_par2=146.339;

//mesh calibration linear fit using calibrated anode
double anodemesh_linearfit_par0=25.1687;
double anodemesh_linearfit_par1=0.620001;

//ccd calibration linear fit using calibrated anode
double anodeccd_linearfit_par0=-88.4827;
double anodeccd_linearfit_par1=0.977151;

//ccd calibration polynomial fit using calibrated mesh
double meshccd_polyfit_par0=-59.7901;
double meshccd_polyfit_par1=1.1158;
double meshccd_polyfit_par2=-5.04164e-05;

//reconstructing energy from length
double El_0=-148.822;
double El_1=79.9262;
double El_2=-0.26157;

//efficiency fits

//passwfandccd/passccd length energy

double efflength0=0.273859;
double efflength1=-3.18967e-05;

//passwfandccd/passccd ccd energy

double effccd0=0.248027;
double effccd1=9.17716e-06;

