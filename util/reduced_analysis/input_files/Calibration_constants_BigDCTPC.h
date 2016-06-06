//calibration constants for Big DCTPC

//order of calibrations
//position of detector center (changes every set)
//energy calibration (changes every sequence)
//radial calibration (never changes...but be aware that detector center does change)

const int TOTAL_SET_NUMBER = 4;
const int TOTAL_SEQ_NUMBER = 17;

//set-by-set position calibration. Each number denotes the pixel corresponding to the center of the detector. 
double POSITION_X_CALIB[TOTAL_SET_NUMBER] = {9.3,5.8,5.7,4.0};
double POSITION_Y_CALIB[TOTAL_SET_NUMBER] = {14.3,22.2,25.3,20.9};
double POSITION_RADIUS_CALIB[TOTAL_SET_NUMBER] = {580.1,582.1,582.9,581.9};	//radius in pixels (radius is 20.0 cm)
double POSITION_ANGLE[TOTAL_SET_NUMBER] = {-33.78,-33.23,-33.02,-33.19}; 	//degrees, slope of spacers in xy plane
double POSITION_IMAGING_AREA[TOTAL_SET_NUMBER] = {1142.26,1136.61,1134.13,1137.44}; 	//cm^2, used fit radius
double POSITION_IMAGING_AREA_ERR[TOTAL_SET_NUMBER] = {6.46,6.65,6.40,6.58};		//cm^2

//sequence-by-sequence energy calibration
//note that the waveform and CCD calibration constants can be a bit different because the charge/light ratio is dependent on gain
double ENERGY_CALIB_wf[TOTAL_SEQ_NUMBER] =   {1.,1.,1.,1.,1.,1.,1.,   1.01,1.05,0.83,0.72,1.38,0.65,0.56,0.56,0.78,.65};
double ENERGY_CALIB_ccd[TOTAL_SEQ_NUMBER] =   {1.,1.,1.,1.,1.,1.,1.,   1.01,1.05,0.80,0.60,1.97,0.65,0.56,0.56,0.78,.65};


//radial calibration, according to the equation: [0] / ([1] + radius)^[2]
double RADIAL_PARAM0_CALIB = 7.56874;
double RADIAL_PARAM1_CALIB = 146.405;
double RADIAL_PARAM2_CALIB = 0.388667;


//drift velocity
double DRIFT_VELOCITY = -1.; //in units of mm/microsec


//diffusion constant
double DIFFUSION_CONSTANT = -1.; 


//pixel size; the CCD is 1024x1024 pixels. We bin at 4x4 but a pixel is defined as 1 of the 1024 (rather than as 1 of 256)
double MMPERPIXEL = 0.3447;
double MMPERPIXEL_ERR = 0.0007;


//detector orientation (angle) relative to ??
double DETECTOR_ORIENTATION = -1.;


//detector height
double DETECTOR_HEIGHT = -1.;


//digitizer sample rate
double DIGITIZER_SAMPLERATE= 250.; //250 MegaSamples/second


//The exposure in seconds for each sequence. This can be found with reduced_analysis/BigDCTPC_physics/run_exposure.C
double EXPOSURE_SEQ[TOTAL_SEQ_NUMBER] = {0.,0.,0.,0.,0.,0.,0.,519160.,438419.,331192.,295480.,1567903.,1286022.,299446.,654818.,380117.,1345758.};

