#!/bin/bash
# Set system dependent variables
export DCTPC_TOP_DIR=/cvmfs/cvmfs.cmsaf.mit.edu/t2srv0008/export/app/d-Chooz/spitzj/DCTPC_soft
export DCTPC_ANA_DIR=$DCTPC_TOP_DIR/AnalysisFramework/v4/
#Run log
export BIGDCTPC_RUNLOG=/net/hisrv0001/home/spitzj/runlog
#Run sequence log. Run sequences are separated by gas refills. 
export DCTPC_SEQLOG=/net/hisrv0001/home/spitzj/DCTPC_simplified_log

case `uname -n` in
    (*ccage*)
	if [[ -z ${ROOTSYS} ]]; then
	    export ROOTSYS=/usr/local/root/v5.28.00f
	fi
	export DCTPC_DATA_DIR=/sps/dchooz/spitzj/DCTPC
    ;;
    (*cmsaf*)
	if [[ -z ${ROOTSYS} ]]; then
	    #export ROOTSYS=/app/d-Chooz/Software/root/root_v5.28.00
            export ROOTSYS=/cvmfs/cvmfs.cmsaf.mit.edu/t2srv0008/export/app/d-Chooz/Software/root/root_v5.28.00
            #export ROOTSYS=/app/d-Chooz/Software/root/root_v5.30.06
	fi
	export DCTPC_DATA_DIR=/net/nudsk0001/d00/scratch/dctpc_tmp/
    ;;
    (*dcnode.in2p3.fr*)
	if [[ -z ${ROOTSYS} ]]; then
	    export ROOTSYS=/home/dconline/ROOT/myROOT
	fi
	export DCTPC_DATA_DIR=/popdata/TPC_Test/
    ;;
    (dcpop.in2p3.fr)
	if [[ -z ${ROOTSYS} ]]; then
	    export ROOTSYS=/home/dconline/ROOT/myROOT
	fi
	export DCTPC_DATA_DIR=/popdata/TPC_Test/
    ;;
    (dcpip.in2p3.fr)
	if [[ -z ${ROOTSYS} ]]; then
	    export ROOTSYS=/home/dconline/ROOT/myROOT
	fi
	export DCTPC_DATA_DIR=/popdata/TPC_Test/
    ;;
    (dcftpc1)
        export DCTPC_DATA_DIR=/home/dctpc/
	export PYTHONPATH=/home/DCTPC/util/src:$PYTHONPATH
	export RECEIVER_ADDRESS="kazuhiro@mit.edu,spitzj@mit.edu,dctpc.farlab@gmail.com"
	export PRESSURE_THRES=2.0           # pressure check threshold ... deviation in percent
	export CCD_TEMPERATURE_THRES=-15.0  # ccd temperature check threshold ... degree celcius
	export CPU_TEMPERATURE_THRES=60.0   # cpu temperature check threshold ... degree celcius
	;;
    (*)
       #echo system not recognized!!! 
    ;;
esac

if [[ -z ${ROOTSYS} ]]; then
    echo Aborting...
else
    
    export ROOTSYS=/cvmfs/cvmfs.cmsaf.mit.edu/t2srv0008/export/app/d-Chooz/Software/root/root_v5.30.06
#export ROOTSYS=/app/d-Chooz/Software/root/root_v5.28.00
    export MCTABLES=$DCTPC_TOP_DIR/MaxCam/tables
    export LD_LIBRARY_PATH=$DCTPC_TOP_DIR/MaxCam:$DCTPC_TOP_DIR/cfitsio:$DCTPC_TOP_DIR/MaxCam/waveformtools/lib:$ROOTSYS/lib:/cvmfs/cvmfs.cmsaf.mit.edu/t2srv0008/export/app/d-Chooz/Software/lib:/cvmfs/cvmfs.cmsaf.mit.edu/t2srv0008/export/app/d-Chooz/Software/fftw-3.3.3/lib:/cvmfs/cvmfs.cmsaf.mit.edu/t2srv0008/export/app/d-Chooz/Software/gsl-1.15/lib:$DCTPC_TOP_DIR/cfitsio/lib:$DCTPC_TOP_DIR/MaxCam/Simulations/v1/lib:$LD_LIBRARY_PATH
    export PATH=$ROOTSYS/bin:$PATH
    export PYTHONPATH=$ROOTSYS/lib:$DCTPC_TOP_DIR/util/src:$DCTPC_ANA_DIR/util/src
    echo Configured for `uname -n` with ROOTSYS=${ROOTSYS}
fi

