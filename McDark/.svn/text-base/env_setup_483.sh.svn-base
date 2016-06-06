# Timur Sahin 08 August 2007
# Script to set up environment variables needed for
# running Geant4. Should be executed as inclusion by
# entering the following into the command line:
# source env_setup_483.sh
# Running with sh will NOT set up search path directory
# correctly!


# NOTE: Please change the G4WORKDIR and the ROOTOUT locations 
# to where you want your executable and output to go!


PATH=$PATH:/net/cockroach/data01/tcsahin/wired/bin/
echo Wired has been set up to the network path. See comments for
echo instructions on installing Wired locally.
# Though it is possible  to run wired off of the network, it will run
# relatively slowly. It is recommended you install Wired onto your own
# machine and comment out these above lines of code. Instructions on
# installing wired (as well as necessary setup files) can be found in
# /net/cockroach/data01/tcsahin/tars/wired.  
echo out these lines in the setup script. For instructions on setting 
PATH=$PATH:/net/cockroach/data01/tcsahin/dawn/dawn_3_88a/
PATH=$PATH:/net/cockroach/data01/tcsahin/dawn/gsview/
echo Search path has been set up.

# Please change this variable to where you want things to go!!
export G4WORKDIR=~/g4work
echo Work directory set to $G4WORKDIR.

if [ ! -d $G4WORKDIR/out ]; then 
  mkdir $G4WORKDIR/out
fi

# Please change this as well if you want the output in a different directory
export ROOTOUT=$G4WORKDIR/out
export G4HEPREPFILE_DIR=$ROOTOUT
echo ROOT output and HepRepFile output currently set to $ROOTOUT




export LD_LIBRARY_PATH=/net/cockroach/data01/tcsahin/CLHEP/lib:/export/data01/root/lib
source /net/cockroach/data01/tcsahin/geant4.8.3/env.sh
echo YOU ARE USING VERSION 4.8.3 OF GEANT!
echo Environment setup complete.
