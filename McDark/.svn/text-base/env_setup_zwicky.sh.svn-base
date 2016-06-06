# Timur Sahin 08 August 2007
#  Updated by Cosmin Deaconu on 07 January 2011

# Script to set up environment variables needed for
# running Geant4 on zwicky. Should be executed as inclusion by
# entering the following into the command line:
# source env_setup_zwicky.sh 


# NOTE: Please change the G4WORKDIR and the ROOTOUT locations 
# to where you want your executable and output to go!

PATH=$PATH:/net/zwicky/dmtpc/software/dawn/dawn_3_88a/
PATH=$PATH:/net/zwicky/dmtpc/software/dawn/gsview/
echo Search path has been set up.

# Please change this variable to where you want things to go!!
export G4WORKDIR=/net/zwicky/dmtpc/`whoami`/g4work
echo Work directory set to $G4WORKDIR.

if [ ! -d $G4WORKDIR/out ]; then 
  mkdir $G4WORKDIR/out
fi

# Please change this as well if you want the output in a different directory
export ROOTOUT=$G4WORKDIR/out
export G4HEPREPFILE_DIR=$ROOTOUT
echo ROOT output and HepRepFile output currently set to $ROOTOUT




export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/net/zwicky/dmtpc/software/clhep/lib:/export/data01/root/lib
source /net/zwicky/dmtpc/software/geant4/env.sh
echo Environment setup complete.
