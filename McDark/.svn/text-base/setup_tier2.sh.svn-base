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
export G4WORKDIR=$PWD
echo Work directory set to $G4WORKDIR.

if [ ! -d $G4WORKDIR/out ]; then 
  mkdir $G4WORKDIR/out
fi

# Please change this as well if you want the output in a different directory
export ROOTOUT=$G4WORKDIR/out
export G4HEPREPFILE_DIR=$ROOTOUT
echo ROOT output and HepRepFile output currently set to $ROOTOUT

export ROOTSYS=/app/d-Chooz/Software/root/root_v5.30.06



export LD_LIBRARY_PATH=$G4LIB/plists/Linux-g++/:$G4LIB/Linux-g++:/app/d-Chooz/Software/CLHEP/CLHEP2042/lib:/net/t2srv0008/app/d-Chooz/Software/geant4/test/geant4.9.4.p03/lib/Linux-g++:$ROOTSYS/lib:/usr/lib:/net/hisrv0001/home/spitzj/cern/2006b/x86_64-slc5-gcc34-opt/lib:$ROOTSYS/lib:/net/hisrv0001/home/spitzj/extralibs:/usr/lib:/app/d-Chooz/Software/geant4/test/geant4.9.4.p03/extralibs:/app/d-Chooz/Software/CLHEP/CLHEP2042/lib:/app/d-Chooz/Software/root/root/lib:/usr/lib:/app/d-Chooz/Software/geant4/test/geant4.9.4.p03/extralibs:$LD_LIBRARY_PATH


#source /net/zwicky/dmtpc/software/geant4/env.sh
unset G4WORKDIR
#source /afs/fnal.gov/files/code/e898/code6/MiniBooNE_scripts/MiniBooNE.login
export G4INSTALL=/net/t2srv0008/app/d-Chooz/Software/geant4/test/geant4.9.4.p03/
#export ROOTSYS=/net/t2srv0008/app/d-Chooz/Software/root/root

export PATH=$ROOTSYS/bin:$PATH
export G4LIB=/net/t2srv0008/app/d-Chooz/Software/geant4/test/geant4.9.4.p03/lib
export G4INCLUDE=/net/t2srv0008/app/d-Chooz/Software/geant4/test/geant4.9.4.p03//include
export CLHEP_BASE_DIR=/app/d-Chooz/Software/CLHEP/CLHEP2042
export CLHEP_INCLUDE_DIR=/app/d-Chooz/Software/CLHEP/CLHEP2042/include
export CLHEP_LIB_DIR=/app/d-Chooz/Software/CLHEP/CLHEP2042/lib
export CLHEP_LIB=CLHEP
export G4WORKDIR=$PWD
#source /afs/fnal.gov/files/data/argoneut/d01/lar/external/CRY/setup
#export G4VIS_USE=1
export G4SYSTEM=Linux-g++
export NeutronHPCrossSections=/net/t2srv0008/app/d-Chooz/Software/geant4/test/geant4.9.4.p03/data/G4NDL3.14
export G4NEUTRONHPDATA=/net/t2srv0008/app/d-Chooz/Software/geant4/test/geant4.9.4.p03/data/G4NDL3.14
export G4LEVELGAMMADATA=$G4INSTALL/data/PhotonEvaporation2.1
export G4RADIOACTIVEDATA=$G4INSTALL/data/RadioactiveDecay3.3
export G4LEDATA=$G4INSTALL/data/G4EMLOW6.23
export G4ANALYSIS_USE_ROOT=1
#setenv GCC_DIR /afs/fnal.gov/ups/gcc/v3_4_3/Linux+2
#setenv COMPILER_PATH /afs/fnal.gov/ups/gcc/v3_4_3/Linux+2/bin
#export OGLHOME=/usr/lib
#export G4VIS_BUILD_OPENGLX_DRIVER=1
#export G4VIS_USE_OPENGLX=1
#export G4LIB_BUILD_SHARED=1
echo Environment setup complete.
