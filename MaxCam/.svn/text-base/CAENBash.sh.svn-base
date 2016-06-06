############
# for ROOT #
############

ROOTSYS=/scratch1/darkmatter/dmtpc/software/32bit/root
PATH=$ROOTSYS/bin:$PATH
if [ ! -d $ROOTSYS ]; then
  echo "Warning: ROOTSYS directory $ROOTSYS not found"
fi
if [ -n "$LD_LIBRARY_PATH" ]; then
  LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
else
  LD_LIBRARY_PATH=$ROOTSYS/lib
fi
export ROOTSYS PATH LD_LIBRARY_PATH

################
# for CAEN lib #
################

export PATH=$PATH:/usr/lib
export PATH=$PATH:/usr/include

CAENLIB=/scratch2/pietrog/dmtpc/projects/DarkMatter/CAEN/lib
export LD_LIBRARY_PATH=$CAENLIB:$LD_LIBRARY_PATH

###################
# GSL and cfitsio #
###################

CFITSIOLIB=/scratch1/darkmatter/dmtpc/software/32bit/cfitsio
GSLLIB=/scratch1/darkmatter/dmtpc/software/32bit/lib

export LD_LIBRARY_PATH=$CFITSIOLIB:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$GSLLIB:$LD_LIBRARY_PATH