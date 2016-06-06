############################################################
#							   #
# Makefile for the ScopeCAEN executable CAEN	           #
# This should help you read data from the CAEN board v1720 #
# Straight into a DmtpcDataset root output form            #
#							   #
# Pietro Giampa, RHUL 2013    				   #
#							   #
############################################################

#include the Makefile root file 
include ../Makefile.root

# Basic Definition of cxx and flags
CXX = g++

# Set Flags for MaxCam and WaveFormtoools
LDFLAGS = -L/net/zwicky/dmtpc/pietrog/new_cvs/projects/DarkMatter/MaxCam -lMaxCam
LDFLAGS += -L/net/zwicky/dmtpc/pietrog/new_cvs/projects/DarkMatter/MaxCam/waveformtools/lib -lWaveformTools
LDFLAGS += -L/net/zwicky/dmtpc/pietrog/new_cvs/projects/DarkMatter/fftw/fftw-3.3.2/lib -lfftw3

# Define objs needed for the CAEN executable
OBJS = tmp/ScopeCAENv1720.o tmp/CAENKeyComm.o \

# Set Flags for CAEN
CXXFLAGS += -I/net/zwicky/dmtpc/pietrog/new_cvs/projects/DarkMatter/CAEN/v1720/include
LDFLAGS += -L/net/zwicky/dmtpc/pietrog/new_cvs/projects/DarkMatter/CAEN/v1720/lib -lCAENDigitizer -lCAENVME

# Compiling step
#.cc.o:
tmp/%.o: src/%.cc
	@echo "Compiling Everything for CAEN"
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) -o $@ -c $<

# Clean Up Everuthing
clean:
	@echo "Cleaning Everything for CAEN"
	@rm -f $(OBJS) *tmp/runCAENv1720.o *caen

# Set Up the caentest executable
caen:	tmp/runCAENv1720.o $(OBJS)
	$(CXX) $(LDFLAGS) $(ROOTGLIBS) tmp/runCAENv1720.o $(OBJS) -o caen
