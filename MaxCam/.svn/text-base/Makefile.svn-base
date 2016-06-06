include	Makefile.root

###
### Is blank on purpose
### filled in if statements later; currently only used for the 4-shooter
### helps pick which slow control and pressure control source files
### to use when building DMSlow and DMPressureControl
###
DET_SUFFIX = 


##
## Set to non-empty to enable labjack 
##
#LABJACK = 0
LABJACK = 

###
### the 10L machines
###
SLOW_CONTROL_10L = pippin.lns.mit.edu
DAQ_10L = gimli.lns.mit.edu
#DAQ_10L = noether.lns.mit.edu

###
### the 4sh machines
###
SLOW_CONTROL_4SH = rubin.lns.mit.edu
#DAQ_4SH = dmtpc00.lns.mit.edu
DAQ_4SH = dmtpc00

###
### the RHUL test stand machine
###
#SLOW_CONTROL_RHUL = rubin.lns.mit.edu
DAQ_RHUL = lab-dm-003.pp.rhul.ac.uk

###
### the OWL machines
###
SLOW_CONTROL_OWL = nest.brynmawr.edu
DAQ_OWL = dmowl.brynmawr.edu
BMC_CLUSTER = powerwulf.brymawr.edu

### 
### destination for the web stuff
###
SCRIPTS_DIR = /var/www/cgi-bin/
HTML_DIR = /var/www/html/

###
### Directory where dbaccess.txt file (or link to it) is
###
DB_DIR = /etc/

###
### Check for libraries
###
MATHMORE_ENABLED=yes

ifeq ($(MATHMORE_ENABLED), yes)
  CXXFLAGS+= -DHAVE_MATHMORE
  LDFLAGS+= -lMathMore
endif

###
###   For analysis on Linux/Mac compile files needed for analysis only
###
OBJS_SCOPE_DAQ_ONLY =  ScopeDataChannel.o ScopeDataChannelCint.o \
	               ScopeWaveform.o ScopeWaveformCint.o \
	               ScopeConfig.o ScopeConfigCint.o  \
	               ScopeChannelConfig.o ScopeChannelConfigCint.o  \
	               ScopeTriggerConfig.o ScopeTriggerConfigCint.o  \
	               ScopeClockConfig.o ScopeClockConfigCint.o  \
	               ScopeBoardInfo.o ScopeBoardInfoCint.o \
			ScopeHandler.o ScopeAlazarATS860.o ScopeData.o Scope.o
#	               ScopeHandlerCint.o 

OBJS_CAEN_DAQ_ONLY =  ScopeCAENv1720.o ScopeCAEN.o ScopeHandlerCAEN.o CAENKeyComm.o

OBJS =  MaxCamConfig.o MaxCamConfigCint.o \
	MaxCamRead.o MaxCamReadCint.o \
	MaxCamQE.o MaxCamQECint.o \
	GausConvSpline.o\
	DmtpcPMTQE.o DmtpcPMTQECint.o \
	Dmtpc4ShooterStitcher.o Dmtpc4ShooterStitcherCint.o \
	DmtpcLensCorrection.o DmtpcLensCorrectionCint.o \
	DmtpcRootTools.o DmtpcRootToolsCint.o \
	DmtpcDecayChain.o DmtpcDecayChainCint.o \
	SimpleConfig.o \
	MaxCamSoftCal.o \
        MaxCamTrack.o MaxCamTrackCint.o  \
        MaxCamMC.o MaxCamMCCint.o \
	      MaxCamImageTools_ImageProcessing.o MaxCamImageTools_ImageTransforms.o \
	      MaxCamImageTools_PixelProcessing.o MaxCamImageTools_Interpolation.o \
	      MaxCamImageTools_ImageCalculations.o MaxCamImageTools_ClusterFinding.o \
			 	MaxCamImageTools_IO.o MaxCamImageToolsCint.o \
        MaxCamSRIM.o MaxCamSRIMCint.o \
	MaxCamTwoBodyKinematics.o MaxCamTwoBodyKinematicsCint.o \
        MaxCamDataset.o MaxCamDatasetCint.o \
	MaxCamENDF.o MaxCamENDFCint.o \
	DmtpcPaschen.o DmtpcPaschenCint.o \
        MaxCamWIMP.o MaxCamWIMPCint.o \
				DmtpcMath.o DmtpcMathCint.o\
	MaxCamStat.o MaxCamStatCint.o  \
	MaxCamVertex.o MaxCamVertexCint.o MaxCamUnits.o \
	MaxCamElectricField.o MaxCamElectricFieldCint.o \
	MaxCamSegment.o MaxCamSegmentCint.o \
	MaxCamCluster.o	MaxCamClusterCint.o Tools.o \
	MaxCamClusterImage.o MaxCamClusterImageCint.o \
	MaxCamChannel.o MaxCamChannelCint.o \
	DmtpcProjection.o DmtpcProjectionCint.o \
	DmtpcTrackFitter.o DmtpcTrackFitterCint.o \
	DmtpcSRIMProjection.o DmtpcSRIMProjectionCint.o \
	MaxCamChannelConfig.o\
        MaxCamImage.o MaxCamImageCint.o \
	MaxCamWaveformTools.o MaxCamWaveformToolsCint.o \
	MaxCamWaveform.o MaxCamWaveformCint.o \
	MaxCamPulse.o MaxCamPulseCint.o \
	MaxCamTriggerGroup.o MaxCamTriggerGroupCint.o \
        MaxCamSerial.o MaxCamSerialCint.o \
	MaxCamMFC.o MaxCamMFCCint.o \
	DmtpcDataset.o DmtpcDatasetCint.o \
	DmtpcSkimDataset.o DmtpcSkimDatasetCint.o \
	DmtpcStitchInfo.o DmtpcStitchInfoCint.o \
	DmtpcMCDataset.o DmtpcMCDatasetCint.o \
	DmtpcSkimPlaylist.o DmtpcSkimPlaylistCint.o \
        ScopeDataInfo.o ScopeDataInfoCint.o \
	ScopeWaveformData.o ScopeWaveformDataCint.o \
	DmtpcEvent.o DmtpcEventCint.o \
	DmtpcEventTable.o DmtpcEventTableCint.o DmtpcEventTableEntryCint.o\
	DmtpcSkimEvent.o DmtpcSkimEventCint.o \
	DmtpcDetectorPart.o DmtpcDetectorPartCint.o \
	DmtpcKeys.o DmtpcKeysCint.o\
	DmtpcAbsLocation.o DmtpcAbsLocationCint.o\
	DmtpcLocation.o DmtpcLocationCint.o \
	McDarkTrack.o McDarkTrackCint.o \
	McDarkDigi.o McDarkDigiCint.o \
        DmtpcEventViewer.o DmtpcEventViewerCint.o \
	DmtpcAstro.o DmtpcAstroCint.o \
	DmtpcTheory.o DmtpcTheoryCint.o \
	DmtpcRose.o DmtpcRoseCint.o \
	DmtpcCameraMap.o DmtpcCameraMapCint.o \
	DmtpcStringTools.o DmtpcStringToolsCint.o \
	SPEFile.o SPEFileCint.o \
	DmtpcSkimRunSummary.o DmtpcSkimRunSummaryCint.o    \
	DmtpcGainMap.o DmtpcGainMapCint.o              \
	RobinHoodTriangle.o RobinHoodTriangleCint.o \
	DmtpcDataConverter.o DmtpcDataConverterCint.o \
	RobinHoodTriangleMaker.o RobinHoodTriangleMakerCint.o



# Alazar scope information
CXXFLAGS_SCOPE = -I/usr/local/AlazarTech/include
LDFLAGS_SCOPE  = -L/usr/local/AlazarTech/lib -lPlxApi

# CAEN scope information
CXXFLAGS_CAEN = -I/net/zwicky/dmtpc/pietrog/new_cvs/projects/DarkMatter/CAEN/v1720/include
LDFLAGS_CAEN  = -L/net/zwicky/dmtpc/pietrog/new_cvs/projects/DarkMatter/CAEN/v1720/lib -lCAENDigitizer -lCAENVME

# CAEN scope information for RHUL
CXXFLAGS_CAEN_RHUL = -I/scratch2/pietrog/dmtpc/projects/DarkMatter/CAEN/include
LDFLAGS_CAEN_RHUL  = -L/scratch2/pietrog/dmtpc/projects/DarkMatter/CAEN/lib -lCAENDigitizer -lCAENVME

#
###
###   For DAQ on Linux compile all files
###

ifeq ($(HOSTNAME),$(SLOW_CONTROL_10L))
#LDFLAGS  += -lpng
LDFLAGS  += -lcomedi
LDFLAGS  += -L/usr/lib/mysql -lmysqlclient -lcurl
#LDFLAGS  += -L/usr/lib/mysql -lmysqlclient
CXXFLAGS += -DDM_DAQ -g
#LDFLAGS += -lserial
OBJS +=   DmtpcLoggableParam.o DmtpcLoggableParamCint.o
OBJS +=   DmtpcLogEntry.o DmtpcLogEntryCint.o
#OBJS +=   MaxCamMFC.o MaxCamMFCCint.o
#OBJS +=   MaxCamSerial.o MaxCamSerialCint.o
#CINTFLAGS += -p -I/usr/include/ 

DET_SUFFIX=

endif

ifeq ($(HOSTNAME),$(SLOW_CONTROL_4SH))

$(warning " ------------------- SLOW_CONTROL_4SH -------------------------") 

DET_SUFFIX=4sh

#LDFLAGS  += -lpng
LDFLAGS  += -lcomedi
LDFLAGS  += -L/usr/lib/mysql -lmysqlclient -lcurl
#LDFLAGS  += -L/usr/lib/mysql -lmysqlclient
CXXFLAGS += -DDM_DAQ
#LDFLAGS += -lserial
OBJS +=   DmtpcLoggableParam.o DmtpcLoggableParamCint.o
OBJS +=   DmtpcLogEntry.o DmtpcLogEntryCint.o
#OBJS +=   MaxCamMFC.o MaxCamMFCCint.o
#OBJS +=   MaxCamSerial.o MaxCamSerialCint.o
#CINTFLAGS += -p -I/usr/include/ 
endif

ifeq ($(HOSTNAME),$(DAQ_4SH)) 

$(warning " ------------------- DAQ_4SH -------------------------") 

LDFLAGS  += -lpng
LDFLAGS  += -lcomedi

# For Fingerlakes camera only
OBJS += MaxCam.o MaxCamCint.o
CXXFLAGS += -I../fli-dist-1.71/libfli -DCAM_FLI
LDFLAGS  += -L../fli-dist-1.71/libfli -lfli

# For Apogee Alta camera only
OBJS += MaxCamAlta.o MaxCamAltaCint.o
LDFLAGS  += /opt/apogee/lib/apogee_USB.so -L/opt/apogee/lib -ltcl8.3 -lccd -lfitsTcl -lusb
CXXFLAGS += -I/opt/apogee/include -I/opt/apogee/src/apogee -DLINUX -DCAM_ALTA

CXXFLAGS += -DDM_DAQ
CXXFLAGS += -I../cfitsio/include -I/usr/local/include

OBJS += MaxCamDummyCamera.o MaxCamDummyCameraCint.o \
	MaxCamCamera.o MaxCamCameraCint.o \
	MaxCamTest.o MaxCamTestCint.o MaxCamTestBase.o MaxCamTestBaseCint.o \
	DmtpcDAQ.o DmtpcDAQCint.o DmtpcRun.o DmtpcRunCint.o

# MySQL
LDFLAGS  +=  -L/usr/lib/mysql -lmysqlclient

# Alazar scope
OBJS     += $(OBJS_SCOPE_DAQ_ONLY)
CXXFLAGS += $(CXXFLAGS_SCOPE)
LDFLAGS  += $(LDFLAGS_SCOPE)

DET_SUFFIX=4sh

endif  # DAQ_4SH

ifeq ($(HOSTNAME),$(DAQ_RHUL))

$(warning " ------------------- DAQ_RHUL -------------------------")

CXXFLAGS += -I/scratch1/darkmatter/dmtpc/software/32bit/cfitsio

OBJS += MaxCamDummyCamera.o MaxCamDummyCameraCint.o \
        MaxCamCamera.o MaxCamCameraCint.o \
        MaxCamTest.o MaxCamTestCint.o MaxCamTestBase.o MaxCamTestBaseCint.o \
        DmtpcDAQ.o DmtpcDAQCint.o DmtpcRun.o DmtpcRunCint.o

# MySQL
LDFLAGS  += -L/usr/lib/mysql -lmysqlclient

# CAEN scope
OBJS     += $(OBJS_CAEN_DAQ_ONLY)
CXXFLAGS += $(CXXFLAGS_CAEN_RHUL) -DSCOPE_CAEN
LDFLAGS  += $(LDFLAGS_CAEN_RHUL) -L/scratch1/darkmatter/dmtpc/software/32bit/lib -lcfitsio -lgsl -lgslcblas

DET_SUFFIX=RHUL

endif  # DAQ_RHUL

ifeq ($(HOSTNAME),$(SLOW_CONTROL_OWL))

$(warning " ------------------- SLOW_CONTROL_OWL -------------------------") 

DET_SUFFIX=OWL

CXXFLAGS += -DDM_DAQ
#LDFLAGS  += -lpng
LDFLAGS  += -lcomedi
LDFLAGS  += -L/usr/lib64/mysql -lmysqlclient 
LDFLAGS  += -L/usr/lib64 -lcurl
#LDFLAGS += -lserial
#OBJS +=   MaxCamSerial.o MaxCamSerialCint.o
#CINTFLAGS += -p -I/usr/include/ 
endif

ifeq ($(HOSTNAME),$(DAQ_10L)) 
LDFLAGS  += -L.
LDFLAGS  += -lpng
LDFLAGS  += -lcomedi
LDFLAGS  += -lmysqlclient
CXXFLAGS += -DDM_DAQ

DET_SUFFIX = 10L

endif  #DAQ_10L

ifeq ($(HOSTNAME),mitdm11.mit.edu)

# Don't worry about test for Raytheon Detector
#OBJS +=	MaxCamTest.o MaxCamTestCint.o \
#	MaxCamTestBase.o MaxCamTestBaseCint.o

# For Fingerlakes camera only - NOT in Raytheon Detector
#OBJS += MaxCam.o MaxCamCint.o
#CXXFLAGS += -I../fli-dist-1.71/libfli -DCAM_FLI
#LDFLAGS  += -L../fli-dist-1.71/libfli -lfli

# For Apogee Alta camera only
OBJS += MaxCamAlta.o MaxCamAltaCint.o
LDFLAGS  += /opt/apogee/lib/apogee_USB.so -L/opt/apogee/lib -ltcl8.3 -lccd -lfitsTcl -lusb
CXXFLAGS += -I/opt/apogee/include -I/opt/apogee/src/apogee -DLINUX -DCAM_ALTA

# For Andor iXon camera only
OBJS += MaxCamAndor.o MaxCamAndorCint.o
LDFLAGS += -landor
CXXFLAGS += -DCAM_ANDOR

OBJS +=	MaxCamDummyCamera.o MaxCamDummyCameraCint.o \
	MaxCamCamera.o MaxCamCameraCint.o \
	DmtpcDAQ.o DmtpcDAQCint.o DmtpcRun.o DmtpcRunCint.o

LD=g++
LDFLAGS += -lpng
LDFLAGS += -lcomedi
CXXFLAGS += -I../cfitsio/include -I.

LDFLAGS  +=  -L/usr/lib/mysql -lmysqlclient
CXXFLAGS += -DDM_DAQ
# Alazar scope
OBJS     += $(OBJS_SCOPE_DAQ_ONLY)
CXXFLAGS += $(CXXFLAGS_SCOPE)
LDFLAGS  += $(LDFLAGS_SCOPE)
endif

ifeq ($(HOSTNAME),$(DAQ_OWL)) 
$(warning " ------------------- DAQ_OWL -------------------------") 
# For Fingerlakes camera only
OBJS += MaxCam.o MaxCamCint.o
#CXXFLAGS += -I../fli-dist-1.71/libfli -DCAM_FLI
#LDFLAGS  += -L../fli-dist-1.71/libfli -lfli
CXXFLAGS += -I../libfli-1.104 -DCAM_FLI 
#CXXFLAGS += -DDM_DAQ
LDFLAGS  += -L../libfli-1.104 -lfli
LDFLAGS  += -lpng

OBJS += MaxCamDummyCamera.o MaxCamDummyCameraCint.o \
	MaxCamCamera.o MaxCamCameraCint.o 

#	MaxCamTest.o MaxCamTestCint.o MaxCamTestBase.o MaxCamTestBaseCint.o \
#	DmtpcDAQ.o DmtpcDAQCint.o DmtpcRun.o DmtpcRunCint.o

endif


ifeq ($(HOSTNAME),$(DAQ_OWL)) 
# For Fingerlakes camera only
OBJS += MaxCam.o MaxCamCint.o
#CXXFLAGS += -I../fli-dist-1.71/libfli -DCAM_FLI
#LDFLAGS  += -L../fli-dist-1.71/libfli -lfli
CXXFLAGS += -I../libfli-1.104 -DCAM_FLI
LDFLAGS  += -L../libfli-1.104 -lfli
LDFLAGS  += -lpng

OBJS += MaxCamDummyCamera.o MaxCamDummyCameraCint.o \
	MaxCamCamera.o MaxCamCameraCint.o

#	MaxCamTest.o MaxCamTestCint.o MaxCamTestBase.o MaxCamTestBaseCint.o \
#	DmtpcDAQ.o DmtpcDAQCint.o DmtpcRun.o DmtpcRunCint.o

endif


### Labjack 

ifdef LABJACK

OBJS +=  LabJackU3.o LabJackU3Cint.o u3.o
LDFLAGS += -llabjackusb 

endif 

### ... for all Linux machines and MACOSX compile as analysis

CXXFLAGS += -I/usr/local/lib -I../cfitsio/include -g
LDFLAGS  += -L../cfitsio -lcfitsio -lSpectrum -lGeom -lgsl -lgslcblas -lMinuit 
CXXFLAGS += -I/usr/local/include -L/usr/local/include
LDFLAGS  += -lcfitsio -lSpectrum -lGeom -lgslcblas -lMinuit -L/app/d-Chooz/Software/gsl-1.15/lib -lgsl 
CXXFLAGS += -I/usr/local/include -L/usr/local/include -L/app/d-Chooz/Software/include
CXXFLAGS += -I/app/d-Chooz/Software/include
CXXFLAGS += -I/usr/local/lib -I../cfitsio/include -g

all:	$(OBJS)
	@echo "Making Library for $(ARCH) on $(HOSTNAME)"
	$(LD) $(SOFLAGS) $(LDFLAGS) $(ROOTGLIBS) $(OBJS) $(OutPutOpt) MaxCam_$(ARCH).so
# must also build MaxCam waveformtools
	make -C waveformtools/
ifneq ($(DET_SUFFIX),)
	touch src/DMSlow$(DET_SUFFIX).cc
	touch src/DMPressureControl$(DET_SUFFIX).cc
#	ln -sf $(CURDIR)/daq/dbaccess_$(DET_SUFFIX).txt $(CURDIR)/daq/dbaccess.txt
	cp -f daq/dbaccess_$(DET_SUFFIX).txt daq/dbaccess.txt
endif
	@echo "Changing permissions..."
	#chmod g+w * 
	#chgrp dmtpc * 

shared:	$(OBJS)
	@echo "Making Shared Library for $(ARCH) on $(HOSTNAME)"
#ifeq ($(ARCH),macosx)
#	$(LD) $(SOFLAGS)libMaxCam.so $(LDFLAGS) $(ROOTLIBS) $(OBJS) $(OutPutOpt) libMaxCam.so
#	ln -sf libMaxCam.so libMaxCam.dylib
#else
#	$(LD) $(SOFLAGS) $(LDFLAGS) $(ROOTGLIBS) $(OBJS) $(OutPutOpt) libMaxCam.so
# apparently the LDFLAGS and ROOTGLIBS have to come after the OBJS in order for the 
# ROOT libraries to be explicitly linked to libMaxCam.so (verify with ldd)
	$(LD) $(SOFLAGS) $(OBJS) $(LDFLAGS) $(ROOTGLIBS) $(OutPutOpt) libMaxCam.so
#endif
clean:
	@rm -f $(OBJS) *Cint.cc *Cint.h core *Cint.o
	make -C waveformtools/ clean

tar:
	tar cvf MaxCam.tar MaxCam*cc MaxCam*hh Makefile Makefile.root examples/*.cxx tables/*
	gzip MaxCam.tar


.SUFFIXES: .cc

###

.SECONDARY: 
 
%Cint.cc: %.hh
	@echo "Generating dictionary ..."
#	Add -I directory for Alazar header files
#	$(ROOTSYS)/bin/rootcint -f $@ -c $(CXXFLAGS_SCOPE) $<
#	$(ROOTSYS)/bin/rootcint -f $@ -c $<
#	$(ROOTSYS)/bin/rootcint -f $@ -c $(CINTFLAGS) $<

#Use linkdef file if it exists
	rootcint -f $@ -c $(CINTFLAGS) $< `ls $*Linkdef.h 2>/dev/null`

u3.o:
	gcc -c u3.c

.cc.o:
	@echo "Compiling... for $(ARCH) on $(HOSTNAME)"
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) -c $<

install: all
#	mv MaxCam_$(ARCH).so /usr/local/lib/MaxCam.so
	cp MaxCam_$(ARCH).so MaxCam.so

#	cp bin/DMPressureControl $(SCRIPTS_DIR) 
#	cp bin/DMSlow $(SCRIPTS_DIR) 
#	ln -sf $(CURDIR)/daq/dbaccess_$(DET_SUFFIX).txt $(DB_DIR)/dbaccess.txt
#	cp -f daq/dbaccess_$(DET_SUFFIX).txt $(DB_DIR)/dbaccess.txt

update:
	cvs update *.cc *.hh Makefile

JBOBJS = DmtpcDataset.o DmtpcDatasetCint.o \
         DmtpcEvent.o DmtpcEventCint.o \
         MaxCamChannel.o MaxCamChannelCint.o \
         MaxCamImage.o MaxCamImageCint.o \
         MaxCamConfig.o MaxCamConfigCint.o

caentest: runCAENv1720.o $(OBJS_CAEN_DAQ_ONLY)
	$(CXX) $(LDFLAGS) $(ROOTGLIBS) runCAENv1720.o $(OBJS_CAEN_DAQ_ONLY) -o caentest

scopetest: scopetest.o $(OBJS_SCOPE_DAQ_ONLY)
	$(CXX) $(LDFLAGS) $(ROOTLIBS) scopetest.o $(OBJS_SCOPE_DAQ_ONLY) -o scopetest 

scopeplan: scopeplan.o $(OBJS_SCOPE_DAQ_ONLY) $(JBOBJS)
	$(CXX) $(LDFLAGS) $(ROOTLIBS) $(OBJS_SCOPE_DAQ_ONLY) $(JBOBJS) scopeplan.o -o scopeplan

scopems: scopems.o Utils.o
	$(CXX) $(LDFLAGS)  scopems.o Utils.o -o scopems


#scopeplan: scopeplan.o $(OBJS_SCOPE_ALL_USERS) $(OBJS_SCOPE_DAQ_ONLY) $(JBOBJS)
#	$(CXX) $(LDFLAGS) $(ROOTLIBS) $(OBJS_SCOPE_ALL_USERS) $(OBJS_SCOPE_DAQ_ONLY) $(JBOBJS) scopeplan.o -o scopeplan

slow: bin/DMSlow 
bin/DMSlow: src/DMSlow$(DET_SUFFIX).o
	$(LD)  $(LDFLAGS)  MaxCamSerial.o MaxCamSerialCint.o SimpleConfig.o MaxCamChannel.o MaxCamChannelConfig.o MaxCamSoftCal.o  MaxCamChannelCint.o src/DMSlow$(DET_SUFFIX).o $(ROOTLIBS) -o bin/DMSlow

rslow: bin/RSlow 
bin/RSlow: src/RSlow.o src/u3.o
	$(LD)  $(LDFLAGS)  -llabjackusb src/RSlow.o src/u3.o MaxCamSerial.o MaxCamSerialCint.o SimpleConfig.o MaxCamChannel.o MaxCamChannelConfig.o MaxCamSoftCal.o MaxCamChannelCint.o  $(ROOTLIBS) -o bin/RSlow

bin/DMKillSynaccess: src/DMKillSynaccess.o
	$(LD)  $(LDFLAGS)   MaxCamSerial.o MaxCamSerialCint.o src/DMKillSynaccess.o  $(ROOTLIBS) -o bin/DMKillSynaccess

bin/DMRampDown: src/DMRampDown.o
	$(LD)  $(LDFLAGS)  MaxCamSerial.o MaxCamSerialCint.o SimpleConfig.o MaxCamChannel.o MaxCamChannelConfig.o MaxCamSoftCal.o  MaxCamChannelCint.o src/DMRampDown.o $(ROOTLIBS) -o bin/DMRampDown

alarm: bin/DMAlarms bin/DMKillSynaccess bin/DMRampDown bin/DM_UPS

bin/DMAlarms: src/DMAlarms.o
	$(LD)  $(LDFLAGS)  MaxCamSerial.o MaxCamSerialCint.o SimpleConfig.o MaxCamChannel.o MaxCamChannelConfig.o MaxCamSoftCal.o  MaxCamChannelCint.o src/DMAlarms.o $(ROOTLIBS) -o bin/DMAlarms

bin/DM_UPS: src/DM_UPS.c 
	gcc src/DM_UPS.c -o bin/DM_UPS 

bin/mcc_usb_dacctl: src/mcc_usb_dacctl.c 
	gcc -g src/mcc_usb_dacctl.c -o bin/mcc_usb_dacctl -lhid -lmcchid -lusb -lm 


press: bin/DMPressureControl
bin/DMPressureControl: src/DMPressureControl$(DET_SUFFIX).o
	$(LD)  $(LDFLAGS)  MaxCamChannel.o MaxCamSerial.o MaxCamSerialCint.o  MaxCamChannelCint.o src/DMPressureControl$(DET_SUFFIX).o MaxCamChannelConfig.o SimpleConfig.o MaxCamSoftCal.o $(ROOTLIBS) -o bin/DMPressureControl

deploy_web: 
	cp -r script/* $(SCRIPTS_DIR)
	cp -r html/* $(HTML_DIR)
	cp -r channels $(HTML_DIR)

src/%.o: src/%.cc 
	@echo 'Compiling ...'
	$(CXX) $(CXXFLAGS)  $(ROOTCFLAGS) -I./ -o $@ -c $<

DataConverter: DataConverter.cc libMaxCam.so
	$(CXX) -g $(CXXFLAGS) $(LDFLAGS) $(ROOTLIBS)  -L./ -lMaxCam DataConverter.cc -o bin/DataConverter


DMCounter: examples/DMCounter.cc LabJackU3.o u3.o 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(ROOTLIBS) -L./ -lMaxCam examples/DMCounter.cc -o bin/DMCounter	

check: libMaxCam.so
	make -C test/ -k all 


jbjunk:
	@echo 'jbjunk ...'
	@echo $(CURDIR)
