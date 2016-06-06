ROOTCONFIG	:= "root-config"
CXXFLAGS 	:= $(shell $(ROOTCONFIG) --cflags) -g 
LDFLAGS 	:= $(shell $(ROOTCONFIG) --ldflags) 
LIBS 		:= $(shell $(ROOTCONFIG) --libs) 
INCDIR  	:= $(shell $(ROOTCONFIG) --incdir)
ARCH		:= $(shell $(ROOTCONFIG) --arch)


SOFLAGS		:= -shared
CXX	:= g++ 
LD	:= g++

ifeq ($(ARCH),linux)
  CXXFLAGS	+= -O -Wall -fPIC
  LDFLAGS	+= -O
  SOFLAGS 	= -shared
endif

ifeq ($(ARCH),linuxx8664gcc)
  CXXFLAGS      += -O -Wall -fPIC
  LDFLAGS       += -O
  SOFLAGS       = -shared
endif

ifeq ($(ARCH),macosx)
  MACOSX_MINOR := $(shell sw_vers | sed -n 's/ProductVersion://p' | cut -d . -f 2)
  MACOSXTARGET := MACOSX_DEPLOYMENT_TARGET=10.$(MACOSX_MINOR)
  LD 	= $(MACOSXTARGET) c++
  CXXFLAGS += -O2 -pipe -Wall -W -Woverloaded-virtual -fPIC
  SOFLAGS  = -arch i386 -dynamiclib -single_module -undefined dynamic_lookup
endif

ifeq ($(ARCH),macosx64)
  MACOSX_MINOR := $(shell sw_vers | sed -n 's/ProductVersion://p' | cut -d . -f 2)
  MACOSXTARGET := MACOSX_DEPLOYMENT_TARGET=10.$(MACOSX_MINOR)
  LD	= $(MACOSXTARGET) c++
  CXXFLAGS += -O2 -pipe -Wall -W -Woverloaded-virtual -fPIC
  SOFLAGS  = -arch x86_64 -dynamiclib -single_module -undefined dynamic_lookup
endif

ifeq ($(MY_INSTALL_DIR),)
  INSTALL_DIR=/usr/local
else
  INSTALL_DIR=$(MY_INSTALL_DIR)
endif


#Uncomment below to enable resampling
#CXXFLAGS += -DHAVE_LIBSAMPLERATE
#LDFLAGS += -lsamplerate

#check fft support
FFT_SUPPORT=$(shell $(ROOTCONFIG) --has-fftw3)


OBJS=  tmp/WaveformToolsDict.Cint \
       tmp/DmtpcPulse.o \
       tmp/SkimWaveform.o \
       tmp/WaveformVector.o \
       tmp/CspPulse.o \
       tmp/CspWaveform.o \
       tmp/CspWfVector.o \
       tmp/PMTPulse.o \
       tmp/PMTWaveform.o \
       tmp/PMTWfVector.o \
       tmp/FastPulse.o \
       tmp/FastWaveform.o \
       tmp/FastWfVector.o \
       tmp/WaveformTools.o \
       tmp/WaveformAnalysis.o \
       tmp/FirFilter.o \
       tmp/FirKernels.o
#       tmp/WaveformT.o \
       tmp/WfVectorT.o\


CINTFLAGS = -I../ -Iinclude
CXXFLAGS += -I./include -I../ 
CXXFLAGS += -I/app/d-Chooz/Software/fftw-3.3.3/include
CXXFLAGS += -L ../ -lMaxCam
CXXFLAGS += -L ./ -lMaxCam
HEADERS = DmtpcPulse.hh \
          SkimWaveform.hh \
          WaveformVector.hh \
          CspPulse.hh \
          CspWaveform.hh \
          CspWfVector.hh\
          PMTPulse.hh \
          PMTWaveform.hh \
          PMTWfVector.hh\
          FastPulse.hh \
          FastWaveform.hh \
          FastWfVector.hh \
          WaveformTools.hh\
          WaveformAnalysis.hh \
          FirFilter.hh \
          FirKernels.hh \


ifeq ($(FFT_SUPPORT),yes)
  OBJS += tmp/WaveformFFT.o tmp/Deconvolution.o 
  HEADERS += WaveformFFT.hh Deconvolution.hh
  LDFLAGS += -L/app/d-Chooz/Software/fftw-3.3.3/lib -lfftw3
endif

.PHONY: all install uninstall clean local_libdir tmpdir install_dir install_incdir install_libdir

all: local_libdir tmpdir $(OBJS)
	$(LD) $(SOFLAGS) $(LDFLAGS) $(LIBS) $(OBJS) -o lib/libWaveformTools.so

install: all install_dir install_incdir install_libdir

	@echo "Installing into $(INSTALL_DIR)"
	cp lib/libWaveformTools.so $(INSTALL_DIR)/lib/
	cp include/*.hh $(INSTALL_DIR)/include/
	cp include/*.h $(INSTALL_DIR)/include/

uninstall:
	@echo "Uninstalling headers and libraries from $(INSTALL_DIR)"
	rm -f $(INSTALL_DIR)/lib/libWaveformTools.so \
	$(INSTALL_DIR)/include/DmtpcPulse.hh \
	$(INSTALL_DIR)/include/SkimWaveform.hh \
	$(INSTALL_DIR)/include/WaveformVector.hh \
	$(INSTALL_DIR)/include/CspPulse.hh \
	$(INSTALL_DIR)/include/CspWaveform.hh \
	$(INSTALL_DIR)/include/CspWfVector.hh \
	$(INSTALL_DIR)/include/PMTPulse.hh \
	$(INSTALL_DIR)/include/PMTWaveform.hh \
	$(INSTALL_DIR)/include/PMTWfVector.hh \
	$(INSTALL_DIR)/include/FastPulse.hh \
	$(INSTALL_DIR)/include/FastWaveform.hh \
	$(INSTALL_DIR)/include/FastWfVector.hh \
	$(INSTALL_DIR)/include/WaveformTools.hh \
	$(INSTALL_DIR)/include/WaveformAnalysis.hh \
	$(INSTALL_DIR)/include/DmtpcWaveformLinkDef.h

clean:  uninstall
	rm -f -r tmp/
	rm -f lib/*.so lib/*.a lib/*.o

.SUFFIXES: .cc .Cint

install_dir:
	if test ! -d $(INSTALL_DIR); \
	then echo "Creating $(INSTALL_DIR)"; \
	mkdir $(INSTALL_DIR);\
	fi

install_incdir: install_dir
	if test ! -d $(INSTALL_DIR)/include; \
	then echo "Creating $(INSTALL_DIR)/include/"; \
	mkdir $(INSTALL_DIR)/include;\
	fi

install_libdir: install_dir
	if test ! -d $(INSTALL_DIR)/lib; \
	then echo "Creating $(INSTALL_DIR)/lib/"; \
	mkdir $(INSTALL_DIR)/lib;\
	fi


local_libdir:
	if test ! -d lib; \
	then echo "Creating ./lib/"; \
	mkdir lib;\
	fi

tmpdir:
	if test ! -d tmp; \
	then echo "Creating ./tmp/"; \
	mkdir tmp;\
	fi

tmp/%Cint.cc:
	@echo "Generating dictionary..."
	rootcint -f $@ -c $(CINTFLAGS) $(HEADERS) include/DmtpcWaveformLinkDef.h

tmp/%.Cint: tmp/%Cint.cc
	$(CXX) $(CXXFLAGS) -o $@ -c $<

tmp/%.o: src/%.cc
	$(CXX) $(CXXFLAGS) -o $@ -c $<
