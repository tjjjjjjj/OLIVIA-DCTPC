name := McDark
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
   G4INSTALL = /usr/local/geant4
endif

include $(G4INSTALL)/config/architecture.gmk

ifndef G4LISTS_BASE
   EXTRALIBS += -L$(G4LIB)/$(G4SYSTEM)
   G4LISTS_BASE = $(G4INSTALL)/physics_lists/hadronic
else
   EXTRALIBS += -L$(G4LIB)/$(G4SYSTEM)
endif


.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

# ROOT stuff
INCFLAGS += -I$(ROOTSYS)/include -I../MaxCam
LDFLAGS += -L$(ROOTSYS)/lib
LDLIBS_PREFINAL += $(shell $(ROOTSYS)/bin/root-config --libs) -L ../MaxCam -lMaxCam -L../MaxCam/waveformtools/lib  -lWaveformTools -L../cfitsio -lcfitsio



