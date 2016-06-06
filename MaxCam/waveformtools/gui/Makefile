ROOTLIBS = $(shell root-config --libs) $(shell root-config --glibs)
ROOTCFLAGS = $(shell root-config --cflags)

CC = g++
CINT = rootcint

OBJS = GuiFirFrame.o GuiFirFrameDict.C

CFLAGS = $(ROOTCFLAGS) -Wall -Wextra

LIBS = $(ROOTLIBS)

all: $(OBJS) FirGUI.cc
	$(CC) $(CFLAGS) $(LIBS) -I./ $(OBJS) FirGUI.cc -o FirFilterGui

clean:
	rm -f $(OBJS) FirFilterGui FirFilterGui.Dict *Dict.cxx *Dict.h


.PHONY: all, clean, shared


GuiFirFrameDict.C: GuiFirFrame.hh
	$(CINT) -f $@ -c $^

GuiFirFrame.o: GuiFirFrame.cc
	$(CC) $(CFLAGS) -c $?
