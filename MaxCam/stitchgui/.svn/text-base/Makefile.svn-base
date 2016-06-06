CC=g++
CFLAGS=-Os -g `root-config --cflags` -Wall -Wextra
LDFLAGS = `root-config --glibs` -L../ -lMaxCam -L../waveformtools/lib -lWaveformTools 


DmtpcStitchGUI: DmtpcStitchGUI.o DmtpcStitchGUIFrame.o DmtpcStitchGUIDict.o
	$(CC) -o $@ $(LDFLAGS) $^

DmtpcStitchGUI.o: DmtpcStitchGUI.cc
	$(CC) $(CFLAGS) -c $?

DmtpcStitchGUIFrame.o: DmtpcStitchGUIFrame.cc
	$(CC) $(CFLAGS) -c $?

DmtpcStitchGUIDict.o : DmtpcStitchGUIDict.cxx
	$(CC) $(CFLAGS) -c $?

DmtpcStitchGUIDict.cxx : DmtpcStitchGUIFrame.hh
	rootcint -f $@ -c $^

clean:
	rm -f *.o DmtpcStitchGUI *Dict.cxx *Dict.h
