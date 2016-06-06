#FFT_SUPPORT=no
FFT_SUPPORT=$(shell root-config --has-fftw3)

CC=g++
CFLAGS=-O0 -g `root-config --cflags` -Wall -Wextra
LFLAGS=`root-config --glibs` -lMathMore -L../ -lMaxCam  -L../waveformtools/lib/ -lWaveformTools


DmtpcSkimViewer: DmtpcSkimViewerScopeFrame.o DmtpcSkimViewerDict.o DmtpcSkimViewerFrame.o DmtpcSkimViewerImageTransform.o DmtpcSkimViewerProjection.o DmtpcSkimViewer.o style.o
	$(CC) -o $@ $(LFLAGS) $^

shared: DmtpcSkimViewerScopeFrame.o DmtpcSkimViewerDict.o DmtpcSkimViewerFrame.o DmtpcSkimViewerImageTransform.o style.o
	$(CC) -shared $(LFLAGS) $^ -o viewer.so

install: 
	cp *.png /net/zwicky/dmtpc/software/viewer/
	cp Readme /net/zwicky/dmtpc/software/viewer/
	cp DmtpcSkimViewer /net/zwicky/dmtpc/software/viewer/

DmtpcSkimViewer.o: DmtpcSkimViewer.cc	
	$(CC) $(CFLAGS) -c $?


style.o: style.cc  
	$(CC) $(CFLAGS) -c $?

DmtpcSkimViewerDict.o: DmtpcSkimViewerDict.cxx  
	$(CC) $(CFLAGS) -c $?

DmtpcSkimViewerImageTransform.o: DmtpcSkimViewerImageTransform.cc
	$(CC) $(CFLAGS) -c $?

DmtpcSkimViewerFrame.o: DmtpcSkimViewerFrame.cc  
ifeq ($(FFT_SUPPORT),yes)
	$(CC) $(CFLAGS) -DHAVE_FFT -c $?
else
	$(CC) $(CFLAGS) -c $?
endif

DmtpcSkimViewerScopeFrame.o: DmtpcSkimViewerScopeFrame.cc
	$(CC) $(CFLAGS) -c $?

DmtpcSkimViewerProjection.o: DmtpcSkimViewerProjection.cc
	$(CC) $(CFLAGS) -c $?

DmtpcSkimViewerDict.cxx : DmtpcSkimViewerImageTransform.hh DmtpcSkimViewerScopeFrame.hh DmtpcSkimViewerFrame.hh DmtpcSkimViewerProjection.hh DmtpcSkimViewerLinkDef.h
	rootcint -f $@ -c $^

clean:
	rm -f DmtpcSkimViewer
	rm *.o *Dict.cxx *Dict.h
