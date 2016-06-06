CC=g++

CFLAGS= -O3 -g `root-config --cflags` -Wno-write-strings -I../../../MaxCam -I../../../MaxCam/waveformtools -I../../../MaxCam/waveformtools/include
LFLAGS= `root-config --libs` -lMathMore -lSpectrum -L../../../MaxCam -lMaxCam -L../../../MaxCam/waveformtools/lib -lWaveformTools -L /net/zwicky/dmtpc/software/tmva/lib -lTMVA




OBJS =  cleanSkimFunctions.o cleanSkimConfig.o cleanSkimProcessor.o 

cleanSkim: cleanSkim.o $(OBJS) 
	$(CC) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

testgm: testgm.o $(OBJS)  
	$(CC) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

MCCompare: MCCompare.o
	$(CC) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

countUnique: countUnique.o
	$(CC) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

MCPhiTest: MCPhiTest.o
	$(CC) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

AlphaFitTest: AlphaFitTest.o
	$(CC) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

DataFitTest: DataFitTest.o
	$(CC) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

DataFitTestShawnTree: DataFitTestShawnTree.o
	$(CC) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

MCFitTest: MCFitTest.o
	$(CC) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

MCFitTestTruth: MCFitTestTruth.o
	$(CC) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

MCPhiTestTruth: MCPhiTestTruth.o
	$(CC) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

Fit2DMC: Fit2DMC.o
	$(CC) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

Fit2DPlaylist: Fit2DPlaylist.o
	$(CC) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

%.o: %.cc
	$(CC) $(CFLAGS) -c $< 

clean: 
	rm -f *.o *.so




