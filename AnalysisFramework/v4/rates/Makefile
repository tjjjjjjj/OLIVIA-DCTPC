CC=g++

CFLAGS=-O3 -g `root-config --cflags` -Wall -I../../../MaxCam -I../../../MaxCam/waveformtools -I../../../MaxCam/waveformtools/include -I../src/
LFLAGS=`root-config --libs` -lSpectrum -lGeom -L../../../MaxCam -lMaxCam -L../../../MaxCam/waveformtools/lib -lWaveformTools -L../lib -lAnalysis -L /net/zwicky/dmtpc/software/lib/ -lTMVA

all: rates efficiency narrowSearch

rates: makeRates.o ../src/AnalysisConfig.o ../src/recoilEnergy.o
	$(CC) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

narrowSearch: narrowSearch.o ../src/AnalysisConfig.o ../src/recoilEnergy.o
	$(CC) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

efficiency: MCefficiency.o ../src/AnalysisConfig.o
	$(CC) $? -o $@  $(LFLAGS)
	mv $@ ../bin/	

%.o: %.cc 
	$(CC) $(CFLAGS) -c $? 

clean:
	rm *.o 
