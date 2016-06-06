CC=g++

CFLAGS=-O0 -g `root-config --cflags` -Wall -fstack-protector-all -I../../../MaxCam -I../../../MaxCam/waveformtools -I../../../MaxCam/waveformtools/include -I../v1/code
LFLAGS=`root-config --libs` -lSpectrum -lMathMore -L../../../MaxCam -lMaxCam -L../../../MaxCam/waveformtools/lib -lWaveformTools -L../v1/lib -lDmtpcSim

all: MCCosMix MCGaus MCCutoff

MCCosMix: MCCosMix.o
	$(CC) $? -o $@  $(LFLAGS)

MCCutoff: MCCutoff.o
	$(CC) $? -o $@  $(LFLAGS)

MCGaus: MCGaus.o
	$(CC) $? -o $@  $(LFLAGS)

%.o: %.cc 
	$(CC) $(CFLAGS) -c $? 

clean:
	rm *.o 
