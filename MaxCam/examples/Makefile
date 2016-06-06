include	../Makefile.root
CXX = g++
LDFLAGS = -L../ -lMaxCam 

LDFLAGS  += /opt/apogee/lib/apogee_USB.so -L/opt/apogee/lib -ltcl8.3 -lccd -lfitsTcl -lusb
CXXFLAGS += -I/opt/apogee/include -I/opt/apogee/src/apogee -DLINUX


altatest: altatest.o
	$(LD)  $(LDFLAGS) ../MaxCamAlta.o altatest.o $(ROOTLIBS) -o altatest

serial: testSerial.o
	$(CXX) testSerial.cc -o serial

makeDecayChain: makeDecayChain.cc
	$(CC) -g makeDecayChain.cc `root-config --cflags` -o makeDecayChain  -L ../ -lMaxCam `root-config --glibs` -L../waveformtools/lib -lWaveformTools

teststitch: teststitch.cc
	$(CC) -g teststitch.cc `root-config --cflags` -o teststitch  -L ../ -lMaxCam `root-config --glibs` -L../waveformtools/lib -lWaveformTools

stitchhists: stitchhists.cc
	$(CC) -g stitchhists.cc `root-config --cflags` -o stitchhists  -L ../ -lMaxCam `root-config --glibs` -L../waveformtools/lib -lWaveformTools

sumDataset: sumDataset.cc
	$(CC) -g sumDataset.cc `root-config --cflags` -o sumDataset  -L ../ -lMaxCam `root-config --glibs` -L../waveformtools/lib -lWaveformTools

pixdist: pixdist.cc
	$(CC) -g pixdist.cc `root-config --cflags` -o pixdist  -L ../ -lMaxCam `root-config --glibs` -L../waveformtools/lib -lWaveformTools

testConv: testConv.cc
	$(CC) -g testConv.cc `root-config --cflags` -o testConv  -L ../ -lMaxCam `root-config --glibs` -L../waveformtools/lib -lWaveformTools

ntriggers: ntriggers.cc
	$(CC) -g ntriggers.cc `root-config --cflags` -o ntriggers  -L ../ -lMaxCam `root-config --glibs` -L../waveformtools/lib -lWaveformTools

biasHotPix: biasHotPix.cc
	$(CC) -g biasHotPix.cc `root-config --cflags` -o biasHotPix  -L ../ -lMaxCam `root-config --glibs` -L../waveformtools/lib -lWaveformTools

teststitch_stored: teststitch_stored.cc
	$(CC) -g teststitch_stored.cc `root-config --cflags` -o teststitch_stored  -L ../ -lMaxCam `root-config --glibs` -L../waveformtools/lib -lWaveformTools

teststitch_image: teststitch_image.cc
	$(CC) -g teststitch_image.cc `root-config --cflags` -o teststitch_image  -L ../ -lMaxCam `root-config --glibs` -L../waveformtools/lib -lWaveformTools -lGeom

testStack: testStack.cc
	$(CC) -g testStack.cc `root-config --cflags` -o testStack  -L ../ -lMaxCam `root-config --glibs` -L../waveformtools/lib -lWaveformTools

compareBias: compareBias.cc
	$(CC) -g compareBias.cc `root-config --cflags` -o compareBias  -L ../ -lMaxCam `root-config --glibs` -L../waveformtools/lib -lWaveformTools


gaussianFit: gaussianFit.cc
	$(CC) -g -O3 -fstack-protector gaussianFit.cc `root-config --cflags` -o gaussianFit  -I../ -L ../ -lMaxCam `root-config --glibs` -lMathMore -L../waveformtools/lib -lWaveformTools

rationalFit: rationalFit.cc
	$(CC) -g -fstack-protector rationalFit.cc `root-config --cflags` -o rationalFit  -I../ -L ../ -lMaxCam `root-config --glibs` -lMathMore -L../waveformtools/lib -lWaveformTools

doSRIMFit: doSRIMFit.cc
	$(CC) -g -fstack-protector doSRIMFit.cc `root-config --cflags` -o doSRIMFit  -I../ -L ../ -lMaxCam `root-config --glibs` -lMathMore -L../waveformtools/lib -lWaveformTools

doSRIMLineFit2D: doSRIMLineFit2D.cc
	$(CC) -g -fstack-protector doSRIMLineFit2D.cc `root-config --cflags` -o doSRIMLineFit2D  -I../ -L ../ -lMaxCam `root-config --glibs` -lMathMore -L../waveformtools/lib -lWaveformTools
