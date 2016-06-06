CC=g++

CFLAGS=-O3 -g -I$(TMVADIR)/inc `root-config --cflags` -Wall -Wextra -I../../../MaxCam -I../../MaxCam/waveformtools -I../../../MaxCam/waveformtools/include
LFLAGS= `root-config --libs` -lGeom -lSpectrum -L../../../MaxCam -lMaxCam -L../../../MaxCam/waveformtools/lib -lWaveformTools -L ../../../tmva/lib -lTMVA 


cleanSkim: cleanSkim.o cleanSkimConfig.o cleanSkimProcessor.o
	$(CC) $? -o $@  $(LFLAGS) -lMinuit -lMLP -lXMLIO -lTreePlayer
	mv $@ ../bin/

cleanSkim_neon: cleanSkim_neon.o cleanSkimConfig.o cleanSkimProcessor.o
	$(CC) $? -o $@  $(LFLAGS) -lMinuit -lMLP -lXMLIO -lTreePlayer
	mv $@ ../bin/

quickSkim: quickSkim.o cleanSkimConfig.o cleanSkimProcessor.o
	$(CC) $? -o $@  $(LFLAGS) -lMinuit -lMLP -lXMLIO -lTreePlayer
	mv $@ ../bin/

convert_skim: convert_skim.o cleanSkimConfig.o cleanSkimProcessor.o
	$(CC) $? -o $@  $(LFLAGS) -lMinuit -lMLP -lXMLIO -lTreePlayer
	mv $@ ../bin/


clusterCount: clusterCount.o
	$(CC) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

MCClusterCompare: MCClusterCompare.o
	$(CC) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

compareNClusters: compareNClusters.o
	$(CC) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

multi: multi.o MultiVariate.o MultiVariateResultCint.o
	$(CC) $? -o $@  $(LFLAGS) -lTMVA -lMinuit -lMLP -lXMLIO -lTreePlayer
	mv $@ ../bin/

MultiVariateResultCint.cc: MultiVariate.hh  
	rootcint -f $@ -c $< LinkDef.h

wr2analysis: wr2analysis.cc recoilEnergy.o MultiVariate.o MultiVariateResultCint.o AnalysisCut.o AnalysisCutCint.o AnalysisConfig.o
	$(CC) $(CFLAGS) $^ -o $@  $(LFLAGS) -lTMVA -lMinuit -lMLP -lXMLIO -lTreePlayer	
	mv $@ ../bin/

gainMap: calibTools.o gainMap.o 
	$(CC) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

wr5analysis: wr5analysis.cc recoilEnergy.o MultiVariate.o MultiVariateResultCint.o AnalysisCut.o AnalysisCutCint.o AnalysisConfig.o
	$(CC) $(CFLAGS) $^ -o $@  $(LFLAGS) -lTMVA -lMinuit -lMLP -lXMLIO -lTreePlayer	
	mv $@ ../bin/

4shanalysis: 4shanalysis.cc recoilEnergy.o MultiVariate.o MultiVariateResultCint.o AnalysisCut.o AnalysisCutCint.o AnalysisConfig.o
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS) -lTMVA -lMinuit -lMLP -lXMLIO -lTreePlayer	
	mv $@ ../bin/

raytheonanalysis: raytheonanalysis.cc recoilEnergy.o MultiVariate.o MultiVariateResultCint.o AnalysisCut.o AnalysisCutCint.o AnalysisConfig.o
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS) -lTMVA -lMinuit -lMLP -lXMLIO -lTreePlayer	
	mv $@ ../bin/

AnalysisCutCint.cc: AnalysisCut.hh
	rootcint -f $@ -c $< 

libMultiVariate.so: MultiVariate.o MultiVariateResultCint.o 
	$(CC) $? -o $@ -shared  $(LFLAGS) -lTMVA -lMinuit -lMLP -lXMLIO -lTreePlayer
	mkdir -p ../lib
	mv $@ ../lib/ 
libAnalysis.so: MultiVariate.o MultiVariateResultCint.o AnalysisCut.o AnalysisCutCint.o 
	$(CC) $? -o $@ -shared  $(LFLAGS) -lTMVA -lMinuit -lMLP -lXMLIO -lTreePlayer
	mkdir -p ../lib
	mv $@ ../lib/ 

%.o: %.cc 
	$(CC) $(CFLAGS) -c $? 

clean:
	rm -f *Cint.cc
	rm -f *Cint.h
	rm -f *.o 


