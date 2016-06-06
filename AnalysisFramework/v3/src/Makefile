CC=g++
CFLAGS=-Os -g `root-config --cflags`
LFLAGS=`root-config --libs` -lSpectrum -L../../../MaxCam -lMaxCam -L../../../MaxCam/waveformtools/lib -lWaveformTools

cleanSkim.o: cleanSkim.cc
	$(CC) $(CFLAGS) -c $? 

compareSkim.o: compareSkim.cc
	$(CC) $(CFLAGS) -c $? 

compareSkim: compareSkim.o 
	$(CC) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

gainMap.o: gainMap.cc
	$(CC) $(CFLAGS) -c $? 

contourPicker.o: contourPicker.cc
	$(CC) $(CFLAGS) -c $? 

createSummary.o: createSummary.cc
	$(CC) $(CFLAGS) -c $?

createSummary: createSummary.o
	$(CC) $? -o $@ $(LFLAGS)
	mv $@ ../bin/

cleanSkim: cleanSkim.o cleanSkimConfig.o 
	$(CC) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

gainMap: calibTools.o gainMap.o 
	$(CC) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

isolateNtuple: isolateNtuple.cc
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

combineRates: combineRates.cc
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

combinePass: combinePass.cc recoilEnergy.o
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

dumpPass: dumpPass.cc style.o
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

dumpRates: dumpRates.cc
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

ntupledump: ntupledump.cc
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

efficiency: efficiency.cc recoilEnergy.o
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

cleanSkimV2Equiv: cleanSkimV2Equiv.cc  
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS)
	mv $@ ../bin/

simpleConfig.o: simpleConfig.cc
	$(CC) $(CFLAGS) -c $? 	

cleanSkimConfig.o: cleanSkimConfig.cc 
	$(CC) $(CFLAGS) -c $? 

treeTest: treeTest.cc
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS)	
	mv $@ ../bin/

cleanSkimUtils: cleanSkimUtils.cc style.o
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS)	
	mv $@ ../bin/

testDataset: testDataset.cc 
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS)	
	mv $@ ../bin/

mergeDiet: mergeDiet.cc 
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS)	
	mv $@ ../bin/

testSkimDataset: testSkimDataset.cc 
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS)	
	mv $@ ../bin/

wr2analysis: wr2analysis.cc recoilEnergy.o AnalysisResult.o AnalysisResultCint.o
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS)	
	mv $@ ../bin/

wr3analysis: wr3analysis.cc recoilEnergy.o AnalysisResult.o MultiVariate.o AnalysisResultCint.o MultiVariateResultCint.o
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS) -lTMVA -lMinuit -lMLP -lXMLIO -lTreePlayer	
	mv $@ ../bin/

wr4analysis: wr4analysis.cc recoilEnergy.o AnalysisResult.o MultiVariate.o AnalysisResultCint.o MultiVariateResultCint.o AnalysisCut.o AnalysisCutCint.o
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS) -lTMVA -lMinuit -lMLP -lXMLIO -lTreePlayer	
	mv $@ ../bin/

CutsBasedAnalysis4sh: CutsBasedAnalysis4sh.cc recoilEnergy.o AnalysisResult.o MultiVariate.o AnalysisResultCint.o MultiVariateResultCint.o AnalysisCut.o AnalysisCutCint.o
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS) -lTMVA -lMinuit -lMLP -lXMLIO -lTreePlayer	
	mv $@ ../bin/

cutOptimization: cutOptimization.cc recoilEnergy.o AnalysisResult.o AnalysisResultCint.o style.o contourPicker.o
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS)	
	mv $@ ../bin/

paperPlots: paperPlots.cc recoilEnergy.o AnalysisResult.o AnalysisResultCint.o style.o
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS)	
	mv $@ ../bin/

mc_compare: mc_compare.cc recoilEnergy.o AnalysisResult.o AnalysisResultCint.o style.o
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS)	
	mv $@ ../bin/

makeContours: makeContours.cc recoilEnergy.o AnalysisResult.o AnalysisResultCint.o style.o
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS)	
	mv $@ ../bin/
result_read: result_read.cc AnalysisResult.o AnalysisResultCint.o
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS)	
	mv $@ ../bin/

testConfig: testConfig.cc cleanSkimConfig.o 
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS)	
	mv $@ ../bin/

makePlots: makePlots.cc recoilEnergy.o style.o
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS)	
	mv $@ ../bin/

recoilEnergy.o: recoilEnergy.cc
	$(CC) $(CFLAGS) -c $? 	

calibTools.o: calibTools.cc
	$(CC) $(CFLAGS) -c $? 	

style.o: style.cc
	$(CC) $(CFLAGS) -c $? 	

clean: 
	rm -f *.o *Cint.h *Cint.cc

AnalysisResult.o: AnalysisResult.cc  
	$(CC) $(CFLAGS) -c $? 

AnalysisResultCint.o: AnalysisResultCint.cc
	$(CC) $(CFLAGS) -c $? 

AnalysisResultCint.cc: AnalysisResult.hh  
	rootcint -f $@ -c $<

multi: multi.cc MultiVariate.o MultiVariateResultCint.o style.o
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS) -lTMVA -lMinuit -lMLP -lXMLIO -lTreePlayer
	mv $@ ../bin/

MultiVariate.o: MultiVariate.cc 
	$(CC) $(CFLAGS) -c $?   

MultiVariateResultCint.cc: MultiVariate.hh  
	rootcint -f $@ -c $< LinkDef.h

MultiVariateResultCint.o: MultiVariateResultCint.cc  
	$(CC) $(CFLAGS) -c $? 

AnalysisCut.o: AnalysisCut.cc
	$(CC) $(CFLAGS) -c $?

AnalysisCutCint.cc: AnalysisCut.hh
	rootcint -f $@ -c $< 

AnalysisCutCint.o: AnalysisCutCint.cc
	$(CC) $(CFLAGS) -c $?

testcut: testcut.cc AnalysisCut.o AnalysisCutCint.o
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS) -lTMVA -lMinuit -lMLP -lXMLIO -lTreePlayer
	mv $@ ../bin/
