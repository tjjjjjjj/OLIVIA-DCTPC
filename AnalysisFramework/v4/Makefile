default:
	mkdir -p bin 
	mkdir -p skim
	make -C  src/ -k cleanSkim
#	make -C src/ -k gainMap
#	make -C src/ -k wr5analysis
	make -C src/ -k 4shanalysis
	make -C src/ -k raytheonanalysis
	make -C src/ -k quickSkim

%:
	make -C src/ -k $*

clean:
	make -C src/ -k clean
	rm -f bin/*


