default:
	mkdir -p bin 
	mkdir -p skim
	make -C src/ -k cleanSkim
	make -C src/ -k wr2analysis
	make -C src/ -k wr3analysis
	make -C src/ -k CutsBasedAnalysis4sh
	make -C src/ -k createSummary
	make -C src/ -k gainMap
	make -C src/ -k multi
%:
	make -C src/ -k $*

clean:
	make -C src/ -k clean
	rm -f bin/*


