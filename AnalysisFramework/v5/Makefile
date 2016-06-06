default:
	mkdir -p bin 
	mkdir -p skim
	make -C src/ -k cleanSkim

%:
	make -C src/ -k $*

clean:
	make -C src/ -k clean
	rm -f bin/*


