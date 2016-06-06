import sys
import speFileProcessor

spefilelist = sys.argv[1]
ff = open(spefilelist, 'r')
spefiles = ff.readlines()

for spefile in spefiles:
  spefile = spefile.strip()
  rootfile = spefile.replace(".SPE", ".root")
  print spefile, " ", rootfile
  speFileProcessor.convertSPEToROOT(spefile, rootfile)
