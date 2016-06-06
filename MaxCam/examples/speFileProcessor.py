# program to convert an SPE file into a ROOT file
# for analysis with the DMTPC analysis framework.
# SPE is the proprietary format of Princeton Instruments
import sys
import os
import piUtils # requires pyfits and numpy
import ROOT
import platform

def determineMaxCamSoFile():

  # get the naming correct for linux vs. mac
  osname = os.uname()[0]
  if osname == 'Darwin':
    soFile = "MaxCam_macosx.so"
  elif osname == 'Linux':
    soFile = "MaxCam_linux.so"
  else:
    print 'unknown osname', osname
    print "exiting"
  sys.exit()
  return "../"+soFile

def convertSPEToROOT(speFile="test.SPE", outroot="test.root"):
  ROOT.gROOT.Reset()
  #soFile = "MaxCam_linux.so"
  #soFile = "MaxCam_macosx.so"
  #soFilePath = "../"+soFile
  soFilePath = determineMaxCamSoFile()
  ROOT.gSystem.Load(soFilePath)
  from ROOT import DmtpcDataset, MaxCamChannel
  from ROOT import MaxCamConfig, MaxCamImageTools
  
  speFile = "test.SPE"
  speDict = piUtils.readSpe(speFile)
  
  # should read this from the SPE data!
  npix  = 512
  nbins = 512
  
  # output root filename
  outroot = "test.root"
  darkFitsFile = "dark.fits"
  
  dataset = DmtpcDataset()
  dataset.createRootFile(outroot, "recreate")
  
  # setup data file
  print "setup data file"
  dataset.setComment("COMMENT");
  dataset.setLocation("LOCATION");
  dataset.setKeyword("KEYWORD");
  dataset.comment().Write();
  dataset.keyword().Write();
  dataset.location().Write();

  # set up experimentConfig values
  print "Set up experimentConfig" 
  driftHV  = MaxCamChannel("driftHV",    "Drift Voltage", -1,  -1);
  anodeHV  = MaxCamChannel("anodeHV",    "Anode Voltage", -1,  -1);
  pressure = MaxCamChannel("pressure",   "Gas pressure",  -1,  -1);    
  driftHV.currentValue  = driftHV.setValue  = 1.5;
  anodeHV.currentValue  = anodeHV.setValue  = 0.68;
  pressure.currentValue = pressure.setValue = 75;
  
  # set up MaxCamConfig file for the dataset
  print "Set up MaxCamConfig" 
  
  ccdconfig = MaxCamConfig("ccdConfig","CCD Camera configuration");
  ccdconfig.cameraID = 0;
  ccdconfig.row_width = npix;
  ccdconfig.img_rows = npix;
  ccdconfig.hbin = npix / nbins;
  ccdconfig.vbin = npix / nbins;
  ccdconfig.ul_x = 0;
  ccdconfig.ul_y = 0;
  ccdconfig.lr_x = npix;
  ccdconfig.lr_y = npix;
  ccdconfig.CCDTemp = -20;
  ccdconfig.CCDTempSet = -20;
  ccdconfig.exposureTime = 1000;
  ccdconfig.frameType = 0;
  ccdconfig.nFlushes = -1;
  ccdconfig.bitDepth = 65535;
  ccdconfig.daqTime = -1;
  
  ##################
  ### dark frame ###
  ##################
  bias = MaxCamImageTools.convertFitsIntoROOT(darkFitsFile,"biasFrame1", 0);
  bias.Write();

  ########################################
  ### Save CCD images and fill dataset ###
  ########################################
  nframes = len(speDict['data'])
  print "nframes = ", nframes

  cameraNumber = 0
  tempNameRoot = "temp"
  nx = len(speDict['data'][0])
  ny = len(speDict['data'][0][0])
  print "nx, ny = ", nx, ", ", ny
  
  #ccdData = dataset.event().ccdData()

  # see DmtpcDAQ::beforeEvent()
  dataset.event().setRunNumber(1)
  
  # need this for TH2F, TClonesArray filling to work...
  # see PyROOT forum on ROOTTalk for more info
  ROOT.TH2F.AddDirectory(False)
  
  for imgnum in range(nframes):
  #for imgnum in range(2):
    # the first image is empty for the ProEM...
    #if imgnum == 0:
    #  continue

    ## events start from 1 (not zero)
    dataset.event().increaseEventNumber()
    
    tempName = tempNameRoot+str(dataset.event().eventNumber())+str(imgnum)
    if imgnum % 50 == 0:
      print "on image number: ", imgnum
    
    img = ROOT.TH2F(tempName, tempName, nx, 0, nx, ny, 0, ny)
    for ii in range(nx):
      for jj in range(ny):
        img.SetBinContent(jj+1, ii+1, speDict['data'][imgnum][ii][jj])

    dataset.event().ccdData()[cameraNumber] = img

    dataset.fill()

  print "done with loop"
  print "dataset.write()"
  dataset.write()  
  print "dataset.file().Close()"
  dataset.file().Close()  
  
  ROOT.TH2F.AddDirectory(True)

def help():
  print "syntax:"
  print "python %s <spefile> <rootfile>" % (sys.argv[0])
  print "  converts an SPE file to ROOT"
  print "  SPE files are produced by Princeton Instruments cameras"
  print "  spefile = name of SPE file to be converted"
  print "  rootfile = name of output root file"
  print ""
  print "  e.g."
  print "  > python %s test.SPE test.root" % (sys.argv[0])
  
if __name__=="__main__":
  if len(sys.argv) < 3:
    help()
    sys.exit()
  convertSPEToROOT(sys.argv[1], sys.argv[2])
  
