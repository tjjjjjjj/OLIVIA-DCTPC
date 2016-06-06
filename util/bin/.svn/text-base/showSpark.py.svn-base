#!/usr/bin/python                                                                                                                             
import sys,os,__main__
tmp_argv=sys.argv
sys.argv=[]
from ROOT import gSystem
from ROOT import gStyle
sys.argv=tmp_argv
gSystem.Load("libMaxCam")
gSystem.Load("libWaveformTools")
gStyle.SetPalette(1)
from ROOT import DmtpcSkimDataset,TCanvas,gErrorIgnoreLevel,TPad
gErrorIgnoreLevel=10000

class TrackViewer:

    SKIM_FILE_NAME="%s/bigdctpc_skim/BigDCTPC_run_%%05dskim.root" % (os.environ['DCTPC_DATA_DIR'])
    RAW_FILE_NAME="%s/bigdctpc_data/BigDCTPC_run_%%05d.root" % (os.environ['DCTPC_DATA_DIR'])
    
    fDataHolder=0
    _canvas=0
    _image=0
    _image_dir='%s/DCTPC_IMAGES' % os.environ['HOME']
    
    _cut_var=(N_TRACK,COUNTS,RANGE,N_WF,SPARK)=xrange(5)
    _cut_name={N_TRACK : "n_track",
               COUNTS  : "counts", 
               RANGE   : "range",
               N_WF    : "n_wf",
               #CHARGE  : "charge",
               SPARK   : "spark"}
    _cut_description={N_TRACK : 'number of tracks in a ccd image',
                      COUNTS  : 'ccd charge (energy variable)',
                      RANGE   : 'length of a reconstructed track',
                      N_WF    : 'number of waveforms from charge readout',
                      #CHARGE  : 'integrated charge from the waveforms (energy variable)',
                      SPARK   : 'binary flag to exclude ccd images flagged as a spark'}
    _cut_value={N_TRACK : 1,
                COUNTS  : 0.,
                RANGE   : 0.,
                N_WF    : 0,
                #CHARGE  : 0.,
                SPARK   : 1}
    _cut_type={N_TRACK  : int,
               COUNTS   : float,
               RANGE    : float,
               N_WF     : int,
               #CHARGE   : float,
               SPARK    : bool}
    _break_program=False

    def __init__(self):
        if not self._canvas:
            self._canvas=TCanvas("c","",600,500)

    def SetImageDir(self,path):
        if not os.path.isdir(path):
            out=os.system('mkdir -p %s' % path)
            if not out==0:
                sys.stderr.write('Cannot create a directory: %s' % path)
                return False
        self._image_dir=path
        return True
                
    def GetCanvas(self):
        return self.c

    def GetSkimFileName(self,run=0):
        return self.SKIM_FILE_NAME % run

    def GetRawFileName(self,run=0):
        return self.RAW_FILE_NAME % run

    def StartLoop(self,r):
        # file existence
        fname1=self.GetSkimFileName(r)
        fname2=self.GetRawFileName(r)            
        for fname in [fname1,fname2]:
            if not os.path.isfile(fname):
                sys.stderr.write('Error: input file missing %s\n' % fname)
                return False

        # loop over events
        self.fDataHolder=DmtpcSkimDataset()
        self.fDataHolder.openRootFile(fname1)
        self.fDataHolder.loadDmtpcEvent(True,fname2)
        if not self.fDataHolder.tree():
            sys.stderr.write("Run %d has no tree saved!\n\n")
            return False

        # Check output directory
        if not os.path.isdir(self._image_dir):
            out=os.system('mkdir %s' % self._image_dir)
            if not out==0:
                sys.stderr.write("Error: cannot make an image dir '%s' ... permission issue?" % self._image_dir)
                sys.exit(1)

        num_events=self.fDataHolder.tree().GetEntries()
        fraction=num_events/100.
        fraction_counter=0
        flusher=[' | ',' \\ ','---',' / ',' | ',' \\ ','---',' / ']
        break_event_loop=False
        for x in xrange(self.fDataHolder.tree().GetEntries()):
            if (x-fraction_counter*fraction)>=0:
                sys.stdout.write('Searching in Run %d ... %3s (%3d%% done)\r' % (self.fDataHolder.event().runNumber(),
                                                                                 flusher[fraction_counter%(len(flusher))],
                                                                                 fraction_counter))
                fraction_counter+=1
                sys.stdout.flush()

            self.fDataHolder.getEvent(x)

            # Apply cuts here
            if self._cut_value[self.N_TRACK] and self._cut_value[self.N_TRACK] > self.fDataHolder.event().ntracks(0):
                continue
            if self._cut_value[self.N_WF] and self._cut_value[self.N_WF] > int(self.fDataHolder.orig_event().scopeData.GetEntries()/4.):
                continue
            show_it = False
            for i in xrange(self.fDataHolder.event().ntracks(0)):
                if self.fDataHolder.event().range(0,i)>self._cut_value[self.RANGE] and self.fDataHolder.event().E(0,x)>self._cut_value[self.COUNTS]:
                    show_it=True
                    break
            if not show_it:
                continue

            print
            print
            print "########################################################################"
            print "#"
            print "# Run: %d ... Event ID: %d " % (self.fDataHolder.event().runNumber(),
                                                    self.fDataHolder.event().eventNumber())
            print '# Found %d waveforms!' % int(self.fDataHolder.orig_event().scopeData().GetEntries()/4)
            self._canvas.Clear()
            self._canvas.cd()
            self._image=self.fDataHolder.event().cluster(0).getImage()

            break_track_loop=False
            for i in xrange(self.fDataHolder.event().ntracks(0)):
                hx=self.fDataHolder.event().x(0,i)
                hy=self.fDataHolder.event().y(0,i)
                length=self.fDataHolder.event().range(0,i)
                if length < self._cut_value[self.RANGE] or self.fDataHolder.event().E(0,x) < self._cut_value[self.COUNTS]:
                    continue
                print "#"
                print '# Event %d ... Track %d: ... (x,y)=(%g,%g) ... length=%g' % (x,i+1,hx,hy,length)
                print "#"
                print "########################################################################"
                self._image.GetXaxis().SetRangeUser(0,1024)
                self._image.GetYaxis().SetRangeUser(0,1024)
                self._image.SetMinimum(-50)
                self._image.SetMaximum(250)
                self._image.Draw("COLZ")
                self._canvas.Update()                                                                     
                
                while 1:
                    sys.stdout.write("\nTrack Image Loop Commands:\n")
                    sys.stdout.write("(*) Hit return ... next track/event.\n")
                    sys.stdout.write("(*) 'n' ... next event.\n")
                    sys.stdout.write("(*) 'r' ... next run.\n")
                    sys.stdout.write("(*) 'c  ... change cut condition.\n")
                    sys.stdout.write("(*) 'q' ... quit the program.\n")
                    sys.stdout.write("(*) 's' ... save the image under %s\n" % self._image_dir)
                    sys.stdout.write(">")
                    sys.stdout.flush()
                    out=sys.stdin.readline().replace('\n','')
                    out=out.split(None)
                    if len(out)==0:
                        break
                    elif out[0].lower() in ['n']:
                        break_track_loop=True
                        break
                    elif out[0].lower() in ['c']:
                        self.ChangeCut()
                        print
                        print '#########################################'
                        print "# Back to track loop command options... #"
                        print '#########################################'
                        print
                        break
                        #try:
                        #    self._range_thres=float(out[1])
                        #except ValueError,IndexError:
                        #    sys.stderr.write("Error: wrong input option usage with 'c'\n")
                    elif out[0].lower() in ['q']:
                        del self.fDataHolder
                        sys.exit(0)
                    elif out[0].lower() in ['r']:
                        break_track_loop=True
                        break_event_loop=True
                        break
                    elif out[0].lower() in ['s']:
                        self._canvas.SaveAs("%s/Image_Run%05d_Event%05d.C" % (self._image_dir,
                                                                              self.fDataHolder.event().runNumber(),
                                                                              x))
                        self._canvas.SaveAs("%s/Image_Run%05d_Event%05d.eps" % (self._image_dir,
                                                                                self.fDataHolder.event().runNumber(),
                                                                                x))
                    else:
                        sys.stderr.write('Error: unrecognized input!\n')
                if break_track_loop:
                    break
            if break_event_loop:
                break
        del self.fDataHolder
        return True

    def ChangeCut(self):
        print
        print "########################################################################"
        sys.stdout.write("# Change Cut Commands:\n")
        sys.stdout.write("# Change cut values by typing info in the fomrat: $VARIABLE $VALUE\n")
        sys.stdout.write("# Here are options for $VARIABLE and curent $VALUE...\n")
        for x in xrange(len(self._cut_var)):
            sys.stdout.write("# (*) %10s = %6s ... %s\n" % (self._cut_name[x],
                                                          self._cut_value[x],
                                                          self._cut_description[x]))
        print "########################################################################"
        sys.stdout.write("\nHit return with no input if you would like to exit this prompt.\n")
        while 1:
            sys.stdout.write("\nYour Input:")
            sys.stdout.flush()
            out=sys.stdin.readline().replace('\n','')
            out=out.split(None)
            if len(out) == 0:
                break
            elif len(out) == 2 and out[0] in self._cut_name.values():
                var_type=0
                for x in xrange(len(self._cut_var)):
                    if out[0]==self._cut_name[x]:
                        var_type=x
                        break
                try:
                    old_value=self._cut_value[var_type]
                    new_value=self._cut_type[var_type](out[1])
                    sys.stdout.write('Changed Variable: %s ... %s => %s\n' % (self._cut_name[var_type],
                                                                            old_value,
                                                                            new_value))
                    self._cut_value[var_type]=new_value
                except ValueError:
                    sys.stderr.write('Error: wrong type value provided for %s\n' % out[0])
            else:
                sys.stderr.write('Wrong number of arguments provided!\n')
        sys.stdout.write('Exiting the command prompt to change cuts...\n\n')
        sys.stdout.flush()

if __name__=='__main__':

    if not len(sys.argv) in [2,3]:
        sys.stderr.write('usage: %s $RUN_START [$RUN_END]\n' % __main__.__file__)
        sys.stderr.write('Error: wrong number of arguments.\n')
        sys.exit(1)

    run_start=-1
    run_end=-1
    try:
        run_start=int(sys.argv[1])
        if len(sys.argv)==3:
            run_end=int(sys.argv[2])
        else:
            run_end=run_start
    except ValueError:
        sys.stderr.write('usage: %s $RUN_START [$RUN_END]\n' % __main__.__file__)
        sys.stderr.write('Error: wrong argument type provided (must be int)\n')
        sys.exit(1)
    if run_start > 0 and run_end >0 and run_end >= run_start:
        k=TrackViewer()
        for x in xrange(run_end-run_start+1):
            k.StartLoop(run_start+x)
        sys.exit(0)
    else:
        sys.stderr.write('Error: unrealistic run number(s) given...\n')
        sys.exit(1)


