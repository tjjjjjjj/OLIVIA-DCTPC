#!/usr/bin/python                                                                                                                             
import sys,os,__main__
tmp_argv=sys.argv
sys.argv=[]
from ROOT import gSystem
sys.argv=tmp_argv
gSystem.Load("libMaxCam")
gSystem.Load("libWaveformTools")
from ROOT import DmtpcSkimDataset,TCanvas,gErrorIgnoreLevel,TPad
gErrorIgnoreLevel=10000

class EventViewer:

    LOG='%s/event_list.txt' % os.environ['DCTPC_DATA_DIR']
    SKIM_FILE_NAME="%s/skim/dmtpc_DC_%%05dskim.root" % (os.environ['DCTPC_DATA_DIR'])
    RAW_FILE_NAME="%s/data/dmtpc_DC_%%05d.root" % (os.environ['DCTPC_DATA_DIR'])
    
    fEventList={}
    fReadEventList={}
    fDataHolder=0

    _canvas=0
    _image=0

    _histoAnode=0
    _histoCanvas=0

    def __init__(self):
        if not self._canvas:
            self._canvas=TCanvas("c","",600,500)
            self._histoCanvas=TCanvas('c_Histo','',600,500)
            
    def GetCanvas(self):
        return self.c

    def GetSkimFileName(self,run=0):
        return self.SKIM_FILE_NAME % run

    def GetRawFileName(self,run=0):
        return self.RAW_FILE_NAME % run

    def LoadEventList(self):
        self.fEventList={}
        if os.path.isfile(self.LOG):
            contents=open(self.LOG,'r').read().split('\n')
            for line in contents:
                if not len(line.split(None))==2:
                    continue
                run_id=int(line.split(None)[0])
                event_id=int(line.split(None)[1])
                if not run_id in self.fEventList.keys():
                    self.fEventList[run_id]=[]
                self.fEventList[run_id].append(event_id)
        else:
            sys.stderr.write('Error: file not found %s' % self.LOG)
            return False
        if len(self.fEventList)==0:
            sys.stderr.write('Error: no events found in the event list %s' % self.LOG)
            return False

        return True

    def StartLoop(self,run_start=1641,run_end=999999):
        if not len(self.fEventList) and not self.LoadEventList():
            return False

        runs=self.fEventList.keys()
        runs.sort()
        status=True
        for r in runs:
            if r < run_start or r>run_end:
                continue
            # file existence
            fname1=self.GetSkimFileName(r)
            fname2=self.GetRawFileName(r)            
            for fname in [fname1,fname2]:
                if not os.path.isfile(fname):
                    sys.stderr.write('Error: input file missing %s\n' % fname)
                    status=False
            if not status:
                continue

            # loop over events
            events=self.fEventList[r]
            for event in events:
                if r in self.fReadEventList and event in self.fReadEventList[r]:
                    continue
                else:
                    if not r in self.fReadEventList:
                        self.fReadEventList[r]=[]
                    self.fReadEventList[r].append(event)

                if not self.fDataHolder:
                    self.fDataHolder=DmtpcSkimDataset()
                    self.fDataHolder.openRootFile(fname1)
                    self.fDataHolder.loadDmtpcEvent(True,fname2)

                self.fDataHolder.getEvent(event)
                print
                print 'Run: %d ... Event ID: %d' % (self.fDataHolder.event().runNumber(),
                                                    self.fDataHolder.event().eventNumber())
                self._canvas.Clear()
                self._histoCanvas.Clear()
                for i in xrange(self.fDataHolder.event().ntracks(0)):
                    print 'Track %d: ... (x,y)=(%g,%g) ... length=%g' % ((i+1),
                                                                         self.fDataHolder.event().x(0,i),
                                                                         self.fDataHolder.event().y(0,i),
                                                                         self.fDataHolder.event().range(0,i))
                print 'Found %d waveforms!' % int(self.fDataHolder.orig_event().scopeData().GetEntries()/4)
                print "Would you like to see image?"
                print "Type 'i' for ccd image, 'w' for waveform, 'b' for both, enter to move next event:",
                sys.stdout.flush()
                out=sys.stdin.readline().replace('\n','')
                if out=='i' or i=='b':
                    self._canvas.cd()
                    self._image=self.fDataHolder.event().cluster(0).getImage().Draw("COLZ")
                    self._canvas.Update()

                if out=='w' or i=='b':
                    for x in xrange(self.fDataHolder.orig_event().scopeData().GetEntries()/4):
                        print 'Displaying Anode Waveform %d ... (hit enter to move next)\r' % (x+1)
                        sys.stdout.flush()
                        self._histoAnode=self.fDataHolder.orig_event().scopeData(x,0,1)
                        self._histoCanvas.cd()
                        self._histoAnode.Draw()
                        self._histoCanvas.Update()
                        out = sys.stdin.readline()
                        if out.replace('\n','').lower() in ['n','no','exit','quit']:
                            status=False
                            break
                if not status or out.lower() in ['n','no','exit','quit']:
                    status=False
                    break

            del self.fDataHolder
            self.fDataHolder=0
            if not status:
                break
        return status

    def GetWFCanvas(self,c):
        pAnode=TPad("pAnode_%s" % c.GetName(),"",0.02,0.67,0.98,0.99)
        pMesh=TPad("pMesh_%s" % c.GetName(),"",0.02,0.34,0.98,0.66)
        pVeto=TPad("pVeto_%s" % c.GetName(),"",0.02,0.01,0.98,0.33)
        c.cd()
        pAnode.Draw()
        pMesh.Draw()
        pVeto.Draw()
        return [c,pAnode,pMesh,pVeto]
                    

if __name__=='__main__':
    k=EventViewer()
    k.StartLoop()


