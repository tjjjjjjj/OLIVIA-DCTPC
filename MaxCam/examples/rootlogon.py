#! /usr/bin/python

#Remember to include Root directories in $PYTHONPATH
#Run in Python as execfile("rootlogon.py")

#Load Basic Python Classes
from datetime import *
from array import array
from math import *
from cmath import *

#Load User Libraries
from ROOT import gSystem
gSystem.Load("../MaxCam/libMaxCam.so")
gSystem.Setenv("MCTABLES","../MaxCam/tables");

#Load ROOT Classes
import ROOT
#Load important ROOT classes separately
from ROOT import gROOT,gStyle, gPad, gRandom
from ROOT import TH1F, TH2F, TGraph, TGraph2D, TProfile
from ROOT import TF1, TF2, TString
from ROOT import TTree, TChain, TFile
from ROOT import TRandom3, TMath, TCanvas, TObject, TObjArray

#Load User Classes
from ROOT import DmtpcDataset, DmtpcEvent
from ROOT import DmtpcSkimDataset, DmtpcSkimEvent


gROOT.SetStyle("Pub")

#Reset color palette
stops = array('d',[0,0.34,0.61,0.84,1.0])
red = array('d',[0.0,0.0,0.87,1.0,0.51])
green = array('d',[0.0,0.81,1.00,0.2,0.0])
blue = array('d',[0.51,1.0,0.12,0.0,0.0])
ROOT.TColor.CreateGradientColorTable(5,stops,red,green,blue,255)
gStyle.SetNumberContours(255)
stops = None
red = None
green = None
blue = None

gStyle.SetMarkerStyle(20)
gStyle.SetMarkerSize(0.5)
gStyle.SetTextSize(0.04)
gStyle.SetLabelSize(0.03,"x")
gStyle.SetTitleSize(0.04,"x")
gStyle.SetLabelSize(0.03,"y")
gStyle.SetTitleSize(0.04,"y")
gStyle.SetLabelSize(0.03,"z")
gStyle.SetTitleSize(0.04,"z")
gStyle.SetTitleOffset(1,"y")
gStyle.SetFuncColor(ROOT.kBlack)
gStyle.SetHistLineColor(ROOT.kBlack)

