#!usr/bin/python
from optparse import OptionParser
import ROOT
from ROOT import *
from setTDRStyle import setTDRStyle
import math
import sys
import subprocess
import os
import numpy as n
#usage = "usage: %prog [options]"
#parser = optparse.OptionParser(usage)
#parser.add_option("--tag",action="store",type="string",dest="tag",default='PU40bx50')

#(options, args) = parser.parse_args()

#tag = options.tag

#gROOT.Reset()
#setTDRStyle()
#gROOT.ForceStyle()
#gROOT.SetStyle('tdrStyle')


f = TFile.Open("/nfs/data/eepgnnb/EGamma/job_spring14_DYJets_20bx25.root")
t = f.Get("ggNtuplizer/EventTree")
nEntries = t.GetEntries()

can1 = TCanvas("test","test",900,900)
testhisto = TF1("testhisto", "sin(x)", -4, 4)
can1.SetGrid()
testhisto.Draw()
can1.SaveAs("testPlot1.pdf")

h_nTrks = TH1F("h_nTrks", "", 50, 0, 150)
nTrks = []

can2 = TCanvas("test2","test2",900,900)
t.Draw("nTrks>>h_nTrks")
#nTrks.append(h_nTrks)
#nTrks[0].Scale(1/nEntries)
can2.SaveAs("testPlot2.pdf")

can3 = TCanvas("test3","test3",900,900)
gaus = TH1D("gaus", "gaussian fit", 20, -20, 20)
gaus.Fill(gRandom.Gaus())
gaus.Draw()
gausfcn = TF1("gausfcn", "gaus", 0, 20)

can3.SaveAs("testPlot3.pdf")
#print nentries

#for i in xrange(nentries):
#    t.GetEntry(i)
#    h_branchname.Fill(histo)



if __name__ == '__main__':

	usage = "usage: %prog [options]"
	parser = OptionParser()
        parser.add_option("-c", "--check", action="store_false", dest="check", default=True, help="check if files are valid files ")
        (options, args) = parser.parse_args()
#        if len(args) >0:

 #               path = args[0]

	

