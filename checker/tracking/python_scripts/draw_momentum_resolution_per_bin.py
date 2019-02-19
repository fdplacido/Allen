# author: Vladimir Gligorov (vladimir.gligorov@cern.ch)
# date:   02/2019
#

import os, sys
import argparse
import ROOT
from ROOT import *
from ROOT import gStyle
from ROOT import gROOT
from ROOT import TStyle
from ROOT import gPad

from array import array

from LHCbStyle import *


def getTrackers():
    return ["Upstream"]  #, "Forward"]


nbins = 8

f = ROOT.TFile.Open("momentum_resolution.root", "read")

setLHCbStyle()

trackers = getTrackers()

canvas = ROOT.TCanvas("perbinp", "perbinp", 800, 800)
canvas.Divide(2, 4)
canvas.cd()

paves = {}

for tracker in trackers:
    for i in range(1, nbins + 1):
        canvas.cd(i)
        plot = f.Get(tracker + "/dp_vs_p_py;" + str(i))
        plot.DrawCopy("")
        paves[str(i)] = ROOT.TPaveText(0.65, 0.65, 0.9, 0.9, "NDC")
        paves[str(i)].SetFillStyle(0)
        paves[str(i)].SetFillColor(0)
        paves[str(i)].SetBorderSize(0)
        paves[str(i)].AddText("Bin " + str(i))
        paves[str(i)].Draw()

canvas.SaveAs("../../../plotsfornote/" + tracker + "MomResByBinOfP.pdf")
