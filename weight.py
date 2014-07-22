#!/usr/bin/env python
# Python modules
import argparse
import numpy
#import xml.etree.ElementTree as ET
# Homemade modules
from particles import *
# Argument parser setup
parser = argparse.ArgumentParser()
parser.add_argument("--infile", help="File in ROOT format to be read", nargs='?', type=str, default="ggHH2.root")
parser.add_argument("--intree", help="input TTree name", nargs='?', type=str, default="test")
parser.add_argument("--infilew", help="File in ROOT format to be read (weights)", nargs='?', type=str, default="../BJetRegression/weights.root")
parser.add_argument("--outfile", help="File in ROOT format for output", nargs='?', type=str, default="weighted_2_hbb_SM.root")
parser.add_argument("--outtree", help="output TTree name", nargs='?', type=str, default="test")

args = parser.parse_args()
# ROOT setup
import ROOT
from ROOT import gROOT
from ROOT import TTree, TFile, TLorentzVector, TChain, TVector3
# ROOT initialization
gROOT.Reset()
gROOT.SetBatch()
chain = TChain(args.intree)
chain.Add(args.infile)

weightFile = TFile.Open(args.infilew)
w_ggHH2_hbb_pt = weightFile.Get("hFrac_ggHH2_8TeV_SM_gr_hbb_p4_pt")
w_ggHH_SM_hbb_pt = weightFile.Get("hFrac_ggHH_8TeV_SM_gr_hbb_p4_pt")
#ibin = w_ggHH0_hbb_pt.FindBin(200)
#print w_ggHH0_hbb_pt.GetBinContent(ibin)


########################################################################################
# WRITE THE ROOT FILE
########################################################################################
outfile = TFile(args.outfile, "recreate")
outtree = TTree(args.outtree, args.outtree)

hbb_SM_pt = numpy.zeros(1, dtype=float)
hbb_new_SM_pt = numpy.zeros(1, dtype=float)

outtree.Branch('gr_hbb_SM_w_pt', hbb_SM_pt, 'gr_hbb_SM_w_pt/D')
outtree.Branch('gr_hbb_new_SM_w_pt', hbb_new_SM_pt, 'gr_hbb_new_SM_w_pt/D')

gr_hbb_p4_pt = numpy.zeros(1, dtype=float)

print chain.GetEntries()
for ievt in xrange(chain.GetEntries()):
#    if ievt > 100:
#        break
    chain.GetEntry(ievt)
    b1 = TLorentzVector()
    b2 = TLorentzVector()
    hbb_new = TLorentzVector()
    hbb = TLorentzVector()
    b1.SetPtEtaPhiM(chain.gr_b1_p4_pt, chain.gr_b1_p4_eta, chain.gr_b1_p4_phi, chain.gr_b1_p4_mass)
    b2.SetPtEtaPhiM(chain.gr_b2_p4_pt, chain.gr_b2_p4_eta, chain.gr_b2_p4_phi, chain.gr_b2_p4_mass)
    hbb.SetPtEtaPhiM(chain.gr_hbb_p4_pt, chain.gr_hbb_p4_eta, chain.gr_hbb_p4_phi, chain.gr_hbb_p4_mass)
    
#    b1.Boost(10*b1.Px(), 10*b1.Py(), b1.Pz())
#    b1.Boost(0.5, 0.5, 0)
#    print b1.Px()
#    b2.Boost(0.5, 0.5, 0)
#    print b2.Px()

    ibin = w_ggHH2_hbb_pt.FindBin(chain.gr_hbb_p4_pt)
    w = w_ggHH2_hbb_pt.GetBinContent(ibin)

    ibin_SM = w_ggHH_SM_hbb_pt.FindBin(chain.gr_hbb_p4_pt)
    w_SM = w_ggHH2_hbb_pt.GetBinContent(ibin_SM)
#    w = h1D_pt_data -> GetBinContent(h1D_pt_data -> FindBin(gr_hbb_p4_pt)
#    print 

    spatial_b1 = TVector3(b1.Px()*w*w_SM, b1.Py()*w*w_SM, b1.Pz())
#    b1_new = b1.SetVectM(spatial_b1, b1.M())
    b1.SetVectM(spatial_b1, b1.M())
#    print b1.Px()
    spatial_b2 = TVector3(b2.Px()*w*w_SM, b2.Py()*w*w_SM, b2.Pz())
    b2.SetVectM(spatial_b2, b2.M())

    spatial_h = TVector3(hbb.Px()*w*w_SM, hbb.Py()*w*w_SM, hbb.Pz())
    hbb.SetVectM(spatial_h, hbb.M())
 
    hbb_SM_new = b1 + b2 

    hbb_SM_pt[0] = hbb.Pt()
    hbb_new_SM_pt[0] = hbb_SM_new.Pt()
#    if ievt % 1000 == 0:
#    print hbb_pt[0] - hbb_new_pt[0], hbb_pt[0], hbb_new_pt[0], w, gr_hbb_p4_pt, ibin
    outtree.Fill()
   
# write the tree into the output file and close the file
outfile.Write()
outfile.Close()

