#!/usr/bin/env python
# Python modules
import argparse
import numpy
import math
from os import path
#import xml.etree.ElementTree as ET
# Homemade modules
from particles import *
# Argument parser setup

afs_plottree = "/afs/cern.ch/user/o/obondu/public/forAkanksha/2014-06-25_selection_noRegression_noMassCut_v15"

parser = argparse.ArgumentParser()
parser.add_argument("--infile", help="File in ROOT format to be read", nargs='?', type=str, default="ggHH_8TeV_noRegression_noMassCut_v15.root")

parser.add_argument("--intree", help="input TTree name", nargs='?', type=str, default="ggHH_8TeV")
parser.add_argument("--infilew", help="File in ROOT format to be read (weights)", nargs='?', type=str, default="../BJetRegression/weights_2D.root")
parser.add_argument("--outfile", help="File in ROOT format for output", nargs='?', type=str, default="weighted_reco_by_lambda2.root")
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
chain.Add(path.join(afs_plottree, args.infile))


weightFile = TFile.Open(args.infilew)
w_hbb_pt_costhetastar_CS = weightFile.Get("hFrac_ggHH2_8TeV_gr_costhetastar_CS_pt_2D")
#ibin = w_ggHH0_hbb_pt.FindBin(200)
#print w_ggHH0_hbb_pt.GetBinContent(ibin)


########################################################################################
# WRITE THE ROOT FILE
########################################################################################
outfile = TFile(args.outfile, "recreate")
outtree = TTree(args.outtree, args.outtree)

# Prepare the output variables
reco_w_pho1_pt = numpy.zeros(1, dtype=float)
reco_w_pho1_eta = numpy.zeros(1, dtype=float)
reco_w_pho1_phi = numpy.zeros(1, dtype=float)
reco_w_pho1_mass = numpy.zeros(1, dtype=float)
reco_w_pho2_pt = numpy.zeros(1, dtype=float)
reco_w_pho2_eta = numpy.zeros(1, dtype=float)
reco_w_pho2_phi = numpy.zeros(1, dtype=float)
reco_w_pho2_mass = numpy.zeros(1, dtype=float)

reco_w_jet1_pt = numpy.zeros(1, dtype=float)
reco_w_jet1_eta = numpy.zeros(1, dtype=float)
reco_w_jet1_phi = numpy.zeros(1, dtype=float)
reco_w_jet1_mass = numpy.zeros(1, dtype=float)
reco_w_jet2_pt = numpy.zeros(1, dtype=float)
reco_w_jet2_eta = numpy.zeros(1, dtype=float)
reco_w_jet2_phi = numpy.zeros(1, dtype=float)
reco_w_jet2_mass = numpy.zeros(1, dtype=float)

reco_w_gg_pt = numpy.zeros(1, dtype=float)
reco_w_gg_eta = numpy.zeros(1, dtype=float)
reco_w_gg_phi = numpy.zeros(1, dtype=float)
reco_w_gg_mass = numpy.zeros(1, dtype=float)

reco_w_jj_pt = numpy.zeros(1, dtype=float)
reco_w_jj_eta = numpy.zeros(1, dtype=float)
reco_w_jj_phi = numpy.zeros(1, dtype=float)
reco_w_jj_mass = numpy.zeros(1, dtype=float)

weight_lambda = numpy.zeros(1, dtype=float)
costhetastar_CS = numpy.zeros(1, dtype= float)
#dEta_gg_jj = numpy.zeros(1, dtype= float)
#dPhi_gg_jj = numpy.zeros(1, dtype= float)

# Create the branch to store the variable in the tree
outtree.Branch('pho1_pt', reco_w_pho1_pt, 'pho1_pt/D')
outtree.Branch('pho1_eta', reco_w_pho1_eta, 'pho1_eta/D')
outtree.Branch('pho1_phi', reco_w_pho1_phi, 'pho1_phi/D')
outtree.Branch('pho1_mass', reco_w_pho1_mass, 'pho1_mass/D')
outtree.Branch('pho2_pt', reco_w_pho2_pt, 'pho2_pt/D')
outtree.Branch('pho2_eta', reco_w_pho2_eta, 'pho2_eta/D')
outtree.Branch('pho2_phi', reco_w_pho2_phi, 'pho2_phi/D')
outtree.Branch('pho2_mass', reco_w_pho2_mass, 'pho2_mass/D')

outtree.Branch('jet1_pt', reco_w_jet1_pt, 'jet1_pt/D')
outtree.Branch('jet1_eta', reco_w_jet1_eta, 'jet1_eta/D')
outtree.Branch('jet1_phi', reco_w_jet1_phi, 'jet1_phi/D')
outtree.Branch('jet1_mass', reco_w_jet1_mass, 'jet1_mass/D')
outtree.Branch('jet2_pt', reco_w_jet2_pt, 'jet2_pt/D')
outtree.Branch('jet2_eta', reco_w_jet2_eta, 'jet2_eta/D')
outtree.Branch('jet2_phi', reco_w_jet2_phi, 'jet2_phi/D')
outtree.Branch('jet2_mass', reco_w_jet2_mass, 'jet2_mass/D')

outtree.Branch('gg_pt', reco_w_gg_pt, 'gg_pt/D')
outtree.Branch('gg_eta', reco_w_gg_eta, 'gg_eta/D')
outtree.Branch('gg_phi', reco_w_gg_phi, 'gg_phi/D')
outtree.Branch('gg_mass', reco_w_gg_mass, 'gg_mass/D')

outtree.Branch('jj_pt', reco_w_jj_pt, 'jj_pt/D')
outtree.Branch('jj_eta', reco_w_jj_eta, 'jj_eta/D')
outtree.Branch('jj_phi', reco_w_jj_phi, 'jj_phi/D')
outtree.Branch('jj_mass', reco_w_jj_mass, 'jj_mass/D')

outtree.Branch('weight_lambda', weight_lambda, 'weight_lambda/D')
outtree.Branch('costhetastar_CS', costhetastar_CS, 'costhetastar_CS/D')
#outtree.Branch('dEta_gg_jj', dEta_gg_jj, 'dEta_gg_jj/D')
#outtree.Branch('dPhi_gg_jj', dPhi_gg_jj, 'dPhi_gg_jj/D')

print chain.GetEntries()
for ievt in xrange(chain.GetEntries()):
#    if ievt > 100:
#        break
    chain.GetEntry(ievt)
    p1 = TLorentzVector()
    p2 = TLorentzVector()
    pho1 = TLorentzVector()
    pho2 = TLorentzVector()
    jet1 = TLorentzVector()
    jet2 = TLorentzVector()
    gg = TLorentzVector()
    jj = TLorentzVector()
    hh = TLorentzVector()
    hbb_CS = TLorentzVector()
    pho1.SetPtEtaPhiM(chain.pho1_pt, chain.pho1_eta, chain.pho1_phi, chain.pho1_mass)
    pho2.SetPtEtaPhiM(chain.pho2_pt, chain.pho2_eta, chain.pho2_phi, chain.pho2_mass)
    jet1.SetPtEtaPhiM(chain.jet1_pt, chain.jet1_eta, chain.jet1_phi, chain.jet1_mass)
    jet2.SetPtEtaPhiM(chain.jet2_pt, chain.jet2_eta, chain.jet2_phi, chain.jet2_mass)
    jj.SetPtEtaPhiM(chain.jj_pt, chain.jj_eta, chain.jj_phi, chain.jj_mass)
    gg.SetPtEtaPhiM(chain.gg_pt, chain.gg_eta, chain.gg_phi, chain.gg_mass)
    hbb_CS.SetPtEtaPhiM(chain.gr_hbb_p4_pt, chain.gr_hbb_p4_eta, chain.gr_hbb_p4_phi, chain.gr_hbb_p4_mass)
    
#    b1.Boost(10*b1.Px(), 10*b1.Py(), b1.Pz())
#    b1.Boost(0.5, 0.5, 0)
#    print b1.Px()
#    b2.Boost(0.5, 0.5, 0)
#    print b2.Px()
    p1.SetPxPyPzE(0,0,4000,4000)
    p2.SetPxPyPzE(0,0,-4000,4000)
#    print p1.Px()
    hh = jj + gg
#    px = hbb.Px() + hgg.Px()
#    py = hbb.Py() + hgg.Py()
#    pz = hbb.Pz() + hgg.Pz()
#    e = hbb.E() + hgg.E()
    boost = -hh.BoostVector()
    x = boost.X()
    y = boost.Y()
    z = boost.Z()
    p1.Boost(x,y,z)
    p2.Boost(x,y,z)
#    print hbb_CS.Px(), hbb.Px(), hbb_CS.Py(), hbb.Py()
    hbb_CS.Boost(x,y,z)
#    print hbb_CS.Px(), hbb.Px(), hbb_CS.Py(), hbb.Py()

#    print 4+5
#    print p1.Px(), math.cos(0.1)
#    hbb.Boost(px/e, py/e, pz/e)
#    hgg.Boost(px/e, py/e, pz/e)
    CS_axis = p2.Vect() - p1.Vect()
    costhetastar = abs(math.cos(CS_axis.Angle(hbb_CS.Vect())))        
#    print abs(3), abs(-3)
    d_eta = jj.Eta() - gg.Eta()
    d_phi = jj.Phi() - gg.Phi()

    ibin = w_hbb_pt_costhetastar_CS.FindBin(abs(chain.gr_hbbhgg_costhetastar_CS), chain.gr_hbb_p4_pt)
    w = w_hbb_pt_costhetastar_CS.GetBinContent(ibin)


#FIXME    print  chain.gr_hgg_p4_pt, chain.gg_pt
#FIXME    w = w * chain.gr_hgg_p4_pt / chain.gg_pt
 
#    spatial_pho1 = TVector3(pho1.Px()*w, pho1.Py()*w, pho1.Pz())
#    pho1.SetVectM(spatial_pho1, pho1.M())

#    spatial_pho2 = TVector3(pho2.Px()*w, pho2.Py()*w, pho2.Pz())
#    pho2.SetVectM(spatial_pho2, pho2.M())

#    spatial_gg = TVector3(gg.Px()*w, gg.Py()*w, gg.Pz())
#    gg.SetVectM(spatial_gg, gg.M())

# jets
#    w = w_ggHH2_hbb_pt.GetBinContent(ibin)
#FIXME    print  chain.gr_hbb_p4_pt, chain.jj_pt
#FIXME    w = w * chain.gr_hbb_p4_pt / chain.jj_pt

#    spatial_jet1 = TVector3(jet1.Px()*w, jet1.Py()*w, jet1.Pz())
#    jet1.SetVectM(spatial_jet1, jet1.M())

#    spatial_jet2 = TVector3(jet2.Px()*w, jet2.Py()*w, jet2.Pz())
#    jet2.SetVectM(spatial_jet2, jet2.M())

#    spatial_jj = TVector3(jj.Px()*w, jj.Py()*w, jj.Pz())
#    jj.SetVectM(spatial_jj, jj.M())
 
# Storing the info for output
    reco_w_pho1_pt[0] = pho1.Pt()
    reco_w_pho1_eta[0] = pho1.Eta()
    reco_w_pho1_phi[0] = pho1.Phi()
    reco_w_pho1_mass[0] = pho1.M()
    reco_w_pho2_pt[0] = pho2.Pt()
    reco_w_pho2_eta[0] = pho2.Eta()
    reco_w_pho2_phi[0] = pho2.Phi()
    reco_w_pho2_mass[0] = pho2.M()

    reco_w_jet1_pt[0] = jet1.Pt()
    reco_w_jet1_eta[0] = jet1.Eta()
    reco_w_jet1_phi[0] = jet1.Phi()
    reco_w_jet1_mass[0] = jet1.M()
    reco_w_jet2_pt[0] = jet2.Pt()
    reco_w_jet2_eta[0] = jet2.Eta()
    reco_w_jet2_phi[0] = jet2.Phi()
    reco_w_jet2_mass[0] = jet2.M()

    reco_w_jj_pt[0] = jj.Pt()
    reco_w_jj_eta[0] = jj.Eta()
    reco_w_jj_phi[0] = jj.Phi()
    reco_w_jj_mass[0] = jj.M()

    reco_w_gg_pt[0] = gg.Pt()
    reco_w_gg_eta[0] = gg.Eta()
    reco_w_gg_phi[0] = gg.Phi()
    reco_w_gg_mass[0] = gg.M()
    
    weight_lambda[0] = w
    costhetastar_CS[0] = costhetastar
#    dEta_gg_jj[0] = d_eta
#    dPhi_gg_jj[0] = d_phi
#    if ievt % 1000 == 0:
#    print hbb_pt[0] - hbb_new_pt[0], hbb_pt[0], hbb_new_pt[0], w, gr_hbb_p4_pt, ibin
    outtree.Fill()
   
# write the tree into the output file and close the file
outfile.Write()
outfile.Close()

