#!/usr/bin/env python
# Python modules
import argparse
import numpy
import xml.etree.ElementTree as ET
# Homemade modules
from particles import *
# Argument parser setup
parser = argparse.ArgumentParser()
parser.add_argument("--infile", help="File in LHE format to be read", nargs='?', type=str, default="GluGluToHHTo2B2G_2hhh_M-125_8TeV_madgraph.lhe")
parser.add_argument("--outfile", help="File in ROOT format for output", nargs='?', type=str, default="ggHH2.root")
parser.add_argument("--outtree", help="output TTree name", nargs='?', type=str, default="test")

args = parser.parse_args()
# ROOT setup
import ROOT
from ROOT import gROOT
from ROOT import TTree, TFile, TLorentzVector
# ROOT initialization
gROOT.Reset()
gROOT.SetBatch()

########################################################################################
######
# Let's go inefficient but human-readable :
#   - first get the relevant data out of the file (ie particle gen info)
#   - then loop over the data to store stuff into ROOT format
#####
########################################################################################

########################################################################################
##### READ THE LHE FILE
########################################################################################
LHEeventlist = ET.parse(args.infile)
root = LHEeventlist.getroot()
eventlist = []

#print root, root.tag, root.attrib, len(root)
for ichild, child in enumerate(root):
    if child.tag != "event":
        continue
#    if ichild > 10:
#        break
    if ichild % 10000 == 0:
        print ichild
    #print ichild, child.tag, child.attrib, child.text
#    print type(child.text)
    # Get the different lines
    LHEparticlelist = child.text.split("\n")
    # List cleanup
    LHEparticlelist = filter(None, LHEparticlelist)
    LHEparticlelist.pop(0)
    # Get the different particles
    particlelist = []
    for iLHEp, LHEp in enumerate(LHEparticlelist):
        # cleanup the particles and prepare the TLorentzVector construction
        LHEp = LHEp.strip()
        LHEp = LHEp.split()
        LHEp = map(int, LHEp[:6]) + map(float, LHEp[6:])
        #print LHEp
        particlelist.append( particle(iLHEp, *LHEp) )
    eventlist.append(particlelist)

########################################################################################
# WRITE THE ROOT FILE
########################################################################################
outfile = TFile(args.outfile, "recreate")
outtree = TTree(args.outtree, args.outtree)

# Prepare the output variables
gr_g1_p4_pt = numpy.zeros(1, dtype=float)
gr_g1_p4_eta = numpy.zeros(1, dtype=float)
gr_g1_p4_phi = numpy.zeros(1, dtype=float)
gr_g1_p4_mass = numpy.zeros(1, dtype=float)
gr_g2_p4_pt = numpy.zeros(1, dtype=float)
gr_g2_p4_eta = numpy.zeros(1, dtype=float)
gr_g2_p4_phi = numpy.zeros(1, dtype=float)
gr_g2_p4_mass = numpy.zeros(1, dtype=float)

gr_b1_p4_pt = numpy.zeros(1, dtype=float)
gr_b1_p4_eta = numpy.zeros(1, dtype=float)
gr_b1_p4_phi = numpy.zeros(1, dtype=float)
gr_b1_p4_mass = numpy.zeros(1, dtype=float)
gr_b2_p4_pt = numpy.zeros(1, dtype=float)
gr_b2_p4_eta = numpy.zeros(1, dtype=float)
gr_b2_p4_phi = numpy.zeros(1, dtype=float)
gr_b2_p4_mass = numpy.zeros(1, dtype=float)

gr_hgg_p4_pt = numpy.zeros(1, dtype=float)
gr_hgg_p4_eta = numpy.zeros(1, dtype=float)
gr_hgg_p4_phi = numpy.zeros(1, dtype=float)
gr_hgg_p4_mass = numpy.zeros(1, dtype=float)

gr_hbb_p4_pt = numpy.zeros(1, dtype=float)
gr_hbb_p4_eta = numpy.zeros(1, dtype=float)
gr_hbb_p4_phi = numpy.zeros(1, dtype=float)
gr_hbb_p4_mass = numpy.zeros(1, dtype=float)

#gr_radion_p4_pt = numpy.zeros(1, dtype=float)
#gr_radion_p4_eta = numpy.zeros(1, dtype=float)
#gr_radion_p4_phi = numpy.zeros(1, dtype=float)
#gr_radion_p4_mass = numpy.zeros(1, dtype=float)


# Create the branch to store the variable in the tree
outtree.Branch('gr_g1_p4_pt', gr_g1_p4_pt, 'gr_g1_p4_pt/D')
outtree.Branch('gr_g1_p4_eta', gr_g1_p4_eta, 'gr_g1_p4_eta/D')
outtree.Branch('gr_g1_p4_phi', gr_g1_p4_phi, 'gr_g1_p4_phi/D')
outtree.Branch('gr_g1_p4_mass', gr_g1_p4_mass, 'gr_g1_p4_mass/D')
outtree.Branch('gr_g2_p4_pt', gr_g2_p4_pt, 'gr_g2_p4_pt/D')
outtree.Branch('gr_g2_p4_eta', gr_g2_p4_eta, 'gr_g2_p4_eta/D')
outtree.Branch('gr_g2_p4_phi', gr_g2_p4_phi, 'gr_g2_p4_phi/D')
outtree.Branch('gr_g2_p4_mass', gr_g2_p4_mass, 'gr_g2_p4_mass/D')

outtree.Branch('gr_b1_p4_pt', gr_b1_p4_pt, 'gr_b1_p4_pt/D')
outtree.Branch('gr_b1_p4_eta', gr_b1_p4_eta, 'gr_b1_p4_eta/D')
outtree.Branch('gr_b1_p4_phi', gr_b1_p4_phi, 'gr_b1_p4_phi/D')
outtree.Branch('gr_b1_p4_mass', gr_b1_p4_mass, 'gr_b1_p4_mass/D')
outtree.Branch('gr_b2_p4_pt', gr_b2_p4_pt, 'gr_b2_p4_pt/D')
outtree.Branch('gr_b2_p4_eta', gr_b2_p4_eta, 'gr_b2_p4_eta/D')
outtree.Branch('gr_b2_p4_phi', gr_b2_p4_phi, 'gr_b2_p4_phi/D')
outtree.Branch('gr_b2_p4_mass', gr_b2_p4_mass, 'gr_b2_p4_mass/D')

outtree.Branch('gr_hgg_p4_pt', gr_hgg_p4_pt, 'gr_hgg_p4_pt/D')
outtree.Branch('gr_hgg_p4_eta', gr_hgg_p4_eta, 'gr_hgg_p4_eta/D')
outtree.Branch('gr_hgg_p4_phi', gr_hgg_p4_phi, 'gr_hgg_p4_phi/D')
outtree.Branch('gr_hgg_p4_mass', gr_hgg_p4_mass, 'gr_hgg_p4_mass/D')

outtree.Branch('gr_hbb_p4_pt', gr_hbb_p4_pt, 'gr_hbb_p4_pt/D')
outtree.Branch('gr_hbb_p4_eta', gr_hbb_p4_eta, 'gr_hbb_p4_eta/D')
outtree.Branch('gr_hbb_p4_phi', gr_hbb_p4_phi, 'gr_hbb_p4_phi/D')
outtree.Branch('gr_hbb_p4_mass', gr_hbb_p4_mass, 'gr_hbb_p4_mass/D')

#outtree.Branch('gr_radion_p4_pt', gr_radion_p4_pt, 'gr_radion_p4_pt/D')
#outtree.Branch('gr_radion_p4_eta', gr_radion_p4_eta, 'gr_radion_p4_eta/D')
#outtree.Branch('gr_radion_p4_phi', gr_radion_p4_phi, 'gr_radion_p4_phi/D')
#outtree.Branch('gr_radion_p4_mass', gr_radion_p4_mass, 'gr_radion_p4_mass/D')


# Loop over events
for event in eventlist:
    pt_tempb1 = 0
    eta_tempb1 = 0
    phi_tempb1 = 0
    mass_tempb1 = 0
    pt_tempb2 = 0
    eta_tempb2 = 0
    phi_tempb2 = 0
    mass_tempb2 = 0
    pt_tempg1 = 0
    eta_tempg1 = 0
    phi_tempg1 = 0
    mass_tempg1 = 0
    pt_tempg2 = 0
    eta_tempg2 = 0
    phi_tempg2 = 0
    mass_tempg2 = 0
    isFirst = True
    for ip, p in enumerate(event):
        #print p.index, p.pid, p.pt
        #gr_radion_p4_pt[0] = p.pt
        #gr_radion_p4_eta[0] = p.eta
        #gr_radion_p4_phi[0] = p.phi
        #gr_radion_p4_mass[0] = p.mass
        if ip == 0 or ip == 1:
            # This is a gluon, we don't care
            continue
        if ip == 2:
            # This is hbb
            gr_hbb_p4_pt[0] = p.pt
            gr_hbb_p4_eta[0] = p.eta
            gr_hbb_p4_phi[0] = p.phi
            gr_hbb_p4_mass[0] = p.mass
            
        
        if p.pid == 5 or p.pid == -5:          
            if isFirst == True :
                pt_tempb1 = p.pt
                eta_tempb1 = p.eta
                phi_tempb1 = p.phi
                mass_tempb1 = p.mass
                isFirst = False
            else:
                pt_tempb2 = p.pt
                eta_tempb2 = p.eta
                phi_tempb2 = p.phi
                mass_tempb2 = p.mass
                isFirst = True   
        
        if ip == 3:
            # This is hgg
            gr_hgg_p4_pt[0] = p.pt
            gr_hgg_p4_eta[0] = p.eta
            gr_hgg_p4_phi[0] = p.phi
            gr_hgg_p4_mass[0] = p.mass
        if p.pid == 22 :
            if isFirst == True :
                pt_tempg1 = p.pt
                eta_tempg1 = p.eta
                phi_tempg1 = p.phi
                mass_tempg1 = p.mass
                isFirst = False
            else:
                pt_tempg2 = p.pt
                eta_tempg2 = p.eta
                phi_tempg2 = p.phi
                mass_tempg2 = p.mass    

        
    if pt_tempb1 >= pt_tempb2 :
        gr_b1_p4_pt[0] = pt_tempb1
        gr_b1_p4_eta[0] = eta_tempb1
        gr_b1_p4_phi[0] = phi_tempb1
        gr_b1_p4_mass[0] = mass_tempb1
        gr_b2_p4_pt[0] = pt_tempb2
        gr_b2_p4_eta[0] = eta_tempb2
        gr_b2_p4_phi[0] = phi_tempb2
        gr_b2_p4_mass[0] = mass_tempb2
    else :
        gr_b2_p4_pt[0] = pt_tempb1
        gr_b2_p4_eta[0] = eta_tempb1
        gr_b2_p4_phi[0] = phi_tempb1
        gr_b2_p4_mass[0] = mass_tempb1
        gr_b1_p4_pt[0] = pt_tempb2
        gr_b1_p4_eta[0] = eta_tempb2
        gr_b1_p4_phi[0] = phi_tempb2
        gr_b1_p4_mass[0] = mass_tempb2

    if pt_tempg1 >= pt_tempg2 :
        gr_g1_p4_pt[0] = pt_tempg1
        gr_g1_p4_eta[0] = eta_tempg1
        gr_g1_p4_phi[0] = phi_tempg1
        gr_g1_p4_mass[0] = mass_tempg1
        gr_g2_p4_pt[0] = pt_tempg2
        gr_g2_p4_eta[0] = eta_tempg2
        gr_g2_p4_phi[0] = phi_tempg2
        gr_g2_p4_mass[0] = mass_tempg2
    else :
        gr_g2_p4_pt[0] = pt_tempg1
        gr_g2_p4_eta[0] = eta_tempg1
        gr_g2_p4_phi[0] = phi_tempg1
        gr_g2_p4_mass[0] = mass_tempg1
        gr_g1_p4_pt[0] = pt_tempg2
        gr_g1_p4_eta[0] = eta_tempg2
        gr_g1_p4_phi[0] = phi_tempg2
        gr_g1_p4_mass[0] = mass_tempg2

   

# Since we are dealing with LO, hbb + hgg is exactly 0
#    hgg = TLorentzVector()
#    hgg.SetPtEtaPhiM(gr_hgg_p4_pt[0], gr_hgg_p4_eta[0], gr_hgg_p4_phi[0], gr_hgg_p4_mass[0])
#    hbb = TLorentzVector()
#    hbb.SetPtEtaPhiM(gr_hbb_p4_pt[0], gr_hbb_p4_eta[0], gr_hbb_p4_phi[0], gr_hbb_p4_mass[0])
#    radion = hbb + hgg
    #print type(radion), radion.Pt(), radion.Pz(), radion.M()
#    gr_radion_p4_pt[0] = radion.Pt()
#    gr_radion_p4_eta[0] = radion.Eta()
#    gr_radion_p4_phi[0] = radion.Phi()
#    gr_radion_p4_mass[0] = radion.M()

    outtree.Fill()
   
# write the tree into the output file and close the file
outfile.Write()
outfile.Close()

