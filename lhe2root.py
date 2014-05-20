#!/usr/bin/env python
# Python modules
import argparse
import numpy
import xml.etree.ElementTree as ET
# Homemade modules
from particles import *
# Argument parser setup
parser = argparse.ArgumentParser()
parser.add_argument("--infile", help="File in LHE format to be read", nargs='?', type=str, default="GluGluToHHTo2B2G_0hhh_M-125_8TeV_madgraph.lhe")
parser.add_argument("--outfile", help="File in ROOT format for output", nargs='?', type=str, default="test.root")
parser.add_argument("--outtree", help="output TTree name", nargs='?', type=str, default="test")
args = parser.parse_args()
# ROOT setup
import ROOT
from ROOT import gROOT
from ROOT import TTree, TFile
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

print root, root.tag, root.attrib, len(root)
for ichild, child in enumerate(root):
    if child.tag != "event":
        continue
    if ichild > 1:
        break
#    print ichild, child.tag, child.attrib, child.text
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
        particlelist.append( particle(iLHEp, *LHEp) )
    eventlist.append(particlelist)

########################################################################################
# WRITE THE ROOT FILE
########################################################################################
outfile = TFile(args.outfile, "recreate")
outtree = TTree(args.outtree, args.outtree)

# Prepare the output variables
gr_g1_p4_pt = numpy.zeros(1, dtype=float)
# Create the branch to store the variable in the tree
outtree.Branch('gr_g1_p4_pt', gr_g1_p4_pt, 'gr_g1_p4_pt/D')
# Loop over events
for event in eventlist:
    for p in event:
        print p.index, p.pid, p.pt
        if p.pid == 22:
            gr_g1_p4_pt[0] = p.pt
    outtree.Fill()

# write the tree into the output file and close the file
outfile.Write()
outfile.Close()

