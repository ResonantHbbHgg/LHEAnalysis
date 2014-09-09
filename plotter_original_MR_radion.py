#!/usr/bin/env python
#// Small dumb code to play with trees
#// O.Bondu, F. Bojarski (May 2014)
# Various python imports
from os import path
from math import log10, pow
# ROOT setup
import ROOT
from ROOT import TFile, TTree, TLine, TChain, TCanvas, TH1D, TLatex, TLegend, TLorentzVector, gStyle
ROOT.gROOT.Reset()
ROOT.gROOT.SetBatch()
ROOT.gROOT.ProcessLine(".x setTDRStyle.C")
ROOT.TGaxis.SetMaxDigits(3);

c1 = TCanvas()
afs_plottree = "/afs/cern.ch/work/o/obondu/public/forRadion/plotTrees/v10/"
#eos_tree = "root://eoscms//eos/cms/store/cmst3/group/hbbhgg/H2GGLOBE/Radion/trees/radion_tree_v09"

intL = 19706.
samples = []
# samples.append([ name, dirpath, subdir, file, tree, weight, color, style, linestyle, label , sigma , N])

samples.append(["MR = 300", "", "../LHEAnalysis", "original_MR_300.root", "test", "1", ROOT.kBlue, 0, 1, "ggHH, MR = 300", 0.0000212, 100000])
samples.append(["MR = 500", "", "../LHEAnalysis", "original_MR_500.root", "test", "1", ROOT.kGreen, 0, 1, "ggHH, MR = 500", 0.0000212, 100000])
samples.append(["MR = 700", "", "../LHEAnalysis", "original_MR_700.root", "test", "1", ROOT.kRed, 0, 1, "ggHH, MR = 700", 0.0000212, 100000])
samples.append(["MR = 1100", "", "../LHEAnalysis", "original_MR_1100.root", "test", "1", ROOT.kYellow, 0, 1, "ggHH, MR = 1100", 0.0000212, 100000])
samples.append(["MR = 2000", "", "../LHEAnalysis", "original_MR_2000.root", "test", "1", ROOT.kCyan, 0, 1, "ggHH, MR = 2000", 0.0000212, 100000])

gStyle.SetPalette(1)

#####plots.append([ name2, variable, cut, norm, binning, title, additional_info, cutline, cutline2 ])
plots = []

hDivided = []
#gamma1 and gamma2
#plots.append(["pho1_pt", "pho1_pt", "", 1., "(100, 0, 500)", "p_{T}^{#gamma1} (GeV)", "", "", ""])
#plots.append(["pho2_pt", "pho2_pt", "", 1., "(100, 0, 500)", "p_{T}^{#gamma2} (GeV)", "", "", ""])
#plots.append(["pho1_eta", "pho1_eta", "", 1., "(50, -4, 4)", "eta^{#gamma1} ", "", "", ""])
#plots.append(["pho2_eta", "pho2_eta", "", 1., "(50, -4, 4)", "eta^{#gamma2} ", "", "", ""])
#plots.append(["pho1_phi", "pho1_phi", "", 1., "(50, -4, 4)", "phi^{#gamma1} (radian)", "", "", ""])
#plots.append(["pho2_phi", "pho2_phi", "", 1., "(50, -4, 4)", "phi^{#gamma2} (radian)", "", "", ""])
#plots.append(["pho1_e", "pho1_e", "", 1., "(100, 0, 500)", "e^{#gamma1} (GeV)", "", "", ""])
#plots.append(["pho2_e", "pho2_e", "", 1., "(100, 0, 500)", "e^{#gamma2} (GeV)", "", "", ""])

#jet1 and jet2
#plots.append(["jet1_pt", "jet1_pt", "", 1., "(100, 0, 500)", "p_{T}^{#jet1} (GeV)", "", "", ""])
#plots.append(["jet2_pt", "jet2_pt", "", 1., "(100, 0, 500)", "p_{T}^{#jet2} (GeV)", "", "", ""])
#plots.append(["jet1_eta", "jet1_eta", "", 1., "(50, -4, 4)", "eta^{#jet1} ", "", "", ""])
#plots.append(["jet2_eta", "jet2_eta", "", 1., "(50, -4, 4)", "eta^{#jet2} ", "", "", ""])
#plots.append(["jet1_phi", "jet1_phi", "", 1., "(50, -4, 4)", "phi^{#jet1} (radian)", "", "", ""])
#plots.append(["jet2_phi", "jet2_phi", "", 1., "(50, -4, 4)", "phi^{#jet2} (radian)", "", "", ""])
#plots.append(["jet1_e", "jet1_e", "", 1., "(100, 0, 500)", "e^{#jet1} (GeV)", "", "", ""])
#plots.append(["jet2_e", "jet2_e", "", 1., "(100, 0, 500)", "e^{#jet2} (GeV)", "", "", ""])

#two gammas
#plots.append(["gg_pt", "gg_pt", "", 1., "(100, 0, 500)", "p_{T}^{#gg} (GeV)", "", "", ""])
#plots.append(["gg_eta", "gg_eta", "", 1., "(50, -6, 6)", "eta^{#gg} ", "", "", ""])
#plots.append(["gg_phi", "gg_phi", "", 1., "(50, -4, 4)", "phi^{#gg} (radian)", "", "", ""])
#plots.append(["gg_e", "gg_e", "", 1., "(100, 0, 1000)", "e^{#gg} (GeV)", "", "", ""])
#plots.append(["gg_mass", "gg_mass", "", 1., "(100, 0, 250)", "mass^{#gg} (GeV)", "", "", ""])

#two jets
#plots.append(["jj_pt", "jj_pt", "", 1., "(100, 0, 500)", "p_{T}^{#jj} (GeV)", "", "", ""])
#plots.append(["jj_eta", "jj_eta", "", 1., "(50, -6, 6)", "eta^{#jj} ", "", "", ""])
#plots.append(["jj_phi", "jj_phi", "", 1., "(50, -4, 4)", "phi^{#jj} (radian)", "", "", ""])
#plots.append(["jj_e", "jj_e", "", 1., "(100, 0, 1000)", "e^{#jj} (GeV)", "", "", ""])
#plots.append(["jj_mass", "jj_mass", "", 1., "(100, 0, 250)", "mass^{#jj} (GeV)", "", "", ""])

#two jets and two gammas
#plots.append(["ggjj_pt", "ggjj_pt", "", 1., "(100, 0, 500)", "p_{T}^{#ggjj} (GeV)", "", "", ""])
#plots.append(["ggjj_eta", "ggjj_eta", "", 1., "(50, -6, 6)", "eta^{#ggjj} ", "", "", ""])
#plots.append(["ggjj_phi", "ggjj_phi", "", 1., "(50, -4, 4)", "phi^{#ggjj} (radian)", "", "", ""])
#plots.append(["ggjj_e", "ggjj_e", "", 1., "(500, 0, 1500)", "e^{#ggjj} (GeV)", "", "", ""])
#plots.append(["ggjj_mass", "ggjj_mass", "", 1., "(500, 0, 1000)", "mass^{#ggjj} (GeV)", "", "", ""])

#gamma1 and gamma2 (lhe2root.py)
#plots.append(["gr_g1_p4_pt_lam_10", "gr_g1_p4_pt", "", 1., "(50, 0, 500)", "p_{T}^{#gamma1} (GeV)", "", "", ""])
#plots.append(["gr_g1_p4_eta_lam_10", "gr_g1_p4_eta", "", 1., "(50, -4, 4)", "eta^{#gamma1} ", "", "", ""])
#plots.append(["gr_g1_p4_phi_lam_10", "gr_g1_p4_phi", "", 1., "(50, -4, 4)", "phi^{#gamma1} (radian)", "", "", ""])
#plots.append(["gr_g1_p4_mass_lam_10", "gr_g1_p4_mass", "", 1., "(100, 0, 200)", "mass^{#gamma1} (GeV", "", "", ""])
#plots.append(["gr_g2_p4_pt_lam_10", "gr_g2_p4_pt", "", 1., "(50, 0, 500)", "p_{T}^{#gamma2} (GeV)", "", "", ""])
#plots.append(["gr_g2_p4_eta_lam_10", "gr_g2_p4_eta", "", 1., "(50, -4, 4)", "eta^{#gamma2} ", "", "", ""])
#plots.append(["gr_g2_p4_phi_lam_10", "gr_g2_p4_phi", "", 1., "(50, -4, 4)", "phi^{#gamma2} (radian)", "", "", ""])
#plots.append(["gr_g2_p4_mass_lam_10", "gr_g2_p4_mass", "", 1., "(100, 0, 200)", "mass^{#gamma2} (GeV)", "", "", ""])

# b quark and anti b quark (lhe2root.py)
#plots.append(["gr_b1_p4_pt_lam_10", "gr_b1_p4_pt", "", 1., "(50, 0, 500)", "p_{T}^{#b1} (GeV)", "", "", ""])
#plots.append(["gr_b1_p4_eta_lam_10", "gr_b1_p4_eta", "", 1., "(50, -4, 4)", "eta^{#b1} ", "", "", ""])
#plots.append(["gr_b1_p4_phi_lam_10", "gr_b1_p4_phi", "", 1., "(50, -4, 4)", "phi^{#b1} (radian)", "", "", ""])
#plots.append(["gr_b1_p4_mass_lam_10", "gr_b1_p4_mass", "", 1., "(100, 0, 200)", "mass^{#b1} (GeV", "", "", ""])
#plots.append(["gr_b2_p4_pt_lam_10", "gr_b2_p4_pt", "", 1., "(50, 0, 500)", "p_{T}^{#b2} (GeV)", "", "", ""])
#plots.append(["gr_b2_p4_eta_lam_10", "gr_b2_p4_eta", "", 1., "(50, -4, 4)", "eta^{#b2} ", "", "", ""])
#plots.append(["gr_b2_p4_phi_lam_10", "gr_b2_p4_phi", "", 1., "(50, -4, 4)", "phi^{#b2} (radian)", "", "", ""])
#plots.append(["gr_b2_p4_mass_lam_10", "gr_b2_p4_mass", "", 1., "(100, 0, 200)", "mass^{#b2} (GeV)", "", "", ""])

# hgg and hbb (lhe2root.py)
#plots.append(["gr_hgg_p4_pt_lam_10", "gr_hgg_p4_pt", "", 1., "(50, 0, 500)", "p_{T}^{#hgg} (GeV)", "", "", ""])
#plots.append(["gr_hgg_p4_eta_lam_10", "gr_hgg_p4_eta", "", 1., "(50, -4, 4)", "eta^{#hgg} ", "", "", ""])
#plots.append(["gr_hgg_p4_phi_lam_10", "gr_hgg_p4_phi", "", 1., "(50, -4, 4)", "phi^{#hgg} (radian)", "", "", ""])
#plots.append(["gr_hgg_p4_mass_lam_10", "gr_hgg_p4_mass", "", 1., "(100, 0, 200)", "mass^{#hgg} (GeV", "", "", ""])
#plots.append(["gr_hbb_p4_pt_lam_10", "gr_hbb_p4_pt", "", 1., "(50, 0, 500)", "p_{T}^{#hbb} (GeV)", "", "", ""])
#plots.append(["gr_hbb_p4_pt_costhetastar_0_0.3", "gr_hbb_p4_pt", "costhetastar_CS < 0.3", 1., "(50, 0, 500)", "p_{T}^{#hbb} (GeV)", "|cos#theta*|^{#hh} < 0.3", "", ""])
##plots.append(["gr_hbb_p4_pt_costhetastar_0.3_0.6", "gr_hbb_p4_pt", "costhetastar_CS > 0.3 && costhetastar_CS < 0.6", 1., "(50, 0, 500)", "p_{T}^{#hbb} (GeV)", "0.3 < |cos#theta*|^{#hh} < 0.6", "", ""])
#plots.append(["gr_hbb_p4_pt_costhetastar_0.6_0.9", "gr_hbb_p4_pt", "costhetastar_CS > 0.6 && costhetastar_CS < 0.9", 1., "(50, 0, 500)", "p_{T}^{#hbb} (GeV)", "0.6 < |cos#theta*|^{#hh} < 0.9", "", ""])
#plots.append(["gr_hbb_p4_pt_costhetastar_0.9_1", "gr_hbb_p4_pt", "costhetastar_CS > 0.9 && costhetastar_CS < 1", 1., "(50, 0, 500)", "p_{T}^{#hbb} (GeV)", "0.9 < |cos#theta*|^{#hh} < 1", "", ""])
#plots.append(["gr_hbb_p4_eta_lam_10", "gr_hbb_p4_eta", "", 1., "(50, -4, 4)", "eta^{#hbb} ", "", "", ""])
#plots.append(["gr_hbb_p4_phi_lam_10", "gr_hbb_p4_phi", "", 1., "(50, -4, 4)", "phi^{#hbb} (radian)", "", "", ""])
#plots.append(["gr_hbb_p4_mass_lam_10", "gr_hbb_p4_mass", "", 1., "(100, 0, 200)", "mass^{#hbb} (GeV)", "", "", ""])
##plots.append(["gr_costhetastar_CS_pt_2D", "gr_hbb_p4_pt:costhetastar_CS_eta_limit", "", 1., "(20, 0, 1, 250, 0, 500)", "|cos#theta*|^{#hh}", "", "", ""])
plots.append(["gr_costhetastar_CS_5", "abs(gr_hbbhgg_costhetastar_CS)", "", 1., "(20, 0, 1)", "|cos#theta*|^{#hh}", "", "", ""])
#plots.append([gr_costhetastar_pt_0_50", "costhetastar_CS", "gr_hbb_p4_pt < 50", 1., "(50, 0, 1)", "|cos#theta*|^{#hh}", "p_{T}^{h} < 50 GeV", "", ""])
#plots.append(["gr_costhetastar_pt_50_100", "costhetastar_CS", "gr_hbb_p4_pt > 50 && gr_hbb_p4_pt < 100", 1., "(50, 0, 1)", "|cos#theta*|^{#hh}", "50 < p_{T}^{h} < 100 GeV", "", ""])
#plots.append(["gr_costhetastar_pt_100_150", "costhetastar_CS", " gr_hbb_p4_pt > 100 && gr_hbb_p4_pt < 150", 1., "(50, 0, 1)", "|cos#theta*|^{#hh}", "100 < p_{T}^{h} < 150 GeV", "", ""])
#plots.append(["gr_costhetastar_pt_150_200", "costhetastar_CS", "gr_hbb_p4_pt > 150 && gr_hbb_p4_pt < 200", 1., "(50, 0, 1)", "|cos#theta*|^{#hh}", "150 < p_{T}^{h} < 200 GeV", "", ""])
#plots.append(["gr_costhetastar_pt_200_500", "costhetastar_CS", "gr_hbb_p4_pt > 200 && gr_hbb_p4_pt < 1000", 1., "(50, 0, 1)", "|cos#theta*|^{#hh}", "200 < p_{T}^{h} < 500 GeV", "", ""])
#plots.append(["2D_gr_delta_eta", "dEta_gg_jj", "", 1., "(50, -6, 6)", "#Delta#eta^{#hh}", "", "", ""])
#plots.append(["2D_gr_delta_phi", "dPhi_gg_jj", "", 1., "(50, -6, 6)", "#Delta#phi^{#hh}", "", "", ""])

#plots.append(["gr_radion_p4_pt_lam_10", "gr_radion_p4_pt", "", 1., "(50, 0, 500)", "p_{T}^{#radion} (GeV)", "", "", ""])
#plots.append(["gr_radion_p4_eta_lam_10", "gr_radion_p4_eta", "", 1., "(50, -4, 4)", "eta^{#radion} ", "", "", ""])
#plots.append(["gr_radion_p4_phi_lam_10", "gr_radion_p4_phi", "", 1., "(50, -4, 4)", "phi^{#radion} (radian)", "", "", ""])
plots.append(["gr_radion_p4_mass_5", "gr_radion_p4_mass", "", 1., "(40, 200, 1000)", "mass^{#radion} (TeV)", "", "", ""])

# dividing the plots

# increased range of pt
#plots.append(["gr_g1_p4_pt", "gr_g1_p4_pt", "", 1., "(500, 0, 1000)", "p_{T}^{#gamma1} (GeV)", "", "", ""])
#plots.append(["gr_g2_p4_pt", "gr_g2_p4_pt", "", 1., "(500, 0, 1000)", "p_{T}^{#gamma2} (GeV)", "", "", ""])
#plots.append(["gr_b1_p4_pt", "gr_b1_p4_pt", "", 1., "(500, 0, 1000)", "p_{T}^{#b1} (GeV)", "", "", ""])
#plots.append(["gr_b2_p4_pt", "gr_b2_p4_pt", "", 1., "(500, 0, 1000)", "p_{T}^{#b2} (GeV)", "", "", ""])
#plots.append(["gr_hgg_p4_pt", "gr_hgg_p4_pt", "", 1., "(500, 0, 1000)", "p_{T}^{#hgg} (GeV)", "", "", ""])
#plots.append(["gr_hbb_p4_pt", "gr_hbb_p4_pt", "", 1., "(500, 0, 1000)", "p_{T}^{#hbb} (GeV)", "", "", ""])

for iplot, [name2, variable, cut, norm, binning, title, additional_info, cutline, cutline2] in enumerate(plots):
    c1 = TCanvas()
    legend = TLegend(0.45, 0.82, 0.90, 0.93, "")
    legend.SetTextSize(0.025)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetLineColor(ROOT.kWhite)
    legend.SetShadowColor(ROOT.kWhite)
    if len(samples) > 1:
        legend.SetNColumns(2)
    xnbin, xlow, xhigh = map(float, binning.strip().strip("()").split(","))
    ymax = -1
    ymin = 10000000
    firsthistname = ""
    if cut == "": cut = "1"

    hDivided = []
    

    for ifile, [ name, dirpath, subdir, file, tree, sample_weight, color, style, linestyle, label , sigma , N] in enumerate(samples):
#        print ifile, file, color, style, label
        chain = TChain(tree)
        chain.Add( path.join(dirpath, subdir, file) )
        sample_cut = cut
        if norm == 1.:
            sample_cut = "((" + sample_cut + ") * (" + sample_weight + "))/" + str( chain.GetEntries() )
        else:
            sample_cut = "((" + sample_cut + ") * (" + sample_weight + ")) * (" + str(sigma) + " * " + str(intL) + ")/" + str(N)
        option = ""
        print sample_cut
        if ifile != 0:
            option = "same"
        chain.Draw(variable + ">>h_tmp" + binning, sample_cut, option)
  
        # Clsosmetics
        h = ROOT.gDirectory.Get("h_tmp")
        h.SetName(name + "_" + name2 + "_" + str(ifile))
        if ifile == 0:
            firsthistname = name + "_" + name2 + "_" + str(ifile)
        h.SetLineWidth(3)
        h.SetLineColor(color)
        h.SetFillColor(color)
        h.SetFillStyle(style)
        h.SetLineStyle(linestyle)
        h.GetXaxis().SetTitle( title )       

        unit = ""
        if title.find("(") != -1:
            unit = title[title.find("(")+1:title.find(")")]
        if norm == 1.:
            h.GetYaxis().SetTitle( "Norm. to unity / ( " + str(((xhigh - xlow) / xnbin)) + " " + unit + " )")
        else:
            h.GetYaxis().SetTitle( "# events / ( " + str(((xhigh - xlow) / xnbin)) + " " + unit + " )")
        legend.AddEntry(h.GetName(), label)
        ymax = max(ymax, h.GetMaximum())
        ymin = min(ymin, h.GetMinimum(0.0))
        del chain, h



    ymin_lin = ymin / 10.
    yrange_lin = ymax - ymin_lin
    ymax_lin = .25 * yrange_lin + ymax
    yrange_log = (log10(ymax) - log10(ymin)) / .77
    ymax_log = pow(10., .25*yrange_log + log10(ymax))
    ymin_log = pow(10., log10(ymin) - .03*yrange_log)

    latexLabel = TLatex()
    latexLabel.SetTextSize(.03)
    latexLabel.SetNDC()
    latexLabel.DrawLatex(.25, .96, "CMS Internal     L = 19.7 fb^{-1}     #sqrt{s} = 8 TeV")
    latexLabel.DrawLatex(.28, .79, "" )
    ROOT.gPad.RedrawAxis()
    legend.Draw()
    c1.Update()

    line = TLine()
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line2 = TLine()
    line2.SetLineStyle(2)
    line2.SetLineWidth(2)

    h = ROOT.gDirectory.Get(firsthistname)
    h.SetMaximum(ymax_lin)
    h.SetMinimum(ymin_lin)
    if cutline != "":
        line.SetX1(cutline); line.SetY1(ymin_lin); line.SetX2(cutline); line.SetY2(ymax)
        line.Draw("same")
    if cutline2 != "":
        line2.SetX1(cutline2); line2.SetY1(ymin_lin); line2.SetX2(cutline2); line2.SetY2(ymax)
        line2.Draw("same")
    c1.Update()
    c1.Print("plots_MR_radion_pdf/" + name2 + ".pdf")
    #c1.Print("plots_Yt_gif/" + name2 + ".gif")
    #c1.Print("plots_Yt_root/" + name2 + ".root")


    c1.SetLogy(1)
    h.SetMaximum(ymax_log)
    h.SetMinimum(ymin_log)
    h.GetYaxis().SetRangeUser(ymin_log, ymax_log)
    if cutline != "":
        line.SetX1(cutline); line.SetY1(ymin_log); line.SetX2(cutline); line.SetY2(ymax)
        line.Draw("same")
    if cutline2 != "":
        line2.SetX1(cutline2); line2.SetY1(ymin_log); line2.SetX2(cutline2); line2.SetY2(ymax)
        line2.Draw("same")
    c1.Update()
    c1.Print("plots_MR_radion_pdf/" + name2 + "_log.pdf")
  #  c1.Print("plots_Yt_gif/" + name2 + "_log.gif")
  #  c1.Print("plots_Yt_root/" + name2 + "_log.root")
    c1.SetLogy(0)

    c1.Clear()

    del c1
    
