#### user inputs ####

filename = 'run_2025_02_26'

#directories in root file
vars = [
        # "JetPt",
        # "JetEta",
        # "JetPhi",
        # "BJetPt",
        # "BJetEta",
        # "BJetPhi",
        # "LeadElePt",
        # "LeadEleEta",
        # "NJet",
        # "NLooseBTag",
        # "NMediumBTag",
        "MET",
        # "NTau",
        "MTSelect",
        # "LeadMuPt",
        # "LeadMuEta"
]

#files within each root directory
channels = [
            "Pre_Lep",
            "Pre_e",
            "Pre_u",
            "Pre_MB_Lep",
            "Pre_MB_e",
            "Pre_MB_u",
            "Lep",
            "Lep_e",
            "Lep_u",
            "Lep_MT",
            "Lep_MT_e",
            "Lep_MT_u",
            "Lep_MB",
            "Lep_MB_e",
            "Lep_MB_u"            
]

year = "2018"

signals = [
           "TT_powheg_pythia8"
]

backgrounds = [
               "Multijet",
               "DYToLL",
               "WToLNu"
]



#### main programs ####

import os
import ROOT
from array import array

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()

def DrawROC(sig, bkg, pname='sync.pdf', isLegend=False):
    # Remove specific strings from the canvas name
    strings_to_remove = ["/", ".png"]
    cname = pname
    for string in strings_to_remove:
        cname = cname.replace(string, "_")
    
    # Create a canvas
    c1 = ROOT.TCanvas(cname, "", 700, 700)
    
    # Initialize lists for ROC curve points
    x1 = []
    y1 = []
    
    # Calculate the ROC curve points
    for i in range(1, sig.GetNbinsX() + 1):
        x1.append(sig.Integral(i, sig.GetNbinsX()) / sig.Integral(1, sig.GetNbinsX()))
        y1.append(1-bkg.Integral(i, bkg.GetNbinsX()) / bkg.Integral(1, bkg.GetNbinsX()))
    
    # Create TGraph objects for the ROC curves
    g1 = ROOT.TGraph(len(x1), array('d', x1), array('d', y1))
    
    # Set graph styles
    g1.SetLineColor(2)
    g1.SetMarkerColor(2)
    g1.SetMarkerStyle(22)
    g1.SetMarkerSize(1.5)
    g1.GetXaxis().SetTitle("sig efficiency")
    g1.GetYaxis().SetTitle("1-bkg efficiency")
    g1.GetYaxis().SetNdivisions(505)
    g1.GetXaxis().SetLimits(0, 1)
    g1.GetYaxis().SetRangeUser(0.00, 1)
    
    # Draw the ROC curves
    g1.Draw('ALP')
            
    cms_lat=ROOT.TLatex()
    cms_lat.SetTextSize(0.05)
    cms_lat.DrawLatex(0,1.01,"#bf{CMS} #scale[0.8]{#font[52]{Run2 Preliminary}}")
    cms_lat.DrawLatex(0.8,1.01,"#scale[0.9]{(13 TeV)}")
    
    for i in range(6):
        cms_lat.SetTextColor(1)
        cms_lat.SetTextSize(0.03)
#        cms_lat.DrawLatex(x1[i+1], y1[i+1]+0.02, f'{10*(i+1):.0f}') python3
#        cms_lat.DrawLatex(x1[i+1], y1[i+1]+0.02, '{}'.format(10*(i+1))) #python2

    # Add legend if specified
    if isLegend:
        leg = ROOT.TLegend(0.15, 0.50, 0.15, 0.30)
        leg.SetBorderSize(0)
        leg.SetFillColor(10)
        leg.SetLineColor(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.04)
        leg.SetHeader('Semi leptonic final states')
        leg.AddEntry(None, channel, "")
        leg.AddEntry(None, "sig: " + sigNames2, "")
        leg.AddEntry(None, "bkg: " + bkgNames2, "")
        leg.AddEntry(g1, var, 'lp')
        leg.Draw()
    
    # Save the canvas
    c1.Print(pname)

def get_histogram(file, hist_name):
    hist = file.Get(hist_name)
    if type(hist) == ROOT.TObject:
        print("no file for " + hist_name)
        return
    else:
        return hist


#append signal names together
#for file names
sigNames = signals[0]
for signal in signals[1:]:
    sigNames = sigNames + '_' + signal
#for plot labels
sigNames2 = signals[0]
for signal in signals[1:]:
    sigNames2 = sigNames2 + ', ' + signal

#append background names together
#for file names
bkgNames = backgrounds[0]
for bkg in backgrounds[1:]:
    bkgNames = bkgNames + '_' + bkg
#for plot labels
bkgNames2 = backgrounds[0]
for bkg in backgrounds[1:]:
    bkgNames2 = bkgNames2 + ', ' + bkg

file = ROOT.TFile.Open(filename + '.root')

# make each ROC curve file
for var in vars:
    for channel in channels:
        
        # get first signal histogram
        hist_name = var + "/" + signals[0] + "_" + year + "_" + channel
        sig = get_histogram(file, hist_name)

        # add rest of the signal histograms to the first
        for signal in signals[1:]:
            hist_name = var + "/" + signal + "_" + year + "_" + channel
            sig.Add(get_histogram(file, hist_name))

        # get first background histogram
        hist_name = var + "/" + backgrounds[0] + "_" + year + "_" + channel
        bkg = get_histogram(file, hist_name)
        
        # add rest of background histograms to the first
        for background in backgrounds[1:]:
            hist_name = var + "/" + background + "_" + year + "_" + channel
            bkg.Add(get_histogram(file, hist_name))

    #make the ROC curve file
        if sig and bkg:
            dir_name = 'ROC_Plots/' + filename + '/' + var + '/'
            try:
                os.makedirs(dir_name)
            except:
                None
            outfile = filename + '_' + channel + '_' + var + '_' + year + '_' + sigNames + '_VS_' + bkgNames + '.png'
            DrawROC(sig,bkg,dir_name + outfile,True)
