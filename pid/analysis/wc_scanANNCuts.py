from __future__ import division
import itertools
import ROOT
import math
from cmath import sqrt

def main():
    ROOT.gStyle.SetOptStat(0)
    ScanCuts()

# This method is designed to work out the selection efficiency and purity for selecting
# nue CC events from the general NuMI beam spectrum.  It does the event weights itself, so
# it should work even if you provide it with events that aren't distributed according to how
# they are in the NuMI beam (i.e. you can give it separate CCQE, NC samples etc.)
# It works by dividing the sample into four non-intersecting sets which cover all possible events
# (not considering taus): nue CCQE events, numu CCQE events, nue CC-non-QE events, numu CC-non-QE
# events, nue NC events, and numu NC events
# Then it works out how many of each sample are selected, renormalises to have 95%:5% numu:nue, 
# 70%:30% CC:NC, and within the CC sample, 20%:80% QE:non-QE.  Finally it uses these numbers
# to plot the efficiency, purity, and efficiency * purity figure of merit as a function of the 
# ANN variable cuts
def ScanCuts():

    relativeElectronNorm = 0.05
    pidChain = ROOT.TChain("PIDTree_ann")
    pidChain.Add("/unix/chips/jtingey/CHIPS/data/geometries/CHIPS_10kton_3inch_HQE_6perCent/PID/output/CHIPS_10kton_3inch_HQE_6perCent_*_readOutput.root")

    # Number of pure CC and NC events
    numNueCCEvents = pidChain.GetEntries("trueBeamPDG==12 && trueCCEvent");
    numNueNCEvents = pidChain.GetEntries("trueBeamPDG==12 && trueNCEvent");
    numNumuCCEvents = pidChain.GetEntries("trueBeamPDG==14 && trueCCEvent");
    numNumuNCEvents = pidChain.GetEntries("trueBeamPDG==14 && trueNCEvent");

    # Number of each type of CC event 
    numNueCCQEEvents = pidChain.GetEntries("trueBeamPDG==12 && trueCCEvent && trueQEEvent");
    numNumuCCQEEvents = pidChain.GetEntries("trueBeamPDG==14 && trueCCEvent && trueQEEvent");
    numNueCCnonQEEvents = pidChain.GetEntries("trueBeamPDG==12 && trueCCEvent && !trueQEEvent");
    numNumuCCnonQEEvents = pidChain.GetEntries("trueBeamPDG==14 && trueCCEvent && !trueQEEvent");

    # Make the plot as if there were a thousand total events 
    normNumEvents = 1000

    # Work out weights - what fraction of the thousand events are each type?
    # This will break down the 1000 events into the ratios given by the numbers below,
    # accounting for the fact that we might have different-sized samples (e.g. a pure CCQE one
    # in addition to all nue events one)

    weightNueNCEvents  = normNumEvents * (0.05 * 0.3)/numNueNCEvents;
    weightNumuNCEvents = normNumEvents * (0.95 * 0.3)/numNumuNCEvents;
    weightNueCCQEEvents  = normNumEvents * (0.05 * 0.7 * 0.2)/numNueCCQEEvents;
    weightNumuCCQEEvents = normNumEvents * (0.95 * 0.7 * 0.2)/numNumuCCQEEvents;
    weightNueCCnonQEEvents  = normNumEvents * (0.05 * 0.7 * 0.8)/numNueCCnonQEEvents;
    weightNumuCCnonQEEvents = normNumEvents * (0.95 * 0.7 * 0.8)/numNumuCCnonQEEvents;
    print("sum of weights = %f" % sum([weightNueNCEvents*numNueNCEvents, weightNumuNCEvents*numNumuNCEvents, weightNueCCQEEvents*numNueCCQEEvents, weightNumuCCQEEvents*numNumuCCQEEvents,weightNueCCnonQEEvents*numNueCCnonQEEvents, weightNumuCCnonQEEvents*numNumuCCnonQEEvents]))

    print("Total numu NC events = %f and nue NC events = %f" % (weightNumuNCEvents*numNumuNCEvents, weightNueNCEvents*numNueNCEvents))

    # The range and step size to sweep out the ANN variable cuts
    min = 0.7
    max = 1.0
    nBins = 30
    width = (max - min)/nBins

    # Define our histograms: the selection efficiency for each sample of signal or background
    histNumuCCAll = ROOT.TH2D("histNumuCCAll","#nu_{#mu} CC rejection efficiency;CCQE e vs. #mu cut;CCQE vs. NC cut", nBins, min, max, nBins, min, max)
    histNueCCQE = ROOT.TH2D("histNueCCQE","#nu_{e} CCQE selection efficiency;CCQE e vs. #mu cut;CCQE vs. NC cut", nBins, min, max, nBins, min, max)
    histNueCCAll = ROOT.TH2D("histNueCCAll","#nu_{e} CC selection efficiency;CCQE e vs. #mu cut;CCQE vs. NC cut", nBins, min, max, nBins, min, max)
    histNuNCAll = ROOT.TH2D("histNuNCAll","NC rejection efficiency;CCQE e vs. #mu cut;CCQE vs. NC cut", nBins, min, max, nBins, min, max)
    
    # Also the purity and eff * purity
    histNueCCAllPurity = ROOT.TH2D("histNueCCAllPurity","#nu_{e} CC purity;CCQE e vs. #mu cut;CCQE vs. NC cut", nBins, min, max, nBins, min, max)
    
    histNueCCAllEffPurity_CC = ROOT.TH2D("histNueCCAllEffPurity_CC","#nu_{e} CCQE efficiency #times purity;CCQE e vs. #mu cut;CCQE vs. NC cut", nBins, min, max, nBins, min, max)
    histFOM_CC = ROOT.TH2D("histFOM_CC","FOM;CCQE e vs. #mu cut;CCQE vs. NC cut", nBins, min, max, nBins, min, max)
    
    histNueCCAllEffPurity_CCQE = ROOT.TH2D("histNueCCAllEffPurity_CCQE","#nu_{e} CCQE efficiency #times purity;CCQE e vs. #mu cut;CCQE vs. NC cut", nBins, min, max, nBins, min, max)
    histFOM_CCQE = ROOT.TH2D("histFOM_CCQE","FOM;CCQE e vs. #mu cut;CCQE vs. NC cut", nBins, min, max, nBins, min, max)


    # Define the cuts used with TTree::GetEntries to select each type of event
    preselectionCut = "preselected"
    isNumuCCQE      = "trueBeamPDG==14 && trueCCEvent && trueQEEvent"
    isNueCCQE       = "trueBeamPDG==12 && trueCCEvent && trueQEEvent"
    isNumuNC        = "trueBeamPDG==14 && trueNCEvent"
    isNueNC         = "trueBeamPDG==12 && trueNCEvent"
    isNumuCCnonQE   = "trueBeamPDG==14 && trueCCEvent && !trueQEEvent"
    isNueCCnonQE    = "trueBeamPDG==12 && trueCCEvent && !trueQEEvent"
    
    for xBin in xrange(1, nBins+1):
        # Define the cut on the numu vs nue ANN for this step
        print "Cuts applied at the upper edge of each bin"
        xCut = min + xBin * width
        muVsElANNCut = "annNueCCQEvsNumuCCQE >= %f" % xCut

        for yBin in xrange(1, nBins+1):
            # Defince the cut on the nue CCQE vs NC ANN for this step
            yCut = min + yBin * width
            CCVsNCANNCut = "annNueCCQEvsNC >= %f" % yCut

            # Our cut to select each sample is comes from joining together the preselection and ANN 
            # cuts, and then adding the cut for this specific sample

            # First the nue CCQEs
            nueCCQECutString = '&&'.join([isNueCCQE, preselectionCut, muVsElANNCut, CCVsNCANNCut])
            numNueCCQECut    = pidChain.GetEntries(nueCCQECutString)


            # Then the numu CCQEs
            numuCCQECutString = '&&'.join([isNumuCCQE, preselectionCut, muVsElANNCut, CCVsNCANNCut])
            numNumuCCQECut    = pidChain.GetEntries(numuCCQECutString)


            # Then the nue NCs
            nueNCCutString = '&&'.join([isNueNC, preselectionCut, muVsElANNCut, CCVsNCANNCut])
            numNueNCCut    = pidChain.GetEntries(nueNCCutString)
            #print("Get %f nue NC events with cut %f, %f" % (numNueNCCut, xCut, yCut))


            # Then the numu NCs
            numuNCCutString = '&&'.join([isNumuNC, preselectionCut, muVsElANNCut, CCVsNCANNCut])
            numNumuNCCut    = pidChain.GetEntries(numuNCCutString)
            #print("Get %f numu NC events with cut %f, %f" % (numNumuNCCut, xCut, yCut))
            #print("So NC eff = %f / %f = %f\n\n" % (numNumuNCCut*weightNumuNCEvents + numNueNCCut*weightNueNCEvents, (numNumuNCEvents*weightNumuNCEvents + numNueNCEvents*weightNueNCEvents), (numNumuNCCut*weightNumuNCEvents + numNueNCCut*weightNueNCEvents)/(numNumuNCEvents*weightNumuNCEvents + numNueNCEvents*weightNueNCEvents)))

            # The the nue CC but non-CCQEs
            nueCCnonQECutString = '&&'.join([isNueCCnonQE, preselectionCut, muVsElANNCut, CCVsNCANNCut])
            numNueCCnonQECut    = pidChain.GetEntries(nueCCnonQECutString)


            # And finally the numu CC but non-CCQEs
            numuCCnonQECutString = '&&'.join([isNumuCCnonQE, preselectionCut, muVsElANNCut, CCVsNCANNCut])
            numNumuCCnonQECut    = pidChain.GetEntries(numuCCnonQECutString)


            # Now we work out the efficiencies: the numerator is the weighted number of selected events, 
            # and the denominator is the weighted number of total events in this sample.  The error is the
            # error on the mean of a binomially distributed variable: err = sqrt(p * (1-p) / n)

            # First we'll do the efficiency for selecting nue CCQE events
            nueCCQEEff   =  (numNueCCQECut*weightNueCCQEEvents)\
                           /(numNueCCQEEvents*weightNueCCQEEvents)
            errorNueCCQE = math.sqrt(nueCCQEEff * (1-nueCCQEEff) / (numNueCCQEEvents*weightNueCCQEEvents))
            histNueCCQE.SetBinContent(xBin, yBin, nueCCQEEff)
            histNueCCQE.SetBinError(xBin, yBin, errorNueCCQE)


            # And now the combined efficiency for nue CCQE and non-CCQE events
            nueCCAllEff   = (numNueCCQECut*weightNueCCQEEvents + numNueCCnonQECut*weightNueCCnonQEEvents)\
                           /(numNueCCQEEvents*weightNueCCQEEvents + numNueCCnonQEEvents*weightNueCCnonQEEvents)
            errorNueCCAll = math.sqrt(nueCCAllEff * (1-nueCCAllEff) / (numNueCCQEEvents*weightNueCCQEEvents + numNueCCnonQEEvents*weightNueCCnonQEEvents))
            histNueCCAll.SetBinContent(xBin, yBin, nueCCAllEff)
            histNueCCAll.SetBinError(xBin, yBin, errorNueCCAll)


            # Also the combined efficiency for selecting numu CCQE and non-QE events (want this to be low)
            numuCCAllEff = (numNumuCCQECut*weightNumuCCQEEvents + numNumuCCnonQECut*weightNumuCCnonQEEvents)\
                         /(numNumuCCQEEvents*weightNumuCCQEEvents + numNumuCCnonQEEvents*weightNumuCCnonQEEvents)
            errorNumuCCAll = math.sqrt(numuCCAllEff * (1-numuCCAllEff) / (numNumuCCQEEvents*weightNumuCCQEEvents + numNumuCCnonQEEvents*weightNumuCCnonQEEvents))
            histNumuCCAll.SetBinContent(xBin, yBin, 1-numuCCAllEff)
            histNumuCCAll.SetBinError(xBin, yBin, errorNumuCCAll)


            # And finally the combined efficiency with which nue NC and numu NC events are selected (also want this low)
            nuNCAllEff = (numNumuNCCut*weightNumuNCEvents + numNueNCCut*weightNueNCEvents)\
                         /(numNumuNCEvents*weightNumuNCEvents + numNueNCEvents*weightNueNCEvents)
            errorNuNCAll = math.sqrt(nuNCAllEff * (1-nuNCAllEff) / (numNumuNCEvents*weightNumuNCEvents + numNueNCEvents*weightNueNCEvents))
            histNuNCAll.SetBinContent(xBin, yBin, 1-nuNCAllEff)
            histNuNCAll.SetBinError(xBin, yBin, errorNuNCAll)


            # Now work out the purity of the nue CCQE selected sample
            nueCCAll  = ((numNueCCQECut*weightNueCCQEEvents) + (numNueCCnonQECut * weightNueCCnonQEEvents))
            nueCCQEAll = (numNueCCQECut*weightNueCCQEEvents)
            numuCCAll = ((numNumuCCQECut*weightNumuCCQEEvents) + (numNumuCCnonQECut * weightNumuCCnonQEEvents))
            nuNCAll   = (numNumuNCCut*weightNumuNCEvents + numNueNCCut*weightNueNCEvents)
            
            purity_cc = nueCCAll/(nueCCAll+numuCCAll+nuNCAll)
            purity_ccqe = nueCCQEAll/(nueCCAll+numuCCAll+nuNCAll)
            fom_cc = nueCCAll/math.sqrt(nueCCAll+numuCCAll+nuNCAll)
            fom_ccqe = nueCCQEAll/math.sqrt(nueCCAll+numuCCAll+nuNCAll)
            
            histNueCCAllPurity.SetBinContent(xBin, yBin, purity_cc)
            
            histNueCCAllEffPurity_CC.SetBinContent(xBin, yBin, nueCCAllEff*purity_cc)
            histNueCCAllEffPurity_CCQE.SetBinContent(xBin, yBin, nueCCQEEff*purity_ccqe)
            histFOM_CC.SetBinContent(xBin, yBin, fom_cc)
            histFOM_CCQE.SetBinContent(xBin, yBin, fom_ccqe)
            
            print "xCut = %s, yCut = %s, effPur_cc = %f, effPur_ccqe = %f, FOM_cc = %f, FOM_ccqe = %f" % (muVsElANNCut, CCVsNCANNCut, nueCCAllEff*purity_cc, nueCCQEEff*purity_ccqe, fom_cc, fom_ccqe)


    # Some formatting
    for hist in [histNueCCQE, histNueCCAll, histNumuCCAll, histNuNCAll, histNueCCAllPurity, histNueCCAllEffPurity_CC, histNueCCAllEffPurity_CCQE, histFOM_CC, histFOM_CCQE]:
        hist.GetYaxis().CenterTitle()
        hist.GetXaxis().CenterTitle()
        hist.GetXaxis().SetTitleSize(0.05)
        hist.GetYaxis().SetTitleSize(0.05)
        hist.GetXaxis().SetTitleOffset(1.30)
        hist.GetYaxis().SetTitleOffset(1.30)
    
    
    # Draw and save the plots
    can = ROOT.TCanvas("can","",800,600)
    histNueCCQE.SetMaximum(1)
    histNueCCQE.SetMinimum(0.0)
    histNueCCQE.Draw("colz")
    can.SaveAs("histNueCCQE.png")
    histNueCCAll.SetMaximum(1)
    histNueCCAll.SetMinimum(0.0)
    histNueCCAll.Draw("colz")
    can.SaveAs("histNueCCAll.png")
    histNumuCCAll.SetMaximum(1)
    histNumuCCAll.SetMinimum(0.9)
    histNumuCCAll.Draw("colz")
    can.SaveAs("histNumuCCAll.png")
    histNuNCAll.SetMaximum(1)
    histNuNCAll.SetMinimum(0.9)
    histNuNCAll.Draw("colz")
    can.SaveAs("histNuNCAll.png")
    histNueCCAllPurity.SetMaximum(1.0)
    histNueCCAllPurity.SetMinimum(0.0)
    histNueCCAllPurity.Draw("colz")
    can.SaveAs("histNueAllPurity.png")
    
    histNueCCAllEffPurity_CC.Draw("colz")
    binX_cc = ROOT.Long(0)
    binY_cc = ROOT.Long(0)
    binZ_cc = ROOT.Long(0)
    histNueCCAllEffPurity_CC.GetMaximumBin(binX_cc, binY_cc, binZ_cc)
    marker = ROOT.TMarker(hist.GetXaxis().GetBinUpEdge(binX_cc), hist.GetYaxis().GetBinUpEdge(binY_cc), 30)
    marker.SetMarkerSize(0.8)
    marker.SetMarkerColor(ROOT.kGreen)
    marker.Draw()
    can.SaveAs("histNueCCQEEffPurity_cc.png")
    can.SaveAs("histNueCCQEEffPurity_cc.pdf")
    can.SaveAs("histNueCCQEEffPurity_cc.root")
    can.SaveAs("histNueCCQEEffPurity_cc.C")
    
    histNueCCAllEffPurity_CCQE.Draw("colz")
    binX_ccqe = ROOT.Long(0)
    binY_ccqe = ROOT.Long(0)
    binZ_ccqe = ROOT.Long(0)
    histNueCCAllEffPurity_CCQE.GetMaximumBin(binX_ccqe, binY_ccqe, binZ_ccqe)
    marker = ROOT.TMarker(hist.GetXaxis().GetBinUpEdge(binX_ccqe), hist.GetYaxis().GetBinUpEdge(binY_ccqe), 30)
    marker.SetMarkerSize(0.8)
    marker.SetMarkerColor(ROOT.kGreen)
    marker.Draw()
    can.SaveAs("histNueCCQEEffPurity_ccqe.png")
    can.SaveAs("histNueCCQEEffPurity_ccqe.pdf")
    can.SaveAs("histNueCCQEEffPurity_ccqe.root")
    can.SaveAs("histNueCCQEEffPurity_ccqe.C")
    
    histFOM_CC.Draw("colz")
    binXFom_cc = ROOT.Long(0)
    binYFom_cc = ROOT.Long(0)
    binZFom_cc = ROOT.Long(0)
    histFOM_CC.GetMaximumBin(binXFom_cc, binYFom_cc, binZFom_cc)
    marker = ROOT.TMarker(hist.GetXaxis().GetBinUpEdge(binXFom_cc), hist.GetYaxis().GetBinUpEdge(binYFom_cc), 30)
    marker.SetMarkerSize(0.8)
    marker.SetMarkerColor(ROOT.kGreen)
    marker.Draw()
    can.SaveAs("histNueCCQEFOM_cc.png")
    can.SaveAs("histNueCCQEFOM_cc.pdf")
    can.SaveAs("histNueCCQEFOM_cc.root")
    can.SaveAs("histNueCCQEFOM_cc.C")
    
    histFOM_CCQE.Draw("colz")
    binXFom_ccqe = ROOT.Long(0)
    binYFom_ccqe = ROOT.Long(0)
    binZFom_ccqe = ROOT.Long(0)
    histFOM_CCQE.GetMaximumBin(binXFom_ccqe, binYFom_ccqe, binZFom_ccqe)
    marker = ROOT.TMarker(hist.GetXaxis().GetBinUpEdge(binXFom_ccqe), hist.GetYaxis().GetBinUpEdge(binYFom_ccqe), 30)
    marker.SetMarkerSize(0.8)
    marker.SetMarkerColor(ROOT.kGreen)
    marker.Draw()
    can.SaveAs("histNueCCQEFOM_ccqe.png")
    can.SaveAs("histNueCCQEFOM_ccqe.pdf")
    can.SaveAs("histNueCCQEFOM_ccqe.root")
    can.SaveAs("histNueCCQEFOM_ccqe.C")

    # Sweep out the efficiency vs. purity
    PlotEffVsPurity(histNueCCAll, histNueCCAllPurity, "efficiencyVsPurityNueCCAll.png")

    #print("Here, eff = %.04f and purity = %.04f" % (histNueCCAll.GetBinContent(binX, binY), histNueCCAllPurity.GetBinContent(binX, binY)))
    #print("Here, numu rejection eff = %.04f and NC rejection eff = %.04f" % (histNumuCCAll.GetBinContent(binX, binY), histNuNCAll.GetBinContent(binX, binY)))
    
    print("Efficiency times purity CC is maximised for nue/numu cut of %.02f, and CCQE/NC cut of %.02f" % (histNueCCAllEffPurity_CC.GetXaxis().GetBinUpEdge(binX_cc), histNueCCAllEffPurity_CC.GetYaxis().GetBinUpEdge(binY_cc)))
    print("Efficiency times purity CCQE is maximised for nue/numu cut of %.02f, and CCQE/NC cut of %.02f" % (histNueCCAllEffPurity_CCQE.GetXaxis().GetBinUpEdge(binX_ccqe), histNueCCAllEffPurity_CCQE.GetYaxis().GetBinUpEdge(binY_ccqe)))
    
    print("FOM CC is maximised for nue/numu cut of %.02f, and CCQE/NC cut of %.02f" % (histFOM_CC.GetXaxis().GetBinUpEdge(binXFom_cc), histFOM_CC.GetYaxis().GetBinUpEdge(binYFom_cc)))
    print("FOM CCQE is maximised for nue/numu cut of %.02f, and CCQE/NC cut of %.02f" % (histFOM_CCQE.GetXaxis().GetBinUpEdge(binXFom_ccqe), histFOM_CCQE.GetYaxis().GetBinUpEdge(binYFom_ccqe)))
    return

def PlotEffVsPurity(effHist, purHist, name):
    hist = ROOT.TH1D("hist",";Efficiency;Purity",100,0,1)
    hist.SetMinimum(0)
    hist.SetMaximum(1)
    for xBin, yBin in itertools.product(range(1, effHist.GetXaxis().GetNbins()+1), range(1, effHist.GetYaxis().GetNbins()+1)):
        eff = effHist.GetBinContent(xBin, yBin)
        pur = purHist.GetBinContent(xBin, yBin)
        bin = hist.GetXaxis().FindBin(eff)
        if(hist.GetBinContent(bin) < pur):
            hist.SetBinContent(bin, pur)
    can = ROOT.TCanvas("can","",1200, 800)
    hist.Draw()
    name = ''.join(name.split(".")[:-1])
    can.SaveAs(name + ".png")
    can.SaveAs(name + ".pdf")
    can.SaveAs(name + ".C")

if __name__ == '__main__':
    main()
