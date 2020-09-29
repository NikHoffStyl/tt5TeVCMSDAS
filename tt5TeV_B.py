from bamboo.analysismodules import NanoAODModule, NanoAODHistoModule, NanoAODSkimmerModule
from bamboo.analysisutils import configureRochesterCorrection
from bamboo import treefunctions as op
import os.path
import os
from bamboo.logging import getLogger
logger = getLogger(__name__)

from numpy import mean, sqrt, square
import numpy as np
from itertools import count, chain, product, repeat
from functools import partial
import math

import bamboo.plots
class Plot(bamboo.plots.Plot):
    def produceResults(self, bareResults, fbe):
        newresults = []
        if (any("__qcdScale" in h.GetName() for h in bareResults)): 

            print("starting syst unc envelope computation")
            hNom = next(h for h in bareResults if "__" not in h.GetName())
            prefix_qcd = f"{hNom.GetName()}__qcdScale"
            hVar_qcdScale = [ h for h in bareResults if h.GetName().startswith(prefix_qcd) ]
            if not all(hv.GetNcells() == hNom.GetNcells() for hv in hVar_qcdScale):
                logger.error("Variation histograms do not have the same binning as the nominal histogram")
            elif len(hVar_qcdScale) < 2:
                logger.error("At least two variations histograms must be provided")
            else: ## make an envelope from maximum deviations
                vars_cont = np.array([ [ hv.GetBinContent(i) for i in range(hv.GetNcells()) ] for hv in hVar_qcdScale ])
                hVar_up_qcd = hNom.Clone(f"{prefix}up")
                hVar_down_qcd = hNom.Clone(f"{prefix}down")
                
                for i,vl,vh in zip(count(), np.amin(vars_cont, axis=0), np.amax(vars_cont, axis=0)):
                    hVar_down_qcd.SetBinContent(i, vl)
                    hVar_up_qcd.SetBinContent(i, vh)
                    print("high:",vh," low:",vl)

                newresults+= hVar_up_qcd
                newresults+= hVar_down_qcd

        if (any("__PDFWeights" in h.GetName() for h in bareResults)):
           
            hNom = next(h for h in bareResults if "__" not in h.GetName())
            prefix_pdf = f"{hNom.GetName()}__PDFWeights"
            hVar_PDFWeights = [ h for h in bareResults if h.GetName().startswith(prefix_pdf) ]
            if not all(hv.GetNcells() == hNom.GetNcells() for hv in hVar_PDFWeights):
                 logger.error("Variation histograms do not have the same binning as the nominal histogram")
            elif len(hVar_PDFWeights) < 2:
                 logger.error("At least two variations histograms must be provided")
            else: ## make an envelope from rms of deviations
                print("pdf")
                pdf_vars=1
                vars_pdf = np.array([ [ hv.GetBinContent(i) for i in range(hv.GetNcells()-2) ] for hv in hVar_PDFWeights ])
                vars_alphaS = np.array([ [ hv.GetBinContent(i) for i in range(hv.GetNcells()-2,hv.GetNcells()) ] for hv in hVar_PDFWeights ])
                hVar_up_pdf = hNom.Clone(f"{prefix}up")
                hVar_down_pdf = hNom.Clone(f"{prefix}down")
                for i,rms in zip(count(), sqrt(mean(square(vars_pdf, axis=0)))):
                    #Xuan's recipe: sqrt(RMS^2 + ((alphas var 1 - alphas var 2)*0.75/2)^2 )
                    unc = sqrt(square(rms) + square((vars_alphaS[i,0]-alphaS[i,1])*0.75/2) ) 
                    print("pdf rms unc:",unc)
                    hVar_down_pdf.SetBinContent(i, 1-rms)
                    hVar_up_pdf.SetBinContent(i, 1+rms)

                newresults+= hVar_up_pdf
                newresults+= hVar_down_pdf
                
        return bareResults+newresults

import bamboo.treedecorators as btd
class ReadVariableVarWithSuffix(btd.ReadVariableVarWithSuffix):
    def getVarName(self, branchName, collgrpname=None):
        variNm = btd.normVarName(branchName[len(self.prefix):].lstrip(self.sep))
        return self.prefix, f"{self.systName}{variNm}" if variNm else self.nomName
class ReadMuScaleCorrection(btd.NanoSystematicVarSpec):
    def appliesTo(self, name):
        return name == "Muon"
    def nomName(self, name):
        return ""
    def getVarName(self, name, collgrpname=None):
        if name.startswith("corrected"):
            vari, v = name[len("corrected"):].split("_")
            return v, f"{self.systName}{btd.normVarName(vari)}" if vari else self.nomName(name)
class ReadElScaleCorrection(btd.NanoSystematicVarSpec):
    def appliesTo(self, name):
        return name == "Electron"
    def getVarName(self, name, collgrpname=None):
        if name.split("_")[0] in ("pt", "mass") and len(name.split("_")) >= 2:
            v = name.split("_")[0]
            vari = "_".join(name.split("_")[1:])
            return v, "nom" if vari == "nom" else f"{self.systName}{vari}"
    def nomName(self, name):
        return "nom"

_descr_5TeV_removedGroups = ["CaloMET_", "ChsMET_", "RawMET_", "TkMET_"]
_descr_5TeV_removedGroups_MC = _descr_5TeV_removedGroups + ["Generator_", "HTXS_"]
_descr_5TeV_removedCollections = ["nFatJet", "nIsoTrack", "nOtherPV", "nSV", "nSoftActivityJet", "nSubJet", "nTau", "nTrigObj"]
_descr_5TeV_removedCollections_MC = _descr_5TeV_removedCollections + ["nGenJetAK8", "nGenVisTau", "nSubGenJetAK8"]
description_nanov5_5TeV_data = btd.NanoAODDescription.get("v5", year="2018", isMC=False,
        removeGroups=_descr_5TeV_removedGroups, removeCollections=_descr_5TeV_removedCollections,
        systVariations=[
            ReadMuScaleCorrection("muScale"),
            ReadElScaleCorrection("elScale"), ## FIXME implement on the fly
            btd.ReadJetMETVar("Jet", "MET", jetsNomName="nom", metNomName="nom")

            ])
description_nanov5_5TeV_MC = btd.NanoAODDescription.get("v5", year="2018", isMC=True,
        removeGroups=_descr_5TeV_removedGroups_MC, removeCollections=_descr_5TeV_removedCollections_MC,
        addGroups=["Pileup_"], addCollections=["nLHEPart"],
        systVariations=[
            btd.nanoRochesterCalc,
            ReadVariableVarWithSuffix("puWeight"),
            ReadVariableVarWithSuffix("PrefireWeight"),
            ReadMuScaleCorrection("muScale"),
            ReadElScaleCorrection("elScale"), ## FIXME implement on the fly
            btd.ReadJetMETVar("Jet", "MET", jetsExclVars=["raw"], metNomName="nom", metExclVars=["raw", "nom"])
            ])

## no need to propagate the scale uncertainty into the scalefactor, always use nominal
import bamboo.scalefactors
binningVariables_nano_noScaleSyst = dict(bamboo.scalefactors.binningVariables_nano)
binningVariables_nano_noScaleSyst.update({ "Pt" : lambda obj : obj.pt.op.arg.wrapped.result[obj.idx] })
## for use with bamboo.scalefactors.get_scalefactor
## e.g. get_scalefactor("lepton", "Muon_RecoToLoose", sfLib=scalefactors_lepMVA, paramDefs=binningVariables_nano_noScaleSyst, systName="muLoose")
scalefactors_lepMVA = {
    f"{lFlav}_{ratioSel}" : os.path.join(os.path.dirname(__file__), "data", f"{lFlav}_{ratioSel}SF.json")
        for ratioSel in ("RecoToLoose", "LooseToTight")
        for lFlav in ("Electron", "Muon")
    }

class Nano5TeVBase(NanoAODModule):
    """ Base module for postprocessed 5TeV NanoAODv5 samples """
    def prepareTree(self, tree, sample=None, sampleCfg=None):
        ## Decorate the tree
        isMC = self.isMC(sample)
        tree,noSel,be,lumiArgs = super(Nano5TeVBase, self).prepareTree(tree, sample=sample, sampleCfg=sampleCfg,
                                description=(description_nanov5_5TeV_MC if isMC else description_nanov5_5TeV_data))
        #rochesterFile = None

        if sampleCfg.get("alt-syst"):
            noSel = noSel.refine("withoutsyst", autoSyst=False)

        if isMC:
                                    
            #psWeight systematics for varying muR and muF
            qcdScale=op.systematic(op.c_float(1.), name="qcdScale", **{f"qcdScale{i:d}" : tree.LHEScaleWeight[i] for i in (0, 1, 3, 5, 7, 8)})
            PDFWeights=op.systematic(op.c_float(1.), name="PDFWeights", **{f"PDFWeights{i:d}" : tree.LHEPdfWeight[i] for i in range(0, 102)})
            psISR = op.systematic(op.c_float(1.), name="psISR", up=tree.PSWeight[2], down=tree.PSWeight[0])
            psFSR = op.systematic(op.c_float(1.), name="psFSR", up=tree.PSWeight[3], down=tree.PSWeight[1])

            noSel = noSel.refine("mcWeight", weight=[ tree.genWeight , tree.puWeight , tree.PrefireWeight , psISR , psFSR, qcdScale, PDFWeights])
            
            #configure rochester correction and variations
            #try:
            #    configureRochesterCorrection(tree._Muon, rochesterFile, isMC=isMC, backend=be, uName=sample)
            #except Exception as ex:
            #    logger.exception("Problem while configuring the Rochester correction")

        return tree,noSel,be,lumiArgs

class Nano5TeVHistoModule(Nano5TeVBase, NanoAODHistoModule):
    pass

class Nano5TeVAnalyzer(Nano5TeVHistoModule):
   
    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        config_noaltsyst = dict(config) ## warning: shallow copy
        config_noaltsyst["samples"] = { smpNm: smpCfg for smpNm, smpCfg in config["samples"].items() if "alt-syst" not in smpCfg }
        super(Nano5TeVAnalyzer, self).postProcess(taskList, config=config_noaltsyst, workdir=workdir, resultsdir=resultsdir)


    def definePlots(self, t, noSel, sample=None, sampleCfg=None):
        from bamboo.plots import Plot, CutFlowReport, SummedPlot
        from bamboo.plots import EquidistantBinning as EqB
        from bamboo import treefunctions as op

        isMC = self.isMC(sample)
        trigCut, trigWeight = None, None
        if isMC:
            trigCut = op.OR(t.HLT.HIEle20_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, t.HLT.HIL3DoubleMu0, t.HLT.HIL3Mu20, t.HLT.HIEle20_WPLoose_Gsf)
            trigWeight = op.switch(op.OR(t.HLT.HIEle20_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, t.HLT.HIL3DoubleMu0), op.c_float(1.),
                    op.switch(t.HLT.HIL3Mu20, op.c_float(306.913/308.545), op.c_float(264.410/308.545))) ## FIXME these are wrong - you will get the final values from team A
        else:
            ## trigger order: dielectron, dimuon or single muon, single electron
            pd = sample.split("_")[0]
            if pd == "SingleMuon":
                trigCut = op.AND(op.NOT(t.HLT.HIEle20_Ele12_CaloIdL_TrackIdL_IsoVL_DZ),
                    op.OR(t.HLT.HIL3DoubleMu0, t.HLT.HIL3Mu20))
            elif pd == "HighEGJet":
                trigCut = op.OR(t.HLT.HIEle20_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,
                    op.AND(op.NOT(op.OR(t.HLT.HIL3DoubleMu0, t.HLT.HIL3Mu20)),
                        t.HLT.HIEle20_WPLoose_Gsf))
        noSel = noSel.refine("trig", cut=trigCut, weight=trigWeight)

        plots = []

        def isGoodElectron(el, ptCut=10.):
            return op.AND(
                el.pt > ptCut,
                op.abs(el.eta) < 2.5,
                el.lostHits == 0, ## do you want this?
                op.abs(el.sip3d) < 8.,
                op.abs(el.dxy) < .05,
                op.abs(el.dz ) < .1,
                el.miniPFRelIso_all < 0.085,
                el.mvaTTH > 0.125,
                op.NOT(op.AND(el.jet.isValid, op.OR(el.jet.btagDeepB > .1522, el.jet.btagDeepB <= -999.)))
                )
        def isGoodMuon(mu, ptCut=10.):
            return op.AND(
                mu.pt > ptCut,
                op.abs(mu.eta) < 2.4,
                mu.mediumPromptId,
                op.abs(mu.sip3d) < 8.,
                op.abs(mu.dxy) < .05,
                op.abs(mu.dz ) < .1,
                mu.miniPFRelIso_all < 0.325,
                mu.mvaTTH > 0.55,
                op.NOT(op.AND(mu.jet.isValid, op.OR(mu.jet.btagDeepB > .1522, mu.jet.btagDeepB <= -999.)))
                )

        goodLeptons = {
            "el" : op.select(t.Electron, partial(isGoodElectron, ptCut=15.)),
            "mu" : op.select(t.Muon, partial(isGoodMuon, ptCut=15.))
            }
        plots += [
            Plot.make1D("trig_nLeptons15", op.rng_len(goodLeptons["el"])+op.rng_len(goodLeptons["mu"]), noSel, EqB(15, 0., 15.)),
            Plot.make1D("trig_nEl15", op.rng_len(goodLeptons["el"]), noSel, EqB(15, 0., 15.)),
            Plot.make1D("trig_nMu15", op.rng_len(goodLeptons["mu"]), noSel, EqB(15, 0., 15.)) 
            ]
        from bamboo.scalefactors import get_scalefactor
        sf_loose = {
            "mu": get_scalefactor("lepton", "Muon_RecoToLoose", sfLib=scalefactors_lepMVA, paramDefs=binningVariables_nano_noScaleSyst, systName="muLoose"),
            "el": get_scalefactor("lepton", "Electron_RecoToLoose", sfLib=scalefactors_lepMVA, paramDefs=binningVariables_nano_noScaleSyst, systName="elLoose")
            }
        sf_tight = {
            "mu": get_scalefactor("lepton", "Muon_LooseToTight", sfLib=scalefactors_lepMVA, paramDefs=binningVariables_nano_noScaleSyst, systName="muTight"),
            "el": get_scalefactor("lepton", "Electron_LooseToTight", sfLib=scalefactors_lepMVA, paramDefs=binningVariables_nano_noScaleSyst, systName="elTight")
            }

        nGoodLeptons = op.rng_len(goodLeptons["el"])+op.rng_len(goodLeptons["mu"])
        hasTwoGoodLeptons = noSel.refine("has2Lep", cut=(nGoodLeptons > 1)) # avoid overlap with 1l
        jets = op.sort(op.select(t.Jet, lambda j : op.AND(
            j.pt > 25.,
            op.abs(j.eta) < 2.4,
            j.jetId & 0x2,
            op.AND(
                op.NOT(op.rng_any(goodLeptons["el"], lambda l : op.deltaR(l.p4, j.p4) < 0.4)),
                op.NOT(op.rng_any(goodLeptons["mu"], lambda l : op.deltaR(l.p4, j.p4) < 0.4)))
            )), lambda j : -j.pt)
        ## WP: see https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
        loosebjets = op.select(jets, lambda j : j.btagDeepB > 0.1522)
        mediumbjets = op.select(jets, lambda j : j.btagDeepB > 0.4941)
        for fl1,fl2 in product(*repeat(goodLeptons.keys(), 2)):
            dilepSel = lambda l1,l2 : op.AND(
                    l1.charge != l2.charge,
                    (l1.p4+l2.p4).M() > 12.
                    )
            if fl1 == fl2:
                lGood = op.sort(goodLeptons[fl1], lambda l : -l.pt)
                dilep = op.combine(lGood, N=2, pred=dilepSel)
            else:
                l1Good = op.sort(goodLeptons[fl1], lambda l : -l.pt)
                l2Good = op.sort(goodLeptons[fl2], lambda l : -l.pt)
                dilep = op.combine((l1Good, l2Good), pred=dilepSel)
            ll = dilep[0]
            hasDilep = hasTwoGoodLeptons.refine(f"hasDilep{fl1}{fl2}", cut=(op.rng_len(dilep) > 0, ll[0].pt > 25.),
                    weight=([ sf_loose[fl1](ll[0]), sf_loose[fl2](ll[1]), sf_tight[fl1](ll[0]), sf_tight[fl2](ll[1]) ] if isMC else None))
            plots += [
                Plot.make1D(f"dilepton_{fl1}{fl2}_Mll", (ll[0].p4+ll[1].p4).M(), hasDilep, EqB(50, 70, 120.), title="Dilepton mass"),
                ]
#            for il,ifl in enumerate((fl1, fl2)):
##                plots += [
#                    Plot.make1D(f"dilepton_{fl1}{fl2}_L{il:d}PT", ll[il].pt, hasDilep, EqB(50, 0., 100.), title=f"Lepton {il:d} PT"),
#                    Plot.make1D(f"dilepton_{fl1}{fl2}_L{il:d}ETA", ll[il].eta, hasDilep, EqB(50, -2.5, 2.5), title=f"Lepton {il:d} ETA"),
#                    ]
#            plots += [
#                Plot.make1D(f"dilepton_{fl1}{fl2}_nJets", op.rng_len(jets), hasDilep, EqB(15, 0, 15.), title="Jet multiplicity"),
#                Plot.make1D(f"dilepton_{fl1}{fl2}_nLooseBJets", op.rng_len(loosebjets), hasDilep, EqB(15, 0, 15.), title="Loose b-jet multiplicity"),
#                Plot.make1D(f"dilepton_{fl1}{fl2}_nMediumBJets", op.rng_len(mediumbjets), hasDilep, EqB(15, 0, 15.), title="Medium b-jet multiplicity")
#                ]

        return plots
