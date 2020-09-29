from functools import partial
from bamboo.analysismodules import NanoAODModule, NanoAODHistoModule, NanoAODSkimmerModule
import os.path

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
            # return name.split("_")[0], "_".join(name.split("_")[1:])
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
        tree,noSel,be,lumiArgs = super(Nano5TeVBase, self).prepareTree(tree, sample=sample, sampleCfg=sampleCfg,
                description=(description_nanov5_5TeV_MC if self.isMC(sample) else description_nanov5_5TeV_data))
        return tree,noSel,be,lumiArgs

class Nano5TeVHistoModule(Nano5TeVBase, NanoAODHistoModule):
    pass

from itertools import chain, product, repeat
import ROOT
class DileptonBGroup(Nano5TeVHistoModule):
    
    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        config_noaltsyst = dict(config) ## warning: shallow copy
        config_noaltsyst["samples"] = { smpNm: smpCfg for smpNm, smpCfg in config["samples"].items() if "alt-syst" not in smpCfg }
        super(DileptonBGroup, self).postProcess(taskList, config=config_noaltsyst, workdir=workdir, resultsdir=resultsdir)
        from bamboo.analysisutils import loadPlotIt
        p_config, samples, plots, systematics, legend = loadPlotIt(config, self.plotList, eras=self.args.eras[1], workdir=workdir, resultsdir=resultsdir, readCounters=self.readCounters, vetoFileAttributes=self.__class__.CustomSampleAttributes, plotDefaults=self.plotDefaults)

        acceptance = {}
        effic = {}
        acceptance1Val = {}
        effic1Val = {}

        set1Dxsec = True

        xsecFile = ROOT.TFile("xsec.root", "RECREATE")

        from plotit.plotit import Stack
        from bamboo.root import gbl

        for plot in plots:
            print(plot.name + " sample name"+ str([smp.cfg.name for smp in samples]))
            obsHists = [smp.getHist(plot) for smp in samples if smp.cfg.type == "DATA"]
            print(plot.name + " Data sample name"+ str([smp.cfg.name for smp in samples if smp.cfg.type == "DATA"]))
            obsStack = Stack(obsHists) 
            bkgHists = [smp.getHist(plot) for smp in samples if (smp.cfg.type == "MC" and ("TT" not in smp.cfg.name) and ("down" not in smp.cfg.name) and ("up" not in smp.cfg.name) )]
            print(plot.name + " Bkgsample name"+ str([smp.cfg.name for smp in samples if (smp.cfg.type == "MC" and ("TT" not in smp.cfg.name) and ("down" not in smp.cfg.name) and ("up" not in smp.cfg.name) )]))
            bkgStack = Stack(bkgHists)
            sigPlot = [smp.getHist(plot) for smp in samples if smp.cfg.name == "TT"]
            print(plot.name + " TT sample name"+ str([smp.cfg.name for smp in samples if smp.cfg.name == "TT"]))

            if "nTotalJets" in plot.name:
                for smp in samples:
                    if smp.cfg.type != "MC": continue
                    accden = smp.getHist(plot).obj.Clone(plot.name+"accden")
                    effden = smp.getHist(plot).obj.Clone(plot.name+"effden")

            if "nJets" in plot.name:
                for smp in samples:
                    if smp.cfg.type != "MC": continue
                    if "elel" in plot.name: branchRat = 0.0123
                    elif "elmu" in plot.name: branchRat = 0.0253
                    elif "mumu" in plot.name: branchRat = 0.0130

                    accTmp = smp.getHist(plot).obj.Clone(plot.name+"acc")
                    accTmpVal = accTmp.Integral()
                    accdenVal = (accden.Integral())
                    acceptance1Val[plot.name] = accTmpVal / (accdenVal * branchRat)
                    effic1Val[plot.name] = accTmpVal / accdenVal
                    print(plot.name+ " A ==  " + str(acceptance1Val[plot.name]))

                    print(plot.name+ " epsilon ==  " + str(effic1Val[plot.name]))
                    accden.Scale(branchRat)
                    accden.Print()
                    accTmp.Divide(accden)
                    accTmp.Print()
                    acceptance[plot.name] = accTmp 

                    effTmp  = smp.getHist(plot).obj.Clone(plot.name+"eff")
                    effTmp.Divide(effden)
                    effic[plot.name] = effTmp

                bkgMerge = bkgHists[0].obj.Clone(plot.name+"merge")
                for nbh, bkghist in enumerate(bkgHists):
                    if nbh !=0:
                        print(plot.name)
                        bkgMerge.Add(bkghist.obj)
                bkgMergeEntries = bkgMerge.Integral()

                obsMerge = obsHists[0].obj.Clone(plot.name+"merge")
                for nh, obshist in enumerate(obsHists):
                    if nh !=0:
                        obsMerge.Add(obshist.obj)
                obsMergeEntries = obsMerge.Integral()

                sigtt = sigPlot[0].obj.Clone(plot.name+"sigtt")
                sigttEntries = sigtt.Integral()
                xsec_ttbar1Val =( (obsMergeEntries - bkgMergeEntries) / sigttEntries)*68.9
                xsec_ttbar1ValError = xsec_ttbar1Val /10  # FIXME: this may not be needed, error should be propagated-used from up-down plots
                print("  >>>>>>>>  " + plot.name+ " xSec ==  " + str(xsec_ttbar1Val))


                if set1Dxsec:
                    xsec_ttbar = ROOT.TH1D("plot.name", "xsec; ; xsec", 60., 80., 100)
                    xsec_ttbar.Fill(xsec_ttbar1Val)
                    xsec_ttbar.SetBinError(xsec_ttbar.FindBin(h.GetBinCenter(1)), xsec_ttbar1ValError)
                    xsec_ttbar.Print()
                else:
                    xsec_ttbar = obsMerge.Clone(plot.name+"xsec")
                    xsec_ttbar.Add(bkgMerge, -1)
                    xsec_ttbar.Print()
                    sigtt.Print()
                    xsec_ttbar.Divide(sigtt)
                cv = gbl.TCanvas(f"c{plot.name}")
                cv.cd(1)
                xsec_ttbar.Draw()
                cv.Update()
                xsecFile.Write(xsec_ttbar)
                print(os.path.join(resultsdir, f"{plot.name}ttbar.png"))
                cv.SaveAs(os.path.join(resultsdir, f"{plot.name}ttbar.pdf"))

        import IPython
        plots = { p.name: p for p in plots } # for convenience: turn list into dictionary
        IPython.embed()

    def definePlots(self, t, noSel, sample=None, sampleCfg=None):
        from bamboo.plots import Plot, CutFlowReport, SummedPlot
        from bamboo.plots import EquidistantBinning as EqB
        from bamboo import treefunctions as op

        isMC = self.isMC(sample)

        if sampleCfg.get("alt-syst"):
            noSel = noSel.refine("withoutsyst", autoSyst=False)

        plots = []

        trigCut, trigWeight = None, None
        if isMC:
            muR = op.systematic(op.c_float(1.), name="muR", up=t.PSWeight[2], down=t.PSWeight[0])
            muF = op.systematic(op.c_float(1.), name="muF", up=t.PSWeight[3], down=t.PSWeight[1])

            noSel = noSel.refine("mcWeight", weight=[ t.genWeight , t.puWeight , t.PrefireWeight , muR , muF])

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
        #plots += [Plot.make1D("nTotalEvents", op.rng_len([1]), noSel , EqB(1, 0, 1.), title="nTotalEvents")]
        plots.append(Plot.make1D("nTotalJets", op.rng_len(t.Jet), noSel, EqB(15, 0, 15.), title="Initial Jet multiplicity"))
        #noSel = noSel.refine("trig", cut=op.OR(t.HLT.HIL3DoubleMu0, t.HLT.HIEle20_Ele12_CaloIdL_TrackIdL_IsoVL_DZ))

        # plots = []
        goodLeptons = {
            "el" : op.select(t.Electron, lambda el : op.AND(el.pt > 15. , op.abs(el.p4.Eta()) < 2.5)),  # op.select(t.Electron, partial(isGoodElectron, ptCut=15.)),
            "mu" : op.select(t.Muon, lambda mu : mu.pt > 20.)  # op.select(t.Muon, partial(isGoodMuon, ptCut=15.))
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
            for il,ifl in enumerate((fl1, fl2)):
                plots += [
                    Plot.make1D(f"dilepton_{fl1}{fl2}_L{il:d}PT", ll[il].pt, hasDilep, EqB(50, 0., 100.), title=f"Lepton {il:d} PT"),
                    Plot.make1D(f"dilepton_{fl1}{fl2}_L{il:d}ETA", ll[il].eta, hasDilep, EqB(50, -2.5, 2.5), title=f"Lepton {il:d} ETA"),
                ]
            plots += [
                Plot.make1D(f"dilepton_{fl1}{fl2}_nJets", op.rng_len(jets), hasDilep, EqB(15, 0, 15.), title="Jet multiplicity"),
                Plot.make1D(f"dilepton_{fl1}{fl2}_nLooseBJets", op.rng_len(loosebjets), hasDilep, EqB(15, 0, 15.), title="Loose b-jet multiplicity"),
                Plot.make1D(f"dilepton_{fl1}{fl2}_nMediumBJets", op.rng_len(mediumbjets), hasDilep, EqB(15, 0, 15.), title="Medium b-jet multiplicity"),
                #Plot.make1D(f"dilepton_{fl1}{fl2}_nSelectedEvents", 1, hasDilep, EqB(1, 0, 1.), title="nSelectedEvents")
            ]


        #muons = op.select(t.Muon, lambda mu : mu.pt > 20.)
        #twoMuSel = noSel.refine("twoMuons", cut=[ op.rng_len(muons) > 1 ])

        #electrons = op.select(t.Electron, lambda el : op.AND(el.pt > 15. , op.abs(el.p4.Eta()) < 2.5))
        #twoElSel = noSel.refine("twoElectrons", cut=[ op.rng_len(electrons) > 1 ])

        #oselmu = op.combine((electrons, muons))
        #leptons = oselmu[0] 
        #twoLepSel = noSel.refine("twoLeptons", cut=[ op.rng_len(electrons) == 1 , op.rng_len(muons) == 1 ])

        #jets = op.select(t.Jet, lambda j : j.pt > 30.)

        #bjets = op.select(jets, lambda j : j.btagDeepFlavB > 0.2217)

        #plots.append(Plot.make1D("dimu_M",
        #    op.invariant_mass(muons[0].p4, muons[1].p4), twoMuSel, EqB(100, 20., 120.),
        #    title="Dimuon invariant mass", plotopts={"show-overflow":False,
        #    "legend-position": [0.2, 0.6, 0.5, 0.9]}))
        #plots.append(Plot.make1D("diel_M",
        #    op.invariant_mass(electrons[0].p4, electrons[1].p4), twoElSel, EqB(100, 20., 120.),
        #    title="Dielectron invariant mass", plotopts={"show-overflow":False,
        #    "legend-position": [0.2, 0.6, 0.5, 0.9]}))
        #plots.append(Plot.make1D("dilep_M",
        #    op.invariant_mass(leptons[0].p4, leptons[1].p4) , twoLepSel, EqB(100, 20., 120.),
        #    title="Dimuon invariant mass", plotopts={"show-overflow":False,
        #    "legend-position": [0.2, 0.6, 0.5, 0.9]}))
        #plots.append(SummedPlot("Mjj", plots, title="m(jj)"))

        #plots.append(Plot.make1D("nJets_dimu",op.rng_len(jets), twoMuSel, EqB(10, -0.5, 9.5),
        #    title="Jet multiplicity", plotopts={"show-overflow":False,
        #    "legend-position": [0.2, 0.6, 0.5, 0.9]}))

        #plots.append(Plot.make1D("nBJets_dimu",op.rng_len(bjets), twoMuSel, EqB(10, -0.5, 9.5),
        #    title="Jet multiplicity", plotopts={"show-overflow":False,
        #    "legend-position": [0.2, 0.6, 0.5, 0.9]}))

        #plots.append(Plot.make1D("nJets_diel",op.rng_len(jets), twoElSel, EqB(10, -0.5, 9.5),
        #    title="Jet multiplicity", plotopts={"show-overflow":False,
        #    "legend-position": [0.2, 0.6, 0.5, 0.9]}))

        #plots.append(Plot.make1D("nBJets_diel",op.rng_len(bjets), twoElSel, EqB(10, -0.5, 9.5),
        #    title="Jet multiplicity", plotopts={"show-overflow":False,
        #    "legend-position": [0.2, 0.6, 0.5, 0.9]}))

        #plots.append(Plot.make1D("nJets_elmu",op.rng_len(jets), twoLepSel, EqB(10, -0.5, 9.5),
        #    title="Jet multiplicity", plotopts={"show-overflow":False,
        #    "legend-position": [0.2, 0.6, 0.5, 0.9]}))

        #plots.append(Plot.make1D("nBJets_elmu",op.rng_len(bjets), twoLepSel, EqB(10, -0.5, 9.5),
        #    title="Jet multiplicity", plotopts={"show-overflow":False,
        #    "legend-position": [0.2, 0.6, 0.5, 0.9]}))

        return plots
