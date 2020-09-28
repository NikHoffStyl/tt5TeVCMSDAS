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
            return name.split("_")[0], "_".join(name.split("_")[1:])
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

class DileptonBGroup(Nano5TeVHistoModule):
    
    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        config_noaltsyst = dict(config) ## warning: shallow copy
        config_noaltsyst["samples"] = { smpNm: smpCfg for smpNm, smpCfg in config["samples"].items() if "alt-syst" not in smpCfg }
        super(DileptonBGroup, self).postProcess(taskList, config=config_noaltsyst, workdir=workdir, resultsdir=resultsdir)

    def definePlots(self, t, noSel, sample=None, sampleCfg=None):
        from bamboo.plots import Plot, CutFlowReport, SummedPlot
        from bamboo.plots import EquidistantBinning as EqB
        from bamboo import treefunctions as op

        isMC = self.isMC(sample)

        if sampleCfg.get("alt-syst"):
            noSel = noSel.refine("withoutsyst", autoSyst=False)

        if isMC:
            muR = op.systematic(op.c_float(1.), name="muR", up=t.PSWeight[2], down=t.PSWeight[0])
            muF = op.systematic(op.c_float(1.), name="muF", up=t.PSWeight[3], down=t.PSWeight[1])

            noSel = noSel.refine("mcWeight", weight=[ t.genWeight , t.puWeight , t.PrefireWeight , muR , muF])
        noSel = noSel.refine("trig", cut=op.OR(t.HLT.HIL3DoubleMu0, t.HLT.HIEle20_Ele12_CaloIdL_TrackIdL_IsoVL_DZ))

        plots = []

        muons = op.select(t.Muon, lambda mu : mu.pt > 20.)
        twoMuSel = noSel.refine("twoMuons", cut=[ op.rng_len(muons) > 1 ])

        electrons = op.select(t.Electron, lambda el : op.AND(el.pt > 15. , op.abs(el.p4.Eta()) < 2.5))
        twoElSel = noSel.refine("twoElectrons", cut=[ op.rng_len(electrons) > 1 ])

        oselmu = op.combine((electrons, muons))
        leptons = oselmu[0] 
        #leptons = muons + electrons
        twoLepSel = noSel.refine("twoLeptons", cut=[ op.rng_len(electrons) == 1 , op.rng_len(muons) == 1 ])

        jets = op.select(t.Jet, lambda j : j.pt > 30.)

        bJets = op.select(jets, lambda j : j.btagDeepFlavB > 0.2217)

        plots.append(Plot.make1D("dimu_M",
            op.invariant_mass(muons[0].p4, muons[1].p4), twoMuSel, EqB(100, 20., 120.),
            title="Dimuon invariant mass", plotopts={"show-overflow":False,
            "legend-position": [0.2, 0.6, 0.5, 0.9]}))
        plots.append(Plot.make1D("diel_M",
            op.invariant_mass(electrons[0].p4, electrons[1].p4), twoElSel, EqB(100, 20., 120.),
            title="Dielectron invariant mass", plotopts={"show-overflow":False,
            "legend-position": [0.2, 0.6, 0.5, 0.9]}))
        plots.append(Plot.make1D("dilep_M",
            op.invariant_mass(leptons[0].p4, leptons[1].p4) , twoLepSel, EqB(100, 20., 120.),
            title="Dimuon invariant mass", plotopts={"show-overflow":False,
            "legend-position": [0.2, 0.6, 0.5, 0.9]}))
        #plots.append(SummedPlot("Mjj", plots, title="m(jj)"))

        print("LALALALAL")

        plots.append(Plot.make1D("nJets",op.rng_len(jets), twoMuSel, EqB(10, -0.5, 9.5),
            title="Jet multiplicity", plotopts={"show-overflow":False,
            "legend-position": [0.2, 0.6, 0.5, 0.9]}))

        return plots
