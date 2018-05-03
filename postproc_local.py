#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor

from PhysicsTools.NanoAODTools.postprocessing.analysis.higgs.vhbb.VHbbProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jecUncertainties import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetRecalib import *
from  PhysicsTools.NanoAODTools.postprocessing.modules.jme.mht import *
#from  PhysicsTools.NanoAODTools.postprocessing.examples.puWeightProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.muonScaleResProducer import *
import argparse

print "args are: ",sys.argv

isMC = True
era = "2017"
#if len(sys.argv) > 1:
#    if sys.argv[1] == "0":
#        isMC = False
#if len(sys.argv) > 2:
#    era = sys.argv[2]
#if era!="2016" and era!="2017":
#    print "Run era must be 2016 or 2017, exiting.."
#    sys.exit(1)
#btagger = "deepcsv"
#if era == "2016":
#    btagger = "cmva"
dataRun = ""
#if len(sys.argv) > 3:
#    dataRun = sys.argv[3]
parser = argparse.ArgumentParser("")
parser.add_argument('-isMC', '--isMC', type=int, default=1, help="")
parser.add_argument('-jobNum', '--jobNum', type=int, default=1, help="")
parser.add_argument('-era', '--era', type=str, default="2017", help="")
parser.add_argument('-dataRun', '--dataRun', type=str, default="X", help="")
args = parser.parse_args()
print "args = ",args
isMC = args.isMC
era = args.era
dataRun = args.dataRun

print "isMC = ",isMC,"era = ",era, "dataRun = ",dataRun

#files=["root://cms-xrd-global.cern.ch//store/user/arizzi/NanoTestProd006/QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer17MiniAOD-92X-NanoCrabProd006/171006_144159/0000/nanolzma_1.root"]
#files=["lzma_1.root"]
#files=["root://cms-xrd-global.cern.ch://store/user/arizzi/NanoCrabProdXmas/WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_pythia8/NanoCrabXmasRunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asym___heIV_v6-v1__215/171221_111154/0000/nano_1.root"]
#files=["root://cms-xrd-global.cern.ch://store/user/arizzi/Nano01Fall17/ZH_HToBB_ZToLL_M120_13TeV_powheg_pythia8/RunIIFall17MiniAOD-94X-Nano01Fall17/180205_183348/0000/test94X_NANO_1.root"]
#files=["root://cms-xrd-global.cern.ch://store/user/arizzi/Nano01_17Nov17/SingleMuon/RunII2017ReReco17Nov17-94X-Nano01/180205_181424/0000/test_data_94X_NANO_1.root"]
#files=["root://cms-xrd-global.cern.ch://store/user/arizzi/Nano01Fall17/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAOD-94X-Nano01Fall17/180214_110649/0000/test94X_NANO_1.root"]
files=["root://cms-xrd-global.cern.ch://store/mc/RunIIFall17NanoAOD/WminusH_HToBB_WToLNu_M125_13TeV_powheg_herwigpp/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/164DBD28-A842-E811-80B3-0CC47A78A418.root"]
#files=["root://cms-xrd-global.cern.ch://store/mc/RunIISummer16NanoAOD/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/NANOAODSIM/PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/00000/E894B8DE-1816-E811-B138-008CFA1660F8.root"]
#files=["root://cmsxrootd.fnal.gov://store/data/Run2016H/SingleElectron/NANOAOD/05Feb2018_ver2-v2/00000/1051F540-4B0B-E811-BB45-002590D9D990.root"]
#files=["root://cms-xrd-global.cern.ch://store/data/Run2016D/SingleElectron/NANOAOD/05Feb2018-v2/00000/2827562B-3F0B-E811-90AE-E0071B7AE500.root"]
#files=["root://cms-xrd-global.cern.ch://store/user/arizzi/NanoTestProd004/WminusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8/NanoCrabProd004/171002_120520/0000/lzma_1.root"]
#files=["root://cms-xrd-global.cern.ch://store/user/arizzi/NanoTestProd004/WplusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8/NanoCrabProd004/171002_120552/0000/lzma_1.root"]
#files=["root://cms-xrd-global.cern.ch://store/user/arizzi/NanoTestProd004/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/NanoCrabProd004/171002_120644/0000/lzma_1.root"]
#files=["root://cms-xrd-global.cern.ch://store/user/arizzi/NanoTestProd004/ggZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/NanoCrabProd004/171002_122256/0000/lzma_1.root"]
#files=["root://cms-xrd-global.cern.ch://store/user/arizzi/NanoTestProd004/ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8/NanoCrabProd004/171002_122221/0000/lzma_1.root"]
filesTTbar= [
'root://cms-xrd-global.cern.ch//store/group/cmst3/group/nanoAOD/NanoTestProd006/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer17MiniAOD-92X-NanoCrabProd006/171006_155430/0000/nanolzma_1.root',
]

#selection='''(Sum$(Electron_pt > 20 && Electron_mvaSpring16GP_WP90) >= 2  ||
# Sum$(Muon_pt > 20) >= 2 ||
# Sum$(Electron_pt > 20 && Electron_mvaSpring16GP_WP80) >= 1   ||
# Sum$(Muon_pt > 20 && Muon_tightId) >= 1 ||
# (Sum$(Muon_pt > 20) == 0 && Sum$(Electron_pt > 20 && Electron_mvaSpring16GP_WP90) == 0 && MET_pt > 150 ) ) 
# &&  Sum$((abs(Jet_eta)<2.5 && Jet_pt > 20 && Jet_jetId)) >= 2 && Entry$ < 1000 
#'''
selection='''(Sum$(Electron_pt > 20 && Electron_mvaSpring16GP_WP90) >= 2  ||
 Sum$(Muon_pt > 20) >= 2 ||
 Sum$(Electron_pt > 20 && Electron_mvaSpring16GP_WP80) >= 1   ||
 Sum$(Muon_pt > 20 && Muon_tightId) >= 1 ||
 (Sum$(Muon_pt > 20) == 0 && Sum$(Electron_pt > 20 && Electron_mvaSpring16GP_WP90) == 0 && MET_pt > 150 ) ) 
 &&   ((Sum$((abs(Jet_eta)<2.5 && Jet_pt > 20 && Jet_jetId)) >= 2)||(Sum$((abs(FatJet_eta)<2.5 && FatJet_pt > 200 && FatJet_jetId)) >= 1)) && Entry$ < 1000 
'''

if era == "2017":
    selection = selection.replace("Electron_mvaSpring16GP","Electron_mvaFall17Iso")

selectionALL='''Sum$(Electron_pt > 20 && Electron_mvaSpring16GP_WP90) >= 2  ||
 Sum$(Electron_pt > 20 && Electron_mvaSpring16GP_WP80) >= 1   ||
 Sum$(Jet_pt > 40 && Jet_jetId) >= 4   || 
Sum$(Jet_pt *(abs(Jet_eta)<2.5 && Jet_pt > 20 && Jet_jetId)) > 160  || 
MET_pt > 100  || Sum$(Muon_pt > 20 && Muon_tightId) >= 1
'''
mhtVHbb = lambda : mhtProducer( lambda j : j.pt > 30,
                            lambda mu : mu.pt > 5 and mu.pfRelIso04_all < 0.4,
                            lambda el : el.pt > 5 and el.pfRelIso03_all < 0.4 )

#p=PostProcessor(".",files,selection.replace('\n',' '),"keep_and_drop.txt",[jetmetUncertainties(),vhbb()],provenance=True)
##p=PostProcessor(".",files,selection.replace('\n',' '),"keep_and_drop.txt",[jetmetUncertaintiesAll(),btagSFProducer("cmva"),vhbb()],provenance=True)
#p=PostProcessor(".",files,selection.replace('\n',' '),"keep_and_drop.txt",[jecUncertAll_cppOut(),jetmetUncertainties(),btagSFProducer("cmva"),vhbb()],provenance=True)
#p=PostProcessor(".",files,selection.replace('\n',' '),"keep_and_drop.txt",[jecUncertAll_cppOut(),jetmetUncertaintiesAll(),btagSFProducer("cmva"),vhbb()],provenance=True)
#p.run()


if isMC:
    if era == "2016":
        p=PostProcessor(".",files,selection.replace('\n',' '),"keep_and_drop.txt",modules=[puWeight(),jetmetUncertainties2016All(),jetmetUncertainties2016AK8PuppiAllNoGroom(),muonScaleRes2016(),mhtVHbb(),btagSFProducer("2016","cmva"),vhbb2016()],provenance=True)
    elif era == "2017":
        p=PostProcessor(".",files,selection.replace('\n',' '),"keep_and_drop.txt",[puAutoWeight(),jetmetUncertainties2017All(),jetmetUncertainties2017AK8PuppiAll(),muonScaleRes2017(),mhtVHbb(),btagSFProducer("2017","deepcsv"),vhbb2017()],provenance=True)
        #p=PostProcessor(".",files,selection.replace('\n',' '),"keep_and_drop.txt",[puAutoWeight(),jetmetUncertainties2017All(),jetmetUncertainties2016AK8PuppiAll(),muonScaleRes2017(),mhtVHbb(),btagSFProducer("2017","deepcsv"),vhbb2017()],provenance=True)
else:
    if era == "2016":
        p=PostProcessor(".",files,selection.replace('\n',' '),"keep_and_drop.txt",modules=[mhtVHbb(),vhbb2016_data()],provenance=True)
    elif era == "2017":
        if dataRun == "B":
            p=PostProcessor(".",files,selection.replace('\n',' '),"keep_and_drop.txt",modules=[jetRecalib2017B(),mhtVHbb(),vhbb2017_data()],provenance=True)
        if dataRun == "C":
            p=PostProcessor(".",files,selection.replace('\n',' '),"keep_and_drop.txt",modules=[jetRecalib2017C(),mhtVHbb(),vhbb2017_data()],provenance=True)
        if dataRun == "D":
            p=PostProcessor(".",files,selection.replace('\n',' '),"keep_and_drop.txt",modules=[jetRecalib2017D(),mhtVHbb(),vhbb2017_data()],provenance=True)
        if dataRun == "E":
            p=PostProcessor(".",files,selection.replace('\n',' '),"keep_and_drop.txt",modules=[jetRecalib2017E(),mhtVHbb(),vhbb2017_data()],provenance=True)
        if dataRun == "F":
            p=PostProcessor(".",files,selection.replace('\n',' '),"keep_and_drop.txt",modules=[jetRecalib2017F(),mhtVHbb(),vhbb2017_data()],provenance=True)
p.run()

print "DONE"
os.system("ls -lR")
