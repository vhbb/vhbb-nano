#!/usr/bin/env python
import os, sys
import ROOT
import argparse
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

print "args are: ",sys.argv

isMC = True
era = "2016"
dataRun = "X"
#if len(sys.argv) > 2:
#    if float(sys.argv[2]) < 0.:
#        isMC = False
#if len(sys.argv) > 3:
#    era = sys.argv[3]
#if era!="2016" and era!="2017":
#    print "Run era must be 2016 or 2017, exiting.."
#    sys.exit(1)
#btagger = "deepcsv"
#if era == "2016":
#    btagger = "cmva"
dataRun = ""
#if len(sys.argv) > 4:
#    dataRun = sys.argv[4]

parser = argparse.ArgumentParser("")
parser.add_argument('-jobNum', '--jobNum', type=int, default=1, help="")
parser.add_argument('-isMC', '--isMC', type=int, default=1, help="")
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
files=["root://cms-xrd-global.cern.ch://store/user/arizzi/NanoCrabProdXmas/WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_pythia8/NanoCrabXmasRunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asym___heIV_v6-v1__215/171221_111154/0000/nano_1.root"]
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
# &&  Sum$((abs(Jet_eta)<2.5 && Jet_pt > 20 && Jet_jetId)) >= 2 
#'''
selection='''(Sum$(Electron_pt > 20 && Electron_mvaSpring16GP_WP90) >= 2  ||
 Sum$(Muon_pt > 20) >= 2 ||
 Sum$(Electron_pt > 20 && Electron_mvaSpring16GP_WP80) >= 1   ||
 Sum$(Muon_pt > 20 && Muon_tightId) >= 1 ||
 (Sum$(Muon_pt > 20) == 0 && Sum$(Electron_pt > 20 && Electron_mvaSpring16GP_WP90) == 0 && MET_pt > 150 ) ) 
 &&   ((Sum$((abs(Jet_eta)<2.5 && Jet_pt > 20 && Jet_jetId)) >= 2)||(Sum$((abs(FatJet_eta)<2.5 && FatJet_pt > 200 && FatJet_jetId)) >= 1)) 
'''

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

#this takes care of converting the input files from CRAB
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles,runsAndLumis

if isMC:
    if era == "2016":
        p=PostProcessor(".",inputFiles(),selection.replace('\n',' '),"keep_and_drop.txt",modules=[puWeight(),jetmetUncertainties2016All(),jetmetUncertainties2016AK8PuppiAll(),muonScaleRes2016(),mhtVHbb(),btagSFProducer("2016","cmva"),vhbb2016()],provenance=True,fwkJobReport=True,jsonInput=runsAndLumis())
    elif era == "2017":
        p=PostProcessor(".",inputFiles(),selection.replace('\n',' '),"keep_and_drop.txt",[puAutoWeight(),jetmetUncertainties2017All(),jetmetUncertainties2016AK8PuppiAll(),muonScaleRes2017(),mhtVHbb(),btagSFProducer("2017","deepcsv"),vhbb2017()],provenance=True,fwkJobReport=True,jsonInput=runsAndLumis())
else:
    if era == "2016":
        p=PostProcessor(".",inputFiles(),selection.replace('\n',' '),"keep_and_drop.txt",modules=[mhtVHbb(),vhbb2016_data()],provenance=True,fwkJobReport=True,jsonInput=runsAndLumis())
    elif era == "2017":
        if dataRun == "B":
            p=PostProcessor(".",inputFiles(),selection.replace('\n',' '),"keep_and_drop.txt",modules=[jetRecalib2017B(),mhtVHbb(),vhbb2017_data()],provenance=True,fwkJobReport=True,jsonInput=runsAndLumis())
        if dataRun == "C":
            p=PostProcessor(".",inputFiles(),selection.replace('\n',' '),"keep_and_drop.txt",modules=[jetRecalib2017C(),mhtVHbb(),vhbb2017_data()],provenance=True,fwkJobReport=True,jsonInput=runsAndLumis())
        if dataRun == "D":
            p=PostProcessor(".",inputFiles(),selection.replace('\n',' '),"keep_and_drop.txt",modules=[jetRecalib2017D(),mhtVHbb(),vhbb2017_data()],provenance=True,fwkJobReport=True,jsonInput=runsAndLumis())
        if dataRun == "E":
            p=PostProcessor(".",inputFiles(),selection.replace('\n',' '),"keep_and_drop.txt",modules=[jetRecalib2017E(),mhtVHbb(),vhbb2017_data()],provenance=True,fwkJobReport=True,jsonInput=runsAndLumis())
        if dataRun == "F":
            p=PostProcessor(".",inputFiles(),selection.replace('\n',' '),"keep_and_drop.txt",modules=[jetRecalib2017F(),mhtVHbb(),vhbb2017_data()],provenance=True,fwkJobReport=True,jsonInput=runsAndLumis())
p.run()

print "DONE"
os.system("ls -lR")
