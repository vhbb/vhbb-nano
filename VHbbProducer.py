import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import sys

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import * #deltaR, matching etc..

class VHbbProducer(Module):
    def __init__(self, isMC, era):
        self.era = era
        self.isMC = isMC
        pass
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Vtype",  "I");
        self.out.branch("V_pt",  "F");
        self.out.branch("V_eta",  "F");
        self.out.branch("V_phi",  "F");
        self.out.branch("V_mass",  "F");
        self.out.branch("V_mt",  "F");
        self.out.branch("Jet_lepFilter",  "O", 1, "nJet");
        self.out.branch("hJidx",  "I", 2);
        self.out.branch("hJidxCMVA",  "I", 2);
        self.out.branch("H_pt",  "F");
        self.out.branch("H_eta",  "F");
        self.out.branch("H_phi",  "F");
        self.out.branch("H_mass",  "F");
        self.out.branch("HCMVA_pt",  "F");
        self.out.branch("HCMVA_eta",  "F");
        self.out.branch("HCMVA_phi",  "F");
        self.out.branch("HCMVA_mass",  "F");
        self.out.branch("SA_Ht",  "F");
        self.out.branch("SA5",  "F");
        self.out.branch("Jet_Pt", "F", 1, "nJet");
        self.out.branch("MET_Pt","F");
        self.out.branch("MET_Phi","F");
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def matchSoftActivity(self,jets,saJets) :
	matched=set()
	for saj in saJets:
	    for j in jets :
	        if deltaR(saj,j) < 0.4 :
                    matched.add(saj)
	return matched
			
    def pt(self, jet, isMC):
        ## the MC has JER smearing applied which has output branch Jet_pt_nom which should be compared 
        ## with data branch Jet_pt. This essentially aliases the two branches to one common jet pt variable.
        if isMC:
            return jet.pt_nom
        else:
            return jet.pt		    
    
    def met(self, met, isMC):
        ## the MC has JER smearing applied which has output branch met_[pt/phi]_nom which should be compared 
        ## with data branch MET_[pt/phi]. This essentially aliases the two branches to one common variable.
        if isMC:
            return (met.pt_nom,met.phi_nom)
        else:
            return (met.pt,met.phi)
	  
    def btag(self, jet):
        ## make the btagging discriminator used a function of era. CMVA for 2016 and deepCSV for 2017
        if (self.era == "2016"):
            return jet.btagCMVA
	elif (self.era == "2017"):
            return jet.btagDeepB
        else:
            print "tried to call btag() for unknown era: %s, quitting..." % self.era
            sys.exit(1) 

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        jets = list(Collection(event, "Jet"))
        met = Object(event, "MET")
        sa = Collection(event, "SoftActivityJet")

        metPt,metPhi = self.met(met,self.isMC)
        self.out.fillBranch("MET_Pt",metPt)
        self.out.fillBranch("MET_Phi",metPhi) 
      
        Vtype = -1

        wElectrons = [x for x in electrons if x.mvaSpring16GP_WP80 and x.pt > 25 and x.pfRelIso03_all < 0.12]      
        wMuons = [x for x in muons if x.pt > 25 and x.tightId >= 1 and x.pfRelIso04_all < 0.15 and x.dxy < 0.05 and x.dz < 0.2]
        zElectrons = [x for x in electrons if x.pt > 20 and x.mvaSpring16GP_WP90 and x.pfRelIso03_all < 0.15]
        zMuons = [x for x in muons if x.pt > 20 and x.pfRelIso04_all < 0.25 and x.dxy < 0.05 and x.dz < 0.2] # muons already preselected with looseId requirement

        zMuons.sort(key=lambda x:x.pt,reverse=True)
        zElectrons.sort(key=lambda x:x.pt,reverse=True)

        vLeptons = [] # decay products of V
        if len(zMuons) >= 2:
            if zMuons[0].pt > 20:
                for i in xrange(1,len(zMuons)):
                    if zMuons[0].charge * zMuons[i].charge < 0:
                        Vtype = 0
                        vLeptons = [zMuons[0],zMuons[i]]
                        break
        elif len(zElectrons) >= 2:
            if zElectrons[0].pt > 20:
                for i in xrange(1,len(zElectrons)):
                    if zElectrons[0].charge * zElectrons[i].charge < 0:
                        Vtype = 1
                        vLeptons = [zElectrons[0],zElectrons[i]]
                        break
        elif len(wElectrons) + len(wMuons) == 1:
            if len(wMuons) == 1:
                Vtype = 2
                vLeptons = [wMuons[0]]
            if len(wElectrons) == 1:
                Vtype=3
                vLeptons = [wElectrons[0]]
        elif len(zElectrons) + len(zMuons) > 0:
            Vtype = 5
        else:
            Vtype = 4
            if metPt < 150:
                Vtype = -1
        self.out.fillBranch("Vtype",Vtype)

        ## add branches for some basic V kinematics
        V = ROOT.TLorentzVector()
        for vLepton in vLeptons:
            vLepton_4vec = ROOT.TLorentzVector()
            vLepton_4vec.SetPtEtaPhiM(vLepton.pt,vLepton.eta,vLepton.phi,vLepton.mass)
            V = V + vLepton_4vec
        if Vtype >=2 and Vtype<=4:
            met_4vec = ROOT.TLorentzVector()
            met_4vec.SetPtEtaPhiM(met.pt,0.,met.phi,0.) # only use met vector to derive transverse quantities
            V = V + met_4vec
        self.out.fillBranch("V_pt",V.Pt())
        self.out.fillBranch("V_eta",V.Eta())
        self.out.fillBranch("V_phi",V.Phi())
        self.out.fillBranch("V_mass",V.M())
        self.out.fillBranch("V_mt",V.Mt())

        ## filter jets that overlap with any of the selected leptons
        allLeptons = zElectrons[:]
        allLeptons.extend(zMuons)
        allLeptons.extend(wElectrons)
        allLeptons.extend(wMuons)
        jetFilterFlags = [True]*len(jets)
        #for jet in jets:
        #    jet.jetFilter = True
        for lepton in allLeptons:
            jetInd = lepton.jetIdx
            if jetInd >= 0:
                jetFilterFlags[jetInd] = False
                #jets[jetInd].jetFilter = False
        self.out.fillBranch("Jet_lepFilter",jetFilterFlags)

        ## alias JER-smeared MC jet pT and data jet pT to the same
        ## branch name
        jetPts = [-99.]*len(jets)
        for i in xrange(len(jets)):
            jetPts[i] = self.pt(jets[i],self.isMC)

        self.out.fillBranch("Jet_Pt",jetPts)

        ## Add explicit indices for selected H(bb) candidate jets
        jetsForHiggs = [x for x in jets if x.lepFilter and x.puId>0 and x.jetId>0 and x.Pt>20 and abs(x.eta)<2.5]
        if (len(jetsForHiggs) < 2): return False
        hJets = sorted(jetsForHiggs, key = lambda jet : self.btag(jet), reverse=True)[0:2]
        hJidx = [jets.index(x) for x in hJets]
        self.out.fillBranch("hJidx",hJidx)
        
        ## Save a few basic reco. H kinematics
        hj1 = ROOT.TLorentzVector()
        hj2 = ROOT.TLorentzVector()
        hj1.SetPtEtaPhiM(jets[hJidx[0]].Pt,jets[hJidx[0]].eta,jets[hJidx[0]].phi,jets[hJidx[0]].mass)
        hj2.SetPtEtaPhiM(jets[hJidx[1]].Pt,jets[hJidx[1]].eta,jets[hJidx[1]].phi,jets[hJidx[1]].mass)
        hbb = hj1 + hj2
        self.out.fillBranch("H_pt",hbb.Pt())
        self.out.fillBranch("H_phi",hbb.Phi())
        self.out.fillBranch("H_eta",hbb.Eta())
        self.out.fillBranch("H_mass",hbb.M())
        
        ## calculate separately selected indices using CMVA, although keep in mind this is already the
        ## default for 2016
        hJetsCMVA = sorted(jetsForHiggs, key = lambda jet : jet.btagCMVA, reverse=True)[0:2]
        hJidxCMVA = [jets.index(x) for x in hJetsCMVA]
        self.out.fillBranch("hJidxCMVA",hJidxCMVA)
        
        ## Save a few basic reco. H kinematics (from CMVA)
        hj1cmva = ROOT.TLorentzVector()
        hj2cmva = ROOT.TLorentzVector()
        hj1cmva.SetPtEtaPhiM(jets[hJidxCMVA[0]].Pt,jets[hJidxCMVA[0]].eta,jets[hJidxCMVA[0]].phi,jets[hJidxCMVA[0]].mass)
        hj2cmva.SetPtEtaPhiM(jets[hJidxCMVA[1]].Pt,jets[hJidxCMVA[1]].eta,jets[hJidxCMVA[1]].phi,jets[hJidxCMVA[1]].mass)
        hbbcmva = hj1cmva + hj2cmva
        self.out.fillBranch("HCMVA_pt",hbbcmva.Pt())
        self.out.fillBranch("HCMVA_phi",hbbcmva.Phi())
        self.out.fillBranch("HCMVA_eta",hbbcmva.Eta())
        self.out.fillBranch("HCMVA_mass",hbbcmva.M())


	## Compute soft activity vetoing Higgs jets
	#find signal footprint
	matchedSAJets=self.matchSoftActivity(hJets,sa)
	# update SA variables 


	softActivityJetHT=event.SoftActivityJetHT2-sum([x.pt for x in matchedSAJets])
        self.out.fillBranch("SA_Ht",softActivityJetHT)

	matchedSAJetsPt5=[x for x in matchedSAJets if x.pt>5]
        softActivityJetNjets5=event.SoftActivityJetNjets5-len(matchedSAJetsPt5)
        self.out.fillBranch("SA5",softActivityJetNjets5)

        return True
                

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

vhbb2016 = lambda : VHbbProducer(True,"2016") 
vhbb2016_data = lambda : VHbbProducer(False,"2016") 
vhbb2017 = lambda : VHbbProducer(True,"2017") 
vhbb2017_data = lambda : VHbbProducer(False,"2017") 
