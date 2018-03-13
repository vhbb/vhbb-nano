import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import * #deltaR, matching etc..

class VHbbBoostedProducer(Module):
    def __init__(self, isMC):
        self.isMC = isMC
        pass
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        
        self.out.branch("Msd_fjidx",  "I");
        self.out.branch("Hbb_fjidx",  "I");
        
        self.out.branch("SAptfj_HT",  "F");
        self.out.branch("SAptfj5",  "F");
        self.out.branch("SAmfj_HT",  "F");
        self.out.branch("SAmfj5",  "F");
        self.out.branch("SAhbbfj_HT",  "F");
        self.out.branch("SAhbbfj5",  "F");
        
        self.out.branch("FatJet_lepFilter",  "O", 1, "nFatJet");
        self.out.branch("FatJet_Pt", "F", 1, "nFatJet");
        self.out.branch("FatJet_Msoftdrop", "F", 1, "nFatJet");
        
        self.out.branch("FatJet_FlavourComposition", "I", 1, "nFatJet"); #55 bb, #54 bc, #5 b, #4 c, #44 cc, #1 other
        
        self.out.branch("FatJet_HiggsProducts", "O", 1, "nFatJet");
        self.out.branch("FatJet_WProducts", "O", 1, "nFatJet");
        self.out.branch("FatJet_ZProducts", "O", 1, "nFatJet");


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def matchSoftActivity(self,jets,saJets) :
	matched=set()
	for saj in saJets:
	    for j in jets :
	        if deltaR(saj,j) < 0.8 :
                    matched.add(saj)
	return matched
    
    def pt(self, jet, isMC):
        ## the MC has JER smearing applied which has output branch Jet_pt_smeared which should be compared 
        ## with data branch Jet_pt. This essentially aliases the two branches to one common jet pt variable.
        if isMC:
            return jet.pt_nom
        else:
            return jet.pt
        
    def msoftdrop(self, jet, isMC):
        ## the MC has JER smearing applied which has output branch Jet_pt_smeared which should be compared 
        ## with data branch Jet_pt. This essentially aliases the two branches to one common jet pt variable.
        if isMC:
            return jet.msoftdrop_nom
        else:
            return jet.msoftdrop
        
    def statusFlags_dict(self, bit):
        
        Dict={0 : "isPrompt", 1 : "isDecayedLeptonHadron", 2 : "isTauDecayProduct", 3 : "isPromptTauDecayProduct", 4 : "isDirectTauDecayProduct", 5 : "isDirectPromptTauDecayProduct", 6 : "isDirectHadronDecayProduct", 7 : "isHardProcess", 8 : "fromHardProcess", 9 : "isHardProcessTauDecayProduct", 10 : "isDirectHardProcessTauDecayProduct", 11 : "fromHardProcessBeforeFSR", 12 : "isFirstCopy", 13 : "isLastCopy", 14 : "isLastCopyBeforeFSR" }
        return Dict[bit]
    
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        fatjets = list(Collection(event, "FatJet"))
        subjets = Collection(event, "SubJet")
        sa = Collection(event, "SoftActivityJet")
        genParticles = Collection(event, "GenPart")
                
        #do we need lepfilter here?
        
        wElectrons = [x for x in electrons if x.mvaSpring16GP_WP80 and x.pt > 25 and x.pfRelIso03_all < 0.12]      
        wMuons = [x for x in muons if x.pt > 25 and x.tightId >= 1 and x.pfRelIso04_all < 0.15 and x.dxy < 0.05 and x.dz < 0.2]
        zElectrons = [x for x in electrons if x.pt > 20 and x.mvaSpring16GP_WP90 and x.pfRelIso03_all < 0.15]
        zMuons = [x for x in muons if x.pt > 20 and x.pfRelIso04_all < 0.25 and x.dxy < 0.05 and x.dz < 0.2] # muons already preselected with looseId requirement

        ## filter jets that overlap with any of the selected leptons
        allLeptons = zElectrons[:]
        allLeptons.extend(zMuons)
        allLeptons.extend(wElectrons)
        allLeptons.extend(wMuons)
        jetFilterFlags = [True]*len(fatjets)
        for fatjet in fatjets:
            fatjet.jetFilter = True
            for lepton in allLeptons:
               if deltaR(fatjet,lepton) < 0.8:
                  jetFilterFlags[fatjets.index(fatjet)] = False

        self.out.fillBranch("FatJet_lepFilter",jetFilterFlags)    
        

        ## alias JER-smeared MC jet pT and data jet pT to the same
        ## branch name
        fatjetPts = [-99.]*len(fatjets)
        for i in xrange(len(fatjets)):
            fatjetPts[i] = self.pt(fatjets[i],self.isMC)
            
        fatjetMSD = [-99.]*len(fatjets)
        for i in xrange(len(fatjets)):
            fatjetMSD[i] = self.msoftdrop(fatjets[i],self.isMC)

        self.out.fillBranch("FatJet_Pt",fatjetPts)
        self.out.fillBranch("FatJet_Msoftdrop",fatjetMSD)
        
        FatJet_FlavourComposition=[1]*len(fatjets)
        FatJet_HiggsProducts=[False]*len(fatjets)
        FatJet_ZProducts=[False]*len(fatjets)
        FatJet_WProducts=[False]*len(fatjets)
        
        #do flavour composition here
        for fatjet in fatjets:
            
            numBHadrons=0
            numCHadrons=0
            H_decay=False
            W_decay=False
            Z_decay=False
            
            for genP in genParticles:
                if deltaR(fatjet,genP)<0.8:                    
                    id_at=max((abs(genP.pdgId)/1000) % 10,(abs(genP.pdgId)/100) % 10)
                    
                    if (id_at==4 or id_at==5):                        
                        lastinChain=True
                        if genP.genPartIdxMother<0: id_mother=-1
                        else: id_mother=max((abs(genParticles[genP.genPartIdxMother].pdgId)/1000) % 10,(abs(genParticles[genP.genPartIdxMother].pdgId)/100) % 10)
                        
                        if (id_mother==4 or id_mother==5):                             
                            lastinChain=False
                            
                        if lastinChain: 
                            if id_at==4: numCHadrons=numCHadrons+1
                            if id_at==5: numBHadrons=numBHadrons+1
                            
                            #subjet matching?
                            
                            #s1=fatjet.subJetIdx1
                            #s2=fatjet.subJetIdx2
   
                            #if s1>-1:
                                #if deltaR(subjets[s1],genP)<0.3: print "it's in subjet 1"
                            #if s2>-1:
                                #if deltaR(subjets[s2],genP)<0.3: print "it's in subjet 2"
                            
                                
                    if (abs(genP.pdgId)==24): W_decay=True
                    if (genP.pdgId==23): Z_decay=True
                    if (genP.pdgId==25): H_decay=True
                    
                    #status flags: maybe we want ot check a few bits in the status before assigning True?
                    
                    #if (genP.pdgId==25): 
                        #print "Higgs boson",genP.pdgId, genP.status, genP.statusFlags,genP.genPartIdxMother
                        #for bit in range(15):
                            #if ((genP.statusFlags>>bit)&1) : print self.statusFlags_dict(bit)
                    
            if numBHadrons>=2:
             FatJet_FlavourComposition[fatjets.index(fatjet)]=55
            elif numBHadrons==1 and numCHadrons>=1:
             FatJet_FlavourComposition[fatjets.index(fatjet)]=54
            elif numBHadrons==1 and numCHadrons==0:
             FatJet_FlavourComposition[fatjets.index(fatjet)]=5
            elif numBHadrons==0 and numCHadrons==2:
             FatJet_FlavourComposition[fatjets.index(fatjet)]=44
            elif numBHadrons==0 and numCHadrons==1:
             FatJet_FlavourComposition[fatjets.index(fatjet)]=4
            
            FatJet_HiggsProducts[fatjets.index(fatjet)]=H_decay
            FatJet_ZProducts[fatjets.index(fatjet)]=Z_decay
            FatJet_WProducts[fatjets.index(fatjet)]=W_decay
                 
        self.out.fillBranch("FatJet_FlavourComposition",FatJet_FlavourComposition)   
        self.out.fillBranch("FatJet_HiggsProducts",FatJet_HiggsProducts)
        self.out.fillBranch("FatJet_ZProducts",FatJet_ZProducts)
        self.out.fillBranch("FatJet_WProducts",FatJet_WProducts)
        
                    

        ## indices for Hjets and soft activity (?)
        fatjetsForHiggs = [x for x in fatjets if x.lepFilter and x.jetId>0 and x.Pt>250 and x.Msoftdrop>40 and (x.eta)<2.5]
        if (len(fatjetsForHiggs) > 1): 
            
            jh = sorted(fatjetsForHiggs, key = lambda jet : jet.Msoftdrop, reverse=True)
            msd_idx = fatjetsForHiggs.index(jh[0])
            jh = sorted(fatjetsForHiggs, key = lambda jet : jet.btagHbb, reverse=True) 
            hbb_idx = fatjetsForHiggs.index(jh[0])
            self.out.fillBranch("Msd_fjidx",msd_idx)
            self.out.fillBranch("Hbb_fjidx",hbb_idx)

            ## SA leading pt
            matchedSAJets=self.matchSoftActivity([fatjetsForHiggs[0]],sa)
            matchedSAJetsPt5=[x for x in matchedSAJets if x.pt>5]
            softActivityJetHT=event.SoftActivityJetHT2-sum([x.pt for x in matchedSAJets])
            self.out.fillBranch("SAptfj_HT",softActivityJetHT)
            softActivityJetNjets5=event.SoftActivityJetNjets5-len(matchedSAJetsPt5)
            self.out.fillBranch("SAptfj5",softActivityJetNjets5)

            ## SA leading mass
            matchedSAJets=self.matchSoftActivity([fatjetsForHiggs[msd_idx]],sa)
            matchedSAJetsPt5=[x for x in matchedSAJets if x.pt>5]
            softActivityJetHT=event.SoftActivityJetHT2-sum([x.pt for x in matchedSAJets])
            self.out.fillBranch("SAmfj_HT",softActivityJetHT)
            softActivityJetNjets5=event.SoftActivityJetNjets5-len(matchedSAJetsPt5)
            self.out.fillBranch("SAmfj5",softActivityJetNjets5)

            ## SA leading mass
            matchedSAJets=self.matchSoftActivity([fatjetsForHiggs[hbb_idx]],sa)
            matchedSAJetsPt5=[x for x in matchedSAJets if x.pt>5]
            softActivityJetHT=event.SoftActivityJetHT2-sum([x.pt for x in matchedSAJets])
            self.out.fillBranch("SAhbbfj_HT",softActivityJetHT)
            softActivityJetNjets5=event.SoftActivityJetNjets5-len(matchedSAJetsPt5)
            self.out.fillBranch("SAhbbfj5",softActivityJetNjets5)
        
        else:
            
            self.out.fillBranch("Msd_fjidx",-1)
            self.out.fillBranch("Hbb_fjidx",-1)
            self.out.fillBranch("SAptfj_HT",-1)
            self.out.fillBranch("SAptfj5",-1)
            self.out.fillBranch("SAmfj_HT",-1)
            self.out.fillBranch("SAmfj5",-1)
            self.out.fillBranch("SAhbbfj_HT",-1)
            self.out.fillBranch("SAhbbfj5",-1)

        
        


        return True
                

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

vhbb_boosted = lambda : VHbbBoostedProducer(True) 
vhbb_boosted_data = lambda : VHbbBoostedProducer(False) 
