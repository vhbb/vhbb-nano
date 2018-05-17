import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import sys

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import * #deltaR, matching etc..

class VHbbProducer(Module):
    def __init__(self, isMC, era, useCMVA=False):
        self.era = era
        self.isMC = isMC
        self.useCMVA = useCMVA
        self.genJets = []
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
        self.out.branch("HFSR_pt",  "F");
        self.out.branch("HFSR_eta",  "F");
        self.out.branch("HFSR_phi",  "F");
        self.out.branch("HFSR_mass",  "F");
        self.out.branch("SA_Ht",  "F");
        self.out.branch("SA5",  "F");
        self.out.branch("Jet_Pt", "F", 1, "nJet");
        self.out.branch("Jet_PtReg", "F", 1, "nJet");
        self.out.branch("MET_Pt","F");
        self.out.branch("MET_Phi","F");

        ## for the boosted analysis
        self.out.branch("Pt_fjidx",  "I");        
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

    def matchSoftActivity(self,jets,saJets,dR=0.4) :
	matched=set()
	for saj in saJets:
	    for j in jets :
	        if deltaR(saj,j) < dR :
                    matched.add(saj)
	return matched
    
    def matchSoftActivityFSR(self,jet1,jet2,saJets,dR=0.4) :
	matched=set()
        drjj = deltaR(jet1,jet2)
        sumDeltaRMin = drjj + 2*dR
	for saj in saJets:
            dr1 = deltaR(saj,jet1)
            dr2 = deltaR(saj,jet2)
            if ((dr1+dr2) < sumDeltaRMin):
                matched.add(saj)
	return matched
			
    def pt(self, jet, isMC, noReg=False):
        ## the MC has JER smearing applied which has output branch Jet_pt_nom which should be compared 
        ## with data branch Jet_pt. This essentially aliases the two branches to one common jet pt variable.
        if noReg:
            if isMC:
                return jet.pt_nom
            else:
                return jet.pt		 
        else:
            if isMC:
                # until we have final post-regression smearing factors we assume a flat 10%
                smearedPt = jet.pt
                if jet.genJetIdx >=0 and  jet.genJetIdx < len(self.genJets) :
                    genJet=self.genJets[jet.genJetIdx]
                    dPt = jet.pt - genJet.pt
                    smearedPt=genJet.pt+1.1*dPt
                return jet.bReg*smearedPt
            else:
                return jet.bReg*jet.pt    
    
    def met(self, met, isMC):
        ## the MC has JER smearing applied which has output branch met_[pt/phi]_nom which should be compared 
        ## with data branch MET_[pt/phi]. This essentially aliases the two branches to one common variable.
        if isMC:
            return (met.pt_nom,met.phi_nom)
        else:
            return (met.pt,met.phi)
	 
    def msoftdrop(self, jet, isMC):
        ## the MC has JER smearing applied which has output branch Jet_pt_smeared which should be compared 
        ## with data branch Jet_pt. This essentially aliases the two branches to one common jet pt variable.
        if isMC and self.era!="2016":
            return jet.msoftdrop_nom
        else:
            return jet.msoftdrop
 
    def btag(self, jet):
        if (self.useCMVA):
            return jet.btagCMVA
	else:
            return jet.btagDeepB

    def elid(self, el, wp):
        if (self.era == "2016" and wp == "80"):
            return el.mvaSpring16GP_WP80
        elif (self.era == "2016" and wp == "90"):
            return el.mvaSpring16GP_WP90
        elif (self.era == "2017" and wp == "80"):
            return el.mvaFall17Iso_WP80
        elif (self.era == "2017" and wp == "90"):
            return el.mvaFall17Iso_WP90
    def statusFlags_dict(self, bit):

        Dict={0 : "isPrompt", 1 : "isDecayedLeptonHadron", 2 : "isTauDecayProduct", 3 : "isPromptTauDecayProduct", 4 : "isDirectTauDecayProduct", 5 : "isDirectPromptTauDecayProduct", 6 : "isDirectHadronDecayProduct", 7 : "isHardProcess", 8 : "fromHardProcess", 9 : "isHardProcessTauDecayProduct", 10 : "isDirectHardProcessTauDecayProduct", 11 : "fromHardProcessBeforeFSR", 12 : "isFirstCopy", 13 : "isLastCopy", 14 : "isLastCopyBeforeFSR" }
        return Dict[bit] 

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        jets = list(Collection(event, "Jet"))
        met = Object(event, "MET")
        sa = Collection(event, "SoftActivityJet")
        fatjets = list(Collection(event, "FatJet"))
        subjets = Collection(event, "SubJet")
        if self.isMC:
            genParticles = Collection(event, "GenPart")
            self.genJets = Collection(event, "GenJet")

        metPt,metPhi = self.met(met,self.isMC)
        self.out.fillBranch("MET_Pt",metPt)
        self.out.fillBranch("MET_Phi",metPhi) 
      
        Vtype = -1

        wElectrons = [x for x in electrons if self.elid(x,"80") and x.pt > 25 and x.pfRelIso03_all < 0.12]      
        wMuons = [x for x in muons if x.pt > 25 and x.tightId >= 1 and x.pfRelIso04_all < 0.15 and x.dxy < 0.05 and x.dz < 0.2]
        zElectrons = [x for x in electrons if x.pt > 20 and self.elid(x,"90") and x.pfRelIso03_all < 0.15]
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
            met_4vec.SetPtEtaPhiM(metPt,0.,metPhi,0.) # only use met vector to derive transverse quantities
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
        fatjetFilterFlags = [True]*len(fatjets)
        #for jet in jets:
        #    jet.jetFilter = True
        for lepton in allLeptons:
            jetInd = lepton.jetIdx
            if jetInd >= 0:
                jetFilterFlags[jetInd] = False
                #jets[jetInd].jetFilter = False
        self.out.fillBranch("Jet_lepFilter",jetFilterFlags)
 
        for fatjet in fatjets:
            fatjet.jetFilter = True
            for lepton in allLeptons:
               if deltaR(fatjet,lepton) < 0.8:
                  jetFilterFlags[fatjets.index(fatjet)] = False

        self.out.fillBranch("FatJet_lepFilter",jetFilterFlags)

        ## alias JER-smeared MC jet pT and data jet pT to the same
        ## branch name
        jetPts = [-99.]*len(jets)
        jetPtRegs = [-99.]*len(jets)
        for i in xrange(len(jets)):
            jetPts[i] = self.pt(jets[i],self.isMC, True)
            jetPtRegs[i] = self.pt(jets[i],self.isMC)

        self.out.fillBranch("Jet_Pt",jetPts)
        self.out.fillBranch("Jet_PtReg",jetPtRegs)

        fatjetPts = [-99.]*len(fatjets)
        for i in xrange(len(fatjets)):
            fatjetPts[i] = self.pt(fatjets[i],self.isMC,True)

        fatjetMSD = [-99.]*len(fatjets)
        for i in xrange(len(fatjets)):
            fatjetMSD[i] = self.msoftdrop(fatjets[i],self.isMC)

        self.out.fillBranch("FatJet_Pt",fatjetPts)
        self.out.fillBranch("FatJet_Msoftdrop",fatjetMSD)

        ## Add explicit indices for selected H(bb) candidate jets
        jetsForHiggs = [x for x in jets if x.lepFilter and x.puId>0 and x.jetId>0 and self.pt(x,self.isMC)>20 and abs(x.eta)<2.5]
        if (len(jetsForHiggs) >= 2): 
            hJets = sorted(jetsForHiggs, key = lambda jet : self.btag(jet), reverse=True)[0:2]
            hJidx = [jets.index(x) for x in hJets]
            self.out.fillBranch("hJidx",hJidx)

            ## Save a few basic reco. H kinematics
            hj1 = ROOT.TLorentzVector()
            hj2 = ROOT.TLorentzVector()
            hj1.SetPtEtaPhiM(self.pt(jets[hJidx[0]],self.isMC),jets[hJidx[0]].eta,jets[hJidx[0]].phi,jets[hJidx[0]].mass)
            hj2.SetPtEtaPhiM(self.pt(jets[hJidx[1]],self.isMC),jets[hJidx[1]].eta,jets[hJidx[1]].phi,jets[hJidx[1]].mass)
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
            hj1cmva.SetPtEtaPhiM(self.pt(jets[hJidxCMVA[0]],self.isMC),jets[hJidxCMVA[0]].eta,jets[hJidxCMVA[0]].phi,jets[hJidxCMVA[0]].mass)
            hj2cmva.SetPtEtaPhiM(self.pt(jets[hJidxCMVA[1]],self.isMC),jets[hJidxCMVA[1]].eta,jets[hJidxCMVA[1]].phi,jets[hJidxCMVA[1]].mass)
            hbbcmva = hj1cmva + hj2cmva
            self.out.fillBranch("HCMVA_pt",hbbcmva.Pt())
            self.out.fillBranch("HCMVA_phi",hbbcmva.Phi())
            self.out.fillBranch("HCMVA_eta",hbbcmva.Eta())
            self.out.fillBranch("HCMVA_mass",hbbcmva.M())

            ## try to recover FSR
            jetsFromFSR = []
            for ijet in xrange(len(jets)):
                if ijet == hJidx[0] or ijet == hJidx[1]: continue
                jet = jets[ijet]
                if self.pt(jet,self.isMC,noReg=True)>20 and abs(jet.eta)<3.0 and jet.puId>0 and jet.jetId>0 and jet.lepFilter:
                   if min(deltaR(jet,jets[hJidx[0]]),deltaR(jet,jets[hJidx[1]])) < 0.8:
                       jetsFromFSR.append(jet)
            HFSR = hbb
            for jet in jetsFromFSR:
                fsrJetToAdd = ROOT.TLorentzVector()
                fsrJetToAdd.SetPtEtaPhiM(jet.Pt,jet.eta,jet.phi,jet.mass)
                HFSR = HFSR + fsrJetToAdd
            self.out.fillBranch("HFSR_pt",HFSR.Pt())
            self.out.fillBranch("HFSR_phi",HFSR.Phi())
            self.out.fillBranch("HFSR_eta",HFSR.Eta())
            self.out.fillBranch("HFSR_mass",HFSR.M())
            


            ## Compute soft activity vetoing Higgs jets
            #find signal footprint
            matchedSAJets=self.matchSoftActivity(hJets,sa)
            #matchedSAJets=self.matchSoftActivityFSR(hJets[0],hJets[1],sa)
            # update SA variables 


            softActivityJetHT=event.SoftActivityJetHT2-sum([x.pt for x in matchedSAJets])
            self.out.fillBranch("SA_Ht",softActivityJetHT)

            matchedSAJetsPt5=[x for x in matchedSAJets if x.pt>5]
            softActivityJetNjets5=event.SoftActivityJetNjets5-len(matchedSAJetsPt5)
            self.out.fillBranch("SA5",softActivityJetNjets5)
            
            if (event.SoftActivityJetNjets5-len(matchedSAJetsPt5)) < 0:
                print "available soft activity jet collection"
                for sajet in sa:
                    print "pt =",sajet.pt,"eta = ",sajet.eta
                print "these jets were removed from the count"
                for sajet in matchedSAJetsPt5:
                    print "pt =",sajet.pt,"eta = ",sajet.eta
                print "in Nano the number of SA5 jets is event.SoftActivityJetNjets5 = ",event.SoftActivityJetNjets5
                print "after removing the matched jets SA5 is ",(event.SoftActivityJetNjets5-len(matchedSAJetsPt5))
        
        else:
            self.out.fillBranch("hJidx",[-1,-1])
            self.out.fillBranch("H_pt",-1)
            self.out.fillBranch("H_phi",-1)
            self.out.fillBranch("H_eta",-1)
            self.out.fillBranch("H_mass",-1)
            self.out.fillBranch("hJidxCMVA",[-1,-1])
            self.out.fillBranch("HCMVA_pt",-1)
            self.out.fillBranch("HCMVA_phi",-1)
            self.out.fillBranch("HCMVA_eta",-1)
            self.out.fillBranch("HCMVA_mass",-1)
            self.out.fillBranch("HFSR_pt",-1)
            self.out.fillBranch("HFSR_phi",-1)
            self.out.fillBranch("HFSR_eta",-1)
            self.out.fillBranch("HFSR_mass",-1)
            self.out.fillBranch("SA_Ht",-1)
            self.out.fillBranch("SA5",-1)
   
        ## indices for Hjets and soft activity (?)
        fatjetsForHiggs = [x for x in fatjets if x.lepFilter and x.jetId>0 and x.Pt>250 and x.Msoftdrop>40 and (x.eta)<2.5]
        if (len(fatjetsForHiggs) >= 1):

            jh = sorted(fatjetsForHiggs, key = lambda jet : jet.Pt, reverse=True)
            pt_idx = fatjets.index(jh[0])
            jh = sorted(fatjetsForHiggs, key = lambda jet : jet.Msoftdrop, reverse=True)
            msd_idx = fatjets.index(jh[0])
            jh = sorted(fatjetsForHiggs, key = lambda jet : jet.btagHbb, reverse=True)
            hbb_idx = fatjets.index(jh[0])
            self.out.fillBranch("Pt_fjidx",pt_idx)
            self.out.fillBranch("Msd_fjidx",msd_idx)
            self.out.fillBranch("Hbb_fjidx",hbb_idx)

            ## SA leading pt
            matchedSAJets=self.matchSoftActivity([fatjets[pt_idx]],sa,0.8)
            matchedSAJetsPt5=[x for x in matchedSAJets if x.pt>5]
            softActivityJetHT=event.SoftActivityJetHT2-sum([x.pt for x in matchedSAJets])
            self.out.fillBranch("SAptfj_HT",softActivityJetHT)
            softActivityJetNjets5=event.SoftActivityJetNjets5-len(matchedSAJetsPt5)
            self.out.fillBranch("SAptfj5",softActivityJetNjets5)

            ## SA leading mass
            matchedSAJets=self.matchSoftActivity([fatjets[msd_idx]],sa,0.8)
            matchedSAJetsPt5=[x for x in matchedSAJets if x.pt>5]
            softActivityJetHT=event.SoftActivityJetHT2-sum([x.pt for x in matchedSAJets])
            self.out.fillBranch("SAmfj_HT",softActivityJetHT)
            softActivityJetNjets5=event.SoftActivityJetNjets5-len(matchedSAJetsPt5)
            self.out.fillBranch("SAmfj5",softActivityJetNjets5)

            ## SA leading mass
            matchedSAJets=self.matchSoftActivity([fatjets[hbb_idx]],sa,0.8)
            matchedSAJetsPt5=[x for x in matchedSAJets if x.pt>5]
            softActivityJetHT=event.SoftActivityJetHT2-sum([x.pt for x in matchedSAJets])
            self.out.fillBranch("SAhbbfj_HT",softActivityJetHT)
            softActivityJetNjets5=event.SoftActivityJetNjets5-len(matchedSAJetsPt5)
            self.out.fillBranch("SAhbbfj5",softActivityJetNjets5)

        else:
            self.out.fillBranch("Pt_fjidx",-1)
            self.out.fillBranch("Msd_fjidx",-1)
            self.out.fillBranch("Hbb_fjidx",-1)
            self.out.fillBranch("SAptfj_HT",-1)
            self.out.fillBranch("SAptfj5",-1)
            self.out.fillBranch("SAmfj_HT",-1)
            self.out.fillBranch("SAmfj5",-1)
            self.out.fillBranch("SAhbbfj_HT",-1)
            self.out.fillBranch("SAhbbfj5",-1)

        FatJet_FlavourComposition=[1]*len(fatjets)
        FatJet_HiggsProducts=[False]*len(fatjets)
        FatJet_ZProducts=[False]*len(fatjets)
        FatJet_WProducts=[False]*len(fatjets)

        if self.isMC:
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
        else:
            # data event, fill with default values
            self.out.fillBranch("FatJet_FlavourComposition",FatJet_FlavourComposition)
            self.out.fillBranch("FatJet_HiggsProducts",FatJet_HiggsProducts)
            self.out.fillBranch("FatJet_ZProducts",FatJet_ZProducts)
            self.out.fillBranch("FatJet_WProducts",FatJet_WProducts)

        return True
                

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

vhbb2016 = lambda : VHbbProducer(True,"2016") 
vhbb2016_data = lambda : VHbbProducer(False,"2016") 
vhbb2017 = lambda : VHbbProducer(True,"2017") 
vhbb2017_data = lambda : VHbbProducer(False,"2017") 
