import sys
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'VHbbPostNano2017_V2'
config.General.workArea = '/afs/cern.ch/work/s/scoopers/private/crabspace/crab_projects/2017/V2/'
#config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'crab_script.sh'
config.JobType.scriptArgs = ['isMC=1','era=2017','dataRun=""']
config.JobType.inputFiles = ['../keep_and_drop.txt','../postproc.py','../../../../../../scripts/haddnano.py'] #hadd nano will not be needed once nano tools are in cmssw
config.JobType.sendPythonFolder	 = True

config.Data.inputDataset = '/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
config.Data.inputDBS = 'phys03'
#config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
#config.Data.outLFNDirBase = '/store/user/%s/VHbbPostNano2017_V1/' % (getUsernameFromSiteDB())
config.Data.outLFNDirBase = '/store/group/phys_higgs/hbb/ntuples/VHbbPostNano/2017/V2/'
config.Data.publication = True
config.Data.outputDatasetTag = 'RunIISummer17MiniAOD-92X-VHbbPostNano2017_V2'
#config.Data.outputDatasetTag = 'RunIISummer16MiniAODv2-PUMoriond17-80X-VHbbPostNano2016_V1'
config.Data.allowNonValidInputDataset = True
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt'
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
config.Site.storageSite = 'T2_CH_CERN'

#sites=['T2_IT_Legnaro','T2_IT_Bari','T2_IT_Pisa','T2_CH_CERN']
sites=['T2_CH_CERN']

if __name__ == '__main__':
    f=open(sys.argv[1]) 
    content = f.readlines()
    content = [x.strip() for x in content] 
    from CRABAPI.RawCommand import crabCommand
    n=79
    for dataset in content :
	#site=sites[n%4]
	#config.Site.storageSite=site
	#if site=='T2_CH_CERN' :
	#	config.Data.outLFNDirBase=  '/store/group/cmst3/group/nanoAOD/NanoTestProd006'
	#else :
   #		config.Data.outLFNDirBase = '/store/user/%s/NanoTestProd006/' % (getUsernameFromSiteDB())

        config.Data.inputDataset = dataset
	config.Data.unitsPerJob = 2000000
	n+=1
	nnn="%s"%n
        config.General.requestName = "VHbbPostNano2017_V2_"+dataset.split('/')[1][:30]+dataset.split('/')[2][:30]+nnn
        config.Data.outputDatasetTag = dataset.split('/')[2][:30]+nnn
        crabCommand('submit', config = config)
