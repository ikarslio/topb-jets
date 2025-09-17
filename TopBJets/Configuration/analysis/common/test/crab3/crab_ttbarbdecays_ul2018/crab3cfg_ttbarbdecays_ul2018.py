from WMCore.Configuration import Configuration

config = Configuration()

config.section_('General')
config.General.requestName     = 'ttbarbdecays_ul2018'
config.General.transferOutputs = True
config.General.transferLogs    = False

config.section_('JobType')
config.JobType.pluginName  = 'Analysis'
config.JobType.allowUndistributedCMSSW = True 
config.JobType.psetName = '/depot/cms/top/chawla19/Topbjets/CMSSW_10_6_32/src/TopBJets/Configuration/analysis/common/test/ConfFile_cfg_2018_MC.py'
config.JobType.maxMemoryMB = 800
config.JobType.pyCfgParams = ['era=2018', 'samplename=ttbarbdecays', 'datasetName=TTTo2L2Nu_CustomBdecays_TuneCP5_13TeV-powheg-pythia8', 'outputFile=ttbarbdecays_ul2018.root']
config.JobType.allowUndistributedCMSSW = True

config.section_('User')

config.section_('Site')
config.Site.storageSite = 'T2_US_Purdue'

config.section_('Data')
config.Data.publication    = False
config.Data.ignoreLocality = False
config.Data.splitting = 'EventAwareLumiBased'
config.Data.inputDataset = '/TTTo2L2Nu_CustomBdecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'
config.Data.outLFNDirBase = '/store/user/rchawla/crab/20230731/ttbarbdecays_ul2018'
config.Data.unitsPerJob = 100
