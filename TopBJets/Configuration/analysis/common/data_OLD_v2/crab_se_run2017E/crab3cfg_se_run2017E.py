from WMCore.Configuration import Configuration

config = Configuration()

config.section_('General')
config.General.requestName     = 'se_run2017E'
config.General.transferOutputs = True
config.General.transferLogs    = False

config.section_('JobType')
config.JobType.pluginName  = 'Analysis'
config.JobType.allowUndistributedCMSSW = True 
config.JobType.psetName = '/depot/cms/top/chawla19/Topbjets/CMSSW_10_6_32/src/TopBJets/Configuration/analysis/common/test/ConfFile_cfg_2017_Data.py'
config.JobType.maxMemoryMB = 2500
config.JobType.pyCfgParams = ['era=2017', 'samplename=data', 'mode=se', 'datasetName=SingleElectron', 'outputFile=se_run2017E.root']
config.JobType.allowUndistributedCMSSW = True

config.section_('User')

config.section_('Site')
config.Site.storageSite = 'T2_US_Purdue'

config.section_('Data')
config.Data.publication    = False
config.Data.ignoreLocality = False
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = '/depot/cms/top/chawla19/Topbjets/CMSSW_10_6_32/src/TopBJets/Configuration/python/Data/Run2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
config.Data.inputDataset = '/SingleElectron/Run2017E-UL2017_MiniAODv2-v1/MINIAOD'
config.Data.outLFNDirBase = '/store/user/rchawla/crab/20230731/se_run2017E'
config.Data.unitsPerJob = 100
