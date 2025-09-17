import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("Demo")

# Load the standard set of configuration modules
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
src = cms.InputTag("prunedGenParticles"),
)

process.printList = cms.EDAnalyzer("ParticleListDrawer",
    src = cms.InputTag("packedGenParticles"),
    maxEventsToPrint = cms.untracked.int32(-1),
    printVertex = cms.untracked.bool(True)
)

# Global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v35')

# Message logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

# Set input to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:0AC879DF-FE75-1E46-BFE8-10CC99EB5628.root'),
    #dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
    #inputCommands=cms.untracked.vstring(
    #        'keep *', 
    #        'drop recoTrackExtrasedmAssociation_muonReducedTrackExtras__RECO')
)

#define the default IDs to produce in VID
_defaultEleIDModules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff']

from EgammaUser.EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runVID=True,
                       applyEnergyCorrections=False,
                       isMiniAOD=True,
                       era='2017-UL',
                       eleIDModules=_defaultEleIDModules,
                       runEnergyCorrections=False) 

process.demo = cms.EDAnalyzer('TopAnalyzer',
    isMC = cms.untracked.bool(False),
    trigger = cms.InputTag("TriggerResults::HLT"),
    egmElectronVetoID = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-veto"),
    egmElectronLooseID = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-loose"),
    egmElectronMediumID = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-medium"),
    egmElectronTightID = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-tight"),
    prunedgenparticles = cms.InputTag("prunedGenParticles"),
    packedgenparticles = cms.InputTag("packedGenParticles"),
    offlinePrimaryVertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
    pfCandidates = cms.InputTag("packedPFCandidates"),
    jets = cms.InputTag("slimmedJets"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    met = cms.InputTag("slimmedMETs")
)

process.primaryVertexFilter  = cms.EDFilter("VertexSelector",
    src = cms.InputTag('offlineSlimmedPrimaryVertices'),
    cut = cms.string('!isFake && ndof > 4.0 && position.Rho < 2.0 && abs(z) < 24'),
    filter = cms.bool(True)
)

# Tfile service
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('data.root')
)

process.p = cms.Path(process.egammaPostRecoSeq+process.demo)
