import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

# process.dump=cms.EDAnalyzer('EventContentAnalyzer')
process.load("FWCore.MessageService.MessageLogger_cfi")



# loads from CJLST
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

import FWCore.ParameterSet.VarParsing as VarParsing
opts = VarParsing.VarParsing('python')
vpbool = VarParsing.VarParsing.varType.bool
opts.register('data'    , False  , mytype=vpbool)

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
fileFormat = 'miniAOD'

if fileFormat == 'AOD' :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)
switchOnVIDPhotonIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = [
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff',
                ]

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

switchOnVIDPhotonIdProducer(process, dataFormat)
# define which IDs we want to produce
my_id_modules = [
                 'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring15_25ns_nonTrig_V2_cff',
                 # 'RecoEgamma.PhotonIdentification.Identification.mvaTLEID_Fall15_V1_cff',
               ]

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)


# Setup ISO calculation
from RecoEgamma.EgammaIsolationAlgos.egmPhotonIsolationMiniAOD_cff import egmPhotonIsolationMiniAOD

process.egmPhotonIsolationMiniAOD = egmPhotonIsolationMiniAOD.clone()


#process.TFileService = cms.Service("TFileService", fileName = cms.string('/data/DATA/temp_pigard/eID/') )

# tle_ntuple = cms.EDAnalyzer('electronML',
tle_ntuple = cms.EDAnalyzer('Ntuplizer',
                                   inputFileFormat = cms.untracked.string(fileFormat),
                                   outputPath = cms.string(""), #process.TFileService.fileName,
                                   outputFile = cms.string('TLE'),
                                   beamSpot = cms.InputTag('offlineBeamSpot'),
                                   # input collection names AOD
                                   electronsAOD = cms.InputTag('gedGsfElectrons'),
                                   verticesAOD = cms.InputTag('offlinePrimaryVertices'),
                                   conversionsAOD = cms.InputTag('allConversions'),
                                   genParticlesAOD = cms.InputTag('genParticles'), 
                                   genEventInfoProductAOD = cms.InputTag('generator'),
                                   PFMETAOD = cms.InputTag('pfMet'),                                
                                   ebReducedRecHitCollectionAOD = cms.InputTag("reducedEcalRecHitsEB"),
                                   eeReducedRecHitCollectionAOD = cms.InputTag("reducedEcalRecHitsEE"),
                                   esReducedRecHitCollectionAOD = cms.InputTag("reducedEcalRecHitsES"),
 
                                   electronsMiniAOD = cms.InputTag('slimmedElectrons'),
                                   photonsMiniAOD = cms.InputTag('slimmedPhotons'),
                                   verticesMiniAOD = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                   conversionsMiniAOD = cms.InputTag('reducedEgamma:reducedConversions'),
                                   genParticlesMiniAOD = cms.InputTag('prunedGenParticles'), 
                                   genEventInfoProductMiniAOD = cms.InputTag('generator'),
                                   PFMETMiniAOD = cms.InputTag('slimmedMETs'),                                
                                   ebReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma:reducedEBRecHits"),
                                   eeReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma:reducedEERecHits"),
                                   esReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma:reducedESRecHits"),

                                   HLTTag          = cms.InputTag('TriggerResults','','HLT'),
                                   isMC = cms.bool(not opts.data),
                                   # do_signal = cms.bool(True),
                                   do_TLE = cms.bool(True),
                                   do_min_dR_cut = cms.bool(True),
                                   MVAId  = cms.VInputTag(),
                                   ID1_use_userFloat = cms.bool(False),
                                   # electronID1 = cms.string("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1"),
                                   electronID1 = cms.string("ElectronMVAEstimatorRun2Spring16GeneralPurposeV1"),
                                   electronID2 = cms.string("ElectronMVAEstimatorRun2Spring15Trig25nsV1"),
                                   electronID1_pass = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
                                   electronID2_pass = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),

                                   photonID1 = cms.string("PhotonMVAEstimatorRun2Spring15NonTrig25nsV2"),
                                   photonID2 = cms.string("TLEMVAEstimatorRun2Fall15V1"),

)

# process.tle_signal = tle_ntuple.clone()
# process.tle_signal.do_signal = cms.bool(True)

# process.tle_background = tle_ntuple.clone()
# process.tle_background.do_signal = cms.bool(False)

tle_ntuple.outputFile = cms.string("output")
tle_ntuple.do_TLE = cms.bool(False)

process.reg = tle_ntuple.clone()

# process.reg_signal = tle_ntuple.clone()
# process.reg_signal.do_signal = cms.bool(True)

# process.reg_background = tle_ntuple.clone()
# process.reg_background.do_signal = cms.bool(False)



# fileNameForSample = 'ntuple'


#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring('file:/data_CMS/cms/davignon/Trigger_WithThomas/CMSSW_7_6_0_pre7/src/L1Trigger/L1TNtuples/00BEC5EF-1472-E511-806C-02163E0141EA.root')
fileNames = cms.untracked.vstring(
    # DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/MINIAODSIM
    # '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/50000/024D0A89-34C4-E611-AEF9-008CFA111348.root',
    'file:../024D0A89-34C4-E611-AEF9-008CFA111348.root', # above  file

# '/store/data/Run2016H/DoubleEG/MINIAOD/03Feb2017_ver2-v1/50000/246AAF40-B3EA-E611-BC2F-0CC47A7DFEC4.root',

)

)

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
 
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v6', '') # MCRUN2_74_V8
# process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v16', '') # MCRUN2_74_V8

# process.TFileService = cms.Service("TFileService", fileName = cms.string(fileNameForSample + '.root') )

# process.p = cms.Path(process.egmPhotonIDSequence * process.electronMVAValueMapProducer * process.egmPhotonIsolationMiniAOD 
# * process.reg_signal * process.reg_background
# * process.tle_signal * process.tle_background)#*process.dump)
# process.p = cms.Path(process.egmPhotonIDSequence * process.electronMVAValueMapProducer * process.egmPhotonIsolationMiniAOD 
# * process.reg_signal * process.reg_background)
process.p = cms.Path(process.egmPhotonIDSequence * process.electronMVAValueMapProducer * process.egmPhotonIsolationMiniAOD * process.reg)
