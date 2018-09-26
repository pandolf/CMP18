import FWCore.ParameterSet.Config as cms

process = cms.Process("CODAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 20


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.TFileService = cms.Service("TFileService",
      #fileName = cms.string("root://eoscms.cern.ch//eos/cms/store/user/pandolf/prova/ganjaTree.root"),
      fileName = cms.string("codTree.root"),
      closeFileFast = cms.untracked.bool(True)
)


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/p/pandolf/public/cmsOpenData7TeV_QCD_Pt_1400to1800_CMSSW_5_3_32.root'
        #'file:/afs/cern.ch/work/p/pandolf/public/cmsOpenData7TeV_QCD_Pt_50to80_CMSSW_5_3_32.root'
    )
)

process.codana = cms.EDAnalyzer('CODAnalyzer'
)


process.p = cms.Path(process.codana)
