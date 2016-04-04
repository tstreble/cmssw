import FWCore.ParameterSet.Config as cms

l1JetRecoTree = cms.EDAnalyzer("L1JetRecoTreeProducer",
  pfJetToken              = cms.untracked.InputTag("ak4PFJetsCHS"),
  jecToken                = cms.untracked.InputTag("ak4PFCHSL1FastL2L3ResidualCorrector"),
  maxJet                  = cms.uint32(20),
  jetptThreshold          = cms.double(30),
  jetetaMax               = cms.double(3.),
  pfMetToken              = cms.untracked.InputTag("pfMetT1"),
<<<<<<< HEAD
  caloMetToken            = cms.untracked.InputTag("caloMet"),
  caloMetBEToken          = cms.untracked.InputTag("caloMetBE")
=======
  caloMetToken            = cms.untracked.InputTag("caloMet")
>>>>>>> 3bfd3dc... pfMet Type-1 correction applied and the addition of caloMet in the sums branch of ntuples
)


