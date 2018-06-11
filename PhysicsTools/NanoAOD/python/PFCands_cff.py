from PhysicsTools.NanoAOD.common_cff import *
import FWCore.ParameterSet.Config as cms

PFCand=cms.EDProducer("PFCandProducer",
                      beamSpot=cms.InputTag("offlineBeamSpot"),
                      PFCandCollection=cms.InputTag("packedPFCandidates")
                      )


PFCandTable=cms.EDProducer("SimpleCandidateFlatTableProducer",
                             src=cms.InputTag("PFCand"),
                             cut=cms.string(""),
                             name=cms.string("PFCand"),
                             doc=cms.string("PF candidates Variable"),
                             singleton=cms.bool(False),
                             extension=cms.bool(False),
                             variables=cms.PSet(pt = Var("daughter(0).pt()",float,doc="pt"),
                                                eta = Var("daughter(0).eta()",float,doc="eta"),
                                                phi = Var("daughter(0).phi()",float,doc="phi"),
                                                mass = Var("daughter(0).mass()",float,doc="mass"),
                                                pdgId = Var("daughter(0).pdgId()",int,doc="PF pdgID"),
                                                DCASig=Var("userFloat('DCASig')", float,doc="significance of xy-distance of closest approach PFCand-beamspot"),
                                                dz = Var("daughter(0).dz()",float,doc="dz (with sign) wrt first PV, in cm",precision=10),
                                                dxy = Var("daughter(0).dxy()",float,doc="dxy (with sign) wrt first PV, in cm",precision=10),
                                                ),
                             )


PFCandSequence=cms.Sequence(PFCand)
PFCandTables=cms.Sequence(PFCandTable)
