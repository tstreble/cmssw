import FWCore.ParameterSet.Config as cms

from PhysicsTools.NanoAOD.common_cff import *

BToKeTrk=cms.EDProducer("BToKeTrkProducer",
                        beamSpot=cms.InputTag("offlineBeamSpot"),
                        electronCollection=cms.InputTag("slimmedElectronsWithUserData"),
                        generalTrackCollection=cms.InputTag("generalTracks"),
                        PFCandCollection=cms.InputTag("packedPFCandidates"),
                        ElectronMinPt=cms.double(1.),
                        ElectronMaxEta=cms.double(2.4),
                        TrackMinPt=cms.double(0.),
                        TrackMaxEta=cms.double(2.4),
                        KaonMinPt=cms.double(1.),
                        KaonMaxEta=cms.double(2.4),
                        KaonMinDCASig=cms.double(3.3),
                        DiElectronChargeCheck=cms.bool(True)
                        )

BToKeTrkTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
                             src=cms.InputTag("BToKeTrk"),
                             cut=cms.string(""),
                             name=cms.string("BToKeTrk"),
                             doc=cms.string("BToKeTrk Variable"),
                             singleton=cms.bool(False),
                             extension=cms.bool(False),
                             variables=cms.PSet(
                                ele_pt=Var("userFloat('ele_pt')", float,doc="pt of electron (refitted)"),
                                ele_eta=Var("userFloat('ele_eta')", float,doc="eta of electron (refitted)"),
                                ele_phi=Var("userFloat('ele_phi')", float,doc="phi of electron (refitted)"),
                                ele_charge=Var("userFloat('ele_charge')", int,doc="charge of electron"),
                                ele_dxy=Var("daughter('ele').dB('PV2D')", float,doc="dxy of electron (with sign) wrt first PV, in cm"),
                                ele_dz=Var("daughter('ele').dB('PVDZ')", float,doc="dz of electron (with sign) wrt first PV, in cm"),
                                trk_pt=Var("userFloat('trk_pt')", float,doc="pt of track (refitted)"),
                                trk_eta=Var("userFloat('trk_eta')", float,doc="eta of track (refitted)"),
                                trk_phi=Var("userFloat('trk_phi')", float,doc="phi of track (refitted)"),
                                trk_charge=Var("userFloat('trk_charge')", int,doc="charge of track"),
                                kaon_pt=Var("userFloat('kaon_pt')", float,doc="pt of kaon (refitted)"),
                                kaon_eta=Var("userFloat('kaon_eta')", float,doc="eta of kaon (refitted)"),
                                kaon_phi=Var("userFloat('kaon_phi')", float,doc="phi of kaon (refitted)"),
                                kaon_charge=Var("userFloat('kaon_charge')", int,doc="charge of kaon"),
                                kaon_DCASig=Var("userFloat('kaon_DCASig')", float,doc="significance of xy-distance of closest approach kaon-beamspot"),
                                kaon_dxy=Var("daughter('kaon').dxy()", float,doc="dxy of kaon (not refitted)"),
                                kaon_dz=Var("daughter('kaon').dz()", float,doc="dz of kaon (not refitted)"),
                                ee_pt=Var("userFloat('ee_pt')", float,doc="dielectron pt (refitted)"),
                                ee_eta=Var("userFloat('ee_eta')", float,doc="dielectron eta (refitted)"),
                                ee_phi=Var("userFloat('ee_phi')", float,doc="dielectron phi (refitted)"),
                                ee_mass=Var("userFloat('ee_mass')", float,doc="dielectron mass (refitted)"),
                                ee_mass_err=Var("userFloat('ee_mass_err')", float,doc="error on dielectron mass"),
                                ee_Lxy=Var("userFloat('ee_Lxy')", float,doc="significance of dielectron vertex-beamspot xy-separation"),
                                ee_CL_vtx=Var("userFloat('ee_CL_vtx')", float,doc="dielectron chi2 vertex probability"),
                                pt=Var("userFloat('pt')", float,doc="pt of BToKee candidate (refitted)"),
                                eta=Var("userFloat('eta')", float,doc="eta of BToKee candidate (refitted)"),
                                phi=Var("userFloat('phi')", float,doc="phi of BToKee candidate (refitted)"),
                                mass=Var("userFloat('mass')", float,doc="mass of BToKee candidate (refitted)"),
                                mass_err=Var("userFloat('mass_err')", float,doc="error on mass of BToKee candidate"),                                
                                Lxy=Var("userFloat('Lxy')", float,doc="significance of BToKee vertex-beamspot xy-separation"),
                                CL_vtx=Var("userFloat('CL_vtx')", float,doc="BToKee chi2 vertex probability"),
                                cosAlpha=Var("userFloat('cosAlpha')", float,doc="cosine of angle between BToKmumu momentum and vertex-beamspot separation"),
                                )
                             )

BToKeTrkSequence=cms.Sequence(BToKeTrk)
BToKeTrkTables=cms.Sequence(BToKeTrkTable)
