#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include <TLorentzVector.h>
#include <TVector.h>
#include <TMatrix.h>


//
// class declaration
//

using namespace std;

class BToKstmumuProducer : public edm::EDProducer {
    
public:
    
    explicit BToKstmumuProducer(const edm::ParameterSet &iConfig);
    
    ~BToKstmumuProducer() override {};
    
    
private:
    
    virtual void produce(edm::Event&, const edm::EventSetup&);
    
    bool MuMuVertexRefitting(const pat::Muon & muon1,
			     const pat::Muon & muon2,
			     edm::ESHandle<TransientTrackBuilder> theTTBuilder,
			     RefCountedKinematicVertex &refitVertex,
			     RefCountedKinematicParticle &refitMuMu,
			     RefCountedKinematicParticle &refitMu1,
			     RefCountedKinematicParticle &refitMu2);

    bool KstVertexRefitting(const pat::PackedCandidate &kaon,
			    const pat::PackedCandidate &pion,
			    edm::ESHandle<TransientTrackBuilder> theTTBuilder,
			    RefCountedKinematicVertex &refitVertex,
			    RefCountedKinematicParticle &refitKst,
			    RefCountedKinematicParticle &refitKaon,
			    RefCountedKinematicParticle &refitPion);

    bool BToKstMuMuVertexRefitting(const pat::Muon &muon1,
                                   const pat::Muon &muon2,
                                   const RefCountedKinematicParticle refitKPi,
                                   edm::ESHandle<TransientTrackBuilder> theTTBuilder,
                                   RefCountedKinematicVertex &refitVertex,
                                   RefCountedKinematicParticle &refitBToKstMuMu,
                                   RefCountedKinematicParticle &refitMu1,
                                   RefCountedKinematicParticle &refitMu2,
                                   RefCountedKinematicParticle &refitKst);

    bool BToKPiMuMuVertexRefitting(const pat::Muon &muon1,
                                   const pat::Muon &muon2,
                                   const pat::PackedCandidate &kaon,
                                   const pat::PackedCandidate &pion,
                                   edm::ESHandle<TransientTrackBuilder> theTTBuilder,
                                   RefCountedKinematicVertex &refitVertex,
                                   RefCountedKinematicParticle &refitBToKstMuMu,
                                   RefCountedKinematicParticle &refitMu1,
                                   RefCountedKinematicParticle &refitMu2,
                                   RefCountedKinematicParticle &refitKaon,
                                   RefCountedKinematicParticle &refitPion);

    bool BToKstJPsiMuMuVertexRefitting(const RefCountedKinematicParticle refitMuMu,
                                       const RefCountedKinematicParticle refitKPi,
                                       RefCountedKinematicVertex &refitVertex,
                                       RefCountedKinematicParticle &refitBToKstJPsiMuMu,
                                       RefCountedKinematicParticle &refitJPsi,
                                       RefCountedKinematicParticle &refitKst);

    
    pair<double,double> computeLS(RefCountedKinematicVertex refitVertex,
                                  reco::BeamSpot beamSpot);
    
    double computeCosAlpha(RefCountedKinematicParticle refitBToKstMuMu,
                           RefCountedKinematicVertex vertexFitTree,
                           reco::BeamSpot beamSpot);

    pair<double,double> computeDCA(const pat::PackedCandidate &kaon,
                                   edm::ESHandle<MagneticField> bFieldHandle,
                                   reco::BeamSpot beamSpot);
    
    // ----------member data ---------------------------
    
    edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;
    edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
    edm::EDGetTokenT<std::vector<pat::Muon>> muonSrc_;
    edm::EDGetTokenT<edm::View<pat::PackedCandidate>> PFCandSrc_;
    
    double ptMinMu_;
    double etaMaxMu_;
    double ptMinKaon_;
    double etaMaxKaon_;
    double DCASigMinKaon_;
    double ptMinPion_;
    double etaMaxPion_;
    double DCASigMinPion_;
    bool diMuonCharge_;
    bool KstCharge_;
    double JPsiMassConstraint_;
    double KstMassConstraint_;
    bool save2TrkRefit_;
    bool save4TrkRefit_;
    
    float MuonMass_ = 0.10565837;
    float MuonMassErr_ = 3.5*1e-9;
    float KaonMass_ = 0.493677;
    float KaonMassErr_ = 1.6e-5;
    float PionMass_ = 0.139570;
    float PionMassErr_ = 3.5e-7;
    //float JPsiMass_ = 3.096916;  //Configurable parameter
    float JPsiMassErr_ = 0.011;
    //float KstMass_ = 0.89176;  //Configurable parameter
    float KstMassErr_ = 0.25e-3;
  
};



BToKstmumuProducer::BToKstmumuProducer(const edm::ParameterSet &iConfig):
beamSpotSrc_( consumes<reco::BeamSpot> ( iConfig.getParameter<edm::InputTag>( "beamSpot" ) ) ),
vertexSrc_( consumes<reco::VertexCollection> ( iConfig.getParameter<edm::InputTag>( "vertexCollection" ) ) ),
muonSrc_( consumes<std::vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "muonCollection" ) ) ),
PFCandSrc_( consumes<edm::View<pat::PackedCandidate>> ( iConfig.getParameter<edm::InputTag>( "PFCandCollection" ) ) ),
ptMinMu_( iConfig.getParameter<double>( "MuonMinPt" ) ),
etaMaxMu_( iConfig.getParameter<double>( "MuonMaxEta" ) ),
ptMinKaon_( iConfig.getParameter<double>( "KaonMinPt" ) ),
etaMaxKaon_( iConfig.getParameter<double>( "KaonMaxEta" ) ),
DCASigMinKaon_( iConfig.getParameter<double>( "KaonMinDCASig" ) ),
ptMinPion_( iConfig.getParameter<double>( "PionMinPt" ) ),
etaMaxPion_( iConfig.getParameter<double>( "PionMaxEta" ) ),
DCASigMinPion_( iConfig.getParameter<double>( "PionMinDCASig" ) ),
diMuonCharge_( iConfig.getParameter<bool>( "DiMuonChargeCheck" ) ),
KstCharge_( iConfig.getParameter<bool>( "KstarChargeCheck" ) ),
JPsiMassConstraint_( iConfig.getParameter<double>( "JPsiMassConstraint" ) ),
KstMassConstraint_( iConfig.getParameter<double>( "KstMassConstraint" ) ),
save2TrkRefit_( iConfig.getParameter<bool>( "save2TrackRefit" ) ),
save4TrkRefit_( iConfig.getParameter<bool>( "save4TrackRefit" ) )
{
    produces<pat::CompositeCandidateCollection>();
}


void BToKstmumuProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::ESHandle<MagneticField> bFieldHandle;
    edm::ESHandle<TransientTrackBuilder> theTTBuilder;

    iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);

    edm::Handle<reco::BeamSpot> beamSpotHandle;
    iEvent.getByToken(beamSpotSrc_, beamSpotHandle);
    if ( ! beamSpotHandle.isValid() ) {
        edm::LogError("BToKstmumuProducer") << "No beam spot available from EventSetup" ;
    }
    reco::BeamSpot beamSpot = *beamSpotHandle;

    edm::Handle<reco::VertexCollection> vertexHandle;
    iEvent.getByToken(vertexSrc_, vertexHandle);
    const reco::Vertex & PV = vertexHandle->front();

    edm::Handle<std::vector<pat::Muon>> muonHandle;
    edm::Handle<edm::View<pat::PackedCandidate>> pfCandHandle;

    iEvent.getByToken(muonSrc_, muonHandle);
    iEvent.getByToken(PFCandSrc_, pfCandHandle);

    unsigned int muonNumber = muonHandle->size();
    unsigned int pfCandNumber = pfCandHandle->size();

    // Output collection
    std::unique_ptr<pat::CompositeCandidateCollection> result( new pat::CompositeCandidateCollection );
    
    if(muonNumber>1){

        // loop on all the mumuKpi quadruplets
        for (unsigned int i = 0; i < muonNumber; ++i) {

            const pat::Muon & muon1 = (*muonHandle)[i];

            if(!(muon1.isLooseMuon() && muon1.isSoftMuon(PV))) continue;
            if(muon1.pt()<ptMinMu_ || abs(muon1.eta())>etaMaxMu_) continue;

            for (unsigned int j = 0; j < muonNumber; ++j) {

                if(i==j) continue;

		const pat::Muon & muon2 = (*muonHandle)[j];

		if(muon1.pt()<muon2.pt()) continue; //Muon 1 is always saved as the leading one
                if(!(muon2.isLooseMuon() && muon2.isSoftMuon(PV))) continue;
                if(muon2.pt()<ptMinMu_ || abs(muon2.eta())>etaMaxMu_) continue;
                if(diMuonCharge_ && muon1.charge()*muon2.charge()>0) continue;

                bool passedDiMuon = false;

                double MuMuLSBS = -1.;
                double MuMuLSBSErr = -1.;
                double MuMuVtx_CL = -1.;
                double MuMu_mass_err = -1.;

                RefCountedKinematicParticle refitMuMu;
                math::XYZVector refitMuMuV3D;

		if(save2TrkRefit_){

                  RefCountedKinematicVertex refitVertexMuMu;
                  RefCountedKinematicParticle refitMu1_MuMu;
                  RefCountedKinematicParticle refitMu2_MuMu;

                  passedDiMuon = MuMuVertexRefitting(muon1, muon2,
						     theTTBuilder,
						     refitVertexMuMu,
						     refitMuMu,
						     refitMu1_MuMu,
						     refitMu2_MuMu);

                  if (passedDiMuon){

		    pair<double,double> MuMuLS = computeLS(refitVertexMuMu,beamSpot);
		    MuMuLSBS = MuMuLS.first;
		    MuMuLSBSErr = MuMuLS.second;

		    MuMuVtx_CL = TMath::Prob((double)refitVertexMuMu->chiSquared(),
					     int(rint(refitVertexMuMu->degreesOfFreedom())));

		    MuMu_mass_err = sqrt(refitMuMu->currentState().kinematicParametersError().matrix()(6,6));

		    math::XYZVector refitMu1V3D_MuMu = refitMu1_MuMu->refittedTransientTrack().track().momentum();
		    math::XYZVector refitMu2V3D_MuMu = refitMu2_MuMu->refittedTransientTrack().track().momentum();
		    refitMuMuV3D = refitMu1V3D_MuMu + refitMu2V3D_MuMu;

		  }

		}


                //Kaon
                for (unsigned int k = 0; k < pfCandNumber; ++k) {

                    const pat::PackedCandidate & kaon = (*pfCandHandle)[k];
                    if(abs(kaon.pdgId())!=211) continue; //Charged hadrons
                    if(!kaon.hasTrackDetails()) continue;
                    if(kaon.pt()<ptMinKaon_ || abs(kaon.eta())>etaMaxKaon_) continue;

                    pair<double,double> DCA_kaon = computeDCA(kaon,
							      bFieldHandle,
							      beamSpot);
                    double DCABS_kaon = DCA_kaon.first;
                    double DCABSErr_kaon = DCA_kaon.second;

                    if(fabs(DCABS_kaon/DCABSErr_kaon)<DCASigMinKaon_) continue;
		    
                    for (unsigned int l = 0; l < pfCandNumber; ++l) {

		      if(k==l) continue;

		      const pat::PackedCandidate & pion = (*pfCandHandle)[l];
		      if(abs(pion.pdgId())!=211) continue; //Charged hadrons
		      if(!pion.hasTrackDetails()) continue;
		      if(pion.pt()<ptMinPion_ || abs(pion.eta())>etaMaxPion_) continue;
		      if(KstCharge_ && kaon.charge()*pion.charge()>0) continue;

		      pair<double,double> DCA_pion = computeDCA(pion,
								bFieldHandle,
								beamSpot);
		      double DCABS_pion = DCA_pion.first;
		      double DCABSErr_pion = DCA_pion.second;
		      
		      if(fabs(DCABS_pion/DCABSErr_pion)<DCASigMinPion_) continue;

		      RefCountedKinematicParticle refitKst;
		      RefCountedKinematicVertex refitVertexKst;
		      RefCountedKinematicParticle refitKaon_Kst;
		      RefCountedKinematicParticle refitPion_Kst;

		      bool passed = KstVertexRefitting(kaon, pion,
						       theTTBuilder,
						       refitVertexKst,
						       refitKst,
						       refitKaon_Kst,
						       refitPion_Kst);

		      if(!passed) continue;

		      pair<double,double> KstLS = computeLS(refitVertexKst,beamSpot);
		      double KstLSBS = KstLS.first;
		      double KstLSBSErr = KstLS.second;

		      double KstVtx_CL = TMath::Prob((double)refitVertexKst->chiSquared(),
						     int(rint(refitVertexKst->degreesOfFreedom())));

		      double Kst_mass_err = sqrt(refitKst->currentState().kinematicParametersError().matrix()(6,6));

		      math::XYZVector refitKaonV3D_Kst = refitKaon_Kst->refittedTransientTrack().track().momentum();
		      math::XYZVector refitPionV3D_Kst = refitPion_Kst->refittedTransientTrack().track().momentum();
		      math::XYZVector refitKstV3D = refitKaonV3D_Kst + refitPionV3D_Kst;

		      RefCountedKinematicVertex refitVertexBToKstMuMu;
		      RefCountedKinematicParticle refitBToKstMuMu;
		      RefCountedKinematicParticle refitMu1;
		      RefCountedKinematicParticle refitMu2;
		      RefCountedKinematicParticle refitKst_BToKstMuMu;

		      passed = BToKstMuMuVertexRefitting(muon1, muon2, refitKst,
							 theTTBuilder,
							 refitVertexBToKstMuMu,
							 refitBToKstMuMu,
							 refitMu1,
							 refitMu2,
							 refitKst_BToKstMuMu);

		      if (!passed) continue;

		      pair<double,double> BToKstMuMuLS = computeLS(refitVertexBToKstMuMu,beamSpot);
		      double LSBS = BToKstMuMuLS.first;
		      double LSBSErr = BToKstMuMuLS.second;

		      double BToKstMuMuVtx_CL = TMath::Prob((double)refitVertexBToKstMuMu->chiSquared(),
							  int(rint(refitVertexBToKstMuMu->degreesOfFreedom())));

		      double cosAlpha = computeCosAlpha(refitBToKstMuMu,refitVertexBToKstMuMu,beamSpot);

		      double mass_err = sqrt(refitBToKstMuMu->currentState().kinematicParametersError().matrix()(6,6));

		      math::XYZVector refitMuon1V3D = refitMu1->refittedTransientTrack().track().momentum();
		      math::XYZVector refitMuon2V3D = refitMu2->refittedTransientTrack().track().momentum();
		      math::XYZVector refitKst_BToKstMuMu_V3D = refitKst_BToKstMuMu->refittedTransientTrack().track().momentum();
		      math::XYZVector refitBToKstMuMuV3D = refitMuon1V3D + refitMuon2V3D + refitKst_BToKstMuMu_V3D;

		      pat::CompositeCandidate BToKstMuMuCand;
		      BToKstMuMuCand.addDaughter( muon1, "muon1");
		      BToKstMuMuCand.addDaughter( muon2, "muon2");
		      BToKstMuMuCand.addDaughter( kaon, "kaon");
		      BToKstMuMuCand.addDaughter( pion, "pion");

		      BToKstMuMuCand.addUserInt("mu1_index", i);
		      BToKstMuMuCand.addUserInt("mu2_index", j);
		      BToKstMuMuCand.addUserInt("kaon_index", k);
		      BToKstMuMuCand.addUserInt("pion_index", l);
		      
		      BToKstMuMuCand.addUserFloat("mu1_pt",     sqrt(refitMuon1V3D.perp2()));
		      BToKstMuMuCand.addUserFloat("mu1_eta",    refitMuon1V3D.eta());
		      BToKstMuMuCand.addUserFloat("mu1_phi",    refitMuon1V3D.phi());
		      BToKstMuMuCand.addUserFloat("mu1_charge", refitMu1->currentState().particleCharge());

		      BToKstMuMuCand.addUserFloat("mu2_pt",     sqrt(refitMuon2V3D.perp2()));
		      BToKstMuMuCand.addUserFloat("mu2_eta",    refitMuon2V3D.eta());
		      BToKstMuMuCand.addUserFloat("mu2_phi",    refitMuon2V3D.phi());
		      BToKstMuMuCand.addUserFloat("mu2_charge", refitMu2->currentState().particleCharge());		    

		      TLorentzVector muon1cand;
		      muon1cand.SetPtEtaPhiM(sqrt(refitMuon1V3D.perp2()), refitMuon1V3D.eta(), refitMuon1V3D.phi(), MuonMass_);
		      TLorentzVector muon2cand;
		      muon2cand.SetPtEtaPhiM(sqrt(refitMuon2V3D.perp2()), refitMuon2V3D.eta(), refitMuon2V3D.phi(), MuonMass_);
		      BToKstMuMuCand.addUserFloat("mumuKPiFit_mumu_mass", (muon1cand+muon2cand).Mag());

		      BToKstMuMuCand.addUserFloat("kaon_pt",    sqrt(refitKaonV3D_Kst.perp2()));
		      BToKstMuMuCand.addUserFloat("kaon_eta",   refitKaonV3D_Kst.eta());
		      BToKstMuMuCand.addUserFloat("kaon_phi",   refitKaonV3D_Kst.phi());
		      BToKstMuMuCand.addUserFloat("kaon_charge",refitKaon_Kst->currentState().particleCharge());
		      BToKstMuMuCand.addUserFloat("kaon_DCASig", DCABS_kaon/DCABSErr_kaon);
		      
		      BToKstMuMuCand.addUserFloat("pion_pt",    sqrt(refitPionV3D_Kst.perp2()));
		      BToKstMuMuCand.addUserFloat("pion_eta",   refitPionV3D_Kst.eta());
		      BToKstMuMuCand.addUserFloat("pion_phi",   refitPionV3D_Kst.phi());
		      BToKstMuMuCand.addUserFloat("pion_charge",refitPion_Kst->currentState().particleCharge());
		      BToKstMuMuCand.addUserFloat("pion_DCASig", DCABS_pion/DCABSErr_pion);

		      BToKstMuMuCand.addUserFloat("Kst_pt", sqrt(refitKstV3D.perp2()));
		      BToKstMuMuCand.addUserFloat("Kst_eta", refitKstV3D.eta());
		      BToKstMuMuCand.addUserFloat("Kst_phi", refitKstV3D.phi());
		      BToKstMuMuCand.addUserFloat("Kst_mass", refitKst->currentState().mass());
		      BToKstMuMuCand.addUserFloat("Kst_mass_err", Kst_mass_err);
		      BToKstMuMuCand.addUserFloat("Kst_Lxy", (float) KstLSBS/KstLSBSErr);
		      BToKstMuMuCand.addUserFloat("Kst_ctxy", (float) KstLSBS/sqrt(refitKstV3D.perp2()));
		      BToKstMuMuCand.addUserFloat("Kst_CL_vtx", (float) KstVtx_CL);

		      BToKstMuMuCand.addUserFloat("pt",     sqrt(refitBToKstMuMuV3D.perp2()));
		      BToKstMuMuCand.addUserFloat("eta",    refitBToKstMuMuV3D.eta());
		      BToKstMuMuCand.addUserFloat("phi",    refitBToKstMuMuV3D.phi());
		      BToKstMuMuCand.addUserFloat("mass",   refitBToKstMuMu->currentState().mass());
		      BToKstMuMuCand.addUserFloat("mass_err", mass_err);
		      BToKstMuMuCand.addUserFloat("Lxy", (float) LSBS/LSBSErr);
		      BToKstMuMuCand.addUserFloat("ctxy", (float) LSBS/sqrt(refitBToKstMuMuV3D.perp2()));
		      BToKstMuMuCand.addUserFloat("CL_vtx", (float) BToKstMuMuVtx_CL);
		      BToKstMuMuCand.addUserFloat("cosAlpha", (float) cosAlpha);
                    
		      BToKstMuMuCand.addUserInt("mumuRefit", (int)passedDiMuon);
		      BToKstMuMuCand.addUserFloat("mumu_pt", (passedDiMuon)? sqrt(refitMuMuV3D.perp2()) : -1.);
		      BToKstMuMuCand.addUserFloat("mumu_eta", (passedDiMuon)? refitMuMuV3D.eta() : -9.);
		      BToKstMuMuCand.addUserFloat("mumu_phi", (passedDiMuon)? refitMuMuV3D.phi() : -9.);
		      BToKstMuMuCand.addUserFloat("mumu_mass", (passedDiMuon)? refitMuMu->currentState().mass() : -1.);
		      BToKstMuMuCand.addUserFloat("mumu_mass_err", (passedDiMuon)?  MuMu_mass_err : -1.);
		      BToKstMuMuCand.addUserFloat("mumu_Lxy", (passedDiMuon)? (float) MuMuLSBS/MuMuLSBSErr : -1.);
		      BToKstMuMuCand.addUserFloat("mumu_ctxy", (passedDiMuon)? (float) MuMuLSBS/sqrt(refitMuMuV3D.perp2()) : -1.);
		      BToKstMuMuCand.addUserFloat("mumu_CL_vtx", (passedDiMuon)? (float) MuMuVtx_CL : -1.);

		      bool passed_2trk = false;
		      float pt_2trk = -9999.;
		      float eta_2trk = -9999.;
		      float phi_2trk = -9999.;
		      float mass_2trk = -9999.;
		      float mass_err_2trk = -9999.;
		      float Lxy_2trk = -9999.;
		      float ctxy_2trk = -9999.;
		      float CL_vtx_2trk = -9999.;
		      float cosAlpha_2trk = -9999.;

		      //This is crashing for some unexplained reason
		      //https://hypernews.cern.ch/HyperNews/CMS/get/physTools/2746/1/2/1/1.html
		      //Disabled for now

		      //if(save2TrkRefit_ && passedDiMuon){
		      if(0){

			RefCountedKinematicVertex refitVertexBToKstJPsiMuMu;
			RefCountedKinematicParticle refitBToKstJPsiMuMu;
			RefCountedKinematicParticle refitJPsiMuMu;
			RefCountedKinematicParticle refitKst_KstJPsi;

			passed_2trk = BToKstJPsiMuMuVertexRefitting(refitMuMu, refitKst,
								    refitVertexBToKstJPsiMuMu,
								    refitBToKstJPsiMuMu,
								    refitJPsiMuMu,
								    refitKst_KstJPsi);

			if(passed_2trk){

			  math::XYZVector refitJPsiMuMuV3D = refitJPsiMuMu->refittedTransientTrack().track().momentum();
			  math::XYZVector refitKstV3D_KJPsi = refitKst_KstJPsi->refittedTransientTrack().track().momentum();
			  math::XYZVector refitBToKstJPsiMuMuV3D = refitJPsiMuMuV3D + refitKstV3D_KJPsi;
			  
			  pt_2trk = sqrt(refitBToKstJPsiMuMuV3D.perp2());
			  eta_2trk = refitBToKstJPsiMuMuV3D.eta();
			  phi_2trk = refitBToKstJPsiMuMuV3D.phi();
			  mass_2trk = refitBToKstJPsiMuMu->currentState().mass();

			  mass_err_2trk = sqrt(refitBToKstJPsiMuMu->currentState().kinematicParametersError().matrix()(6,6));

			  pair<double,double> BToKstJPsiMuMuLS = computeLS(refitVertexBToKstJPsiMuMu,beamSpot);
			  double LSBS_2trk = BToKstJPsiMuMuLS.first;
			  double LSBSErr_2trk = BToKstJPsiMuMuLS.second;
			  Lxy_2trk = LSBS_2trk/LSBSErr_2trk;
			  ctxy_2trk = LSBS_2trk/pt_2trk;
			  CL_vtx_2trk = TMath::Prob((double)refitVertexBToKstJPsiMuMu->chiSquared(),
						    int(rint(refitVertexBToKstJPsiMuMu->degreesOfFreedom())));
			  cosAlpha_2trk = computeCosAlpha(refitBToKstJPsiMuMu,refitVertexBToKstJPsiMuMu,beamSpot);

			}

		      }

                     BToKstMuMuCand.addUserInt("2trkRefit", (int)passed_2trk);
                     BToKstMuMuCand.addUserFloat("pt_2trk", pt_2trk);
                     BToKstMuMuCand.addUserFloat("eta_2trk", eta_2trk);
                     BToKstMuMuCand.addUserFloat("phi_2trk", phi_2trk);
                     BToKstMuMuCand.addUserFloat("mass_2trk", mass_2trk);
                     BToKstMuMuCand.addUserFloat("mass_err_2trk", mass_err_2trk);
                     BToKstMuMuCand.addUserFloat("Lxy_2trk", Lxy_2trk);
                     BToKstMuMuCand.addUserFloat("ctxy_2trk", ctxy_2trk);
                     BToKstMuMuCand.addUserFloat("CL_vtx_2trk", CL_vtx_2trk);
                     BToKstMuMuCand.addUserFloat("cosAlpha_2trk", cosAlpha_2trk);

                     bool passed_4trk = false;
                     float pt_4trk = -9999.;
                     float eta_4trk = -9999.;
                     float phi_4trk = -9999.;
                     float mass_4trk = -9999.;
                     float mass_err_4trk = -9999.;
                     float Lxy_4trk = -9999.;
                     float ctxy_4trk = -9999.;
                     float CL_vtx_4trk = -9999.;
                     float cosAlpha_4trk = -9999.;

                     if(save4TrkRefit_){

		       RefCountedKinematicVertex refitVertexBToKPiMuMu;
		       RefCountedKinematicParticle refitBToKPiMuMu;
		       RefCountedKinematicParticle refitMu1_BToKPiMuMu;
		       RefCountedKinematicParticle refitMu2_BToKPiMuMu;
		       RefCountedKinematicParticle refitKaon_BToKPiMuMu;
		       RefCountedKinematicParticle refitPion_BToKPiMuMu;
		       
		       passed_4trk = BToKPiMuMuVertexRefitting(muon1, muon2, kaon, pion,
							       theTTBuilder,
							       refitVertexBToKPiMuMu,
							       refitBToKPiMuMu,
							       refitMu1_BToKPiMuMu,
							       refitMu2_BToKPiMuMu,
							       refitKaon_BToKPiMuMu,
							       refitPion_BToKPiMuMu);

                       if(passed_4trk){

                         math::XYZVector refitMu1_BToKPiMuMu_V3D = refitMu1_BToKPiMuMu->refittedTransientTrack().track().momentum();
                         math::XYZVector refitMu2_BToKPiMuMu_V3D = refitMu2_BToKPiMuMu->refittedTransientTrack().track().momentum();
			 math::XYZVector refitKaon_BToKPiMuMu_V3D = refitKaon_BToKPiMuMu->refittedTransientTrack().track().momentum();
			 math::XYZVector refitPion_BToKPiMuMu_V3D = refitPion_BToKPiMuMu->refittedTransientTrack().track().momentum();
			 math::XYZVector refitBToKPiMuMuV3D = refitMu1_BToKPiMuMu_V3D + refitMu2_BToKPiMuMu_V3D + refitKaon_BToKPiMuMu_V3D + refitPion_BToKPiMuMu_V3D;

			 pt_4trk = sqrt(refitBToKPiMuMuV3D.perp2());
			 eta_4trk = refitBToKPiMuMuV3D.eta();
			 phi_4trk = refitBToKPiMuMuV3D.phi();
			 mass_4trk = refitBToKPiMuMu->currentState().mass();
			 mass_err_4trk = sqrt(refitBToKPiMuMu->currentState().kinematicParametersError().matrix()(6,6));

			 pair<double,double> BToKPiMuMuLS = computeLS(refitVertexBToKPiMuMu,beamSpot);
			 double LSBS_4trk = BToKPiMuMuLS.first;
			 double LSBSErr_4trk = BToKPiMuMuLS.second;
			 Lxy_4trk = LSBS_4trk/LSBSErr_4trk;
			 ctxy_4trk = LSBS_4trk/pt_4trk;
			 CL_vtx_4trk = TMath::Prob((double)refitVertexBToKPiMuMu->chiSquared(),
						   int(rint(refitVertexBToKPiMuMu->degreesOfFreedom())));
			 cosAlpha_4trk = computeCosAlpha(refitBToKPiMuMu,refitVertexBToKPiMuMu,beamSpot);

		       }

                     }

                     BToKstMuMuCand.addUserInt("4trkRefit", (int)passed_4trk);
                     BToKstMuMuCand.addUserFloat("pt_4trk", pt_4trk);
                     BToKstMuMuCand.addUserFloat("eta_4trk", eta_4trk);
                     BToKstMuMuCand.addUserFloat("phi_4trk", phi_4trk);
                     BToKstMuMuCand.addUserFloat("mass_4trk", mass_4trk);
                     BToKstMuMuCand.addUserFloat("mass_err_4trk", mass_err_4trk);
                     BToKstMuMuCand.addUserFloat("Lxy_4trk", Lxy_4trk);
                     BToKstMuMuCand.addUserFloat("ctxy_4trk", ctxy_4trk);
                     BToKstMuMuCand.addUserFloat("CL_vtx_4trk", CL_vtx_4trk);
                     BToKstMuMuCand.addUserFloat("cosAlpha_4trk", cosAlpha_4trk);

		     result->push_back(BToKstMuMuCand);
                    
		    }
                    
                }
                
            }
            
        }
        
    }
    
    iEvent.put(std::move(result));
    
}



bool BToKstmumuProducer::MuMuVertexRefitting(const pat::Muon & muon1,
					     const pat::Muon & muon2,
					     edm::ESHandle<TransientTrackBuilder> theTTBuilder,
					     RefCountedKinematicVertex &refitVertex,
					     RefCountedKinematicParticle &refitMuMu,
					     RefCountedKinematicParticle &refitMuon1,
					     RefCountedKinematicParticle &refitMuon2){
    
    const reco::TransientTrack muon1TT = theTTBuilder->build(muon1.innerTrack());
    const reco::TransientTrack muon2TT = theTTBuilder->build(muon2.innerTrack());
    
    KinematicParticleFactoryFromTransientTrack partFactory;
    KinematicParticleVertexFitter PartVtxFitter;
    
    std::vector<RefCountedKinematicParticle> muParticles;
    double chi = 0.;
    double ndf = 0.;
    muParticles.push_back(partFactory.particle(muon1TT,MuonMass_,chi,ndf,MuonMassErr_));
    muParticles.push_back(partFactory.particle(muon2TT,MuonMass_,chi,ndf,MuonMassErr_));
    RefCountedKinematicTree mumuVertexFitTree = PartVtxFitter.fit(muParticles);
    
    if ( !mumuVertexFitTree->isValid()) return false;
    
    mumuVertexFitTree->movePointerToTheTop();
    refitVertex = mumuVertexFitTree->currentDecayVertex();
    refitMuMu = mumuVertexFitTree->currentParticle();
    
    if ( !refitVertex->vertexIsValid()) return false;

    // extract the re-fitted tracks
    mumuVertexFitTree->movePointerToTheTop();
    
    mumuVertexFitTree->movePointerToTheFirstChild();
    refitMuon1 = mumuVertexFitTree->currentParticle();
    
    mumuVertexFitTree->movePointerToTheNextChild();
    refitMuon2 = mumuVertexFitTree->currentParticle();
    
    return true;
    
}




bool BToKstmumuProducer::KstVertexRefitting(const pat::PackedCandidate &kaon,
					    const pat::PackedCandidate &pion,
					    edm::ESHandle<TransientTrackBuilder> theTTBuilder,
					    RefCountedKinematicVertex &refitVertex,
					    RefCountedKinematicParticle &refitKst,
					    RefCountedKinematicParticle &refitKaon,
					    RefCountedKinematicParticle &refitPion){
    
    const reco::TransientTrack kaonTT = theTTBuilder->build(kaon.bestTrack());
    const reco::TransientTrack pionTT = theTTBuilder->build(pion.bestTrack());

    KinematicParticleFactoryFromTransientTrack partFactory;
    KinematicParticleVertexFitter PartVtxFitter;
    
    std::vector<RefCountedKinematicParticle> KstParticles;
    double chi = 0.;
    double ndf = 0.;
    KstParticles.push_back(partFactory.particle(kaonTT,KaonMass_,chi,ndf,KaonMassErr_));
    KstParticles.push_back(partFactory.particle(pionTT,PionMass_,chi,ndf,PionMassErr_));
    RefCountedKinematicTree KstVertexFitTree = PartVtxFitter.fit(KstParticles);
    
    if ( !KstVertexFitTree->isValid()) return false;
    
    KstVertexFitTree->movePointerToTheTop();
    refitVertex = KstVertexFitTree->currentDecayVertex();
    refitKst = KstVertexFitTree->currentParticle();
    
    if ( !refitVertex->vertexIsValid()) return false;

    // extract the re-fitted tracks
    KstVertexFitTree->movePointerToTheTop();
    
    KstVertexFitTree->movePointerToTheFirstChild();
    refitKaon = KstVertexFitTree->currentParticle();
    
    KstVertexFitTree->movePointerToTheNextChild();
    refitPion = KstVertexFitTree->currentParticle();
    
    return true;
    
}




bool BToKstmumuProducer::BToKstMuMuVertexRefitting(const pat::Muon &muon1,
						   const pat::Muon &muon2,
						   const RefCountedKinematicParticle refitKPi,
						   edm::ESHandle<TransientTrackBuilder> theTTBuilder,
						   RefCountedKinematicVertex &refitVertex,
						   RefCountedKinematicParticle &refitBToKstMuMu,
						   RefCountedKinematicParticle &refitMuon1,
						   RefCountedKinematicParticle &refitMuon2,
						   RefCountedKinematicParticle &refitKst){

    const reco::TransientTrack muon1TT = theTTBuilder->build(muon1.innerTrack());
    const reco::TransientTrack muon2TT = theTTBuilder->build(muon2.innerTrack());
    const reco::TransientTrack KPiTT = refitKPi->refittedTransientTrack();

    KinematicParticleFactoryFromTransientTrack partFactory;
    KinematicParticleVertexFitter PartVtxFitter;

    float Kst_mass = refitKPi->currentState().mass();
    float Kst_mass_err = sqrt(refitKPi->currentState().kinematicParametersError().matrix()(6,6));
    if(KstMassConstraint_ > 0){
      Kst_mass = KstMassConstraint_;
      Kst_mass_err = KstMassErr_;
    }

    std::vector<RefCountedKinematicParticle> BToKstMuMuParticles;
    double chi = 0.;
    double ndf = 0.;
    BToKstMuMuParticles.push_back(partFactory.particle(muon1TT,MuonMass_,chi,ndf,MuonMassErr_));
    BToKstMuMuParticles.push_back(partFactory.particle(muon2TT,MuonMass_,chi,ndf,MuonMassErr_));
    BToKstMuMuParticles.push_back(partFactory.particle(KPiTT,Kst_mass,chi,ndf,Kst_mass_err));

    RefCountedKinematicTree BToKstMuMuVertexFitTree = PartVtxFitter.fit(BToKstMuMuParticles);

    if ( !BToKstMuMuVertexFitTree->isValid()) return false;

    BToKstMuMuVertexFitTree->movePointerToTheTop();
    refitVertex = BToKstMuMuVertexFitTree->currentDecayVertex();
    refitBToKstMuMu = BToKstMuMuVertexFitTree->currentParticle();

    if ( !refitVertex->vertexIsValid()) return false;

    // extract the re-fitted tracks
    BToKstMuMuVertexFitTree->movePointerToTheTop();

    BToKstMuMuVertexFitTree->movePointerToTheFirstChild();
    refitMuon1 = BToKstMuMuVertexFitTree->currentParticle();

    BToKstMuMuVertexFitTree->movePointerToTheNextChild();
    refitMuon2 = BToKstMuMuVertexFitTree->currentParticle();

    BToKstMuMuVertexFitTree->movePointerToTheNextChild();
    refitKst = BToKstMuMuVertexFitTree->currentParticle();

    return true;

}




bool BToKstmumuProducer::BToKPiMuMuVertexRefitting(const pat::Muon &muon1,
						   const pat::Muon &muon2,
						   const pat::PackedCandidate &kaon,
						   const pat::PackedCandidate &pion,
						   edm::ESHandle<TransientTrackBuilder> theTTBuilder,					   
						   RefCountedKinematicVertex &refitVertex,
						   RefCountedKinematicParticle &refitBToKstMuMu,
						   RefCountedKinematicParticle &refitMuon1,
						   RefCountedKinematicParticle &refitMuon2,
						   RefCountedKinematicParticle &refitKaon,
						   RefCountedKinematicParticle &refitPion){

    const reco::TransientTrack muon1TT = theTTBuilder->build(muon1.innerTrack());
    const reco::TransientTrack muon2TT = theTTBuilder->build(muon2.innerTrack());
    const reco::TransientTrack kaonTT = theTTBuilder->build(kaon.bestTrack());
    const reco::TransientTrack pionTT = theTTBuilder->build(pion.bestTrack());

    KinematicParticleFactoryFromTransientTrack partFactory;
    KinematicParticleVertexFitter PartVtxFitter;

    std::vector<RefCountedKinematicParticle> BToKstMuMuParticles;
    double chi = 0.;
    double ndf = 0.;
    BToKstMuMuParticles.push_back(partFactory.particle(muon1TT,MuonMass_,chi,ndf,MuonMassErr_));
    BToKstMuMuParticles.push_back(partFactory.particle(muon2TT,MuonMass_,chi,ndf,MuonMassErr_));
    BToKstMuMuParticles.push_back(partFactory.particle(kaonTT,KaonMass_,chi,ndf,KaonMassErr_));
    BToKstMuMuParticles.push_back(partFactory.particle(pionTT,PionMass_,chi,ndf,PionMassErr_));
 
    RefCountedKinematicTree BToKstMuMuVertexFitTree = PartVtxFitter.fit(BToKstMuMuParticles);
    
    if ( !BToKstMuMuVertexFitTree->isValid()) return false;
    
    BToKstMuMuVertexFitTree->movePointerToTheTop();
    refitVertex = BToKstMuMuVertexFitTree->currentDecayVertex();
    refitBToKstMuMu = BToKstMuMuVertexFitTree->currentParticle();
    
    if ( !refitVertex->vertexIsValid()) return false;
    
    // extract the re-fitted tracks
    BToKstMuMuVertexFitTree->movePointerToTheTop();
    
    BToKstMuMuVertexFitTree->movePointerToTheFirstChild();
    refitMuon1 = BToKstMuMuVertexFitTree->currentParticle();
    
    BToKstMuMuVertexFitTree->movePointerToTheNextChild();
    refitMuon2 = BToKstMuMuVertexFitTree->currentParticle();
    
    BToKstMuMuVertexFitTree->movePointerToTheNextChild();
    refitKaon = BToKstMuMuVertexFitTree->currentParticle();
    
    BToKstMuMuVertexFitTree->movePointerToTheNextChild();
    refitPion = BToKstMuMuVertexFitTree->currentParticle();

    return true;



}



bool BToKstmumuProducer::BToKstJPsiMuMuVertexRefitting(const RefCountedKinematicParticle refitMuMu,
						       const RefCountedKinematicParticle refitKPi,
						       RefCountedKinematicVertex &refitVertex,
						       RefCountedKinematicParticle &refitBToKstJPsiMuMu,
						       RefCountedKinematicParticle &refitJPsi,
						       RefCountedKinematicParticle &refitKst){

  const reco::TransientTrack MuMuTT = refitMuMu->refittedTransientTrack();
  const reco::TransientTrack KPiTT = refitKPi->refittedTransientTrack();

  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter PartVtxFitter;

  std::vector<RefCountedKinematicParticle> BToKstMuMuParticles;
  double chi = 0.;
  double ndf = 0.;

  float MuMu_mass = refitMuMu->currentState().mass();
  float MuMu_mass_err = sqrt(refitMuMu->currentState().kinematicParametersError().matrix()(6,6));
  if(JPsiMassConstraint_ > 0){
    MuMu_mass = JPsiMassConstraint_;
    MuMu_mass_err = JPsiMassErr_;
  }

  float Kst_mass = refitKPi->currentState().mass();
  float Kst_mass_err = sqrt(refitKPi->currentState().kinematicParametersError().matrix()(6,6));
  if(KstMassConstraint_ > 0){
    Kst_mass = KstMassConstraint_;
    Kst_mass_err = KstMassErr_;
  }

  BToKstMuMuParticles.push_back(partFactory.particle(MuMuTT,MuMu_mass,chi,ndf,MuMu_mass_err));
  BToKstMuMuParticles.push_back(partFactory.particle(KPiTT,Kst_mass,chi,ndf,Kst_mass_err));

  RefCountedKinematicTree BToKstMuMuVertexFitTree = PartVtxFitter.fit(BToKstMuMuParticles);

  if ( !BToKstMuMuVertexFitTree->isValid()) return false;

  BToKstMuMuVertexFitTree->movePointerToTheTop();
  refitVertex = BToKstMuMuVertexFitTree->currentDecayVertex();
  refitBToKstJPsiMuMu = BToKstMuMuVertexFitTree->currentParticle();

  if ( !refitVertex->vertexIsValid()) return false;

  // extract the re-fitted tracks
  BToKstMuMuVertexFitTree->movePointerToTheTop();

  BToKstMuMuVertexFitTree->movePointerToTheFirstChild();
  refitJPsi = BToKstMuMuVertexFitTree->currentParticle();

  BToKstMuMuVertexFitTree->movePointerToTheNextChild();
  refitKst = BToKstMuMuVertexFitTree->currentParticle();

  return true;



}




pair<double,double> BToKstmumuProducer::computeLS(RefCountedKinematicVertex refitVertex,
						  reco::BeamSpot beamSpot){
    
  TVector v(2);
  v[0] = refitVertex->position().x()-beamSpot.position().x();
  v[1] = refitVertex->position().y()-beamSpot.position().y();

  TMatrix errVtx(2,2);
  errVtx(0,0) = refitVertex->error().cxx();
  errVtx(0,1) = refitVertex->error().matrix()(0,1);
  errVtx(1,0) = errVtx(0,1);
  errVtx(1,1) = refitVertex->error().cyy();

  TMatrix errBS(2,2);
  errBS(0,0) = beamSpot.covariance()(0,0);
  errBS(0,1) = beamSpot.covariance()(0,1);
  errBS(1,0) = beamSpot.covariance()(1,0);
  errBS(1,1) = beamSpot.covariance()(1,1);
    
  double LSBS = sqrt(v.Norm2Sqr());
  double LSBSErr = sqrt( v*(errVtx*v) + v*(errBS*v) ) / LSBS;
    
  pair<double,double> LS = make_pair(LSBS,LSBSErr);
    
  return LS;
}



double BToKstmumuProducer::computeCosAlpha(RefCountedKinematicParticle refitBToKstMuMu,
					   RefCountedKinematicVertex refitVertex,
					   reco::BeamSpot beamSpot){
    
  TVector v(2);
  v[0] = refitVertex->position().x()-beamSpot.position().x();
  v[1] = refitVertex->position().y()-beamSpot.position().y();

  TVector w(2);
  w[0] = refitBToKstMuMu->currentState().globalMomentum().x();
  w[1] = refitBToKstMuMu->currentState().globalMomentum().y();
    
  double cosAlpha = v*w/sqrt(v.Norm2Sqr()*w.Norm2Sqr());
  return cosAlpha;
}






pair<double,double> BToKstmumuProducer::computeDCA(const pat::PackedCandidate &kaon,
						   edm::ESHandle<MagneticField> bFieldHandle,
						   reco::BeamSpot beamSpot){
  
  const reco::TransientTrack trackTT((*(kaon.bestTrack())), &(*bFieldHandle));

  TrajectoryStateClosestToPoint theDCAXBS = trackTT.trajectoryStateClosestToPoint( GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()) );  
  
  double DCABS = theDCAXBS.perigeeParameters().transverseImpactParameter();
  double DCABSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
    
  pair<double,double> DCA = make_pair(DCABS,DCABSErr);
    
  return DCA;
}




DEFINE_FWK_MODULE(BToKstmumuProducer);
