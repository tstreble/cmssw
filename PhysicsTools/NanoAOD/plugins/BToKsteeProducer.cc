#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

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

class BToKsteeProducer : public edm::EDProducer {
    
public:
    
    explicit BToKsteeProducer(const edm::ParameterSet &iConfig);
    
    ~BToKsteeProducer() override {};
    
    
private:
    
    virtual void produce(edm::Event&, const edm::EventSetup&);
    
    bool EEVertexRefitting(const pat::Electron & ele1,
			   const pat::Electron & ele2,
			   edm::ESHandle<TransientTrackBuilder> theTTBuilder,
			   RefCountedKinematicVertex &refitVertex,
			   RefCountedKinematicParticle &refitEE,
			   RefCountedKinematicParticle &refitEle1,
			   RefCountedKinematicParticle &refitEle2);

    bool KstVertexRefitting(const pat::PackedCandidate &kaon,
			    const pat::PackedCandidate &pion,
			    edm::ESHandle<TransientTrackBuilder> theTTBuilder,
			    RefCountedKinematicVertex &refitVertex,
			    RefCountedKinematicParticle &refitKst,
			    RefCountedKinematicParticle &refitKaon,
			    RefCountedKinematicParticle &refitPion);

    bool BToKstEEVertexRefitting(const pat::Electron &ele1,
                                 const pat::Electron &ele2,
                                 const RefCountedKinematicParticle refitKPi,
                                 edm::ESHandle<TransientTrackBuilder> theTTBuilder,
                                 RefCountedKinematicVertex &refitVertex,
                                 RefCountedKinematicParticle &refitBToKstEE,
                                 RefCountedKinematicParticle &refitEle1,
                                 RefCountedKinematicParticle &refitEle2,
                                 RefCountedKinematicParticle &refitKst);

    bool BToKPiEEVertexRefitting(const pat::Electron &ele1,
                                 const pat::Electron &ele2,
                                 const pat::PackedCandidate &kaon,
                                 const pat::PackedCandidate &pion,
                                 edm::ESHandle<TransientTrackBuilder> theTTBuilder,
                                 RefCountedKinematicVertex &refitVertex,
                                 RefCountedKinematicParticle &refitBToKPiEE,
                                 RefCountedKinematicParticle &refitEle1,
                                 RefCountedKinematicParticle &refitEle2,
                                 RefCountedKinematicParticle &refitKaon,
                                 RefCountedKinematicParticle &refitPion);

    bool BToKstJPsiEEVertexRefitting(const RefCountedKinematicParticle refitEE,
                                     const RefCountedKinematicParticle refitKPi,
                                     RefCountedKinematicVertex &refitVertex,
                                     RefCountedKinematicParticle &refitBToKstJPsiEE,
                                     RefCountedKinematicParticle &refitJPsi,
                                     RefCountedKinematicParticle &refitKst);
    
    pair<double,double> computeLS(RefCountedKinematicVertex refitVertex,
                                  reco::BeamSpot beamSpot);
    
    double computeCosAlpha(RefCountedKinematicParticle refitBToKstEE,
                           RefCountedKinematicVertex vertexFitTree,
                           reco::BeamSpot beamSpot);

    pair<double,double> computeDCA(const pat::PackedCandidate &kaon,
				   edm::ESHandle<MagneticField> bFieldHandle,
				   reco::BeamSpot beamSpot);
    
    // ----------member data ---------------------------
    
    edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;
    edm::EDGetTokenT<std::vector<pat::Electron>> electronSrc_;
    edm::EDGetTokenT<edm::View<pat::PackedCandidate>> PFCandSrc_;
    
    double ptMinEle_;
    double etaMaxEle_;
    double ptMinKaon_;
    double etaMaxKaon_;
    double DCASigMinKaon_;
    double ptMinPion_;
    double etaMaxPion_;
    double DCASigMinPion_;
    bool diEleCharge_;
    bool KstCharge_;
    double JPsiMassConstraint_;
    double KstMassConstraint_;
    bool save2TrkRefit_;
    bool save4TrkRefit_;
    
    float ElectronMass_ = 0.5109989e-3;
    float ElectronMassErr_ = 3.1*1e-12;
    float KaonMass_ = 0.493677;
    float KaonMassErr_ = 1.6e-5;
    float PionMass_ = 0.139570;
    float PionMassErr_ = 3.5e-7;
    //float JPsiMass_ = 3.096916;  //Configurable parameter
    float JPsiMassErr_ = 0.011;
    //float KstMass_ = 0.89176;  //Configurable parameter
    float KstMassErr_ = 0.25e-3;

};



BToKsteeProducer::BToKsteeProducer(const edm::ParameterSet &iConfig):
beamSpotSrc_( consumes<reco::BeamSpot> ( iConfig.getParameter<edm::InputTag>( "beamSpot" ) ) ),
electronSrc_( consumes<std::vector<pat::Electron>> ( iConfig.getParameter<edm::InputTag>( "electronCollection" ) ) ),
PFCandSrc_( consumes<edm::View<pat::PackedCandidate>> ( iConfig.getParameter<edm::InputTag>( "PFCandCollection" ) ) ),
ptMinEle_( iConfig.getParameter<double>( "ElectronMinPt" ) ),
etaMaxEle_( iConfig.getParameter<double>( "ElectronMaxEta" ) ),
ptMinKaon_( iConfig.getParameter<double>( "KaonMinPt" ) ),
etaMaxKaon_( iConfig.getParameter<double>( "KaonMaxEta" ) ),
DCASigMinKaon_( iConfig.getParameter<double>( "KaonMinDCASig" ) ),
ptMinPion_( iConfig.getParameter<double>( "PionMinPt" ) ),
etaMaxPion_( iConfig.getParameter<double>( "PionMaxEta" ) ),
DCASigMinPion_( iConfig.getParameter<double>( "PionMinDCASig" ) ),
diEleCharge_( iConfig.getParameter<bool>( "DiElectronChargeCheck" ) ),
KstCharge_( iConfig.getParameter<bool>( "KstarChargeCheck" ) ),
JPsiMassConstraint_( iConfig.getParameter<double>( "JPsiMassConstraint" ) ),
KstMassConstraint_( iConfig.getParameter<double>( "KstMassConstraint" ) ),
save2TrkRefit_( iConfig.getParameter<bool>( "save2TrackRefit" ) ),
save4TrkRefit_( iConfig.getParameter<bool>( "save4TrackRefit" ) )
{
    produces<pat::CompositeCandidateCollection>();
}


void BToKsteeProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    
    edm::ESHandle<MagneticField> bFieldHandle;
    edm::ESHandle<TransientTrackBuilder> theTTBuilder;
 
    iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);

    edm::Handle<reco::BeamSpot> beamSpotHandle;
    iEvent.getByToken(beamSpotSrc_, beamSpotHandle);
    
    if ( ! beamSpotHandle.isValid() ) {
        edm::LogError("BToKsteeProducer") << "No beam spot available from EventSetup" ;
    }
    
    reco::BeamSpot beamSpot = *beamSpotHandle;

    edm::Handle<std::vector<pat::Electron>> electronHandle;
    edm::Handle<edm::View<pat::PackedCandidate>> pfCandHandle;

    iEvent.getByToken(electronSrc_, electronHandle);
    iEvent.getByToken(PFCandSrc_, pfCandHandle);

    unsigned int electronNumber = electronHandle->size();
    unsigned int pfCandNumber = pfCandHandle->size();

    // Output collection
    std::unique_ptr<pat::CompositeCandidateCollection> result( new pat::CompositeCandidateCollection );

    if(electronNumber>1){

        // loop on all the eeKPi quadruplets
        for (unsigned int i = 0; i < electronNumber; ++i) {

            const pat::Electron & ele1 = (*electronHandle)[i];

	    //could implement ele ID criteria here
            if(ele1.pt()<ptMinEle_ || abs(ele1.eta())>etaMaxEle_) continue;
            
            for (unsigned int j = 0; j < electronNumber; ++j) {

                if(i==j) continue;

                const pat::Electron & ele2 = (*electronHandle)[j];

		if(ele1.pt()<ele2.pt()) continue; //Electron 1 is always saved as the leading one
		//could implement ele ID criteria here
                if(ele2.pt()<ptMinEle_ || abs(ele2.eta())>etaMaxEle_) continue;
                if(diEleCharge_ && ele1.charge()*ele2.charge()>0) continue;

                bool passedDiEle = false;

		double EELSBS = -1.;
                double EELSBSErr = -1.;
                double EEVtx_CL = -1.;
                double EE_mass_err = -1.;

                RefCountedKinematicParticle refitEE;
                math::XYZVector refitEEV3D;

                if(save2TrkRefit_){

                  RefCountedKinematicVertex refitVertexEE;
                  RefCountedKinematicParticle refitEle1_EE;
                  RefCountedKinematicParticle refitEle2_EE;

                  passedDiEle = EEVertexRefitting(ele1, ele2,
						  theTTBuilder,
						  refitVertexEE,
						  refitEE,
						  refitEle1_EE,
						  refitEle2_EE);

		  if(passedDiEle){

		    math::XYZVector refitEle1V3D_EE = refitEle1_EE->refittedTransientTrack().track().momentum();
		    math::XYZVector refitEle2V3D_EE = refitEle2_EE->refittedTransientTrack().track().momentum();
		    refitEEV3D = refitEle1V3D_EE + refitEle2V3D_EE;

		    pair<double,double> EELS = computeLS(refitVertexEE,beamSpot);
		    EELSBS = EELS.first;
		    EELSBSErr = EELS.second;

		    EEVtx_CL = TMath::Prob((double)refitVertexEE->chiSquared(),
					   int(rint(refitVertexEE->degreesOfFreedom())));

		    EE_mass_err = sqrt(refitEE->currentState().kinematicParametersError().matrix()(6,6));

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

		      
		      RefCountedKinematicVertex refitVertexBToKstEE;
		      RefCountedKinematicParticle refitBToKstEE;
		      RefCountedKinematicParticle refitEle1;
		      RefCountedKinematicParticle refitEle2;
		      RefCountedKinematicParticle refitKst_BToKstEE;

		      passed = BToKstEEVertexRefitting(ele1, ele2, refitKst,
						       theTTBuilder,
						       refitVertexBToKstEE,
						       refitBToKstEE,
						       refitEle1,
						       refitEle2,
						       refitKst_BToKstEE);
                    
		      if (!passed) continue;
                    
		      pair<double,double> BToKstEELS = computeLS(refitVertexBToKstEE,beamSpot);
		      double LSBS = BToKstEELS.first;
		      double LSBSErr = BToKstEELS.second;
                    
		      double BToKstEEVtx_CL = TMath::Prob((double)refitVertexBToKstEE->chiSquared(),
							  int(rint(refitVertexBToKstEE->degreesOfFreedom())));

                    
		      double cosAlpha = computeCosAlpha(refitBToKstEE,refitVertexBToKstEE,beamSpot);
                    
		      double mass_err = sqrt(refitBToKstEE->currentState().kinematicParametersError().matrix()(6,6));

		      math::XYZVector refitEle1V3D = refitEle1->refittedTransientTrack().track().momentum();
		      math::XYZVector refitEle2V3D = refitEle2->refittedTransientTrack().track().momentum();
		      math::XYZVector refitKst_BToKstEE_V3D = refitKst_BToKstEE->refittedTransientTrack().track().momentum();
		      math::XYZVector refitBToKstEEV3D = refitEle1V3D + refitEle2V3D + refitKst_BToKstEE_V3D;

      
		      pat::CompositeCandidate BToKstEECand;
		      BToKstEECand.addDaughter( ele1, "ele1");
		      BToKstEECand.addDaughter( ele2, "ele2");
		      BToKstEECand.addDaughter( kaon, "kaon");
		      BToKstEECand.addDaughter( pion, "pion");

		      BToKstEECand.addUserInt("ele1_index", i);
		      BToKstEECand.addUserInt("ele2_index", j);
		      BToKstEECand.addUserInt("kaon_index", k);
		      BToKstEECand.addUserInt("pion_index", l);
		      
		      BToKstEECand.addUserFloat("ele1_pt",     sqrt(refitEle1V3D.perp2()));
		      BToKstEECand.addUserFloat("ele1_eta",    refitEle1V3D.eta());
		      BToKstEECand.addUserFloat("ele1_phi",    refitEle1V3D.phi());
		      BToKstEECand.addUserFloat("ele1_charge", refitEle1->currentState().particleCharge());

		      BToKstEECand.addUserFloat("ele2_pt",     sqrt(refitEle2V3D.perp2()));
		      BToKstEECand.addUserFloat("ele2_eta",    refitEle2V3D.eta());
		      BToKstEECand.addUserFloat("ele2_phi",    refitEle2V3D.phi());
		      BToKstEECand.addUserFloat("ele2_charge", refitEle2->currentState().particleCharge());		    

		      TLorentzVector ele1cand;
		      ele1cand.SetPtEtaPhiM(sqrt(refitEle1V3D.perp2()), refitEle1V3D.eta(), refitEle1V3D.phi(), ElectronMass_);
		      TLorentzVector ele2cand;
		      ele2cand.SetPtEtaPhiM(sqrt(refitEle2V3D.perp2()), refitEle2V3D.eta(), refitEle2V3D.phi(), ElectronMass_);
		      BToKstEECand.addUserFloat("eeKPiFit_e_mass", (ele1cand+ele2cand).Mag());

		      BToKstEECand.addUserFloat("kaon_pt",     sqrt(refitKaonV3D_Kst.perp2()));
		      BToKstEECand.addUserFloat("kaon_eta",    refitKaonV3D_Kst.eta());
		      BToKstEECand.addUserFloat("kaon_phi",    refitKaonV3D_Kst.phi());
		      BToKstEECand.addUserFloat("kaon_charge", refitKaon_Kst->currentState().particleCharge());
		      BToKstEECand.addUserFloat("kaon_DCASig", DCABS_kaon/DCABSErr_kaon);
		      
		      BToKstEECand.addUserFloat("pion_pt",     sqrt(refitPionV3D_Kst.perp2()));
		      BToKstEECand.addUserFloat("pion_eta",    refitPionV3D_Kst.eta());
		      BToKstEECand.addUserFloat("pion_phi",    refitPionV3D_Kst.phi());
		      BToKstEECand.addUserFloat("pion_charge", refitPion_Kst->currentState().particleCharge());
		      BToKstEECand.addUserFloat("pion_DCASig", DCABS_pion/DCABSErr_pion);

		      BToKstEECand.addUserFloat("Kst_pt", sqrt(refitKstV3D.perp2()));
		      BToKstEECand.addUserFloat("Kst_eta", refitKstV3D.eta());
		      BToKstEECand.addUserFloat("Kst_phi", refitKstV3D.phi());
		      BToKstEECand.addUserFloat("Kst_mass", refitKst->currentState().mass());
		      BToKstEECand.addUserFloat("Kst_mass_err",  Kst_mass_err);
		      BToKstEECand.addUserFloat("Kst_Lxy", (float) KstLSBS/KstLSBSErr);
		      BToKstEECand.addUserFloat("Kst_ctxy", (float) KstLSBS/sqrt(refitKstV3D.perp2()));
		      BToKstEECand.addUserFloat("Kst_CL_vtx", (float) KstVtx_CL);

		      BToKstEECand.addUserFloat("pt",     sqrt(refitBToKstEEV3D.perp2()));
		      BToKstEECand.addUserFloat("eta",    refitBToKstEEV3D.eta());
		      BToKstEECand.addUserFloat("phi",    refitBToKstEEV3D.phi());
		      BToKstEECand.addUserFloat("mass",   refitBToKstEE->currentState().mass());
		      BToKstEECand.addUserFloat("mass_err", mass_err);
		      BToKstEECand.addUserFloat("Lxy", (float) LSBS/LSBSErr);
		      BToKstEECand.addUserFloat("ctxy", (float) LSBS/sqrt(refitBToKstEEV3D.perp2()));
		      BToKstEECand.addUserFloat("CL_vtx", (float) BToKstEEVtx_CL);
		      BToKstEECand.addUserFloat("cosAlpha", (float) cosAlpha);

                    
		      BToKstEECand.addUserInt("eeRefit", (int)passedDiEle);
		      BToKstEECand.addUserFloat("ee_pt", (passedDiEle)? sqrt(refitEEV3D.perp2()) : -1.);
		      BToKstEECand.addUserFloat("ee_eta", (passedDiEle)? refitEEV3D.eta() : -9.);
		      BToKstEECand.addUserFloat("ee_phi", (passedDiEle)? refitEEV3D.phi() : -9.);
		      BToKstEECand.addUserFloat("ee_mass", (passedDiEle)? refitEE->currentState().mass() : -1.);
		      BToKstEECand.addUserFloat("ee_mass_err", (passedDiEle)?  EE_mass_err : -1.);
		      BToKstEECand.addUserFloat("ee_Lxy", (passedDiEle)? (float) EELSBS/EELSBSErr : -1.);
		      BToKstEECand.addUserFloat("ee_ctxy", (passedDiEle)? (float) EELSBS/sqrt(refitEEV3D.perp2()) : -1.);
		      BToKstEECand.addUserFloat("ee_CL_vtx", (passedDiEle)? (float) EEVtx_CL : -1.);


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

		      //if(save2TrkRefit_ && passedDiEle){
		      if(0){

			RefCountedKinematicVertex refitVertexBToKstJPsiEE;
			RefCountedKinematicParticle refitBToKstJPsiEE;
			RefCountedKinematicParticle refitJPsiEE;
			RefCountedKinematicParticle refitKst_KstJPsi;

			bool passed_2trk = BToKstJPsiEEVertexRefitting(refitEE, refitKst,
								       refitVertexBToKstJPsiEE,
								       refitBToKstJPsiEE,
								       refitJPsiEE,
								       refitKst_KstJPsi);
			if(passed_2trk){

			  math::XYZVector refitJPsiEEV3D = refitJPsiEE->refittedTransientTrack().track().momentum();
			  math::XYZVector refitKstV3D_KJPsi = refitKst_KstJPsi->refittedTransientTrack().track().momentum();
			  math::XYZVector refitBToKstJPsiEEV3D = refitJPsiEEV3D + refitKstV3D_KJPsi;
			  
			  pt_2trk = sqrt(refitBToKstJPsiEEV3D.perp2());
			  eta_2trk = refitBToKstJPsiEEV3D.eta();
			  phi_2trk = refitBToKstJPsiEEV3D.phi();
			  mass_2trk = refitBToKstJPsiEE->currentState().mass();

			  mass_err_2trk = sqrt(refitBToKstJPsiEE->currentState().kinematicParametersError().matrix()(6,6));

			  pair<double,double> BToKstJPsiEELS = computeLS(refitVertexBToKstJPsiEE,beamSpot);
			  double LSBS_2trk = BToKstJPsiEELS.first;
			  double LSBSErr_2trk = BToKstJPsiEELS.second;
			  Lxy_2trk = LSBS_2trk/LSBSErr_2trk;
			  ctxy_2trk = LSBS_2trk/pt_2trk;
			  CL_vtx_2trk = TMath::Prob((double)refitVertexBToKstJPsiEE->chiSquared(),
						    int(rint(refitVertexBToKstJPsiEE->degreesOfFreedom())));
			  cosAlpha_2trk = computeCosAlpha(refitBToKstJPsiEE,refitVertexBToKstJPsiEE,beamSpot);

			}

		      }

		      BToKstEECand.addUserInt("2trkRefit", (int)passed_2trk);
		      BToKstEECand.addUserFloat("pt_2trk", pt_2trk);
		      BToKstEECand.addUserFloat("eta_2trk", eta_2trk);
		      BToKstEECand.addUserFloat("phi_2trk", phi_2trk);
		      BToKstEECand.addUserFloat("mass_2trk", mass_2trk);
		      BToKstEECand.addUserFloat("mass_err_2trk", mass_err_2trk);
		      BToKstEECand.addUserFloat("Lxy_2trk", Lxy_2trk);
		      BToKstEECand.addUserFloat("ctxy_2trk", ctxy_2trk);
		      BToKstEECand.addUserFloat("CL_vtx_2trk", CL_vtx_2trk);
		      BToKstEECand.addUserFloat("cosAlpha_2trk", cosAlpha_2trk);



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

			RefCountedKinematicVertex refitVertexBToKPiEE;
			RefCountedKinematicParticle refitBToKPiEE;
			RefCountedKinematicParticle refitEle1_BToKPiEE;
			RefCountedKinematicParticle refitEle2_BToKPiEE;
			RefCountedKinematicParticle refitKaon_BToKPiEE;
			RefCountedKinematicParticle refitPion_BToKPiEE;

			passed_4trk = BToKPiEEVertexRefitting(ele1, ele2, kaon, pion,
							      theTTBuilder,
							      refitVertexBToKPiEE,
							      refitBToKPiEE,
							      refitEle1_BToKPiEE,
							      refitEle2_BToKPiEE,
							      refitKaon_BToKPiEE,
							      refitPion_BToKPiEE);


			if(passed_4trk){

			  math::XYZVector refitEle1_BToKPiEE_V3D = refitEle1_BToKPiEE->refittedTransientTrack().track().momentum();
			  math::XYZVector refitEle2_BToKPiEE_V3D = refitEle2_BToKPiEE->refittedTransientTrack().track().momentum();
			  math::XYZVector refitKaon_BToKPiEE_V3D = refitKaon_BToKPiEE->refittedTransientTrack().track().momentum();
			  math::XYZVector refitPion_BToKPiEE_V3D = refitPion_BToKPiEE->refittedTransientTrack().track().momentum();

			  math::XYZVector refitBToKPiEEV3D = refitEle1_BToKPiEE_V3D + refitEle2_BToKPiEE_V3D + refitKaon_BToKPiEE_V3D + refitPion_BToKPiEE_V3D;

			  pt_4trk = sqrt(refitBToKPiEEV3D.perp2());
			  eta_4trk = refitBToKPiEEV3D.eta();
			  phi_4trk = refitBToKPiEEV3D.phi();
			  mass_4trk = refitBToKPiEE->currentState().mass();

			  mass_err_4trk = sqrt(refitBToKPiEE->currentState().kinematicParametersError().matrix()(6,6));

			  pair<double,double> BToKPiEELS = computeLS(refitVertexBToKPiEE,beamSpot);
			  double LSBS_4trk = BToKPiEELS.first;
			  double LSBSErr_4trk = BToKPiEELS.second;
			  Lxy_4trk = LSBS_4trk/LSBSErr_4trk;
			  ctxy_4trk = LSBS_4trk/pt_4trk;
			  CL_vtx_4trk = TMath::Prob((double)refitVertexBToKPiEE->chiSquared(),
						    int(rint(refitVertexBToKPiEE->degreesOfFreedom())));
			  cosAlpha_4trk = computeCosAlpha(refitBToKPiEE,refitVertexBToKPiEE,beamSpot);

			}

		      }

		      BToKstEECand.addUserInt("4trkRefit", (int)passed_4trk);
		      BToKstEECand.addUserFloat("pt_4trk", pt_4trk);
		      BToKstEECand.addUserFloat("eta_4trk", eta_4trk);
		      BToKstEECand.addUserFloat("phi_4trk", phi_4trk);
		      BToKstEECand.addUserFloat("mass_4trk", mass_4trk);
		      BToKstEECand.addUserFloat("mass_err_4trk", mass_err_4trk);
		      BToKstEECand.addUserFloat("Lxy_4trk", Lxy_4trk);
		      BToKstEECand.addUserFloat("ctxy_4trk", ctxy_4trk);
		      BToKstEECand.addUserFloat("CL_vtx_4trk", CL_vtx_4trk);
		      BToKstEECand.addUserFloat("cosAlpha_4trk", cosAlpha_4trk);


		      result->push_back(BToKstEECand);
                    
		    }
                    
                }
                
            }
            
        }
        
    }
    
    iEvent.put(std::move(result));
    
}



bool BToKsteeProducer::EEVertexRefitting(const pat::Electron & ele1,
				       const pat::Electron & ele2,
				       edm::ESHandle<TransientTrackBuilder> theTTBuilder,
				       RefCountedKinematicVertex &refitVertex,
				       RefCountedKinematicParticle &refitEE,
				       RefCountedKinematicParticle &refitEle1,
				       RefCountedKinematicParticle &refitEle2){
    
    const reco::TransientTrack ele1TT = theTTBuilder->build(ele1.gsfTrack()); //closestCtfTrackRef ?
    const reco::TransientTrack ele2TT = theTTBuilder->build(ele2.gsfTrack());
    
    KinematicParticleFactoryFromTransientTrack partFactory;
    KinematicParticleVertexFitter PartVtxFitter;
    
    std::vector<RefCountedKinematicParticle> eleParticles;
    double chi = 0.;
    double ndf = 0.;
    eleParticles.push_back(partFactory.particle(ele1TT,ElectronMass_,chi,ndf,ElectronMassErr_));
    eleParticles.push_back(partFactory.particle(ele2TT,ElectronMass_,chi,ndf,ElectronMassErr_));
    RefCountedKinematicTree eeVertexFitTree = PartVtxFitter.fit(eleParticles);
    
    if ( !eeVertexFitTree->isValid()) return false;
    
    eeVertexFitTree->movePointerToTheTop();
    refitVertex = eeVertexFitTree->currentDecayVertex();
    refitEE = eeVertexFitTree->currentParticle();
    
    if ( !refitVertex->vertexIsValid()) return false;

    // extract the re-fitted tracks
    eeVertexFitTree->movePointerToTheTop();
    
    eeVertexFitTree->movePointerToTheFirstChild();
    refitEle1 = eeVertexFitTree->currentParticle();
    
    eeVertexFitTree->movePointerToTheNextChild();
    refitEle2 = eeVertexFitTree->currentParticle();
    
    return true;
    
}




bool BToKsteeProducer::KstVertexRefitting(const pat::PackedCandidate &kaon,
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





bool BToKsteeProducer::BToKstEEVertexRefitting(const pat::Electron &ele1,
					       const pat::Electron &ele2,
					       const RefCountedKinematicParticle refitKPi,
					       edm::ESHandle<TransientTrackBuilder> theTTBuilder,					   
					       RefCountedKinematicVertex &refitVertex,
					       RefCountedKinematicParticle &refitBToKstEE,
					       RefCountedKinematicParticle &refitEle1,
					       RefCountedKinematicParticle &refitEle2,
					       RefCountedKinematicParticle &refitKst){

    const reco::TransientTrack ele1TT = theTTBuilder->build(ele1.gsfTrack());
    const reco::TransientTrack ele2TT = theTTBuilder->build(ele2.gsfTrack());
    const reco::TransientTrack KPiTT = refitKPi->refittedTransientTrack();

    KinematicParticleFactoryFromTransientTrack partFactory;
    KinematicParticleVertexFitter PartVtxFitter;

    float Kst_mass = refitKPi->currentState().mass();
    float Kst_mass_err = sqrt(refitKPi->currentState().kinematicParametersError().matrix()(6,6));
    if(KstMassConstraint_ > 0){
      Kst_mass = KstMassConstraint_;
      Kst_mass_err = KstMassErr_;
    }

    std::vector<RefCountedKinematicParticle> BToKstEEParticles;
    double chi = 0.;
    double ndf = 0.;
    BToKstEEParticles.push_back(partFactory.particle(ele1TT,ElectronMass_,chi,ndf,ElectronMassErr_));
    BToKstEEParticles.push_back(partFactory.particle(ele2TT,ElectronMass_,chi,ndf,ElectronMassErr_));
    BToKstEEParticles.push_back(partFactory.particle(KPiTT,Kst_mass,chi,ndf,Kst_mass_err));

    RefCountedKinematicTree BToKstEEVertexFitTree = PartVtxFitter.fit(BToKstEEParticles);

    if ( !BToKstEEVertexFitTree->isValid()) return false;

    BToKstEEVertexFitTree->movePointerToTheTop();
    refitVertex = BToKstEEVertexFitTree->currentDecayVertex();
    refitBToKstEE = BToKstEEVertexFitTree->currentParticle();

    if ( !refitVertex->vertexIsValid()) return false;

    // extract the re-fitted tracks
    BToKstEEVertexFitTree->movePointerToTheTop();

    BToKstEEVertexFitTree->movePointerToTheFirstChild();
    refitEle1 = BToKstEEVertexFitTree->currentParticle();

    BToKstEEVertexFitTree->movePointerToTheNextChild();
    refitEle2 = BToKstEEVertexFitTree->currentParticle();

    BToKstEEVertexFitTree->movePointerToTheNextChild();
    refitKst = BToKstEEVertexFitTree->currentParticle();

    return true;

}





bool BToKsteeProducer::BToKPiEEVertexRefitting(const pat::Electron &ele1,
					       const pat::Electron &ele2,
					       const pat::PackedCandidate &kaon,
					       const pat::PackedCandidate &pion,
					       edm::ESHandle<TransientTrackBuilder> theTTBuilder,
					       RefCountedKinematicVertex &refitVertex,
					       RefCountedKinematicParticle &refitBToKPiEE,
					       RefCountedKinematicParticle &refitEle1,
					       RefCountedKinematicParticle &refitEle2,
					       RefCountedKinematicParticle &refitKaon,
					       RefCountedKinematicParticle &refitPion){

    const reco::TransientTrack ele1TT = theTTBuilder->build(ele1.gsfTrack());
    const reco::TransientTrack ele2TT = theTTBuilder->build(ele2.gsfTrack());
    const reco::TransientTrack kaonTT = theTTBuilder->build(kaon.bestTrack());
    const reco::TransientTrack pionTT = theTTBuilder->build(pion.bestTrack());

    KinematicParticleFactoryFromTransientTrack partFactory;
    KinematicParticleVertexFitter PartVtxFitter;

    std::vector<RefCountedKinematicParticle> BToKPiEEParticles;
    double chi = 0.;
    double ndf = 0.;
    BToKPiEEParticles.push_back(partFactory.particle(ele1TT,ElectronMass_,chi,ndf,ElectronMassErr_));
    BToKPiEEParticles.push_back(partFactory.particle(ele2TT,ElectronMass_,chi,ndf,ElectronMassErr_));
    BToKPiEEParticles.push_back(partFactory.particle(kaonTT,KaonMass_,chi,ndf,KaonMassErr_));
    BToKPiEEParticles.push_back(partFactory.particle(pionTT,PionMass_,chi,ndf,PionMassErr_));

    RefCountedKinematicTree BToKPiEEVertexFitTree = PartVtxFitter.fit(BToKPiEEParticles);
    
    if ( !BToKPiEEVertexFitTree->isValid()) return false;

    BToKPiEEVertexFitTree->movePointerToTheTop();
    refitVertex = BToKPiEEVertexFitTree->currentDecayVertex();
    refitBToKPiEE = BToKPiEEVertexFitTree->currentParticle();

    if ( !refitVertex->vertexIsValid()) return false;

    // extract the re-fitted tracks
    BToKPiEEVertexFitTree->movePointerToTheTop();

    BToKPiEEVertexFitTree->movePointerToTheFirstChild();
    refitEle1 = BToKPiEEVertexFitTree->currentParticle();

    BToKPiEEVertexFitTree->movePointerToTheNextChild();
    refitEle2 = BToKPiEEVertexFitTree->currentParticle();

    BToKPiEEVertexFitTree->movePointerToTheNextChild();
    refitKaon = BToKPiEEVertexFitTree->currentParticle();

    BToKPiEEVertexFitTree->movePointerToTheNextChild();
    refitPion = BToKPiEEVertexFitTree->currentParticle();

    return true;



}



bool BToKsteeProducer::BToKstJPsiEEVertexRefitting(const RefCountedKinematicParticle refitEE,
						   const RefCountedKinematicParticle refitKPi,						  
						   RefCountedKinematicVertex &refitVertex,
						   RefCountedKinematicParticle &refitBToKstJPsiEE,
						   RefCountedKinematicParticle &refitJPsi,
						   RefCountedKinematicParticle &refitKst){

  const reco::TransientTrack EETT = refitEE->refittedTransientTrack();
  const reco::TransientTrack KPiTT = refitKPi->refittedTransientTrack();

  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter PartVtxFitter;

  std::vector<RefCountedKinematicParticle> BToKstEEParticles;
  double chi = 0.;
  double ndf = 0.;

  float EE_mass = refitEE->currentState().mass();
  float EE_mass_err = sqrt(refitEE->currentState().kinematicParametersError().matrix()(6,6));
  if(JPsiMassConstraint_ > 0){
    EE_mass = JPsiMassConstraint_;
    EE_mass_err = JPsiMassErr_;
  }

  float Kst_mass = refitKPi->currentState().mass();
  float Kst_mass_err = sqrt(refitKPi->currentState().kinematicParametersError().matrix()(6,6));
  if(KstMassConstraint_ > 0){
    Kst_mass = KstMassConstraint_;
    Kst_mass_err = KstMassErr_;
  }

  BToKstEEParticles.push_back(partFactory.particle(EETT,EE_mass,chi,ndf,EE_mass_err));
  BToKstEEParticles.push_back(partFactory.particle(KPiTT,Kst_mass,chi,ndf,Kst_mass_err));

  RefCountedKinematicTree BToKstEEVertexFitTree = PartVtxFitter.fit(BToKstEEParticles);

  if ( !BToKstEEVertexFitTree->isValid()) return false;

  BToKstEEVertexFitTree->movePointerToTheTop();
  refitVertex = BToKstEEVertexFitTree->currentDecayVertex();
  refitBToKstJPsiEE = BToKstEEVertexFitTree->currentParticle();

  if ( !refitVertex->vertexIsValid()) return false;

  // extract the re-fitted tracks
  BToKstEEVertexFitTree->movePointerToTheTop();

  BToKstEEVertexFitTree->movePointerToTheFirstChild();
  refitJPsi = BToKstEEVertexFitTree->currentParticle();

  BToKstEEVertexFitTree->movePointerToTheNextChild();
  refitKst = BToKstEEVertexFitTree->currentParticle();

  return true;



}




pair<double,double> BToKsteeProducer::computeLS(RefCountedKinematicVertex refitVertex,
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



double BToKsteeProducer::computeCosAlpha(RefCountedKinematicParticle refitBToKstEE,
				       RefCountedKinematicVertex refitVertex,
				       reco::BeamSpot beamSpot){
    
  TVector v(2);
  v[0] = refitVertex->position().x()-beamSpot.position().x();
  v[1] = refitVertex->position().y()-beamSpot.position().y();

  TVector w(2);
  w[0] = refitBToKstEE->currentState().globalMomentum().x();
  w[1] = refitBToKstEE->currentState().globalMomentum().y();
    
  double cosAlpha = v*w/sqrt(v.Norm2Sqr()*w.Norm2Sqr());
  return cosAlpha;
}






pair<double,double> BToKsteeProducer::computeDCA(const pat::PackedCandidate &kaon,
					       edm::ESHandle<MagneticField> bFieldHandle,
					       reco::BeamSpot beamSpot){
  
  const reco::TransientTrack trackTT((*(kaon.bestTrack())), &(*bFieldHandle));

  TrajectoryStateClosestToPoint theDCAXBS = trackTT.trajectoryStateClosestToPoint( GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()) );  
  
  double DCABS = theDCAXBS.perigeeParameters().transverseImpactParameter();
  double DCABSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
    
  pair<double,double> DCA = make_pair(DCABS,DCABSErr);
    
  return DCA;
}




DEFINE_FWK_MODULE(BToKsteeProducer);
