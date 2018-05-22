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

class BToKmumuProducer : public edm::EDProducer {
    
public:
    
    explicit BToKmumuProducer(const edm::ParameterSet &iConfig);
    
    ~BToKmumuProducer() override {};
    
    
private:
    
    virtual void produce(edm::Event&, const edm::EventSetup&);
    
    bool MuMuVertexRefitting(const pat::Muon & muon1,
                             const pat::Muon & muon2,
                             edm::ESHandle<MagneticField> bFieldHandle,
                             RefCountedKinematicVertex &refitVertex,
                             RefCountedKinematicParticle &refitMuMu,
			     RefCountedKinematicParticle &refitMu1,
			     RefCountedKinematicParticle &refitMu2);
    
    bool BToKMuMuVertexRefitting(const pat::Muon &muon1,
                                 const pat::Muon &muon2,
                                 const pat::PackedCandidate &kaon,
                                 edm::ESHandle<MagneticField> bFieldHandle,
                                 RefCountedKinematicVertex &refitVertex,
                                 RefCountedKinematicParticle &refitBToKMuMu,
                                 RefCountedKinematicParticle &refitMu1,
                                 RefCountedKinematicParticle &refitMu2,
                                 RefCountedKinematicParticle &refitKaon);
    
    pair<double,double> computeLS(RefCountedKinematicVertex refitVertex,
				  reco::BeamSpot beamSpot);
    
    double computeCosAlpha(RefCountedKinematicParticle refitBToKMuMu,
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
    bool diMuonCharge_;
    
    float MuonMass_ = 0.10565837;
    float MuonMassErr_ = 3.5*1e-9;
    float KaonMass_ = 0.493677;
    float KaonMassErr_ = 1.6e-5;
    
};



BToKmumuProducer::BToKmumuProducer(const edm::ParameterSet &iConfig):
beamSpotSrc_( consumes<reco::BeamSpot>                  ( iConfig.getParameter<edm::InputTag>( "beamSpot" ) ) ),
vertexSrc_( consumes<reco::VertexCollection>            ( iConfig.getParameter<edm::InputTag>( "vertexCollection" ) ) ),
muonSrc_(     consumes<std::vector<pat::Muon>>          ( iConfig.getParameter<edm::InputTag>( "muonCollection" ) ) ),
PFCandSrc_(   consumes<edm::View<pat::PackedCandidate>> ( iConfig.getParameter<edm::InputTag>( "PFCandCollection" ) ) ),
ptMinMu_(     iConfig.getParameter<double>( "MuonMinPt" ) ),
etaMaxMu_(    iConfig.getParameter<double>( "MuonMaxEta" ) ),
ptMinKaon_(   iConfig.getParameter<double>( "KaonMinPt" ) ),
etaMaxKaon_(  iConfig.getParameter<double>( "KaonMaxEta" ) ),
DCASigMinKaon_(   iConfig.getParameter<double>( "KaonMinDCASig" ) ),
diMuonCharge_( iConfig.getParameter<bool>( "DiMuonChargeCheck" ) )
{
    produces<pat::CompositeCandidateCollection>();
}


void BToKmumuProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    
    edm::ESHandle<MagneticField> bFieldHandle;
    edm::Handle<reco::BeamSpot> beamSpotHandle;
    edm::Handle<reco::VertexCollection> vertexHandle;
    
    iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

    iEvent.getByToken(beamSpotSrc_, beamSpotHandle);
    
    if ( ! beamSpotHandle.isValid() ) {
        edm::LogError("BToKMuMuProducer") << "No beam spot available from EventSetup" ;
    }
    
    reco::BeamSpot beamSpot = *beamSpotHandle;
    
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
        
        // loop on all the mumuK triplets
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
                
                RefCountedKinematicVertex refitVertexMuMu;
                RefCountedKinematicParticle refitMuMu;
                RefCountedKinematicParticle refitMu1_MuMu;
		RefCountedKinematicParticle refitMu2_MuMu;
                
                bool passed = MuMuVertexRefitting(muon1, muon2,
                                                  bFieldHandle,
                                                  refitVertexMuMu,
                                                  refitMuMu,
						  refitMu1_MuMu,
						  refitMu2_MuMu);
                
                if ( !passed) continue;
                
                pair<double,double> MuMuLS = computeLS(refitVertexMuMu,beamSpot);
                double MuMuLSBS = MuMuLS.first;
                double MuMuLSBSErr = MuMuLS.second;

                
                double MuMuVtx_CL = TMath::Prob((double)refitVertexMuMu->chiSquared(),
                                                int(rint(refitVertexMuMu->degreesOfFreedom())));
                
                double MuMu_mass_err = sqrt(refitMuMu->currentState().kinematicParametersError().matrix()(6,6));

		math::XYZVector refitMu1V3D_MuMu = refitMu1_MuMu->refittedTransientTrack().track().momentum();
		math::XYZVector refitMu2V3D_MuMu = refitMu2_MuMu->refittedTransientTrack().track().momentum();		
		math::XYZVector refitMuMuV3D = refitMu1V3D_MuMu + refitMu2V3D_MuMu;
              
                //Kaon
                for (unsigned int k = 0; k < pfCandNumber; ++k) {
                    
                    const pat::PackedCandidate & pfCand = (*pfCandHandle)[k];
                    if(abs(pfCand.pdgId())!=211) continue; //Charged hadrons
                    if(!pfCand.hasTrackDetails()) continue;
                    if(pfCand.pt()<ptMinKaon_ || abs(pfCand.eta())>etaMaxKaon_) continue;

		    pair<double,double> DCA = computeDCA(pfCand,
							 bFieldHandle,
							 beamSpot);
		    double DCABS = DCA.first;
		    double DCABSErr = DCA.second;

		    if(DCABS/DCABSErr<DCASigMinKaon_) continue;
                    
                    RefCountedKinematicVertex refitVertexBToKMuMu;
                    RefCountedKinematicParticle refitBToKMuMu;
                    RefCountedKinematicParticle refitMuon1;
                    RefCountedKinematicParticle refitMuon2;
                    RefCountedKinematicParticle refitKaon;
                    
                    
                    passed = BToKMuMuVertexRefitting(muon1, muon2, pfCand,
                                                     bFieldHandle,
                                                     refitVertexBToKMuMu,
                                                     refitBToKMuMu,
                                                     refitMuon1,
                                                     refitMuon2,
                                                     refitKaon);
                    
                    if (!passed) continue;
                    
                    pair<double,double> BToKMuMuLS = computeLS(refitVertexBToKMuMu,beamSpot);
                    double LSBS = BToKMuMuLS.first;
                    double LSBSErr = BToKMuMuLS.second;
                    
                    double BToKMuMuVtx_CL = TMath::Prob((double)refitVertexBToKMuMu->chiSquared(),
                                                    int(rint(refitVertexBToKMuMu->degreesOfFreedom())));

                    
                    double cosAlpha = computeCosAlpha(refitBToKMuMu,refitVertexBToKMuMu,beamSpot);
                    
                    double mass_err = sqrt(refitBToKMuMu->currentState().kinematicParametersError().matrix()(6,6));
                    
                    pat::CompositeCandidate BToKMuMuCand;
                    BToKMuMuCand.addDaughter( muon1 , "muon1");
                    BToKMuMuCand.addDaughter( muon2 , "muon2");
                    BToKMuMuCand.addDaughter( pfCand, "kaon");
                    
		    math::XYZVector refitMu1V3D = refitMuon1->refittedTransientTrack().track().momentum();
                    BToKMuMuCand.addUserFloat("mu1_pt",     sqrt(refitMu1V3D.perp2()));
                    BToKMuMuCand.addUserFloat("mu1_eta",    refitMu1V3D.eta());
                    BToKMuMuCand.addUserFloat("mu1_phi",    refitMu1V3D.phi());
                    BToKMuMuCand.addUserFloat("mu1_charge", refitMuon1->currentState().particleCharge());

		    math::XYZVector refitMu2V3D = refitMuon2->refittedTransientTrack().track().momentum();
                    BToKMuMuCand.addUserFloat("mu2_pt",     sqrt(refitMu2V3D.perp2()));
                    BToKMuMuCand.addUserFloat("mu2_eta",    refitMu2V3D.eta());
                    BToKMuMuCand.addUserFloat("mu2_phi",    refitMu2V3D.phi());
		    BToKMuMuCand.addUserFloat("mu2_charge", refitMuon2->currentState().particleCharge());		    

		    math::XYZVector refitKaonV3D = refitKaon->refittedTransientTrack().track().momentum();
                    BToKMuMuCand.addUserFloat("kaon_pt",    sqrt(refitKaonV3D.perp2()));
                    BToKMuMuCand.addUserFloat("kaon_eta",   refitKaonV3D.eta());
                    BToKMuMuCand.addUserFloat("kaon_phi",   refitKaonV3D.phi());
                    BToKMuMuCand.addUserFloat("kaon_charge",refitKaon->currentState().particleCharge());
		    BToKMuMuCand.addUserFloat("kaon_DCASig", DCABS/DCABSErr);
                    		    
                    BToKMuMuCand.addUserFloat("mumu_pt",     sqrt(refitMuMuV3D.perp2()));
                    BToKMuMuCand.addUserFloat("mumu_eta",    refitMuMuV3D.eta());
                    BToKMuMuCand.addUserFloat("mumu_phi",    refitMuMuV3D.phi());
                    BToKMuMuCand.addUserFloat("mumu_mass",   refitMuMu->currentState().mass());
                    BToKMuMuCand.addUserFloat("mumu_mass_err",   MuMu_mass_err);
                    
                    BToKMuMuCand.addUserFloat("mumu_Lxy", (float) MuMuLSBS/MuMuLSBSErr);
                    BToKMuMuCand.addUserFloat("mumu_CL_vtx", (float) MuMuVtx_CL);
                    
		    math::XYZVector refitBToKMuMuV3D = refitMu1V3D + refitMu2V3D + refitKaonV3D;
                    BToKMuMuCand.addUserFloat("pt",     sqrt(refitBToKMuMuV3D.perp2()));
                    BToKMuMuCand.addUserFloat("eta",    refitBToKMuMuV3D.eta());
                    BToKMuMuCand.addUserFloat("phi",    refitBToKMuMuV3D.phi());
                    BToKMuMuCand.addUserFloat("mass",   refitBToKMuMu->currentState().mass());
                    BToKMuMuCand.addUserFloat("mass_err", mass_err);
                    
                    BToKMuMuCand.addUserFloat("Lxy", (float) LSBS/LSBSErr);
                    BToKMuMuCand.addUserFloat("CL_vtx", (float) BToKMuMuVtx_CL);
                    BToKMuMuCand.addUserFloat("cosAlpha", (float) cosAlpha);
                    
                    result->push_back(BToKMuMuCand);
                    
                    
                }
                
            }
            
        }
        
    }
    
    //iEvent.put(result);
    iEvent.put(std::move(result));
    
}



bool BToKmumuProducer::MuMuVertexRefitting(const pat::Muon & muon1,
                                           const pat::Muon & muon2,
                                           edm::ESHandle<MagneticField> bFieldHandle,
                                           RefCountedKinematicVertex &refitVertex,
                                           RefCountedKinematicParticle &refitMuMu,
					   RefCountedKinematicParticle &refitMu1,
					   RefCountedKinematicParticle &refitMu2){
    
    const reco::TransientTrack muon1TT(muon1.innerTrack(), &(*bFieldHandle));
    const reco::TransientTrack muon2TT(muon2.innerTrack(), &(*bFieldHandle));
    
    KinematicParticleFactoryFromTransientTrack partFactory;
    KinematicParticleVertexFitter PartVtxFitter;
    
    std::vector<RefCountedKinematicParticle> muonParticles;
    double chi = 0.;
    double ndf = 0.;
    muonParticles.push_back(partFactory.particle(muon1TT,MuonMass_,chi,ndf,MuonMassErr_));
    muonParticles.push_back(partFactory.particle(muon2TT,MuonMass_,chi,ndf,MuonMassErr_));
    RefCountedKinematicTree mumuVertexFitTree = PartVtxFitter.fit(muonParticles);
    
    if ( !mumuVertexFitTree->isValid()) return false;
    
    mumuVertexFitTree->movePointerToTheTop();
    refitVertex = mumuVertexFitTree->currentDecayVertex();
    refitMuMu = mumuVertexFitTree->currentParticle();
    
    if ( !refitVertex->vertexIsValid()) return false;

    // extract the re-fitted tracks
    mumuVertexFitTree->movePointerToTheTop();
    
    mumuVertexFitTree->movePointerToTheFirstChild();
    refitMu1 = mumuVertexFitTree->currentParticle();
    
    mumuVertexFitTree->movePointerToTheNextChild();
    refitMu2 = mumuVertexFitTree->currentParticle();
    
    return true;
    
}




bool BToKmumuProducer::BToKMuMuVertexRefitting(const pat::Muon &muon1,
                                               const pat::Muon &muon2,
                                               const pat::PackedCandidate &kaon,
                                               edm::ESHandle<MagneticField> bFieldHandle,
                                               RefCountedKinematicVertex &refitVertex,
                                               RefCountedKinematicParticle &refitBToKMuMu,
                                               RefCountedKinematicParticle &refitMu1,
                                               RefCountedKinematicParticle &refitMu2,
                                               RefCountedKinematicParticle &refitKaon){

    const reco::TransientTrack muon1TT(muon1.innerTrack(), &(*bFieldHandle));
    const reco::TransientTrack muon2TT(muon2.innerTrack(), &(*bFieldHandle));
    const reco::TransientTrack kaonTT((*(kaon.bestTrack())), &(*bFieldHandle));
    
    KinematicParticleFactoryFromTransientTrack partFactory;
    KinematicParticleVertexFitter PartVtxFitter;
    
    std::vector<RefCountedKinematicParticle> BToKMuMuParticles;
    double chi = 0.;
    double ndf = 0.;
    BToKMuMuParticles.push_back(partFactory.particle(muon1TT,MuonMass_,chi,ndf,MuonMassErr_));
    BToKMuMuParticles.push_back(partFactory.particle(muon2TT,MuonMass_,chi,ndf,MuonMassErr_));
    BToKMuMuParticles.push_back(partFactory.particle(kaonTT,KaonMass_,chi,ndf,KaonMassErr_));

    RefCountedKinematicTree BToKMuMuVertexFitTree = PartVtxFitter.fit(BToKMuMuParticles);
    
    if ( !BToKMuMuVertexFitTree->isValid()) return false;
    
    BToKMuMuVertexFitTree->movePointerToTheTop();
    refitVertex = BToKMuMuVertexFitTree->currentDecayVertex();
    refitBToKMuMu = BToKMuMuVertexFitTree->currentParticle();
    
    if ( !refitVertex->vertexIsValid()) return false;
    
    // extract the re-fitted tracks
    BToKMuMuVertexFitTree->movePointerToTheTop();
    
    BToKMuMuVertexFitTree->movePointerToTheFirstChild();
    refitMu1 = BToKMuMuVertexFitTree->currentParticle();
    
    BToKMuMuVertexFitTree->movePointerToTheNextChild();
    refitMu2 = BToKMuMuVertexFitTree->currentParticle();
    
    BToKMuMuVertexFitTree->movePointerToTheNextChild();
    refitKaon = BToKMuMuVertexFitTree->currentParticle();
    
    return true;



}






pair<double,double> BToKmumuProducer::computeLS(RefCountedKinematicVertex refitVertex,
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



double BToKmumuProducer::computeCosAlpha(RefCountedKinematicParticle refitBToKMuMu,
                                         RefCountedKinematicVertex refitVertex,
                                         reco::BeamSpot beamSpot){
    
  TVector v(2);
  v[0] = refitVertex->position().x()-beamSpot.position().x();
  v[1] = refitVertex->position().y()-beamSpot.position().y();

  TVector w(2);
  w[0] = refitBToKMuMu->currentState().globalMomentum().x();
  w[1] = refitBToKMuMu->currentState().globalMomentum().y();
    
  double cosAlpha = v*w/sqrt(v.Norm2Sqr()*w.Norm2Sqr());
  return cosAlpha;
}






pair<double,double> BToKmumuProducer::computeDCA(const pat::PackedCandidate &kaon,
						 edm::ESHandle<MagneticField> bFieldHandle,
						 reco::BeamSpot beamSpot){

  const reco::TransientTrack trackTT((*(kaon.bestTrack())), &(*bFieldHandle));

  TrajectoryStateClosestToPoint theDCAXBS = trackTT.trajectoryStateClosestToPoint( GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()) );  
  
  double DCABS = theDCAXBS.perigeeParameters().transverseImpactParameter();
  double DCABSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
    
  pair<double,double> DCA = make_pair(DCABS,DCABSErr);
    
  return DCA;
}




DEFINE_FWK_MODULE(BToKmumuProducer);
