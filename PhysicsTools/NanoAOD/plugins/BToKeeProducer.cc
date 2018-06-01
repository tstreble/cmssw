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

class BToKeeProducer : public edm::EDProducer {
    
public:
    
    explicit BToKeeProducer(const edm::ParameterSet &iConfig);
    
    ~BToKeeProducer() override {};
    
    
private:
    
    virtual void produce(edm::Event&, const edm::EventSetup&);
    
    bool EEVertexRefitting(const pat::Electron & ele1,
			   const pat::Electron & ele2,
			   edm::ESHandle<TransientTrackBuilder> theTTBuilder,
			   RefCountedKinematicVertex &refitVertex,
			   RefCountedKinematicParticle &refitEE,
			   RefCountedKinematicParticle &refitEle1,
			   RefCountedKinematicParticle &refitEle2);
    
    bool BToKEEVertexRefitting(const pat::Electron &ele1,
			       const pat::Electron &ele2,
			       const pat::PackedCandidate &kaon,
			       edm::ESHandle<TransientTrackBuilder> theTTBuilder,
			       RefCountedKinematicVertex &refitVertex,
			       RefCountedKinematicParticle &refitBToKEE,
			       RefCountedKinematicParticle &refitEle1,
			       RefCountedKinematicParticle &refitEle2,
			       RefCountedKinematicParticle &refitKaon);
    
    pair<double,double> computeLS(RefCountedKinematicVertex refitVertex,
				  reco::BeamSpot beamSpot);
    
    double computeCosAlpha(RefCountedKinematicParticle refitBToKEE,
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
    bool diEleCharge_;
    
    float ElectronMass_ = 0.5109989e-3;
    float ElectronMassErr_ = 3.1*1e-12;
    float KaonMass_ = 0.493677;
    float KaonMassErr_ = 1.6e-5;
    
};



BToKeeProducer::BToKeeProducer(const edm::ParameterSet &iConfig):
beamSpotSrc_( consumes<reco::BeamSpot>                  ( iConfig.getParameter<edm::InputTag>( "beamSpot" ) ) ),
electronSrc_(     consumes<std::vector<pat::Electron>>          ( iConfig.getParameter<edm::InputTag>( "electronCollection" ) ) ),
PFCandSrc_(   consumes<edm::View<pat::PackedCandidate>> ( iConfig.getParameter<edm::InputTag>( "PFCandCollection" ) ) ),
ptMinEle_(     iConfig.getParameter<double>( "ElectronMinPt" ) ),
etaMaxEle_(    iConfig.getParameter<double>( "ElectronMaxEta" ) ),
ptMinKaon_(   iConfig.getParameter<double>( "KaonMinPt" ) ),
etaMaxKaon_(  iConfig.getParameter<double>( "KaonMaxEta" ) ),
DCASigMinKaon_(   iConfig.getParameter<double>( "KaonMinDCASig" ) ),
diEleCharge_( iConfig.getParameter<bool>( "DiElectronChargeCheck" ) )
{
    produces<pat::CompositeCandidateCollection>();
}


void BToKeeProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    
    edm::ESHandle<MagneticField> bFieldHandle;
    edm::ESHandle<TransientTrackBuilder> theTTBuilder;
 
    iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);

    edm::Handle<reco::BeamSpot> beamSpotHandle;
    edm::Handle<reco::VertexCollection> vertexHandle;   

    iEvent.getByToken(beamSpotSrc_, beamSpotHandle);
    
    if ( ! beamSpotHandle.isValid() ) {
        edm::LogError("BToKeeProducer") << "No beam spot available from EventSetup" ;
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
        
        // loop on all the eeK triplets
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
                
                RefCountedKinematicVertex refitVertexEE;
                RefCountedKinematicParticle refitEE;
                RefCountedKinematicParticle refitEle1_EE;
		RefCountedKinematicParticle refitEle2_EE;
                
                bool passed = EEVertexRefitting(ele1, ele2,
						theTTBuilder,
						refitVertexEE,
						refitEE,
						refitEle1_EE,
						refitEle2_EE);
                
                if ( !passed) continue;
                
                pair<double,double> EELS = computeLS(refitVertexEE,beamSpot);
                double EELSBS = EELS.first;
                double EELSBSErr = EELS.second;

                
                double EEVtx_CL = TMath::Prob((double)refitVertexEE->chiSquared(),
                                                int(rint(refitVertexEE->degreesOfFreedom())));
                
                double EE_mass_err = sqrt(refitEE->currentState().kinematicParametersError().matrix()(6,6));

		math::XYZVector refitEle1V3D_EE = refitEle1_EE->refittedTransientTrack().track().momentum();
		math::XYZVector refitEle2V3D_EE = refitEle2_EE->refittedTransientTrack().track().momentum();		
		math::XYZVector refitEEV3D = refitEle1V3D_EE + refitEle2V3D_EE;
              
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
                    
                    RefCountedKinematicVertex refitVertexBToKEE;
                    RefCountedKinematicParticle refitBToKEE;
                    RefCountedKinematicParticle refitEle1;
                    RefCountedKinematicParticle refitEle2;
                    RefCountedKinematicParticle refitKaon;
                    
                    
                    passed = BToKEEVertexRefitting(ele1, ele2, pfCand,
						   theTTBuilder,
						   refitVertexBToKEE,
						   refitBToKEE,
						   refitEle1,
						   refitEle2,
						   refitKaon);
                    
                    if (!passed) continue;
                    
                    pair<double,double> BToKEELS = computeLS(refitVertexBToKEE,beamSpot);
                    double LSBS = BToKEELS.first;
                    double LSBSErr = BToKEELS.second;
                    
                    double BToKEEVtx_CL = TMath::Prob((double)refitVertexBToKEE->chiSquared(),
                                                    int(rint(refitVertexBToKEE->degreesOfFreedom())));

                    
                    double cosAlpha = computeCosAlpha(refitBToKEE,refitVertexBToKEE,beamSpot);
                    
                    double mass_err = sqrt(refitBToKEE->currentState().kinematicParametersError().matrix()(6,6));
                    
                    pat::CompositeCandidate BToKEECand;
                    BToKEECand.addDaughter( ele1 , "ele1");
                    BToKEECand.addDaughter( ele2 , "ele2");
                    BToKEECand.addDaughter( pfCand, "kaon");
                    
		    math::XYZVector refitEle1V3D = refitEle1->refittedTransientTrack().track().momentum();
                    BToKEECand.addUserFloat("ele1_pt",     sqrt(refitEle1V3D.perp2()));
                    BToKEECand.addUserFloat("ele1_eta",    refitEle1V3D.eta());
                    BToKEECand.addUserFloat("ele1_phi",    refitEle1V3D.phi());
                    BToKEECand.addUserFloat("ele1_charge", refitEle1->currentState().particleCharge());

		    math::XYZVector refitEle2V3D = refitEle2->refittedTransientTrack().track().momentum();
                    BToKEECand.addUserFloat("ele2_pt",     sqrt(refitEle2V3D.perp2()));
                    BToKEECand.addUserFloat("ele2_eta",    refitEle2V3D.eta());
                    BToKEECand.addUserFloat("ele2_phi",    refitEle2V3D.phi());
		    BToKEECand.addUserFloat("ele2_charge", refitEle2->currentState().particleCharge());		    

		    math::XYZVector refitKaonV3D = refitKaon->refittedTransientTrack().track().momentum();
                    BToKEECand.addUserFloat("kaon_pt",    sqrt(refitKaonV3D.perp2()));
                    BToKEECand.addUserFloat("kaon_eta",   refitKaonV3D.eta());
                    BToKEECand.addUserFloat("kaon_phi",   refitKaonV3D.phi());
                    BToKEECand.addUserFloat("kaon_charge",refitKaon->currentState().particleCharge());
		    BToKEECand.addUserFloat("kaon_DCASig", DCABS/DCABSErr);

                    BToKEECand.addUserFloat("ee_pt",     sqrt(refitEEV3D.perp2()));
                    BToKEECand.addUserFloat("ee_eta",    refitEEV3D.eta());
                    BToKEECand.addUserFloat("ee_phi",    refitEEV3D.phi());
                    BToKEECand.addUserFloat("ee_mass",   refitEE->currentState().mass());
                    BToKEECand.addUserFloat("ee_mass_err",   EE_mass_err);

                    BToKEECand.addUserFloat("ee_Lxy", (float) EELSBS/EELSBSErr);
                    BToKEECand.addUserFloat("ee_CL_vtx", (float) EEVtx_CL);
                    
		    math::XYZVector refitBToKEEV3D = refitEle1V3D + refitEle2V3D + refitKaonV3D;
                    BToKEECand.addUserFloat("pt",     sqrt(refitBToKEEV3D.perp2()));
                    BToKEECand.addUserFloat("eta",    refitBToKEEV3D.eta());
                    BToKEECand.addUserFloat("phi",    refitBToKEEV3D.phi());
                    BToKEECand.addUserFloat("mass",   refitBToKEE->currentState().mass());
                    BToKEECand.addUserFloat("mass_err", mass_err);

                    BToKEECand.addUserFloat("Lxy", (float) LSBS/LSBSErr);
                    BToKEECand.addUserFloat("CL_vtx", (float) BToKEEVtx_CL);
                    BToKEECand.addUserFloat("cosAlpha", (float) cosAlpha);
                    
                    result->push_back(BToKEECand);
                    
                    
                }
                
            }
            
        }
        
    }
    
    //iEvent.put(result);
    iEvent.put(std::move(result));
    
}



bool BToKeeProducer::EEVertexRefitting(const pat::Electron & ele1,
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




bool BToKeeProducer::BToKEEVertexRefitting(const pat::Electron &ele1,
					   const pat::Electron &ele2,
					   const pat::PackedCandidate &kaon,
					   edm::ESHandle<TransientTrackBuilder> theTTBuilder,					   
					   RefCountedKinematicVertex &refitVertex,
					   RefCountedKinematicParticle &refitBToKEE,
					   RefCountedKinematicParticle &refitEle1,
					   RefCountedKinematicParticle &refitEle2,
					   RefCountedKinematicParticle &refitKaon){

    const reco::TransientTrack ele1TT = theTTBuilder->build(ele1.gsfTrack()); //closestCtfTrackRef ?
    const reco::TransientTrack ele2TT = theTTBuilder->build(ele2.gsfTrack());
    const reco::TransientTrack kaonTT = theTTBuilder->build(kaon.bestTrack());
    
    KinematicParticleFactoryFromTransientTrack partFactory;
    KinematicParticleVertexFitter PartVtxFitter;
    
    std::vector<RefCountedKinematicParticle> BToKEEParticles;
    double chi = 0.;
    double ndf = 0.;
    BToKEEParticles.push_back(partFactory.particle(ele1TT,ElectronMass_,chi,ndf,ElectronMassErr_));
    BToKEEParticles.push_back(partFactory.particle(ele2TT,ElectronMass_,chi,ndf,ElectronMassErr_));
    BToKEEParticles.push_back(partFactory.particle(kaonTT,KaonMass_,chi,ndf,KaonMassErr_));

    RefCountedKinematicTree BToKEEVertexFitTree = PartVtxFitter.fit(BToKEEParticles);
    
    if ( !BToKEEVertexFitTree->isValid()) return false;
    
    BToKEEVertexFitTree->movePointerToTheTop();
    refitVertex = BToKEEVertexFitTree->currentDecayVertex();
    refitBToKEE = BToKEEVertexFitTree->currentParticle();
    
    if ( !refitVertex->vertexIsValid()) return false;
    
    // extract the re-fitted tracks
    BToKEEVertexFitTree->movePointerToTheTop();
    
    BToKEEVertexFitTree->movePointerToTheFirstChild();
    refitEle1 = BToKEEVertexFitTree->currentParticle();
    
    BToKEEVertexFitTree->movePointerToTheNextChild();
    refitEle2 = BToKEEVertexFitTree->currentParticle();
    
    BToKEEVertexFitTree->movePointerToTheNextChild();
    refitKaon = BToKEEVertexFitTree->currentParticle();
    
    return true;



}






pair<double,double> BToKeeProducer::computeLS(RefCountedKinematicVertex refitVertex,
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



double BToKeeProducer::computeCosAlpha(RefCountedKinematicParticle refitBToKEE,
				       RefCountedKinematicVertex refitVertex,
				       reco::BeamSpot beamSpot){
    
  TVector v(2);
  v[0] = refitVertex->position().x()-beamSpot.position().x();
  v[1] = refitVertex->position().y()-beamSpot.position().y();

  TVector w(2);
  w[0] = refitBToKEE->currentState().globalMomentum().x();
  w[1] = refitBToKEE->currentState().globalMomentum().y();
    
  double cosAlpha = v*w/sqrt(v.Norm2Sqr()*w.Norm2Sqr());
  return cosAlpha;
}






pair<double,double> BToKeeProducer::computeDCA(const pat::PackedCandidate &kaon,
					       edm::ESHandle<MagneticField> bFieldHandle,
					       reco::BeamSpot beamSpot){
  
  const reco::TransientTrack trackTT((*(kaon.bestTrack())), &(*bFieldHandle));

  TrajectoryStateClosestToPoint theDCAXBS = trackTT.trajectoryStateClosestToPoint( GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()) );  
  
  double DCABS = theDCAXBS.perigeeParameters().transverseImpactParameter();
  double DCABSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
    
  pair<double,double> DCA = make_pair(DCABS,DCABSErr);
    
  return DCA;
}




DEFINE_FWK_MODULE(BToKeeProducer);
